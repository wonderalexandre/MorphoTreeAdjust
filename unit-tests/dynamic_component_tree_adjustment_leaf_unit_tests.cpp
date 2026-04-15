#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../morphoTreeAdjust/include/DynamicComponentTreeAdjustmentLeaf.hpp"

#include "../dev-tools/external/stb/stb_image.h"

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr load_gray_image(const std::string &path) {
    int numCols = 0;
    int numRows = 0;
    int numChannels = 0;
    uint8_t *data = stbi_load(path.c_str(), &numCols, &numRows, &numChannels, 1);
    if (data == nullptr) {
        throw std::runtime_error("failed to load image: " + path);
    }

    ImageUInt8Ptr image = ImageUInt8::create(numRows, numCols);
    std::copy(data, data + (numRows * numCols), image->rawData());
    stbi_image_free(data);
    return image;
}

ImageUInt8Ptr make_demo_image() {
    auto input = ImageUInt8::create(4, 4);
    const std::vector<uint8_t> values = {
        4, 4, 2, 1,
        4, 3, 2, 1,
        5, 5, 2, 0,
        5, 6, 6, 0,
    };
    for (size_t i = 0; i < values.size(); ++i) {
        (*input)[static_cast<int>(i)] = values[i];
    }
    return input;
}

std::string image_to_string(ImageUInt8Ptr image) {
    std::ostringstream out;
    out << image->getNumRows() << "x" << image->getNumCols() << ":";
    for (int row = 0; row < image->getNumRows(); ++row) {
        out << " [";
        for (int col = 0; col < image->getNumCols(); ++col) {
            if (col != 0) {
                out << ' ';
            }
            out << static_cast<int>((*image)[row * image->getNumCols() + col]);
        }
        out << "]";
    }
    return out.str();
}

ImageUInt8Ptr make_random_image(int numRows,
                                int numCols,
                                std::mt19937 &rng,
                                int maxValueInclusive) {
    auto image = ImageUInt8::create(numRows, numCols);
    for (int pixelId = 0; pixelId < numRows * numCols; ++pixelId) {
        (*image)[pixelId] = static_cast<uint8_t>(rng() % static_cast<std::mt19937::result_type>(maxValueInclusive + 1));
    }
    return image;
}

int deterministic_int(std::mt19937 &rng, int minValue, int maxValue) {
    const auto span = static_cast<std::mt19937::result_type>(maxValue - minValue + 1);
    return minValue + static_cast<int>(rng() % span);
}

std::vector<NodeId> getNodesToPrune(const DynamicComponentTree &tree,
                                    const std::vector<float> &attribute,
                                    int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree.getRoot());
    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (attribute[static_cast<size_t>(nodeId)] > threshold) {
            for (NodeId childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }
    return out;
}

void collectPostOrderNodes(const DynamicComponentTree &tree, NodeId nodeId, std::vector<NodeId> &out) {
    for (NodeId childId : tree.getChildren(nodeId)) {
        collectPostOrderNodes(tree, childId, out);
    }
    out.push_back(nodeId);
}

NodeId firstNonRootLeaf(const DynamicComponentTree &tree) {
    std::vector<NodeId> postOrderNodes;
    collectPostOrderNodes(tree, tree.getRoot(), postOrderNodes);
    for (NodeId nodeId : postOrderNodes) {
        if (nodeId != tree.getRoot() && tree.isAlive(nodeId) && tree.isLeaf(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

void require_tree_has_no_empty_live_nodes(const DynamicComponentTree &tree,
                                          const std::string &treeName) {
    require(tree.isNode(tree.getRoot()) && tree.isAlive(tree.getRoot()),
            treeName + ": the root must remain alive");
    for (NodeId nodeId = 0; nodeId < tree.getGlobalIdSpaceSize(); ++nodeId) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        require(tree.getNumProperParts(nodeId) > 0,
                treeName + ": live node without proper parts at id " + std::to_string(nodeId));
        if (!tree.isRoot(nodeId)) {
            const NodeId parentId = tree.getNodeParent(nodeId);
            require(parentId != InvalidNode && tree.isNode(parentId) && tree.isAlive(parentId),
                    treeName + ": live node with invalid parent at id " + std::to_string(nodeId));
            require(tree.hasChild(parentId, nodeId),
                    treeName + ": parent-child inconsistency at id " + std::to_string(nodeId));
        }
    }
}

void require_trees_stay_consistent(const DynamicComponentTree &minTree,
                                   const DynamicComponentTree &maxTree,
                                   const std::string &stepLabel) {
    require_tree_has_no_empty_live_nodes(minTree, stepLabel + " / minTree");
    require_tree_has_no_empty_live_nodes(maxTree, stepLabel + " / maxTree");
    require(minTree.reconstructionImage()->isEqual(maxTree.reconstructionImage()),
            stepLabel + ": min-tree and max-tree must reconstruct the same image");
}

ImageUInt8Ptr applyNaiveLeafThreshold(ImageUInt8Ptr input,
                                      AdjacencyRelationPtr adj,
                                      int threshold) {
    DynamicComponentTree maxTree(input, true, adj);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    std::vector<float> maxArea = maxAreaComputer.compute();
    auto nodesToPrune = getNodesToPrune(maxTree, maxArea, threshold);
    std::vector<NodeId> postOrderNodes;
    for (NodeId nodeId : nodesToPrune) {
        if (nodeId == maxTree.getRoot()) {
            continue;
        }
        postOrderNodes.clear();
        collectPostOrderNodes(maxTree, nodeId, postOrderNodes);
        for (NodeId leafId : postOrderNodes) {
            if (leafId != maxTree.getRoot() && maxTree.isAlive(leafId) && maxTree.isLeaf(leafId)) {
                maxTree.pruneNode(leafId);
            }
        }
    }

    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adj);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> minArea = minAreaComputer.compute();
    nodesToPrune = getNodesToPrune(minTree, minArea, threshold);
    for (NodeId nodeId : nodesToPrune) {
        if (nodeId == minTree.getRoot()) {
            continue;
        }
        postOrderNodes.clear();
        collectPostOrderNodes(minTree, nodeId, postOrderNodes);
        for (NodeId leafId : postOrderNodes) {
            if (leafId != minTree.getRoot() && minTree.isAlive(leafId) && minTree.isLeaf(leafId)) {
                minTree.pruneNode(leafId);
            }
        }
    }

    return minTree.reconstructionImage();
}

void testDynamicLeafAdjustmentMatchesNaiveLeafBaseline() {
    constexpr int kThreshold = 4;
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DynamicComponentTreeAdjustmentLeaf<AltitudeType> adjust(&minTree, &maxTree, *adj);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> maxArea(static_cast<size_t>(maxTree.getGlobalIdSpaceSize()), 0.0f);
    std::vector<float> minArea(static_cast<size_t>(minTree.getGlobalIdSpaceSize()), 0.0f);
    maxAreaComputer.compute(std::span<float>(maxArea));
    minAreaComputer.compute(std::span<float>(minArea));
    adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                std::span<float>(minArea), std::span<float>(maxArea));

    auto nodesToPrune = getNodesToPrune(maxTree, maxArea, kThreshold);
    require(!nodesToPrune.empty(), "expected candidates in the max-tree");
    adjust.pruneMaxTreeAndUpdateMinTree(nodesToPrune);

    nodesToPrune = getNodesToPrune(minTree, minArea, kThreshold);
    adjust.pruneMinTreeAndUpdateMaxTree(nodesToPrune);

    const auto dynamicImage = minTree.reconstructionImage();
    const auto baselineImage = applyNaiveLeafThreshold(input, adj, kThreshold);
    require(dynamicImage->isEqual(baselineImage), "the dynamic leaf result must match the naive leaf baseline");
}

void testDirectLeafAdjustOverloadsKeepTreesConsistent() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DynamicComponentTreeAdjustmentLeaf<AltitudeType> adjust(&minTree, &maxTree, *adj);

    NodeId maxLeafId = firstNonRootLeaf(maxTree);
    require(maxLeafId != InvalidNode, "expected a non-root leaf in the max-tree");
    adjust.pruneMaxTreeAndUpdateMinTree(maxLeafId);
    require(minTree.reconstructionImage()->isEqual(maxTree.reconstructionImage()),
            "the trees must remain consistent after leaf-based adjustMinTree");

    NodeId minLeafId = firstNonRootLeaf(minTree);
    require(minLeafId != InvalidNode, "expected a non-root leaf in the min-tree");
    adjust.pruneMinTreeAndUpdateMaxTree(minLeafId);
    require(minTree.reconstructionImage()->isEqual(maxTree.reconstructionImage()),
            "the trees must remain consistent after leaf-based adjustMaxTree");
}

void testDynamicLeafRandomStressMatchesNaiveLeaf() {
    constexpr int kNumCases = 80;
    std::mt19937 rng(789123u);
    const std::array<double, 2> radii = {1.0, 1.5};

    for (int caseIndex = 0; caseIndex < kNumCases; ++caseIndex) {
        const int numRows = deterministic_int(rng, 2, 6);
        const int numCols = deterministic_int(rng, 2, 6);
        auto image = make_random_image(numRows, numCols, rng, 7);
        const double radius = radii[static_cast<size_t>(deterministic_int(rng, 0, 1))];
        auto adj = std::make_shared<AdjacencyRelation>(numRows, numCols, radius);

        DynamicComponentTree maxTree(image, true, adj);
        DynamicComponentTree minTree(image, false, adj);
        DynamicComponentTreeAdjustmentLeaf<AltitudeType> adjust(&minTree, &maxTree, *adj);
        DynamicAreaComputer maxAreaComputer(&maxTree);
        DynamicAreaComputer minAreaComputer(&minTree);
        std::vector<float> maxArea(maxTree.getGlobalIdSpaceSize(), 0.0f);
        std::vector<float> minArea(minTree.getGlobalIdSpaceSize(), 0.0f);
        maxAreaComputer.compute(std::span<float>(maxArea));
        minAreaComputer.compute(std::span<float>(minArea));
        adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                    std::span<float>(minArea), std::span<float>(maxArea));

        std::vector<int> thresholds;
        const int maxThreshold = std::min(numRows * numCols, 6);
        for (int threshold = 1; threshold <= maxThreshold; threshold += 2) {
            thresholds.push_back(threshold);
        }

        for (int threshold : thresholds) {
            auto nodesToPrune = getNodesToPrune(maxTree, maxArea, threshold);
            adjust.pruneMaxTreeAndUpdateMinTree(nodesToPrune);

            std::ostringstream phase1Label;
            phase1Label << "leaf stress case " << caseIndex << " threshold " << threshold << " phase max->min";
            require_trees_stay_consistent(minTree, maxTree, phase1Label.str());

            nodesToPrune = getNodesToPrune(minTree, minArea, threshold);
            adjust.pruneMinTreeAndUpdateMaxTree(nodesToPrune);

            std::ostringstream phase2Label;
            phase2Label << "leaf stress case " << caseIndex << " threshold " << threshold << " phase min->max";
            require_trees_stay_consistent(minTree, maxTree, phase2Label.str());

            auto baselineImage = applyNaiveLeafThreshold(image, adj, threshold);
            require(minTree.reconstructionImage()->isEqual(baselineImage),
                    "leaf stress must match the naive leaf baseline; "
                    "case=" + std::to_string(caseIndex) +
                    " threshold=" + std::to_string(threshold) +
                    " radius=" + std::to_string(radius) +
                    " image=" + image_to_string(image));
            image = baselineImage;
            adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radius);
        }
    }
}

void testDynamicLeafHouseThreshold50Regression() {
    constexpr int kThreshold = 50;
    const auto housePath =
        (std::filesystem::path(__FILE__).parent_path().parent_path() / "dat/misc256/house.png").string();
    auto image = load_gray_image(housePath);
    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);

    DynamicComponentTree maxTree(image, true, adj);
    DynamicComponentTree minTree(image, false, adj);
    DynamicComponentTreeAdjustmentLeaf<AltitudeType> adjust(&minTree, &maxTree, *adj);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> maxArea(maxTree.getGlobalIdSpaceSize(), 0.0f);
    std::vector<float> minArea(minTree.getGlobalIdSpaceSize(), 0.0f);
    maxAreaComputer.compute(std::span<float>(maxArea));
    minAreaComputer.compute(std::span<float>(minArea));
    adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                std::span<float>(minArea), std::span<float>(maxArea));

    auto nodesToPrune = getNodesToPrune(maxTree, maxArea, kThreshold);
    adjust.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
    require_trees_stay_consistent(minTree, maxTree, "house threshold 50 apos max->min");

    nodesToPrune = getNodesToPrune(minTree, minArea, kThreshold);
    adjust.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
    require_trees_stay_consistent(minTree, maxTree, "house threshold 50 apos min->max");

    const auto dynamicImage = minTree.reconstructionImage();
    const auto baselineImage = applyNaiveLeafThreshold(image, adj, kThreshold);
    require(dynamicImage->isEqual(baselineImage),
            "dynamic-leaf must match the naive leaf baseline on house threshold 50");
}

} // namespace

int main() {
    try {
        testDynamicLeafAdjustmentMatchesNaiveLeafBaseline();
        testDirectLeafAdjustOverloadsKeepTreesConsistent();
        testDynamicLeafRandomStressMatchesNaiveLeaf();
        testDynamicLeafHouseThreshold50Regression();
    } catch (const std::exception &e) {
        std::cerr << "dynamic_component_tree_adjustment_leaf_unit_tests: " << e.what() << "\n";
        return 1;
    }

    std::cout << "dynamic_component_tree_adjustment_leaf_unit_tests: OK\n";
    return 0;
}
