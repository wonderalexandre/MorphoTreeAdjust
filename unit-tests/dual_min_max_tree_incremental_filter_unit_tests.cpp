#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#define PRINT_LOG 1
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilter.hpp"

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
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
        (*input)[(int) i] = values[i];
    }
    return input;
}

ImageUInt8Ptr make_large_demo_image() {
    auto image = ImageUInt8::create(12, 12);
    const std::vector<uint8_t> raw = {
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
        2, 5, 5, 5, 3, 3, 3, 2, 2, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 5, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 6, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 5, 2, 3, 1,
        2, 5, 5, 5, 3, 3, 3, 2, 2, 2, 3, 1,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1,
        2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 1,
        6, 6, 6, 6, 6, 2, 4, 4, 4, 4, 4, 1,
        7, 3, 6, 4, 6, 2, 4, 8, 4, 8, 4, 1,
        7, 8, 6, 6, 6, 2, 4, 4, 4, 4, 4, 1,
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
    };
    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[(int) i] = raw[i];
    }
    return image;
}

ImageUInt8Ptr loadPgmImage(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    require(in.good(), "failed to open PGM image: " + path);

    auto readToken = [&]() -> std::string {
        std::string token;
        while (in >> token) {
            if (!token.empty() && token[0] == '#') {
                std::string comment;
                std::getline(in, comment);
                continue;
            }
            return token;
        }
        return {};
    };

    const std::string magic = readToken();
    require(magic == "P2" || magic == "P5", "unsupported PGM format: " + magic);

    const int numCols = std::stoi(readToken());
    const int numRows = std::stoi(readToken());
    const int maxValue = std::stoi(readToken());
    require(maxValue == 255, "the test expects PGM with maxValue=255");

    auto image = ImageUInt8::create(numRows, numCols);
    if (magic == "P2") {
        for (int i = 0; i < numRows * numCols; ++i) {
            (*image)[i] = static_cast<uint8_t>(std::stoi(readToken()));
        }
        return image;
    }

    in.get();
    in.read(reinterpret_cast<char *>(image->rawData()), static_cast<std::streamsize>(numRows * numCols));
    require(in.good(), "failed to read PGM raster: " + path);
    return image;
}

NodeId firstNonRootNode(const DynamicComponentTree &tree) {
    for (int nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isRoot(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

NodeId firstStrictNonRootNode(const DynamicComponentTree &tree) {
    for (int nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (tree.isRoot(nodeId)) {
            continue;
        }

        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
            continue;
        }

        const int parentLevel = tree.getAltitude(parentId);
        const int nodeLevel = tree.getAltitude(nodeId);
        if ((tree.isMaxtree() && parentLevel < nodeLevel) ||
            (!tree.isMaxtree() && parentLevel > nodeLevel)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

template <typename Fn>
void requireThrows(Fn &&fn, const std::string &message) {
    bool threw = false;
    try {
        fn();
    } catch (const std::exception &) {
        threw = true;
    }
    require(threw, message);
}

std::vector<NodeId> getNodesToPrune(const DynamicComponentTree &tree,
                                    const std::vector<float> &attribute,
                                    int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree.getRoot());
    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (attribute[(size_t) nodeId] > threshold) {
            for (NodeId childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }
    return out;
}

std::vector<NodeId> getNodesToPrune(DynamicComponentTree *tree,
                                    const std::vector<float> &attribute,
                                    int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree->getRoot());
    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (attribute[(size_t) nodeId] > threshold) {
            for (NodeId childId : tree->getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }
    return out;
}

ImageUInt8Ptr applyNaiveThreshold(ImageUInt8Ptr input,
                                  AdjacencyRelationPtr adj,
                                  int threshold) {
    DynamicComponentTree maxTree(input, true, adj);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    std::vector<float> maxArea = maxAreaComputer.compute();
    for (NodeId nodeId : getNodesToPrune(&maxTree, maxArea, threshold)) {
        if (nodeId != maxTree.getRoot()) {
            maxTree.pruneNode(nodeId);
        }
    }

    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adj);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> minArea = minAreaComputer.compute();
    for (NodeId nodeId : getNodesToPrune(&minTree, minArea, threshold)) {
        if (nodeId != minTree.getRoot()) {
            minTree.pruneNode(nodeId);
        }
    }
    return minTree.reconstructionImage();
}

int count_alive_nodes(const DynamicComponentTree &tree) {
    int count = 0;
    for (int nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId)) {
            ++count;
        }
    }
    return count;
}

bool has_unreachable_alive_nodes(const DynamicComponentTree &tree) {
    std::vector<bool> reachable((size_t) tree.getNumInternalNodeSlots(), false);
    if (tree.getRoot() != InvalidNode && tree.isAlive(tree.getRoot())) {
        for (int nodeId : tree.getNodeSubtree(tree.getRoot())) {
            reachable[(size_t) nodeId] = true;
        }
    }

    for (int nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId) && !reachable[(size_t) nodeId]) {
            return true;
        }
    }
    return false;
}

void require_tree_consistency(const DynamicComponentTree &tree) {
    require(tree.getRoot() != InvalidNode, "root must be valid");
    require(tree.isAlive(tree.getRoot()), "root must be alive");
    require(!has_unreachable_alive_nodes(tree), "all alive nodes must be reachable from root");

    int totalProperParts = 0;
    std::vector<bool> seen((size_t) tree.getNumTotalProperParts(), false);
    for (int nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }

        for (int childId : tree.getChildren(nodeId)) {
            require(tree.isAlive(childId), "child must be alive");
            require(tree.getNodeParent(childId) == nodeId, "child parent must match");
            require(tree.hasChild(nodeId, childId), "parent-child relation must be symmetric");
        }

        require(tree.getNumProperParts(nodeId) > 0, "alive nodes must keep at least one proper part");
        int countedProperParts = 0;
        for (int pixelId : tree.getProperParts(nodeId)) {
            require(pixelId >= 0 && pixelId < tree.getNumTotalProperParts(), "proper part must be in bounds");
            require(tree.getSmallestComponent(pixelId) == nodeId, "proper part owner mismatch");
            require(!seen[(size_t) pixelId], "proper part visited twice");
            seen[(size_t) pixelId] = true;
            ++countedProperParts;
        }
        require(countedProperParts == tree.getNumProperParts(nodeId), "proper part count mismatch");
        totalProperParts += countedProperParts;
    }

    require(totalProperParts == tree.getNumTotalProperParts(), "global proper part count mismatch");
    require(std::all_of(seen.begin(), seen.end(), [](bool v) { return v; }), "proper part partition is incomplete");
}

class RecordingAttributeComputer final : public DynamicAttributeComputer {
private:
    DynamicComponentTree *tree_ = nullptr;

    void computeNode(NodeId nodeId, std::span<float> buffer) const {
        preProcessing(nodeId, buffer);
        for (NodeId childId : tree_->getChildren(nodeId)) {
            computeNode(childId, buffer);
            mergeProcessing(nodeId, childId, buffer);
        }
        postProcessing(nodeId, buffer);
    }

public:
    mutable int moveProperPartsCount = 0;
    mutable int nodeRemovedCount = 0;
    mutable int postProcessingCount = 0;
    mutable NodeId lastMoveTarget = InvalidNode;
    mutable NodeId lastMoveSource = InvalidNode;
    mutable NodeId lastRemovedNode = InvalidNode;
    mutable NodeId lastPostProcessedNode = InvalidNode;

    explicit RecordingAttributeComputer(DynamicComponentTree *tree)
        : tree_(tree) {}

    void resetCounters() const {
        moveProperPartsCount = 0;
        nodeRemovedCount = 0;
        postProcessingCount = 0;
        lastMoveTarget = InvalidNode;
        lastMoveSource = InvalidNode;
        lastRemovedNode = InvalidNode;
        lastPostProcessedNode = InvalidNode;
    }

    void preProcessing(NodeId nodeId, std::span<float> buffer) const override {
        buffer[static_cast<std::size_t>(nodeId)] =
            static_cast<float>(tree_->getNumProperParts(nodeId));
    }

    void mergeProcessing(NodeId parentId, NodeId childId, std::span<float> buffer) const override {
        buffer[static_cast<std::size_t>(parentId)] +=
            buffer[static_cast<std::size_t>(childId)];
    }

    void postProcessing(NodeId nodeId, std::span<float>) const override {
        ++postProcessingCount;
        lastPostProcessedNode = nodeId;
    }

    void onMoveProperParts(NodeId targetId, NodeId sourceId) const override {
        ++moveProperPartsCount;
        lastMoveTarget = targetId;
        lastMoveSource = sourceId;
    }

    void onNodeRemoved(NodeId nodeId) const override {
        ++nodeRemovedCount;
        lastRemovedNode = nodeId;
    }

    void compute(std::span<float> buffer) const {
        if (tree_ == nullptr || tree_->getRoot() == InvalidNode) {
            return;
        }
        computeNode(tree_->getRoot(), buffer);
    }
};

void testUpdateTreeEmitsAdjustmentLog() {
    // With PRINT_LOG enabled, the adjuster should log the main
    // elements of the incremental step: `tau_S`, `F_lambda`, and `nodeCa`.
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
    adjust.setRuntimePostConditionValidationEnabled(true);

    const NodeId rootSubtree = firstNonRootNode(maxTree);
    require(rootSubtree != InvalidNode, "expected a non-root subtree in the max-tree");

    adjust.updateTree(&minTree, rootSubtree);
    const std::string log = adjust.getOutputLog();
    require(!log.empty(), "PRINT_LOG must produce a non-empty log");
    require(log.find("Proper parts (tau_S): [") != std::string::npos, "log must include the tau_S section");
    require(log.find("F_λ = {") != std::string::npos, "log must include the F_lambda section");
    require(log.find("nodeCa: Id:") != std::string::npos, "log must include the nodeCa description");
}

void testDynamicAdjustmentMatchesNaiveBaseline() {
    // In a small case, the incremental min/max flow should produce the same
    // reconstruction as the naive prune-then-rebuild strategy.
    constexpr int kThreshold = 4;
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
    adjust.setRuntimePostConditionValidationEnabled(true);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> maxArea((size_t) maxTree.getNumInternalNodeSlots(), 0.0f);
    std::vector<float> minArea((size_t) minTree.getNumInternalNodeSlots(), 0.0f);
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
    const auto baselineImage = applyNaiveThreshold(input, adj, kThreshold);
    require(dynamicImage->isEqual(baselineImage), "the dynamic CASF result must match the naive baseline");
}

void testDynamicAdjustmentMatchesNaiveOnRecordedSubtreeRegression() {
    const auto inputPath =
        (std::filesystem::path(__FILE__).parent_path().parent_path() / "dat/cameraman.pgm").string();
    const std::vector<int> thresholds = {32, 64, 96, 128, 160, 192};

    auto input = loadPgmImage(inputPath);
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
    adjust.setRuntimePostConditionValidationEnabled(true);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> maxArea((size_t) maxTree.getNumInternalNodeSlots(), 0.0f);
    std::vector<float> minArea((size_t) minTree.getNumInternalNodeSlots(), 0.0f);
    maxAreaComputer.compute(std::span<float>(maxArea));
    minAreaComputer.compute(std::span<float>(minArea));
    adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                std::span<float>(minArea), std::span<float>(maxArea));

    auto baseline = input->clone();
    for (int threshold : thresholds) {
        auto nodesToPrune = getNodesToPrune(maxTree, maxArea, threshold);
        adjust.pruneMaxTreeAndUpdateMinTree(nodesToPrune);

        nodesToPrune = getNodesToPrune(minTree, minArea, threshold);
        adjust.pruneMinTreeAndUpdateMaxTree(nodesToPrune);

        baseline = applyNaiveThreshold(baseline, adj, threshold);
    }

    const auto dynamicImage = minTree.reconstructionImage();
    require(dynamicImage->isEqual(baseline),
            "DualMinMaxTreeIncrementalFilter diverged from the naive baseline in the documented regression case");
}

void testSequentialMintreePrunesMatchDualReconstructionOnLargeFixture() {
    auto input = make_large_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);

    for (int pixelId : {32, 82}) {
        const NodeId rootSubtree = minTree.getSmallestComponent(pixelId);
        require(rootSubtree != InvalidNode, "subtree root must be valid");
        require(minTree.isAlive(rootSubtree), "subtree root must be alive");
        require(!minTree.isRoot(rootSubtree), "subtree root must be non-root");

        adjust.updateTree(&maxTree, rootSubtree);
        minTree.pruneNode(rootSubtree);

        require(maxTree.reconstructionImage()->isEqual(minTree.reconstructionImage()),
                "updated dual tree must reconstruct the same filtered image as the pruned primal tree");
    }
}

void testPublicPruneSubtreeAndUpdateDualKeepsTreesSynchronized() {
    auto input = make_large_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    {
        DynamicComponentTree maxTree(input, true, adj);
        DynamicComponentTree minTree(input, false, adj);
        DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);

        const NodeId rootSubtree = firstStrictNonRootNode(maxTree);
        require(rootSubtree != InvalidNode, "expected a strict max-tree subtree root");

        adjust.pruneSubtreeAndUpdateDual(&maxTree, &minTree, rootSubtree);

        require_tree_consistency(minTree);
        require(minTree.reconstructionImage()->isEqual(maxTree.reconstructionImage()),
                "public pruneSubtreeAndUpdateDual must synchronize min-tree after a max-tree prune");
    }

    {
        DynamicComponentTree maxTree(input, true, adj);
        DynamicComponentTree minTree(input, false, adj);
        DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);

        std::vector<NodeId> nodesToPrune = {firstStrictNonRootNode(minTree)};
        require(nodesToPrune.front() != InvalidNode, "expected a strict min-tree subtree root");

        adjust.pruneSubtreesAndUpdateDual(&minTree, &maxTree,
                                          std::span<const NodeId>(nodesToPrune.data(), nodesToPrune.size()));

        require_tree_consistency(maxTree);
        require(maxTree.reconstructionImage()->isEqual(minTree.reconstructionImage()),
                "public pruneSubtreesAndUpdateDual must synchronize max-tree after a min-tree prune");
    }
}

void testPublicCleanApiRejectsInvalidNodes() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input, true, adj);
    DynamicComponentTree minTree(input, false, adj);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);

    requireThrows([&]() {
        adjust.pruneSubtreeAndUpdateDual(&maxTree, &minTree, maxTree.getRoot());
    }, "public pruneSubtreeAndUpdateDual must reject the primal root");

    requireThrows([&]() {
        adjust.mergeNodeIntoParentAndUpdateDual(&minTree, &maxTree, minTree.getRoot());
    }, "public mergeNodeIntoParentAndUpdateDual must reject the primal root");
}

void requirePublicMergeNodeIntoParentAndUpdateDualSuccess(bool primalIsMinTree) {
    auto input = make_large_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    DynamicComponentTree referenceMaxTree(input, true, adj);
    DynamicComponentTree referenceMinTree(input, false, adj);
    const DynamicComponentTree &referencePrimalTree = primalIsMinTree ? referenceMinTree : referenceMaxTree;

    std::vector<NodeId> candidateNodes;
    for (NodeId nodeId = 0; nodeId < referencePrimalTree.getNumInternalNodeSlots(); ++nodeId) {
        if (referencePrimalTree.isNode(nodeId) &&
            referencePrimalTree.isAlive(nodeId) &&
            !referencePrimalTree.isRoot(nodeId) &&
            referencePrimalTree.getNumProperParts(nodeId) > 0) {
            candidateNodes.push_back(nodeId);
        }
    }
    require(!candidateNodes.empty(), "expected merge-node candidates");

    for (NodeId nodeId : candidateNodes) {
        DynamicComponentTree maxTree(input, true, adj);
        DynamicComponentTree minTree(input, false, adj);
        DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
        DynamicComponentTree *primalTree = primalIsMinTree ? &minTree : &maxTree;
        DynamicComponentTree *dualTree = primalIsMinTree ? &maxTree : &minTree;

        try {
            adjust.mergeNodeIntoParentAndUpdateDual(primalTree, dualTree, nodeId);
        } catch (const std::exception &) {
            continue;
        }

        require_tree_consistency(*primalTree);
        require_tree_consistency(*dualTree);
        require(primalTree->reconstructionImage()->isEqual(dualTree->reconstructionImage()),
                "public mergeNodeIntoParentAndUpdateDual must keep primal and dual reconstructions synchronized");
        return;
    }

    require(false, primalIsMinTree
        ? "expected at least one successful min-tree merge-node event"
        : "expected at least one successful max-tree merge-node event");
}

void testPublicMergeNodeIntoParentAndUpdateDualKeepsTreesSynchronized() {
    requirePublicMergeNodeIntoParentAndUpdateDualSuccess(true);
    requirePublicMergeNodeIntoParentAndUpdateDualSuccess(false);
}

void requirePublicMergeNodeIntoParentAndUpdateDualUpdatesPrimalAttributeState(bool primalIsMinTree) {
    auto input = make_large_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    DynamicComponentTree referenceMaxTree(input, true, adj);
    DynamicComponentTree referenceMinTree(input, false, adj);
    const DynamicComponentTree &referencePrimalTree = primalIsMinTree ? referenceMinTree : referenceMaxTree;

    std::vector<NodeId> candidateNodes;
    for (NodeId nodeId = 0; nodeId < referencePrimalTree.getNumInternalNodeSlots(); ++nodeId) {
        if (referencePrimalTree.isNode(nodeId) &&
            referencePrimalTree.isAlive(nodeId) &&
            !referencePrimalTree.isRoot(nodeId) &&
            referencePrimalTree.getNumProperParts(nodeId) > 0) {
            candidateNodes.push_back(nodeId);
        }
    }
    require(!candidateNodes.empty(), "expected merge-node candidates for attribute hook test");

    for (NodeId nodeId : candidateNodes) {
        DynamicComponentTree maxTree(input, true, adj);
        DynamicComponentTree minTree(input, false, adj);
        DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
        RecordingAttributeComputer maxComputer(&maxTree);
        RecordingAttributeComputer minComputer(&minTree);
        std::vector<float> maxBuffer(static_cast<std::size_t>(maxTree.getNumInternalNodeSlots()), 0.0f);
        std::vector<float> minBuffer(static_cast<std::size_t>(minTree.getNumInternalNodeSlots()), 0.0f);
        maxComputer.compute(std::span<float>(maxBuffer));
        minComputer.compute(std::span<float>(minBuffer));
        maxComputer.resetCounters();
        minComputer.resetCounters();
        adjust.setAttributeComputer(minComputer, maxComputer,
                                    std::span<float>(minBuffer), std::span<float>(maxBuffer));

        DynamicComponentTree *primalTree = primalIsMinTree ? &minTree : &maxTree;
        DynamicComponentTree *dualTree = primalIsMinTree ? &maxTree : &minTree;
        RecordingAttributeComputer &primalComputer = primalIsMinTree ? minComputer : maxComputer;
        std::vector<float> &primalBuffer = primalIsMinTree ? minBuffer : maxBuffer;

        const NodeId parentId = primalTree->getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !primalTree->isAlive(parentId)) {
            continue;
        }

        try {
            adjust.mergeNodeIntoParentAndUpdateDual(primalTree, dualTree, nodeId);
        } catch (const std::exception &) {
            continue;
        }

        require(primalComputer.moveProperPartsCount == 1,
                "public merge-node update must notify the primal attribute about the direct proper-part merge");
        require(primalComputer.lastMoveTarget == parentId,
                "public merge-node update must notify the primal merge target");
        require(primalComputer.lastMoveSource == nodeId,
                "public merge-node update must notify the primal merge source");
        require(primalComputer.nodeRemovedCount == 1,
                "public merge-node update must notify the primal attribute about the removed node");
        require(primalComputer.lastRemovedNode == nodeId,
                "public merge-node update must notify the removed primal node id");
        require(primalComputer.postProcessingCount >= 1 &&
                    primalComputer.lastPostProcessedNode == parentId,
                "public merge-node update must recompute the impacted primal parent");

        DynamicAreaComputer expectedComputer(primalTree);
        const auto expected = expectedComputer.compute();
        require(primalBuffer[static_cast<std::size_t>(parentId)] ==
                    expected[static_cast<std::size_t>(parentId)],
                "public merge-node update must leave the impacted primal parent attribute exact");
        return;
    }

    require(false, primalIsMinTree
        ? "expected at least one successful min-tree merge-node attribute update"
        : "expected at least one successful max-tree merge-node attribute update");
}

void testPublicMergeNodeIntoParentAndUpdateDualUpdatesPrimalAttributeState() {
    requirePublicMergeNodeIntoParentAndUpdateDualUpdatesPrimalAttributeState(true);
    requirePublicMergeNodeIntoParentAndUpdateDualUpdatesPrimalAttributeState(false);
}

void testUpdateTreeKeepsFinalTreeConnectedAndAreaConsistent() {
    auto input = make_large_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    for (int pixelId : {32, 37, 63}) {
        DynamicComponentTree maxTree(input, true, adj);
        DynamicComponentTree minTree(input, false, adj);
        DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
        DynamicAreaComputer maxAreaComputer(&maxTree);
        DynamicAreaComputer minAreaComputer(&minTree);
        std::vector<float> maxArea((size_t) maxTree.getNumInternalNodeSlots(), 0.0f);
        std::vector<float> minArea((size_t) minTree.getNumInternalNodeSlots(), 0.0f);
        maxAreaComputer.compute(std::span<float>(maxArea));
        minAreaComputer.compute(std::span<float>(minArea));
        adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                    std::span<float>(minArea), std::span<float>(maxArea));

        const NodeId rootSubtree = minTree.getSmallestComponent(pixelId);
        require(rootSubtree != InvalidNode, "subtree root must be valid");
        const int oldNumNodes = count_alive_nodes(maxTree);

        std::vector<NodeId> nodesToPrune = {rootSubtree};
        adjust.pruneMinTreeAndUpdateMaxTree(nodesToPrune);

        require(count_alive_nodes(maxTree) < oldNumNodes, "update should contract the affected dual tree");
        require_tree_consistency(maxTree);
        require(maxTree.reconstructionImage()->isEqual(minTree.reconstructionImage()),
                "updated dual tree and pruned primal tree must reconstruct the same image");
    }
}

void testUpdateTreeRemainsStructurallyValidOnAllSharedRoots() {
    auto input = make_large_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);

    for (bool updateMaxTree : {true, false}) {
        DynamicComponentTree referenceMax(input, true, adj);
        DynamicComponentTree referenceMin(input, false, adj);

        std::vector<NodeId> sharedRoots;
        for (int nodeId = 0; nodeId < referenceMax.getNumInternalNodeSlots(); ++nodeId) {
            if (referenceMax.isAlive(nodeId) && referenceMin.isAlive(nodeId) &&
                !referenceMax.isRoot(nodeId) && !referenceMin.isRoot(nodeId)) {
                sharedRoots.push_back(nodeId);
            }
        }

        for (NodeId rootSubtree : sharedRoots) {
            DynamicComponentTree maxTree(input, true, adj);
            DynamicComponentTree minTree(input, false, adj);
            DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adj);
            DynamicAreaComputer maxAreaComputer(&maxTree);
            DynamicAreaComputer minAreaComputer(&minTree);
            std::vector<float> maxArea((size_t) maxTree.getNumInternalNodeSlots(), 0.0f);
            std::vector<float> minArea((size_t) minTree.getNumInternalNodeSlots(), 0.0f);
            maxAreaComputer.compute(std::span<float>(maxArea));
            minAreaComputer.compute(std::span<float>(minArea));
            adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                        std::span<float>(minArea), std::span<float>(maxArea));

            std::vector<NodeId> nodesToPrune = {rootSubtree};
            if (updateMaxTree) {
                adjust.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
                require_tree_consistency(maxTree);
                require(maxTree.reconstructionImage()->isEqual(minTree.reconstructionImage()),
                        "max-tree update must stay synchronized with the pruned mintree");
            } else {
                adjust.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
                require_tree_consistency(minTree);
                require(minTree.reconstructionImage()->isEqual(maxTree.reconstructionImage()),
                        "min-tree update must stay synchronized with the pruned maxtree");
            }
        }
    }
}

} // namespace

int main() {
    try {
        testUpdateTreeEmitsAdjustmentLog();
        testDynamicAdjustmentMatchesNaiveBaseline();
        testDynamicAdjustmentMatchesNaiveOnRecordedSubtreeRegression();
        testSequentialMintreePrunesMatchDualReconstructionOnLargeFixture();
        testPublicPruneSubtreeAndUpdateDualKeepsTreesSynchronized();
        testPublicCleanApiRejectsInvalidNodes();
        testPublicMergeNodeIntoParentAndUpdateDualKeepsTreesSynchronized();
        testPublicMergeNodeIntoParentAndUpdateDualUpdatesPrimalAttributeState();
        testUpdateTreeKeepsFinalTreeConnectedAndAreaConsistent();
        testUpdateTreeRemainsStructurallyValidOnAllSharedRoots();
    } catch (const std::exception &e) {
        std::cerr << "dual_min_max_tree_incremental_filter_unit_tests: " << e.what() << "\n";
        return 1;
    }

    std::cout << "dual_min_max_tree_incremental_filter_unit_tests: OK\n";
    return 0;
}
