#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeCasf.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeCasfNonIncreasing.hpp"
#include "../morphoTreeAdjust/include/DynamicComponentTree.hpp"

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_near(float actual, float expected, float epsilon, const std::string &message) {
    if (std::fabs(actual - expected) > epsilon) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr make_demo_image() {
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

ImageUInt8Ptr make_structured_benchmark_image(int numRows, int numCols) {
    auto image = ImageUInt8::create(numRows, numCols);
    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            const int p = r * numCols + c;
            const int value = (37 * r + 19 * c + 11 * ((r / 4) % 7) + 23 * ((c / 8) % 5) + ((r * c) % 29)) % 256;
            (*image)[p] = static_cast<uint8_t>(value);
        }
    }
    return image;
}

std::vector<int> make_area_thresholds(int numPixels, int numThresholds) {
    std::vector<int> thresholds;
    thresholds.reserve((size_t) numThresholds);
    for (int step = 1; step <= numThresholds; ++step) {
        const double ratio = (double) step / (double) (numThresholds + 1);
        const int threshold = std::max(1, (int) std::llround(ratio * (double) numPixels));
        if (thresholds.empty() || threshold > thresholds.back()) {
            thresholds.push_back(threshold);
        }
    }
    return thresholds;
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

struct NonIncreasingCandidateOrderEntry {
    NodeId nodeId = InvalidNode;
    int altitude = 0;
    PixelId supportMinPixel = std::numeric_limits<PixelId>::max();
    int supportSize = 0;
};

std::vector<NodeId> getNodesToMergeForNonIncreasingAttribute(const DynamicComponentTree &tree,
                                                             const std::vector<float> &attribute,
                                                             float threshold) {
    FastQueue<NodeId> queue;
    std::vector<NodeId> traversalOrder;
    std::vector<NodeId> scanOrder;
    queue.push(tree.getRoot());
    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (nodeId == InvalidNode ||
            !tree.isNode(nodeId) ||
            !tree.isAlive(nodeId)) {
            continue;
        }

        traversalOrder.push_back(nodeId);
        if (nodeId != tree.getRoot() &&
            static_cast<std::size_t>(nodeId) < attribute.size() &&
            attribute[static_cast<std::size_t>(nodeId)] <= threshold) {
            scanOrder.push_back(nodeId);
        }

        for (NodeId childId : tree.getChildren(nodeId)) {
            queue.push(childId);
        }
    }

    const auto nodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    std::vector<int> supportSizeByNode(nodeSlots, 0);
    std::vector<PixelId> supportMinPixelByNode(
        nodeSlots,
        std::numeric_limits<PixelId>::max());

    for (auto it = traversalOrder.rbegin(); it != traversalOrder.rend(); ++it) {
        const NodeId nodeId = *it;
        const auto index = static_cast<std::size_t>(nodeId);
        int supportSize = tree.getNumProperParts(nodeId);
        PixelId supportMinPixel = std::numeric_limits<PixelId>::max();

        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            supportMinPixel = std::min(supportMinPixel, pixelId);
        }

        for (NodeId childId : tree.getChildren(nodeId)) {
            const auto childIndex = static_cast<std::size_t>(childId);
            supportSize += supportSizeByNode[childIndex];
            supportMinPixel = std::min(supportMinPixel, supportMinPixelByNode[childIndex]);
        }

        supportSizeByNode[index] = supportSize;
        supportMinPixelByNode[index] = supportMinPixel;
    }

    std::vector<NonIncreasingCandidateOrderEntry> ordered(scanOrder.size());
    for (std::size_t index = 0; index < scanOrder.size(); ++index) {
        const NodeId nodeId = scanOrder[index];
        auto &entry = ordered[index];
        entry.nodeId = nodeId;
        entry.altitude = tree.getAltitude(nodeId);
        const auto nodeIndex = static_cast<std::size_t>(nodeId);
        entry.supportMinPixel = supportMinPixelByNode[nodeIndex];
        entry.supportSize = supportSizeByNode[nodeIndex];
    }

    std::sort(ordered.begin(), ordered.end(), [](const auto &lhs, const auto &rhs) {
        return std::tie(lhs.altitude, lhs.supportMinPixel, lhs.supportSize, lhs.nodeId) <
               std::tie(rhs.altitude, rhs.supportMinPixel, rhs.supportSize, rhs.nodeId);
    });

    std::vector<NodeId> out;
    out.reserve(ordered.size());
    for (const auto &entry : ordered) {
        out.push_back(entry.nodeId);
    }
    return out;
}

std::vector<float> compute_attribute(DynamicComponentTree &tree, Attribute attribute) {
    if (attribute == AREA) {
        DynamicAreaComputer computer(&tree);
        return computer.compute();
    }

    DynamicBoundingBoxComputer computer(&tree, attribute);
    return computer.compute();
}

float compute_bruteforce_bbox_attribute(const DynamicComponentTree &tree,
                                        NodeId nodeId,
                                        Attribute attribute) {
    int xmin = std::numeric_limits<int>::max();
    int xmax = std::numeric_limits<int>::min();
    int ymin = std::numeric_limits<int>::max();
    int ymax = std::numeric_limits<int>::min();
    int area = 0;

    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            const auto [row, col] = ImageUtils::to2D(pixelId, tree.getNumColsOfImage());
            xmin = std::min(xmin, col);
            xmax = std::max(xmax, col);
            ymin = std::min(ymin, row);
            ymax = std::max(ymax, row);
            ++area;
        }
    }

    if (area == 0) {
        return 0.0f;
    }

    const float width = static_cast<float>(xmax - xmin + 1);
    const float height = static_cast<float>(ymax - ymin + 1);
    const float areaValue = static_cast<float>(area);
    if (attribute == BOX_WIDTH) {
        return width;
    }
    if (attribute == BOX_HEIGHT) {
        return height;
    }
    if (attribute == DIAGONAL_LENGTH) {
        return std::sqrt(width * width + height * height);
    }
    if (attribute == AREA_PER_BOX_PERIMETER) {
        return areaValue / (2.0f * width + 2.0f * height);
    }
    if (attribute == RECTANGULARITY) {
        return areaValue / (width * height);
    }
    if (attribute == BOX_ASPECT_RATIO) {
        return std::max(width, height) / std::min(width, height);
    }
    if (attribute == BOX_ELONGATION) {
        return 1.0f - (std::min(width, height) / std::max(width, height));
    }
    if (attribute == AREA_PER_BOX_DIAGONAL_SQUARED) {
        return areaValue / (width * width + height * height);
    }
    throw std::runtime_error("unexpected attribute in brute-force bbox computation");
}

ImageUInt8Ptr apply_naive_threshold(ImageUInt8Ptr input,
                                    AdjacencyRelationPtr adj,
                                    int threshold,
                                    Attribute attribute) {
    DynamicComponentTree maxTree(input, true, adj);
    const auto maxAttribute = compute_attribute(maxTree, attribute);
    for (NodeId nodeId : getNodesToPrune(maxTree, maxAttribute, threshold)) {
        if (nodeId != maxTree.getRoot()) {
            maxTree.pruneNode(nodeId);
        }
    }

    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adj);
    const auto minAttribute = compute_attribute(minTree, attribute);
    for (NodeId nodeId : getNodesToPrune(minTree, minAttribute, threshold)) {
        if (nodeId != minTree.getRoot()) {
            minTree.pruneNode(nodeId);
        }
    }
    return minTree.reconstructionImage();
}

ImageUInt8Ptr apply_naive_non_increasing_threshold(ImageUInt8Ptr input,
                                                   AdjacencyRelationPtr adj,
                                                   float threshold,
                                                   Attribute attribute) {
    DynamicComponentTree maxTree(input, true, adj);
    const auto maxAttribute = compute_attribute(maxTree, attribute);
    for (NodeId nodeId : getNodesToMergeForNonIncreasingAttribute(maxTree, maxAttribute, threshold)) {
        if (nodeId != maxTree.getRoot() &&
            maxTree.isNode(nodeId) &&
            maxTree.isAlive(nodeId) &&
            maxTree.getNumProperParts(nodeId) > 0) {
            maxTree.mergeNodeIntoParent(nodeId);
        }
    }

    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adj);
    const auto minAttribute = compute_attribute(minTree, attribute);
    for (NodeId nodeId : getNodesToMergeForNonIncreasingAttribute(minTree, minAttribute, threshold)) {
        if (nodeId != minTree.getRoot() &&
            minTree.isNode(nodeId) &&
            minTree.isAlive(nodeId) &&
            minTree.getNumProperParts(nodeId) > 0) {
            minTree.mergeNodeIntoParent(nodeId);
        }
    }

    return minTree.reconstructionImage();
}

ImageUInt8Ptr run_naive_sequence(ImageUInt8Ptr input,
                                 AdjacencyRelationPtr adj,
                                 const std::vector<int> &thresholds,
                                 Attribute attribute) {
    auto current = input->clone();
    for (int threshold : thresholds) {
        current = apply_naive_threshold(current, adj, threshold, attribute);
    }
    return current;
}

ImageUInt8Ptr run_naive_non_increasing_sequence(ImageUInt8Ptr input,
                                                AdjacencyRelationPtr adj,
                                                const std::vector<float> &thresholds,
                                                Attribute attribute) {
    auto current = input->clone();
    for (float threshold : thresholds) {
        current = apply_naive_non_increasing_threshold(current, adj, threshold, attribute);
    }
    return current;
}

void test_area_sequence_matches_naive_baseline() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    const std::vector<int> thresholds = {1, 14, 15, 100};

    ComponentTreeCasf<AltitudeType> runner(input, 1.0, AREA);
    const auto filtered = runner.filter(thresholds, ComponentTreeCasf<AltitudeType>::Mode::Updating);

    const auto baseline = run_naive_sequence(input, adj, thresholds, AREA);
    require(filtered->isEqual(baseline),
            "dynamic subtree CASF by area must match the naive rebuild baseline");
}

void test_bounding_box_sequences_match_naive_baseline() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    const std::vector<int> thresholds = {2, 4, 6};

    for (Attribute attribute : {BOX_WIDTH, BOX_HEIGHT, DIAGONAL_LENGTH}) {
        ComponentTreeCasf<AltitudeType> runner(input, 1.0, attribute);
        const auto filtered = runner.filter(thresholds);

        const auto baseline = run_naive_sequence(input, adj, thresholds, attribute);
        require(filtered->isEqual(baseline),
                "dynamic subtree CASF by bounding box must match the naive rebuild baseline");
    }
}

void test_bounding_box_derived_attributes_match_bruteforce() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree tree(input, true, adj);

    for (Attribute attribute : {
             BOX_WIDTH,
             BOX_HEIGHT,
             DIAGONAL_LENGTH,
             AREA_PER_BOX_PERIMETER,
             RECTANGULARITY,
             BOX_ASPECT_RATIO,
             BOX_ELONGATION,
             AREA_PER_BOX_DIAGONAL_SQUARED,
         }) {
        DynamicBoundingBoxComputer computer(&tree, attribute);
        const auto values = computer.compute();
        for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
            if (!tree.isAlive(nodeId)) {
                continue;
            }
            const float expected = compute_bruteforce_bbox_attribute(tree, nodeId, attribute);
            require_near(
                values[static_cast<std::size_t>(nodeId)],
                expected,
                1e-5f,
                "bbox-derived attribute must match brute-force subtree geometry");
        }
    }
}

void test_bounding_box_dirty_target_move_matches_fresh_recompute() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree tree(input, true, adj);
    DynamicBoundingBoxComputer cachedComputer(&tree, BOX_HEIGHT);

    const NodeId donor = tree.getSmallestComponent(32);
    const NodeId singlePixelChild = tree.getSmallestComponent(44);
    const NodeId receiver = tree.getSmallestComponent(96);

    tree.moveProperPart(receiver, donor, 32);
    cachedComputer.onMoveProperPart(receiver, donor, 32);

    tree.moveProperPart(donor, singlePixelChild, 44);
    cachedComputer.onMoveProperPart(donor, singlePixelChild, 44);

    tree.moveProperPart(receiver, donor, 44);
    cachedComputer.onMoveProperPart(receiver, donor, 44);

    const auto cached = cachedComputer.compute();
    DynamicBoundingBoxComputer freshComputer(&tree, BOX_HEIGHT);
    const auto fresh = freshComputer.compute();

    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId)) {
            require(cached[(size_t) nodeId] == fresh[(size_t) nodeId],
                    "bounding box cache must not duplicate border counts when a dirty target receives a pixel");
        }
    }
}

void test_casf_is_deterministic_and_empty_sequence_is_noop() {
    auto input = make_demo_image();
    const std::vector<int> thresholds = {1, 14, 15, 100};

    ComponentTreeCasf<AltitudeType> runnerA(input, 1.0, AREA);
    ComponentTreeCasf<AltitudeType> runnerB(input, 1.0, AREA);
    const auto imageA = runnerA.filter(thresholds);
    const auto imageB = runnerB.filter(thresholds);

    require(imageA->isEqual(imageB), "independent CASF runners must be deterministic");

    ComponentTreeCasf<AltitudeType> noopRunner(input, 1.0, AREA);
    require(noopRunner.filter({})->isEqual(input), "empty threshold sequence must preserve the original image");
}

void test_area_stress_matches_naive_sequence_on_structured_images() {
    for (int size : {32, 64, 128, 256}) {
        auto input = make_structured_benchmark_image(size, size);
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
        const auto thresholds = make_area_thresholds(size * size, 10);

        ComponentTreeCasf<AltitudeType> runner(input, 1.0, AREA);
        const auto filtered = runner.filter(thresholds);

        const auto baseline = run_naive_sequence(input, adj, thresholds, AREA);
        require(filtered->isEqual(baseline),
                "dynamic subtree CASF stress sequence must match the naive rebuild baseline");
    }
}

void test_naive_and_hybrid_modes_match_baseline() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    const std::vector<int> thresholds = {1, 14, 15, 100};

    ComponentTreeCasf<AltitudeType> naiveRunner(input, 1.0, AREA);
    ComponentTreeCasf<AltitudeType> hybridRunner(input, 1.0, AREA);

    const auto naiveImage = naiveRunner.filter(thresholds, ComponentTreeCasf<AltitudeType>::Mode::Naive);
    const auto hybridImage = hybridRunner.filter(thresholds, "hybrid");
    const auto baseline = run_naive_sequence(input, adj, thresholds, AREA);

    require(naiveImage->isEqual(baseline), "ComponentTreeCasf naive mode must match the naive rebuild baseline");
    require(hybridImage->isEqual(baseline), "ComponentTreeCasf hybrid mode must match the naive rebuild baseline");
}

void test_tree_accessors_expose_internal_trees() {
    auto input = make_demo_image();
    const std::vector<int> thresholds = {1, 14};

    ComponentTreeCasf<AltitudeType> runner(input, 1.0, AREA);
    const auto filtered = runner.filter(thresholds);

    require(runner.getMaxTree().reconstructionImage()->isEqual(filtered),
            "ComponentTreeCasf max-tree accessor must expose the internal tree state");
    require(runner.getMinTree().reconstructionImage()->isEqual(filtered),
            "ComponentTreeCasf min-tree accessor must expose the internal tree state");
}

void test_non_increasing_bbox_derived_attributes_match_naive_baseline() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    struct Case {
        Attribute attribute;
        std::vector<float> thresholds;
    };
    const std::vector<Case> cases = {
        {AREA_PER_BOX_PERIMETER, {0.25f, 0.4f, 0.75f}},
        {RECTANGULARITY, {0.35f, 0.6f, 0.85f}},
        {BOX_ASPECT_RATIO, {1.0f, 1.25f, 1.75f}},
        {BOX_ELONGATION, {0.05f, 0.25f, 0.5f}},
        {AREA_PER_BOX_DIAGONAL_SQUARED, {0.15f, 0.3f, 0.45f}},
    };

    for (const auto &testCase : cases) {
        ComponentTreeCasfNonIncreasing<AltitudeType> runner(
            input,
            adj,
            testCase.attribute);
        const auto filtered = runner.filter(testCase.thresholds);

        const auto baseline = run_naive_non_increasing_sequence(
            input,
            adj,
            testCase.thresholds,
            testCase.attribute);
        require(filtered->isEqual(baseline),
                "non-increasing CASF by bbox-derived attribute must match the naive rebuild baseline");
    }
}

} // namespace

int main() {
    try {
        test_area_sequence_matches_naive_baseline();
        test_bounding_box_sequences_match_naive_baseline();
        test_bounding_box_derived_attributes_match_bruteforce();
        test_bounding_box_dirty_target_move_matches_fresh_recompute();
        test_casf_is_deterministic_and_empty_sequence_is_noop();
        test_area_stress_matches_naive_sequence_on_structured_images();
        test_naive_and_hybrid_modes_match_baseline();
        test_tree_accessors_expose_internal_trees();
        test_non_increasing_bbox_derived_attributes_match_naive_baseline();
    } catch (const std::exception &e) {
        std::cerr << "dynamic_component_tree_casf_unit_tests: FAIL\n" << e.what() << "\n";
        return 1;
    }

    std::cout << "dynamic_component_tree_casf_unit_tests: OK\n";
    return 0;
}
