#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeCasf.hpp"
#include "../morphoTreeAdjust/include/DynamicComponentTree.hpp"

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
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

std::vector<float> compute_attribute(DynamicComponentTree &tree, Attribute attribute) {
    if (attribute == AREA) {
        DynamicAreaComputer computer(&tree);
        return computer.compute();
    }

    DynamicBoundingBoxComputer computer(&tree, attribute);
    return computer.compute();
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

} // namespace

int main() {
    try {
        test_area_sequence_matches_naive_baseline();
        test_bounding_box_sequences_match_naive_baseline();
        test_casf_is_deterministic_and_empty_sequence_is_noop();
        test_area_stress_matches_naive_sequence_on_structured_images();
        test_naive_and_hybrid_modes_match_baseline();
        test_tree_accessors_expose_internal_trees();
    } catch (const std::exception &e) {
        std::cerr << "dynamic_component_tree_casf_unit_tests: FAIL\n" << e.what() << "\n";
        return 1;
    }

    std::cout << "dynamic_component_tree_casf_unit_tests: OK\n";
    return 0;
}
