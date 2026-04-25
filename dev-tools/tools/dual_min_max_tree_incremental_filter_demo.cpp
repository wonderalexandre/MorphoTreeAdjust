#include <iostream>
#include <memory>
#include <vector>

#include "../../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../../morphoTreeAdjust/include/Common.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilter.hpp"

static std::vector<int> getNodesThreshold(const DynamicComponentTree &tree,
                                          const std::vector<float> &area,
                                          int threshold) {
    std::vector<int> out;
    FastQueue<int> queue;
    queue.push(tree.getRoot());

    while (!queue.empty()) {
        const int nodeId = queue.pop();
        if (area[(size_t) nodeId] > threshold) {
            for (int childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }
    return out;
}

static ImageUInt8Ptr runBaselineAdjustMin(ImageUInt8Ptr input,
                                          AdjacencyRelationPtr adj,
                                          int threshold) {
    DynamicComponentTree maxTree(input, true, adj);
    const auto nodesToPrune = DynamicComponentTree::getNodesThreshold(&maxTree, threshold);
    for (int nodeId : nodesToPrune) {
        if (nodeId != maxTree.getRoot()) {
            maxTree.pruneNode(nodeId);
        }
    }
    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adj);
    return minTree.reconstructionImage();
}

int main() {
    constexpr int kThreshold = 4;
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

    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree dynamicMaxTree(input, true, adj);
    DynamicComponentTree dynamicMinTree(input, false, adj);

    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&dynamicMinTree, &dynamicMaxTree, *adj);
    DynamicAreaComputer maxAreaComputer(&dynamicMaxTree);
    DynamicAreaComputer minAreaComputer(&dynamicMinTree);
    std::vector<float> areaMax((size_t) dynamicMaxTree.getNumInternalNodeSlots(), 0.0f);
    std::vector<float> areaMin((size_t) dynamicMinTree.getNumInternalNodeSlots(), 0.0f);
    maxAreaComputer.compute(std::span<float>(areaMax));
    minAreaComputer.compute(std::span<float>(areaMin));
    adjust.setAttributeComputer(minAreaComputer, maxAreaComputer, std::span<float>(areaMin), std::span<float>(areaMax));
    auto nodesToPrune = getNodesThreshold(dynamicMaxTree, areaMax, kThreshold);
    if (nodesToPrune.empty()) {
        std::cerr << "dual_min_max_tree_incremental_filter_demo: no nodes to prune\n";
        return 1;
    }

    adjust.pruneMaxTreeAndUpdateMinTree(nodesToPrune);

    if (!dynamicMinTree.isAlive(dynamicMinTree.getRoot()) || !dynamicMaxTree.isAlive(dynamicMaxTree.getRoot())) {
        std::cerr << "dual_min_max_tree_incremental_filter_demo: root became invalid after adjustMinTree\n";
        return 1;
    }
    if (dynamicMinTree.getNumInternalNodeSlots() <= 0 || dynamicMaxTree.getNumInternalNodeSlots() <= 0) {
        std::cerr << "dual_min_max_tree_incremental_filter_demo: invalid internal node slot count after adjustMinTree\n";
        return 1;
    }

    const auto dynamicImage = dynamicMinTree.reconstructionImage();
    const auto baselineImage = runBaselineAdjustMin(input, adj, kThreshold);
    if (!dynamicImage->isEqual(baselineImage)) {
        std::cerr << "dual_min_max_tree_incremental_filter_demo: dynamic result diverges from rebuild baseline\n";
        return 1;
    }

    std::cout << "DualMinMaxTreeIncrementalFilter demo: OK\n";
    return 0;
}
