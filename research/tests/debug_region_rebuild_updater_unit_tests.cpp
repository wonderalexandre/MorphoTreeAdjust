#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChangeRebuildUpdater.hpp"

namespace region_debug = debug::region;

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr make_demo_image() {
    auto image = ImageUInt8::create(4, 4);
    const std::vector<uint8_t> raw = {
        4, 4, 2, 1,
        4, 3, 2, 1,
        5, 5, 2, 0,
        5, 6, 6, 0,
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

NodeId firstNonRootNode(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isRoot(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

void requireRebuildPairMatchesImage(const DynamicComponentTree &minTree,
                                    const DynamicComponentTree &maxTree,
                                    const ImageUInt8Ptr &expected) {
    require(minTree.reconstructionImage()->isEqual(expected),
            "rebuilt min-tree reconstruction must match expected image");
    require(maxTree.reconstructionImage()->isEqual(expected),
            "rebuilt max-tree reconstruction must match expected image");
    require(minTree.reconstructionImage()->isEqual(maxTree.reconstructionImage()),
            "rebuilt min-tree and max-tree reconstructions must match");
}

void requireThrowsDisconnectedRegion(DynamicComponentTree &minTree,
                                     DynamicComponentTree &maxTree,
                                     const AdjacencyRelationPtr &adj,
                                     const region_debug::RegionLevelChange<AltitudeType> &change) {
    try {
        (void) region_debug::RegionLevelChangeRebuildUpdater::updatePair(minTree, maxTree, adj, change);
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error("disconnected region update should throw");
}

void testExplicitNonSubtreeRaiseByRebuild() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsC = {0, 1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        static_cast<AltitudeType>(7),
        false,
        true);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionLevelChangeRebuildUpdater::updatePair(minTree, maxTree, adj, change);

    require(updated->isEqual(expected), "updater must return the globally changed image");
    requireRebuildPairMatchesImage(minTree, maxTree, expected);
}

void testSubtreeDescriptorByRebuild() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const NodeId subtreeRoot = firstNonRootNode(maxTree);
    require(subtreeRoot != InvalidNode, "test image must have a non-root subtree");

    auto anchored = region_debug::buildSubtreeRegionChange<AltitudeType>(maxTree, minTree, subtreeRoot);
    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, anchored.change);
    auto updated = region_debug::RegionLevelChangeRebuildUpdater::updatePair(minTree, maxTree, adj, anchored.change);

    require(updated->isEqual(expected), "subtree descriptor updater must return the globally changed image");
    requireRebuildPairMatchesImage(minTree, maxTree, expected);
}

void testDisconnectedRegionValidation() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> disconnectedC = {0, 15};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(disconnectedC),
        static_cast<AltitudeType>(4),
        static_cast<AltitudeType>(7),
        false,
        true);

    requireThrowsDisconnectedRegion(minTree, maxTree, adj, change);
    requireRebuildPairMatchesImage(minTree, maxTree, input);
}

} // namespace

int main() {
    testExplicitNonSubtreeRaiseByRebuild();
    testSubtreeDescriptorByRebuild();
    testDisconnectedRegionValidation();
    return 0;
}
