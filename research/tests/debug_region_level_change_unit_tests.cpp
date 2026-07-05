#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"

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

std::vector<PixelId> collectSubtreePixels(const DynamicComponentTree &tree, NodeId subtreeRoot) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(pixelId);
        }
    }
    std::sort(pixels.begin(), pixels.end());
    return pixels;
}

NodeId firstNonRootNode(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isRoot(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

template<typename PixelType>
region_debug::RegionUpdateDirection expectedDirection(PixelType a, PixelType b) {
    if (b > a) {
        return region_debug::RegionUpdateDirection::Raise;
    }
    if (b < a) {
        return region_debug::RegionUpdateDirection::Lower;
    }
    return region_debug::RegionUpdateDirection::Unchanged;
}

void testExplicitRegionChangeCopiesPixelsAndClassifiesDirection() {
    std::vector<PixelId> pixels = {1, 2, 5};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixels),
        static_cast<AltitudeType>(2),
        static_cast<AltitudeType>(5),
        true,
        false);

    pixels[0] = 99;
    require(change.pixelsC == std::vector<PixelId>({1, 2, 5}), "explicit region must copy pixels");
    require(change.a == static_cast<AltitudeType>(2), "explicit region must keep level a");
    require(change.b == static_cast<AltitudeType>(5), "explicit region must keep level b");
    require(change.direction == region_debug::RegionUpdateDirection::Raise, "explicit region must classify raise updates");
    require(change.assumeConnected, "explicit region must keep connectedness assumption");
    require(!change.assumeMaterializedAsProperParts, "explicit region must keep materialization assumption");

    auto lower = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(change.pixelsC),
        static_cast<AltitudeType>(5),
        static_cast<AltitudeType>(2));
    require(lower.direction == region_debug::RegionUpdateDirection::Lower, "explicit region must classify lower updates");
}

void testConnectedRegionValidation() {
    AdjacencyRelation graph(2, 2, 1.0);
    const std::vector<PixelId> connected = {0, 1, 3};
    const std::vector<PixelId> disconnected = {0, 3};
    const std::vector<PixelId> duplicated = {0, 1, 1};

    require(region_debug::isConnectedRegion(graph, std::span<const PixelId>(connected)), "region should be connected");
    require(!region_debug::isConnectedRegion(graph, std::span<const PixelId>(disconnected)), "region should be disconnected");
    require(!region_debug::isConnectedRegion(graph, std::span<const PixelId>(duplicated)), "duplicated region should be rejected");
}

void testSubtreeBuilderMatchesCurrentImplicitCollection() {
    auto image = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adj);
    DynamicComponentTree minTree(image->clone(), false, adj);

    const NodeId subtreeRoot = firstNonRootNode(maxTree);
    require(subtreeRoot != InvalidNode, "test image must have a non-root node");

    auto anchored = region_debug::buildSubtreeRegionChange<AltitudeType>(maxTree, minTree, subtreeRoot);
    auto expectedPixels = collectSubtreePixels(maxTree, subtreeRoot);
    auto actualPixels = anchored.change.pixelsC;
    std::sort(actualPixels.begin(), actualPixels.end());

    const NodeId subtreeParent = maxTree.getNodeParent(subtreeRoot);
    require(actualPixels == expectedPixels, "subtree builder must collect the current subtree proper parts");
    require(anchored.nodeCa != InvalidNode, "subtree builder must locate nodeCa");
    require(anchored.change.a == static_cast<AltitudeType>(minTree.getAltitude(anchored.nodeCa)),
            "subtree builder level a must match nodeCa altitude in the dual tree");
    require(anchored.change.b == static_cast<AltitudeType>(maxTree.getAltitude(subtreeParent)),
            "subtree builder level b must match the subtree parent altitude");
    require(anchored.change.direction == expectedDirection(anchored.change.a, anchored.change.b),
            "subtree builder must classify the update direction");
    require(anchored.change.assumeConnected, "subtree builder may assume connected C");
    require(anchored.change.assumeMaterializedAsProperParts,
            "subtree builder may assume C is materialized as current proper parts");

    auto [nodeCa, altitudeCa] = region_debug::findNodeCa<AltitudeType>(
        minTree,
        std::span<const PixelId>(anchored.change.pixelsC));
    require(nodeCa == anchored.nodeCa, "findNodeCa must match the subtree builder anchor");
    require(altitudeCa == anchored.change.a, "findNodeCa must match the subtree builder level a");
}

} // namespace

int main() {
    testExplicitRegionChangeCopiesPixelsAndClassifiesDirection();
    testConnectedRegionValidation();
    testSubtreeBuilderMatchesCurrentImplicitCollection();
    return 0;
}
