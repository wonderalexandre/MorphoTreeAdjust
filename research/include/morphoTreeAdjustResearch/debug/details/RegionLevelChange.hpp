#pragma once

#include <algorithm>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"
#include "RegionUpdateDirection.hpp"

namespace debug::region {

template<typename PixelType = AltitudeType>
struct RegionLevelChange {
    std::vector<PixelId> pixelsC;
    PixelType a{};
    PixelType b{};
    RegionUpdateDirection direction = RegionUpdateDirection::Unchanged;
    bool assumeConnected = false;
    bool assumeMaterializedAsProperParts = false;
};

template<typename PixelType = AltitudeType>
struct AnchoredRegionLevelChange {
    RegionLevelChange<PixelType> change;
    NodeId nodeCa = InvalidNode;
};

inline bool isConnectedRegion(AdjacencyRelation &graph, std::span<const PixelId> pixelsC) {
    if (pixelsC.empty()) {
        return false;
    }

    std::unordered_set<PixelId> region;
    region.reserve(pixelsC.size());
    for (PixelId pixelId : pixelsC) {
        region.insert(pixelId);
    }
    if (region.size() != pixelsC.size()) {
        return false;
    }

    std::vector<PixelId> stack;
    std::unordered_set<PixelId> visited;
    stack.reserve(pixelsC.size());
    visited.reserve(pixelsC.size());
    stack.push_back(pixelsC.front());
    visited.insert(pixelsC.front());

    while (!stack.empty()) {
        const PixelId current = stack.back();
        stack.pop_back();
        for (PixelId neighbor : graph.getNeighborPixels(current)) {
            if (region.find(neighbor) == region.end()) {
                continue;
            }
            if (visited.insert(neighbor).second) {
                stack.push_back(neighbor);
            }
        }
    }

    return visited.size() == region.size();
}

template<typename PixelType = AltitudeType>
std::pair<NodeId, PixelType> findNodeCa(const DynamicComponentTree &dualTree,
                                        std::span<const PixelId> pixelsC) {
    NodeId nodeCa = InvalidNode;
    PixelType altitudeCa{};
    const bool isMaxtree = dualTree.isMaxtree();

    for (PixelId pixelId : pixelsC) {
        const NodeId ownerNodeId = dualTree.getSmallestComponent(pixelId);
        if (ownerNodeId == InvalidNode) {
            continue;
        }

        const PixelType altitude = static_cast<PixelType>(dualTree.getAltitude(ownerNodeId));
        if (nodeCa == InvalidNode ||
            (isMaxtree && altitude < altitudeCa) ||
            (!isMaxtree && altitude > altitudeCa)) {
            nodeCa = ownerNodeId;
            altitudeCa = altitude;
        }
    }

    if (nodeCa == InvalidNode) {
        throw std::runtime_error("RegionLevelChange cannot locate nodeCa for an empty or ownerless region");
    }

    return {nodeCa, altitudeCa};
}

template<typename PixelType = AltitudeType>
RegionLevelChange<PixelType> buildExplicitRegionChange(std::span<const PixelId> pixelsC,
                                                       PixelType a,
                                                       PixelType b,
                                                       bool assumeConnected = false,
                                                       bool assumeMaterializedAsProperParts = false) {
    RegionLevelChange<PixelType> change;
    change.pixelsC.assign(pixelsC.begin(), pixelsC.end());
    change.a = a;
    change.b = b;
    change.direction = directionFromLevels(a, b);
    change.assumeConnected = assumeConnected;
    change.assumeMaterializedAsProperParts = assumeMaterializedAsProperParts;
    return change;
}

template<typename PixelType = AltitudeType>
AnchoredRegionLevelChange<PixelType> buildExplicitAnchoredRegionChange(const DynamicComponentTree &dualTree,
                                                                       std::span<const PixelId> pixelsC,
                                                                       PixelType b,
                                                                       bool assumeConnected = false,
                                                                       bool assumeMaterializedAsProperParts = false) {
    auto [nodeCa, a] = findNodeCa<PixelType>(dualTree, pixelsC);
    AnchoredRegionLevelChange<PixelType> out;
    out.change = buildExplicitRegionChange<PixelType>(pixelsC, a, b, assumeConnected, assumeMaterializedAsProperParts);
    out.nodeCa = nodeCa;
    return out;
}

template<typename PixelType = AltitudeType>
AnchoredRegionLevelChange<PixelType> buildSubtreeRegionChange(const DynamicComponentTree &primalTree,
                                                             const DynamicComponentTree &dualTree,
                                                             NodeId subtreeRoot) {
    if (subtreeRoot == InvalidNode ||
        !primalTree.isNode(subtreeRoot) ||
        !primalTree.isAlive(subtreeRoot) ||
        primalTree.isRoot(subtreeRoot)) {
        throw std::runtime_error("buildSubtreeRegionChange expects a live non-root subtree root");
    }

    const NodeId subtreeParentId = primalTree.getNodeParent(subtreeRoot);
    if (subtreeParentId == InvalidNode || subtreeParentId == subtreeRoot) {
        throw std::runtime_error("buildSubtreeRegionChange cannot locate the subtree parent");
    }

    std::vector<PixelId> pixelsC;
    for (NodeId subtreeNodeId : primalTree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : primalTree.getProperParts(subtreeNodeId)) {
            pixelsC.push_back(pixelId);
        }
    }

    const PixelType b = static_cast<PixelType>(primalTree.getAltitude(subtreeParentId));
    return buildExplicitAnchoredRegionChange<PixelType>(dualTree, std::span<const PixelId>(pixelsC), b, true, true);
}

} // namespace debug::region
