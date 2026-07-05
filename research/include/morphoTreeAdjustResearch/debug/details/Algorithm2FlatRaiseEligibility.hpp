#pragma once

#include <span>
#include <string>
#include <utility>
#include <vector>

#include "RegionLevelChange.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem2 {

using ::debug::region::buildExplicitRegionChange;
using ::debug::region::isConnectedRegion;
using ::debug::region::RegionLevelChange;

struct Algorithm2FlatRaiseEligibility {
    bool accepted = false;
    std::string reason;
    NodeId sourceNodeQa = InvalidNode;
    AltitudeType a = 0;
    AltitudeType b = 0;
    int pixelsP = 0;
    int directProperPartsOfQa = 0;
    int childrenOfQa = 0;
    int higherBoundaryNeighbors = 0;
    bool releasesQa = false;
};

inline Algorithm2FlatRaiseEligibility rejectAlgorithm2FlatRaise(
    std::string reason,
    std::span<const PixelId> pixelsP,
    AltitudeType a = 0,
    AltitudeType b = 0) {
    Algorithm2FlatRaiseEligibility out;
    out.accepted = false;
    out.reason = std::move(reason);
    out.a = a;
    out.b = b;
    out.pixelsP = static_cast<int>(pixelsP.size());
    return out;
}

inline bool isDirectProperPartOfNode(const DynamicComponentTree &tree,
                                     NodeId nodeId,
                                     PixelId pixelId) {
    for (PixelId properPart : tree.getProperParts(nodeId)) {
        if (properPart == pixelId) {
            return true;
        }
    }
    return false;
}

inline int countHigherBoundaryNeighbors(const DynamicComponentTree &maxTree,
                                        AdjacencyRelation &adjacency,
                                        std::span<const PixelId> pixelsP,
                                        AltitudeType a) {
    std::vector<unsigned char> inRegion(
        static_cast<std::size_t>(maxTree.getNumTotalProperParts()),
        0);
    for (PixelId pixelId : pixelsP) {
        inRegion[static_cast<std::size_t>(pixelId)] = 1;
    }

    int count = 0;
    for (PixelId pixelId : pixelsP) {
        for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
            if (inRegion[static_cast<std::size_t>(neighbor)] != 0) {
                continue;
            }
            const NodeId ownerNode = maxTree.getSmallestComponent(neighbor);
            if (ownerNode != InvalidNode &&
                maxTree.isAlive(ownerNode) &&
                maxTree.getAltitude(ownerNode) > static_cast<int>(a)) {
                ++count;
            }
        }
    }
    return count;
}

template<typename PixelType = AltitudeType>
Algorithm2FlatRaiseEligibility checkAlgorithm2FlatRaiseEligibility(
    const DynamicComponentTree &maxTree,
    AdjacencyRelation &adjacency,
    const RegionLevelChange<PixelType> &change) {
    if (!maxTree.isMaxtree()) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise requires a max-tree",
            std::span<const PixelId>(change.pixelsC),
            static_cast<AltitudeType>(change.a),
            static_cast<AltitudeType>(change.b));
    }
    if (change.pixelsC.empty()) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise requires a non-empty flat region P",
            std::span<const PixelId>(change.pixelsC),
            static_cast<AltitudeType>(change.a),
            static_cast<AltitudeType>(change.b));
    }
    if (change.b <= change.a) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise currently supports only raise updates b > a",
            std::span<const PixelId>(change.pixelsC),
            static_cast<AltitudeType>(change.a),
            static_cast<AltitudeType>(change.b));
    }

    const int domainSize = maxTree.getNumTotalProperParts();
    std::vector<unsigned char> seen(static_cast<std::size_t>(domainSize), 0);
    for (PixelId pixelId : change.pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            return rejectAlgorithm2FlatRaise(
                "Algorithm 2 flat raise received an out-of-bounds pixel",
                std::span<const PixelId>(change.pixelsC),
                static_cast<AltitudeType>(change.a),
                static_cast<AltitudeType>(change.b));
        }
        auto index = static_cast<std::size_t>(pixelId);
        if (seen[index] != 0) {
            return rejectAlgorithm2FlatRaise(
                "Algorithm 2 flat raise received duplicated pixels in P",
                std::span<const PixelId>(change.pixelsC),
                static_cast<AltitudeType>(change.a),
                static_cast<AltitudeType>(change.b));
        }
        seen[index] = 1;
    }

    if (!change.assumeConnected &&
        !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise requires connected region P",
            std::span<const PixelId>(change.pixelsC),
            static_cast<AltitudeType>(change.a),
            static_cast<AltitudeType>(change.b));
    }

    NodeId sourceNode = InvalidNode;
    for (PixelId pixelId : change.pixelsC) {
        const NodeId ownerNode = maxTree.getSmallestComponent(pixelId);
        if (ownerNode == InvalidNode || !maxTree.isAlive(ownerNode)) {
            return rejectAlgorithm2FlatRaise(
                "Algorithm 2 flat raise cannot locate a live owner for P",
                std::span<const PixelId>(change.pixelsC),
                static_cast<AltitudeType>(change.a),
                static_cast<AltitudeType>(change.b));
        }
        if (maxTree.getAltitude(ownerNode) != static_cast<int>(change.a)) {
            return rejectAlgorithm2FlatRaise(
                "Algorithm 2 flat raise requires every pixel in P to have level a",
                std::span<const PixelId>(change.pixelsC),
                static_cast<AltitudeType>(change.a),
                static_cast<AltitudeType>(change.b));
        }
        if (sourceNode == InvalidNode) {
            sourceNode = ownerNode;
        } else if (sourceNode != ownerNode) {
            return rejectAlgorithm2FlatRaise(
                "Algorithm 2 flat raise currently requires P to be owned by a single source node Q_a",
                std::span<const PixelId>(change.pixelsC),
                static_cast<AltitudeType>(change.a),
                static_cast<AltitudeType>(change.b));
        }
    }

    if (sourceNode == InvalidNode) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise cannot locate Q_a",
            std::span<const PixelId>(change.pixelsC),
            static_cast<AltitudeType>(change.a),
            static_cast<AltitudeType>(change.b));
    }
    if (maxTree.getAltitude(sourceNode) != static_cast<int>(change.a)) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise requires Q_a altitude to be a",
            std::span<const PixelId>(change.pixelsC),
            static_cast<AltitudeType>(change.a),
            static_cast<AltitudeType>(change.b));
    }
    for (PixelId pixelId : change.pixelsC) {
        if (!isDirectProperPartOfNode(maxTree, sourceNode, pixelId)) {
            return rejectAlgorithm2FlatRaise(
                "Algorithm 2 flat raise requires every pixel in P to be a direct proper part of Q_a",
                std::span<const PixelId>(change.pixelsC),
                static_cast<AltitudeType>(change.a),
                static_cast<AltitudeType>(change.b));
        }
    }

    Algorithm2FlatRaiseEligibility out;
    out.accepted = true;
    out.reason = "accepted";
    out.sourceNodeQa = sourceNode;
    out.a = static_cast<AltitudeType>(change.a);
    out.b = static_cast<AltitudeType>(change.b);
    out.pixelsP = static_cast<int>(change.pixelsC.size());
    out.directProperPartsOfQa = maxTree.getNumProperParts(sourceNode);
    out.childrenOfQa = maxTree.getNumChildren(sourceNode);
    out.higherBoundaryNeighbors = countHigherBoundaryNeighbors(
        maxTree,
        adjacency,
        std::span<const PixelId>(change.pixelsC),
        out.a);
    out.releasesQa = out.pixelsP == out.directProperPartsOfQa;
    return out;
}

inline Algorithm2FlatRaiseEligibility checkAlgorithm2FlatRaiseEligibility(
    const DynamicComponentTree &maxTree,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsP,
    AltitudeType b,
    bool assumeConnected = false,
    bool assumeMaterializedAsProperParts = true) {
    if (pixelsP.empty()) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise requires a non-empty flat region P",
            pixelsP,
            0,
            b);
    }
    if (pixelsP.front() < 0 || pixelsP.front() >= maxTree.getNumTotalProperParts()) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise received an out-of-bounds pixel",
            pixelsP,
            0,
            b);
    }

    const NodeId firstOwner = maxTree.getSmallestComponent(pixelsP.front());
    if (firstOwner == InvalidNode || !maxTree.isAlive(firstOwner)) {
        return rejectAlgorithm2FlatRaise(
            "Algorithm 2 flat raise cannot locate a live owner for P",
            pixelsP,
            0,
            b);
    }

    auto change = buildExplicitRegionChange<AltitudeType>(
        pixelsP,
        static_cast<AltitudeType>(maxTree.getAltitude(firstOwner)),
        b,
        assumeConnected,
        assumeMaterializedAsProperParts);
    return checkAlgorithm2FlatRaiseEligibility(maxTree, adjacency, change);
}

template<typename PixelType = AltitudeType>
Algorithm2FlatRaiseEligibility requireAlgorithm2FlatRaiseEligibility(
    const DynamicComponentTree &maxTree,
    AdjacencyRelation &adjacency,
    const RegionLevelChange<PixelType> &change) {
    auto eligibility = checkAlgorithm2FlatRaiseEligibility(maxTree, adjacency, change);
    if (!eligibility.accepted) {
        throw std::runtime_error(eligibility.reason);
    }
    return eligibility;
}

} // namespace debug::problem2
