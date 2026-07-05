#pragma once

#include <span>
#include <unordered_set>
#include <vector>

#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"
#include "RegionUpdateDirection.hpp"

namespace debug::problem1 {

using debug::region::RegionUpdateDirection;
using debug::region::directionFromLevels;

enum class Problem1PreflightStatus {
    Accepted,
    EmptyRegion,
    OutOfBoundsPixel,
    DuplicatePixel,
    MissingLiveOwner,
    DisconnectedRegion,
    InvalidTargetLevel,
    RequiresLocalMaterialization,
    RequiresTargetNodeCreation,
    UnsupportedLocalTopology,
};

inline bool problem1PreflightStatusAcceptsPublishedProblem(Problem1PreflightStatus status) {
    switch (status) {
        case Problem1PreflightStatus::Accepted:
        case Problem1PreflightStatus::RequiresLocalMaterialization:
        case Problem1PreflightStatus::RequiresTargetNodeCreation:
        case Problem1PreflightStatus::UnsupportedLocalTopology:
            return true;
        case Problem1PreflightStatus::EmptyRegion:
        case Problem1PreflightStatus::OutOfBoundsPixel:
        case Problem1PreflightStatus::DuplicatePixel:
        case Problem1PreflightStatus::MissingLiveOwner:
        case Problem1PreflightStatus::DisconnectedRegion:
        case Problem1PreflightStatus::InvalidTargetLevel:
            return false;
    }
    return false;
}

inline bool problem1PreflightStatusIsUpdaterReady(
    Problem1PreflightStatus status,
    bool allowTargetNodeCreation = false) {
    if (status == Problem1PreflightStatus::Accepted) {
        return true;
    }
    return status == Problem1PreflightStatus::RequiresTargetNodeCreation &&
        allowTargetNodeCreation;
}

inline bool problem1PreflightStatusRequiresLocalMaterialization(
    Problem1PreflightStatus status) {
    return status == Problem1PreflightStatus::RequiresLocalMaterialization;
}

inline bool problem1PreflightStatusRequiresTargetNodeCreation(
    Problem1PreflightStatus status) {
    return status == Problem1PreflightStatus::RequiresTargetNodeCreation;
}

inline const char *problem1PreflightStatusName(Problem1PreflightStatus status) {
    switch (status) {
        case Problem1PreflightStatus::Accepted:
            return "accepted";
        case Problem1PreflightStatus::EmptyRegion:
            return "empty_region";
        case Problem1PreflightStatus::OutOfBoundsPixel:
            return "out_of_bounds_pixel";
        case Problem1PreflightStatus::DuplicatePixel:
            return "duplicate_pixel";
        case Problem1PreflightStatus::MissingLiveOwner:
            return "missing_live_owner";
        case Problem1PreflightStatus::DisconnectedRegion:
            return "disconnected_region";
        case Problem1PreflightStatus::InvalidTargetLevel:
            return "invalid_target_level";
        case Problem1PreflightStatus::RequiresLocalMaterialization:
            return "requires_local_materialization";
        case Problem1PreflightStatus::RequiresTargetNodeCreation:
            return "requires_target_node_creation";
        case Problem1PreflightStatus::UnsupportedLocalTopology:
            return "unsupported_local_topology";
    }
    return "unknown";
}

template<typename PixelType = AltitudeType>
struct Problem1PreflightResult {
    Problem1PreflightStatus status = Problem1PreflightStatus::EmptyRegion;
    bool isMaxtree = false;
    RegionUpdateDirection direction = RegionUpdateDirection::Unchanged;
    NodeId nodeCa = InvalidNode;
    PixelType a{};
    PixelType minC{};
    PixelType maxC{};
    PixelType b{};
    int regionSize = 0;

    bool acceptedByPublishedProblem() const {
        return problem1PreflightStatusAcceptsPublishedProblem(status);
    }

    bool readyForCurrentUpdater(bool allowTargetNodeCreation = false) const {
        return problem1PreflightStatusIsUpdaterReady(
            status,
            allowTargetNodeCreation);
    }

    bool requiresLocalMaterialization() const {
        return problem1PreflightStatusRequiresLocalMaterialization(status);
    }

    bool requiresTargetNodeCreation() const {
        return problem1PreflightStatusRequiresTargetNodeCreation(status);
    }
};

inline bool isConnectedExplicitRegionProblem1Region(AdjacencyRelation &graph,
                                                    std::span<const PixelId> pixelsC) {
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

namespace detail {

template<typename PixelType = AltitudeType>
Problem1PreflightResult<PixelType> rejectPreflight(
    Problem1PreflightStatus status) {
    Problem1PreflightResult<PixelType> out;
    out.status = status;
    return out;
}

} // namespace detail

template<typename PixelType = AltitudeType>
bool hasExplicitRegionProblem1TargetLevelUnionCandidate(
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsC,
    const std::vector<unsigned char> &inRegion,
    PixelType a,
    PixelType b) {
    const int domainSize = dualTree.getNumTotalProperParts();
    const int nodeSlots = dualTree.getNumInternalNodeSlots();
    std::vector<unsigned char> adjacentSeedMarks(static_cast<size_t>(nodeSlots), 0);
    std::vector<unsigned char> climbedNodeMarks(static_cast<size_t>(nodeSlots), 0);

    auto markNode = [nodeSlots](std::vector<unsigned char> &marks, NodeId nodeId) -> bool {
        if (nodeId < 0 || nodeId >= nodeSlots) {
            return false;
        }
        const size_t index = static_cast<size_t>(nodeId);
        if (marks[index] != 0) {
            return false;
        }
        marks[index] = 1;
        return true;
    };

    const bool isMaxtree = dualTree.isMaxtree();
    for (PixelId pixelInC : pixelsC) {
        for (PixelId neighbor : adjacency.getNeighborPixels(pixelInC)) {
            if (neighbor < 0 || neighbor >= domainSize) {
                continue;
            }
            if (inRegion[static_cast<size_t>(neighbor)] != 0) {
                continue;
            }

            const NodeId seed = dualTree.getSmallestComponent(neighbor);
            if (seed == InvalidNode || !dualTree.isAlive(seed)) {
                continue;
            }

            const PixelType seedLevel =
                static_cast<PixelType>(dualTree.getAltitude(seed));
            const bool validSeed =
                (isMaxtree && seedLevel >= a) ||
                (!isMaxtree && seedLevel <= a);
            if (!validSeed || !markNode(adjacentSeedMarks, seed)) {
                continue;
            }

            NodeId current = seed;
            while (current != InvalidNode && dualTree.isAlive(current)) {
                if (!markNode(climbedNodeMarks, current)) {
                    break;
                }

                const PixelType currentLevel =
                    static_cast<PixelType>(dualTree.getAltitude(current));
                const bool insideNodeCaInterval =
                    (isMaxtree && currentLevel >= a) ||
                    (!isMaxtree && currentLevel <= a);
                if (!insideNodeCaInterval) {
                    break;
                }

                const bool insideMergeInterval =
                    (isMaxtree && currentLevel <= b) ||
                    (!isMaxtree && currentLevel >= b);
                if (insideMergeInterval && currentLevel == b) {
                    return true;
                }

                const NodeId parentId = dualTree.getNodeParent(current);
                if (parentId == current) {
                    break;
                }
                current = parentId;
            }
        }
    }

    return false;
}

template<typename PixelType = AltitudeType>
Problem1PreflightResult<PixelType> checkExplicitRegionProblem1Preflight(
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false) {
    if (pixelsC.empty()) {
        return detail::rejectPreflight<PixelType>(
            Problem1PreflightStatus::EmptyRegion);
    }

    const int domainSize = dualTree.getNumTotalProperParts();
    const int nodeSlots = dualTree.getNumInternalNodeSlots();
    std::vector<unsigned char> inRegion(static_cast<size_t>(domainSize), 0);
    std::vector<int> ownerInsideCounts(static_cast<size_t>(nodeSlots), 0);
    std::vector<NodeId> affectedOwners;

    bool firstPixel = true;
    NodeId nodeCa = InvalidNode;
    PixelType minC{};
    PixelType maxC{};

    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            return detail::rejectPreflight<PixelType>(
                Problem1PreflightStatus::OutOfBoundsPixel);
        }

        const size_t pixelIndex = static_cast<size_t>(pixelId);
        if (inRegion[pixelIndex] != 0) {
            return detail::rejectPreflight<PixelType>(
                Problem1PreflightStatus::DuplicatePixel);
        }
        inRegion[pixelIndex] = 1;

        const NodeId ownerNode = dualTree.getSmallestComponent(pixelId);
        if (ownerNode == InvalidNode ||
            ownerNode < 0 ||
            ownerNode >= nodeSlots ||
            !dualTree.isAlive(ownerNode)) {
            return detail::rejectPreflight<PixelType>(
                Problem1PreflightStatus::MissingLiveOwner);
        }

        if (ownerInsideCounts[static_cast<size_t>(ownerNode)] == 0) {
            affectedOwners.push_back(ownerNode);
        }
        ++ownerInsideCounts[static_cast<size_t>(ownerNode)];

        const PixelType level =
            static_cast<PixelType>(dualTree.getAltitude(ownerNode));
        if (firstPixel) {
            minC = level;
            maxC = level;
            nodeCa = ownerNode;
            firstPixel = false;
        } else {
            if (level < minC) {
                minC = level;
            }
            if (level > maxC) {
                maxC = level;
            }
            if ((dualTree.isMaxtree() && level == minC) ||
                (!dualTree.isMaxtree() && level == maxC)) {
                nodeCa = ownerNode;
            }
        }
    }

    if (!assumeConnected &&
        !isConnectedExplicitRegionProblem1Region(adjacency, pixelsC)) {
        return detail::rejectPreflight<PixelType>(
            Problem1PreflightStatus::DisconnectedRegion);
    }

    Problem1PreflightResult<PixelType> out;
    out.status = Problem1PreflightStatus::Accepted;
    out.isMaxtree = dualTree.isMaxtree();
    out.nodeCa = nodeCa;
    out.minC = minC;
    out.maxC = maxC;
    out.b = b;
    out.regionSize = static_cast<int>(pixelsC.size());

    if (dualTree.isMaxtree()) {
        out.a = minC;
        out.direction = directionFromLevels(maxC, b);
        if (b <= maxC) {
            out.status = Problem1PreflightStatus::InvalidTargetLevel;
            return out;
        }
    } else {
        out.a = maxC;
        out.direction = directionFromLevels(minC, b);
        if (b >= minC) {
            out.status = Problem1PreflightStatus::InvalidTargetLevel;
            return out;
        }
    }

    for (NodeId ownerNode : affectedOwners) {
        const int insideCount = ownerInsideCounts[static_cast<size_t>(ownerNode)];
        if (insideCount > 0 && insideCount < dualTree.getNumProperParts(ownerNode)) {
            out.status = Problem1PreflightStatus::RequiresLocalMaterialization;
            return out;
        }
    }

    if (!hasExplicitRegionProblem1TargetLevelUnionCandidate<PixelType>(
            dualTree,
            adjacency,
            pixelsC,
            inRegion,
            out.a,
            b)) {
        out.status = Problem1PreflightStatus::RequiresTargetNodeCreation;
        return out;
    }

    return out;
}

} // namespace debug::problem1
