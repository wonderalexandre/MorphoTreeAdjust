#pragma once

#include <span>
#include <string>
#include <vector>

#include "../Region.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"
#include "ExplicitRegionProblem1Preflight.hpp"

namespace debug::problem1 {

using ::debug::region::AnchoredRegionLevelChange;
using ::debug::region::RegionLevelChange;
using ::debug::region::RegionLevelChangeRebuildUpdater;
using ::debug::region::RegionMaterializationPlan;
using ::debug::region::RegionMergeWorkspace;
using ::debug::region::RegionUpdateDirection;
using ::debug::region::analyzeRegionMaterialization;
using ::debug::region::buildExplicitAnchoredRegionChange;
using ::debug::region::buildExplicitRegionChange;
using ::debug::region::buildRegionMergeWorkspace;
using ::debug::region::buildSubtreeRegionChange;
using ::debug::region::isConnectedRegion;
using ::debug::region::materializePartialOwners;
using ::debug::region::materializePartialOwnersInsideOnly;

struct ExplicitRegionProblem1ContractDiagnostics {
    int affectedOwners = 0;
    int partialOwners = 0;
    int fullyCoveredOwners = 0;
    int affectedLevels = 0;
    int affectedIntervalLevels = 0;
    int missingLevelUnionCandidates = 0;
    int missingIntermediateLevelUnionCandidates = 0;
    bool hasTargetLevelUnionCandidate = false;
    std::string currentUpdaterLimitation;
};

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1Contract {
    Problem1PreflightResult<PixelType> preflight;
    ExplicitRegionProblem1ContractDiagnostics diagnostics;

    // Legacy diagnostic mirrors kept for current debug callers.
    bool acceptedByPublishedProblem = false;
    bool compatibleWithCurrentUpdater = false;
    bool isMaxtree = false;
    RegionUpdateDirection direction = RegionUpdateDirection::Unchanged;
    NodeId nodeCa = InvalidNode;
    PixelType a{};
    PixelType minC{};
    PixelType maxC{};
    PixelType b{};
    int regionSize = 0;
    int affectedOwners = 0;
    int partialOwners = 0;
    int fullyCoveredOwners = 0;
    int affectedLevels = 0;
    int affectedIntervalLevels = 0;
    int missingLevelUnionCandidates = 0;
    int missingIntermediateLevelUnionCandidates = 0;
    bool hasTargetLevelUnionCandidate = false;
    std::string rejectionReason;
    std::string currentUpdaterLimitation;

    bool acceptsPublishedProblem() const {
        return preflight.acceptedByPublishedProblem();
    }

    bool readyForLocalUpdater() const {
        return preflight.readyForCurrentUpdater();
    }

    bool requiresLocalMaterialization() const {
        return preflight.requiresLocalMaterialization();
    }

    bool requiresTargetNodeCreation() const {
        return preflight.requiresTargetNodeCreation();
    }
};

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1Contract<PixelType> rejectExplicitRegionProblem1Contract(
    Problem1PreflightStatus status,
    const std::string &reason) {
    ExplicitRegionProblem1Contract<PixelType> out;
    out.preflight.status = status;
    out.acceptedByPublishedProblem =
        problem1PreflightStatusAcceptsPublishedProblem(status);
    out.compatibleWithCurrentUpdater =
        problem1PreflightStatusIsUpdaterReady(status);
    out.rejectionReason = reason;
    return out;
}

template<typename PixelType = AltitudeType>
void syncExplicitRegionProblem1ContractMirrors(
    ExplicitRegionProblem1Contract<PixelType> &contract) {
    contract.acceptedByPublishedProblem =
        contract.preflight.acceptedByPublishedProblem();
    contract.compatibleWithCurrentUpdater =
        contract.preflight.readyForCurrentUpdater();
    contract.isMaxtree = contract.preflight.isMaxtree;
    contract.direction = contract.preflight.direction;
    contract.nodeCa = contract.preflight.nodeCa;
    contract.a = contract.preflight.a;
    contract.minC = contract.preflight.minC;
    contract.maxC = contract.preflight.maxC;
    contract.b = contract.preflight.b;
    contract.regionSize = contract.preflight.regionSize;

    contract.affectedOwners = contract.diagnostics.affectedOwners;
    contract.partialOwners = contract.diagnostics.partialOwners;
    contract.fullyCoveredOwners = contract.diagnostics.fullyCoveredOwners;
    contract.affectedLevels = contract.diagnostics.affectedLevels;
    contract.affectedIntervalLevels = contract.diagnostics.affectedIntervalLevels;
    contract.missingLevelUnionCandidates =
        contract.diagnostics.missingLevelUnionCandidates;
    contract.missingIntermediateLevelUnionCandidates =
        contract.diagnostics.missingIntermediateLevelUnionCandidates;
    contract.hasTargetLevelUnionCandidate =
        contract.diagnostics.hasTargetLevelUnionCandidate;
    contract.currentUpdaterLimitation =
        contract.diagnostics.currentUpdaterLimitation;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1Contract<PixelType> checkExplicitRegionProblem1Contract(
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false) {
    if (pixelsC.empty()) {
        return rejectExplicitRegionProblem1Contract<PixelType>(
            Problem1PreflightStatus::EmptyRegion,
            "Problem 1 requires a non-empty region C");
    }

    const int domainSize = dualTree.getNumTotalProperParts();
    std::vector<unsigned char> inRegion(static_cast<size_t>(domainSize), 0);

    bool firstPixel = true;
    NodeId nodeCa = InvalidNode;
    PixelType minC{};
    PixelType maxC{};

    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            return rejectExplicitRegionProblem1Contract<PixelType>(
                Problem1PreflightStatus::OutOfBoundsPixel,
                "Problem 1 region C contains an out-of-bounds pixel");
        }

        const size_t pixelIndex = static_cast<size_t>(pixelId);
        if (inRegion[pixelIndex] != 0) {
            return rejectExplicitRegionProblem1Contract<PixelType>(
                Problem1PreflightStatus::DuplicatePixel,
                "Problem 1 region C contains duplicated pixels");
        }
        inRegion[pixelIndex] = 1;

        const NodeId ownerNode = dualTree.getSmallestComponent(pixelId);
        if (ownerNode == InvalidNode || !dualTree.isAlive(ownerNode)) {
            return rejectExplicitRegionProblem1Contract<PixelType>(
                Problem1PreflightStatus::MissingLiveOwner,
                "Problem 1 region C contains a pixel without a live owner");
        }

        const PixelType level = static_cast<PixelType>(dualTree.getAltitude(ownerNode));
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

    if (!assumeConnected && !isConnectedRegion(adjacency, pixelsC)) {
        return rejectExplicitRegionProblem1Contract<PixelType>(
            Problem1PreflightStatus::DisconnectedRegion,
            "Problem 1 requires connected region C");
    }

    ExplicitRegionProblem1Contract<PixelType> out;
    out.preflight.status = Problem1PreflightStatus::Accepted;
    out.preflight.isMaxtree = dualTree.isMaxtree();
    out.preflight.nodeCa = nodeCa;
    out.preflight.minC = minC;
    out.preflight.maxC = maxC;
    out.preflight.b = b;
    out.preflight.regionSize = static_cast<int>(pixelsC.size());

    if (dualTree.isMaxtree()) {
        out.preflight.a = minC;
        out.preflight.direction = directionFromLevels(maxC, b);
        if (b <= maxC) {
            return rejectExplicitRegionProblem1Contract<PixelType>(
                Problem1PreflightStatus::InvalidTargetLevel,
                "Problem 1 on a max-tree requires b > max I-(C)");
        }
    } else {
        out.preflight.a = maxC;
        out.preflight.direction = directionFromLevels(minC, b);
        if (b >= minC) {
            return rejectExplicitRegionProblem1Contract<PixelType>(
                Problem1PreflightStatus::InvalidTargetLevel,
                "The min-tree dual of Problem 1 requires b < min I-(C)");
        }
    }

    const RegionMaterializationPlan materialization =
        analyzeRegionMaterialization(dualTree, adjacency, pixelsC);
    out.diagnostics.affectedOwners =
        static_cast<int>(materialization.affectedOwners.size());
    out.diagnostics.partialOwners = materialization.numPartialOwners();
    out.diagnostics.fullyCoveredOwners = materialization.numFullyCoveredOwners();

    RegionLevelChange<PixelType> change;
    change.pixelsC.assign(pixelsC.begin(), pixelsC.end());
    change.a = out.preflight.a;
    change.b = b;
    change.direction = dualTree.isMaxtree()
        ? RegionUpdateDirection::Raise
        : RegionUpdateDirection::Lower;
    change.assumeConnected = true;
    change.assumeMaterializedAsProperParts = !materialization.requiresLocalMaterialization();

    const RegionMergeWorkspace<PixelType> workspace =
        buildRegionMergeWorkspace(dualTree, adjacency, change);
    out.diagnostics.affectedLevels =
        static_cast<int>(workspace.affectedThresholds.size());

    const int firstAffectedLevel = dualTree.isMaxtree()
        ? static_cast<int>(out.preflight.a) + 1
        : static_cast<int>(out.preflight.a) - 1;
    const int lastAffectedLevel = static_cast<int>(b);
    const int step = dualTree.isMaxtree() ? 1 : -1;
    for (int level = firstAffectedLevel;
         dualTree.isMaxtree() ? level <= lastAffectedLevel : level >= lastAffectedLevel;
         level += step) {
        ++out.diagnostics.affectedIntervalLevels;
        const auto *bucket = workspace.mergeNodesVMinus(static_cast<PixelType>(level));
        if (bucket == nullptr || bucket->empty()) {
            ++out.diagnostics.missingLevelUnionCandidates;
            if (level != lastAffectedLevel) {
                ++out.diagnostics.missingIntermediateLevelUnionCandidates;
            }
        }
    }

    const auto *targetBucket = workspace.mergeNodesVMinus(b);
    out.diagnostics.hasTargetLevelUnionCandidate =
        targetBucket != nullptr && !targetBucket->empty();

    if (materialization.requiresLocalMaterialization()) {
        out.preflight.status = Problem1PreflightStatus::RequiresLocalMaterialization;
        out.diagnostics.currentUpdaterLimitation =
            "C cuts direct proper parts; local materialization is required first";
    } else if (!out.diagnostics.hasTargetLevelUnionCandidate) {
        out.preflight.status = Problem1PreflightStatus::RequiresTargetNodeCreation;
        out.diagnostics.currentUpdaterLimitation =
            "current updater has no level-b union candidate and must create a local node";
    } else {
        out.preflight.status = Problem1PreflightStatus::Accepted;
    }

    syncExplicitRegionProblem1ContractMirrors(out);

    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1Contract<PixelType> checkExplicitRegionProblem1Contract(
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    const RegionLevelChange<PixelType> &change) {
    return checkExplicitRegionProblem1Contract<PixelType>(
        dualTree,
        adjacency,
        std::span<const PixelId>(change.pixelsC),
        change.b,
        change.assumeConnected);
}

} // namespace debug::problem1
