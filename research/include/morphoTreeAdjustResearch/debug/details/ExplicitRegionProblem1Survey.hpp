#pragma once

#include <algorithm>
#include <set>
#include <span>
#include <string>
#include <vector>

#include "ExplicitRegionProblem1CheckedUpdate.hpp"
#include "ExplicitRegionProblem1LocalCorrectnessValidator.hpp"
#include "ExplicitRegionProblem1MultiAnchor.hpp"
#include "ExplicitRegionProblem1Oracle.hpp"
#include "ExplicitRegionProblem1Updater.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"

namespace debug::problem1 {

using ::debug::structural::compareDynamicComponentTreeStructure;

struct ExplicitRegionProblem1SurveySummary {
    int connectedRegions = 0;
    int totalEvents = 0;
    int rejectedByContract = 0;
    int acceptedByPublishedProblem = 0;
    int compatibleWithCurrentUpdater = 0;
    int knownUpdaterLimitations = 0;
    int requiresMaterialization = 0;
    int requiresTargetNodeCreation = 0;
    int otherKnownLimitations = 0;
    int materializationAttempted = 0;
    int materializationApplied = 0;
    int localCorrectnessValidationPerformed = 0;
    int localCorrectnessValidationPassed = 0;
    int rejectedByLocalCorrectnessValidation = 0;
    int localCorrectnessRejectedH2 = 0;
    int localCorrectnessRejectedH5 = 0;
    int localCorrectnessRejectedH6 = 0;
    int localCorrectnessRejectedH7 = 0;
    int localCorrectnessRejectedH8 = 0;
    int materializedLocalUpdateMatchedRebuild = 0;
    int targetNodeCreationLocalUpdateMatchedRebuild = 0;
    int materializedTargetNodeCreationLocalUpdateMatchedRebuild = 0;
    int localUpdateMatchedRebuild = 0;
    int localUpdateMismatchedRebuild = 0;
    int localUpdateThrew = 0;
    int localMaterializationThrew = 0;
    int safeUpdateAttempted = 0;
    int safeUpdateApplied = 0;
    int safeUpdateMatchedRebuild = 0;
    int safeUpdateMismatchedRebuild = 0;
    int directUpdateAttempted = 0;
    int directUpdateApplied = 0;
    int directUpdateMatchedRebuild = 0;
    int directUpdateMismatchedRebuild = 0;
    int directUpdateThrew = 0;
    int rebuildFallbackAttempted = 0;
    int rebuildFallbackApplied = 0;
    int rebuildFallbackMatchedRebuild = 0;
    int rebuildFallbackMismatchedRebuild = 0;
    int rebuildFallbackThrew = 0;
    int structuralValidationPerformed = 0;
    int structuralMatchedRebuild = 0;
    int structuralMismatchedRebuild = 0;
    int sameLevelEdgesCollapsed = 0;
    long long affectedIntervalLevels = 0;
    long long missingLevelUnionCandidates = 0;
    long long missingIntermediateLevelUnionCandidates = 0;

    bool hasLocalUpdateFailure() const {
        return localUpdateMismatchedRebuild != 0 ||
               localUpdateThrew != 0 ||
               localMaterializationThrew != 0 ||
               safeUpdateMismatchedRebuild != 0 ||
               directUpdateMismatchedRebuild != 0 ||
               directUpdateThrew != 0 ||
               rebuildFallbackMismatchedRebuild != 0 ||
               rebuildFallbackThrew != 0;
    }
};

struct ExplicitRegionProblem1MultiAnchorSurveySummary {
    int connectedRegions = 0;
    int totalEvents = 0;
    int rejectedByContract = 0;
    int acceptedByPublishedProblem = 0;
    int analyzedEvents = 0;
    int materializationAttempted = 0;
    int materializationApplied = 0;
    int multipleOwnerEvents = 0;
    int multiOwnerSingleWorkspaceEvents = 0;
    int multiAnchorCandidateEvents = 0;
    int multipleExtremalOwnerEvents = 0;
    int disconnectedOwnerPartitionEvents = 0;
    int multipleSparseWorkspaceAnchorEvents = 0;
    int multipleGeneralWorkspaceAnchorEvents = 0;
    int closureRootEvents = 0;
    int closureRootMatchesNodeCaEvents = 0;
    int closureRootDiffersFromNodeCaEvents = 0;
    int multiOwnerClosureRootDiffersFromNodeCaEvents = 0;
    int nodeCaReachesClosureRootEvents = 0;
    int nodeCaToClosureRootSameLevelPathEvents = 0;
    int closureRootDiffersButSameLevelPathEvents = 0;
    int ownerIntervalBucketGapEvents = 0;
    int ownerIntervalBucketUnexplainedGapEvents = 0;
    int localCorrectnessCandidateEvents = 0;
    int h6UniqueGeneralWorkspaceRootEvents = 0;
    int h7SameLevelClosureEvents = 0;
    int h8BucketCompletenessEvents = 0;
    int targetLevelPresentOrCreatedEvents = 0;
    long long nodeCaToClosureRootPathLengthSum = 0;
    long long owners = 0;
    long long extremalOwners = 0;
    long long ownerPartitionComponents = 0;
    long long disconnectedOwnerPartitions = 0;
    long long ownerIncidenceEdges = 0;
    long long sparseWorkspaceNodes = 0;
    long long sparseWorkspaceAnchors = 0;
    long long ownerIntervalWorkspaceNodes = 0;
    long long ownerIntervalWorkspaceAnchors = 0;
    long long generalWorkspaceNodes = 0;
    long long generalWorkspaceAnchors = 0;
    long long ownerIntervalBucketMissingNodes = 0;
    long long ownerIntervalBucketClosurePathMissingNodes = 0;
    long long ownerIntervalBucketCOnlyProperPartMissingNodes = 0;
    long long ownerIntervalBucketUnexplainedMissingNodes = 0;
    long long frontierNodes = 0;
    int maxOwners = 0;
    int maxExtremalOwners = 0;
    int maxOwnerPartitionComponents = 0;
    int maxOwnerIncidenceEdges = 0;
    int maxSparseWorkspaceAnchors = 0;
    int maxOwnerIntervalWorkspaceAnchors = 0;
    int maxGeneralWorkspaceAnchors = 0;
    int maxNodeCaToClosureRootPathLength = 0;
    int maxOwnerIntervalBucketMissingNodes = 0;
    int maxOwnerIntervalBucketUnexplainedMissingNodes = 0;

    bool sawMultiAnchorCandidate() const {
        return multiAnchorCandidateEvents > 0;
    }
};

namespace detail {

template<typename PixelType = AltitudeType>
bool localCorrectnessFailedHypothesis(
    const ExplicitRegionProblem1LocalCorrectnessValidation<PixelType> &validation,
    const std::string &hypothesis) {
    return std::find(
               validation.failedHypotheses.begin(),
               validation.failedHypotheses.end(),
               hypothesis) != validation.failedHypotheses.end();
}

inline void enumeratePixelCombinations(int domainSize,
                                       int targetSize,
                                       int start,
                                       std::vector<PixelId> &current,
                                       std::vector<std::vector<PixelId>> &out) {
    if (static_cast<int>(current.size()) == targetSize) {
        out.push_back(current);
        return;
    }

    const int remaining = targetSize - static_cast<int>(current.size());
    for (int pixelId = start; pixelId <= domainSize - remaining; ++pixelId) {
        current.push_back(static_cast<PixelId>(pixelId));
        enumeratePixelCombinations(domainSize, targetSize, pixelId + 1, current, out);
        current.pop_back();
    }
}

inline std::vector<std::vector<PixelId>> enumerateConnectedRegions(AdjacencyRelation &adjacency,
                                                                   int domainSize,
                                                                   int maxRegionSize) {
    std::vector<std::vector<PixelId>> connectedRegions;
    if (domainSize <= 0 || maxRegionSize <= 0) {
        return connectedRegions;
    }

    const int boundedMaxRegionSize = std::min(domainSize, maxRegionSize);
    for (int regionSize = 1; regionSize <= boundedMaxRegionSize; ++regionSize) {
        std::vector<std::vector<PixelId>> candidates;
        std::vector<PixelId> current;
        current.reserve(static_cast<size_t>(regionSize));
        enumeratePixelCombinations(domainSize, regionSize, 0, current, candidates);

        for (auto &candidate : candidates) {
            if (isConnectedRegion(adjacency, std::span<const PixelId>(candidate))) {
                connectedRegions.push_back(std::move(candidate));
            }
        }
    }

    return connectedRegions;
}

inline void appendUniqueConnectedRegion(AdjacencyRelation &adjacency,
                                        std::vector<PixelId> region,
                                        int maxRegionSize,
                                        std::set<std::vector<PixelId>> &seen,
                                        std::vector<std::vector<PixelId>> &out) {
    std::sort(region.begin(), region.end());
    region.erase(std::unique(region.begin(), region.end()), region.end());
    if (region.empty() || static_cast<int>(region.size()) > maxRegionSize) {
        return;
    }
    if (!isConnectedRegion(adjacency, std::span<const PixelId>(region))) {
        return;
    }
    if (seen.insert(region).second) {
        out.push_back(std::move(region));
    }
}

inline std::vector<std::vector<PixelId>> enumerateSeededConnectedRegions(
    AdjacencyRelation &adjacency,
    int domainSize,
    int maxRegionSize,
    int maxRegions) {
    std::vector<std::vector<PixelId>> connectedRegions;
    if (domainSize <= 0 || maxRegionSize <= 0 || maxRegions <= 0) {
        return connectedRegions;
    }

    std::set<std::vector<PixelId>> seen;
    const int boundedMaxRegionSize = std::min(domainSize, maxRegionSize);
    for (PixelId seed = 0;
         seed < domainSize && static_cast<int>(connectedRegions.size()) < maxRegions;
         ++seed) {
        appendUniqueConnectedRegion(
            adjacency,
            {seed},
            boundedMaxRegionSize,
            seen,
            connectedRegions);

        std::vector<PixelId> grownRegion = {seed};
        for (int step = 1; step < boundedMaxRegionSize; ++step) {
            std::vector<PixelId> candidates;
            for (PixelId pixelId : grownRegion) {
                for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
                    if (neighbor < 0 || neighbor >= domainSize) {
                        continue;
                    }
                    if (std::find(grownRegion.begin(),
                                  grownRegion.end(),
                                  neighbor) != grownRegion.end()) {
                        continue;
                    }
                    candidates.push_back(neighbor);
                }
            }

            std::sort(candidates.begin(), candidates.end());
            candidates.erase(std::unique(candidates.begin(), candidates.end()),
                             candidates.end());
            if (candidates.empty()) {
                break;
            }

            const PixelId chosen =
                candidates[static_cast<size_t>(seed + step * 3) % candidates.size()];
            grownRegion.push_back(chosen);
            appendUniqueConnectedRegion(
                adjacency,
                grownRegion,
                boundedMaxRegionSize,
                seen,
                connectedRegions);
        }

        for (PixelId neighbor : adjacency.getNeighborPixels(seed)) {
            if (neighbor < 0 || neighbor >= domainSize) {
                continue;
            }
            appendUniqueConnectedRegion(
                adjacency,
                {seed, neighbor},
                boundedMaxRegionSize,
                seen,
                connectedRegions);

            for (PixelId secondNeighbor : adjacency.getNeighborPixels(neighbor)) {
                if (secondNeighbor < 0 ||
                    secondNeighbor >= domainSize ||
                    secondNeighbor == seed) {
                    continue;
                }
                appendUniqueConnectedRegion(
                    adjacency,
                    {seed, neighbor, secondNeighbor},
                    boundedMaxRegionSize,
                    seen,
                    connectedRegions);
            }
        }
    }

    if (static_cast<int>(connectedRegions.size()) > maxRegions) {
        connectedRegions.resize(static_cast<size_t>(maxRegions));
    }
    return connectedRegions;
}

template<typename PixelType = AltitudeType>
void accumulateExplicitRegionProblem1SurveyResult(
    ExplicitRegionProblem1SurveySummary &summary,
    const ExplicitRegionProblem1OracleResult<PixelType> &result,
    bool countDirectUpdateMethod = false) {
    if (result.materializationAttempted) {
        ++summary.materializationAttempted;
    }
    if (result.materializationApplied) {
        ++summary.materializationApplied;
    }
    if (result.localCorrectnessValidationPerformed) {
        ++summary.localCorrectnessValidationPerformed;
        if (result.localCorrectnessValidationPassed) {
            ++summary.localCorrectnessValidationPassed;
        }
        if (result.rejectedByLocalCorrectnessValidation) {
            ++summary.rejectedByLocalCorrectnessValidation;
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H2")) {
                ++summary.localCorrectnessRejectedH2;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H5")) {
                ++summary.localCorrectnessRejectedH5;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H6")) {
                ++summary.localCorrectnessRejectedH6;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H7")) {
                ++summary.localCorrectnessRejectedH7;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H8")) {
                ++summary.localCorrectnessRejectedH8;
            }
        }
    }
    summary.sameLevelEdgesCollapsed += result.sameLevelEdgesCollapsed;
    if (result.structuralValidationPerformed) {
        ++summary.structuralValidationPerformed;
        if (result.structuralMatchedRebuild) {
            ++summary.structuralMatchedRebuild;
        } else {
            ++summary.structuralMismatchedRebuild;
        }
    }

    if (!result.contract.acceptsPublishedProblem()) {
        ++summary.rejectedByContract;
        return;
    }

    ++summary.acceptedByPublishedProblem;
    if (countDirectUpdateMethod &&
        result.status != ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation) {
        ++summary.directUpdateAttempted;
    }
    if (result.contract.readyForLocalUpdater()) {
        ++summary.compatibleWithCurrentUpdater;
    }
    const auto &activeContract = result.hasPostMaterializationContract
        ? result.postMaterializationContract
        : result.contract;
    summary.affectedIntervalLevels += activeContract.affectedIntervalLevels;
    summary.missingLevelUnionCandidates += activeContract.missingLevelUnionCandidates;
    summary.missingIntermediateLevelUnionCandidates +=
        activeContract.missingIntermediateLevelUnionCandidates;

    switch (result.status) {
        case ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation:
            ++summary.knownUpdaterLimitations;
            if (result.contract.requiresLocalMaterialization()) {
                ++summary.requiresMaterialization;
            } else if (result.contract.requiresTargetNodeCreation()) {
                ++summary.requiresTargetNodeCreation;
            } else if (result.rejectedByLocalCorrectnessValidation) {
                ++summary.otherKnownLimitations;
            } else {
                ++summary.otherKnownLimitations;
            }
            break;
        case ExplicitRegionProblem1OracleStatus::LocalMaterializationThrew:
            ++summary.localMaterializationThrew;
            break;
        case ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild:
            ++summary.localUpdateMatchedRebuild;
            if (countDirectUpdateMethod) {
                ++summary.directUpdateApplied;
                ++summary.directUpdateMatchedRebuild;
            }
            if (result.materializationApplied) {
                ++summary.materializedLocalUpdateMatchedRebuild;
            }
            if (result.targetNodeCreationAllowed) {
                ++summary.targetNodeCreationLocalUpdateMatchedRebuild;
            }
            if (result.materializationApplied &&
                result.targetNodeCreationAllowed) {
                ++summary.materializedTargetNodeCreationLocalUpdateMatchedRebuild;
            }
            break;
        case ExplicitRegionProblem1OracleStatus::LocalUpdateMismatchedRebuild:
            ++summary.localUpdateMismatchedRebuild;
            if (countDirectUpdateMethod) {
                ++summary.directUpdateApplied;
                ++summary.directUpdateMismatchedRebuild;
            }
            break;
        case ExplicitRegionProblem1OracleStatus::LocalUpdateThrew:
            ++summary.localUpdateThrew;
            if (countDirectUpdateMethod) {
                ++summary.directUpdateThrew;
            }
            break;
        case ExplicitRegionProblem1OracleStatus::RejectedByContract:
            ++summary.rejectedByContract;
            --summary.acceptedByPublishedProblem;
            break;
    }
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1OracleResult<PixelType> runDirectUpdateTreeWithExplicitRegionSurveyEvent(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false,
    bool allowLocalMaterialization = true,
    bool allowTargetNodeCreation = true,
    bool validateStructureAgainstRebuild = false) {
    ExplicitRegionProblem1OracleResult<PixelType> out;

    if (!input) {
        out.message = "direct update survey requires a valid input image";
        return out;
    }
    if (!adjacency) {
        out.message = "direct update survey requires a valid adjacency relation";
        return out;
    }

    DynamicComponentTree minTree(input->clone(), false, adjacency);
    DynamicComponentTree maxTree(input->clone(), true, adjacency);
    DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;

    out.contract = checkExplicitRegionProblem1Contract<PixelType>(
        dualTree,
        *adjacency,
        pixelsC,
        b,
        assumeConnected);

    if (!out.contract.acceptsPublishedProblem()) {
        out.status = ExplicitRegionProblem1OracleStatus::RejectedByContract;
        out.message = out.contract.rejectionReason;
        return out;
    }

    RegionLevelChange<PixelType> change;
    change.pixelsC.assign(pixelsC.begin(), pixelsC.end());
    change.a = out.contract.a;
    change.b = b;
    change.direction = out.contract.direction;
    change.assumeConnected = true;
    change.assumeMaterializedAsProperParts =
        !out.contract.requiresLocalMaterialization();
    out.expectedImage = RegionLevelChangeRebuildUpdater::applyToImage(input, change);

    ExplicitRegionProblem1LocalCorrectnessValidationOptions validationOptions;
    validationOptions.assumeConnected = assumeConnected;
    validationOptions.allowLocalMaterialization = allowLocalMaterialization;
    validationOptions.allowTargetNodeCreation = allowTargetNodeCreation;
    const auto diagnostic =
        diagnoseExplicitRegionProblem1LocalCorrectness<PixelType>(
            input,
            adjacency,
            updateMaxTree,
            pixelsC,
            b,
            validationOptions);
    out.localCorrectnessValidationPerformed = true;
    out.localCorrectnessValidation = diagnostic.validation;
    out.localCorrectnessValidationPassed =
        diagnostic.validation.satisfiesIntegratedReleaseCriteria();
    out.rejectedByLocalCorrectnessValidation =
        !out.localCorrectnessValidationPassed;

    if (diagnostic.analysis.hasPostMaterializationContract) {
        out.postMaterializationContract =
            diagnostic.analysis.postMaterializationContract;
        out.hasPostMaterializationContract = true;
    }

    const ExplicitRegionProblem1Contract<PixelType> &activeContract =
        diagnostic.analysis.activeContract();
    out.materializationAttempted =
        allowLocalMaterialization && out.contract.requiresLocalMaterialization();
    out.materializationApplied =
        out.materializationAttempted &&
        diagnostic.analysis.materializationApplied;
    out.targetNodeCreationRequired =
        explicitRegionProblem1OracleNeedsOnlyTargetNodeCreation(activeContract);
    out.targetNodeCreationAllowed =
        out.targetNodeCreationRequired && allowTargetNodeCreation;

    if (!out.localCorrectnessValidationPassed) {
        out.status = ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation;
        out.message = diagnostic.validation.message;
        return out;
    }

    try {
        DualMinMaxTreeIncrementalFilter<PixelType> adjust(&minTree, &maxTree, *adjacency);
        adjust.setRuntimePostConditionValidationEnabled(true);
        updateTreeWithExplicitRegionChecked(
            adjust,
            &dualTree,
            *adjacency,
            pixelsC,
            b,
            allowTargetNodeCreation);
        out.sameLevelEdgesCollapsed = adjust.getLastSameLevelEdgesCollapsed();
        out.localReconstruction = dualTree.reconstructionImage();
    } catch (const std::exception &error) {
        out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateThrew;
        out.message = error.what();
        return out;
    }

    if (out.localReconstruction->isEqual(out.expectedImage)) {
        if (validateStructureAgainstRebuild) {
            DynamicComponentTree expectedTree(
                out.expectedImage->clone(),
                updateMaxTree,
                adjacency);
            const auto structuralComparison =
                compareDynamicComponentTreeStructure(dualTree, expectedTree);
            out.structuralValidationPerformed = true;
            out.structuralMatchedRebuild = structuralComparison.equivalent;
            out.structuralValidationMessage = structuralComparison.message;
            if (!structuralComparison.equivalent) {
                out.status =
                    ExplicitRegionProblem1OracleStatus::LocalUpdateMismatchedRebuild;
                out.message =
                    "direct update reconstruction matches rebuild, but structure differs: " +
                    structuralComparison.message;
                return out;
            }
        }

        out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild;
        out.message = "direct updateTreeWithExplicitRegion reconstruction matches rebuild oracle";
    } else {
        out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateMismatchedRebuild;
        out.message = "direct updateTreeWithExplicitRegion reconstruction differs from rebuild oracle";
    }

    return out;
}

template<typename PixelType = AltitudeType>
void accumulateExplicitRegionProblem1SafeUpdateDiagnostics(
    ExplicitRegionProblem1SurveySummary &summary,
    const ExplicitRegionProblem1UpdateDiagnostics<PixelType> &result,
    const ImageUInt8Ptr &input,
    DynamicComponentTree &minTree,
    DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool validateStructureAgainstRebuild) {
    if (result.localCorrectnessValidationPerformed) {
        ++summary.localCorrectnessValidationPerformed;
        if (result.localCorrectnessValidationPassed) {
            ++summary.localCorrectnessValidationPassed;
        }
        if (result.rejectedByLocalCorrectnessValidation) {
            ++summary.rejectedByLocalCorrectnessValidation;
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H2")) {
                ++summary.localCorrectnessRejectedH2;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H5")) {
                ++summary.localCorrectnessRejectedH5;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H6")) {
                ++summary.localCorrectnessRejectedH6;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H7")) {
                ++summary.localCorrectnessRejectedH7;
            }
            if (localCorrectnessFailedHypothesis(
                    result.localCorrectnessValidation,
                    "H8")) {
                ++summary.localCorrectnessRejectedH8;
            }
        }
    }

    summary.sameLevelEdgesCollapsed += result.sameLevelEdgesCollapsed;
    if (result.structuralValidationPerformed) {
        ++summary.structuralValidationPerformed;
        if (result.structuralMatchedRebuild) {
            ++summary.structuralMatchedRebuild;
        } else {
            ++summary.structuralMismatchedRebuild;
        }
    }

    if (!result.contract.acceptsPublishedProblem()) {
        ++summary.rejectedByContract;
        return;
    }

    ++summary.acceptedByPublishedProblem;
    ++summary.safeUpdateAttempted;
    if (result.contract.readyForLocalUpdater()) {
        ++summary.compatibleWithCurrentUpdater;
    }
    const auto &activeContract = result.hasPostMaterializationContract
        ? result.postMaterializationContract
        : result.contract;
    summary.affectedIntervalLevels += activeContract.affectedIntervalLevels;
    summary.missingLevelUnionCandidates += activeContract.missingLevelUnionCandidates;
    summary.missingIntermediateLevelUnionCandidates +=
        activeContract.missingIntermediateLevelUnionCandidates;

    if (result.materializationAttempted) {
        ++summary.materializationAttempted;
    }
    if (result.materializationApplied) {
        ++summary.materializationApplied;
    }
    if (result.fallbackAttempted) {
        ++summary.rebuildFallbackAttempted;
    }

    switch (result.status) {
        case ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation:
            ++summary.knownUpdaterLimitations;
            if (result.contract.requiresLocalMaterialization()) {
                ++summary.requiresMaterialization;
            } else if (activeContract.requiresTargetNodeCreation()) {
                ++summary.requiresTargetNodeCreation;
            } else {
                ++summary.otherKnownLimitations;
            }
            break;
        case ExplicitRegionProblem1UpdateStatus::LocalMaterializationThrew:
            ++summary.localMaterializationThrew;
            break;
        case ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied:
            ++summary.safeUpdateApplied;
            ++summary.safeUpdateMatchedRebuild;
            ++summary.localUpdateMatchedRebuild;
            if (result.materializationApplied) {
                ++summary.materializedLocalUpdateMatchedRebuild;
            }
            if (result.targetNodeCreationAllowed) {
                ++summary.targetNodeCreationLocalUpdateMatchedRebuild;
            }
            if (result.materializationApplied &&
                result.targetNodeCreationAllowed) {
                ++summary.materializedTargetNodeCreationLocalUpdateMatchedRebuild;
            }
            break;
        case ExplicitRegionProblem1UpdateStatus::LocalUpdateMismatchedRebuild:
            ++summary.safeUpdateApplied;
            ++summary.safeUpdateMismatchedRebuild;
            ++summary.localUpdateMismatchedRebuild;
            break;
        case ExplicitRegionProblem1UpdateStatus::LocalUpdateThrew:
            ++summary.localUpdateThrew;
            break;
        case ExplicitRegionProblem1UpdateStatus::RebuildFallbackApplied: {
            ++summary.safeUpdateApplied;
            ++summary.rebuildFallbackApplied;

            const auto change = makeExplicitRegionProblem1ChangeForFallback<PixelType>(
                result.contract,
                pixelsC,
                b);
            const auto expectedImage = RegionLevelChangeRebuildUpdater::applyToImage(
                input,
                change);

            const bool reconstructionMatched =
                expectedImage &&
                result.fallbackReconstruction &&
                result.fallbackReconstruction->isEqual(expectedImage) &&
                minTree.reconstructionImage()->isEqual(expectedImage) &&
                maxTree.reconstructionImage()->isEqual(expectedImage);
            bool structuralMatched = true;
            if (validateStructureAgainstRebuild && expectedImage) {
                DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;
                DynamicComponentTree expectedTree(
                    expectedImage->clone(),
                    updateMaxTree,
                    adjacency);
                const auto structuralComparison =
                    compareDynamicComponentTreeStructure(dualTree, expectedTree);
                ++summary.structuralValidationPerformed;
                structuralMatched = structuralComparison.equivalent;
                if (structuralMatched) {
                    ++summary.structuralMatchedRebuild;
                } else {
                    ++summary.structuralMismatchedRebuild;
                }
            }

            if (reconstructionMatched && structuralMatched) {
                ++summary.safeUpdateMatchedRebuild;
                ++summary.rebuildFallbackMatchedRebuild;
            } else {
                ++summary.safeUpdateMismatchedRebuild;
                ++summary.rebuildFallbackMismatchedRebuild;
            }
            break;
        }
        case ExplicitRegionProblem1UpdateStatus::RebuildFallbackThrew:
            ++summary.rebuildFallbackThrew;
            break;
        case ExplicitRegionProblem1UpdateStatus::RejectedByContract:
            ++summary.rejectedByContract;
            --summary.acceptedByPublishedProblem;
            --summary.safeUpdateAttempted;
            break;
    }
}

template<typename PixelType = AltitudeType>
void accumulateExplicitRegionProblem1MultiAnchorSurveyResult(
    ExplicitRegionProblem1MultiAnchorSurveySummary &summary,
    const ExplicitRegionProblem1MultiAnchorAnalysis<PixelType> &analysis,
    bool allowTargetNodeCreation = true) {
    if (!analysis.acceptedByPublishedProblem) {
        ++summary.rejectedByContract;
        return;
    }

    ++summary.acceptedByPublishedProblem;
    if (analysis.materializationAttempted) {
        ++summary.materializationAttempted;
    }
    if (analysis.materializationApplied) {
        ++summary.materializationApplied;
    }
    if (!analysis.analyzed) {
        return;
    }

    const auto validation =
        validateExplicitRegionProblem1LocalCorrectnessFromAnalysis<PixelType>(
            analysis,
            allowTargetNodeCreation);

    ++summary.analyzedEvents;
    summary.owners += analysis.ownerCount;
    summary.extremalOwners += analysis.extremalOwnerCount;
    summary.ownerPartitionComponents += analysis.ownerPartitionComponents;
    summary.disconnectedOwnerPartitions += analysis.disconnectedOwnerPartitions;
    summary.ownerIncidenceEdges += analysis.ownerIncidenceEdges;
    summary.sparseWorkspaceNodes += analysis.sparseWorkspaceNodeCount;
    summary.sparseWorkspaceAnchors += analysis.sparseWorkspaceAnchorCount;
    summary.ownerIntervalWorkspaceNodes +=
        analysis.ownerIntervalWorkspaceNodeCount;
    summary.ownerIntervalWorkspaceAnchors +=
        analysis.ownerIntervalWorkspaceAnchorCount;
    summary.generalWorkspaceNodes += analysis.generalWorkspaceNodeCount;
    summary.generalWorkspaceAnchors += analysis.generalWorkspaceAnchorCount;
    summary.ownerIntervalBucketMissingNodes +=
        analysis.ownerIntervalBucketMissingNodeCount;
    summary.ownerIntervalBucketClosurePathMissingNodes +=
        analysis.ownerIntervalBucketClosurePathMissingNodeCount;
    summary.ownerIntervalBucketCOnlyProperPartMissingNodes +=
        analysis.ownerIntervalBucketCOnlyProperPartMissingNodeCount;
    summary.ownerIntervalBucketUnexplainedMissingNodes +=
        analysis.ownerIntervalBucketUnexplainedMissingNodeCount;
    summary.frontierNodes += analysis.frontierNodeCount;

    summary.maxOwners = std::max(summary.maxOwners, analysis.ownerCount);
    summary.maxExtremalOwners =
        std::max(summary.maxExtremalOwners, analysis.extremalOwnerCount);
    summary.maxOwnerPartitionComponents =
        std::max(summary.maxOwnerPartitionComponents,
                 analysis.ownerPartitionComponents);
    summary.maxOwnerIncidenceEdges =
        std::max(summary.maxOwnerIncidenceEdges, analysis.ownerIncidenceEdges);
    summary.maxSparseWorkspaceAnchors =
        std::max(summary.maxSparseWorkspaceAnchors,
                 analysis.sparseWorkspaceAnchorCount);
    summary.maxOwnerIntervalWorkspaceAnchors =
        std::max(summary.maxOwnerIntervalWorkspaceAnchors,
                 analysis.ownerIntervalWorkspaceAnchorCount);
    summary.maxGeneralWorkspaceAnchors =
        std::max(summary.maxGeneralWorkspaceAnchors,
                 analysis.generalWorkspaceAnchorCount);
    if (analysis.nodeCaReachesClosureRoot()) {
        summary.maxNodeCaToClosureRootPathLength =
            std::max(summary.maxNodeCaToClosureRootPathLength,
                     analysis.nodeCaToClosureRootPathLength);
    }
    summary.maxOwnerIntervalBucketMissingNodes =
        std::max(summary.maxOwnerIntervalBucketMissingNodes,
                 analysis.ownerIntervalBucketMissingNodeCount);
    summary.maxOwnerIntervalBucketUnexplainedMissingNodes =
        std::max(summary.maxOwnerIntervalBucketUnexplainedMissingNodes,
                 analysis.ownerIntervalBucketUnexplainedMissingNodeCount);

    if (analysis.hasMultipleExtremalOwners()) {
        ++summary.multipleExtremalOwnerEvents;
    }
    if (analysis.hasMultipleOwners()) {
        ++summary.multipleOwnerEvents;
    }
    if (analysis.isMultiOwnerSingleWorkspaceCandidate()) {
        ++summary.multiOwnerSingleWorkspaceEvents;
    }
    if (analysis.hasDisconnectedOwnerPartition()) {
        ++summary.disconnectedOwnerPartitionEvents;
    }
    if (analysis.hasMultipleSparseWorkspaceAnchors()) {
        ++summary.multipleSparseWorkspaceAnchorEvents;
    }
    if (analysis.hasMultipleGeneralWorkspaceAnchors()) {
        ++summary.multipleGeneralWorkspaceAnchorEvents;
    }
    if (analysis.hasClosureRoot()) {
        ++summary.closureRootEvents;
        if (analysis.nodeCaReachesClosureRoot()) {
            ++summary.nodeCaReachesClosureRootEvents;
            summary.nodeCaToClosureRootPathLengthSum +=
                analysis.nodeCaToClosureRootPathLength;
            if (analysis.nodeCaToClosureRootSameLevelPath) {
                ++summary.nodeCaToClosureRootSameLevelPathEvents;
            }
        }
        if (analysis.closureRootMatchesActiveNodeCa()) {
            ++summary.closureRootMatchesNodeCaEvents;
        } else if (analysis.closureRootDiffersFromActiveNodeCa()) {
            ++summary.closureRootDiffersFromNodeCaEvents;
            if (analysis.hasMultipleOwners()) {
                ++summary.multiOwnerClosureRootDiffersFromNodeCaEvents;
            }
            if (analysis.nodeCaToClosureRootSameLevelPath) {
                ++summary.closureRootDiffersButSameLevelPathEvents;
            }
        }
    }
    if (analysis.ownerIntervalBucketMissingNodeCount > 0) {
        ++summary.ownerIntervalBucketGapEvents;
    }
    if (analysis.hasUnexplainedOwnerIntervalBucketGap()) {
        ++summary.ownerIntervalBucketUnexplainedGapEvents;
    }
    if (validation.h6UniqueGeneralWorkspaceRoot) {
        ++summary.h6UniqueGeneralWorkspaceRootEvents;
    }
    if (validation.h7SameLevelClosure) {
        ++summary.h7SameLevelClosureEvents;
    }
    if (validation.h8BucketCompleteness) {
        ++summary.h8BucketCompletenessEvents;
    }
    if (validation.targetLevelPresentOrCreated) {
        ++summary.targetLevelPresentOrCreatedEvents;
    }
    if (validation.satisfiesIntegratedReleaseCriteria()) {
        ++summary.localCorrectnessCandidateEvents;
    }
    if (analysis.isMultiAnchorCandidate()) {
        ++summary.multiAnchorCandidateEvents;
    }
}

} // namespace detail

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1SurveySummary surveyExplicitRegionProblem1Events(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    int maxRegionSize,
    PixelType minB,
    PixelType maxB,
    bool allowLocalMaterialization = false,
    bool allowTargetNodeCreation = false,
    bool validateStructureAgainstRebuild = false,
    bool enforceLocalCorrectnessValidation = false,
    bool useSafeUpdateWrapper = false,
    bool enableRebuildFallback = false,
    bool useDirectUpdateMethod = false) {
    ExplicitRegionProblem1SurveySummary summary;
    if (!input || !adjacency || minB > maxB) {
        return summary;
    }

    auto connectedRegions = detail::enumerateConnectedRegions(
        *adjacency,
        input->getSize(),
        maxRegionSize);
    summary.connectedRegions = static_cast<int>(connectedRegions.size());

    for (const auto &region : connectedRegions) {
        for (int bValue = static_cast<int>(minB);
             bValue <= static_cast<int>(maxB);
             ++bValue) {
            const PixelType b = static_cast<PixelType>(bValue);
            ++summary.totalEvents;
            if (useSafeUpdateWrapper) {
                DynamicComponentTree minTree(input->clone(), false, adjacency);
                DynamicComponentTree maxTree(input->clone(), true, adjacency);
                ExplicitRegionProblem1UpdateOptions updateOptions;
                updateOptions.assumeConnected = true;
                updateOptions.allowLocalMaterialization = allowLocalMaterialization;
                updateOptions.allowTargetNodeCreation = allowTargetNodeCreation;
                updateOptions.enforceLocalCorrectnessValidation =
                    enforceLocalCorrectnessValidation;
                updateOptions.validateAgainstRebuild = true;
                updateOptions.validateStructureAgainstRebuild =
                    validateStructureAgainstRebuild;
                updateOptions.fallbackPolicy = enableRebuildFallback
                    ? ExplicitRegionProblem1FallbackPolicy::RebuildPair
                    : ExplicitRegionProblem1FallbackPolicy::NoFallback;

                const auto result = safeUpdateExplicitRegionProblem1WithDiagnostics<PixelType>(
                    minTree,
                    maxTree,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    updateOptions);
                detail::accumulateExplicitRegionProblem1SafeUpdateDiagnostics(
                    summary,
                    result,
                    input,
                    minTree,
                    maxTree,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    validateStructureAgainstRebuild);
            } else if (useDirectUpdateMethod) {
                const auto result =
                    detail::runDirectUpdateTreeWithExplicitRegionSurveyEvent<PixelType>(
                        input,
                        adjacency,
                        updateMaxTree,
                        std::span<const PixelId>(region),
                        b,
                        true,
                        allowLocalMaterialization,
                        allowTargetNodeCreation,
                        validateStructureAgainstRebuild);
                detail::accumulateExplicitRegionProblem1SurveyResult(
                    summary,
                    result,
                    true);
            } else {
                const auto result = runExplicitRegionProblem1Oracle<PixelType>(
                    input,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    true,
                    allowLocalMaterialization,
                    allowTargetNodeCreation,
                    validateStructureAgainstRebuild,
                    enforceLocalCorrectnessValidation);

                detail::accumulateExplicitRegionProblem1SurveyResult(summary, result);
            }
        }
    }

    return summary;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1SurveySummary surveySeededExplicitRegionProblem1Events(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    int maxRegionSize,
    int maxRegions,
    PixelType minB,
    PixelType maxB,
    bool allowLocalMaterialization = false,
    bool allowTargetNodeCreation = false,
    bool validateStructureAgainstRebuild = false,
    bool enforceLocalCorrectnessValidation = false,
    bool useSafeUpdateWrapper = false,
    bool enableRebuildFallback = false,
    bool useDirectUpdateMethod = false) {
    ExplicitRegionProblem1SurveySummary summary;
    if (!input || !adjacency || minB > maxB) {
        return summary;
    }

    auto connectedRegions = detail::enumerateSeededConnectedRegions(
        *adjacency,
        input->getSize(),
        maxRegionSize,
        maxRegions);
    summary.connectedRegions = static_cast<int>(connectedRegions.size());

    for (const auto &region : connectedRegions) {
        for (int bValue = static_cast<int>(minB);
             bValue <= static_cast<int>(maxB);
             ++bValue) {
            const PixelType b = static_cast<PixelType>(bValue);
            ++summary.totalEvents;
            if (useSafeUpdateWrapper) {
                DynamicComponentTree minTree(input->clone(), false, adjacency);
                DynamicComponentTree maxTree(input->clone(), true, adjacency);
                ExplicitRegionProblem1UpdateOptions updateOptions;
                updateOptions.assumeConnected = true;
                updateOptions.allowLocalMaterialization = allowLocalMaterialization;
                updateOptions.allowTargetNodeCreation = allowTargetNodeCreation;
                updateOptions.enforceLocalCorrectnessValidation =
                    enforceLocalCorrectnessValidation;
                updateOptions.validateAgainstRebuild = true;
                updateOptions.validateStructureAgainstRebuild =
                    validateStructureAgainstRebuild;
                updateOptions.fallbackPolicy = enableRebuildFallback
                    ? ExplicitRegionProblem1FallbackPolicy::RebuildPair
                    : ExplicitRegionProblem1FallbackPolicy::NoFallback;

                const auto result = safeUpdateExplicitRegionProblem1WithDiagnostics<PixelType>(
                    minTree,
                    maxTree,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    updateOptions);
                detail::accumulateExplicitRegionProblem1SafeUpdateDiagnostics(
                    summary,
                    result,
                    input,
                    minTree,
                    maxTree,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    validateStructureAgainstRebuild);
            } else if (useDirectUpdateMethod) {
                const auto result =
                    detail::runDirectUpdateTreeWithExplicitRegionSurveyEvent<PixelType>(
                        input,
                        adjacency,
                        updateMaxTree,
                        std::span<const PixelId>(region),
                        b,
                        true,
                        allowLocalMaterialization,
                        allowTargetNodeCreation,
                        validateStructureAgainstRebuild);
                detail::accumulateExplicitRegionProblem1SurveyResult(
                    summary,
                    result,
                    true);
            } else {
                const auto result = runExplicitRegionProblem1Oracle<PixelType>(
                    input,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    true,
                    allowLocalMaterialization,
                    allowTargetNodeCreation,
                    validateStructureAgainstRebuild,
                    enforceLocalCorrectnessValidation);
                detail::accumulateExplicitRegionProblem1SurveyResult(summary, result);
            }
        }
    }

    return summary;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1MultiAnchorSurveySummary
surveyExplicitRegionProblem1MultiAnchorEvents(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    int maxRegionSize,
    PixelType minB,
    PixelType maxB,
    bool allowLocalMaterialization = true,
    bool allowTargetNodeCreation = true) {
    ExplicitRegionProblem1MultiAnchorSurveySummary summary;
    if (!input || !adjacency || minB > maxB) {
        return summary;
    }

    auto connectedRegions = detail::enumerateConnectedRegions(
        *adjacency,
        input->getSize(),
        maxRegionSize);
    summary.connectedRegions = static_cast<int>(connectedRegions.size());

    for (const auto &region : connectedRegions) {
        for (int bValue = static_cast<int>(minB);
             bValue <= static_cast<int>(maxB);
             ++bValue) {
            const PixelType b = static_cast<PixelType>(bValue);
            ++summary.totalEvents;
            const auto analysis =
                analyzeExplicitRegionProblem1MultiAnchor<PixelType>(
                    input,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    true,
                    allowLocalMaterialization);
            detail::accumulateExplicitRegionProblem1MultiAnchorSurveyResult(
                summary,
                analysis,
                allowTargetNodeCreation);
        }
    }

    return summary;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1MultiAnchorSurveySummary
surveySeededExplicitRegionProblem1MultiAnchorEvents(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    int maxRegionSize,
    int maxRegions,
    PixelType minB,
    PixelType maxB,
    bool allowLocalMaterialization = true,
    bool allowTargetNodeCreation = true) {
    ExplicitRegionProblem1MultiAnchorSurveySummary summary;
    if (!input || !adjacency || minB > maxB) {
        return summary;
    }

    auto connectedRegions = detail::enumerateSeededConnectedRegions(
        *adjacency,
        input->getSize(),
        maxRegionSize,
        maxRegions);
    summary.connectedRegions = static_cast<int>(connectedRegions.size());

    for (const auto &region : connectedRegions) {
        for (int bValue = static_cast<int>(minB);
             bValue <= static_cast<int>(maxB);
             ++bValue) {
            const PixelType b = static_cast<PixelType>(bValue);
            ++summary.totalEvents;
            const auto analysis =
                analyzeExplicitRegionProblem1MultiAnchor<PixelType>(
                    input,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    b,
                    true,
                    allowLocalMaterialization);
            detail::accumulateExplicitRegionProblem1MultiAnchorSurveyResult(
                summary,
                analysis,
                allowTargetNodeCreation);
        }
    }

    return summary;
}

} // namespace debug::problem1
