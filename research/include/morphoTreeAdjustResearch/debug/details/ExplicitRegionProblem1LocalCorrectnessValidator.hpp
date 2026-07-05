#pragma once

#include <span>
#include <string>
#include <vector>

#include "ExplicitRegionProblem1MultiAnchor.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"

namespace debug::problem1 {

struct ExplicitRegionProblem1LocalCorrectnessValidationOptions {
    bool assumeConnected = false;
    bool allowLocalMaterialization = true;
    bool allowTargetNodeCreation = true;
};

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1LocalCorrectnessValidation {
    bool acceptedByPublishedProblem = false;
    bool analyzed = false;
    bool materializedOwnersReady = false;
    bool targetLevelPresent = false;
    bool targetLevelCreationAllowed = false;
    bool targetLevelPresentOrCreated = false;
    bool uniqueGeneralWorkspaceRoot = false;
    bool nodeCaReachesClosureRoot = false;
    // H7: nodeCa and closureRoot must represent the same level-a threshold
    // component, possibly split by same-level dynamic-tree nodes.
    bool nodeCaToClosureRootSameLevelPath = false;
    // H8: every owner-interval node omitted from concrete merge buckets must
    // be explained by the same-level closure path or by direct C-only support.
    bool noUnexplainedOwnerIntervalBucketGap = false;
    bool h6UniqueGeneralWorkspaceRoot = false;
    bool h7SameLevelClosure = false;
    bool h8BucketCompleteness = false;
    NodeId closureRoot = InvalidNode;
    NodeId nodeCa = InvalidNode;
    int generalWorkspaceAnchorCount = 0;
    int ownerIntervalBucketUnexplainedMissingNodeCount = 0;
    int ownerIntervalBucketMissingNodeCount = 0;
    std::vector<std::string> failedHypotheses;
    std::string message;

    bool satisfiesH6ToH8() const {
        return h6UniqueGeneralWorkspaceRoot &&
               h7SameLevelClosure &&
               h8BucketCompleteness;
    }

    bool satisfiesIntegratedReleaseCriteria() const {
        return acceptedByPublishedProblem &&
               analyzed &&
               materializedOwnersReady &&
               targetLevelPresentOrCreated &&
               satisfiesH6ToH8();
    }
};

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1LocalCorrectnessDiagnostic {
    ExplicitRegionProblem1MultiAnchorAnalysis<PixelType> analysis;
    ExplicitRegionProblem1LocalCorrectnessValidation<PixelType> validation;
};

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1LocalCorrectnessValidation<PixelType>
validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
    const ExplicitRegionProblem1MultiAnchorAnalysis<PixelType> &analysis,
    bool allowTargetNodeCreation = true) {
    ExplicitRegionProblem1LocalCorrectnessValidation<PixelType> out;

    out.acceptedByPublishedProblem = analysis.acceptedByPublishedProblem;
    out.analyzed = analysis.analyzed;
    out.closureRoot = analysis.closureRoot;
    out.nodeCa = analysis.activeContract().nodeCa;
    out.generalWorkspaceAnchorCount = analysis.generalWorkspaceAnchorCount;
    out.ownerIntervalBucketMissingNodeCount =
        analysis.ownerIntervalBucketMissingNodeCount;
    out.ownerIntervalBucketUnexplainedMissingNodeCount =
        analysis.ownerIntervalBucketUnexplainedMissingNodeCount;

    if (!analysis.acceptedByPublishedProblem) {
        out.failedHypotheses.push_back("H1");
        out.message = analysis.message.empty()
            ? "event rejected by the published Problem 1 contract"
            : analysis.message;
        return out;
    }

    if (!analysis.analyzed) {
        out.failedHypotheses.push_back("analysis");
        out.message = analysis.message.empty()
            ? "local correctness analysis was not completed"
            : analysis.message;
        return out;
    }

    const auto &contract = analysis.activeContract();
    out.materializedOwnersReady = !contract.requiresLocalMaterialization();
    out.targetLevelPresent = contract.diagnostics.hasTargetLevelUnionCandidate;
    out.targetLevelCreationAllowed =
        !out.targetLevelPresent && allowTargetNodeCreation;
    out.targetLevelPresentOrCreated =
        out.targetLevelPresent || out.targetLevelCreationAllowed;
    out.uniqueGeneralWorkspaceRoot =
        analysis.generalWorkspaceAnchorCount == 1 && analysis.hasClosureRoot();
    out.nodeCaReachesClosureRoot = analysis.nodeCaReachesClosureRoot();
    out.nodeCaToClosureRootSameLevelPath =
        analysis.nodeCaToClosureRootSameLevelPath;
    out.noUnexplainedOwnerIntervalBucketGap =
        !analysis.hasUnexplainedOwnerIntervalBucketGap();

    out.h6UniqueGeneralWorkspaceRoot = out.uniqueGeneralWorkspaceRoot;
    out.h7SameLevelClosure =
        out.nodeCaReachesClosureRoot &&
        out.nodeCaToClosureRootSameLevelPath;
    out.h8BucketCompleteness = out.noUnexplainedOwnerIntervalBucketGap;

    if (!out.materializedOwnersReady) {
        out.failedHypotheses.push_back("H2");
    }
    if (!out.targetLevelPresentOrCreated) {
        out.failedHypotheses.push_back("H5");
    }
    if (!out.h6UniqueGeneralWorkspaceRoot) {
        out.failedHypotheses.push_back("H6");
    }
    if (!out.h7SameLevelClosure) {
        out.failedHypotheses.push_back("H7");
    }
    if (!out.h8BucketCompleteness) {
        out.failedHypotheses.push_back("H8");
    }

    out.message = out.satisfiesIntegratedReleaseCriteria()
        ? "event satisfies the integrated Problem 1 local correctness release criteria"
        : "event violates one or more integrated Problem 1 local correctness release criteria";
    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1LocalCorrectnessDiagnostic<PixelType>
diagnoseExplicitRegionProblem1LocalCorrectness(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    const ExplicitRegionProblem1LocalCorrectnessValidationOptions &options = {}) {
    ExplicitRegionProblem1LocalCorrectnessDiagnostic<PixelType> out;
    out.analysis = analyzeExplicitRegionProblem1MultiAnchor<PixelType>(
        input,
        adjacency,
        updateMaxTree,
        pixelsC,
        b,
        options.assumeConnected,
        options.allowLocalMaterialization);
    out.validation = validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
        out.analysis,
        options.allowTargetNodeCreation);
    return out;
}

} // namespace debug::problem1
