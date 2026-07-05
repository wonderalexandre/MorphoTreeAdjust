#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/debug/NonIncreasingCasfDiagnostics.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingCasfKernel.hpp"
#include "dynamic_casf_apply_common.hpp"

namespace common = dynamic_casf_apply_common;
namespace noninc = morphotree::details::non_increasing;
namespace noninc_debug = debug::non_increasing;

namespace {

enum class Mode {
    ExplicitUpdate,
    Naive,
    Compare,
};

const char *phaseName(noninc::NonIncreasingPhase phase);
void printIntVectorPreview(const std::vector<int> &values);
void printCandidateContext(const char *prefix,
                           const noninc::NonIncreasingCandidateHierarchyContext &context);

struct MaxCanonicalizationEventStats {
    bool present = false;
    long long micros = 0;
    float threshold = 0.0f;
    noninc::NonIncreasingPhase phase = noninc::NonIncreasingPhase::MaxToMin;
    NodeId nodeId = InvalidNode;
    NodeId parentId = InvalidNode;
    int nodeLevel = 0;
    int parentLevel = 0;
    int numChildren = 0;
    int componentCount = 0;
    int totalComponentPixels = 0;
    int strippedLocalPixelCount = 0;
    int strippedLocalCreatedNodes = 0;
    int strippedLocalRecompositionNodes = 0;
    int strippedCreatedConnectorNodes = 0;
    int strippedEmptyCreatedNodes = 0;
    int strippedEmptyConnectorNodes = 0;
    long long materializationMicros = 0;
    long long releaseMicros = 0;
    long long targetCheckMicros = 0;
};

struct MaxPrimalMergeEventStats {
    bool present = false;
    long long micros = 0;
    float threshold = 0.0f;
    noninc::NonIncreasingPhase phase = noninc::NonIncreasingPhase::MaxToMin;
    NodeId nodeId = InvalidNode;
    NodeId parentId = InvalidNode;
    int nodeLevel = 0;
    int parentLevel = 0;
    int directProperParts = 0;
    int numChildren = 0;
    int componentCount = 0;
    int totalComponentPixels = 0;
    int strippedLocalPixelCount = 0;
    int strippedLocalCreatedNodes = 0;
    int strippedLocalRecompositionNodes = 0;
    int strippedCreatedConnectorNodes = 0;
    long long materializationMicros = 0;
    long long canonicalizationMicros = 0;
};

struct MergeEventStats {
    std::size_t events = 0;
    std::size_t maxToMinEvents = 0;
    std::size_t minToMaxEvents = 0;
    std::size_t multiComponentEvents = 0;
    std::size_t childEvents = 0;
    std::size_t dualRebuildEvents = 0;
    std::size_t phaseDualRebuilds = 0;
    std::size_t phaseCandidateGuardRebuilds = 0;
    std::size_t phaseCandidateGuardCanonicalizationChanges = 0;
    std::size_t strippedLocalAttempts = 0;
    std::size_t strippedLocalAccepted = 0;
    std::size_t strippedLocalApplied = 0;
    std::size_t strippedDirectCommits = 0;
    std::size_t strippedLocalPlanAuditChecks = 0;
    std::size_t strippedLocalTargetMatches = 0;
    std::size_t strippedLocalTargetChecks = 0;
    std::size_t strippedLocalOracleFreePlans = 0;
    std::size_t strippedMaterializationAuditChecks = 0;
    std::size_t strippedLocalOracleFreeMaterializations = 0;
    std::size_t strippedFinalReconstructionChecks = 0;
    std::size_t strippedLocalFallbacks = 0;
    std::size_t strippedLocalRebuildFallbacks = 0;
    std::size_t strippedLocalPhaseRebuildFallbacks = 0;
    std::size_t strippedLocalPhaseRebuildSkips = 0;
    std::size_t partialFAttempts = 0;
    std::size_t partialFAccepted = 0;
    std::size_t partialFRejected = 0;
    std::size_t partialFMutationAttempted = 0;
    std::size_t partialFMcpcLiveAccepted = 0;
    std::size_t partialFMcpcHiddenSupportRedundantLiveAccepted = 0;
    std::size_t partialFMcpcParentLinkBridgeLiveAccepted = 0;
    std::size_t partialFMcpcStableChildBridgeLiveAccepted = 0;
    std::size_t partialFMcpcChildLinkBridgeLiveAccepted = 0;
    std::size_t partialFMcpcEmptyOwnerReleaseLiveAccepted = 0;
    std::size_t partialFMcpcStableChildMoveLiveAccepted = 0;
    std::size_t partialFRejectInvalidEvent = 0;
    std::size_t partialFRejectEmptyC = 0;
    std::size_t partialFRejectNonConnectedC = 0;
    std::size_t partialFRejectInvalidPixel = 0;
    std::size_t partialFRejectDuplicatedPixels = 0;
    std::size_t partialFRejectDeadOwnerPixel = 0;
    std::size_t partialFRejectInvalidTargetLevel = 0;
    std::size_t partialFRejectPublicProblem1Preflight = 0;
    std::size_t partialFRejectMutationSafetyUnproven = 0;
    std::size_t partialFRejectD12BasePreflightFailed = 0;
    std::size_t partialFRejectD12PlanIncomplete = 0;
    std::size_t partialFRejectD12PlanOracleNotFree = 0;
    std::size_t partialFRejectD12MaterializationContextMissing = 0;
    std::size_t partialFRejectD12MaterializationOracleNotFree = 0;
    std::size_t partialFRejectD12MaterializationInputsInvalid = 0;
    std::size_t partialFRejectD12RollbackSummaryMissing = 0;
    std::size_t partialFRejectD12InvalidStableChildRef = 0;
    std::size_t partialFRejectD12StableChildCycleRisk = 0;
    std::size_t partialFRejectD12DirectV1PrecheckMissing = 0;
    std::size_t partialFRejectD12HiddenSupportNotRedundant = 0;
    std::size_t partialFRejectD12ParentLinkNotRepresented = 0;
    std::size_t partialFRejectD12ParentLinkInvalidLevel = 0;
    std::size_t partialFRejectD12DirectV1PrecheckRejected = 0;
    std::size_t partialFRejectD12ExceptionFreeUnproven = 0;
    std::size_t partialFRejectD12DryRunEmptyPlan = 0;
    std::size_t partialFRejectD12DryRunAuxiliaryNodesPresent = 0;
    std::size_t partialFRejectD12DryRunStableChildBridgePresent = 0;
    std::size_t partialFRejectD12DryRunNonOracleHiddenSupportConnectorPresent = 0;
    std::size_t partialFRejectD12DryRunParentLinkBridgePresent = 0;
    std::size_t partialFRejectD12DryRunChildLinkBridgePresent = 0;
    std::size_t partialFRejectD12DryRunExternalAttachmentInvalid = 0;
    std::size_t partialFRejectD12DryRunExternalAttachmentEmpty = 0;
    std::size_t partialFRejectD12DryRunStableChildMovePresent = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerPresent = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerInvalidDirectProperPart = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerNonLiveSourceOwner = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerInvalidStableChild = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerRootSourceOwner = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerSourceOwnerNotLeaf = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidParent = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidChild = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidChildLevel = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenCycleRisk = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenDuplicateChild = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidPlannedLevel = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildPresent = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildSnapshotMissing = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildTouchesModified = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildSupportOverlapsPlan = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerWithChildrenExternalDependency = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerSourceCountMismatch = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerRootPreviousParent = 0;
    std::size_t partialFRejectD12DryRunEmptyOwnerPreviousParentCountMismatch = 0;
    std::size_t partialFRejectD12DryRunInvalidParentPlan = 0;
    std::size_t partialFRejectD12DryRunInvalidParentLevel = 0;
    std::size_t partialFRejectD12DryRunEmptyPlannedNode = 0;
    std::size_t partialFRejectD12DryRunChildMovePresent = 0;
    std::size_t partialFRejectD12DryRunInvalidPixel = 0;
    std::size_t partialFRejectD12DryRunDuplicatedPixel = 0;
    std::size_t partialFRejectCurrentUpdaterLimitation = 0;
    std::size_t partialFRejectOther = 0;
    std::size_t mergeNodeLocalAttempts = 0;
    std::size_t mergeNodeLocalAccepted = 0;
    std::size_t mergeNodeLocalApplied = 0;
    std::size_t mergeNodeLocalCertificationChecks = 0;
    std::size_t mergeNodeLocalCertificationPassed = 0;
    std::size_t mergeNodeLocalDirectCommits = 0;
    std::size_t mergeNodeLocalTargetCreationAllowedEvents = 0;
    std::size_t mergeNodeLocalTargetCreationRequiredEvents = 0;
    std::size_t mergeNodeLocalPhaseRebuildFallbacks = 0;
    std::size_t mergeNodeLocalPhaseRebuildSkips = 0;
    std::size_t mergeNodeLocalPreDualMatches = 0;
    std::size_t mergeNodeLocalDualMatches = 0;
    long long mergeNodeLocalRejectedComponents = 0;
    long long mergeNodeLocalProblem1PreflightRejectedComponents = 0;
    long long mergeNodeLocalProblem1PreflightRejectedByPublishedProblem = 0;
    long long mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation = 0;
    long long mergeNodeLocalPartialOwnerComponents = 0;
    long long mergeNodeLocalRawMissingTargetComponents = 0;
    long long mergeNodeLocalMissingTargetComponents = 0;
    long long mergeNodeLocalPartialOwnerAndMissingTargetComponents = 0;
    long long mergeNodeLocalPixels = 0;
    long long mergeNodeLocalCopyMicros = 0;
    long long mergeNodeLocalCertificationMicros = 0;
    long long mergeNodeComponentCollectMicros = 0;
    long long mergeNodeLocalUpdateMicros = 0;
    long long mergeNodePrimalMergeMicros = 0;
    long long mergeNodeEmptyReleaseMicros = 0;
    long long strippedLocalCreatedNodes = 0;
    long long strippedLocalRecompositionNodes = 0;
    long long strippedLocalPixels = 0;
    long long strippedLocalReleasedOwners = 0;
    long long strippedEmptyCreatedNodes = 0;
    long long strippedEmptyConnectorNodes = 0;
    long long strippedLocalCopyMicros = 0;
    long long strippedLocalUpdateMicros = 0;
    long long strippedOracleFallbackMicros = 0;
    long long strippedRebuildFallbackMicros = 0;
    long long strippedInputExtractionMicros = 0;
    long long strippedPartitionPlanMicros = 0;
    long long strippedRecompositionPlanMicros = 0;
    long long strippedPlanAuditMicros = 0;
    long long strippedMaterializationContextMicros = 0;
    long long strippedMaterializationAuditMicros = 0;
    long long strippedMaterializationMicros = 0;
    long long strippedPostLocalSnapshotMicros = 0;
    long long strippedPrimalMergeReleaseMicros = 0;
    long long strippedFinalReconstructionCheckMicros = 0;
    long long strippedMaterializationValidationMicros = 0;
    long long strippedMaterializationCreateNodesMicros = 0;
    long long strippedMaterializationAttachNodesMicros = 0;
    long long strippedMaterializationMoveProperPartsMicros = 0;
    long long strippedMaterializationChildBridgeMicros = 0;
    long long strippedMaterializationStructuralRepairMicros = 0;
    long long strippedMaterializationReleaseMicros = 0;
    long long strippedMaterializationRootPromotionMicros = 0;
    long long strippedMaterializationCanonicalizationMicros = 0;
    long long strippedMaterializationTargetCheckMicros = 0;
    std::size_t oracleLocalEvents = 0;
    std::size_t oracleLocalTargetMatches = 0;
    std::size_t oracleLocalStructureMatches = 0;
    std::size_t oracleLocalStructureMismatches = 0;
    long long oracleLocalCreatedNodes = 0;
    long long oracleLocalRecompositionNodes = 0;
    long long oracleLocalReleasedOwners = 0;
    long long totalComponentPixels = 0;
    int maxComponents = 0;
    int maxC = 0;
    int maxStrippedLocalPixels = 0;
    MaxCanonicalizationEventStats maxCanonicalizationEvent;
    MaxPrimalMergeEventStats maxPrimalMergeEvent;
};

struct PhaseReleaseStats {
    std::size_t comparisons = 0;
    std::size_t releaseEvaluated = 0;
    std::size_t releaseSatisfied = 0;
    std::size_t releaseRejected = 0;
    std::size_t releaseNoFocus = 0;
    std::size_t h1Failures = 0;
    std::size_t h2Failures = 0;
    std::size_t h3Failures = 0;
    std::size_t h4BoundaryFailures = 0;
    std::size_t h4IncidenceFailures = 0;
    std::size_t h5OrderFailures = 0;
    int firstOrderMismatchMin = -1;
};

void updateMaxCanonicalizationEvent(MergeEventStats &stats,
                                    const noninc::NonIncreasingMergeEvent &event) {
    if (!event.strippedLocalAttempted ||
        event.strippedMaterializationCanonicalizationMicros <= 0 ||
        event.strippedMaterializationCanonicalizationMicros <=
            stats.maxCanonicalizationEvent.micros) {
        return;
    }

    auto &maxEvent = stats.maxCanonicalizationEvent;
    maxEvent.present = true;
    maxEvent.micros = event.strippedMaterializationCanonicalizationMicros;
    maxEvent.threshold = event.threshold;
    maxEvent.phase = event.phase;
    maxEvent.nodeId = event.nodeId;
    maxEvent.parentId = event.parentId;
    maxEvent.nodeLevel = event.nodeLevel;
    maxEvent.parentLevel = event.parentLevel;
    maxEvent.numChildren = event.numChildren;
    maxEvent.componentCount = event.componentCount;
    maxEvent.totalComponentPixels = event.totalComponentPixels;
    maxEvent.strippedLocalPixelCount = event.strippedLocalPixelCount;
    maxEvent.strippedLocalCreatedNodes = event.strippedLocalCreatedNodes;
    maxEvent.strippedLocalRecompositionNodes = event.strippedLocalRecompositionNodes;
    maxEvent.strippedCreatedConnectorNodes = event.strippedCreatedConnectorNodes;
    maxEvent.strippedEmptyCreatedNodes =
        event.strippedEmptyCreatedNodesAfterCanonicalization;
    maxEvent.strippedEmptyConnectorNodes =
        event.strippedEmptyConnectorNodesAfterCanonicalization;
    maxEvent.materializationMicros = event.strippedMaterializationMicros;
    maxEvent.releaseMicros = event.strippedMaterializationReleaseMicros;
    maxEvent.targetCheckMicros = event.strippedMaterializationTargetCheckMicros;
}

void updateMaxPrimalMergeEvent(MergeEventStats &stats,
                               const noninc::NonIncreasingMergeEvent &event) {
    if (!event.strippedLocalAttempted ||
        event.strippedPrimalMergeReleaseMicros <= 0 ||
        event.strippedPrimalMergeReleaseMicros <= stats.maxPrimalMergeEvent.micros) {
        return;
    }

    auto &maxEvent = stats.maxPrimalMergeEvent;
    maxEvent.present = true;
    maxEvent.micros = event.strippedPrimalMergeReleaseMicros;
    maxEvent.threshold = event.threshold;
    maxEvent.phase = event.phase;
    maxEvent.nodeId = event.nodeId;
    maxEvent.parentId = event.parentId;
    maxEvent.nodeLevel = event.nodeLevel;
    maxEvent.parentLevel = event.parentLevel;
    maxEvent.directProperParts = event.directProperParts;
    maxEvent.numChildren = event.numChildren;
    maxEvent.componentCount = event.componentCount;
    maxEvent.totalComponentPixels = event.totalComponentPixels;
    maxEvent.strippedLocalPixelCount = event.strippedLocalPixelCount;
    maxEvent.strippedLocalCreatedNodes = event.strippedLocalCreatedNodes;
    maxEvent.strippedLocalRecompositionNodes = event.strippedLocalRecompositionNodes;
    maxEvent.strippedCreatedConnectorNodes = event.strippedCreatedConnectorNodes;
    maxEvent.materializationMicros = event.strippedMaterializationMicros;
    maxEvent.canonicalizationMicros =
        event.strippedMaterializationCanonicalizationMicros;
}

void printMaxCanonicalizationEventFields(const MergeEventStats &stats) {
    const auto &maxEvent = stats.maxCanonicalizationEvent;
    std::cout << " max_stripped_mat_canonicalize_ms="
              << common::microsToMillis(maxEvent.micros);
    if (!maxEvent.present) {
        return;
    }
    std::cout << " max_stripped_mat_canonicalize_threshold=" << maxEvent.threshold
              << " max_stripped_mat_canonicalize_phase=" << phaseName(maxEvent.phase)
              << " max_stripped_mat_canonicalize_node=" << maxEvent.nodeId
              << " max_stripped_mat_canonicalize_parent=" << maxEvent.parentId
              << " max_stripped_mat_canonicalize_node_level=" << maxEvent.nodeLevel
              << " max_stripped_mat_canonicalize_parent_level=" << maxEvent.parentLevel
              << " max_stripped_mat_canonicalize_children=" << maxEvent.numChildren
              << " max_stripped_mat_canonicalize_components="
              << maxEvent.componentCount
              << " max_stripped_mat_canonicalize_C="
              << maxEvent.totalComponentPixels
              << " max_stripped_mat_canonicalize_local_pixels="
              << maxEvent.strippedLocalPixelCount
              << " max_stripped_mat_canonicalize_created_nodes="
              << maxEvent.strippedLocalCreatedNodes
              << " max_stripped_mat_canonicalize_recomposition_nodes="
              << maxEvent.strippedLocalRecompositionNodes
              << " max_stripped_mat_canonicalize_connector_nodes="
              << maxEvent.strippedCreatedConnectorNodes
              << " max_stripped_mat_canonicalize_empty_created_nodes="
              << maxEvent.strippedEmptyCreatedNodes
              << " max_stripped_mat_canonicalize_empty_connector_nodes="
              << maxEvent.strippedEmptyConnectorNodes
              << " max_stripped_mat_canonicalize_materialization_ms="
              << common::microsToMillis(maxEvent.materializationMicros)
              << " max_stripped_mat_canonicalize_release_ms="
              << common::microsToMillis(maxEvent.releaseMicros)
              << " max_stripped_mat_canonicalize_target_check_ms="
              << common::microsToMillis(maxEvent.targetCheckMicros);
}

void printMaxPrimalMergeEventFields(const MergeEventStats &stats) {
    const auto &maxEvent = stats.maxPrimalMergeEvent;
    std::cout << " max_stripped_primal_merge_ms="
              << common::microsToMillis(maxEvent.micros);
    if (!maxEvent.present) {
        return;
    }
    std::cout << " max_stripped_primal_merge_threshold=" << maxEvent.threshold
              << " max_stripped_primal_merge_phase=" << phaseName(maxEvent.phase)
              << " max_stripped_primal_merge_node=" << maxEvent.nodeId
              << " max_stripped_primal_merge_parent=" << maxEvent.parentId
              << " max_stripped_primal_merge_node_level=" << maxEvent.nodeLevel
              << " max_stripped_primal_merge_parent_level=" << maxEvent.parentLevel
              << " max_stripped_primal_merge_direct_pp="
              << maxEvent.directProperParts
              << " max_stripped_primal_merge_children=" << maxEvent.numChildren
              << " max_stripped_primal_merge_components="
              << maxEvent.componentCount
              << " max_stripped_primal_merge_C="
              << maxEvent.totalComponentPixels
              << " max_stripped_primal_merge_local_pixels="
              << maxEvent.strippedLocalPixelCount
              << " max_stripped_primal_merge_created_nodes="
              << maxEvent.strippedLocalCreatedNodes
              << " max_stripped_primal_merge_recomposition_nodes="
              << maxEvent.strippedLocalRecompositionNodes
              << " max_stripped_primal_merge_connector_nodes="
              << maxEvent.strippedCreatedConnectorNodes
              << " max_stripped_primal_merge_materialization_ms="
              << common::microsToMillis(maxEvent.materializationMicros)
              << " max_stripped_primal_merge_canonicalize_ms="
              << common::microsToMillis(maxEvent.canonicalizationMicros);
}

void addPartialFRejectReason(MergeEventStats &stats,
                             const std::string &reasonCode) {
    if (reasonCode == "invalid_event") {
        ++stats.partialFRejectInvalidEvent;
    } else if (reasonCode == "empty_C") {
        ++stats.partialFRejectEmptyC;
    } else if (reasonCode == "non_connected_C") {
        ++stats.partialFRejectNonConnectedC;
    } else if (reasonCode == "invalid_pixel") {
        ++stats.partialFRejectInvalidPixel;
    } else if (reasonCode == "duplicated_pixels") {
        ++stats.partialFRejectDuplicatedPixels;
    } else if (reasonCode == "dead_owner_pixel") {
        ++stats.partialFRejectDeadOwnerPixel;
    } else if (reasonCode == "invalid_target_level") {
        ++stats.partialFRejectInvalidTargetLevel;
    } else if (reasonCode == "public_problem1_preflight_rejected" ||
               reasonCode == "published_problem_1_contract_violation") {
        ++stats.partialFRejectPublicProblem1Preflight;
    } else if (reasonCode == "mutation_safety_unproven") {
        ++stats.partialFRejectMutationSafetyUnproven;
    } else if (reasonCode == "d12_base_preflight_failed") {
        ++stats.partialFRejectD12BasePreflightFailed;
    } else if (reasonCode == "d12_plan_incomplete") {
        ++stats.partialFRejectD12PlanIncomplete;
    } else if (reasonCode == "d12_plan_oracle_not_free") {
        ++stats.partialFRejectD12PlanOracleNotFree;
    } else if (reasonCode == "d12_materialization_context_missing") {
        ++stats.partialFRejectD12MaterializationContextMissing;
    } else if (reasonCode == "d12_materialization_oracle_not_free") {
        ++stats.partialFRejectD12MaterializationOracleNotFree;
    } else if (reasonCode == "d12_materialization_inputs_invalid") {
        ++stats.partialFRejectD12MaterializationInputsInvalid;
    } else if (reasonCode == "d12_rollback_summary_missing") {
        ++stats.partialFRejectD12RollbackSummaryMissing;
    } else if (reasonCode == "d12_invalid_stable_child_ref") {
        ++stats.partialFRejectD12InvalidStableChildRef;
    } else if (reasonCode == "d12_stable_child_cycle_risk") {
        ++stats.partialFRejectD12StableChildCycleRisk;
    } else if (reasonCode == "d12_direct_v1_precheck_missing") {
        ++stats.partialFRejectD12DirectV1PrecheckMissing;
    } else if (reasonCode == "d12_hidden_support_not_redundant") {
        ++stats.partialFRejectD12HiddenSupportNotRedundant;
    } else if (reasonCode == "d12_parent_link_not_represented") {
        ++stats.partialFRejectD12ParentLinkNotRepresented;
    } else if (reasonCode == "d12_parent_link_invalid_level") {
        ++stats.partialFRejectD12ParentLinkInvalidLevel;
    } else if (reasonCode == "d12_direct_v1_precheck_rejected") {
        ++stats.partialFRejectD12DirectV1PrecheckRejected;
    } else if (reasonCode == "d12_exception_free_unproven") {
        ++stats.partialFRejectD12ExceptionFreeUnproven;
    } else if (reasonCode == "d12_dry_run_empty_plan") {
        ++stats.partialFRejectD12DryRunEmptyPlan;
    } else if (reasonCode == "d12_dry_run_auxiliary_nodes_present") {
        ++stats.partialFRejectD12DryRunAuxiliaryNodesPresent;
    } else if (reasonCode == "d12_dry_run_stable_child_bridge_present") {
        ++stats.partialFRejectD12DryRunStableChildBridgePresent;
    } else if (reasonCode == "d12_dry_run_non_oracle_hidden_support_connector_present") {
        ++stats.partialFRejectD12DryRunNonOracleHiddenSupportConnectorPresent;
    } else if (reasonCode == "d12_dry_run_parent_link_bridge_present") {
        ++stats.partialFRejectD12DryRunParentLinkBridgePresent;
    } else if (reasonCode == "d12_dry_run_child_link_bridge_present") {
        ++stats.partialFRejectD12DryRunChildLinkBridgePresent;
    } else if (reasonCode == "d12_dry_run_external_attachment_invalid") {
        ++stats.partialFRejectD12DryRunExternalAttachmentInvalid;
    } else if (reasonCode == "d12_dry_run_external_attachment_empty") {
        ++stats.partialFRejectD12DryRunExternalAttachmentEmpty;
    } else if (reasonCode == "d12_dry_run_stable_child_move_present") {
        ++stats.partialFRejectD12DryRunStableChildMovePresent;
    } else if (reasonCode == "d12_dry_run_empty_owner_present") {
        ++stats.partialFRejectD12DryRunEmptyOwnerPresent;
    } else if (reasonCode == "d12_dry_run_empty_owner_invalid_direct_proper_part") {
        ++stats.partialFRejectD12DryRunEmptyOwnerInvalidDirectProperPart;
    } else if (reasonCode == "d12_dry_run_empty_owner_non_live_source_owner") {
        ++stats.partialFRejectD12DryRunEmptyOwnerNonLiveSourceOwner;
    } else if (reasonCode == "d12_dry_run_empty_owner_invalid_stable_child") {
        ++stats.partialFRejectD12DryRunEmptyOwnerInvalidStableChild;
    } else if (reasonCode == "d12_dry_run_empty_owner_root_source_owner") {
        ++stats.partialFRejectD12DryRunEmptyOwnerRootSourceOwner;
    } else if (reasonCode == "d12_dry_run_empty_owner_source_owner_not_leaf") {
        ++stats.partialFRejectD12DryRunEmptyOwnerSourceOwnerNotLeaf;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_invalid_parent") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidParent;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_invalid_child") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidChild;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_invalid_child_level") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidChildLevel;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_cycle_risk") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenCycleRisk;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_duplicate_child") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenDuplicateChild;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_invalid_planned_level") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidPlannedLevel;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_residual_child_present") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildPresent;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_residual_child_snapshot_missing") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildSnapshotMissing;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_residual_child_touches_modified") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildTouchesModified;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_residual_child_support_overlaps_plan") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildSupportOverlapsPlan;
    } else if (reasonCode == "d12_dry_run_empty_owner_with_children_external_dependency") {
        ++stats.partialFRejectD12DryRunEmptyOwnerWithChildrenExternalDependency;
    } else if (reasonCode == "d12_dry_run_empty_owner_source_count_mismatch") {
        ++stats.partialFRejectD12DryRunEmptyOwnerSourceCountMismatch;
    } else if (reasonCode == "d12_dry_run_empty_owner_root_previous_parent") {
        ++stats.partialFRejectD12DryRunEmptyOwnerRootPreviousParent;
    } else if (reasonCode == "d12_dry_run_empty_owner_previous_parent_count_mismatch") {
        ++stats.partialFRejectD12DryRunEmptyOwnerPreviousParentCountMismatch;
    } else if (reasonCode == "d12_dry_run_invalid_parent_plan") {
        ++stats.partialFRejectD12DryRunInvalidParentPlan;
    } else if (reasonCode == "d12_dry_run_invalid_parent_level") {
        ++stats.partialFRejectD12DryRunInvalidParentLevel;
    } else if (reasonCode == "d12_dry_run_empty_planned_node") {
        ++stats.partialFRejectD12DryRunEmptyPlannedNode;
    } else if (reasonCode == "d12_dry_run_child_move_present") {
        ++stats.partialFRejectD12DryRunChildMovePresent;
    } else if (reasonCode == "d12_dry_run_invalid_pixel") {
        ++stats.partialFRejectD12DryRunInvalidPixel;
    } else if (reasonCode == "d12_dry_run_duplicated_pixel") {
        ++stats.partialFRejectD12DryRunDuplicatedPixel;
    } else if (reasonCode == "current_updater_limitation") {
        ++stats.partialFRejectCurrentUpdaterLimitation;
    } else {
        ++stats.partialFRejectOther;
    }
}

std::string partialFRejectStatsText(const MergeEventStats &stats) {
    std::ostringstream out;
    out << " partial_f_reject_invalid_event=" << stats.partialFRejectInvalidEvent
        << " partial_f_reject_empty_c=" << stats.partialFRejectEmptyC
        << " partial_f_reject_non_connected_c=" << stats.partialFRejectNonConnectedC
        << " partial_f_reject_invalid_pixel=" << stats.partialFRejectInvalidPixel
        << " partial_f_reject_duplicated_pixels=" << stats.partialFRejectDuplicatedPixels
        << " partial_f_reject_dead_owner_pixel=" << stats.partialFRejectDeadOwnerPixel
        << " partial_f_reject_invalid_target_level=" << stats.partialFRejectInvalidTargetLevel
        << " partial_f_reject_public_problem1_preflight="
        << stats.partialFRejectPublicProblem1Preflight
        << " partial_f_reject_mutation_safety_unproven="
        << stats.partialFRejectMutationSafetyUnproven
        << " partial_f_reject_d12_base_preflight_failed="
        << stats.partialFRejectD12BasePreflightFailed
        << " partial_f_reject_d12_plan_incomplete="
        << stats.partialFRejectD12PlanIncomplete
        << " partial_f_reject_d12_plan_oracle_not_free="
        << stats.partialFRejectD12PlanOracleNotFree
        << " partial_f_reject_d12_materialization_context_missing="
        << stats.partialFRejectD12MaterializationContextMissing
        << " partial_f_reject_d12_materialization_oracle_not_free="
        << stats.partialFRejectD12MaterializationOracleNotFree
        << " partial_f_reject_d12_materialization_inputs_invalid="
        << stats.partialFRejectD12MaterializationInputsInvalid
        << " partial_f_reject_d12_rollback_summary_missing="
        << stats.partialFRejectD12RollbackSummaryMissing
        << " partial_f_reject_d12_invalid_stable_child_ref="
        << stats.partialFRejectD12InvalidStableChildRef
        << " partial_f_reject_d12_stable_child_cycle_risk="
        << stats.partialFRejectD12StableChildCycleRisk
        << " partial_f_reject_d12_direct_v1_precheck_missing="
        << stats.partialFRejectD12DirectV1PrecheckMissing
        << " partial_f_reject_d12_hidden_support_not_redundant="
        << stats.partialFRejectD12HiddenSupportNotRedundant
        << " partial_f_reject_d12_parent_link_not_represented="
        << stats.partialFRejectD12ParentLinkNotRepresented
        << " partial_f_reject_d12_parent_link_invalid_level="
        << stats.partialFRejectD12ParentLinkInvalidLevel
        << " partial_f_reject_d12_direct_v1_precheck_rejected="
        << stats.partialFRejectD12DirectV1PrecheckRejected
        << " partial_f_reject_d12_exception_free_unproven="
        << stats.partialFRejectD12ExceptionFreeUnproven
        << " partial_f_reject_d12_dry_run_empty_plan="
        << stats.partialFRejectD12DryRunEmptyPlan
        << " partial_f_reject_d12_dry_run_auxiliary_nodes_present="
        << stats.partialFRejectD12DryRunAuxiliaryNodesPresent
        << " partial_f_reject_d12_dry_run_stable_child_bridge_present="
        << stats.partialFRejectD12DryRunStableChildBridgePresent
        << " partial_f_reject_d12_dry_run_non_oracle_hidden_support_connector_present="
        << stats.partialFRejectD12DryRunNonOracleHiddenSupportConnectorPresent
        << " partial_f_reject_d12_dry_run_parent_link_bridge_present="
        << stats.partialFRejectD12DryRunParentLinkBridgePresent
        << " partial_f_reject_d12_dry_run_child_link_bridge_present="
        << stats.partialFRejectD12DryRunChildLinkBridgePresent
        << " partial_f_reject_d12_dry_run_external_attachment_invalid="
        << stats.partialFRejectD12DryRunExternalAttachmentInvalid
        << " partial_f_reject_d12_dry_run_external_attachment_empty="
        << stats.partialFRejectD12DryRunExternalAttachmentEmpty
        << " partial_f_reject_d12_dry_run_stable_child_move_present="
        << stats.partialFRejectD12DryRunStableChildMovePresent
        << " partial_f_reject_d12_dry_run_empty_owner_present="
        << stats.partialFRejectD12DryRunEmptyOwnerPresent
        << " partial_f_reject_d12_dry_run_empty_owner_invalid_direct_proper_part="
        << stats.partialFRejectD12DryRunEmptyOwnerInvalidDirectProperPart
        << " partial_f_reject_d12_dry_run_empty_owner_non_live_source_owner="
        << stats.partialFRejectD12DryRunEmptyOwnerNonLiveSourceOwner
        << " partial_f_reject_d12_dry_run_empty_owner_invalid_stable_child="
        << stats.partialFRejectD12DryRunEmptyOwnerInvalidStableChild
        << " partial_f_reject_d12_dry_run_empty_owner_root_source_owner="
        << stats.partialFRejectD12DryRunEmptyOwnerRootSourceOwner
        << " partial_f_reject_d12_dry_run_empty_owner_source_owner_not_leaf="
        << stats.partialFRejectD12DryRunEmptyOwnerSourceOwnerNotLeaf
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_invalid_parent="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidParent
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_invalid_child="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidChild
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_invalid_child_level="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidChildLevel
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_cycle_risk="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenCycleRisk
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_duplicate_child="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenDuplicateChild
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_invalid_planned_level="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenInvalidPlannedLevel
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_residual_child_present="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildPresent
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_residual_child_snapshot_missing="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildSnapshotMissing
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_residual_child_touches_modified="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildTouchesModified
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_residual_child_support_overlaps_plan="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenResidualChildSupportOverlapsPlan
        << " partial_f_reject_d12_dry_run_empty_owner_with_children_external_dependency="
        << stats.partialFRejectD12DryRunEmptyOwnerWithChildrenExternalDependency
        << " partial_f_reject_d12_dry_run_empty_owner_source_count_mismatch="
        << stats.partialFRejectD12DryRunEmptyOwnerSourceCountMismatch
        << " partial_f_reject_d12_dry_run_empty_owner_root_previous_parent="
        << stats.partialFRejectD12DryRunEmptyOwnerRootPreviousParent
        << " partial_f_reject_d12_dry_run_empty_owner_previous_parent_count_mismatch="
        << stats.partialFRejectD12DryRunEmptyOwnerPreviousParentCountMismatch
        << " partial_f_reject_d12_dry_run_invalid_parent_plan="
        << stats.partialFRejectD12DryRunInvalidParentPlan
        << " partial_f_reject_d12_dry_run_invalid_parent_level="
        << stats.partialFRejectD12DryRunInvalidParentLevel
        << " partial_f_reject_d12_dry_run_empty_planned_node="
        << stats.partialFRejectD12DryRunEmptyPlannedNode
        << " partial_f_reject_d12_dry_run_child_move_present="
        << stats.partialFRejectD12DryRunChildMovePresent
        << " partial_f_reject_d12_dry_run_invalid_pixel="
        << stats.partialFRejectD12DryRunInvalidPixel
        << " partial_f_reject_d12_dry_run_duplicated_pixel="
        << stats.partialFRejectD12DryRunDuplicatedPixel
        << " partial_f_reject_current_updater_limitation="
        << stats.partialFRejectCurrentUpdaterLimitation
        << " partial_f_reject_other=" << stats.partialFRejectOther;
    return out.str();
}

void addEventStats(MergeEventStats &stats,
                   const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (const auto &event : events) {
        ++stats.events;
        if (event.phase == noninc::NonIncreasingPhase::MaxToMin) {
            ++stats.maxToMinEvents;
        } else {
            ++stats.minToMaxEvents;
        }
        if (event.componentCount > 1) {
            ++stats.multiComponentEvents;
        }
        if (event.numChildren > 0) {
            ++stats.childEvents;
        }
        if (event.dualRebuiltAfterLocalMismatch) {
            ++stats.dualRebuildEvents;
        }
        if (event.strippedLocalAttempted) {
            ++stats.strippedLocalAttempts;
            if (event.strippedLocalAccepted) {
                ++stats.strippedLocalAccepted;
            }
            if (event.strippedLocalApplied) {
                ++stats.strippedLocalApplied;
            }
            if (event.strippedDirectCommit) {
                ++stats.strippedDirectCommits;
            }
            if (event.strippedLocalPlanAuditChecked) {
                ++stats.strippedLocalPlanAuditChecks;
            }
            if (event.strippedLocalMatchesTarget) {
                ++stats.strippedLocalTargetMatches;
            }
            if (event.strippedLocalTargetChecked) {
                ++stats.strippedLocalTargetChecks;
            }
            if (event.strippedLocalPlanOracleFree) {
                ++stats.strippedLocalOracleFreePlans;
            }
            if (event.strippedMaterializationAuditChecked) {
                ++stats.strippedMaterializationAuditChecks;
            }
            if (event.strippedMaterializationInputOracleFree) {
                ++stats.strippedLocalOracleFreeMaterializations;
            }
            if (event.strippedFinalReconstructionChecked) {
                ++stats.strippedFinalReconstructionChecks;
            }
            if (event.strippedLocalFallbackToOracle) {
                ++stats.strippedLocalFallbacks;
            }
            if (event.strippedLocalFallbackToRebuild) {
                ++stats.strippedLocalRebuildFallbacks;
            }
            stats.strippedLocalCreatedNodes += event.strippedLocalCreatedNodes;
            stats.strippedLocalRecompositionNodes += event.strippedLocalRecompositionNodes;
            stats.strippedLocalPixels += event.strippedLocalPixelCount;
            stats.maxStrippedLocalPixels =
                std::max(stats.maxStrippedLocalPixels, event.strippedLocalPixelCount);
            stats.strippedLocalReleasedOwners += event.strippedLocalReleasedSourceOwners;
            stats.strippedEmptyCreatedNodes +=
                event.strippedEmptyCreatedNodesAfterCanonicalization;
            stats.strippedEmptyConnectorNodes +=
                event.strippedEmptyConnectorNodesAfterCanonicalization;
            stats.strippedLocalCopyMicros += event.strippedLocalCopyMicros;
            stats.strippedLocalUpdateMicros += event.strippedLocalUpdateMicros;
            stats.strippedOracleFallbackMicros += event.strippedOracleFallbackMicros;
            stats.strippedRebuildFallbackMicros += event.strippedRebuildFallbackMicros;
            stats.strippedInputExtractionMicros += event.strippedInputExtractionMicros;
            stats.strippedPartitionPlanMicros += event.strippedPartitionPlanMicros;
            stats.strippedRecompositionPlanMicros += event.strippedRecompositionPlanMicros;
            stats.strippedPlanAuditMicros += event.strippedPlanAuditMicros;
            stats.strippedMaterializationContextMicros += event.strippedMaterializationContextMicros;
            stats.strippedMaterializationAuditMicros += event.strippedMaterializationAuditMicros;
            stats.strippedMaterializationMicros += event.strippedMaterializationMicros;
            stats.strippedPostLocalSnapshotMicros += event.strippedPostLocalSnapshotMicros;
            stats.strippedPrimalMergeReleaseMicros += event.strippedPrimalMergeReleaseMicros;
            stats.strippedFinalReconstructionCheckMicros +=
                event.strippedFinalReconstructionCheckMicros;
            updateMaxPrimalMergeEvent(stats, event);
            stats.strippedMaterializationValidationMicros +=
                event.strippedMaterializationValidationMicros;
            stats.strippedMaterializationCreateNodesMicros +=
                event.strippedMaterializationCreateNodesMicros;
            stats.strippedMaterializationAttachNodesMicros +=
                event.strippedMaterializationAttachNodesMicros;
            stats.strippedMaterializationMoveProperPartsMicros +=
                event.strippedMaterializationMoveProperPartsMicros;
            stats.strippedMaterializationChildBridgeMicros +=
                event.strippedMaterializationChildBridgeMicros;
            stats.strippedMaterializationStructuralRepairMicros +=
                event.strippedMaterializationStructuralRepairMicros;
            stats.strippedMaterializationReleaseMicros +=
                event.strippedMaterializationReleaseMicros;
            stats.strippedMaterializationRootPromotionMicros +=
                event.strippedMaterializationRootPromotionMicros;
            stats.strippedMaterializationCanonicalizationMicros +=
                event.strippedMaterializationCanonicalizationMicros;
            stats.strippedMaterializationTargetCheckMicros +=
                event.strippedMaterializationTargetCheckMicros;
            updateMaxCanonicalizationEvent(stats, event);
        }
        if (event.strippedLocalFallbackToPhaseRebuild) {
            ++stats.strippedLocalPhaseRebuildFallbacks;
        }
        if (event.strippedLocalSkippedByPhaseRebuild) {
            ++stats.strippedLocalPhaseRebuildSkips;
        }
        if (event.mergeNodeLocalAttempted) {
            ++stats.mergeNodeLocalAttempts;
            if (event.partialFAttempted) {
                ++stats.partialFAttempts;
            }
            if (event.partialFAccepted) {
                ++stats.partialFAccepted;
            }
            if (event.partialFRejected) {
                ++stats.partialFRejected;
                addPartialFRejectReason(stats, event.partialFRejectReasonCode);
            }
            if (event.partialFMutationAttempted) {
                ++stats.partialFMutationAttempted;
            }
            if (event.partialFMcpcLiveUpdate) {
                ++stats.partialFMcpcLiveAccepted;
            }
            if (event.partialFMcpcHiddenSupportRedundantLiveUpdate) {
                ++stats.partialFMcpcHiddenSupportRedundantLiveAccepted;
            }
            if (event.partialFMcpcParentLinkBridgeLiveUpdate) {
                ++stats.partialFMcpcParentLinkBridgeLiveAccepted;
            }
            if (event.partialFMcpcStableChildBridgeLiveUpdate) {
                ++stats.partialFMcpcStableChildBridgeLiveAccepted;
            }
            if (event.partialFMcpcChildLinkBridgeLiveUpdate) {
                ++stats.partialFMcpcChildLinkBridgeLiveAccepted;
            }
            if (event.partialFMcpcEmptyOwnerReleaseLiveUpdate) {
                ++stats.partialFMcpcEmptyOwnerReleaseLiveAccepted;
            }
            if (event.partialFMcpcStableChildMoveLiveUpdate) {
                ++stats.partialFMcpcStableChildMoveLiveAccepted;
            }
            if (event.mergeNodeLocalAccepted) {
                ++stats.mergeNodeLocalAccepted;
            }
            if (event.mergeNodeLocalApplied) {
                ++stats.mergeNodeLocalApplied;
            }
            if (event.mergeNodeLocalCertificationChecked) {
                ++stats.mergeNodeLocalCertificationChecks;
            }
            if (event.mergeNodeLocalCertificationPassed) {
                ++stats.mergeNodeLocalCertificationPassed;
            }
            if (event.mergeNodeLocalDirectCommit) {
                ++stats.mergeNodeLocalDirectCommits;
            }
            if (event.mergeNodeLocalTargetCreationAllowed) {
                ++stats.mergeNodeLocalTargetCreationAllowedEvents;
            }
            if (event.mergeNodeLocalTargetCreationRequired) {
                ++stats.mergeNodeLocalTargetCreationRequiredEvents;
            }
            if (event.mergeNodeLocalDualMatchesPrimalBeforeMerge) {
                ++stats.mergeNodeLocalPreDualMatches;
            }
            if (event.mergeNodeLocalDualMatchesPrimalAfterMerge) {
                ++stats.mergeNodeLocalDualMatches;
            }
            stats.mergeNodeLocalRejectedComponents +=
                event.mergeNodeLocalRejectedComponents;
            stats.mergeNodeLocalProblem1PreflightRejectedComponents +=
                event.mergeNodeLocalProblem1PreflightRejectedComponents;
            stats.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem +=
                event.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem;
            stats.mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation +=
                event.mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation;
            stats.mergeNodeLocalPartialOwnerComponents +=
                event.mergeNodeLocalPartialOwnerComponents;
            stats.mergeNodeLocalRawMissingTargetComponents +=
                event.mergeNodeLocalRawMissingTargetComponents;
            stats.mergeNodeLocalMissingTargetComponents +=
                event.mergeNodeLocalMissingTargetComponents;
            stats.mergeNodeLocalPartialOwnerAndMissingTargetComponents +=
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents;
            stats.mergeNodeLocalPixels += event.mergeNodeLocalPixelCount;
            stats.mergeNodeLocalCopyMicros += event.mergeNodeLocalCopyMicros;
            stats.mergeNodeLocalCertificationMicros +=
                event.mergeNodeLocalCertificationMicros;
            stats.mergeNodeComponentCollectMicros +=
                event.mergeNodeComponentCollectMicros;
            stats.mergeNodeLocalUpdateMicros += event.mergeNodeLocalUpdateMicros;
            stats.mergeNodePrimalMergeMicros += event.mergeNodePrimalMergeMicros;
            stats.mergeNodeEmptyReleaseMicros += event.mergeNodeEmptyReleaseMicros;
        }
        if (event.mergeNodeLocalFallbackToPhaseRebuild) {
            ++stats.mergeNodeLocalPhaseRebuildFallbacks;
        }
        if (event.mergeNodeLocalSkippedByPhaseRebuild) {
            ++stats.mergeNodeLocalPhaseRebuildSkips;
        }
        if (event.oracleGuidedLocalApplied) {
            ++stats.oracleLocalEvents;
            if (event.oracleGuidedLocalMatchesTarget) {
                ++stats.oracleLocalTargetMatches;
            }
            if (event.oracleGuidedLocalStructureMatchesOracle) {
                ++stats.oracleLocalStructureMatches;
            } else {
                ++stats.oracleLocalStructureMismatches;
            }
            stats.oracleLocalCreatedNodes += event.oracleGuidedCreatedNodes;
            stats.oracleLocalRecompositionNodes += event.oracleGuidedRecompositionNodes;
            stats.oracleLocalReleasedOwners += event.oracleGuidedReleasedSourceOwners;
        }
        stats.totalComponentPixels += event.totalComponentPixels;
        stats.maxComponents = std::max(stats.maxComponents, event.componentCount);
        stats.maxC = std::max(stats.maxC, event.totalComponentPixels);
    }
}

std::vector<float> parseThresholds(int argc, char **argv, int argi) {
    std::vector<float> thresholds;
    thresholds.reserve(static_cast<std::size_t>(argc - argi));
    for (int i = argi; i < argc; ++i) {
        char *endPtr = nullptr;
        const float value = std::strtof(argv[i], &endPtr);
        if (endPtr == argv[i] || *endPtr != '\0') {
            throw std::runtime_error(std::string("Invalid threshold: ") + argv[i]);
        }
        thresholds.push_back(value);
    }
    return thresholds;
}

void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0
              << " [--mode explicit|naive|compare]"
              << " [--attribute gray-variance|box-area-perimeter-ratio]"
              << " [--dual-maintenance partial-f|problem1-restricted|p1-restricted|problem1-close-recompose|explicit-only|materialize-explicit|oracle-local|stripped-oracle-local|stripped-rebuild-local|stripped-phase-rebuild-local|merge-node-phase-rebuild|merge-node-certified-phase-rebuild|merge-node-certified-create-b-phase-rebuild|merge-node-certified-create-b-mcpc-phase-rebuild|merge-node-clean-unchecked|phase-rebuild]"
              << " [--stripped-validation safe|no-reconstruction|unchecked]"
              << " [--stripped-commit safe-copy|direct]"
              << " [--strict-direct-v1-precheck]"
              << " [--mcpc-no-copy-direct-commit]"
              << " [--mcpc-validated-direct-commit]"
              << " [--diagnose-first-divergence]"
              << " [--stop-on-first-phase-fallback]"
              << " [--find-divergent-prefix]"
              << " [--dump-prefix-fixture <n>]"
              << " [--prefix-fixture-padding <n>]"
              << " [--diagnose-events]"
              << " [--diagnose-phase-oracle]"
              << " [--phase-candidate-guard]"
              << " [--dump-phase-oracle-fixture]"
              << " [--dump-smallest-phase-oracle-fixture]"
              << " [--dump-incident-closure-fixture]"
              << " [--dump-first-local-oracle]"
              << " [--dump-local-oracle <threshold> <max->min|min->max> <nodeId>]"
              << " [--max-merge-events <n>]"
              << " [--time-budget-ms <n>]"
              << " [--disable-leaf-fast-path]"
              << " [--radio-adj <radius>]"
              << " [--no-output]"
              << " <input.png> [<output.png>] <threshold1> [threshold2 ...]\n";
}

noninc::NonIncreasingAttributeKind parseAttributeKind(const std::string &value) {
    if (value == "gray-variance" ||
        value == "gray_level_variance" ||
        value == "gray-variance-attribute") {
        return noninc::NonIncreasingAttributeKind::GrayVariance;
    }
    if (value == "box-area-perimeter-ratio" ||
        value == "bbox-area-perimeter-ratio" ||
        value == "area-per-box-perimeter") {
        return noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio;
    }
    throw std::runtime_error("Invalid attribute: " + value);
}

std::string attributeKindName(noninc::NonIncreasingAttributeKind attributeKind) {
    switch (attributeKind) {
        case noninc::NonIncreasingAttributeKind::GrayVariance:
            return "gray_level_variance";
        case noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio:
            return "box_area_perimeter_ratio";
    }
    return "unknown";
}

void printPhaseOracleComparison(const noninc::NonIncreasingPhaseOracleComparison &comparison) {
    std::cout << "  phase_oracle"
              << " threshold=" << comparison.threshold
              << " point="
              << (comparison.beforeMaxToMin ? "before_max->min" : "after_phase")
              << " phase=" << phaseName(comparison.completedPhase)
              << " reconstruction_matches="
              << (comparison.localReconstructionMatchesRebuild ? "true" : "false")
              << " candidate_set_matches="
              << (comparison.candidateSetMatchesRebuild ? "true" : "false")
              << " candidate_order_matches="
              << (comparison.candidateScanOrderMatchesRebuild ? "true" : "false")
              << " local_candidates=" << comparison.localCandidateCount
              << " rebuild_candidates=" << comparison.rebuildCandidateCount
              << " local_only=" << comparison.localOnlyCandidateCount
              << " rebuild_only=" << comparison.rebuildOnlyCandidateCount
              << " first_candidate_order_mismatch="
              << comparison.firstCandidateScanOrderMismatchIndex
              << " envelope_structure_compared="
              << (comparison.envelopeStructureCompared ? "true" : "false")
              << " envelope_structure_matches="
              << (comparison.envelopeStructureMatchesRebuild ? "true" : "false")
              << " local_envelope_nodes=" << comparison.localEnvelopeNodeCount
              << " rebuild_envelope_nodes=" << comparison.rebuildEnvelopeNodeCount
              << " boundary_structure_compared="
              << (comparison.boundaryStructureCompared ? "true" : "false")
              << " boundary_structure_matches="
              << (comparison.boundaryStructureMatchesRebuild ? "true" : "false")
              << " boundary_incidence_matches="
              << (comparison.boundaryIncidenceMatchesRebuild ? "true" : "false")
              << " local_boundary_nodes=" << comparison.localBoundaryNodeCount
              << " rebuild_boundary_nodes=" << comparison.rebuildBoundaryNodeCount
              << " local_boundary_incidences=" << comparison.localBoundaryIncidenceCount
              << " rebuild_boundary_incidences=" << comparison.rebuildBoundaryIncidenceCount
              << " release_evaluated="
              << (comparison.releasePredicateEvaluated ? "true" : "false")
              << " release_satisfied="
              << (comparison.releasePredicateSatisfied ? "true" : "false")
              << " release_focus_available="
              << (comparison.releaseFocusAvailable ? "true" : "false")
              << " release_h1="
              << (comparison.releaseH1ReconstructionMatches ? "true" : "false")
              << " release_h2="
              << (comparison.releaseH2CandidateSetMatches ? "true" : "false")
              << " release_h3="
              << (comparison.releaseH3EnvelopeStructureMatches ? "true" : "false")
              << " release_h4_boundary="
              << (comparison.releaseH4BoundaryStructureMatches ? "true" : "false")
              << " release_h4_incidence="
              << (comparison.releaseH4BoundaryIncidenceMatches ? "true" : "false")
              << " release_h5_order="
              << (comparison.releaseH5CandidateScanOrderMatches ? "true" : "false")
              << " release_focus_nodes=" << comparison.releaseFocusNodeCount
              << " release_focus_pixels=" << comparison.releaseFocusPixelCount
              << " release_closure_nodes=" << comparison.releaseClosureNodeCount
              << " release_closure_pixels=" << comparison.releaseClosurePixelCount
              << " first_local_only_level=" << comparison.firstLocalOnlyAltitude
              << " first_local_only_support=" << comparison.firstLocalOnlySupportSize
              << " first_rebuild_only_level=" << comparison.firstRebuildOnlyAltitude
              << " first_rebuild_only_support=" << comparison.firstRebuildOnlySupportSize;
    if (comparison.firstLocalOnlyAltitude >= 0) {
        printCandidateContext("local_only", comparison.firstLocalOnlyContext);
    }
    if (comparison.firstRebuildOnlyAltitude >= 0) {
        printCandidateContext("rebuild_only", comparison.firstRebuildOnlyContext);
    }
    std::cout << "\n";
}

void printPhaseOracleComparisons(
    const std::vector<noninc::NonIncreasingPhaseOracleComparison> &comparisons) {
    PhaseReleaseStats stats;
    for (const auto &comparison : comparisons) {
        printPhaseOracleComparison(comparison);
        ++stats.comparisons;
        if (!comparison.releasePredicateEvaluated) {
            ++stats.releaseNoFocus;
        } else {
            ++stats.releaseEvaluated;
            if (comparison.releasePredicateSatisfied) {
                ++stats.releaseSatisfied;
            } else {
                ++stats.releaseRejected;
            }
        }
        if (!comparison.releaseH1ReconstructionMatches) {
            ++stats.h1Failures;
        }
        if (!comparison.releaseH2CandidateSetMatches) {
            ++stats.h2Failures;
        }
        if (comparison.releasePredicateEvaluated &&
            !comparison.releaseH3EnvelopeStructureMatches) {
            ++stats.h3Failures;
        }
        if (comparison.releasePredicateEvaluated &&
            !comparison.releaseH4BoundaryStructureMatches) {
            ++stats.h4BoundaryFailures;
        }
        if (comparison.releasePredicateEvaluated &&
            !comparison.releaseH4BoundaryIncidenceMatches) {
            ++stats.h4IncidenceFailures;
        }
        if (!comparison.releaseH5CandidateScanOrderMatches) {
            ++stats.h5OrderFailures;
        }
        if (comparison.firstCandidateScanOrderMismatchIndex >= 0 &&
            (stats.firstOrderMismatchMin < 0 ||
             comparison.firstCandidateScanOrderMismatchIndex < stats.firstOrderMismatchMin)) {
            stats.firstOrderMismatchMin =
                comparison.firstCandidateScanOrderMismatchIndex;
        }
    }
    std::cout << "  phase_release_summary"
              << " comparisons=" << stats.comparisons
              << " release_evaluated=" << stats.releaseEvaluated
              << " release_satisfied=" << stats.releaseSatisfied
              << " release_rejected=" << stats.releaseRejected
              << " release_no_focus=" << stats.releaseNoFocus
              << " h1_failures=" << stats.h1Failures
              << " h2_failures=" << stats.h2Failures
              << " h3_failures=" << stats.h3Failures
              << " h4_boundary_failures=" << stats.h4BoundaryFailures
              << " h4_incidence_failures=" << stats.h4IncidenceFailures
              << " h5_order_failures=" << stats.h5OrderFailures
              << " first_order_mismatch_min=" << stats.firstOrderMismatchMin
              << "\n";
}

void printPhaseOracleMismatchFixture(
    const noninc::NonIncreasingPhaseOracleMismatchFixture &fixture,
    const char *label = "First") {
    std::cout << "\n" << label << " phase-oracle mismatch fixture\n"
              << "  threshold=" << fixture.threshold
              << " point=" << (fixture.beforeMaxToMin ? "before_max->min" : "after_phase")
              << " phase=" << phaseName(fixture.completedPhase)
              << " bbox=" << fixture.minRow << "," << fixture.minCol
              << "-" << fixture.maxRow << "," << fixture.maxCol
              << "\n";
    if (fixture.croppedImage) {
        std::cout << "  crop_size=" << fixture.croppedImage->getNumRows()
                  << "x" << fixture.croppedImage->getNumCols()
                  << "\n";
    }
    std::cout << "  local_only_support_in_crop=";
    printIntVectorPreview(fixture.localOnlySupportInCrop);
    std::cout << "\n  rebuild_only_support_in_crop=";
    printIntVectorPreview(fixture.rebuildOnlySupportInCrop);
    std::cout << "\n";
    if (fixture.localOnlyContext.supportSize > 0) {
        std::cout << "  local_only_context";
        printCandidateContext("local_only", fixture.localOnlyContext);
        std::cout << "\n";
    }
    if (fixture.rebuildOnlyContext.supportSize > 0) {
        std::cout << "  rebuild_only_context";
        printCandidateContext("rebuild_only", fixture.rebuildOnlyContext);
        std::cout << "\n";
    }
    if (!fixture.croppedImage) {
        return;
    }
    constexpr int maxPrintedCropArea = 400;
    if (fixture.croppedImage->getSize() > maxPrintedCropArea) {
        std::cout << "  crop_values=omitted"
                  << " area=" << fixture.croppedImage->getSize()
                  << " max_printed_area=" << maxPrintedCropArea
                  << "\n";
        return;
    }
    std::cout << "  crop_values:\n";
    for (int row = 0; row < fixture.croppedImage->getNumRows(); ++row) {
        std::cout << "    ";
        for (int col = 0; col < fixture.croppedImage->getNumCols(); ++col) {
            if (col != 0) {
                std::cout << " ";
            }
            const int pixelId = row * fixture.croppedImage->getNumCols() + col;
            std::cout << static_cast<int>((*fixture.croppedImage)[pixelId]);
        }
        std::cout << "\n";
    }
}

struct ImageWindow {
    bool valid = false;
    int minRow = 0;
    int minCol = 0;
    int maxRow = -1;
    int maxCol = -1;

    int numRows() const {
        return valid ? maxRow - minRow + 1 : 0;
    }

    int numCols() const {
        return valid ? maxCol - minCol + 1 : 0;
    }

    int area() const {
        return numRows() * numCols();
    }
};

struct ImageDiffSummary {
    bool differs = false;
    int differingPixels = 0;
    ImageWindow bbox;
    std::vector<PixelId> previewPixels;
};

void includePoint(ImageWindow &window, int row, int col) {
    if (!window.valid) {
        window.valid = true;
        window.minRow = row;
        window.maxRow = row;
        window.minCol = col;
        window.maxCol = col;
        return;
    }
    window.minRow = std::min(window.minRow, row);
    window.maxRow = std::max(window.maxRow, row);
    window.minCol = std::min(window.minCol, col);
    window.maxCol = std::max(window.maxCol, col);
}

ImageWindow expandedWindow(ImageWindow window,
                           int padding,
                           int numRows,
                           int numCols) {
    if (!window.valid) {
        return window;
    }
    window.minRow = std::max(0, window.minRow - padding);
    window.minCol = std::max(0, window.minCol - padding);
    window.maxRow = std::min(numRows - 1, window.maxRow + padding);
    window.maxCol = std::min(numCols - 1, window.maxCol + padding);
    return window;
}

ImageWindow mergeWindows(ImageWindow a,
                         ImageWindow b,
                         int numRows,
                         int numCols) {
    if (!a.valid) {
        return b;
    }
    if (!b.valid) {
        return a;
    }
    includePoint(a, b.minRow, b.minCol);
    includePoint(a, b.maxRow, b.maxCol);
    a.minRow = std::max(0, a.minRow);
    a.minCol = std::max(0, a.minCol);
    a.maxRow = std::min(numRows - 1, a.maxRow);
    a.maxCol = std::min(numCols - 1, a.maxCol);
    return a;
}

ImageWindow eventWindow(const noninc::NonIncreasingMergeEvent &event) {
    ImageWindow window;
    if (event.minRow < 0 ||
        event.minCol < 0 ||
        event.maxRow < event.minRow ||
        event.maxCol < event.minCol) {
        return window;
    }
    window.valid = true;
    window.minRow = event.minRow;
    window.minCol = event.minCol;
    window.maxRow = event.maxRow;
    window.maxCol = event.maxCol;
    return window;
}

ImageDiffSummary summarizeImageDiff(const ImageUInt8Ptr &a,
                                    const ImageUInt8Ptr &b) {
    if (!a || !b ||
        a->getNumRows() != b->getNumRows() ||
        a->getNumCols() != b->getNumCols()) {
        throw std::runtime_error("image diff requires same-size images");
    }

    ImageDiffSummary summary;
    constexpr std::size_t maxPreviewPixels = 12;
    const int numCols = a->getNumCols();
    for (PixelId pixelId = 0; pixelId < a->getSize(); ++pixelId) {
        if ((*a)[pixelId] == (*b)[pixelId]) {
            continue;
        }
        summary.differs = true;
        ++summary.differingPixels;
        if (summary.previewPixels.size() < maxPreviewPixels) {
            summary.previewPixels.push_back(pixelId);
        }
        includePoint(summary.bbox, pixelId / numCols, pixelId % numCols);
    }
    return summary;
}

void printWindowSummary(const char *label, const ImageWindow &window) {
    std::cout << "  " << label << "_bbox=";
    if (!window.valid) {
        std::cout << "none\n";
        return;
    }
    std::cout << window.minRow << "," << window.minCol
              << "-" << window.maxRow << "," << window.maxCol
              << " size=" << window.numRows() << "x" << window.numCols()
              << " area=" << window.area() << "\n";
}

void printPixelCoordPreview(const char *label,
                            const std::vector<PixelId> &pixels,
                            int numCols) {
    std::cout << "  " << label << "=[";
    const std::size_t limit = std::min<std::size_t>(pixels.size(), 12);
    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << ",";
        }
        const PixelId pixelId = pixels[i];
        std::cout << pixelId << "(" << (pixelId / numCols)
                  << "," << (pixelId % numCols) << ")";
    }
    if (pixels.size() > limit) {
        std::cout << ",...";
    }
    std::cout << "]\n";
}

void printPixelValueTriples(const char *label,
                            const ImageUInt8Ptr &previous,
                            const ImageUInt8Ptr &current,
                            const ImageUInt8Ptr &reference,
                            const std::vector<PixelId> &pixels) {
    std::cout << "  " << label << "=[";
    const std::size_t limit = std::min<std::size_t>(pixels.size(), 12);
    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << ",";
        }
        const PixelId pixelId = pixels[i];
        std::cout << "{id=" << pixelId
                  << ",prev=" << static_cast<int>((*previous)[pixelId])
                  << ",explicit=" << static_cast<int>((*current)[pixelId])
                  << ",naive=" << static_cast<int>((*reference)[pixelId])
                  << "}";
    }
    if (pixels.size() > limit) {
        std::cout << ",...";
    }
    std::cout << "]\n";
}

void printImageWindowValues(const char *label,
                            const ImageUInt8Ptr &image,
                            const ImageWindow &window,
                            int maxPrintedArea = 400) {
    if (!image || !window.valid) {
        std::cout << "  " << label << "_crop=none\n";
        return;
    }
    std::cout << "  " << label << "_crop_bbox="
              << window.minRow << "," << window.minCol
              << "-" << window.maxRow << "," << window.maxCol
              << " size=" << window.numRows() << "x" << window.numCols()
              << "\n";
    if (window.area() > maxPrintedArea) {
        std::cout << "  " << label << "_crop_values=omitted"
                  << " area=" << window.area()
                  << " max_printed_area=" << maxPrintedArea
                  << "\n";
        return;
    }
    std::cout << "  " << label << "_crop_values:\n";
    for (int row = window.minRow; row <= window.maxRow; ++row) {
        std::cout << "    ";
        for (int col = window.minCol; col <= window.maxCol; ++col) {
            if (col != window.minCol) {
                std::cout << " ";
            }
            const PixelId pixelId = row * image->getNumCols() + col;
            std::cout << static_cast<int>((*image)[pixelId]);
        }
        std::cout << "\n";
    }
}

void printDiffPointCrops(const ImageUInt8Ptr &previous,
                         const ImageUInt8Ptr &current,
                         const ImageUInt8Ptr &reference,
                         const std::vector<PixelId> &pixels,
                         int padding) {
    constexpr std::size_t maxDiffPointCrops = 6;
    if (!previous || pixels.empty()) {
        return;
    }
    const int numCols = previous->getNumCols();
    const int numRows = previous->getNumRows();
    const std::size_t limit = std::min<std::size_t>(pixels.size(), maxDiffPointCrops);
    std::cout << "\n  Diff-point crops"
              << " count=" << limit
              << " padding=" << padding
              << "\n";
    for (std::size_t i = 0; i < limit; ++i) {
        const PixelId pixelId = pixels[i];
        ImageWindow pointWindow;
        includePoint(pointWindow, pixelId / numCols, pixelId % numCols);
        pointWindow = expandedWindow(pointWindow, padding, numRows, numCols);
        std::cout << "  diff_point[" << i << "]"
                  << " id=" << pixelId
                  << " row=" << (pixelId / numCols)
                  << " col=" << (pixelId % numCols)
                  << "\n";
        printImageWindowValues("previous_explicit", previous, pointWindow);
        printImageWindowValues("current_explicit", current, pointWindow);
        printImageWindowValues("current_naive", reference, pointWindow);
    }
}

const char *phaseName(noninc::NonIncreasingPhase phase) {
    return phase == noninc::NonIncreasingPhase::MaxToMin ? "max->min" : "min->max";
}

noninc::NonIncreasingPhase parsePhase(const std::string &value) {
    if (value == "max->min" || value == "max-to-min" || value == "max2min") {
        return noninc::NonIncreasingPhase::MaxToMin;
    }
    if (value == "min->max" || value == "min-to-max" || value == "min2max") {
        return noninc::NonIncreasingPhase::MinToMax;
    }
    throw std::runtime_error("Invalid phase: " + value);
}

NodeId parseNodeId(const std::string &value) {
    char *endPtr = nullptr;
    const long parsed = std::strtol(value.c_str(), &endPtr, 10);
    if (endPtr == value.c_str() || *endPtr != '\0' || parsed < 0) {
        throw std::runtime_error("Invalid node id: " + value);
    }
    return static_cast<NodeId>(parsed);
}

const char *dualMaintenanceName(noninc::NonIncreasingDualMaintenance mode) {
    switch (mode) {
    case noninc::NonIncreasingDualMaintenance::PartialF:
        return "partial-f";
    case noninc::NonIncreasingDualMaintenance::ExplicitOnly:
        return "problem1-restricted";
    case noninc::NonIncreasingDualMaintenance::MaterializeThenExplicit:
        return "materialize-explicit";
    case noninc::NonIncreasingDualMaintenance::Problem1CloseRecompose:
        return "problem1-close-recompose";
    case noninc::NonIncreasingDualMaintenance::RebuildAfterPhase:
        return "phase-rebuild";
    case noninc::NonIncreasingDualMaintenance::OracleGuidedLocalRecomposition:
        return "oracle-local";
    case noninc::NonIncreasingDualMaintenance::StrippedThenOracleLocalRecomposition:
        return "stripped-oracle-local";
    case noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition:
        return "stripped-rebuild-local";
    case noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition:
        return "stripped-phase-rebuild-local";
    case noninc::NonIncreasingDualMaintenance::MergeNodeUpdateThenPhaseRebuild:
        return "merge-node-phase-rebuild";
    case noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild:
        return "merge-node-certified-phase-rebuild";
    case noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild:
        return "merge-node-certified-create-b-phase-rebuild";
    case noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild:
        return "merge-node-certified-create-b-mcpc-phase-rebuild";
    case noninc::NonIncreasingDualMaintenance::MergeNodeCleanUnchecked:
        return "merge-node-clean-unchecked";
    }
    return "unknown";
}

noninc::NonIncreasingStrippedLocalUpdateOptions strippedValidationOptions(const std::string &value) {
    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    if (value == "safe") {
        return options;
    }
    if (value == "no-reconstruction" || value == "no-reconstruction-checks") {
        options.checkInputPreReconstruction = false;
        options.buildFullUpdatedImage = false;
        options.checkPreMaterializationReconstruction = false;
        options.checkMaterializationTargetReconstruction = false;
        options.checkNoEmptyLiveNodesAfterCanonicalization = false;
        options.checkFinalDualPrimalReconstruction = false;
        options.releaseEmptyNodesAfterPrimalMerge = false;
        return options;
    }
    if (value == "unchecked") {
        options.checkInputPreReconstruction = false;
        options.buildFullUpdatedImage = false;
        options.auditPlanOracleDependencies = false;
        options.auditMaterializationOracleDependencies = false;
        options.checkPreMaterializationReconstruction = false;
        options.checkMaterializationTargetReconstruction = false;
        options.checkNoEmptyLiveNodesAfterCanonicalization = false;
        options.checkFinalDualPrimalReconstruction = false;
        options.releaseEmptyNodesAfterPrimalMerge = false;
        return options;
    }
    throw std::runtime_error("Invalid stripped validation mode: " + value);
}

noninc::NonIncreasingStrippedCommitMode strippedCommitModeFromString(const std::string &value) {
    if (value == "safe-copy" || value == "copy") {
        return noninc::NonIncreasingStrippedCommitMode::SafeCopy;
    }
    if (value == "direct") {
        return noninc::NonIncreasingStrippedCommitMode::Direct;
    }
    throw std::runtime_error("Invalid stripped commit mode: " + value);
}

void printIntVectorPreview(const std::vector<int> &values) {
    std::cout << "[";
    const std::size_t limit = std::min<std::size_t>(values.size(), 8);
    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << ",";
        }
        std::cout << values[i];
    }
    if (values.size() > limit) {
        std::cout << ",...";
    }
    std::cout << "]";
}

void printEmptyCreatedNodeDiagnosticsPreview(
    const std::vector<noninc::LocalEmptyCreatedNodeDiagnostic> &nodes) {
    std::cout << "[";
    const std::size_t limit = std::min<std::size_t>(nodes.size(), 8);
    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << ",";
        }
        const auto &node = nodes[i];
        std::cout << "{id=" << node.nodeId
                  << ",plan=" << node.planIndex
                  << ",parent=" << node.parentNodeId
                  << ",level=" << node.level
                  << ",children=" << node.childCount
                  << ",moved_pp=" << node.movedDirectProperPartCount
                  << ",moved_children=" << node.movedChildNodeCount
                  << "}";
    }
    if (nodes.size() > limit) {
        std::cout << ",...";
    }
    std::cout << "]";
}

void printCandidateContext(const char *prefix,
                           const noninc::NonIncreasingCandidateHierarchyContext &context) {
    std::cout << " " << prefix << "_node=" << context.nodeId
              << " " << prefix << "_level=" << context.altitude
              << " " << prefix << "_support=" << context.supportSize
              << " " << prefix << "_bbox="
              << context.minRow << "," << context.minCol
              << "-" << context.maxRow << "," << context.maxCol
              << " " << prefix << "_parent_level=" << context.parentAltitude
              << " " << prefix << "_parent_support=" << context.parentSupportSize
              << " " << prefix << "_children=" << context.childCount
              << " " << prefix << "_first_child_level=" << context.firstChildAltitude
              << " " << prefix << "_first_child_support=" << context.firstChildSupportSize
              << " " << prefix << "_pixels=";
    printIntVectorPreview(context.supportPixels);
}

ImageWindow pixelWindow(std::span<const PixelId> pixels, int numCols) {
    ImageWindow window;
    if (numCols <= 0) {
        return window;
    }
    for (PixelId pixelId : pixels) {
        if (pixelId < 0) {
            continue;
        }
        includePoint(window, pixelId / numCols, pixelId % numCols);
    }
    return window;
}

bool isIncidentClosureComparison(
    const noninc::NonIncreasingPhaseOracleComparison &comparison) {
    return comparison.releasePredicateEvaluated &&
           comparison.releaseFocusAvailable &&
           !comparison.releasePredicateSatisfied &&
           comparison.releaseFocusPixelCount > 0 &&
           comparison.releaseClosureNodeCount > 0 &&
           comparison.releaseClosurePixelCount > 0 &&
           comparison.localReconstructionImage;
}

const noninc::NonIncreasingPhaseOracleComparison *firstIncidentClosureComparison(
    const std::vector<noninc::NonIncreasingPhaseOracleComparison> &comparisons) {
    for (const auto &comparison : comparisons) {
        if (isIncidentClosureComparison(comparison)) {
            return &comparison;
        }
    }
    return nullptr;
}

void printIncidentClosureFixture(
    const noninc::NonIncreasingPhaseOracleComparison &comparison) {
    std::cout << "\nIncident closure fixture\n"
              << "  threshold=" << comparison.threshold
              << " point="
              << (comparison.beforeMaxToMin ? "before_max->min" : "after_phase")
              << " phase=" << phaseName(comparison.completedPhase)
              << "\n"
              << "  release_h1="
              << (comparison.releaseH1ReconstructionMatches ? "true" : "false")
              << " release_h2="
              << (comparison.releaseH2CandidateSetMatches ? "true" : "false")
              << " release_h3="
              << (comparison.releaseH3EnvelopeStructureMatches ? "true" : "false")
              << " release_h4_boundary="
              << (comparison.releaseH4BoundaryStructureMatches ? "true" : "false")
              << " release_h4_incidence="
              << (comparison.releaseH4BoundaryIncidenceMatches ? "true" : "false")
              << " release_h5_order="
              << (comparison.releaseH5CandidateScanOrderMatches ? "true" : "false")
              << "\n"
              << "  focus_nodes=" << comparison.releaseFocusNodeCount
              << " focus_pixels=" << comparison.releaseFocusPixelCount
              << " closure_nodes=" << comparison.releaseClosureNodeCount
              << " closure_pixels=" << comparison.releaseClosurePixelCount
              << "\n"
              << "  local_envelope_nodes=" << comparison.localEnvelopeNodeCount
              << " rebuild_envelope_nodes=" << comparison.rebuildEnvelopeNodeCount
              << " local_boundary_nodes=" << comparison.localBoundaryNodeCount
              << " rebuild_boundary_nodes=" << comparison.rebuildBoundaryNodeCount
              << " local_boundary_incidences="
              << comparison.localBoundaryIncidenceCount
              << " rebuild_boundary_incidences="
              << comparison.rebuildBoundaryIncidenceCount
              << "\n";

    if (comparison.firstLocalOnlyContext.supportSize > 0) {
        std::cout << "  local_only_context";
        printCandidateContext("local_only", comparison.firstLocalOnlyContext);
        std::cout << "\n";
    }
    if (comparison.firstRebuildOnlyContext.supportSize > 0) {
        std::cout << "  rebuild_only_context";
        printCandidateContext("rebuild_only", comparison.firstRebuildOnlyContext);
        std::cout << "\n";
    }

    const int numCols = comparison.localReconstructionImage->getNumCols();
    printPixelCoordPreview(
        "focus_pixels",
        comparison.releaseFocusPixels,
        numCols);
    printPixelCoordPreview(
        "closure_pixels",
        comparison.releaseClosurePixels,
        numCols);
    const auto window = pixelWindow(
        std::span<const PixelId>(comparison.releaseClosurePixels),
        numCols);
    printWindowSummary("closure", window);
    printImageWindowValues(
        "closure_local_reconstruction",
        comparison.localReconstructionImage,
        window,
        1200);
}

void printNodeVectorPreview(const std::vector<NodeId> &values) {
    std::cout << "[";
    const std::size_t limit = std::min<std::size_t>(values.size(), 8);
    for (std::size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            std::cout << ",";
        }
        std::cout << values[i];
    }
    if (values.size() > limit) {
        std::cout << ",...";
    }
    std::cout << "]";
}

void printMergeEvent(const noninc::NonIncreasingMergeEvent &event,
                     std::size_t index) {
    std::cout << "  event[" << index << "]"
              << " phase=" << phaseName(event.phase)
              << " threshold=" << event.threshold
              << " node=" << event.nodeId
              << " level=" << event.nodeLevel
              << " parent=" << event.parentId
              << " b=" << event.parentLevel
              << " attr=" << event.attributeValue
              << " children=" << event.numChildren
              << " direct_pp=" << event.directProperParts
              << " components=" << event.componentCount
              << " |C|=" << event.totalComponentPixels
              << " C_bbox=" << event.minRow << "," << event.minCol
              << "-" << event.maxRow << "," << event.maxCol
              << " dual_rebuild=" << (event.dualRebuiltAfterLocalMismatch ? "true" : "false")
              << " stripped_attempted=" << (event.strippedLocalAttempted ? "true" : "false")
              << " stripped_accepted=" << (event.strippedLocalAccepted ? "true" : "false")
              << " stripped_applied=" << (event.strippedLocalApplied ? "true" : "false")
              << " stripped_direct_commit=" << (event.strippedDirectCommit ? "true" : "false")
              << " stripped_plan_audit_checked=" << (event.strippedLocalPlanAuditChecked ? "true" : "false")
              << " stripped_target=" << (event.strippedLocalMatchesTarget ? "true" : "false")
              << " stripped_target_checked=" << (event.strippedLocalTargetChecked ? "true" : "false")
              << " stripped_structure_checked="
              << (event.strippedLocalStructureChecked ? "true" : "false")
              << " stripped_structure_matches_rebuild="
              << (event.strippedLocalStructureMatchesRebuild ? "true" : "false")
              << " stripped_plan_oracle_free=" << (event.strippedLocalPlanOracleFree ? "true" : "false")
              << " stripped_materialization_audit_checked="
              << (event.strippedMaterializationAuditChecked ? "true" : "false")
              << " stripped_materialization_oracle_free="
              << (event.strippedMaterializationInputOracleFree ? "true" : "false")
              << " stripped_final_checked=" << (event.strippedFinalReconstructionChecked ? "true" : "false")
              << " stripped_dual_matches="
              << (event.strippedLocalDualMatchesPrimalAfterMerge ? "true" : "false")
              << " stripped_fallback=" << (event.strippedLocalFallbackToOracle ? "true" : "false")
              << " stripped_rebuild_fallback=" << (event.strippedLocalFallbackToRebuild ? "true" : "false")
              << " stripped_phase_rebuild_fallback=" << (event.strippedLocalFallbackToPhaseRebuild ? "true" : "false")
              << " stripped_phase_rebuild_skip=" << (event.strippedLocalSkippedByPhaseRebuild ? "true" : "false")
              << " partial_f_attempted=" << (event.partialFAttempted ? "true" : "false")
              << " partial_f_accepted=" << (event.partialFAccepted ? "true" : "false")
              << " partial_f_rejected=" << (event.partialFRejected ? "true" : "false")
              << " partial_f_mutation_attempted="
              << (event.partialFMutationAttempted ? "true" : "false")
              << " partial_f_mcpc_live="
              << (event.partialFMcpcLiveUpdate ? "true" : "false")
              << " partial_f_mcpc_hidden_support_redundant_live="
              << (event.partialFMcpcHiddenSupportRedundantLiveUpdate ? "true" : "false")
              << " partial_f_mcpc_parent_link_bridge_live="
              << (event.partialFMcpcParentLinkBridgeLiveUpdate ? "true" : "false")
              << " partial_f_mcpc_stable_child_bridge_live="
              << (event.partialFMcpcStableChildBridgeLiveUpdate ? "true" : "false")
              << " partial_f_mcpc_child_link_bridge_live="
              << (event.partialFMcpcChildLinkBridgeLiveUpdate ? "true" : "false")
              << " partial_f_mcpc_empty_owner_release_live="
              << (event.partialFMcpcEmptyOwnerReleaseLiveUpdate ? "true" : "false")
              << " partial_f_mcpc_stable_child_move_live="
              << (event.partialFMcpcStableChildMoveLiveUpdate ? "true" : "false")
              << " merge_node_attempted=" << (event.mergeNodeLocalAttempted ? "true" : "false")
              << " merge_node_accepted=" << (event.mergeNodeLocalAccepted ? "true" : "false")
              << " merge_node_applied=" << (event.mergeNodeLocalApplied ? "true" : "false")
              << " merge_node_cert_checked="
              << (event.mergeNodeLocalCertificationChecked ? "true" : "false")
              << " merge_node_cert_passed="
              << (event.mergeNodeLocalCertificationPassed ? "true" : "false")
              << " merge_node_direct_commit="
              << (event.mergeNodeLocalDirectCommit ? "true" : "false")
              << " merge_node_target_creation_allowed="
              << (event.mergeNodeLocalTargetCreationAllowed ? "true" : "false")
              << " merge_node_target_creation_required="
              << (event.mergeNodeLocalTargetCreationRequired ? "true" : "false")
              << " merge_node_pre_dual_matches="
              << (event.mergeNodeLocalDualMatchesPrimalBeforeMerge ? "true" : "false")
              << " merge_node_dual_matches="
              << (event.mergeNodeLocalDualMatchesPrimalAfterMerge ? "true" : "false")
              << " merge_node_phase_rebuild_fallback="
              << (event.mergeNodeLocalFallbackToPhaseRebuild ? "true" : "false")
              << " merge_node_phase_rebuild_skip="
              << (event.mergeNodeLocalSkippedByPhaseRebuild ? "true" : "false")
              << " merge_node_components=" << event.mergeNodeLocalComponentCount
              << " merge_node_rejected_components="
              << event.mergeNodeLocalRejectedComponents
              << " merge_node_problem1_preflight_rejected_components="
              << event.mergeNodeLocalProblem1PreflightRejectedComponents
              << " merge_node_problem1_preflight_rejected_by_published="
              << event.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem
              << " merge_node_problem1_preflight_rejected_by_updater="
              << event.mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation
              << " merge_node_partial_owner_components="
              << event.mergeNodeLocalPartialOwnerComponents
              << " merge_node_raw_missing_target_b_components="
              << event.mergeNodeLocalRawMissingTargetComponents
              << " merge_node_certified_missing_target_b_components="
              << event.mergeNodeLocalMissingTargetComponents
              << " merge_node_partial_owner_and_missing_target_b_components="
              << event.mergeNodeLocalPartialOwnerAndMissingTargetComponents
              << " merge_node_missing_target_components="
              << event.mergeNodeLocalMissingTargetComponents
              << " merge_node_pixels=" << event.mergeNodeLocalPixelCount
              << " merge_node_copy_us=" << event.mergeNodeLocalCopyMicros
              << " merge_node_cert_us=" << event.mergeNodeLocalCertificationMicros
              << " merge_node_component_collect_us="
              << event.mergeNodeComponentCollectMicros
              << " merge_node_update_us=" << event.mergeNodeLocalUpdateMicros
              << " merge_node_primal_merge_us=" << event.mergeNodePrimalMergeMicros
              << " merge_node_empty_release_us=" << event.mergeNodeEmptyReleaseMicros
              << " stripped_recomp_nodes=" << event.strippedLocalRecompositionNodes
              << " stripped_created_nodes=" << event.strippedLocalCreatedNodes
              << " stripped_local_pixels=" << event.strippedLocalPixelCount
              << " stripped_created_levels=";
    printIntVectorPreview(event.strippedCreatedNodeLevels);
    std::cout
              << " stripped_created_plan_indices=";
    printIntVectorPreview(event.strippedCreatedNodePlanIndices);
    std::cout
              << " stripped_created_connector_nodes="
              << event.strippedCreatedConnectorNodes
              << " stripped_empty_created_nodes="
              << event.strippedEmptyCreatedNodesAfterCanonicalization
              << " stripped_empty_connector_nodes="
              << event.strippedEmptyConnectorNodesAfterCanonicalization
              << " stripped_released_owners=" << event.strippedLocalReleasedSourceOwners
              << " stripped_copy_us=" << event.strippedLocalCopyMicros
              << " stripped_update_us=" << event.strippedLocalUpdateMicros
              << " stripped_fallback_us=" << event.strippedOracleFallbackMicros
              << " stripped_rebuild_fallback_us=" << event.strippedRebuildFallbackMicros
              << " stripped_input_us=" << event.strippedInputExtractionMicros
              << " stripped_partition_us=" << event.strippedPartitionPlanMicros
              << " stripped_recomposition_us=" << event.strippedRecompositionPlanMicros
              << " stripped_materialization_us=" << event.strippedMaterializationMicros
              << " stripped_mat_canonicalize_us="
              << event.strippedMaterializationCanonicalizationMicros
              << " stripped_mat_release_us=" << event.strippedMaterializationReleaseMicros
              << " stripped_mat_target_check_us="
              << event.strippedMaterializationTargetCheckMicros
              << " stripped_primal_release_us=" << event.strippedPrimalMergeReleaseMicros
              << " stripped_final_check_us=" << event.strippedFinalReconstructionCheckMicros
              << " oracle_local=" << (event.oracleGuidedLocalApplied ? "true" : "false")
              << " oracle_target=" << (event.oracleGuidedLocalMatchesTarget ? "true" : "false")
              << " oracle_structure=" << (event.oracleGuidedLocalStructureMatchesOracle ? "true" : "false")
              << " oracle_recomp_nodes=" << event.oracleGuidedRecompositionNodes
              << " oracle_created_nodes=" << event.oracleGuidedCreatedNodes
              << " oracle_released_owners=" << event.oracleGuidedReleasedSourceOwners
              << " component_sizes=";
    printIntVectorPreview(event.componentSizes);
    std::cout << " nodeCa=";
    printNodeVectorPreview(event.componentNodeCa);
    std::cout << " altitudeCa=";
    printIntVectorPreview(event.componentAltitudeCa);
    std::cout << " C_pixels=";
    printIntVectorPreview(event.affectedPixelsPreview);
    if (!event.strippedEmptyCreatedNodeDiagnostics.empty()) {
        std::cout << " empty_created_nodes=";
        printEmptyCreatedNodeDiagnosticsPreview(event.strippedEmptyCreatedNodeDiagnostics);
    }
    if (!event.strippedLocalFailureSummary.empty()) {
        std::cout << " stripped_failure_summary=\"" << event.strippedLocalFailureSummary << "\"";
    }
    if (!event.mergeNodeLocalFailureSummary.empty()) {
        std::cout << " merge_node_failure_summary=\"" << event.mergeNodeLocalFailureSummary << "\"";
    }
    if (!event.mergeNodeLocalProblem1FailureMode.empty()) {
        std::cout << " merge_node_problem1_failure_mode=\""
                  << event.mergeNodeLocalProblem1FailureMode << "\"";
    }
    if (!event.strippedLocalStructureMismatchSummary.empty()) {
        std::cout << " stripped_structure_mismatch=\""
                  << event.strippedLocalStructureMismatchSummary << "\"";
    }
    std::cout << "\n";
}

bool bboxIntersects(int aMinRow,
                    int aMinCol,
                    int aMaxRow,
                    int aMaxCol,
                    int bMinRow,
                    int bMinCol,
                    int bMaxRow,
                    int bMaxCol) {
    if (aMinRow < 0 || aMinCol < 0 || aMaxRow < aMinRow || aMaxCol < aMinCol ||
        bMinRow < 0 || bMinCol < 0 || bMaxRow < bMinRow || bMaxCol < bMinCol) {
        return false;
    }
    return aMinRow <= bMaxRow &&
           aMaxRow >= bMinRow &&
           aMinCol <= bMaxCol &&
           aMaxCol >= bMinCol;
}

void printPhaseOracleFixtureEventLinks(
    const noninc::NonIncreasingPhaseOracleMismatchFixture &fixture,
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    constexpr std::size_t maxPrintedEvents = 5;
    std::size_t matches = 0;
    std::cout << "  overlapping_C_events:\n";
    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto &event = events[i];
        if (!bboxIntersects(event.minRow,
                            event.minCol,
                            event.maxRow,
                            event.maxCol,
                            fixture.minRow,
                            fixture.minCol,
                            fixture.maxRow,
                            fixture.maxCol)) {
            continue;
        }
        ++matches;
        if (matches <= maxPrintedEvents) {
            printMergeEvent(event, i);
        }
    }
    if (matches == 0) {
        std::cout << "    none\n";
    } else if (matches > maxPrintedEvents) {
        std::cout << "    ... " << (matches - maxPrintedEvents)
                  << " more overlapping events omitted\n";
    }
}

void printMergeEventSummary(const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    MergeEventStats stats;
    addEventStats(stats, events);
    std::size_t multiComponentEvents = 0;
    std::size_t maxComponentsIndex = 0;
    std::size_t maxCIndex = 0;
    int maxComponents = 0;
    int maxC = 0;

    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto &event = events[i];
        if (event.componentCount > 1) {
            ++multiComponentEvents;
        }
        if (event.componentCount > maxComponents) {
            maxComponents = event.componentCount;
            maxComponentsIndex = i;
        }
        if (event.totalComponentPixels > maxC) {
            maxC = event.totalComponentPixels;
            maxCIndex = i;
        }
    }

    std::cout << "\nMerge-event diagnostics\n";
    std::cout << "  events=" << events.size()
              << " multi_component_events=" << multiComponentEvents
              << " child_events=" << stats.childEvents
              << " dual_rebuilds=" << stats.dualRebuildEvents
              << " phase_dual_rebuilds=" << stats.phaseDualRebuilds
              << " phase_candidate_guard_rebuilds=" << stats.phaseCandidateGuardRebuilds
              << " phase_candidate_guard_canonicalize_changes="
              << stats.phaseCandidateGuardCanonicalizationChanges
              << " stripped_attempts=" << stats.strippedLocalAttempts
              << " stripped_accepted=" << stats.strippedLocalAccepted
              << " stripped_applied=" << stats.strippedLocalApplied
              << " stripped_direct_commits=" << stats.strippedDirectCommits
              << " stripped_plan_audits=" << stats.strippedLocalPlanAuditChecks
              << " stripped_target_checks=" << stats.strippedLocalTargetChecks
              << " stripped_target_matches=" << stats.strippedLocalTargetMatches
              << " stripped_plan_oracle_free=" << stats.strippedLocalOracleFreePlans
              << " stripped_materialization_audits="
              << stats.strippedMaterializationAuditChecks
              << " stripped_materialization_oracle_free="
              << stats.strippedLocalOracleFreeMaterializations
              << " stripped_final_checks=" << stats.strippedFinalReconstructionChecks
              << " stripped_fallbacks=" << stats.strippedLocalFallbacks
              << " stripped_rebuild_fallbacks=" << stats.strippedLocalRebuildFallbacks
              << " stripped_phase_rebuild_fallbacks=" << stats.strippedLocalPhaseRebuildFallbacks
              << " stripped_phase_rebuild_skips=" << stats.strippedLocalPhaseRebuildSkips
              << " partial_f_attempts=" << stats.partialFAttempts
              << " partial_f_accepted=" << stats.partialFAccepted
              << " partial_f_rejected=" << stats.partialFRejected
              << " partial_f_mutation_attempted=" << stats.partialFMutationAttempted
              << " partial_f_mcpc_live_accepted=" << stats.partialFMcpcLiveAccepted
              << " partial_f_mcpc_hidden_support_redundant_live_accepted="
              << stats.partialFMcpcHiddenSupportRedundantLiveAccepted
              << " partial_f_mcpc_parent_link_bridge_live_accepted="
              << stats.partialFMcpcParentLinkBridgeLiveAccepted
              << " partial_f_mcpc_stable_child_bridge_live_accepted="
              << stats.partialFMcpcStableChildBridgeLiveAccepted
              << " partial_f_mcpc_child_link_bridge_live_accepted="
              << stats.partialFMcpcChildLinkBridgeLiveAccepted
              << " partial_f_mcpc_empty_owner_release_live_accepted="
              << stats.partialFMcpcEmptyOwnerReleaseLiveAccepted
              << " partial_f_mcpc_stable_child_move_live_accepted="
              << stats.partialFMcpcStableChildMoveLiveAccepted
              << partialFRejectStatsText(stats)
              << " merge_node_attempts=" << stats.mergeNodeLocalAttempts
              << " merge_node_accepted=" << stats.mergeNodeLocalAccepted
              << " merge_node_applied=" << stats.mergeNodeLocalApplied
              << " merge_node_cert_checks=" << stats.mergeNodeLocalCertificationChecks
              << " merge_node_cert_passed=" << stats.mergeNodeLocalCertificationPassed
              << " merge_node_direct_commits=" << stats.mergeNodeLocalDirectCommits
              << " merge_node_target_creation_allowed_events="
              << stats.mergeNodeLocalTargetCreationAllowedEvents
              << " merge_node_target_creation_required_events="
              << stats.mergeNodeLocalTargetCreationRequiredEvents
              << " merge_node_pre_dual_matches=" << stats.mergeNodeLocalPreDualMatches
              << " merge_node_dual_matches=" << stats.mergeNodeLocalDualMatches
              << " merge_node_phase_rebuild_fallbacks="
              << stats.mergeNodeLocalPhaseRebuildFallbacks
              << " merge_node_phase_rebuild_skips="
              << stats.mergeNodeLocalPhaseRebuildSkips
              << " merge_node_rejected_components="
              << stats.mergeNodeLocalRejectedComponents
              << " merge_node_problem1_preflight_rejected_components="
              << stats.mergeNodeLocalProblem1PreflightRejectedComponents
              << " merge_node_problem1_preflight_rejected_by_published="
              << stats.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem
              << " merge_node_problem1_preflight_rejected_by_updater="
              << stats.mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation
              << " merge_node_partial_owner_components="
              << stats.mergeNodeLocalPartialOwnerComponents
              << " merge_node_raw_missing_target_b_components="
              << stats.mergeNodeLocalRawMissingTargetComponents
              << " merge_node_certified_missing_target_b_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_partial_owner_and_missing_target_b_components="
              << stats.mergeNodeLocalPartialOwnerAndMissingTargetComponents
              << " merge_node_missing_target_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_pixels=" << stats.mergeNodeLocalPixels
              << " merge_node_copy_ms=" << common::microsToMillis(stats.mergeNodeLocalCopyMicros)
              << " merge_node_cert_ms="
              << common::microsToMillis(stats.mergeNodeLocalCertificationMicros)
              << " merge_node_component_collect_ms="
              << common::microsToMillis(stats.mergeNodeComponentCollectMicros)
              << " merge_node_update_ms=" << common::microsToMillis(stats.mergeNodeLocalUpdateMicros)
              << " merge_node_primal_merge_ms="
              << common::microsToMillis(stats.mergeNodePrimalMergeMicros)
              << " merge_node_empty_release_ms="
              << common::microsToMillis(stats.mergeNodeEmptyReleaseMicros)
              << " stripped_created_nodes=" << stats.strippedLocalCreatedNodes
              << " stripped_recomposition_nodes=" << stats.strippedLocalRecompositionNodes
              << " stripped_local_pixels=" << stats.strippedLocalPixels
              << " max_stripped_local_pixels=" << stats.maxStrippedLocalPixels
              << " stripped_empty_created_nodes=" << stats.strippedEmptyCreatedNodes
              << " stripped_empty_connector_nodes=" << stats.strippedEmptyConnectorNodes
              << " stripped_released_owners=" << stats.strippedLocalReleasedOwners
              << " stripped_copy_ms=" << common::microsToMillis(stats.strippedLocalCopyMicros)
              << " stripped_update_ms=" << common::microsToMillis(stats.strippedLocalUpdateMicros)
              << " stripped_fallback_ms=" << common::microsToMillis(stats.strippedOracleFallbackMicros)
              << " stripped_rebuild_fallback_ms=" << common::microsToMillis(stats.strippedRebuildFallbackMicros)
              << " stripped_input_ms=" << common::microsToMillis(stats.strippedInputExtractionMicros)
              << " stripped_partition_ms=" << common::microsToMillis(stats.strippedPartitionPlanMicros)
              << " stripped_recomposition_ms=" << common::microsToMillis(stats.strippedRecompositionPlanMicros)
              << " stripped_plan_audit_ms=" << common::microsToMillis(stats.strippedPlanAuditMicros)
              << " stripped_mat_context_ms="
              << common::microsToMillis(stats.strippedMaterializationContextMicros)
              << " stripped_mat_audit_ms="
              << common::microsToMillis(stats.strippedMaterializationAuditMicros)
              << " stripped_materialization_ms="
              << common::microsToMillis(stats.strippedMaterializationMicros)
              << " stripped_mat_validation_ms="
              << common::microsToMillis(stats.strippedMaterializationValidationMicros)
              << " stripped_mat_create_nodes_ms="
              << common::microsToMillis(stats.strippedMaterializationCreateNodesMicros)
              << " stripped_mat_attach_nodes_ms="
              << common::microsToMillis(stats.strippedMaterializationAttachNodesMicros)
              << " stripped_mat_move_pp_ms="
              << common::microsToMillis(stats.strippedMaterializationMoveProperPartsMicros)
              << " stripped_mat_child_bridge_ms="
              << common::microsToMillis(stats.strippedMaterializationChildBridgeMicros)
              << " stripped_mat_structural_repair_ms="
              << common::microsToMillis(stats.strippedMaterializationStructuralRepairMicros)
              << " stripped_mat_release_ms="
              << common::microsToMillis(stats.strippedMaterializationReleaseMicros)
              << " stripped_mat_root_promotion_ms="
              << common::microsToMillis(stats.strippedMaterializationRootPromotionMicros)
              << " stripped_mat_canonicalize_ms="
              << common::microsToMillis(stats.strippedMaterializationCanonicalizationMicros);
    printMaxCanonicalizationEventFields(stats);
    std::cout << " stripped_mat_target_check_ms="
              << common::microsToMillis(stats.strippedMaterializationTargetCheckMicros)
              << " stripped_post_snapshot_ms="
              << common::microsToMillis(stats.strippedPostLocalSnapshotMicros)
              << " stripped_primal_release_ms="
              << common::microsToMillis(stats.strippedPrimalMergeReleaseMicros);
    printMaxPrimalMergeEventFields(stats);
    std::cout << " stripped_final_check_ms="
              << common::microsToMillis(stats.strippedFinalReconstructionCheckMicros)
              << " oracle_local_events=" << stats.oracleLocalEvents
              << " oracle_local_target_matches=" << stats.oracleLocalTargetMatches
              << " oracle_local_structure_matches=" << stats.oracleLocalStructureMatches
              << " oracle_local_structure_mismatches=" << stats.oracleLocalStructureMismatches
              << " oracle_created_nodes=" << stats.oracleLocalCreatedNodes
              << " oracle_recomposition_nodes=" << stats.oracleLocalRecompositionNodes
              << " oracle_released_owners=" << stats.oracleLocalReleasedOwners
              << " total_|C|=" << stats.totalComponentPixels
              << " max_components=" << maxComponents
              << " max_|C|=" << maxC
              << "\n";

    const std::size_t firstLimit = std::min<std::size_t>(events.size(), 10);
    std::cout << "  first_events:\n";
    for (std::size_t i = 0; i < firstLimit; ++i) {
        printMergeEvent(events[i], i);
    }

    if (!events.empty()) {
        std::cout << "  largest_component_count_event:\n";
        printMergeEvent(events[maxComponentsIndex], maxComponentsIndex);
        std::cout << "  largest_C_event:\n";
        printMergeEvent(events[maxCIndex], maxCIndex);
    }

    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].componentCount > 1) {
            std::cout << "  first_multi_component_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }

    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_live_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcHiddenSupportRedundantLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_hidden_support_redundant_live_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcParentLinkBridgeLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_parent_link_bridge_live_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcStableChildBridgeLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_stable_child_bridge_live_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcChildLinkBridgeLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_child_link_bridge_live_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcStableChildMoveLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_stable_child_move_live_event:\n";
            printMergeEvent(events[i], i);
            break;
        }
    }
}

void printFirstOracleStructureMismatchEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto &event = events[i];
        if (event.oracleGuidedLocalApplied &&
            !event.oracleGuidedLocalStructureMatchesOracle) {
            std::cout << "  first_oracle_structure_mismatch_event:\n";
            printMergeEvent(event, i);
            if (!event.oracleGuidedStructureMismatchSummary.empty()) {
                std::cout << "    structure_mismatch_summary:\n"
                          << "    " << event.oracleGuidedStructureMismatchSummary << "\n";
            }
            return;
        }
    }
}

void printFirstPhaseRebuildFallbackEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        const auto &event = events[i];
        if ((event.strippedLocalFallbackToPhaseRebuild &&
             !event.strippedLocalSkippedByPhaseRebuild) ||
            (event.mergeNodeLocalFallbackToPhaseRebuild &&
             !event.mergeNodeLocalSkippedByPhaseRebuild)) {
            std::cout << "  first_phase_rebuild_fallback_event:\n";
            printMergeEvent(event, i);
            return;
        }
    }
}

void printFirstPartialFMcpcLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

void printFirstPartialFMcpcParentLinkBridgeLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcParentLinkBridgeLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_parent_link_bridge_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

void printFirstPartialFMcpcHiddenSupportRedundantLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcHiddenSupportRedundantLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_hidden_support_redundant_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

void printFirstPartialFMcpcStableChildBridgeLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcStableChildBridgeLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_stable_child_bridge_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

void printFirstPartialFMcpcChildLinkBridgeLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcChildLinkBridgeLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_child_link_bridge_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

void printFirstPartialFMcpcEmptyOwnerReleaseLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcEmptyOwnerReleaseLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_empty_owner_release_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

void printFirstPartialFMcpcStableChildMoveLiveEvent(
    const std::vector<noninc::NonIncreasingMergeEvent> &events) {
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcStableChildMoveLiveUpdate) {
            std::cout << "  first_partial_f_mcpc_stable_child_move_live_event:\n";
            printMergeEvent(events[i], i);
            return;
        }
    }
}

struct PrefixCompareResult {
    bool matches = false;
    int requestedEvents = 0;
    int explicitEvents = 0;
    long long explicitMicros = 0;
    long long naiveMicros = 0;
    ImageUInt8Ptr explicitOutput;
    ImageUInt8Ptr naiveOutput;
    MergeEventStats stats;
    std::vector<noninc::NonIncreasingMergeEvent> events;
    bool localUpdateRejectionStopReached = false;
};

noninc_debug::NonIncreasingCasfDiagnosticsOptions makeBaseDiagnosticsOptions(
    const noninc::NonIncreasingStrippedLocalUpdateOptions &strippedOptions,
    noninc::NonIncreasingStrippedCommitMode strippedCommitMode,
    bool leafFastPathEnabled,
    bool phaseCandidateGuardEnabled,
    bool strictDirectV1PrecheckEnabled,
    bool mcpcNoCopyDirectCommitEnabled,
    bool mcpcValidatedDirectCommitEnabled) {
    noninc_debug::NonIncreasingCasfDiagnosticsOptions options;
    options.strippedLocalUpdateOptions = strippedOptions;
    options.strippedCommitMode = strippedCommitMode;
    options.leafExplicitFastPath = leafFastPathEnabled;
    options.phaseCandidateEquivalenceGuard = phaseCandidateGuardEnabled;
    options.mcpcDirectV1StrictPrecheck = strictDirectV1PrecheckEnabled;
    options.mcpcNoCopyDirectCommit = mcpcNoCopyDirectCommitEnabled;
    options.mcpcValidatedDirectCommit = mcpcValidatedDirectCommitEnabled;
    return options;
}

PrefixCompareResult runPrefixCompare(ImageUInt8Ptr input,
                                     const AdjacencyRelationPtr &adjacency,
                                     double radiusAdj,
                                     float threshold,
                                     noninc::NonIncreasingDualMaintenance dualMaintenance,
                                     noninc::NonIncreasingAttributeKind attributeKind,
                                     const noninc::NonIncreasingStrippedLocalUpdateOptions &strippedOptions,
                                     noninc::NonIncreasingStrippedCommitMode strippedCommitMode,
                                     bool leafFastPathEnabled,
                                     bool phaseCandidateGuardEnabled,
                                     bool strictDirectV1PrecheckEnabled,
                                     bool mcpcNoCopyDirectCommitEnabled,
                                     bool mcpcValidatedDirectCommitEnabled,
                                     int maxEvents) {
    PrefixCompareResult result;
    result.requestedEvents = maxEvents;

    noninc::NonIncreasingCasfKernel runner(
        input,
        radiusAdj,
        dualMaintenance,
        attributeKind);
    auto diagnostics = makeBaseDiagnosticsOptions(
        strippedOptions,
        strippedCommitMode,
        leafFastPathEnabled,
        phaseCandidateGuardEnabled,
        strictDirectV1PrecheckEnabled,
        mcpcNoCopyDirectCommitEnabled,
        mcpcValidatedDirectCommitEnabled);
    diagnostics.collectMergeEvents = true;
    diagnostics.maxAppliedMergeEventsPerThreshold = maxEvents;
    noninc_debug::applyNonIncreasingCasfDiagnostics(runner, diagnostics);

    const auto explicitStart = std::chrono::steady_clock::now();
    runner.applyThreshold(threshold);
    const auto explicitEnd = std::chrono::steady_clock::now();
    result.explicitMicros = common::elapsedMicros(explicitStart, explicitEnd);
    result.events = runner.lastMergeEvents();
    result.explicitEvents = static_cast<int>(result.events.size());
    addEventStats(result.stats, result.events);
    result.stats.phaseDualRebuilds = static_cast<std::size_t>(runner.lastPhaseDualRebuilds());
    result.stats.phaseCandidateGuardRebuilds =
        static_cast<std::size_t>(runner.lastPhaseCandidateGuardRebuilds());
    result.stats.phaseCandidateGuardCanonicalizationChanges =
        static_cast<std::size_t>(runner.lastPhaseCandidateGuardCanonicalizationChanges());
    result.localUpdateRejectionStopReached =
        runner.lastLocalUpdateRejectionStopReached();

    const auto naiveStart = std::chrono::steady_clock::now();
    const auto naiveOutput = noninc::applyNaiveAttributeMergeThresholdPrefix(
        input,
        adjacency,
        threshold,
        maxEvents,
        attributeKind);
    const auto naiveEnd = std::chrono::steady_clock::now();
    result.naiveMicros = common::elapsedMicros(naiveStart, naiveEnd);
    result.explicitOutput = runner.reconstructImage();
    result.naiveOutput = naiveOutput;
    result.matches = result.explicitOutput->isEqual(result.naiveOutput);
    return result;
}

void printPrefixProbe(const PrefixCompareResult &result) {
    std::cout << "  prefix_probe"
              << " requested_events=" << result.requestedEvents
              << " explicit_events=" << result.explicitEvents
              << " explicit_ms=" << common::microsToMillis(result.explicitMicros)
              << " naive_ms=" << common::microsToMillis(result.naiveMicros)
              << " matches=" << (result.matches ? "true" : "false")
              << " phase_candidate_guard_rebuilds="
              << result.stats.phaseCandidateGuardRebuilds
              << " phase_candidate_guard_canonicalize_changes="
              << result.stats.phaseCandidateGuardCanonicalizationChanges
              << " local_update_rejection_stop="
              << (result.localUpdateRejectionStopReached ? "true" : "false")
              << " stripped_created_nodes=" << result.stats.strippedLocalCreatedNodes
              << " stripped_empty_created_nodes=" << result.stats.strippedEmptyCreatedNodes
              << " stripped_empty_connector_nodes=" << result.stats.strippedEmptyConnectorNodes
              << "\n";
    std::cout.flush();
}

void printDivergentPrefixDetails(const PrefixCompareResult &result) {
    printPrefixProbe(result);
    if (result.events.empty()) {
        return;
    }
    const std::size_t eventIndex = result.events.size() - 1;
    std::cout << "  divergent_prefix_last_event:\n";
    printMergeEvent(result.events[eventIndex], eventIndex);

    const auto &event = result.events[eventIndex];
    if (event.strippedEmptyCreatedNodeDiagnostics.empty()) {
        std::cout << "  divergent_prefix_empty_connectors=none_on_last_event\n";
        return;
    }
    std::cout << "  divergent_prefix_empty_connectors_on_last_event="
              << event.strippedEmptyConnectorNodesAfterCanonicalization
              << "\n";
}

void findAndPrintDivergentPrefix(ImageUInt8Ptr input,
                                 const AdjacencyRelationPtr &adjacency,
                                 double radiusAdj,
                                 float threshold,
                                 noninc::NonIncreasingDualMaintenance dualMaintenance,
                                 noninc::NonIncreasingAttributeKind attributeKind,
                                 const noninc::NonIncreasingStrippedLocalUpdateOptions &strippedOptions,
                                 noninc::NonIncreasingStrippedCommitMode strippedCommitMode,
                                 bool leafFastPathEnabled,
                                 bool phaseCandidateGuardEnabled,
                                 bool strictDirectV1PrecheckEnabled,
                                 bool mcpcNoCopyDirectCommitEnabled,
                                 bool mcpcValidatedDirectCommitEnabled) {
    std::cout << "\nPrefix divergence finder\n";
    std::cout << "  threshold=" << threshold << "\n";
    std::cout.flush();

    int lowMatchingPrefix = 0;
    int highDivergentPrefix = -1;
    PrefixCompareResult highResult;
    int probePrefix = 1;

    std::cout << "  exponential_prefix_search:\n";
    while (highDivergentPrefix < 0) {
        auto probeResult = runPrefixCompare(
            input,
            adjacency,
            radiusAdj,
            threshold,
            dualMaintenance,
            attributeKind,
            strippedOptions,
            strippedCommitMode,
            leafFastPathEnabled,
            phaseCandidateGuardEnabled,
            strictDirectV1PrecheckEnabled,
            mcpcNoCopyDirectCommitEnabled,
            mcpcValidatedDirectCommitEnabled,
            probePrefix);
        printPrefixProbe(probeResult);

        if (!probeResult.matches) {
            highDivergentPrefix = probeResult.explicitEvents;
            highResult = std::move(probeResult);
            break;
        }

        lowMatchingPrefix = probeResult.explicitEvents;
        if (probeResult.explicitEvents < probePrefix) {
            std::cout << "  no divergent prefix found because the threshold ended before the next prefix and matched the naive prefix reference.\n";
            return;
        }
        if (probePrefix > 500000000) {
            std::cout << "  cannot continue exponential prefix search safely after prefix="
                      << probePrefix << "\n";
            return;
        }
        probePrefix *= 2;
    }

    if (highDivergentPrefix <= 0) {
        std::cout << "  cannot bisect divergent prefix because no explicit events were applied.\n";
        return;
    }

    std::cout << "  binary_prefix_search:\n";
    while (lowMatchingPrefix + 1 < highDivergentPrefix) {
        const int mid = lowMatchingPrefix + (highDivergentPrefix - lowMatchingPrefix) / 2;
        auto midResult = runPrefixCompare(
            input,
            adjacency,
            radiusAdj,
            threshold,
            dualMaintenance,
            attributeKind,
            strippedOptions,
            strippedCommitMode,
            leafFastPathEnabled,
            phaseCandidateGuardEnabled,
            strictDirectV1PrecheckEnabled,
            mcpcNoCopyDirectCommitEnabled,
            mcpcValidatedDirectCommitEnabled,
            mid);
        printPrefixProbe(midResult);
        if (midResult.matches) {
            lowMatchingPrefix = mid;
        } else {
            highDivergentPrefix = mid;
            highResult = std::move(midResult);
        }
    }

    std::cout << "\nFirst divergent prefix candidate\n";
    std::cout << "  last_matching_prefix=" << lowMatchingPrefix
              << " first_divergent_prefix=" << highDivergentPrefix
              << "\n";
    if (lowMatchingPrefix > 0) {
        auto previous = runPrefixCompare(
            input,
            adjacency,
            radiusAdj,
            threshold,
            dualMaintenance,
            attributeKind,
            strippedOptions,
            strippedCommitMode,
            leafFastPathEnabled,
            phaseCandidateGuardEnabled,
            strictDirectV1PrecheckEnabled,
            mcpcNoCopyDirectCommitEnabled,
            mcpcValidatedDirectCommitEnabled,
            lowMatchingPrefix);
        std::cout << "  previous_prefix_check:\n";
        printPrefixProbe(previous);
    }
    std::cout << "  divergent_prefix_check:\n";
    printDivergentPrefixDetails(highResult);
}

void dumpPrefixFixtureDiagnostic(ImageUInt8Ptr input,
                                 const AdjacencyRelationPtr &adjacency,
                                 double radiusAdj,
                                 float threshold,
                                 noninc::NonIncreasingDualMaintenance dualMaintenance,
                                 noninc::NonIncreasingAttributeKind attributeKind,
                                 const noninc::NonIncreasingStrippedLocalUpdateOptions &strippedOptions,
                                 noninc::NonIncreasingStrippedCommitMode strippedCommitMode,
                                 bool leafFastPathEnabled,
                                 bool phaseCandidateGuardEnabled,
                                 bool strictDirectV1PrecheckEnabled,
                                 bool mcpcNoCopyDirectCommitEnabled,
                                 bool mcpcValidatedDirectCommitEnabled,
                                 int prefixEvents,
                                 int padding) {
    if (prefixEvents <= 0) {
        throw std::runtime_error("--dump-prefix-fixture requires a positive prefix length");
    }
    if (padding < 0) {
        throw std::runtime_error("--prefix-fixture-padding must be non-negative");
    }

    std::cout << "\nPrefix fixture dump\n"
              << "  threshold=" << threshold
              << " prefix_events=" << prefixEvents
              << " previous_prefix_events=" << (prefixEvents - 1)
              << " padding=" << padding
              << "\n";
    std::cout.flush();

    PrefixCompareResult previous;
    if (prefixEvents > 1) {
        previous = runPrefixCompare(
            input,
            adjacency,
            radiusAdj,
            threshold,
            dualMaintenance,
            attributeKind,
            strippedOptions,
            strippedCommitMode,
            leafFastPathEnabled,
            phaseCandidateGuardEnabled,
            strictDirectV1PrecheckEnabled,
            mcpcNoCopyDirectCommitEnabled,
            mcpcValidatedDirectCommitEnabled,
            prefixEvents - 1);
    } else {
        previous.requestedEvents = 0;
        previous.explicitOutput = input->clone();
        previous.naiveOutput = input->clone();
        previous.matches = true;
    }
    std::cout << "  previous_prefix_check:\n";
    printPrefixProbe(previous);

    auto current = runPrefixCompare(
        input,
        adjacency,
        radiusAdj,
        threshold,
        dualMaintenance,
        attributeKind,
        strippedOptions,
        strippedCommitMode,
        leafFastPathEnabled,
        phaseCandidateGuardEnabled,
        strictDirectV1PrecheckEnabled,
        mcpcNoCopyDirectCommitEnabled,
        mcpcValidatedDirectCommitEnabled,
        prefixEvents);
    std::cout << "  current_prefix_check:\n";
    printPrefixProbe(current);

    if (current.events.empty()) {
        std::cout << "  prefix_fixture_last_event=none\n";
        return;
    }

    const std::size_t eventIndex = current.events.size() - 1;
    const auto &event = current.events[eventIndex];
    std::cout << "  prefix_fixture_last_event:\n";
    printMergeEvent(event, eventIndex);

    const auto previousDiff = summarizeImageDiff(previous.explicitOutput, previous.naiveOutput);
    const auto currentDiff = summarizeImageDiff(current.explicitOutput, current.naiveOutput);
    std::cout << "  previous_prefix_differs=" << (previousDiff.differs ? "true" : "false")
              << " differing_pixels=" << previousDiff.differingPixels << "\n";
    printWindowSummary("previous_diff", previousDiff.bbox);
    std::cout << "  current_prefix_differs=" << (currentDiff.differs ? "true" : "false")
              << " differing_pixels=" << currentDiff.differingPixels << "\n";
    printWindowSummary("current_diff", currentDiff.bbox);
    printPixelCoordPreview("current_diff_preview", currentDiff.previewPixels, input->getNumCols());
    printPixelValueTriples(
        "current_diff_values",
        previous.explicitOutput,
        current.explicitOutput,
        current.naiveOutput,
        currentDiff.previewPixels);

    const auto rawEventWindow = eventWindow(event);
    const auto paddedEventWindow = expandedWindow(
        rawEventWindow,
        padding,
        input->getNumRows(),
        input->getNumCols());
    const auto paddedDiffWindow = expandedWindow(
        currentDiff.bbox,
        padding,
        input->getNumRows(),
        input->getNumCols());
    const auto combinedWindow = mergeWindows(
        paddedEventWindow,
        paddedDiffWindow,
        input->getNumRows(),
        input->getNumCols());

    printWindowSummary("event", rawEventWindow);
    printWindowSummary("event_padded", paddedEventWindow);
    printWindowSummary("diff_padded", paddedDiffWindow);
    printWindowSummary("combined", combinedWindow);
    printPixelCoordPreview("affected_pixels_preview", event.affectedPixelsPreview, input->getNumCols());

    std::cout << "\n  Event-window crops\n";
    printImageWindowValues("previous_explicit", previous.explicitOutput, paddedEventWindow);
    printImageWindowValues("current_explicit", current.explicitOutput, paddedEventWindow);
    printImageWindowValues("current_naive", current.naiveOutput, paddedEventWindow);

    std::cout << "\n  Diff-window crops\n";
    printImageWindowValues("previous_explicit", previous.explicitOutput, paddedDiffWindow);
    printImageWindowValues("current_explicit", current.explicitOutput, paddedDiffWindow);
    printImageWindowValues("current_naive", current.naiveOutput, paddedDiffWindow);

    if (combinedWindow.valid && combinedWindow.area() <= 400) {
        std::cout << "\n  Combined fixture crop\n";
        printImageWindowValues("previous_explicit", previous.explicitOutput, combinedWindow);
        printImageWindowValues("current_explicit", current.explicitOutput, combinedWindow);
        printImageWindowValues("current_naive", current.naiveOutput, combinedWindow);
    } else if (combinedWindow.valid) {
        std::cout << "\n  Combined fixture crop omitted"
                  << " area=" << combinedWindow.area()
                  << " max_printed_area=400\n";
    }

    printDiffPointCrops(
        previous.explicitOutput,
        current.explicitOutput,
        current.naiveOutput,
        currentDiff.previewPixels,
        padding);
}

void printThresholdEventStats(float threshold,
                              const MergeEventStats &stats,
                              long long explicitTotalUs,
                              long long naiveTotalUs,
                              bool matches) {
    std::cout << "  threshold_summary"
              << " threshold=" << threshold
              << " explicit_ms=" << common::microsToMillis(explicitTotalUs)
              << " naive_ms=" << common::microsToMillis(naiveTotalUs)
              << " matches=" << (matches ? "true" : "false")
              << " events=" << stats.events
              << " max_to_min=" << stats.maxToMinEvents
              << " min_to_max=" << stats.minToMaxEvents
              << " dual_rebuilds=" << stats.dualRebuildEvents
              << " phase_dual_rebuilds=" << stats.phaseDualRebuilds
              << " phase_candidate_guard_rebuilds=" << stats.phaseCandidateGuardRebuilds
              << " phase_candidate_guard_canonicalize_changes="
              << stats.phaseCandidateGuardCanonicalizationChanges
              << " multi_component_events=" << stats.multiComponentEvents
              << " child_events=" << stats.childEvents
              << " stripped_attempts=" << stats.strippedLocalAttempts
              << " stripped_accepted=" << stats.strippedLocalAccepted
              << " stripped_applied=" << stats.strippedLocalApplied
              << " stripped_direct_commits=" << stats.strippedDirectCommits
              << " stripped_plan_audits=" << stats.strippedLocalPlanAuditChecks
              << " stripped_target_checks=" << stats.strippedLocalTargetChecks
              << " stripped_target_matches=" << stats.strippedLocalTargetMatches
              << " stripped_plan_oracle_free=" << stats.strippedLocalOracleFreePlans
              << " stripped_materialization_audits="
              << stats.strippedMaterializationAuditChecks
              << " stripped_materialization_oracle_free="
              << stats.strippedLocalOracleFreeMaterializations
              << " stripped_final_checks=" << stats.strippedFinalReconstructionChecks
              << " stripped_fallbacks=" << stats.strippedLocalFallbacks
              << " stripped_rebuild_fallbacks=" << stats.strippedLocalRebuildFallbacks
              << " stripped_phase_rebuild_fallbacks=" << stats.strippedLocalPhaseRebuildFallbacks
              << " stripped_phase_rebuild_skips=" << stats.strippedLocalPhaseRebuildSkips
              << " partial_f_attempts=" << stats.partialFAttempts
              << " partial_f_accepted=" << stats.partialFAccepted
              << " partial_f_rejected=" << stats.partialFRejected
              << " partial_f_mutation_attempted=" << stats.partialFMutationAttempted
              << " partial_f_mcpc_live_accepted=" << stats.partialFMcpcLiveAccepted
              << " partial_f_mcpc_hidden_support_redundant_live_accepted="
              << stats.partialFMcpcHiddenSupportRedundantLiveAccepted
              << " partial_f_mcpc_parent_link_bridge_live_accepted="
              << stats.partialFMcpcParentLinkBridgeLiveAccepted
              << " partial_f_mcpc_stable_child_bridge_live_accepted="
              << stats.partialFMcpcStableChildBridgeLiveAccepted
              << " partial_f_mcpc_child_link_bridge_live_accepted="
              << stats.partialFMcpcChildLinkBridgeLiveAccepted
              << " partial_f_mcpc_empty_owner_release_live_accepted="
              << stats.partialFMcpcEmptyOwnerReleaseLiveAccepted
              << " partial_f_mcpc_stable_child_move_live_accepted="
              << stats.partialFMcpcStableChildMoveLiveAccepted
              << partialFRejectStatsText(stats)
              << " merge_node_attempts=" << stats.mergeNodeLocalAttempts
              << " merge_node_accepted=" << stats.mergeNodeLocalAccepted
              << " merge_node_applied=" << stats.mergeNodeLocalApplied
              << " merge_node_cert_checks=" << stats.mergeNodeLocalCertificationChecks
              << " merge_node_cert_passed=" << stats.mergeNodeLocalCertificationPassed
              << " merge_node_direct_commits=" << stats.mergeNodeLocalDirectCommits
              << " merge_node_target_creation_allowed_events="
              << stats.mergeNodeLocalTargetCreationAllowedEvents
              << " merge_node_target_creation_required_events="
              << stats.mergeNodeLocalTargetCreationRequiredEvents
              << " merge_node_pre_dual_matches=" << stats.mergeNodeLocalPreDualMatches
              << " merge_node_dual_matches=" << stats.mergeNodeLocalDualMatches
              << " merge_node_phase_rebuild_fallbacks="
              << stats.mergeNodeLocalPhaseRebuildFallbacks
              << " merge_node_phase_rebuild_skips="
              << stats.mergeNodeLocalPhaseRebuildSkips
              << " merge_node_rejected_components="
              << stats.mergeNodeLocalRejectedComponents
              << " merge_node_problem1_preflight_rejected_components="
              << stats.mergeNodeLocalProblem1PreflightRejectedComponents
              << " merge_node_problem1_preflight_rejected_by_published="
              << stats.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem
              << " merge_node_problem1_preflight_rejected_by_updater="
              << stats.mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation
              << " merge_node_partial_owner_components="
              << stats.mergeNodeLocalPartialOwnerComponents
              << " merge_node_raw_missing_target_b_components="
              << stats.mergeNodeLocalRawMissingTargetComponents
              << " merge_node_certified_missing_target_b_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_partial_owner_and_missing_target_b_components="
              << stats.mergeNodeLocalPartialOwnerAndMissingTargetComponents
              << " merge_node_missing_target_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_pixels=" << stats.mergeNodeLocalPixels
              << " merge_node_copy_ms=" << common::microsToMillis(stats.mergeNodeLocalCopyMicros)
              << " merge_node_cert_ms="
              << common::microsToMillis(stats.mergeNodeLocalCertificationMicros)
              << " merge_node_component_collect_ms="
              << common::microsToMillis(stats.mergeNodeComponentCollectMicros)
              << " merge_node_update_ms=" << common::microsToMillis(stats.mergeNodeLocalUpdateMicros)
              << " merge_node_primal_merge_ms="
              << common::microsToMillis(stats.mergeNodePrimalMergeMicros)
              << " merge_node_empty_release_ms="
              << common::microsToMillis(stats.mergeNodeEmptyReleaseMicros)
              << " stripped_created_nodes=" << stats.strippedLocalCreatedNodes
              << " stripped_recomposition_nodes=" << stats.strippedLocalRecompositionNodes
              << " stripped_local_pixels=" << stats.strippedLocalPixels
              << " max_stripped_local_pixels=" << stats.maxStrippedLocalPixels
              << " stripped_empty_created_nodes=" << stats.strippedEmptyCreatedNodes
              << " stripped_empty_connector_nodes=" << stats.strippedEmptyConnectorNodes
              << " stripped_released_owners=" << stats.strippedLocalReleasedOwners
              << " stripped_copy_ms=" << common::microsToMillis(stats.strippedLocalCopyMicros)
              << " stripped_update_ms=" << common::microsToMillis(stats.strippedLocalUpdateMicros)
              << " stripped_fallback_ms=" << common::microsToMillis(stats.strippedOracleFallbackMicros)
              << " stripped_rebuild_fallback_ms=" << common::microsToMillis(stats.strippedRebuildFallbackMicros)
              << " stripped_input_ms=" << common::microsToMillis(stats.strippedInputExtractionMicros)
              << " stripped_partition_ms=" << common::microsToMillis(stats.strippedPartitionPlanMicros)
              << " stripped_recomposition_ms=" << common::microsToMillis(stats.strippedRecompositionPlanMicros)
              << " stripped_materialization_ms="
              << common::microsToMillis(stats.strippedMaterializationMicros)
              << " stripped_mat_canonicalize_ms="
              << common::microsToMillis(stats.strippedMaterializationCanonicalizationMicros);
    printMaxCanonicalizationEventFields(stats);
    std::cout << " stripped_mat_release_ms="
              << common::microsToMillis(stats.strippedMaterializationReleaseMicros)
              << " stripped_mat_target_check_ms="
              << common::microsToMillis(stats.strippedMaterializationTargetCheckMicros)
              << " stripped_primal_release_ms="
              << common::microsToMillis(stats.strippedPrimalMergeReleaseMicros);
    printMaxPrimalMergeEventFields(stats);
    std::cout << " stripped_final_check_ms="
              << common::microsToMillis(stats.strippedFinalReconstructionCheckMicros)
              << " oracle_local_events=" << stats.oracleLocalEvents
              << " oracle_local_target_matches=" << stats.oracleLocalTargetMatches
              << " oracle_local_structure_matches=" << stats.oracleLocalStructureMatches
              << " oracle_local_structure_mismatches=" << stats.oracleLocalStructureMismatches
              << " oracle_created_nodes=" << stats.oracleLocalCreatedNodes
              << " oracle_recomposition_nodes=" << stats.oracleLocalRecompositionNodes
              << " oracle_released_owners=" << stats.oracleLocalReleasedOwners
              << " total_C=" << stats.totalComponentPixels
              << " max_components=" << stats.maxComponents
              << " max_C=" << stats.maxC
              << "\n";
}

void printExplicitThresholdEventStats(float threshold,
                                      const MergeEventStats &stats,
                                      long long explicitTotalUs) {
    std::cout << "  threshold_summary"
              << " threshold=" << threshold
              << " explicit_ms=" << common::microsToMillis(explicitTotalUs)
              << " events=" << stats.events
              << " max_to_min=" << stats.maxToMinEvents
              << " min_to_max=" << stats.minToMaxEvents
              << " dual_rebuilds=" << stats.dualRebuildEvents
              << " phase_dual_rebuilds=" << stats.phaseDualRebuilds
              << " phase_candidate_guard_rebuilds=" << stats.phaseCandidateGuardRebuilds
              << " phase_candidate_guard_canonicalize_changes="
              << stats.phaseCandidateGuardCanonicalizationChanges
              << " multi_component_events=" << stats.multiComponentEvents
              << " child_events=" << stats.childEvents
              << " stripped_attempts=" << stats.strippedLocalAttempts
              << " stripped_accepted=" << stats.strippedLocalAccepted
              << " stripped_applied=" << stats.strippedLocalApplied
              << " stripped_direct_commits=" << stats.strippedDirectCommits
              << " stripped_plan_audits=" << stats.strippedLocalPlanAuditChecks
              << " stripped_target_checks=" << stats.strippedLocalTargetChecks
              << " stripped_target_matches=" << stats.strippedLocalTargetMatches
              << " stripped_plan_oracle_free=" << stats.strippedLocalOracleFreePlans
              << " stripped_materialization_audits="
              << stats.strippedMaterializationAuditChecks
              << " stripped_materialization_oracle_free="
              << stats.strippedLocalOracleFreeMaterializations
              << " stripped_final_checks=" << stats.strippedFinalReconstructionChecks
              << " stripped_fallbacks=" << stats.strippedLocalFallbacks
              << " stripped_rebuild_fallbacks=" << stats.strippedLocalRebuildFallbacks
              << " stripped_phase_rebuild_fallbacks=" << stats.strippedLocalPhaseRebuildFallbacks
              << " stripped_phase_rebuild_skips=" << stats.strippedLocalPhaseRebuildSkips
              << " partial_f_attempts=" << stats.partialFAttempts
              << " partial_f_accepted=" << stats.partialFAccepted
              << " partial_f_rejected=" << stats.partialFRejected
              << " partial_f_mutation_attempted=" << stats.partialFMutationAttempted
              << " partial_f_mcpc_live_accepted=" << stats.partialFMcpcLiveAccepted
              << " partial_f_mcpc_hidden_support_redundant_live_accepted="
              << stats.partialFMcpcHiddenSupportRedundantLiveAccepted
              << " partial_f_mcpc_parent_link_bridge_live_accepted="
              << stats.partialFMcpcParentLinkBridgeLiveAccepted
              << " partial_f_mcpc_stable_child_bridge_live_accepted="
              << stats.partialFMcpcStableChildBridgeLiveAccepted
              << " partial_f_mcpc_child_link_bridge_live_accepted="
              << stats.partialFMcpcChildLinkBridgeLiveAccepted
              << " partial_f_mcpc_empty_owner_release_live_accepted="
              << stats.partialFMcpcEmptyOwnerReleaseLiveAccepted
              << " partial_f_mcpc_stable_child_move_live_accepted="
              << stats.partialFMcpcStableChildMoveLiveAccepted
              << partialFRejectStatsText(stats)
              << " merge_node_attempts=" << stats.mergeNodeLocalAttempts
              << " merge_node_accepted=" << stats.mergeNodeLocalAccepted
              << " merge_node_applied=" << stats.mergeNodeLocalApplied
              << " merge_node_cert_checks=" << stats.mergeNodeLocalCertificationChecks
              << " merge_node_cert_passed=" << stats.mergeNodeLocalCertificationPassed
              << " merge_node_direct_commits=" << stats.mergeNodeLocalDirectCommits
              << " merge_node_target_creation_allowed_events="
              << stats.mergeNodeLocalTargetCreationAllowedEvents
              << " merge_node_target_creation_required_events="
              << stats.mergeNodeLocalTargetCreationRequiredEvents
              << " merge_node_pre_dual_matches=" << stats.mergeNodeLocalPreDualMatches
              << " merge_node_dual_matches=" << stats.mergeNodeLocalDualMatches
              << " merge_node_phase_rebuild_fallbacks="
              << stats.mergeNodeLocalPhaseRebuildFallbacks
              << " merge_node_phase_rebuild_skips="
              << stats.mergeNodeLocalPhaseRebuildSkips
              << " merge_node_rejected_components="
              << stats.mergeNodeLocalRejectedComponents
              << " merge_node_problem1_preflight_rejected_components="
              << stats.mergeNodeLocalProblem1PreflightRejectedComponents
              << " merge_node_problem1_preflight_rejected_by_published="
              << stats.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem
              << " merge_node_problem1_preflight_rejected_by_updater="
              << stats.mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation
              << " merge_node_partial_owner_components="
              << stats.mergeNodeLocalPartialOwnerComponents
              << " merge_node_raw_missing_target_b_components="
              << stats.mergeNodeLocalRawMissingTargetComponents
              << " merge_node_certified_missing_target_b_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_partial_owner_and_missing_target_b_components="
              << stats.mergeNodeLocalPartialOwnerAndMissingTargetComponents
              << " merge_node_missing_target_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_pixels=" << stats.mergeNodeLocalPixels
              << " merge_node_copy_ms=" << common::microsToMillis(stats.mergeNodeLocalCopyMicros)
              << " merge_node_cert_ms="
              << common::microsToMillis(stats.mergeNodeLocalCertificationMicros)
              << " merge_node_component_collect_ms="
              << common::microsToMillis(stats.mergeNodeComponentCollectMicros)
              << " merge_node_update_ms=" << common::microsToMillis(stats.mergeNodeLocalUpdateMicros)
              << " merge_node_primal_merge_ms="
              << common::microsToMillis(stats.mergeNodePrimalMergeMicros)
              << " merge_node_empty_release_ms="
              << common::microsToMillis(stats.mergeNodeEmptyReleaseMicros)
              << " stripped_created_nodes=" << stats.strippedLocalCreatedNodes
              << " stripped_recomposition_nodes=" << stats.strippedLocalRecompositionNodes
              << " stripped_local_pixels=" << stats.strippedLocalPixels
              << " max_stripped_local_pixels=" << stats.maxStrippedLocalPixels
              << " stripped_empty_created_nodes=" << stats.strippedEmptyCreatedNodes
              << " stripped_empty_connector_nodes=" << stats.strippedEmptyConnectorNodes
              << " stripped_released_owners=" << stats.strippedLocalReleasedOwners
              << " stripped_copy_ms=" << common::microsToMillis(stats.strippedLocalCopyMicros)
              << " stripped_update_ms=" << common::microsToMillis(stats.strippedLocalUpdateMicros)
              << " stripped_fallback_ms=" << common::microsToMillis(stats.strippedOracleFallbackMicros)
              << " stripped_rebuild_fallback_ms=" << common::microsToMillis(stats.strippedRebuildFallbackMicros)
              << " stripped_input_ms=" << common::microsToMillis(stats.strippedInputExtractionMicros)
              << " stripped_partition_ms=" << common::microsToMillis(stats.strippedPartitionPlanMicros)
              << " stripped_recomposition_ms=" << common::microsToMillis(stats.strippedRecompositionPlanMicros)
              << " stripped_materialization_ms="
              << common::microsToMillis(stats.strippedMaterializationMicros)
              << " stripped_mat_canonicalize_ms="
              << common::microsToMillis(stats.strippedMaterializationCanonicalizationMicros);
    printMaxCanonicalizationEventFields(stats);
    std::cout << " stripped_mat_release_ms="
              << common::microsToMillis(stats.strippedMaterializationReleaseMicros)
              << " stripped_mat_target_check_ms="
              << common::microsToMillis(stats.strippedMaterializationTargetCheckMicros)
              << " stripped_primal_release_ms="
              << common::microsToMillis(stats.strippedPrimalMergeReleaseMicros);
    printMaxPrimalMergeEventFields(stats);
    std::cout << " stripped_final_check_ms="
              << common::microsToMillis(stats.strippedFinalReconstructionCheckMicros)
              << " oracle_local_events=" << stats.oracleLocalEvents
              << " oracle_local_target_matches=" << stats.oracleLocalTargetMatches
              << " oracle_local_structure_matches=" << stats.oracleLocalStructureMatches
              << " oracle_local_structure_mismatches=" << stats.oracleLocalStructureMismatches
              << " oracle_created_nodes=" << stats.oracleLocalCreatedNodes
              << " oracle_recomposition_nodes=" << stats.oracleLocalRecompositionNodes
              << " oracle_released_owners=" << stats.oracleLocalReleasedOwners
              << " total_C=" << stats.totalComponentPixels
              << " max_components=" << stats.maxComponents
              << " max_C=" << stats.maxC
              << "\n";
}

void printAggregateEventStats(const MergeEventStats &stats) {
    std::cout << "\nAggregate merge-event diagnostics\n";
    std::cout << "  events=" << stats.events
              << " max_to_min=" << stats.maxToMinEvents
              << " min_to_max=" << stats.minToMaxEvents
              << " dual_rebuilds=" << stats.dualRebuildEvents
              << " phase_dual_rebuilds=" << stats.phaseDualRebuilds
              << " phase_candidate_guard_rebuilds=" << stats.phaseCandidateGuardRebuilds
              << " phase_candidate_guard_canonicalize_changes="
              << stats.phaseCandidateGuardCanonicalizationChanges
              << " multi_component_events=" << stats.multiComponentEvents
              << " child_events=" << stats.childEvents
              << " stripped_attempts=" << stats.strippedLocalAttempts
              << " stripped_accepted=" << stats.strippedLocalAccepted
              << " stripped_applied=" << stats.strippedLocalApplied
              << " stripped_direct_commits=" << stats.strippedDirectCommits
              << " stripped_plan_audits=" << stats.strippedLocalPlanAuditChecks
              << " stripped_target_checks=" << stats.strippedLocalTargetChecks
              << " stripped_target_matches=" << stats.strippedLocalTargetMatches
              << " stripped_plan_oracle_free=" << stats.strippedLocalOracleFreePlans
              << " stripped_materialization_audits="
              << stats.strippedMaterializationAuditChecks
              << " stripped_materialization_oracle_free="
              << stats.strippedLocalOracleFreeMaterializations
              << " stripped_final_checks=" << stats.strippedFinalReconstructionChecks
              << " stripped_fallbacks=" << stats.strippedLocalFallbacks
              << " stripped_rebuild_fallbacks=" << stats.strippedLocalRebuildFallbacks
              << " stripped_phase_rebuild_fallbacks=" << stats.strippedLocalPhaseRebuildFallbacks
              << " stripped_phase_rebuild_skips=" << stats.strippedLocalPhaseRebuildSkips
              << " partial_f_attempts=" << stats.partialFAttempts
              << " partial_f_accepted=" << stats.partialFAccepted
              << " partial_f_rejected=" << stats.partialFRejected
              << " partial_f_mutation_attempted=" << stats.partialFMutationAttempted
              << " partial_f_mcpc_live_accepted=" << stats.partialFMcpcLiveAccepted
              << " partial_f_mcpc_hidden_support_redundant_live_accepted="
              << stats.partialFMcpcHiddenSupportRedundantLiveAccepted
              << " partial_f_mcpc_parent_link_bridge_live_accepted="
              << stats.partialFMcpcParentLinkBridgeLiveAccepted
              << " partial_f_mcpc_stable_child_bridge_live_accepted="
              << stats.partialFMcpcStableChildBridgeLiveAccepted
              << " partial_f_mcpc_child_link_bridge_live_accepted="
              << stats.partialFMcpcChildLinkBridgeLiveAccepted
              << " partial_f_mcpc_empty_owner_release_live_accepted="
              << stats.partialFMcpcEmptyOwnerReleaseLiveAccepted
              << " partial_f_mcpc_stable_child_move_live_accepted="
              << stats.partialFMcpcStableChildMoveLiveAccepted
              << partialFRejectStatsText(stats)
              << " merge_node_attempts=" << stats.mergeNodeLocalAttempts
              << " merge_node_accepted=" << stats.mergeNodeLocalAccepted
              << " merge_node_applied=" << stats.mergeNodeLocalApplied
              << " merge_node_cert_checks=" << stats.mergeNodeLocalCertificationChecks
              << " merge_node_cert_passed=" << stats.mergeNodeLocalCertificationPassed
              << " merge_node_direct_commits=" << stats.mergeNodeLocalDirectCommits
              << " merge_node_target_creation_allowed_events="
              << stats.mergeNodeLocalTargetCreationAllowedEvents
              << " merge_node_target_creation_required_events="
              << stats.mergeNodeLocalTargetCreationRequiredEvents
              << " merge_node_pre_dual_matches=" << stats.mergeNodeLocalPreDualMatches
              << " merge_node_dual_matches=" << stats.mergeNodeLocalDualMatches
              << " merge_node_phase_rebuild_fallbacks="
              << stats.mergeNodeLocalPhaseRebuildFallbacks
              << " merge_node_phase_rebuild_skips="
              << stats.mergeNodeLocalPhaseRebuildSkips
              << " merge_node_rejected_components="
              << stats.mergeNodeLocalRejectedComponents
              << " merge_node_partial_owner_components="
              << stats.mergeNodeLocalPartialOwnerComponents
              << " merge_node_raw_missing_target_b_components="
              << stats.mergeNodeLocalRawMissingTargetComponents
              << " merge_node_certified_missing_target_b_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_partial_owner_and_missing_target_b_components="
              << stats.mergeNodeLocalPartialOwnerAndMissingTargetComponents
              << " merge_node_missing_target_components="
              << stats.mergeNodeLocalMissingTargetComponents
              << " merge_node_pixels=" << stats.mergeNodeLocalPixels
              << " merge_node_copy_ms=" << common::microsToMillis(stats.mergeNodeLocalCopyMicros)
              << " merge_node_cert_ms="
              << common::microsToMillis(stats.mergeNodeLocalCertificationMicros)
              << " merge_node_component_collect_ms="
              << common::microsToMillis(stats.mergeNodeComponentCollectMicros)
              << " merge_node_update_ms=" << common::microsToMillis(stats.mergeNodeLocalUpdateMicros)
              << " merge_node_primal_merge_ms="
              << common::microsToMillis(stats.mergeNodePrimalMergeMicros)
              << " merge_node_empty_release_ms="
              << common::microsToMillis(stats.mergeNodeEmptyReleaseMicros)
              << " stripped_created_nodes=" << stats.strippedLocalCreatedNodes
              << " stripped_recomposition_nodes=" << stats.strippedLocalRecompositionNodes
              << " stripped_local_pixels=" << stats.strippedLocalPixels
              << " max_stripped_local_pixels=" << stats.maxStrippedLocalPixels
              << " stripped_empty_created_nodes=" << stats.strippedEmptyCreatedNodes
              << " stripped_empty_connector_nodes=" << stats.strippedEmptyConnectorNodes
              << " stripped_released_owners=" << stats.strippedLocalReleasedOwners
              << " stripped_copy_ms=" << common::microsToMillis(stats.strippedLocalCopyMicros)
              << " stripped_update_ms=" << common::microsToMillis(stats.strippedLocalUpdateMicros)
              << " stripped_fallback_ms=" << common::microsToMillis(stats.strippedOracleFallbackMicros)
              << " stripped_rebuild_fallback_ms=" << common::microsToMillis(stats.strippedRebuildFallbackMicros)
              << " stripped_input_ms=" << common::microsToMillis(stats.strippedInputExtractionMicros)
              << " stripped_partition_ms=" << common::microsToMillis(stats.strippedPartitionPlanMicros)
              << " stripped_recomposition_ms=" << common::microsToMillis(stats.strippedRecompositionPlanMicros)
              << " stripped_materialization_ms="
              << common::microsToMillis(stats.strippedMaterializationMicros)
              << " stripped_mat_canonicalize_ms="
              << common::microsToMillis(stats.strippedMaterializationCanonicalizationMicros);
    printMaxCanonicalizationEventFields(stats);
    std::cout << " stripped_mat_release_ms="
              << common::microsToMillis(stats.strippedMaterializationReleaseMicros)
              << " stripped_mat_target_check_ms="
              << common::microsToMillis(stats.strippedMaterializationTargetCheckMicros)
              << " stripped_primal_release_ms="
              << common::microsToMillis(stats.strippedPrimalMergeReleaseMicros);
    printMaxPrimalMergeEventFields(stats);
    std::cout << " stripped_final_check_ms="
              << common::microsToMillis(stats.strippedFinalReconstructionCheckMicros)
              << " oracle_local_events=" << stats.oracleLocalEvents
              << " oracle_local_target_matches=" << stats.oracleLocalTargetMatches
              << " oracle_local_structure_matches=" << stats.oracleLocalStructureMatches
              << " oracle_local_structure_mismatches=" << stats.oracleLocalStructureMismatches
              << " oracle_created_nodes=" << stats.oracleLocalCreatedNodes
              << " oracle_recomposition_nodes=" << stats.oracleLocalRecompositionNodes
              << " oracle_released_owners=" << stats.oracleLocalReleasedOwners
              << " total_C=" << stats.totalComponentPixels
              << " max_components=" << stats.maxComponents
              << " max_C=" << stats.maxC
              << "\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        if (argc < 4) {
            printUsage(argv[0]);
            return 1;
        }

        int argi = 1;
        Mode mode = Mode::ExplicitUpdate;
        auto dualMaintenance = noninc::NonIncreasingDualMaintenance::PartialF;
        auto attributeKind = noninc::NonIncreasingAttributeKind::GrayVariance;
        std::string strippedValidationMode = "safe";
        auto strippedLocalUpdateOptions = noninc::NonIncreasingStrippedLocalUpdateOptions{};
        std::string strippedCommitModeName = "safe-copy";
        auto strippedCommitMode = noninc::NonIncreasingStrippedCommitMode::SafeCopy;
        double radioAdj = 1.0;
        bool noOutput = false;
        bool diagnoseFirstDivergence = false;
        bool stopOnFirstPhaseFallback = false;
        bool findDivergentPrefix = false;
        bool dumpPrefixFixture = false;
        bool diagnoseEvents = false;
        bool diagnosePhaseOracle = false;
        bool phaseCandidateGuard = false;
        bool dumpPhaseOracleFixture = false;
        bool dumpSmallestPhaseOracleFixture = false;
        bool dumpIncidentClosureFixture = false;
        bool dumpFirstLocalOracle = false;
        bool dumpTargetLocalOracle = false;
        bool strictDirectV1Precheck = false;
        bool mcpcNoCopyDirectCommit = false;
        bool mcpcValidatedDirectCommit = false;
        bool leafFastPathEnabled = true;
        int maxMergeEvents = 0;
        long long timeBudgetMs = 0;
        int prefixFixtureEvents = 0;
        int prefixFixturePadding = 2;
        float targetLocalOracleThreshold = 0.0f;
        noninc::NonIncreasingPhase targetLocalOraclePhase = noninc::NonIncreasingPhase::MaxToMin;
        NodeId targetLocalOracleNodeId = InvalidNode;
        while (argi < argc && std::string(argv[argi]).rfind("--", 0) == 0) {
            const std::string option = argv[argi];
            if (option == "--mode") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                const std::string value = argv[argi + 1];
                if (value == "explicit") {
                    mode = Mode::ExplicitUpdate;
                } else if (value == "naive") {
                    mode = Mode::Naive;
                } else if (value == "compare") {
                    mode = Mode::Compare;
                } else {
                    throw std::runtime_error("Invalid mode: " + value);
                }
                argi += 2;
                continue;
            }
            if (option == "--attribute") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                attributeKind = parseAttributeKind(argv[argi + 1]);
                argi += 2;
                continue;
            }
            if (option == "--dual-maintenance") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                const std::string value = argv[argi + 1];
                if (value == "partial-f" || value == "f") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::PartialF;
                } else if (value == "problem1-restricted" ||
                    value == "p1-restricted" ||
                    value == "problem1-close-recompose" ||
                    value == "problem1-close" ||
                    value == "explicit-only") {
                    dualMaintenance = (value == "problem1-close-recompose" ||
                                      value == "problem1-close")
                                         ? noninc::NonIncreasingDualMaintenance::Problem1CloseRecompose
                                         : noninc::NonIncreasingDualMaintenance::ExplicitOnly;
                } else if (value == "materialize-explicit") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::MaterializeThenExplicit;
                } else if (value == "oracle-local") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::OracleGuidedLocalRecomposition;
                } else if (value == "stripped-oracle-local") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::StrippedThenOracleLocalRecomposition;
                } else if (value == "stripped-rebuild-local") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition;
                } else if (value == "stripped-phase-rebuild-local") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition;
                } else if (value == "merge-node-phase-rebuild") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::MergeNodeUpdateThenPhaseRebuild;
                } else if (value == "merge-node-certified-phase-rebuild") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild;
                } else if (value == "merge-node-certified-create-b-phase-rebuild") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild;
                } else if (value == "merge-node-certified-create-b-mcpc-phase-rebuild") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild;
                } else if (value == "merge-node-clean-unchecked") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::MergeNodeCleanUnchecked;
                } else if (value == "phase-rebuild") {
                    dualMaintenance = noninc::NonIncreasingDualMaintenance::RebuildAfterPhase;
                } else {
                    throw std::runtime_error("Invalid dual maintenance mode: " + value);
                }
                argi += 2;
                continue;
            }
            if (option == "--stripped-validation") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                strippedValidationMode = argv[argi + 1];
                strippedLocalUpdateOptions = strippedValidationOptions(strippedValidationMode);
                argi += 2;
                continue;
            }
            if (option == "--stripped-commit") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                strippedCommitModeName = argv[argi + 1];
                strippedCommitMode = strippedCommitModeFromString(strippedCommitModeName);
                argi += 2;
                continue;
            }
            if (option == "--strict-direct-v1-precheck") {
                strictDirectV1Precheck = true;
                ++argi;
                continue;
            }
            if (option == "--mcpc-no-copy-direct-commit") {
                mcpcNoCopyDirectCommit = true;
                ++argi;
                continue;
            }
            if (option == "--mcpc-validated-direct-commit") {
                mcpcValidatedDirectCommit = true;
                ++argi;
                continue;
            }
            if (option == "--diagnose-first-divergence") {
                diagnoseFirstDivergence = true;
                ++argi;
                continue;
            }
            if (option == "--stop-on-first-phase-fallback") {
                stopOnFirstPhaseFallback = true;
                diagnoseEvents = true;
                ++argi;
                continue;
            }
            if (option == "--find-divergent-prefix") {
                findDivergentPrefix = true;
                ++argi;
                continue;
            }
            if (option == "--dump-prefix-fixture") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char *endPtr = nullptr;
                const long value = std::strtol(argv[argi + 1], &endPtr, 10);
                if (endPtr == argv[argi + 1] || *endPtr != '\0' || value <= 0) {
                    throw std::runtime_error("Invalid --dump-prefix-fixture value. Expected a positive integer.");
                }
                dumpPrefixFixture = true;
                prefixFixtureEvents = static_cast<int>(value);
                argi += 2;
                continue;
            }
            if (option == "--prefix-fixture-padding") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char *endPtr = nullptr;
                const long value = std::strtol(argv[argi + 1], &endPtr, 10);
                if (endPtr == argv[argi + 1] || *endPtr != '\0' || value < 0) {
                    throw std::runtime_error("Invalid --prefix-fixture-padding value. Expected a non-negative integer.");
                }
                prefixFixturePadding = static_cast<int>(value);
                argi += 2;
                continue;
            }
            if (option == "--diagnose-events") {
                diagnoseEvents = true;
                ++argi;
                continue;
            }
            if (option == "--diagnose-phase-oracle") {
                diagnosePhaseOracle = true;
                ++argi;
                continue;
            }
            if (option == "--phase-candidate-guard") {
                phaseCandidateGuard = true;
                ++argi;
                continue;
            }
            if (option == "--dump-phase-oracle-fixture") {
                dumpPhaseOracleFixture = true;
                diagnosePhaseOracle = true;
                ++argi;
                continue;
            }
            if (option == "--dump-smallest-phase-oracle-fixture") {
                dumpSmallestPhaseOracleFixture = true;
                diagnosePhaseOracle = true;
                ++argi;
                continue;
            }
            if (option == "--dump-incident-closure-fixture") {
                dumpIncidentClosureFixture = true;
                diagnosePhaseOracle = true;
                ++argi;
                continue;
            }
            if (option == "--dump-first-local-oracle") {
                dumpFirstLocalOracle = true;
                diagnoseFirstDivergence = true;
                ++argi;
                continue;
            }
            if (option == "--dump-local-oracle") {
                if (argi + 3 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char *endPtr = nullptr;
                targetLocalOracleThreshold = std::strtof(argv[argi + 1], &endPtr);
                if (endPtr == argv[argi + 1] || *endPtr != '\0') {
                    throw std::runtime_error("Invalid local oracle threshold: " + std::string(argv[argi + 1]));
                }
                targetLocalOraclePhase = parsePhase(argv[argi + 2]);
                targetLocalOracleNodeId = parseNodeId(argv[argi + 3]);
                dumpTargetLocalOracle = true;
                argi += 4;
                continue;
            }
            if (option == "--max-merge-events") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char *endPtr = nullptr;
                const long value = std::strtol(argv[argi + 1], &endPtr, 10);
                if (endPtr == argv[argi + 1] || *endPtr != '\0' || value < 0) {
                    throw std::runtime_error("Invalid --max-merge-events value. Expected a non-negative integer.");
                }
                maxMergeEvents = static_cast<int>(value);
                argi += 2;
                continue;
            }
            if (option == "--time-budget-ms") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char *endPtr = nullptr;
                const long long value = std::strtoll(argv[argi + 1], &endPtr, 10);
                if (endPtr == argv[argi + 1] || *endPtr != '\0' || value < 0) {
                    throw std::runtime_error("Invalid --time-budget-ms value. Expected a non-negative integer.");
                }
                timeBudgetMs = value;
                argi += 2;
                continue;
            }
            if (option == "--disable-leaf-fast-path") {
                leafFastPathEnabled = false;
                ++argi;
                continue;
            }
            if (option == "--radio-adj") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char *endPtr = nullptr;
                radioAdj = std::strtod(argv[argi + 1], &endPtr);
                if (endPtr == argv[argi + 1] || *endPtr != '\0' || radioAdj <= 0.0) {
                    throw std::runtime_error("Invalid --radio-adj value. Expected a positive radius.");
                }
                argi += 2;
                continue;
            }
            if (option == "--no-output") {
                noOutput = true;
                ++argi;
                continue;
            }
            throw std::runtime_error("Unknown option: " + option);
        }

        if ((!noOutput && argc - argi < 3) || (noOutput && argc - argi < 2)) {
            printUsage(argv[0]);
            return 1;
        }
        if (maxMergeEvents > 0 && mode != Mode::ExplicitUpdate) {
            throw std::runtime_error("--max-merge-events is only supported with --mode explicit");
        }
        if (timeBudgetMs > 0 && mode != Mode::ExplicitUpdate) {
            throw std::runtime_error("--time-budget-ms is only supported with --mode explicit");
        }
        if (stopOnFirstPhaseFallback && mode != Mode::ExplicitUpdate) {
            throw std::runtime_error("--stop-on-first-phase-fallback is only supported with --mode explicit");
        }
        if (findDivergentPrefix && mode != Mode::Compare) {
            throw std::runtime_error("--find-divergent-prefix is only supported with --mode compare");
        }
        if (dumpPrefixFixture && mode != Mode::Compare) {
            throw std::runtime_error("--dump-prefix-fixture is only supported with --mode compare");
        }
        if (dumpIncidentClosureFixture && mode != Mode::ExplicitUpdate) {
            throw std::runtime_error("--dump-incident-closure-fixture is only supported with --mode explicit");
        }

        const std::string inputPath = argv[argi++];
        const std::string outputPath = noOutput ? std::string() : argv[argi++];
        const auto thresholds = parseThresholds(argc, argv, argi);
        auto input = common::loadGrayImage(inputPath);
        auto adjacency = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), radioAdj);

        std::cout << "Input: " << inputPath << "\n";
        if (noOutput) {
            std::cout << "Output: disabled\n";
        } else {
            std::cout << "Output: " << outputPath << "\n";
        }
        std::cout << "Image size: " << input->getNumRows() << "x" << input->getNumCols() << "\n";
        std::cout << "Mode: "
                  << (mode == Mode::ExplicitUpdate ? "explicit"
                      : mode == Mode::Naive ? "naive"
                      : "compare")
                  << "\n";
        std::cout << "Attribute: " << attributeKindName(attributeKind) << "\n";
        std::cout << "Adjacency radius: " << radioAdj << "\n";
        std::cout << "Dual maintenance: " << dualMaintenanceName(dualMaintenance) << "\n";
        if (dualMaintenance == noninc::NonIncreasingDualMaintenance::OracleGuidedLocalRecomposition ||
            dualMaintenance == noninc::NonIncreasingDualMaintenance::Problem1CloseRecompose) {
            std::cout << "Leaf explicit fast path: " << (leafFastPathEnabled ? "enabled" : "disabled") << "\n";
        } else if (dualMaintenance == noninc::NonIncreasingDualMaintenance::StrippedThenOracleLocalRecomposition ||
                   dualMaintenance == noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition ||
                   dualMaintenance == noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition) {
            std::cout << "Leaf explicit fast path: bypassed by stripped-first mode\n";
            std::cout << "Stripped validation: " << strippedValidationMode << "\n";
            std::cout << "Stripped commit: " << strippedCommitModeName << "\n";
            std::cout << "Stripped fallback: "
                      << (dualMaintenance == noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition
                              ? "phase rebuild"
                              : dualMaintenance == noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition
                                    ? "event rebuild"
                                    : "oracle-guided local")
                      << "\n";
        }
        if (maxMergeEvents > 0) {
            std::cout << "Max merge events per threshold: " << maxMergeEvents << "\n";
        }
        if (timeBudgetMs > 0) {
            std::cout << "Time budget per threshold: " << timeBudgetMs << " ms\n";
        }
        if (diagnoseFirstDivergence) {
            std::cout << "Diagnostics: first divergence\n";
        }
        if (diagnoseEvents) {
            std::cout << "Diagnostics: merge event summary\n";
        }
        if (findDivergentPrefix) {
            std::cout << "Diagnostics: divergent prefix finder\n";
        }
        if (dumpPrefixFixture) {
            std::cout << "Diagnostics: prefix fixture dump"
                      << " prefix_events=" << prefixFixtureEvents
                      << " padding=" << prefixFixturePadding
                      << "\n";
        }
        if (dumpFirstLocalOracle) {
            std::cout << "Diagnostics: first local oracle dump\n";
        }
        if (dumpTargetLocalOracle) {
            std::cout << "Diagnostics: targeted local oracle dump"
                      << " threshold=" << targetLocalOracleThreshold
                      << " phase=" << phaseName(targetLocalOraclePhase)
                      << " node=" << targetLocalOracleNodeId
                      << "\n";
        }
        if (phaseCandidateGuard) {
            std::cout << "Phase candidate guard: enabled\n";
        }
        if (dumpPhaseOracleFixture) {
            std::cout << "Diagnostics: first phase-oracle mismatch fixture\n";
        }
        if (dumpSmallestPhaseOracleFixture) {
            std::cout << "Diagnostics: smallest phase-oracle mismatch fixture\n";
        }
        if (dumpIncidentClosureFixture) {
            std::cout << "Diagnostics: incident closure fixture\n";
        }
        if (strictDirectV1Precheck) {
            std::cout << "MCPC: direct-v1 strict precheck enabled\n";
        }
        if (mcpcNoCopyDirectCommit) {
            std::cout << "MCPC: no-copy direct commit enabled\n";
        }
        if (mcpcValidatedDirectCommit) {
            std::cout << "MCPC: validated direct commit enabled\n";
        }

        ImageUInt8Ptr output;
        if (mode == Mode::ExplicitUpdate) {
            if (maxMergeEvents <= 0 &&
                timeBudgetMs <= 0 &&
                !diagnoseEvents && !diagnoseFirstDivergence && !diagnosePhaseOracle && !dumpTargetLocalOracle) {
                const auto start = std::chrono::steady_clock::now();
                noninc::NonIncreasingCasfKernel runner(input, radioAdj, dualMaintenance, attributeKind);
                auto diagnostics = makeBaseDiagnosticsOptions(
                    strippedLocalUpdateOptions,
                    strippedCommitMode,
                    leafFastPathEnabled,
                    phaseCandidateGuard,
                    strictDirectV1Precheck,
                    mcpcNoCopyDirectCommit,
                    mcpcValidatedDirectCommit);
                noninc_debug::applyNonIncreasingCasfDiagnostics(runner, diagnostics);
                output = runner.filter(thresholds);
                const auto end = std::chrono::steady_clock::now();
                std::cout << "Total time: " << common::microsToMillis(common::elapsedMicros(start, end)) << " ms\n";
            } else {
                noninc::NonIncreasingCasfKernel runner(input, radioAdj, dualMaintenance, attributeKind);
                auto diagnostics = makeBaseDiagnosticsOptions(
                    strippedLocalUpdateOptions,
                    strippedCommitMode,
                    leafFastPathEnabled,
                    phaseCandidateGuard,
                    strictDirectV1Precheck,
                    mcpcNoCopyDirectCommit,
                    mcpcValidatedDirectCommit);
                diagnostics.collectMergeEvents =
                    diagnoseEvents ||
                    diagnoseFirstDivergence ||
                    stopOnFirstPhaseFallback ||
                    maxMergeEvents > 0 ||
                    timeBudgetMs > 0;
                diagnostics.collectPhaseOracleComparisons = diagnosePhaseOracle;
                diagnostics.stopOnFirstInvariantMismatch = diagnoseFirstDivergence;
                diagnostics.stopOnFirstPhaseRebuildFallback = stopOnFirstPhaseFallback;
                diagnostics.stopOnFirstPhaseReleaseRejection = dumpIncidentClosureFixture;
                diagnostics.collectFirstInvariantMismatchOracle = dumpFirstLocalOracle;
                diagnostics.maxAppliedMergeEventsPerThreshold = maxMergeEvents;
                diagnostics.maxThresholdMicros = timeBudgetMs * 1000;
                if (dumpTargetLocalOracle) {
                    diagnostics.captureLocalOracle = true;
                    diagnostics.localOracleThreshold = targetLocalOracleThreshold;
                    diagnostics.localOraclePhase = targetLocalOraclePhase;
                    diagnostics.localOracleNodeId = targetLocalOracleNodeId;
                }
                noninc_debug::applyNonIncreasingCasfDiagnostics(runner, diagnostics);

                long long explicitTotalUs = 0;
                bool invariantMismatch = false;
                bool targetOracleDumped = false;
                bool phaseOracleFixtureDumped = false;
                bool incidentClosureFixtureDumped = false;
                bool mergeEventLimitReached = false;
                bool timeBudgetReached = false;
                bool phaseFallbackStopReached = false;
                bool localUpdateRejectionStopReached = false;
                float mismatchThreshold = 0.0f;
                MergeEventStats aggregateStats;

                for (float threshold : thresholds) {
                    const auto start = std::chrono::steady_clock::now();
                    runner.applyThreshold(threshold);
                    const auto end = std::chrono::steady_clock::now();
                    explicitTotalUs += common::elapsedMicros(start, end);
                    output = runner.reconstructImage();

                    MergeEventStats thresholdStats;
                    addEventStats(thresholdStats, runner.lastMergeEvents());
                    thresholdStats.phaseDualRebuilds = static_cast<std::size_t>(runner.lastPhaseDualRebuilds());
                    thresholdStats.phaseCandidateGuardRebuilds =
                        static_cast<std::size_t>(runner.lastPhaseCandidateGuardRebuilds());
                    thresholdStats.phaseCandidateGuardCanonicalizationChanges =
                        static_cast<std::size_t>(runner.lastPhaseCandidateGuardCanonicalizationChanges());
                    addEventStats(aggregateStats, runner.lastMergeEvents());
                    aggregateStats.phaseDualRebuilds += thresholdStats.phaseDualRebuilds;
                    aggregateStats.phaseCandidateGuardRebuilds +=
                        thresholdStats.phaseCandidateGuardRebuilds;
                    aggregateStats.phaseCandidateGuardCanonicalizationChanges +=
                        thresholdStats.phaseCandidateGuardCanonicalizationChanges;
                    if (diagnoseEvents) {
                        printExplicitThresholdEventStats(threshold, thresholdStats, explicitTotalUs);
                        printFirstOracleStructureMismatchEvent(runner.lastMergeEvents());
                        printFirstPhaseRebuildFallbackEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcLiveEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcHiddenSupportRedundantLiveEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcParentLinkBridgeLiveEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcStableChildBridgeLiveEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcChildLinkBridgeLiveEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcEmptyOwnerReleaseLiveEvent(runner.lastMergeEvents());
                        printFirstPartialFMcpcStableChildMoveLiveEvent(runner.lastMergeEvents());
                    }
                    if (diagnosePhaseOracle) {
                        printPhaseOracleComparisons(runner.lastPhaseOracleComparisons());
                    }
                    if (dumpIncidentClosureFixture &&
                        runner.lastPhaseReleaseRejectionStopReached()) {
                        if (const auto *comparison =
                                firstIncidentClosureComparison(runner.lastPhaseOracleComparisons())) {
                            incidentClosureFixtureDumped = true;
                            printIncidentClosureFixture(*comparison);
                        }
                        break;
                    }
                    if (dumpTargetLocalOracle && runner.hasCapturedLocalOracle()) {
                        targetOracleDumped = true;
                        std::cout << "\nTarget local oracle before merge:\n"
                                  << noninc::debugString(runner.capturedLocalOracle());
                        const auto partitionPlan = noninc::buildLocalPartitionPlan(
                            runner.capturedLocalOracle(),
                            *adjacency);
                        const auto recompositionPlan = noninc::buildLocalRecompositionPlan(
                            runner.capturedLocalOracle(),
                            partitionPlan);
                        std::cout << noninc::debugString(partitionPlan)
                                  << noninc::debugString(recompositionPlan);
                        break;
                    }
                    if (runner.lastMergeEventLimitReached()) {
                        mergeEventLimitReached = true;
                        std::cout << "  merge_event_limit_reached=true"
                                  << " threshold=" << threshold
                                  << " max_events=" << maxMergeEvents
                                  << "\n";
                        break;
                    }
                    if (runner.lastTimeBudgetReached()) {
                        timeBudgetReached = true;
                        std::cout << "  time_budget_reached=true"
                                  << " threshold=" << threshold
                                  << " budget_ms=" << timeBudgetMs
                                  << " events=" << runner.lastMergeEvents().size()
                                  << "\n";
                        break;
                    }
                    if (runner.lastPhaseFallbackStopReached()) {
                        phaseFallbackStopReached = true;
                        std::cout << "  phase_fallback_stop_reached=true"
                                  << " threshold=" << threshold
                                  << " events=" << runner.lastMergeEvents().size()
                                  << "\n";
                        break;
                    }
                    if (runner.lastLocalUpdateRejectionStopReached()) {
                        localUpdateRejectionStopReached = true;
                        std::cout << "  local_update_rejection_stop_reached=true"
                                  << " threshold=" << threshold
                                  << " events=" << runner.lastMergeEvents().size()
                                  << "\n";
                        if (!runner.lastMergeEvents().empty()) {
                            std::cout << "First rejected local-update event:\n";
                            printMergeEvent(runner.lastMergeEvents().back(),
                                            runner.lastMergeEvents().size() - 1);
                        }
                        break;
                    }
                    if ((dumpPhaseOracleFixture || dumpSmallestPhaseOracleFixture) && !phaseOracleFixtureDumped) {
                        const auto comparisonSpan =
                            std::span<const noninc::NonIncreasingPhaseOracleComparison>(
                                runner.lastPhaseOracleComparisons());
                        const auto fixture = dumpSmallestPhaseOracleFixture
                                                 ? noninc::buildSmallestPhaseOracleMismatchFixture(comparisonSpan, 1)
                                                 : noninc::buildFirstPhaseOracleMismatchFixture(comparisonSpan, 1);
                        if (fixture.found) {
                            phaseOracleFixtureDumped = true;
                            printPhaseOracleMismatchFixture(
                                fixture,
                                dumpSmallestPhaseOracleFixture ? "Smallest" : "First");
                            printPhaseOracleFixtureEventLinks(fixture, runner.lastMergeEvents());
                            break;
                        }
                    }
                    if (runner.hasFirstInvariantMismatch()) {
                        invariantMismatch = true;
                        mismatchThreshold = threshold;
                        std::cout << "\nFirst reconstruction-invariant mismatch inside threshold: "
                                  << mismatchThreshold << "\n";
                        std::cout << "First failing merge event:\n";
                        printMergeEvent(runner.firstInvariantMismatchEvent(),
                                        runner.lastMergeEvents().empty() ? 0 : runner.lastMergeEvents().size() - 1);
                        if (dumpFirstLocalOracle && runner.hasFirstInvariantMismatchOracle()) {
                            std::cout << "\nLocal oracle before failing merge:\n"
                                      << noninc::debugString(runner.firstInvariantMismatchOracle());
                        }
                        printMergeEventSummary(runner.lastMergeEvents());
                        break;
                    }
                }

                    if (!invariantMismatch) {
                        if (targetOracleDumped) {
                            std::cout << "\nTarget local oracle captured.\n";
                        } else if (phaseOracleFixtureDumped) {
                            std::cout << "\nPhase-oracle mismatch fixture captured.\n";
                        } else if (incidentClosureFixtureDumped) {
                            std::cout << "\nIncident closure fixture captured.\n";
                        } else if (mergeEventLimitReached) {
                            std::cout << "\nStopped after the configured merge-event limit.\n";
                        } else if (timeBudgetReached) {
                            std::cout << "\nStopped after the configured time budget.\n";
                        } else if (phaseFallbackStopReached) {
                            std::cout << "\nStopped at the first phase-rebuild fallback.\n";
                        } else if (localUpdateRejectionStopReached) {
                            std::cout << "\nStopped at the first local-update domain rejection.\n";
                        } else if (dumpTargetLocalOracle) {
                            std::cout << "\nTarget local oracle was not found in the provided threshold sequence.\n";
                        } else if (dumpIncidentClosureFixture) {
                            std::cout << "\nIncident closure fixture was not found in the provided threshold sequence.\n";
                        } else if (diagnoseFirstDivergence) {
                        std::cout << "\nNo reconstruction-invariant mismatch found in the provided threshold sequence.\n";
                    } else {
                        std::cout << "\nReconstruction-invariant check: disabled.\n";
                    }
                    if (diagnoseEvents) {
                        printAggregateEventStats(aggregateStats);
                    }
                } else {
                    throw std::runtime_error("explicit non-increasing CASF broke the reconstruction invariant");
                }
            }
        } else if (mode == Mode::Naive) {
            const auto start = std::chrono::steady_clock::now();
            output = noninc::applyNaiveAttributeMergeSequence(input, adjacency, thresholds, attributeKind);
            const auto end = std::chrono::steady_clock::now();
            std::cout << "Total time: " << common::microsToMillis(common::elapsedMicros(start, end)) << " ms\n";
        } else if (dumpPrefixFixture) {
            if (thresholds.size() != 1) {
                throw std::runtime_error("--dump-prefix-fixture requires exactly one threshold");
            }
            dumpPrefixFixtureDiagnostic(
                input,
                adjacency,
                radioAdj,
                thresholds.front(),
                dualMaintenance,
                attributeKind,
                strippedLocalUpdateOptions,
                strippedCommitMode,
                leafFastPathEnabled,
                phaseCandidateGuard,
                strictDirectV1Precheck,
                mcpcNoCopyDirectCommit,
                mcpcValidatedDirectCommit,
                prefixFixtureEvents,
                prefixFixturePadding);
            output = input->clone();
        } else if (findDivergentPrefix) {
            if (thresholds.size() != 1) {
                throw std::runtime_error("--find-divergent-prefix requires exactly one threshold");
            }
            findAndPrintDivergentPrefix(
                input,
                adjacency,
                radioAdj,
                thresholds.front(),
                dualMaintenance,
                attributeKind,
                strippedLocalUpdateOptions,
                strippedCommitMode,
                leafFastPathEnabled,
                phaseCandidateGuard,
                strictDirectV1Precheck,
                mcpcNoCopyDirectCommit,
                mcpcValidatedDirectCommit);
            output = input->clone();
        } else if (!diagnoseFirstDivergence && !diagnoseEvents && !diagnosePhaseOracle && !dumpTargetLocalOracle) {
            const auto explicitStart = std::chrono::steady_clock::now();
            noninc::NonIncreasingCasfKernel runner(input, radioAdj, dualMaintenance, attributeKind);
            auto diagnostics = makeBaseDiagnosticsOptions(
                strippedLocalUpdateOptions,
                strippedCommitMode,
                leafFastPathEnabled,
                phaseCandidateGuard,
                strictDirectV1Precheck,
                mcpcNoCopyDirectCommit,
                mcpcValidatedDirectCommit);
            noninc_debug::applyNonIncreasingCasfDiagnostics(runner, diagnostics);
            const auto explicitOutput = runner.filter(thresholds);
            const auto explicitEnd = std::chrono::steady_clock::now();

            const auto naiveStart = std::chrono::steady_clock::now();
            const auto naiveOutput =
                noninc::applyNaiveAttributeMergeSequence(input, adjacency, thresholds, attributeKind);
            const auto naiveEnd = std::chrono::steady_clock::now();

            output = explicitOutput;
            std::cout << "explicit total: "
                      << common::microsToMillis(common::elapsedMicros(explicitStart, explicitEnd))
                      << " ms\n";
            std::cout << "naive total: "
                      << common::microsToMillis(common::elapsedMicros(naiveStart, naiveEnd))
                      << " ms\n";
            std::cout << "explicit == naive: "
                      << (explicitOutput->isEqual(naiveOutput) ? "true" : "false") << "\n";
            if (!explicitOutput->isEqual(naiveOutput)) {
                throw std::runtime_error("explicit non-increasing CASF diverged from naive baseline");
            }
        } else {
            noninc::NonIncreasingCasfKernel runner(input, radioAdj, dualMaintenance, attributeKind);
            auto diagnostics = makeBaseDiagnosticsOptions(
                strippedLocalUpdateOptions,
                strippedCommitMode,
                leafFastPathEnabled,
                phaseCandidateGuard,
                strictDirectV1Precheck,
                mcpcNoCopyDirectCommit,
                mcpcValidatedDirectCommit);
            diagnostics.collectMergeEvents = diagnoseEvents || diagnoseFirstDivergence;
            diagnostics.collectPhaseOracleComparisons = diagnosePhaseOracle;
            diagnostics.stopOnFirstInvariantMismatch = diagnoseFirstDivergence;
            diagnostics.collectFirstInvariantMismatchOracle = dumpFirstLocalOracle;
            if (dumpTargetLocalOracle) {
                diagnostics.captureLocalOracle = true;
                diagnostics.localOracleThreshold = targetLocalOracleThreshold;
                diagnostics.localOraclePhase = targetLocalOraclePhase;
                diagnostics.localOracleNodeId = targetLocalOracleNodeId;
            }
            noninc_debug::applyNonIncreasingCasfDiagnostics(runner, diagnostics);

            auto naiveOutput = input->clone();
            long long explicitTotalUs = 0;
            long long naiveTotalUs = 0;
            bool diverged = false;
            bool targetOracleDumped = false;
            bool phaseOracleFixtureDumped = false;
            bool localUpdateRejectionStopReached = false;
            float divergentThreshold = 0.0f;
            MergeEventStats aggregateStats;

            for (float threshold : thresholds) {
                const auto explicitStart = std::chrono::steady_clock::now();
                runner.applyThreshold(threshold);
                const auto explicitEnd = std::chrono::steady_clock::now();
                explicitTotalUs += common::elapsedMicros(explicitStart, explicitEnd);

                output = runner.reconstructImage();
                if (dumpTargetLocalOracle && runner.hasCapturedLocalOracle()) {
                    targetOracleDumped = true;
                    std::cout << "\nTarget local oracle before merge:\n"
                              << noninc::debugString(runner.capturedLocalOracle());
                    const auto partitionPlan = noninc::buildLocalPartitionPlan(
                        runner.capturedLocalOracle(),
                        *adjacency);
                    const auto recompositionPlan = noninc::buildLocalRecompositionPlan(
                        runner.capturedLocalOracle(),
                        partitionPlan);
                    std::cout << noninc::debugString(partitionPlan)
                              << noninc::debugString(recompositionPlan);
                    break;
                }
                if (runner.hasFirstInvariantMismatch()) {
                    MergeEventStats thresholdStats;
                    addEventStats(thresholdStats, runner.lastMergeEvents());
                    thresholdStats.phaseDualRebuilds = static_cast<std::size_t>(runner.lastPhaseDualRebuilds());
                    thresholdStats.phaseCandidateGuardRebuilds =
                        static_cast<std::size_t>(runner.lastPhaseCandidateGuardRebuilds());
                    thresholdStats.phaseCandidateGuardCanonicalizationChanges =
                        static_cast<std::size_t>(runner.lastPhaseCandidateGuardCanonicalizationChanges());
                    addEventStats(aggregateStats, runner.lastMergeEvents());
                    aggregateStats.phaseDualRebuilds += thresholdStats.phaseDualRebuilds;
                    aggregateStats.phaseCandidateGuardRebuilds +=
                        thresholdStats.phaseCandidateGuardRebuilds;
                    aggregateStats.phaseCandidateGuardCanonicalizationChanges +=
                        thresholdStats.phaseCandidateGuardCanonicalizationChanges;
                    diverged = true;
                    divergentThreshold = threshold;
                    std::cout << "  threshold=" << threshold
                              << " explicit_ms=" << common::microsToMillis(explicitTotalUs)
                              << " invariant_mismatch=true"
                              << "\n";
                    std::cout << "\nFirst reconstruction-invariant mismatch inside threshold: "
                              << divergentThreshold << "\n";
                    std::cout << "First failing merge event:\n";
                    printMergeEvent(runner.firstInvariantMismatchEvent(), runner.lastMergeEvents().empty() ? 0 : runner.lastMergeEvents().size() - 1);
                    if (dumpFirstLocalOracle && runner.hasFirstInvariantMismatchOracle()) {
                        std::cout << "\nLocal oracle before failing merge:\n"
                                  << noninc::debugString(runner.firstInvariantMismatchOracle());
                    }
                    if (diagnosePhaseOracle) {
                        printPhaseOracleComparisons(runner.lastPhaseOracleComparisons());
                    }
                    printThresholdEventStats(threshold,
                                             thresholdStats,
                                             explicitTotalUs,
                                             naiveTotalUs,
                                             false);
                    printMergeEventSummary(runner.lastMergeEvents());
                    break;
                }

                if (runner.lastLocalUpdateRejectionStopReached()) {
                    MergeEventStats thresholdStats;
                    addEventStats(thresholdStats, runner.lastMergeEvents());
                    thresholdStats.phaseDualRebuilds = static_cast<std::size_t>(runner.lastPhaseDualRebuilds());
                    thresholdStats.phaseCandidateGuardRebuilds =
                        static_cast<std::size_t>(runner.lastPhaseCandidateGuardRebuilds());
                    thresholdStats.phaseCandidateGuardCanonicalizationChanges =
                        static_cast<std::size_t>(runner.lastPhaseCandidateGuardCanonicalizationChanges());
                    addEventStats(aggregateStats, runner.lastMergeEvents());
                    aggregateStats.phaseDualRebuilds += thresholdStats.phaseDualRebuilds;
                    aggregateStats.phaseCandidateGuardRebuilds +=
                        thresholdStats.phaseCandidateGuardRebuilds;
                    aggregateStats.phaseCandidateGuardCanonicalizationChanges +=
                        thresholdStats.phaseCandidateGuardCanonicalizationChanges;
                    localUpdateRejectionStopReached = true;
                    divergentThreshold = threshold;
                    std::cout << "  threshold=" << threshold
                              << " explicit_ms=" << common::microsToMillis(explicitTotalUs)
                              << " local_update_rejection_stop=true"
                              << " events=" << runner.lastMergeEvents().size()
                              << "\n";
                    if (!runner.lastMergeEvents().empty()) {
                        std::cout << "\nFirst rejected local-update event:\n";
                        printMergeEvent(runner.lastMergeEvents().back(),
                                        runner.lastMergeEvents().size() - 1);
                    }
                    if (diagnoseEvents) {
                        printThresholdEventStats(threshold,
                                                 thresholdStats,
                                                 explicitTotalUs,
                                                 naiveTotalUs,
                                                 false);
                    }
                    printMergeEventSummary(runner.lastMergeEvents());
                    break;
                }

                const auto naiveStart = std::chrono::steady_clock::now();
                naiveOutput = noninc::applyNaiveAttributeMergeThreshold(
                    naiveOutput,
                    adjacency,
                    threshold,
                    attributeKind);
                const auto naiveEnd = std::chrono::steady_clock::now();
                naiveTotalUs += common::elapsedMicros(naiveStart, naiveEnd);

                const bool matches = output->isEqual(naiveOutput);
                MergeEventStats thresholdStats;
                addEventStats(thresholdStats, runner.lastMergeEvents());
                thresholdStats.phaseDualRebuilds = static_cast<std::size_t>(runner.lastPhaseDualRebuilds());
                thresholdStats.phaseCandidateGuardRebuilds =
                    static_cast<std::size_t>(runner.lastPhaseCandidateGuardRebuilds());
                thresholdStats.phaseCandidateGuardCanonicalizationChanges =
                    static_cast<std::size_t>(runner.lastPhaseCandidateGuardCanonicalizationChanges());
                addEventStats(aggregateStats, runner.lastMergeEvents());
                aggregateStats.phaseDualRebuilds += thresholdStats.phaseDualRebuilds;
                aggregateStats.phaseCandidateGuardRebuilds +=
                    thresholdStats.phaseCandidateGuardRebuilds;
                aggregateStats.phaseCandidateGuardCanonicalizationChanges +=
                    thresholdStats.phaseCandidateGuardCanonicalizationChanges;
                if (diagnoseEvents) {
                    printThresholdEventStats(threshold,
                                             thresholdStats,
                                             explicitTotalUs,
                                             naiveTotalUs,
                                             matches);
                    printFirstOracleStructureMismatchEvent(runner.lastMergeEvents());
                    printFirstPhaseRebuildFallbackEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcLiveEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcHiddenSupportRedundantLiveEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcParentLinkBridgeLiveEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcStableChildBridgeLiveEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcChildLinkBridgeLiveEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcEmptyOwnerReleaseLiveEvent(runner.lastMergeEvents());
                    printFirstPartialFMcpcStableChildMoveLiveEvent(runner.lastMergeEvents());
                } else {
                    std::cout << "  threshold=" << threshold
                              << " explicit_ms=" << common::microsToMillis(explicitTotalUs)
                              << " naive_ms=" << common::microsToMillis(naiveTotalUs)
                              << " matches=" << (matches ? "true" : "false")
                              << "\n";
                }
                if (diagnosePhaseOracle) {
                    printPhaseOracleComparisons(runner.lastPhaseOracleComparisons());
                }
                if ((dumpPhaseOracleFixture || dumpSmallestPhaseOracleFixture) && !phaseOracleFixtureDumped) {
                    const auto comparisonSpan =
                        std::span<const noninc::NonIncreasingPhaseOracleComparison>(
                            runner.lastPhaseOracleComparisons());
                    const auto fixture = dumpSmallestPhaseOracleFixture
                                             ? noninc::buildSmallestPhaseOracleMismatchFixture(comparisonSpan, 1)
                                             : noninc::buildFirstPhaseOracleMismatchFixture(comparisonSpan, 1);
                    if (fixture.found) {
                        phaseOracleFixtureDumped = true;
                        printPhaseOracleMismatchFixture(
                            fixture,
                            dumpSmallestPhaseOracleFixture ? "Smallest" : "First");
                        printPhaseOracleFixtureEventLinks(fixture, runner.lastMergeEvents());
                        break;
                    }
                }

                if (!matches) {
                    diverged = true;
                    divergentThreshold = threshold;
                    std::cout << "\nFirst divergent threshold: " << divergentThreshold << "\n";
                    printMergeEventSummary(runner.lastMergeEvents());
                    break;
                }
            }

            if (!diverged && !localUpdateRejectionStopReached) {
                output = runner.reconstructImage();
                if (targetOracleDumped) {
                    std::cout << "\nTarget local oracle captured.\n";
                } else if (phaseOracleFixtureDumped) {
                    std::cout << "\nPhase-oracle mismatch fixture captured.\n";
                } else if (dumpTargetLocalOracle) {
                    std::cout << "\nTarget local oracle was not found in the provided threshold sequence.\n";
                } else {
                    std::cout << "\nNo divergence found in the provided threshold sequence.\n";
                }
                if (diagnoseEvents) {
                    printAggregateEventStats(aggregateStats);
                }
            } else {
                throw std::runtime_error("explicit non-increasing CASF diverged from naive baseline");
            }
        }

        if (!noOutput) {
            common::saveGrayImage(outputPath, output);
            std::cout << "\nSaved filtered image to " << outputPath << "\n";
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "non_increasing_casf_apply: " << e.what() << "\n";
        return 1;
    }
}
