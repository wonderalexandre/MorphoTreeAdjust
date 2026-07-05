#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../AdjacencyRelation.hpp"
#include "../AttributeComputer.hpp"
#include "../Common.hpp"
#include "../DynamicComponentTree.hpp"
#include "../DualMinMaxTreeIncrementalFilter.hpp"
#include "../debug/details/ExplicitRegionProblem1Preflight.hpp"
#include "../debug/details/RegionMaterializationPreflight.hpp"
#include "../debug/details/AttributeDirtyRecomputer.hpp"
#include "../debug/details/NonIncreasingMcpcDirectPreflight.hpp"
#include "NonIncreasingLocalUpdate.hpp"
#include "NonIncreasingLocalOracle.hpp"
#include "../debug/details/RegionLevelChange.hpp"
#include "../debug/details/RegionMergeUpdater.hpp"

namespace morphotree::details::non_increasing {


using ::debug::problem1::Problem1PreflightResult;
using ::debug::problem1::Problem1PreflightStatus;
using ::debug::problem1::checkExplicitRegionProblem1Preflight;
using ::debug::problem1::hasExplicitRegionProblem1TargetLevelUnionCandidate;
using ::debug::region::RegionMaterializationPlan;
using ::debug::region::analyzeRegionMaterialization;

using ::debug::non_increasing::NonIncreasingMcpcD12DomainDiagnostics;
using ::debug::non_increasing::NonIncreasingMcpcDirectCommitCandidateDiagnostics;
using ::debug::non_increasing::NonIncreasingMcpcDirectPreflightArtifacts;
using ::debug::non_increasing::NonIncreasingMcpcDirectPreflightResult;
using ::debug::non_increasing::NonIncreasingMcpcDirectSummaryComparison;
using ::debug::non_increasing::NonIncreasingMcpcDirectV1AcceptanceDiagnostics;
using ::debug::non_increasing::NonIncreasingMcpcLiveMutationDryRun;
using ::debug::non_increasing::NonIncreasingMcpcNoCopyCandidateDiagnostics;
using ::debug::non_increasing::applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight;
using ::debug::non_increasing::buildMcpcDirectPreflightArtifacts;
using ::debug::non_increasing::compareMcpcDirectPreflightSummaryWithStrippedUpdate;
using ::debug::non_increasing::diagnoseMcpcD12Domain;
using ::debug::non_increasing::diagnoseMcpcDirectCommitCandidate;
using ::debug::non_increasing::diagnoseMcpcDirectPreflightV1Acceptance;
using ::debug::non_increasing::diagnoseMcpcLiveMutationDryRun;
using ::debug::non_increasing::diagnoseMcpcNoCopyCandidate;
using ::debug::non_increasing::evaluateMcpcDirectV1PrecheckAcceptance;
using ::debug::non_increasing::preflightMcpcDirectMergeNodeUpdate;

inline long long nonIncreasingElapsedMicros(std::chrono::steady_clock::time_point start,
                                            std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

enum class NonIncreasingDualMaintenance {
    PartialF,
    ExplicitOnly,
    MaterializeThenExplicit,
    RebuildAfterPhase,
    OracleGuidedLocalRecomposition,
    Problem1CloseRecompose,
    StrippedThenOracleLocalRecomposition,
    StrippedThenRebuildLocalRecomposition,
    StrippedThenPhaseRebuildLocalRecomposition,
    MergeNodeUpdateThenPhaseRebuild,
    MergeNodeCertifiedThenPhaseRebuild,
    MergeNodeCertifiedCreateBThenPhaseRebuild,
    MergeNodeCertifiedCreateBMcpcThenPhaseRebuild,
    MergeNodeCleanUnchecked,
};

enum class NonIncreasingStrippedCommitMode {
    SafeCopy,
    Direct,
};

enum class NonIncreasingAttributeKind {
    GrayVariance,
    BoxAreaPerimeterRatio,
};

enum class NonIncreasingPhase {
    MaxToMin,
    MinToMax,
};

inline bool usesNonIncreasingPhaseRebuildFallback(NonIncreasingDualMaintenance mode) {
    return mode == NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition ||
           mode == NonIncreasingDualMaintenance::MergeNodeUpdateThenPhaseRebuild ||
           mode == NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild ||
           mode == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild ||
           mode == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild;
}

template<typename PixelType = AltitudeType>
struct NonIncreasingProblem1ExplicitRegionPreflight {
    Problem1PreflightResult<PixelType> preflight;
    RegionMaterializationPlan materialization;
    bool materializationAnalyzed = false;
    bool hasTargetLevelUnionCandidate = false;

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

    int affectedOwners() const {
        return materializationAnalyzed
            ? static_cast<int>(materialization.affectedOwners.size())
            : 0;
    }
};

inline std::string problem1ExplicitRegionPreflightMessage(
    Problem1PreflightStatus status) {
    switch (status) {
        case Problem1PreflightStatus::Accepted:
            return "";
        case Problem1PreflightStatus::EmptyRegion:
            return "Problem 1 requires a non-empty region C";
        case Problem1PreflightStatus::OutOfBoundsPixel:
            return "Problem 1 region C contains an out-of-bounds pixel";
        case Problem1PreflightStatus::DuplicatePixel:
            return "Problem 1 region C contains duplicated pixels";
        case Problem1PreflightStatus::MissingLiveOwner:
            return "Problem 1 region C contains a pixel without a live owner";
        case Problem1PreflightStatus::DisconnectedRegion:
            return "Problem 1 requires connected region C";
        case Problem1PreflightStatus::InvalidTargetLevel:
            return "Problem 1 target level is outside the admissible monotone interval";
        case Problem1PreflightStatus::RequiresLocalMaterialization:
            return "C cuts direct proper parts; local materialization is required first";
        case Problem1PreflightStatus::RequiresTargetNodeCreation:
            return "current updater has no level-b union candidate and must create a local node";
        case Problem1PreflightStatus::UnsupportedLocalTopology:
            return "current updater rejected this local topology";
    }
    return "unknown Problem 1 preflight status";
}

template<typename PixelType>
inline bool nonIncreasingProblem1HasTargetLevelUnionCandidate(
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsC,
    const Problem1PreflightResult<PixelType> &preflight) {
    if (!preflight.acceptedByPublishedProblem()) {
        return false;
    }

    std::vector<unsigned char> inRegion(
        static_cast<size_t>(dualTree.getNumTotalProperParts()),
        0);
    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
            return false;
        }
        inRegion[static_cast<size_t>(pixelId)] = 1;
    }

    return hasExplicitRegionProblem1TargetLevelUnionCandidate<PixelType>(
        dualTree,
        adjacency,
        pixelsC,
        inRegion,
        preflight.a,
        preflight.b);
}

template<typename PixelType>
inline NonIncreasingProblem1ExplicitRegionPreflight<PixelType>
analyzeNonIncreasingProblem1ExplicitRegionPreflight(
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false) {
    NonIncreasingProblem1ExplicitRegionPreflight<PixelType> out;
    out.preflight = checkExplicitRegionProblem1Preflight<PixelType>(
        dualTree,
        adjacency,
        pixelsC,
        b,
        assumeConnected);

    if (out.preflight.acceptedByPublishedProblem()) {
        out.materialization = analyzeRegionMaterialization(
            dualTree,
            adjacency,
            pixelsC);
        out.materializationAnalyzed = true;
        out.hasTargetLevelUnionCandidate =
            nonIncreasingProblem1HasTargetLevelUnionCandidate(
                dualTree,
                adjacency,
                pixelsC,
                out.preflight);
    }

    return out;
}

template <typename PixelType>
inline std::string classifyProblem1ExplicitRegionFailureMode(
    const NonIncreasingProblem1ExplicitRegionPreflight<PixelType> &analysis) {
    switch (analysis.preflight.status) {
        case Problem1PreflightStatus::EmptyRegion:
            return "empty_C";
        case Problem1PreflightStatus::OutOfBoundsPixel:
            return "invalid_pixel";
        case Problem1PreflightStatus::DuplicatePixel:
            return "duplicated_pixels";
        case Problem1PreflightStatus::MissingLiveOwner:
            return "dead_owner_pixel";
        case Problem1PreflightStatus::DisconnectedRegion:
            return "non_connected_C";
        case Problem1PreflightStatus::InvalidTargetLevel:
            return "invalid_target_level";
        case Problem1PreflightStatus::RequiresLocalMaterialization:
            return "non_materialized_C";
        case Problem1PreflightStatus::RequiresTargetNodeCreation:
            return "requires_node_creation";
        case Problem1PreflightStatus::UnsupportedLocalTopology:
            return "current_updater_limitation";
        case Problem1PreflightStatus::Accepted:
            break;
    }
    if (analysis.affectedOwners() > 1) {
        return "multiple_nodeCa";
    }
    return "contract_rejection";
}

struct NonIncreasingLocalOracleCaptureTarget {
    float threshold = 0.0f;
    NonIncreasingPhase phase = NonIncreasingPhase::MaxToMin;
    NodeId nodeId = InvalidNode;
    bool enabled = false;
};

struct NonIncreasingMergeEvent {
    bool collectDetailedDiagnostics = true;
    float threshold = 0.0f;
    NonIncreasingPhase phase = NonIncreasingPhase::MaxToMin;
    NodeId nodeId = InvalidNode;
    NodeId parentId = InvalidNode;
    int nodeLevel = 0;
    int parentLevel = 0;
    float attributeValue = 0.0f;
    int numChildren = 0;
    int directProperParts = 0;
    int componentCount = 0;
    int totalComponentPixels = 0;
    int minRow = -1;
    int minCol = -1;
    int maxRow = -1;
    int maxCol = -1;
    std::vector<int> componentSizes;
    std::vector<NodeId> componentNodeCa;
    std::vector<int> componentAltitudeCa;
    std::vector<PixelId> affectedPixelsPreview;
    bool applied = false;
    bool dualRebuiltAfterLocalMismatch = false;
    bool strippedLocalAttempted = false;
    bool strippedLocalAccepted = false;
    bool strippedLocalApplied = false;
    bool strippedDirectCommit = false;
    bool strippedLocalPlanAuditChecked = false;
    bool strippedLocalMatchesTarget = false;
    bool strippedLocalTargetChecked = false;
    bool strippedLocalStructureChecked = false;
    bool strippedLocalStructureMatchesRebuild = false;
    bool strippedLocalPlanOracleFree = false;
    bool strippedMaterializationAuditChecked = false;
    bool strippedMaterializationInputOracleFree = false;
    bool strippedFinalReconstructionChecked = false;
    bool strippedLocalDualMatchesPrimalAfterMerge = false;
    bool strippedLocalFallbackToOracle = false;
    bool strippedLocalFallbackToRebuild = false;
    bool strippedLocalFallbackToPhaseRebuild = false;
    bool strippedLocalSkippedByPhaseRebuild = false;
    bool mergeNodeLocalAttempted = false;
    bool mergeNodeLocalAccepted = false;
    bool mergeNodeLocalApplied = false;
    bool mergeNodeLocalFallbackToPhaseRebuild = false;
    bool mergeNodeLocalSkippedByPhaseRebuild = false;
    bool mergeNodeLocalDualMatchesPrimalBeforeMerge = false;
    bool mergeNodeLocalDualMatchesPrimalAfterMerge = false;
    bool mergeNodeLocalCertificationChecked = false;
    bool mergeNodeLocalCertificationPassed = false;
    bool mergeNodeLocalDirectCommit = false;
    bool mergeNodeLocalTargetCreationAllowed = false;
    bool mergeNodeLocalTargetCreationRequired = false;
    int mergeNodeLocalComponentCount = 0;
    int mergeNodeLocalPixelCount = 0;
    int mergeNodeLocalRejectedComponents = 0;
    int mergeNodeLocalPartialOwnerComponents = 0;
    int mergeNodeLocalRawMissingTargetComponents = 0;
    int mergeNodeLocalMissingTargetComponents = 0;
    int mergeNodeLocalPartialOwnerAndMissingTargetComponents = 0;
    int mergeNodeLocalProblem1PreflightRejectedComponents = 0;
    int mergeNodeLocalProblem1PreflightRejectedByPublishedProblem = 0;
    int mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation = 0;
    std::string mergeNodeLocalProblem1FailureMode;
    long long mergeNodeLocalCopyMicros = 0;
    long long mergeNodeLocalCertificationMicros = 0;
    long long mergeNodeComponentCollectMicros = 0;
    long long mergeNodeLocalUpdateMicros = 0;
    long long mergeNodePrimalMergeMicros = 0;
    long long mergeNodeEmptyReleaseMicros = 0;
    std::string mergeNodeLocalFailureSummary;
    bool mcpcDirectPreflightAttempted = false;
    bool mcpcDirectPreflightChecked = false;
    bool mcpcDirectPreflightPassed = false;
    bool mcpcDirectPreflightDirectV1PrecheckEvaluated = false;
    bool mcpcDirectPreflightDirectV1PrecheckAccepted = false;
    bool mcpcDirectPreflightDirectV1Evaluated = false;
    bool mcpcDirectPreflightDirectV1Accepted = false;
    bool mcpcDirectPreflightDirectV1StrictGateEnabled = false;
    bool mcpcDirectPreflightDirectV1StrictGateRejected = false;
    bool mcpcDirectPreflightPlanOracleFree = false;
    bool mcpcDirectPreflightMaterializationOracleFree = false;
    bool mcpcDirectPreflightMaterializationInputsValid = false;
    bool mcpcDirectPreflightSummaryBuilt = false;
    bool mcpcDirectPreflightSharedWithStrippedUpdate = false;
    bool mcpcDirectPreflightSummaryComparedWithStrippedUpdate = false;
    bool mcpcDirectPreflightSummaryMatchesStrippedUpdate = false;
    bool mcpcNoCopyCandidateEvaluated = false;
    bool mcpcNoCopyCandidateAccepted = false;
    bool mcpcNoCopyDirectCommitEnabled = false;
    bool mcpcNoCopyDirectCommitAttempted = false;
    bool mcpcNoCopyDirectCommitAccepted = false;
    bool mcpcNoCopyDirectCommitFellBackToPhaseRebuild = false;
    bool mcpcNoCopyDirectCommitDualMatchesPrimal = false;
    bool mcpcDirectCommitCandidateEvaluated = false;
    bool mcpcDirectCommitCandidateAccepted = false;
    bool mcpcValidatedDirectCommitEnabled = false;
    bool mcpcValidatedDirectCommitAttempted = false;
    bool mcpcValidatedDirectCommitAccepted = false;
    bool mcpcValidatedDirectCommitFellBackToCopy = false;
    bool mcpcValidatedDirectCommitDualMatchesPrimal = false;
    bool mcpcValidatedDirectCommitMatchesCopiedOutput = false;
    bool partialFAttempted = false;
    bool partialFAccepted = false;
    bool partialFRejected = false;
    bool partialFMutationAttempted = false;
    bool partialFMcpcLiveUpdate = false;
    bool partialFMcpcHiddenSupportRedundantLiveUpdate = false;
    bool partialFMcpcParentLinkBridgeLiveUpdate = false;
    bool partialFMcpcStableChildBridgeLiveUpdate = false;
    bool partialFMcpcChildLinkBridgeLiveUpdate = false;
    bool partialFMcpcEmptyOwnerReleaseLiveUpdate = false;
    bool partialFMcpcStableChildMoveLiveUpdate = false;
    std::string partialFRejectReasonCode;
    int mcpcDirectPreflightAffectedLevels = 0;
    int mcpcDirectPreflightLocalPixels = 0;
    int mcpcDirectPreflightPlannedNodes = 0;
    int mcpcDirectPreflightPlannedDirectProperParts = 0;
    int mcpcDirectPreflightPlannedExternalAttachments = 0;
    int mcpcDirectPreflightSourceOwnersToBeEmptied = 0;
    int mcpcDirectPreflightStableChildrenToMove = 0;
    int mcpcDirectPreflightPreviousParentsPossiblyEmptied = 0;
    int mcpcDirectPreflightInvalidStableChildRefs = 0;
    int mcpcDirectPreflightStableChildCycleRisks = 0;
    int mcpcDirectPreflightPredictedAuxiliaryNodes = 0;
    int mcpcDirectPreflightPredictedHiddenExternalParentConnectorNodes = 0;
    int mcpcDirectPreflightPredictedStableChildBridgeNodes = 0;
    int mcpcDirectPreflightSummaryMismatchCount = 0;
    int mcpcDirectPreflightSummaryPlannedNodeCountMismatch = 0;
    int mcpcDirectPreflightSummaryDirectProperPartCountMismatch = 0;
    int mcpcDirectPreflightSummaryAuxiliaryNodeCountMismatch = 0;
    int mcpcDirectPreflightSummaryHiddenExternalParentConnectorMismatch = 0;
    int mcpcDirectPreflightSummaryStableChildBridgeMismatch = 0;
    int mcpcDirectPreflightSummaryEmptiedLocalNodeMismatch = 0;
    int mcpcDirectPreflightSummaryMovedStableChildMismatch = 0;
    int mcpcDirectPreflightPredictedEmptiedLocalNodes = 0;
    int mcpcDirectPreflightObservedEmptiedLocalNodes = 0;
    int mcpcDirectPreflightCandidateNonOracleHiddenSupportConnectors = 0;
    int mcpcDirectPreflightCandidateNonOracleParentLinkBridges = 0;
    int mcpcDirectPreflightCandidateNonOracleChildLinkBridges = 0;
    int mcpcDirectPreflightV1PrecheckHiddenSupportCandidates = 0;
    int mcpcDirectPreflightV1PrecheckHiddenSupportAlreadyRepresented = 0;
    int mcpcDirectPreflightV1PrecheckHiddenSupportNotRedundant = 0;
    int mcpcDirectPreflightV1PrecheckParentLinksTouchingPlan = 0;
    int mcpcDirectPreflightV1PrecheckParentLinksBothRepresented = 0;
    int mcpcDirectPreflightV1PrecheckParentLinksNotBothRepresented = 0;
    int mcpcDirectPreflightV1PrecheckParentLinksBlockedByModifiedSnapshot = 0;
    int mcpcDirectPreflightV1PrecheckParentLinksInvalidLevel = 0;
    std::vector<int> mcpcDirectPreflightPlannedNodeLevels;
    long long mcpcDirectPreflightMicros = 0;
    long long mcpcNoCopyDirectCommitMicros = 0;
    long long mcpcValidatedDirectCommitMicros = 0;
    std::string mcpcDirectPreflightFailureSummary;
    std::string mcpcDirectPreflightSummaryMismatchSummary;
    std::string mcpcDirectPreflightDirectV1PrecheckRejectionSummary;
    std::string mcpcDirectPreflightDirectV1RejectionSummary;
    std::string mcpcNoCopyCandidateRejectionSummary;
    std::string mcpcNoCopyDirectCommitFailureSummary;
    std::string mcpcDirectCommitCandidateRejectionSummary;
    std::string mcpcValidatedDirectCommitFailureSummary;
    int strippedLocalRecompositionNodes = 0;
    int strippedLocalCreatedNodes = 0;
    int strippedLocalPixelCount = 0;
    std::vector<int> strippedCreatedNodeLevels;
    std::vector<int> strippedCreatedNodePlanIndices;
    int strippedCreatedConnectorNodes = 0;
    int strippedCreatedHiddenExternalParentConnectorNodes = 0;
    int strippedCreatedOracleBridgeNodes = 0;
    int strippedCreatedStableChildBridgeNodes = 0;
    int strippedCreatedHiddenSupportConnectorNodes = 0;
    int strippedCreatedNonOracleHiddenSupportConnectorNodes = 0;
    int strippedCreatedNonOracleParentLinkBridgeNodes = 0;
    int strippedCreatedNonOracleChildLinkBridgeNodes = 0;
    int strippedNonOracleHiddenSupportCandidates = 0;
    int strippedNonOracleHiddenSupportAlreadyRepresented = 0;
    int strippedNonOracleHiddenSupportMaterializable = 0;
    int strippedNonOracleParentLinksTouchingPlan = 0;
    int strippedNonOracleParentLinksBothRepresented = 0;
    int strippedNonOracleParentLinksBlockedByModifiedSnapshot = 0;
    int strippedNonOracleParentLinksInvalidLevel = 0;
    int strippedNonOracleParentLinksAlreadyDirect = 0;
    int strippedNonOracleParentLinksAlreadyIntermediate = 0;
    int strippedNonOracleParentLinksMoved = 0;
    int strippedNonOracleParentLinkBridgeAttempts = 0;
    int strippedNonOracleParentLinkBridgeMaterializable = 0;
    int strippedNonOracleChildLinkBridgeAttempts = 0;
    int strippedNonOracleChildLinkBridgeInvalidLevel = 0;
    int strippedNonOracleChildLinkBridgeMaterializable = 0;
    int strippedLocalReleasedSourceOwners = 0;
    int strippedEmptyCreatedNodesAfterCanonicalization = 0;
    int strippedEmptyConnectorNodesAfterCanonicalization = 0;
    int strippedPostEmptyNodeReleasesRequired = 0;
    int strippedPostSameLevelParentChildCollapsesRequired = 0;
    int strippedPostSameLevelSiblingMergesRequired = 0;
    int strippedPostComparableSiblingMovesRequired = 0;
    int strippedPostCanonicalizationChangesRequired = 0;
    std::vector<LocalEmptyCreatedNodeDiagnostic> strippedEmptyCreatedNodeDiagnostics;
    std::string strippedLocalStructureMismatchSummary;
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
    std::string strippedLocalFailureSummary;
    bool oracleGuidedLocalApplied = false;
    bool oracleGuidedLocalMatchesTarget = false;
    bool oracleGuidedLocalStructureMatchesOracle = false;
    bool oracleGuidedDualMatchesPrimalAfterMerge = false;
    int oracleGuidedRecompositionNodes = 0;
    int oracleGuidedCreatedNodes = 0;
    int oracleGuidedReleasedSourceOwners = 0;
    std::string oracleGuidedStructureMismatchSummary;
};

struct NonIncreasingCandidateSignature {
    NodeId nodeId = InvalidNode;
    int altitude = 0;
    std::vector<PixelId> supportPixels;

    bool operator<(const NonIncreasingCandidateSignature &other) const {
        return std::tie(altitude, supportPixels) <
               std::tie(other.altitude, other.supportPixels);
    }

    bool operator==(const NonIncreasingCandidateSignature &other) const {
        return altitude == other.altitude && supportPixels == other.supportPixels;
    }
};

struct NonIncreasingCandidateOrderEntry {
    NodeId nodeId = InvalidNode;
    int altitude = 0;
    PixelId supportMinPixel = std::numeric_limits<PixelId>::max();
    int supportSize = 0;
};

struct NonIncreasingCandidateCollectionScratch {
    FastQueue<NodeId> queue;
    std::vector<NodeId> traversalOrder;
    std::vector<NodeId> scanOrder;
    std::vector<NonIncreasingCandidateOrderEntry> ordered;
    std::vector<NodeId> output;
    std::vector<int> supportSizeByNode;
    std::vector<PixelId> supportMinPixelByNode;
};

struct NonIncreasingCandidateHierarchyContext {
    NodeId nodeId = InvalidNode;
    int altitude = -1;
    int supportSize = 0;
    std::vector<PixelId> supportPixels;
    int minRow = -1;
    int minCol = -1;
    int maxRow = -1;
    int maxCol = -1;
    int parentAltitude = -1;
    int parentSupportSize = 0;
    int childCount = 0;
    int firstChildAltitude = -1;
    int firstChildSupportSize = 0;
};

struct NonIncreasingPhaseOracleComparison {
    float threshold = 0.0f;
    NonIncreasingPhase completedPhase = NonIncreasingPhase::MaxToMin;
    bool beforeMaxToMin = false;
    ImageUInt8Ptr localReconstructionImage;
    bool localReconstructionMatchesRebuild = false;
    bool candidateSetMatchesRebuild = false;
    bool candidateScanOrderMatchesRebuild = false;
    int localCandidateCount = 0;
    int rebuildCandidateCount = 0;
    int localOnlyCandidateCount = 0;
    int rebuildOnlyCandidateCount = 0;
    int firstCandidateScanOrderMismatchIndex = -1;
    int firstLocalOnlyAltitude = -1;
    int firstLocalOnlySupportSize = 0;
    NonIncreasingCandidateHierarchyContext firstLocalOnlyContext;
    int firstRebuildOnlyAltitude = -1;
    int firstRebuildOnlySupportSize = 0;
    NonIncreasingCandidateHierarchyContext firstRebuildOnlyContext;
    bool envelopeStructureCompared = false;
    bool envelopeStructureMatchesRebuild = false;
    int localEnvelopeNodeCount = 0;
    int rebuildEnvelopeNodeCount = 0;
    bool boundaryStructureCompared = false;
    bool boundaryStructureMatchesRebuild = false;
    bool boundaryIncidenceMatchesRebuild = false;
    int localBoundaryNodeCount = 0;
    int rebuildBoundaryNodeCount = 0;
    int localBoundaryIncidenceCount = 0;
    int rebuildBoundaryIncidenceCount = 0;
    bool releasePredicateEvaluated = false;
    bool releasePredicateSatisfied = false;
    bool releaseFocusAvailable = false;
    bool releaseH1ReconstructionMatches = false;
    bool releaseH2CandidateSetMatches = false;
    bool releaseH3EnvelopeStructureMatches = false;
    bool releaseH4BoundaryStructureMatches = false;
    bool releaseH4BoundaryIncidenceMatches = false;
    bool releaseH5CandidateScanOrderMatches = false;
    int releaseFocusNodeCount = 0;
    int releaseFocusPixelCount = 0;
    int releaseClosureNodeCount = 0;
    int releaseClosurePixelCount = 0;
    std::vector<PixelId> releaseFocusPixels;
    std::vector<PixelId> releaseClosurePixels;
};

struct NonIncreasingCasfDebugState {
    NonIncreasingMergeEvent firstInvariantMismatchEvent;
    LocalMergeOracle firstInvariantMismatchOracle;
    LocalMergeOracle capturedLocalOracle;
    std::vector<NonIncreasingMergeEvent> mergeEvents;
    std::vector<NonIncreasingPhaseOracleComparison> phaseOracleComparisons;
};

struct NonIncreasingCasfDebugOptionsState {
    bool collectMergeEvents = false;
    bool collectPhaseOracleComparisons = false;
    bool phaseCandidateEquivalenceGuardEnabled = false;
    bool stopOnFirstInvariantMismatch = false;
    bool stopOnFirstPhaseRebuildFallback = false;
    bool stopOnFirstPhaseReleaseRejection = false;
    bool collectFirstInvariantMismatchOracle = false;
    NonIncreasingLocalOracleCaptureTarget localOracleCaptureTarget;
    bool leafExplicitFastPathEnabled = true;
    bool mcpcDirectPreflightDiagnosticsEnabled = false;
    bool mcpcDirectV1StrictPrecheckEnabled = false;
    bool mcpcNoCopyDirectCommitEnabled = false;
    bool mcpcValidatedDirectCommitEnabled = false;
    bool operationalMergeEventInstrumentationEnabled = true;
    NonIncreasingStrippedLocalUpdateOptions strippedLocalUpdateOptions;
    NonIncreasingStrippedCommitMode strippedCommitMode = NonIncreasingStrippedCommitMode::SafeCopy;
    int maxAppliedMergeEventsPerThreshold = 0;
    long long maxThresholdMicros = 0;
};

struct NonIncreasingPhaseOracleMismatchFixture {
    bool found = false;
    float threshold = 0.0f;
    NonIncreasingPhase completedPhase = NonIncreasingPhase::MaxToMin;
    bool beforeMaxToMin = false;
    int minRow = -1;
    int minCol = -1;
    int maxRow = -1;
    int maxCol = -1;
    ImageUInt8Ptr croppedImage;
    NonIncreasingCandidateHierarchyContext localOnlyContext;
    NonIncreasingCandidateHierarchyContext rebuildOnlyContext;
    std::vector<PixelId> localOnlySupportInCrop;
    std::vector<PixelId> rebuildOnlySupportInCrop;
};

/**
 * @brief Incremental gray-level variance over the support of each tree node.
 *
 * The gray value associated with a direct proper part is the altitude of its
 * current owner node. A node attribute is therefore the variance of the
 * reconstructed gray levels over the node support.
 */
class DynamicGrayLevelVarianceComputer : public DynamicAttributeComputer {
private:
    DynamicComponentTree *tree_ = nullptr;
    mutable std::vector<double> count_;
    mutable std::vector<double> sum_;
    mutable std::vector<double> sumSquares_;

    void ensureStorage() const {
        const auto slots = static_cast<std::size_t>(tree_ ? tree_->getNumInternalNodeSlots() : 0);
        if (count_.size() < slots) {
            count_.resize(slots, 0.0);
            sum_.resize(slots, 0.0);
            sumSquares_.resize(slots, 0.0);
        }
    }

    void computeNode(NodeId nodeId, std::span<float> buffer) const {
        preProcessing(nodeId, buffer);
        for (NodeId childId : tree_->getChildren(nodeId)) {
            computeNode(childId, buffer);
            mergeProcessing(nodeId, childId, buffer);
        }
        postProcessing(nodeId, buffer);
    }

public:
    explicit DynamicGrayLevelVarianceComputer(DynamicComponentTree *tree)
        : tree_(tree),
          count_(static_cast<std::size_t>(tree ? tree->getNumInternalNodeSlots() : 0), 0.0),
          sum_(static_cast<std::size_t>(tree ? tree->getNumInternalNodeSlots() : 0), 0.0),
          sumSquares_(static_cast<std::size_t>(tree ? tree->getNumInternalNodeSlots() : 0), 0.0) {}

    void preProcessing(NodeId nodeId, std::span<float>) const override {
        ensureStorage();
        if (tree_ == nullptr || nodeId == InvalidNode || !tree_->isNode(nodeId) || !tree_->isAlive(nodeId)) {
            return;
        }

        const auto idx = static_cast<std::size_t>(nodeId);
        const double n = static_cast<double>(tree_->getNumProperParts(nodeId));
        const double level = static_cast<double>(tree_->getAltitude(nodeId));
        count_[idx] = n;
        sum_[idx] = n * level;
        sumSquares_[idx] = n * level * level;
    }

    void mergeProcessing(NodeId parentId, NodeId childId, std::span<float>) const override {
        ensureStorage();
        const auto parentIdx = static_cast<std::size_t>(parentId);
        const auto childIdx = static_cast<std::size_t>(childId);
        count_[parentIdx] += count_[childIdx];
        sum_[parentIdx] += sum_[childIdx];
        sumSquares_[parentIdx] += sumSquares_[childIdx];
    }

    void postProcessing(NodeId nodeId, std::span<float> buffer) const override {
        ensureStorage();
        if (nodeId == InvalidNode || static_cast<std::size_t>(nodeId) >= buffer.size()) {
            return;
        }

        const auto idx = static_cast<std::size_t>(nodeId);
        if (count_[idx] <= 1.0) {
            buffer[idx] = 0.0f;
            return;
        }

        const double mean = sum_[idx] / count_[idx];
        const double variance = std::max(0.0, (sumSquares_[idx] / count_[idx]) - (mean * mean));
        buffer[idx] = static_cast<float>(variance);
    }

    void onNodeRemoved(NodeId nodeId) const override {
        ensureStorage();
        if (nodeId == InvalidNode || static_cast<std::size_t>(nodeId) >= count_.size()) {
            return;
        }
        const auto idx = static_cast<std::size_t>(nodeId);
        count_[idx] = 0.0;
        sum_[idx] = 0.0;
        sumSquares_[idx] = 0.0;
    }

    std::vector<float> compute() const {
        std::vector<float> buffer(static_cast<std::size_t>(tree_ ? tree_->getNumInternalNodeSlots() : 0), 0.0f);
        compute(std::span<float>(buffer));
        return buffer;
    }

    void compute(std::span<float> buffer) const {
        if (tree_ == nullptr || tree_->getRoot() == InvalidNode) {
            return;
        }
        ensureStorage();
        computeNode(tree_->getRoot(), buffer);
    }
};

inline std::vector<float> computeGrayLevelVariance(DynamicComponentTree &tree) {
    DynamicGrayLevelVarianceComputer computer(&tree);
    return computer.compute();
}

inline std::vector<float> computeBoxAreaPerimeterRatio(DynamicComponentTree &tree) {
    DynamicBoundingBoxComputer computer(&tree, Attribute::AREA_PER_BOX_PERIMETER);
    return computer.compute();
}

inline std::vector<float> computeNonIncreasingAttribute(
    DynamicComponentTree &tree,
    NonIncreasingAttributeKind attributeKind) {
    if (attributeKind == NonIncreasingAttributeKind::BoxAreaPerimeterRatio) {
        return computeBoxAreaPerimeterRatio(tree);
    }
    return computeGrayLevelVariance(tree);
}

inline void collectNodeSupportPixelsInto(const DynamicComponentTree &tree,
                                         NodeId nodeId,
                                         std::vector<PixelId> &pixels) {
    pixels.clear();
    if (nodeId == InvalidNode || !tree.isAlive(nodeId)) {
        return;
    }
    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(pixelId);
        }
    }
    std::sort(pixels.begin(), pixels.end());
    pixels.erase(std::unique(pixels.begin(), pixels.end()), pixels.end());
}

inline std::vector<NodeId> collectNonIncreasingMergeCandidates(const DynamicComponentTree &tree,
                                                               const std::vector<float> &attribute,
                                                               float threshold,
                                                               NonIncreasingCandidateCollectionScratch *scratch) {
    NonIncreasingCandidateCollectionScratch localScratch;
    auto &work = scratch == nullptr ? localScratch : *scratch;

    work.queue.clear();
    work.traversalOrder.clear();
    work.scanOrder.clear();
    work.queue.push(tree.getRoot());
    while (!work.queue.empty()) {
        const NodeId nodeId = work.queue.pop();
        if (nodeId == InvalidNode ||
            !tree.isNode(nodeId) ||
            !tree.isAlive(nodeId)) {
            continue;
        }
        work.traversalOrder.push_back(nodeId);
        if (nodeId != tree.getRoot() &&
            static_cast<std::size_t>(nodeId) < attribute.size() &&
            attribute[static_cast<std::size_t>(nodeId)] <= threshold) {
            work.scanOrder.push_back(nodeId);
        }

        for (NodeId childId : tree.getChildren(nodeId)) {
            work.queue.push(childId);
        }
    }

    const std::size_t nodeSlots =
        static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    if (work.supportSizeByNode.size() < nodeSlots) {
        work.supportSizeByNode.resize(nodeSlots, 0);
    }
    if (work.supportMinPixelByNode.size() < nodeSlots) {
        work.supportMinPixelByNode.resize(
            nodeSlots,
            std::numeric_limits<PixelId>::max());
    }

    for (auto it = work.traversalOrder.rbegin(); it != work.traversalOrder.rend(); ++it) {
        const NodeId nodeId = *it;
        const auto index = static_cast<std::size_t>(nodeId);
        int supportSize = tree.getNumProperParts(nodeId);
        PixelId supportMinPixel = std::numeric_limits<PixelId>::max();
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            supportMinPixel = std::min(supportMinPixel, pixelId);
        }
        for (NodeId childId : tree.getChildren(nodeId)) {
            if (childId == InvalidNode ||
                childId < 0 ||
                static_cast<std::size_t>(childId) >= work.supportSizeByNode.size()) {
                continue;
            }
            const auto childIndex = static_cast<std::size_t>(childId);
            supportSize += work.supportSizeByNode[childIndex];
            supportMinPixel =
                std::min(supportMinPixel, work.supportMinPixelByNode[childIndex]);
        }
        work.supportSizeByNode[index] = supportSize;
        work.supportMinPixelByNode[index] = supportMinPixel;
    }

    work.ordered.resize(work.scanOrder.size());
    for (std::size_t index = 0; index < work.scanOrder.size(); ++index) {
        const NodeId nodeId = work.scanOrder[index];
        auto &entry = work.ordered[index];
        entry.nodeId = nodeId;
        entry.altitude = tree.getAltitude(nodeId);
        const auto nodeIndex = static_cast<std::size_t>(nodeId);
        entry.supportMinPixel = work.supportMinPixelByNode[nodeIndex];
        entry.supportSize = work.supportSizeByNode[nodeIndex];
    }
    std::sort(work.ordered.begin(),
              work.ordered.end(),
              [](const NonIncreasingCandidateOrderEntry &lhs,
                 const NonIncreasingCandidateOrderEntry &rhs) {
                  return std::tie(lhs.altitude, lhs.supportMinPixel, lhs.supportSize, lhs.nodeId) <
                         std::tie(rhs.altitude, rhs.supportMinPixel, rhs.supportSize, rhs.nodeId);
              });

    work.output.clear();
    work.output.reserve(work.ordered.size());
    for (const auto &entry : work.ordered) {
        work.output.push_back(entry.nodeId);
    }
    return work.output;
}

inline std::vector<NodeId> collectNonIncreasingMergeCandidates(const DynamicComponentTree &tree,
                                                               const std::vector<float> &attribute,
                                                               float threshold) {
    return collectNonIncreasingMergeCandidates(tree, attribute, threshold, nullptr);
}

inline NonIncreasingCandidateSignature describeNonIncreasingCandidateSignature(
    const DynamicComponentTree &tree,
    NodeId nodeId) {
    NonIncreasingCandidateSignature signature;
    signature.nodeId = nodeId;
    signature.altitude = tree.getAltitude(nodeId);
    signature.supportPixels = collectNodeSupportPixels(tree, nodeId);
    return signature;
}

inline std::vector<NonIncreasingCandidateSignature> collectNonIncreasingCandidateScanSignatures(
    const DynamicComponentTree &tree,
    const std::vector<float> &attribute,
    float threshold) {
    std::vector<NonIncreasingCandidateSignature> signatures;
    const auto candidates = collectNonIncreasingMergeCandidates(tree, attribute, threshold);
    signatures.reserve(candidates.size());
    for (NodeId nodeId : candidates) {
        signatures.push_back(describeNonIncreasingCandidateSignature(tree, nodeId));
    }
    return signatures;
}

inline std::vector<NonIncreasingCandidateSignature> collectNonIncreasingCandidateSignatures(
    const DynamicComponentTree &tree,
    const std::vector<float> &attribute,
    float threshold) {
    auto signatures = collectNonIncreasingCandidateScanSignatures(tree, attribute, threshold);
    std::sort(signatures.begin(), signatures.end());
    return signatures;
}

inline NonIncreasingCandidateHierarchyContext describeNonIncreasingCandidateHierarchy(
    const DynamicComponentTree &tree,
    const NonIncreasingCandidateSignature &signature) {
    NonIncreasingCandidateHierarchyContext context;
    context.nodeId = signature.nodeId;
    context.altitude = signature.altitude;
    context.supportPixels = signature.supportPixels;
    context.supportSize = static_cast<int>(signature.supportPixels.size());
    const int numCols = tree.getNumColsOfImage();
    for (PixelId pixelId : context.supportPixels) {
        if (pixelId < 0 || numCols <= 0) {
            continue;
        }
        const int row = pixelId / numCols;
        const int col = pixelId % numCols;
        if (context.minRow < 0 || row < context.minRow) {
            context.minRow = row;
        }
        if (context.maxRow < 0 || row > context.maxRow) {
            context.maxRow = row;
        }
        if (context.minCol < 0 || col < context.minCol) {
            context.minCol = col;
        }
        if (context.maxCol < 0 || col > context.maxCol) {
            context.maxCol = col;
        }
    }

    const NodeId nodeId = signature.nodeId;
    if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
        return context;
    }

    const NodeId parentId = tree.getNodeParent(nodeId);
    if (parentId != InvalidNode && parentId != nodeId && tree.isNode(parentId) && tree.isAlive(parentId)) {
        context.parentAltitude = tree.getAltitude(parentId);
        context.parentSupportSize = static_cast<int>(collectNodeSupportPixels(tree, parentId).size());
    }

    context.childCount = tree.getNumChildren(nodeId);
    std::vector<std::pair<int, int>> childSummaries;
    childSummaries.reserve(static_cast<std::size_t>(context.childCount));
    for (NodeId childId : tree.getChildren(nodeId)) {
        if (childId == InvalidNode || !tree.isNode(childId) || !tree.isAlive(childId)) {
            continue;
        }
        childSummaries.push_back({
            tree.getAltitude(childId),
            static_cast<int>(collectNodeSupportPixels(tree, childId).size())});
    }
    std::sort(childSummaries.begin(), childSummaries.end());
    if (!childSummaries.empty()) {
        context.firstChildAltitude = childSummaries.front().first;
        context.firstChildSupportSize = childSummaries.front().second;
    }
    return context;
}

inline void addPhaseCandidateGuardFocusNode(const DynamicComponentTree &tree,
                                            std::vector<NodeId> &focusNodes,
                                            NodeId nodeId,
                                            int ancestorLimit = 6) {
    NodeId current = nodeId;
    int depth = 0;
    while (current != InvalidNode &&
           tree.isNode(current) &&
           tree.isAlive(current) &&
           depth <= ancestorLimit) {
        focusNodes.push_back(current);
        if (tree.isRoot(current)) {
            break;
        }
        const NodeId parentId = tree.getNodeParent(current);
        if (parentId == InvalidNode || parentId == current) {
            break;
        }
        current = parentId;
        ++depth;
    }
}

inline void addPhaseCandidateGuardContextFocus(
    const DynamicComponentTree &localTree,
    const NonIncreasingCandidateHierarchyContext &context,
    bool includeContextNodeId,
    std::vector<NodeId> &focusNodes,
    std::vector<PixelId> &localPixels) {
    if (context.supportPixels.empty()) {
        return;
    }
    if (includeContextNodeId) {
        addPhaseCandidateGuardFocusNode(localTree, focusNodes, context.nodeId);
    }
    for (PixelId pixelId : context.supportPixels) {
        if (pixelId < 0 || pixelId >= localTree.getNumTotalProperParts()) {
            continue;
        }
        localPixels.push_back(pixelId);
        addPhaseCandidateGuardFocusNode(
            localTree,
            focusNodes,
            localTree.getSmallestComponent(pixelId));
    }
}

inline std::pair<std::vector<NodeId>, std::vector<PixelId>>
buildPhaseCandidateGuardFocus(const DynamicComponentTree &localTree,
                              const NonIncreasingPhaseOracleComparison &comparison) {
    std::vector<NodeId> focusNodes;
    std::vector<PixelId> localPixels;
    addPhaseCandidateGuardContextFocus(
        localTree,
        comparison.firstLocalOnlyContext,
        true,
        focusNodes,
        localPixels);
    addPhaseCandidateGuardContextFocus(
        localTree,
        comparison.firstRebuildOnlyContext,
        false,
        focusNodes,
        localPixels);
    return {
        sortedUniqueValidNodes(std::move(focusNodes)),
        sortedUniquePixels(std::move(localPixels))};
}

inline std::vector<PixelId> collectSupportPixelsForNodes(
    const DynamicComponentTree &tree,
    std::span<const NodeId> nodes) {
    std::vector<PixelId> pixels;
    for (NodeId nodeId : nodes) {
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        const auto support = collectNodeSupportPixels(tree, nodeId);
        pixels.insert(pixels.end(), support.begin(), support.end());
    }
    return sortedUniquePixels(std::move(pixels));
}

struct NonIncreasingTreeNodeStructureSignature {
    int altitude = 0;
    std::vector<PixelId> supportPixels;
    std::vector<PixelId> directProperParts;
    std::vector<std::pair<int, std::vector<PixelId>>> children;

    bool operator<(const NonIncreasingTreeNodeStructureSignature &other) const {
        return std::tie(altitude, supportPixels, directProperParts, children) <
               std::tie(other.altitude, other.supportPixels, other.directProperParts, other.children);
    }

    bool operator==(const NonIncreasingTreeNodeStructureSignature &other) const {
        return altitude == other.altitude &&
               supportPixels == other.supportPixels &&
               directProperParts == other.directProperParts &&
               children == other.children;
    }
};

enum class NonIncreasingBoundaryRelation {
    Parent = 0,
    Child = 1,
    AdjacentSupport = 2,
};

struct NonIncreasingBoundaryIncidenceSignature {
    NonIncreasingTreeNodeStructureSignature envelopeNode;
    NonIncreasingTreeNodeStructureSignature boundaryNode;
    NonIncreasingBoundaryRelation relation = NonIncreasingBoundaryRelation::Parent;

    bool operator<(const NonIncreasingBoundaryIncidenceSignature &other) const {
        if (envelopeNode < other.envelopeNode) {
            return true;
        }
        if (other.envelopeNode < envelopeNode) {
            return false;
        }
        if (boundaryNode < other.boundaryNode) {
            return true;
        }
        if (other.boundaryNode < boundaryNode) {
            return false;
        }
        return static_cast<int>(relation) < static_cast<int>(other.relation);
    }

    bool operator==(const NonIncreasingBoundaryIncidenceSignature &other) const {
        return envelopeNode == other.envelopeNode &&
               boundaryNode == other.boundaryNode &&
               relation == other.relation;
    }
};

inline NonIncreasingTreeNodeStructureSignature describeNonIncreasingTreeNodeStructure(
    const DynamicComponentTree &tree,
    NodeId nodeId) {
    NonIncreasingTreeNodeStructureSignature signature;
    signature.altitude = tree.getAltitude(nodeId);
    signature.supportPixels = collectNodeSupportPixels(tree, nodeId);
    signature.directProperParts = collectDirectProperParts(tree, nodeId);
    for (NodeId childId : tree.getChildren(nodeId)) {
        if (childId == InvalidNode || !tree.isNode(childId) || !tree.isAlive(childId)) {
            continue;
        }
        signature.children.push_back({
            tree.getAltitude(childId),
            collectNodeSupportPixels(tree, childId)});
    }
    std::sort(signature.children.begin(), signature.children.end());
    return signature;
}

inline std::vector<NonIncreasingTreeNodeStructureSignature>
collectNonIncreasingTreeStructureSignaturesForNodes(const DynamicComponentTree &tree,
                                                    std::span<const NodeId> nodes) {
    std::vector<NonIncreasingTreeNodeStructureSignature> signatures;
    signatures.reserve(nodes.size());
    for (NodeId nodeId : nodes) {
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        signatures.push_back(describeNonIncreasingTreeNodeStructure(tree, nodeId));
    }
    std::sort(signatures.begin(), signatures.end());
    return signatures;
}

inline std::vector<NonIncreasingTreeNodeStructureSignature>
collectNonIncreasingTreeStructureSignatures(const DynamicComponentTree &tree) {
    std::vector<NonIncreasingTreeNodeStructureSignature> signatures;
    const NodeId root = tree.getRoot();
    if (root == InvalidNode || !tree.isNode(root) || !tree.isAlive(root)) {
        return signatures;
    }
    for (NodeId nodeId : tree.getNodeSubtree(root)) {
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        signatures.push_back(describeNonIncreasingTreeNodeStructure(tree, nodeId));
    }
    std::sort(signatures.begin(), signatures.end());
    return signatures;
}

inline bool nonIncreasingTreeStructureSignaturesMatch(const DynamicComponentTree &localTree,
                                                      const DynamicComponentTree &rebuiltTree) {
    return collectNonIncreasingTreeStructureSignatures(localTree) ==
           collectNonIncreasingTreeStructureSignatures(rebuiltTree);
}

inline void addNonIncreasingStructureEnvelopeNode(const DynamicComponentTree &tree,
                                                 std::vector<NodeId> &nodes,
                                                 NodeId nodeId,
                                                 int ancestorLimit = 6) {
    NodeId current = nodeId;
    int depth = 0;
    while (current != InvalidNode &&
           tree.isNode(current) &&
           tree.isAlive(current) &&
           depth <= ancestorLimit) {
        nodes.push_back(current);
        for (NodeId childId : tree.getChildren(current)) {
            if (childId != InvalidNode && tree.isNode(childId) && tree.isAlive(childId)) {
                nodes.push_back(childId);
            }
        }
        if (tree.isRoot(current)) {
            break;
        }
        const NodeId parentId = tree.getNodeParent(current);
        if (parentId == InvalidNode || parentId == current) {
            break;
        }
        current = parentId;
        ++depth;
    }
}

inline std::vector<NodeId> collectNonIncreasingStructureEnvelopeNodes(
    const DynamicComponentTree &tree,
    std::span<const NodeId> focusNodes,
    std::span<const PixelId> localPixels,
    int ancestorLimit = 6) {
    std::vector<NodeId> nodes;
    for (NodeId nodeId : focusNodes) {
        addNonIncreasingStructureEnvelopeNode(tree, nodes, nodeId, ancestorLimit);
    }
    for (PixelId pixelId : localPixels) {
        if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
            continue;
        }
        addNonIncreasingStructureEnvelopeNode(
            tree,
            nodes,
            tree.getSmallestComponent(pixelId),
            ancestorLimit);
    }

    nodes = sortedUniqueValidNodes(std::move(nodes));
    nodes.erase(
        std::remove_if(nodes.begin(),
                       nodes.end(),
                       [&](NodeId nodeId) {
                           return nodeId == InvalidNode ||
                                  !tree.isNode(nodeId) ||
                                  !tree.isAlive(nodeId);
                       }),
        nodes.end());
    return nodes;
}

inline std::vector<NonIncreasingTreeNodeStructureSignature>
collectNonIncreasingStructureEnvelopeSignatures(const DynamicComponentTree &tree,
                                                std::span<const NodeId> focusNodes,
                                                std::span<const PixelId> localPixels,
                                                int ancestorLimit = 6) {
    const auto nodes = collectNonIncreasingStructureEnvelopeNodes(
        tree,
        focusNodes,
        localPixels,
        ancestorLimit);
    return collectNonIncreasingTreeStructureSignaturesForNodes(
        tree,
        std::span<const NodeId>(nodes));
}

inline bool nonIncreasingStructureEnvelopeSignaturesMatch(const DynamicComponentTree &localTree,
                                                         const DynamicComponentTree &rebuiltTree,
                                                         std::span<const NodeId> localFocusNodes,
                                                         std::span<const PixelId> localPixels,
                                                         int ancestorLimit = 6) {
    return collectNonIncreasingStructureEnvelopeSignatures(
               localTree,
               localFocusNodes,
               localPixels,
               ancestorLimit) ==
           collectNonIncreasingStructureEnvelopeSignatures(
               rebuiltTree,
               std::span<const NodeId>{},
               localPixels,
               ancestorLimit);
}

inline bool nonIncreasingPhaseEnvelopeStructureSignaturesMatch(
    const DynamicComponentTree &localTree,
    const DynamicComponentTree &rebuiltTree,
    const NonIncreasingPhaseOracleComparison &comparison,
    int ancestorLimit = 6) {
    const auto [localFocusNodes, localPixels] =
        buildPhaseCandidateGuardFocus(localTree, comparison);
    return nonIncreasingStructureEnvelopeSignaturesMatch(
        localTree,
        rebuiltTree,
        std::span<const NodeId>(localFocusNodes),
        std::span<const PixelId>(localPixels),
        ancestorLimit);
}

inline bool containsSortedNode(std::span<const NodeId> nodes, NodeId nodeId) {
    return nodeId != InvalidNode &&
           std::binary_search(nodes.begin(), nodes.end(), nodeId);
}

inline bool nodeOrAncestorInSortedNodes(const DynamicComponentTree &tree,
                                        NodeId nodeId,
                                        std::span<const NodeId> nodes) {
    NodeId current = nodeId;
    while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
        if (containsSortedNode(nodes, current)) {
            return true;
        }
        const NodeId parentId = tree.getNodeParent(current);
        if (parentId == InvalidNode || parentId == current) {
            break;
        }
        current = parentId;
    }
    return false;
}

inline std::vector<NodeId> collectNonIncreasingStructureBoundaryNodes(
    const DynamicComponentTree &tree,
    std::span<const NodeId> envelopeNodes,
    AdjacencyRelation *adjacency = nullptr) {
    std::vector<NodeId> boundaryNodes;
    for (NodeId nodeId : envelopeNodes) {
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId != InvalidNode &&
            parentId != nodeId &&
            tree.isNode(parentId) &&
            tree.isAlive(parentId) &&
            !containsSortedNode(envelopeNodes, parentId)) {
            boundaryNodes.push_back(parentId);
        }
        for (NodeId childId : tree.getChildren(nodeId)) {
            if (childId != InvalidNode &&
                tree.isNode(childId) &&
                tree.isAlive(childId) &&
                !containsSortedNode(envelopeNodes, childId)) {
                boundaryNodes.push_back(childId);
            }
        }
    }

    if (adjacency != nullptr) {
        for (NodeId nodeId : envelopeNodes) {
            if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
                continue;
            }
            for (PixelId pixelId : collectNodeSupportPixels(tree, nodeId)) {
                for (PixelId neighbor : adjacency->getNeighborPixels(pixelId)) {
                    if (neighbor < 0 || neighbor >= tree.getNumTotalProperParts()) {
                        continue;
                    }
                    const NodeId ownerId = tree.getSmallestComponent(neighbor);
                    if (ownerId == InvalidNode ||
                        !tree.isNode(ownerId) ||
                        !tree.isAlive(ownerId) ||
                        nodeOrAncestorInSortedNodes(tree, ownerId, envelopeNodes)) {
                        continue;
                    }
                    boundaryNodes.push_back(ownerId);
                }
            }
        }
    }

    boundaryNodes = sortedUniqueValidNodes(std::move(boundaryNodes));
    boundaryNodes.erase(
        std::remove_if(boundaryNodes.begin(),
                       boundaryNodes.end(),
                       [&](NodeId nodeId) {
                           return nodeId == InvalidNode ||
                                  !tree.isNode(nodeId) ||
                                  !tree.isAlive(nodeId) ||
                                  containsSortedNode(envelopeNodes, nodeId);
                       }),
        boundaryNodes.end());
    return boundaryNodes;
}

inline std::vector<NonIncreasingTreeNodeStructureSignature>
collectNonIncreasingStructureBoundarySignatures(const DynamicComponentTree &tree,
                                                std::span<const NodeId> envelopeNodes,
                                                AdjacencyRelation *adjacency = nullptr) {
    const auto boundaryNodes = collectNonIncreasingStructureBoundaryNodes(
        tree,
        envelopeNodes,
        adjacency);
    return collectNonIncreasingTreeStructureSignaturesForNodes(
        tree,
        std::span<const NodeId>(boundaryNodes));
}

inline void addNonIncreasingBoundaryIncidence(
    const DynamicComponentTree &tree,
    std::vector<NonIncreasingBoundaryIncidenceSignature> &incidences,
    NodeId envelopeNodeId,
    NodeId boundaryNodeId,
    NonIncreasingBoundaryRelation relation) {
    if (envelopeNodeId == InvalidNode ||
        boundaryNodeId == InvalidNode ||
        !tree.isNode(envelopeNodeId) ||
        !tree.isNode(boundaryNodeId) ||
        !tree.isAlive(envelopeNodeId) ||
        !tree.isAlive(boundaryNodeId)) {
        return;
    }
    NonIncreasingBoundaryIncidenceSignature incidence;
    incidence.envelopeNode = describeNonIncreasingTreeNodeStructure(tree, envelopeNodeId);
    incidence.boundaryNode = describeNonIncreasingTreeNodeStructure(tree, boundaryNodeId);
    incidence.relation = relation;
    incidences.push_back(std::move(incidence));
}

inline std::vector<NonIncreasingBoundaryIncidenceSignature>
collectNonIncreasingBoundaryIncidenceSignatures(const DynamicComponentTree &tree,
                                                std::span<const NodeId> envelopeNodes,
                                                std::span<const NodeId> boundaryNodes,
                                                AdjacencyRelation *adjacency = nullptr) {
    std::vector<NonIncreasingBoundaryIncidenceSignature> incidences;
    for (NodeId nodeId : envelopeNodes) {
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (containsSortedNode(boundaryNodes, parentId)) {
            addNonIncreasingBoundaryIncidence(
                tree,
                incidences,
                nodeId,
                parentId,
                NonIncreasingBoundaryRelation::Parent);
        }
        for (NodeId childId : tree.getChildren(nodeId)) {
            if (containsSortedNode(boundaryNodes, childId)) {
                addNonIncreasingBoundaryIncidence(
                    tree,
                    incidences,
                    nodeId,
                    childId,
                    NonIncreasingBoundaryRelation::Child);
            }
        }
    }

    if (adjacency != nullptr) {
        for (NodeId nodeId : envelopeNodes) {
            if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
                continue;
            }
            for (PixelId pixelId : collectNodeSupportPixels(tree, nodeId)) {
                for (PixelId neighbor : adjacency->getNeighborPixels(pixelId)) {
                    if (neighbor < 0 || neighbor >= tree.getNumTotalProperParts()) {
                        continue;
                    }
                    const NodeId ownerId = tree.getSmallestComponent(neighbor);
                    if (containsSortedNode(boundaryNodes, ownerId)) {
                        addNonIncreasingBoundaryIncidence(
                            tree,
                            incidences,
                            nodeId,
                            ownerId,
                            NonIncreasingBoundaryRelation::AdjacentSupport);
                    }
                }
            }
        }
    }

    std::sort(incidences.begin(), incidences.end());
    incidences.erase(std::unique(incidences.begin(), incidences.end()), incidences.end());
    return incidences;
}

inline NonIncreasingPhaseOracleComparison compareNonIncreasingPhaseAgainstRebuild(
    const DynamicComponentTree &localTree,
    DynamicComponentTree &rebuiltTree,
    float threshold,
    NonIncreasingPhase completedPhase,
    bool beforeMaxToMin = false,
    AdjacencyRelation *adjacency = nullptr,
    NonIncreasingAttributeKind attributeKind = NonIncreasingAttributeKind::GrayVariance) {
    NonIncreasingPhaseOracleComparison comparison;
    comparison.threshold = threshold;
    comparison.completedPhase = completedPhase;
    comparison.beforeMaxToMin = beforeMaxToMin;
    comparison.localReconstructionMatchesRebuild =
        localTree.reconstructionImage()->isEqual(rebuiltTree.reconstructionImage());

    const auto localAttribute = computeNonIncreasingAttribute(
        *const_cast<DynamicComponentTree *>(&localTree),
        attributeKind);
    const auto rebuiltAttribute = computeNonIncreasingAttribute(
        rebuiltTree,
        attributeKind);
    const auto localScanSignatures = collectNonIncreasingCandidateScanSignatures(
        localTree,
        localAttribute,
        threshold);
    const auto rebuiltScanSignatures = collectNonIncreasingCandidateScanSignatures(
        rebuiltTree,
        rebuiltAttribute,
        threshold);
    const auto localSignatures = collectNonIncreasingCandidateSignatures(
        localTree,
        localAttribute,
        threshold);
    const auto rebuiltSignatures = collectNonIncreasingCandidateSignatures(
        rebuiltTree,
        rebuiltAttribute,
        threshold);

    comparison.localCandidateCount = static_cast<int>(localSignatures.size());
    comparison.rebuildCandidateCount = static_cast<int>(rebuiltSignatures.size());
    comparison.candidateSetMatchesRebuild = localSignatures == rebuiltSignatures;
    comparison.candidateScanOrderMatchesRebuild =
        localScanSignatures == rebuiltScanSignatures;
    comparison.releaseH1ReconstructionMatches =
        comparison.localReconstructionMatchesRebuild;
    comparison.releaseH2CandidateSetMatches =
        comparison.candidateSetMatchesRebuild;
    comparison.releaseH5CandidateScanOrderMatches =
        comparison.candidateScanOrderMatchesRebuild;
    const std::size_t maxScanSignatureCount =
        std::max(localScanSignatures.size(), rebuiltScanSignatures.size());
    for (std::size_t i = 0; i < maxScanSignatureCount; ++i) {
        if (i >= localScanSignatures.size() ||
            i >= rebuiltScanSignatures.size() ||
            !(localScanSignatures[i] == rebuiltScanSignatures[i])) {
            comparison.firstCandidateScanOrderMismatchIndex = static_cast<int>(i);
            break;
        }
    }

    std::size_t localIndex = 0;
    std::size_t rebuiltIndex = 0;
    while (localIndex < localSignatures.size() || rebuiltIndex < rebuiltSignatures.size()) {
        if (rebuiltIndex == rebuiltSignatures.size() ||
            (localIndex < localSignatures.size() && localSignatures[localIndex] < rebuiltSignatures[rebuiltIndex])) {
            ++comparison.localOnlyCandidateCount;
            if (comparison.firstLocalOnlyAltitude < 0) {
                comparison.firstLocalOnlyContext =
                    describeNonIncreasingCandidateHierarchy(localTree, localSignatures[localIndex]);
                comparison.firstLocalOnlyAltitude = comparison.firstLocalOnlyContext.altitude;
                comparison.firstLocalOnlySupportSize = comparison.firstLocalOnlyContext.supportSize;
            }
            ++localIndex;
        } else if (localIndex == localSignatures.size() ||
                   rebuiltSignatures[rebuiltIndex] < localSignatures[localIndex]) {
            ++comparison.rebuildOnlyCandidateCount;
            if (comparison.firstRebuildOnlyAltitude < 0) {
                comparison.firstRebuildOnlyContext =
                    describeNonIncreasingCandidateHierarchy(rebuiltTree, rebuiltSignatures[rebuiltIndex]);
                comparison.firstRebuildOnlyAltitude = comparison.firstRebuildOnlyContext.altitude;
                comparison.firstRebuildOnlySupportSize = comparison.firstRebuildOnlyContext.supportSize;
            }
            ++rebuiltIndex;
        } else {
            ++localIndex;
            ++rebuiltIndex;
        }
    }

    if (!comparison.candidateSetMatchesRebuild ||
        !comparison.localReconstructionMatchesRebuild) {
        const auto [localFocusNodes, localPixels] =
            buildPhaseCandidateGuardFocus(localTree, comparison);
        if (!localFocusNodes.empty() || !localPixels.empty()) {
            const auto localEnvelopeNodes =
                collectNonIncreasingStructureEnvelopeNodes(
                    localTree,
                    std::span<const NodeId>(localFocusNodes),
                    std::span<const PixelId>(localPixels));
            const auto rebuildEnvelopeNodes =
                collectNonIncreasingStructureEnvelopeNodes(
                    rebuiltTree,
                    std::span<const NodeId>{},
                    std::span<const PixelId>(localPixels));
            comparison.releaseFocusNodeCount =
                static_cast<int>(localFocusNodes.size());
            comparison.releaseFocusPixels = localPixels;
            comparison.releaseFocusPixelCount =
                static_cast<int>(comparison.releaseFocusPixels.size());
            comparison.releaseClosureNodeCount =
                static_cast<int>(localEnvelopeNodes.size());
            comparison.releaseClosurePixels = collectSupportPixelsForNodes(
                localTree,
                std::span<const NodeId>(localEnvelopeNodes));
            comparison.releaseClosurePixelCount =
                static_cast<int>(comparison.releaseClosurePixels.size());
            const auto localEnvelopeSignatures =
                collectNonIncreasingTreeStructureSignaturesForNodes(
                    localTree,
                    std::span<const NodeId>(localEnvelopeNodes));
            const auto rebuildEnvelopeSignatures =
                collectNonIncreasingTreeStructureSignaturesForNodes(
                    rebuiltTree,
                    std::span<const NodeId>(rebuildEnvelopeNodes));
            comparison.envelopeStructureCompared = true;
            comparison.envelopeStructureMatchesRebuild =
                localEnvelopeSignatures == rebuildEnvelopeSignatures;
            comparison.localEnvelopeNodeCount =
                static_cast<int>(localEnvelopeSignatures.size());
            comparison.rebuildEnvelopeNodeCount =
                static_cast<int>(rebuildEnvelopeSignatures.size());
            const auto localBoundaryNodes =
                collectNonIncreasingStructureBoundaryNodes(
                    localTree,
                    std::span<const NodeId>(localEnvelopeNodes),
                    adjacency);
            const auto rebuildBoundaryNodes =
                collectNonIncreasingStructureBoundaryNodes(
                    rebuiltTree,
                    std::span<const NodeId>(rebuildEnvelopeNodes),
                    adjacency);
            const auto localBoundarySignatures =
                collectNonIncreasingTreeStructureSignaturesForNodes(
                    localTree,
                    std::span<const NodeId>(localBoundaryNodes));
            const auto rebuildBoundarySignatures =
                collectNonIncreasingTreeStructureSignaturesForNodes(
                    rebuiltTree,
                    std::span<const NodeId>(rebuildBoundaryNodes));
            const auto localBoundaryIncidences =
                collectNonIncreasingBoundaryIncidenceSignatures(
                    localTree,
                    std::span<const NodeId>(localEnvelopeNodes),
                    std::span<const NodeId>(localBoundaryNodes),
                    adjacency);
            const auto rebuildBoundaryIncidences =
                collectNonIncreasingBoundaryIncidenceSignatures(
                    rebuiltTree,
                    std::span<const NodeId>(rebuildEnvelopeNodes),
                    std::span<const NodeId>(rebuildBoundaryNodes),
                    adjacency);
            comparison.boundaryStructureCompared = true;
            comparison.boundaryStructureMatchesRebuild =
                localBoundarySignatures == rebuildBoundarySignatures;
            comparison.boundaryIncidenceMatchesRebuild =
                localBoundaryIncidences == rebuildBoundaryIncidences;
            comparison.localBoundaryNodeCount =
                static_cast<int>(localBoundarySignatures.size());
            comparison.rebuildBoundaryNodeCount =
                static_cast<int>(rebuildBoundarySignatures.size());
            comparison.localBoundaryIncidenceCount =
                static_cast<int>(localBoundaryIncidences.size());
            comparison.rebuildBoundaryIncidenceCount =
                static_cast<int>(rebuildBoundaryIncidences.size());
            comparison.releaseFocusAvailable = true;
            comparison.releasePredicateEvaluated = true;
            comparison.releaseH3EnvelopeStructureMatches =
                comparison.envelopeStructureMatchesRebuild;
            comparison.releaseH4BoundaryStructureMatches =
                comparison.boundaryStructureMatchesRebuild;
            comparison.releaseH4BoundaryIncidenceMatches =
                comparison.boundaryIncidenceMatchesRebuild;
            comparison.releasePredicateSatisfied =
                comparison.releaseH1ReconstructionMatches &&
                comparison.releaseH2CandidateSetMatches &&
                comparison.releaseH3EnvelopeStructureMatches &&
                comparison.releaseH4BoundaryStructureMatches &&
                comparison.releaseH4BoundaryIncidenceMatches &&
                comparison.releaseH5CandidateScanOrderMatches;
        }
        comparison.localReconstructionImage = localTree.reconstructionImage();
    }

    return comparison;
}

inline void includeCandidateContextInWindow(const NonIncreasingCandidateHierarchyContext &context,
                                            int &minRow,
                                            int &minCol,
                                            int &maxRow,
                                            int &maxCol) {
    if (context.supportPixels.empty() ||
        context.minRow < 0 ||
        context.minCol < 0 ||
        context.maxRow < context.minRow ||
        context.maxCol < context.minCol) {
        return;
    }

    if (minRow < 0 || context.minRow < minRow) {
        minRow = context.minRow;
    }
    if (minCol < 0 || context.minCol < minCol) {
        minCol = context.minCol;
    }
    if (maxRow < 0 || context.maxRow > maxRow) {
        maxRow = context.maxRow;
    }
    if (maxCol < 0 || context.maxCol > maxCol) {
        maxCol = context.maxCol;
    }
}

inline std::vector<PixelId> remapSupportToCrop(std::span<const PixelId> support,
                                               int originalNumCols,
                                               int cropMinRow,
                                               int cropMinCol,
                                               int cropNumCols) {
    std::vector<PixelId> remapped;
    remapped.reserve(support.size());
    for (PixelId pixelId : support) {
        if (pixelId < 0 || originalNumCols <= 0) {
            continue;
        }
        const int row = pixelId / originalNumCols;
        const int col = pixelId % originalNumCols;
        const int cropRow = row - cropMinRow;
        const int cropCol = col - cropMinCol;
        if (cropRow < 0 || cropCol < 0) {
            continue;
        }
        remapped.push_back(cropRow * cropNumCols + cropCol);
    }
    return sortedUniquePixels(std::move(remapped));
}

inline ImageUInt8Ptr cropImageWindow(const ImageUInt8Ptr &image,
                                     int minRow,
                                     int minCol,
                                     int maxRow,
                                     int maxCol) {
    if (!image || minRow < 0 || minCol < 0 || maxRow < minRow || maxCol < minCol) {
        return nullptr;
    }
    const int cropRows = maxRow - minRow + 1;
    const int cropCols = maxCol - minCol + 1;
    auto crop = ImageUInt8::create(cropRows, cropCols);
    for (int row = 0; row < cropRows; ++row) {
        for (int col = 0; col < cropCols; ++col) {
            const int src = (minRow + row) * image->getNumCols() + (minCol + col);
            const int dst = row * cropCols + col;
            (*crop)[dst] = (*image)[src];
        }
    }
    return crop;
}

inline NonIncreasingPhaseOracleMismatchFixture buildPhaseOracleMismatchFixture(
    const NonIncreasingPhaseOracleComparison &comparison,
    int padding = 1) {
    if (comparison.candidateSetMatchesRebuild || !comparison.localReconstructionImage) {
        return {};
    }

    int minRow = -1;
    int minCol = -1;
    int maxRow = -1;
    int maxCol = -1;
    includeCandidateContextInWindow(comparison.firstLocalOnlyContext, minRow, minCol, maxRow, maxCol);
    includeCandidateContextInWindow(comparison.firstRebuildOnlyContext, minRow, minCol, maxRow, maxCol);
    if (minRow < 0 || minCol < 0 || maxRow < minRow || maxCol < minCol) {
        return {};
    }

    minRow = std::max(0, minRow - padding);
    minCol = std::max(0, minCol - padding);
    maxRow = std::min(comparison.localReconstructionImage->getNumRows() - 1, maxRow + padding);
    maxCol = std::min(comparison.localReconstructionImage->getNumCols() - 1, maxCol + padding);

    NonIncreasingPhaseOracleMismatchFixture fixture;
    fixture.found = true;
    fixture.threshold = comparison.threshold;
    fixture.completedPhase = comparison.completedPhase;
    fixture.beforeMaxToMin = comparison.beforeMaxToMin;
    fixture.minRow = minRow;
    fixture.minCol = minCol;
    fixture.maxRow = maxRow;
    fixture.maxCol = maxCol;
    fixture.croppedImage = cropImageWindow(
        comparison.localReconstructionImage,
        minRow,
        minCol,
        maxRow,
        maxCol);
    fixture.localOnlyContext = comparison.firstLocalOnlyContext;
    fixture.rebuildOnlyContext = comparison.firstRebuildOnlyContext;
    const int originalNumCols = comparison.localReconstructionImage->getNumCols();
    const int cropNumCols = maxCol - minCol + 1;
    fixture.localOnlySupportInCrop = remapSupportToCrop(
        std::span<const PixelId>(comparison.firstLocalOnlyContext.supportPixels),
        originalNumCols,
        minRow,
        minCol,
        cropNumCols);
    fixture.rebuildOnlySupportInCrop = remapSupportToCrop(
        std::span<const PixelId>(comparison.firstRebuildOnlyContext.supportPixels),
        originalNumCols,
        minRow,
        minCol,
        cropNumCols);
    return fixture;
}

inline NonIncreasingPhaseOracleMismatchFixture buildFirstPhaseOracleMismatchFixture(
    std::span<const NonIncreasingPhaseOracleComparison> comparisons,
    int padding = 1) {
    for (const auto &comparison : comparisons) {
        auto fixture = buildPhaseOracleMismatchFixture(comparison, padding);
        if (fixture.found) {
            return fixture;
        }
    }
    return {};
}

inline NonIncreasingPhaseOracleMismatchFixture buildSmallestPhaseOracleMismatchFixture(
    std::span<const NonIncreasingPhaseOracleComparison> comparisons,
    int padding = 1) {
    NonIncreasingPhaseOracleMismatchFixture best;
    int bestArea = -1;
    for (const auto &comparison : comparisons) {
        auto fixture = buildPhaseOracleMismatchFixture(comparison, padding);
        if (!fixture.found || !fixture.croppedImage) {
            continue;
        }
        const int area = fixture.croppedImage->getSize();
        if (bestArea < 0 || area < bestArea) {
            bestArea = area;
            best = std::move(fixture);
        }
    }
    return best;
}

struct NonIncreasingDirectComponentScratch {
    GenerationStampSet inRegionMarks;
    GenerationStampSet visitedMarks;
    GenerationStampSet releaseFocusMarks;
    GenerationStampSet releaseCandidateMarks;
    std::vector<std::vector<PixelId>> components;
    std::size_t componentCount = 0;
    std::vector<PixelId> seeds;
    std::vector<PixelId> stack;
    std::vector<NodeId> dualReleaseFocusNodes;
    std::vector<NodeId> primalReleaseFocusNodes;
    std::vector<NodeId> releaseWorklist;
};

struct NonIncreasingAttributeDirtyRecorder {
    std::vector<NodeId> *primalDirtySeeds = nullptr;
    std::vector<NodeId> *primalReleasedNodes = nullptr;
    std::vector<NodeId> *dualDirtySeeds = nullptr;
    std::vector<NodeId> *dualReleasedNodes = nullptr;
};

inline void appendNonIncreasingDirtyNode(std::vector<NodeId> *nodes,
                                         NodeId nodeId) {
    if (nodes != nullptr && nodeId != InvalidNode) {
        nodes->push_back(nodeId);
    }
}

inline void appendNonIncreasingDirtyNodeAndParent(const DynamicComponentTree &tree,
                                                  std::vector<NodeId> *nodes,
                                                  NodeId nodeId) {
    appendNonIncreasingDirtyNode(nodes, nodeId);
    if (nodes == nullptr ||
        nodeId == InvalidNode ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId) ||
        tree.isRoot(nodeId)) {
        return;
    }
    const NodeId parentId = tree.getNodeParent(nodeId);
    if (parentId != InvalidNode && parentId != nodeId) {
        nodes->push_back(parentId);
    }
}

inline void sortUniqueNonIncreasingNodes(std::vector<NodeId> &nodes) {
    nodes.erase(
        std::remove(nodes.begin(), nodes.end(), InvalidNode),
        nodes.end());
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(
        std::unique(nodes.begin(), nodes.end()),
        nodes.end());
}

inline std::span<const std::vector<PixelId>> directProperPartComponentsInto(
    const DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    NodeId nodeId,
    NonIncreasingDirectComponentScratch &state) {
    auto appendReusableComponent = [](NonIncreasingDirectComponentScratch &state) -> std::vector<PixelId> & {
        if (state.componentCount >= state.components.size()) {
            state.components.emplace_back();
        }
        auto &component = state.components[state.componentCount++];
        component.clear();
        return component;
    };

    auto componentSpan = [](const NonIncreasingDirectComponentScratch &state) {
        return std::span<const std::vector<PixelId>>(state.components.data(), state.componentCount);
    };

    state.componentCount = 0;

    if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
        return componentSpan(state);
    }

    const int numPixels = tree.getNumTotalProperParts();
    const int directProperPartCount = tree.getNumProperParts(nodeId);
    if (directProperPartCount <= 0) {
        return componentSpan(state);
    }
    if (directProperPartCount == 1) {
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            if (pixelId < 0 || pixelId >= numPixels) {
                continue;
            }
            auto &component = appendReusableComponent(state);
            component.push_back(pixelId);
            break;
        }
        return componentSpan(state);
    }

    state.inRegionMarks.ensureSize(static_cast<std::size_t>(numPixels));
    state.visitedMarks.ensureSize(static_cast<std::size_t>(numPixels));
    state.inRegionMarks.resetAll();
    state.visitedMarks.resetAll();
    state.seeds.clear();
    state.seeds.reserve(static_cast<std::size_t>(directProperPartCount));

    for (PixelId pixelId : tree.getProperParts(nodeId)) {
        if (pixelId < 0 || pixelId >= numPixels) {
            continue;
        }
        const auto pixelIndex = static_cast<std::size_t>(pixelId);
        if (state.inRegionMarks.isMarked(pixelIndex)) {
            continue;
        }
        state.inRegionMarks.mark(pixelIndex);
        state.seeds.push_back(pixelId);
    }

    state.stack.clear();
    for (PixelId seed : state.seeds) {
        const auto seedIdx = static_cast<std::size_t>(seed);
        if (state.visitedMarks.isMarked(seedIdx)) {
            continue;
        }

        auto &component = appendReusableComponent(state);
        state.stack.clear();
        state.stack.push_back(seed);
        state.visitedMarks.mark(seedIdx);

        while (!state.stack.empty()) {
            const PixelId pixelId = state.stack.back();
            state.stack.pop_back();
            component.push_back(pixelId);

            for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
                if (neighbor < 0 || neighbor >= numPixels) {
                    continue;
                }
                const auto neighborIdx = static_cast<std::size_t>(neighbor);
                if (!state.inRegionMarks.isMarked(neighborIdx) ||
                    state.visitedMarks.isMarked(neighborIdx)) {
                    continue;
                }
                state.visitedMarks.mark(neighborIdx);
                state.stack.push_back(neighbor);
            }
        }
    }

    return componentSpan(state);
}

inline std::vector<std::vector<PixelId>> directProperPartComponents(
    const DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    NodeId nodeId,
    NonIncreasingDirectComponentScratch *scratch) {
    NonIncreasingDirectComponentScratch localScratch;
    auto &state = scratch == nullptr ? localScratch : *scratch;
    const auto components = directProperPartComponentsInto(
        tree,
        adjacency,
        nodeId,
        state);
    return {components.begin(), components.end()};
}

inline std::vector<std::vector<PixelId>> directProperPartComponents(const DynamicComponentTree &tree,
                                                                    AdjacencyRelation &adjacency,
                                                                    NodeId nodeId) {
    return directProperPartComponents(tree, adjacency, nodeId, nullptr);
}

inline NodeId findNodeCaForExplicitRegion(const DynamicComponentTree &dualTree,
                                          std::span<const PixelId> pixelsC) {
    const bool isMaxtree = dualTree.isMaxtree();
    NodeId nodeCa = InvalidNode;
    int altitudeCa = 0;

    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
            continue;
        }
        const NodeId ownerNodeId = dualTree.getSmallestComponent(pixelId);
        if (ownerNodeId == InvalidNode || !dualTree.isNode(ownerNodeId) || !dualTree.isAlive(ownerNodeId)) {
            continue;
        }

        const int altitude = dualTree.getAltitude(ownerNodeId);
        if (nodeCa == InvalidNode ||
            (isMaxtree && altitude < altitudeCa) ||
            (!isMaxtree && altitude > altitudeCa)) {
            nodeCa = ownerNodeId;
            altitudeCa = altitude;
        }
    }

    return nodeCa;
}

inline void addUniqueReleaseFocusNode(std::vector<NodeId> &focusNodes,
                                      NodeId nodeId,
                                      GenerationStampSet *marks = nullptr) {
    if (nodeId == InvalidNode) {
        return;
    }
    if (marks != nullptr) {
        const auto index = static_cast<std::size_t>(nodeId);
        if (marks->isMarked(index)) {
            return;
        }
        marks->mark(index);
        focusNodes.push_back(nodeId);
        return;
    }
    if (std::find(focusNodes.begin(), focusNodes.end(), nodeId) != focusNodes.end()) {
        return;
    }
    focusNodes.push_back(nodeId);
}

inline NodeId findNodeCaAndAppendExplicitRegionOwnerFocusNodes(
    const DynamicComponentTree &dualTree,
    std::span<const PixelId> pixelsC,
    std::vector<NodeId> &focusNodes,
    GenerationStampSet *focusMarks = nullptr) {
    const bool isMaxtree = dualTree.isMaxtree();
    NodeId nodeCa = InvalidNode;
    int altitudeCa = 0;
    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
            continue;
        }
        const NodeId ownerNodeId = dualTree.getSmallestComponent(pixelId);
        if (ownerNodeId == InvalidNode ||
            !dualTree.isNode(ownerNodeId) ||
            !dualTree.isAlive(ownerNodeId)) {
            continue;
        }
        addUniqueReleaseFocusNode(focusNodes, ownerNodeId, focusMarks);

        const int altitude = dualTree.getAltitude(ownerNodeId);
        if (nodeCa == InvalidNode ||
            (isMaxtree && altitude < altitudeCa) ||
            (!isMaxtree && altitude > altitudeCa)) {
            nodeCa = ownerNodeId;
            altitudeCa = altitude;
        }
    }
    return nodeCa;
}

inline void appendExplicitRegionOwnerFocusNodes(const DynamicComponentTree &dualTree,
                                                std::span<const PixelId> pixelsC,
                                                std::vector<NodeId> &focusNodes) {
    (void) findNodeCaAndAppendExplicitRegionOwnerFocusNodes(
        dualTree,
        pixelsC,
        focusNodes,
        nullptr);
}

struct NonIncreasingMergeNodeCertificationResult {
    bool checked = false;
    bool passed = false;
    bool dualMatchesPrimalBeforeMerge = false;
    int componentCount = 0;
    int pixelCount = 0;
    int rejectedComponents = 0;
    int partialOwnerComponents = 0;
    int rawMissingTargetComponents = 0;
    int missingTargetComponents = 0;
    int partialOwnerAndMissingTargetComponents = 0;
    long long certificationMicros = 0;
    std::string message;
};

struct NonIncreasingPartialFResult {
    bool checked = false;
    bool accepted = false;
    bool rejected = false;
    bool mutationAttempted = false;
    bool localUpdateApplied = false;
    bool primalMergeApplied = false;
    bool mcpcLiveUpdateRequired = false;
    bool mcpcHiddenSupportRedundantLiveUpdate = false;
    bool mcpcParentLinkBridgeLiveUpdate = false;
    bool mcpcStableChildBridgeLiveUpdate = false;
    bool mcpcChildLinkBridgeLiveUpdate = false;
    bool mcpcEmptyOwnerReleaseLiveUpdate = false;
    bool mcpcStableChildMoveLiveUpdate = false;
    AltitudeType b{};
    int componentCount = 0;
    int pixelCount = 0;
    int rejectedComponents = 0;
    int partialOwnerComponents = 0;
    int rawMissingTargetComponents = 0;
    int missingTargetComponents = 0;
    int partialOwnerAndMissingTargetComponents = 0;
    long long preflightMicros = 0;
    long long collectComponentsUs = 0;
    long long localUpdateUs = 0;
    long long primalMergeUs = 0;
    std::string rejectReasonCode;
    std::string message;
};

inline void rejectNonIncreasingPartialF(NonIncreasingPartialFResult &out,
                                        const std::string &reason,
                                        const std::string &message = {}) {
    out.accepted = false;
    out.rejected = true;
    out.rejectReasonCode = reason;
    out.message = message.empty() ? reason : message;
}

inline NonIncreasingPartialFResult preflightNonIncreasingCasfPartialF(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    NodeId nodeId) {
    NonIncreasingPartialFResult out;
    const auto preflightStart = std::chrono::steady_clock::now();
    out.checked = true;

    auto finish = [&]() {
        out.preflightMicros =
            nonIncreasingElapsedMicros(preflightStart, std::chrono::steady_clock::now());
        return out;
    };

    if (!dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage())) {
        rejectNonIncreasingPartialF(
            out,
            "unsynchronized_primal_dual",
            "partial F requires synchronized primal/dual reconstructions");
        return finish();
    }

    if (nodeId == InvalidNode ||
        nodeId == primalTree.getRoot() ||
        !primalTree.isNode(nodeId) ||
        !primalTree.isAlive(nodeId)) {
        rejectNonIncreasingPartialF(
            out,
            "invalid_event",
            "partial F received an invalid primal node");
        return finish();
    }

    const NodeId parentId = primalTree.getNodeParent(nodeId);
    if (parentId == InvalidNode || parentId == nodeId || !primalTree.isAlive(parentId)) {
        rejectNonIncreasingPartialF(
            out,
            "invalid_event",
            "partial F cannot locate a live primal parent");
        return finish();
    }

    if (primalTree.getNumProperParts(nodeId) <= 0) {
        rejectNonIncreasingPartialF(
            out,
            "empty_C",
            "partial F requires a non-empty direct proper-part region");
        return finish();
    }

    out.b = static_cast<AltitudeType>(primalTree.getAltitude(parentId));
    const auto components = directProperPartComponents(primalTree, adjacency, nodeId);
    out.componentCount = static_cast<int>(components.size());
    if (components.empty()) {
        rejectNonIncreasingPartialF(
            out,
            "empty_C",
            "partial F found no connected component in C");
        return finish();
    }

    for (const auto &component : components) {
        if (component.empty()) {
            continue;
        }
        out.pixelCount += static_cast<int>(component.size());
        const auto analysis = analyzeNonIncreasingProblem1ExplicitRegionPreflight<AltitudeType>(
            dualTree,
            adjacency,
            std::span<const PixelId>(component),
            out.b,
            true);
        if (!analysis.acceptsPublishedProblem()) {
            ++out.rejectedComponents;
            rejectNonIncreasingPartialF(
                out,
                classifyProblem1ExplicitRegionFailureMode(analysis),
                problem1ExplicitRegionPreflightMessage(analysis.preflight.status));
            return finish();
        }
        if (!analysis.hasTargetLevelUnionCandidate) {
            ++out.rawMissingTargetComponents;
            ++out.missingTargetComponents;
        }
        if (analysis.requiresLocalMaterialization()) {
            ++out.partialOwnerComponents;
            if (!analysis.hasTargetLevelUnionCandidate) {
                ++out.partialOwnerAndMissingTargetComponents;
            }
        }
        if (!analysis.readyForLocalUpdater() &&
            !analysis.requiresLocalMaterialization() &&
            !analysis.requiresTargetNodeCreation()) {
            ++out.rejectedComponents;
            rejectNonIncreasingPartialF(
                out,
                "current_updater_limitation",
                problem1ExplicitRegionPreflightMessage(analysis.preflight.status).empty()
                    ? "partial F rejects updater limitations not covered by target-b creation"
                    : problem1ExplicitRegionPreflightMessage(analysis.preflight.status));
            return finish();
        }
    }

    if (out.pixelCount <= 0) {
        rejectNonIncreasingPartialF(
            out,
            "empty_C",
            "partial F found no valid pixel in C");
        return finish();
    }

    if (out.partialOwnerComponents > 0) {
        const AdjacencyRelationPtr adjacencyPtr(
            &adjacency,
            [](AdjacencyRelation *) {});
        const auto artifacts = buildMcpcDirectPreflightArtifacts(
            primalTree,
            dualTree,
            adjacencyPtr,
            nodeId);
        const auto d12 = diagnoseMcpcD12Domain(artifacts.preflight);
        if (d12.rejectReasonCode != "d12_exception_free_unproven") {
            rejectNonIncreasingPartialF(
                out,
                d12.rejectReasonCode.empty()
                    ? "d12_exception_free_unproven"
                    : d12.rejectReasonCode,
                d12.rejectionSummary.empty()
                    ? "partial F rejects P1-MCPC/partial-owner events until D12 is proved"
                    : d12.rejectionSummary);
            return finish();
        }

        const auto dryRun = diagnoseMcpcLiveMutationDryRun(dualTree, artifacts);
        if (!dryRun.accepted) {
            rejectNonIncreasingPartialF(
                out,
                dryRun.rejectReasonCode.empty()
                    ? "d12_exception_free_unproven"
                    : dryRun.rejectReasonCode,
                dryRun.rejectionSummary.empty()
                    ? "D12 live mutation dry-run has not accepted this P1-MCPC event"
                    : dryRun.rejectionSummary);
            return finish();
        }

        out.mcpcLiveUpdateRequired = true;
        out.mcpcHiddenSupportRedundantLiveUpdate =
            dryRun.candidateNonOracleHiddenSupportConnectorCount > 0;
        out.mcpcParentLinkBridgeLiveUpdate =
            dryRun.candidateNonOracleParentLinkBridgeCount > 0;
        out.mcpcStableChildBridgeLiveUpdate =
            dryRun.predictedStableChildBridgeNodeCount > 0;
        out.mcpcChildLinkBridgeLiveUpdate =
            dryRun.candidateNonOracleChildLinkBridgeCount > 0;
        out.mcpcEmptyOwnerReleaseLiveUpdate =
            dryRun.sourceOwnerToBeEmptiedCount > 0 ||
            dryRun.previousParentPossiblyEmptiedCount > 0;
        out.mcpcStableChildMoveLiveUpdate =
            dryRun.stableChildToMoveCount > 0;
    }

    out.accepted = true;
    out.rejected = false;
    out.message = "partial F preflight accepted";
    return finish();
}

inline NonIncreasingPartialFResult applyNonIncreasingCasfPartialF(
    DynamicComponentTree &primalTree,
    DynamicComponentTree &dualTree,
    DualMinMaxTreeIncrementalFilter<AltitudeType> &adjust,
    AdjacencyRelation &adjacency,
    NodeId nodeId) {
    auto out = preflightNonIncreasingCasfPartialF(
        primalTree,
        dualTree,
        adjacency,
        nodeId);
    if (!out.accepted) {
        return out;
    }

    out.mutationAttempted = true;
    if (out.mcpcLiveUpdateRequired) {
        const AdjacencyRelationPtr adjacencyPtr(
            &adjacency,
            [](AdjacencyRelation *) {});
        const auto artifacts = buildMcpcDirectPreflightArtifacts(
            primalTree,
            dualTree,
            adjacencyPtr,
            nodeId);
        const auto dryRun = diagnoseMcpcLiveMutationDryRun(dualTree, artifacts);
        if (!dryRun.accepted) {
            out.accepted = false;
            out.rejected = false;
            out.rejectReasonCode = "bug_mcpc_dry_run_rejected_after_preflight_acceptance";
            out.message = dryRun.rejectionSummary.empty()
                ? "P1-MCPC dry-run rejected after accepted partial F preflight"
                : dryRun.rejectionSummary;
            return out;
        }

        NonIncreasingStrippedLocalUpdateOptions options;
        options.checkNoEmptyLiveNodesAfterCanonicalization = true;
        options.checkFinalDualPrimalReconstruction = true;
        options.releaseEmptyNodesAfterPrimalMerge = true;
        options.releaseEmptySourceOwnersWithResidualChildren = true;
        const auto mcpcStart = std::chrono::steady_clock::now();
        const auto update = applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight(
            primalTree,
            dualTree,
            adjacencyPtr,
            nodeId,
            artifacts,
            options);
        out.localUpdateUs =
            nonIncreasingElapsedMicros(mcpcStart, std::chrono::steady_clock::now());
        out.localUpdateApplied = update.localRecompositionApplied;
        out.primalMergeApplied = update.applied;
        out.componentCount = update.recompositionNodeCount;
        out.pixelCount = update.localPixelCount;

        if (!update.applied ||
            !update.localRecompositionApplied ||
            !update.localRecompositionMatchesTarget ||
            !update.dualMatchesPrimalAfterMerge ||
            !update.emptyCreatedNodesAfterCanonicalization.empty()) {
            out.accepted = false;
            out.rejected = false;
            out.rejectReasonCode = "bug_mcpc_live_update_failed_after_preflight_acceptance";
            out.message = "P1-MCPC live update failed proof assertions after accepted preflight";
            return out;
        }

        out.message = "partial F P1-MCPC live update accepted";
        return out;
    }

    const auto directUpdate = adjust.mergeNodeIntoParentAndUpdate(
        &primalTree,
        &dualTree,
        nodeId);
    out.collectComponentsUs = directUpdate.collectComponentsUs;
    out.localUpdateUs = directUpdate.localUpdateUs;
    out.primalMergeUs = directUpdate.primalMergeUs;
    out.localUpdateApplied = directUpdate.localUpdateApplied;
    out.primalMergeApplied = directUpdate.primalMergeApplied;
    out.componentCount = directUpdate.componentCount;
    out.pixelCount = directUpdate.regionSize;

    if (!directUpdate.inputValid ||
        !directUpdate.localUpdateApplied ||
        !directUpdate.primalMergeApplied) {
        out.accepted = false;
        out.rejected = false;
        out.rejectReasonCode = "bug_unexpected_direct_update_rejection";
        out.message = directUpdate.message.empty()
            ? "partial F direct update rejected after accepted preflight"
            : directUpdate.message;
        return out;
    }

    out.message = "partial F applied";
    return out;
}

inline void recordNonIncreasingPartialFResult(
    NonIncreasingMergeEvent &event,
    const NonIncreasingPartialFResult &result,
    const DynamicComponentTree *primalTree = nullptr,
    const DynamicComponentTree *dualTree = nullptr,
    bool acceptedResultImpliesDualMatch = false) {
    event.partialFAccepted = result.accepted;
    event.partialFRejected = result.rejected;
    event.partialFMutationAttempted = result.mutationAttempted;
    event.partialFMcpcLiveUpdate =
        result.accepted && result.mcpcLiveUpdateRequired;
    event.partialFMcpcHiddenSupportRedundantLiveUpdate =
        result.accepted && result.mcpcHiddenSupportRedundantLiveUpdate;
    event.partialFMcpcParentLinkBridgeLiveUpdate =
        result.accepted && result.mcpcParentLinkBridgeLiveUpdate;
    event.partialFMcpcStableChildBridgeLiveUpdate =
        result.accepted && result.mcpcStableChildBridgeLiveUpdate;
    event.partialFMcpcChildLinkBridgeLiveUpdate =
        result.accepted && result.mcpcChildLinkBridgeLiveUpdate;
    event.partialFMcpcEmptyOwnerReleaseLiveUpdate =
        result.accepted && result.mcpcEmptyOwnerReleaseLiveUpdate;
    event.partialFMcpcStableChildMoveLiveUpdate =
        result.accepted && result.mcpcStableChildMoveLiveUpdate;
    event.partialFRejectReasonCode = result.rejectReasonCode;
    event.mergeNodeLocalApplied =
        result.localUpdateApplied && result.primalMergeApplied;
    event.mergeNodeLocalAccepted = result.accepted && event.mergeNodeLocalApplied;
    event.mergeNodeLocalComponentCount = result.componentCount;
    event.mergeNodeLocalPixelCount = result.pixelCount;
    event.mergeNodeLocalRejectedComponents =
        result.rejected ? std::max(1, result.rejectedComponents) : 0;
    event.mergeNodeLocalPartialOwnerComponents =
        result.partialOwnerComponents;
    event.mergeNodeLocalRawMissingTargetComponents =
        result.rawMissingTargetComponents;
    event.mergeNodeLocalMissingTargetComponents =
        result.missingTargetComponents;
    event.mergeNodeLocalPartialOwnerAndMissingTargetComponents =
        result.partialOwnerAndMissingTargetComponents;
    event.mergeNodeLocalTargetCreationRequired =
        result.missingTargetComponents > 0;
    event.mergeNodeLocalUpdateMicros =
        result.preflightMicros + result.collectComponentsUs + result.localUpdateUs;
    event.mergeNodePrimalMergeMicros = result.primalMergeUs;
    if (acceptedResultImpliesDualMatch) {
        event.mergeNodeLocalDualMatchesPrimalAfterMerge =
            event.mergeNodeLocalApplied;
    } else {
        event.mergeNodeLocalDualMatchesPrimalAfterMerge =
            event.mergeNodeLocalApplied &&
            primalTree != nullptr &&
            dualTree != nullptr &&
            dualTree->reconstructionImage()->isEqual(primalTree->reconstructionImage());
    }
    if (result.rejected) {
        event.mergeNodeLocalProblem1FailureMode = result.rejectReasonCode;
        event.mergeNodeLocalFailureSummary = result.message;
    }
    if (event.mergeNodeLocalAccepted) {
        event.applied = true;
    }
}

inline NonIncreasingMergeNodeCertificationResult certifyMergeNodeExplicitComponents(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    AdjacencyRelation &adjacency,
    std::span<const std::vector<PixelId>> components,
    AltitudeType b,
    bool allowTargetNodeCreation = false) {
    NonIncreasingMergeNodeCertificationResult out;
    const auto certificationStart = std::chrono::steady_clock::now();
    out.checked = true;
    out.dualMatchesPrimalBeforeMerge =
        dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
    if (!out.dualMatchesPrimalBeforeMerge) {
        out.message = "merge-node certification requires synchronized primal/dual reconstructions";
        out.certificationMicros =
            nonIncreasingElapsedMicros(certificationStart, std::chrono::steady_clock::now());
        return out;
    }

    for (const auto &component : components) {
        if (component.empty()) {
            continue;
        }
        ++out.componentCount;
        out.pixelCount += static_cast<int>(component.size());
        const auto analysis = analyzeNonIncreasingProblem1ExplicitRegionPreflight<AltitudeType>(
            dualTree,
            adjacency,
            std::span<const PixelId>(component),
            b,
            true);
        if (!analysis.acceptsPublishedProblem()) {
            ++out.rejectedComponents;
            if (out.message.empty()) {
                out.message =
                    problem1ExplicitRegionPreflightMessage(analysis.preflight.status);
            }
            continue;
        }
        const bool missingTarget = !analysis.hasTargetLevelUnionCandidate;
        if (missingTarget) {
            ++out.rawMissingTargetComponents;
        }
        if (analysis.requiresLocalMaterialization()) {
            ++out.partialOwnerComponents;
            if (missingTarget) {
                ++out.partialOwnerAndMissingTargetComponents;
            }
            if (out.message.empty()) {
                out.message = "merge-node certification rejected a component with partial owners";
            }
            continue;
        }
        if (missingTarget) {
            ++out.missingTargetComponents;
            if (!allowTargetNodeCreation && out.message.empty()) {
                out.message = "merge-node certification requires an existing level-b union candidate";
            }
            if (!allowTargetNodeCreation) {
                continue;
            }
        }
        if (!analysis.readyForLocalUpdater() &&
            !(allowTargetNodeCreation && analysis.requiresTargetNodeCreation())) {
            ++out.rejectedComponents;
            if (out.message.empty()) {
                out.message = problem1ExplicitRegionPreflightMessage(analysis.preflight.status).empty()
                    ? "merge-node certification rejected a component incompatible with the direct updater"
                    : problem1ExplicitRegionPreflightMessage(analysis.preflight.status);
            }
        }
    }

    out.passed =
        out.rejectedComponents == 0 &&
        out.partialOwnerComponents == 0 &&
        (allowTargetNodeCreation || out.missingTargetComponents == 0);
    if (out.passed && out.message.empty()) {
        out.message = "merge-node certification passed";
    } else if (!out.passed && out.message.empty()) {
        out.message = "merge-node certification rejected the event";
    }
    out.certificationMicros =
        nonIncreasingElapsedMicros(certificationStart, std::chrono::steady_clock::now());
    return out;
}

inline void mergeNodeIntoParentWithExplicitDualUpdate(DynamicComponentTree &primalTree,
                                                      DynamicComponentTree &dualTree,
                                                      DualMinMaxTreeIncrementalFilter<AltitudeType> &adjust,
                                                      const AdjacencyRelationPtr &adjacency,
                                                      NonIncreasingDualMaintenance dualMaintenance,
                                                      NodeId nodeId,
                                                      NonIncreasingMergeEvent *event = nullptr,
                                                      bool repairDualMismatchAfterMerge = false,
                                                      bool useOracleLocalLeafExplicitFastPath = true,
                                                      const NonIncreasingStrippedLocalUpdateOptions &strippedLocalUpdateOptions = {},
                                                      NonIncreasingStrippedCommitMode strippedCommitMode = NonIncreasingStrippedCommitMode::SafeCopy,
                                                      bool skipDualUpdateAndMergePrimalOnly = false,
                                                      bool mcpcDirectPreflightDiagnosticsEnabled = false,
                                                      bool mcpcDirectV1StrictPrecheckEnabled = false,
                                                      bool mcpcNoCopyDirectCommitEnabled = false,
                                                      bool mcpcValidatedDirectCommitEnabled = false,
                                                      bool allowMergeNodePhaseRebuildFallback = true,
                                                      NonIncreasingDirectComponentScratch *componentScratch = nullptr,
                                                      NonIncreasingAttributeDirtyRecorder *attributeDirtyRecorder = nullptr) {
    if (!adjacency) {
        throw std::runtime_error("mergeNodeIntoParentWithExplicitDualUpdate requires a valid adjacency relation");
    }
    if (nodeId == InvalidNode ||
        nodeId == primalTree.getRoot() ||
        !primalTree.isNode(nodeId) ||
        !primalTree.isAlive(nodeId)) {
        return;
    }

    const NodeId parentId = primalTree.getNodeParent(nodeId);
    if (parentId == InvalidNode || parentId == nodeId || !primalTree.isAlive(parentId)) {
        return;
    }

    const auto b = static_cast<AltitudeType>(primalTree.getAltitude(parentId));
    std::chrono::steady_clock::time_point componentCollectStart;
    if (event != nullptr) {
        componentCollectStart = std::chrono::steady_clock::now();
    }
    NonIncreasingDirectComponentScratch localComponentScratch;
    auto &directComponentScratch =
        componentScratch == nullptr ? localComponentScratch : *componentScratch;
    const auto components = directProperPartComponentsInto(
        primalTree,
        *adjacency,
        nodeId,
        directComponentScratch);
    long long componentCollectMicros = 0;
    if (event != nullptr) {
        componentCollectMicros =
            nonIncreasingElapsedMicros(componentCollectStart, std::chrono::steady_clock::now());
    }

    const bool collectDetailedEventDiagnostics =
        event != nullptr && event->collectDetailedDiagnostics;

    if (event != nullptr) {
        event->nodeId = nodeId;
        event->parentId = parentId;
        event->nodeLevel = primalTree.getAltitude(nodeId);
        event->parentLevel = primalTree.getAltitude(parentId);
        event->numChildren = primalTree.getNumChildren(nodeId);
        event->directProperParts = primalTree.getNumProperParts(nodeId);
        event->componentCount = static_cast<int>(components.size());
        event->mergeNodeComponentCollectMicros = componentCollectMicros;
        event->totalComponentPixels = 0;
        event->minRow = -1;
        event->minCol = -1;
        event->maxRow = -1;
        event->maxCol = -1;
        if (collectDetailedEventDiagnostics) {
            event->componentSizes.clear();
            event->componentNodeCa.clear();
            event->componentAltitudeCa.clear();
            event->affectedPixelsPreview.clear();
            event->componentSizes.reserve(components.size());
            event->componentNodeCa.reserve(components.size());
            event->componentAltitudeCa.reserve(components.size());
            constexpr std::size_t affectedPreviewLimit = 16;
            const int numCols = primalTree.getNumColsOfImage();
            for (const auto &component : components) {
                event->totalComponentPixels += static_cast<int>(component.size());
                event->componentSizes.push_back(static_cast<int>(component.size()));
                for (PixelId pixelId : component) {
                    if (pixelId < 0 || numCols <= 0) {
                        continue;
                    }
                    const int row = pixelId / numCols;
                    const int col = pixelId % numCols;
                    if (event->minRow < 0 || row < event->minRow) {
                        event->minRow = row;
                    }
                    if (event->minCol < 0 || col < event->minCol) {
                        event->minCol = col;
                    }
                    if (event->maxRow < 0 || row > event->maxRow) {
                        event->maxRow = row;
                    }
                    if (event->maxCol < 0 || col > event->maxCol) {
                        event->maxCol = col;
                    }
                    if (event->affectedPixelsPreview.size() < affectedPreviewLimit) {
                        event->affectedPixelsPreview.push_back(pixelId);
                    }
                }
            }
        }
    }

    if (skipDualUpdateAndMergePrimalOnly) {
        const auto mergeStart = std::chrono::steady_clock::now();
        primalTree.mergeNodeIntoParent(nodeId);
        const auto mergeEnd = std::chrono::steady_clock::now();
        const auto releaseStart = std::chrono::steady_clock::now();
        (void) releaseAllEmptyLeafNodes(primalTree);
        const auto releaseEnd = std::chrono::steady_clock::now();
        if (event != nullptr) {
            const auto mergeMicros = nonIncreasingElapsedMicros(mergeStart, mergeEnd);
            const auto releaseMicros = nonIncreasingElapsedMicros(releaseStart, releaseEnd);
            if (dualMaintenance == NonIncreasingDualMaintenance::MergeNodeUpdateThenPhaseRebuild ||
                dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild ||
                dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild ||
                dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild) {
                event->mergeNodeLocalSkippedByPhaseRebuild = true;
                event->mergeNodeLocalFallbackToPhaseRebuild = true;
                event->mergeNodePrimalMergeMicros = mergeMicros;
                event->mergeNodeEmptyReleaseMicros = releaseMicros;
            } else {
                event->strippedLocalSkippedByPhaseRebuild = true;
                event->strippedLocalFallbackToPhaseRebuild = true;
                event->strippedPrimalMergeReleaseMicros = mergeMicros + releaseMicros;
            }
            event->applied = true;
        }
        return;
    }

    const auto fallbackToPrimalMergeWithPhaseRebuild = [&](const std::string &message) {
        const auto fallbackStart = std::chrono::steady_clock::now();
        primalTree.mergeNodeIntoParent(nodeId);
        const auto fallbackMergeEnd = std::chrono::steady_clock::now();
        const auto fallbackReleaseStart = std::chrono::steady_clock::now();
        (void) releaseAllEmptyLeafNodes(primalTree);
        const auto fallbackEnd = std::chrono::steady_clock::now();
        if (event != nullptr) {
            event->mergeNodeLocalFallbackToPhaseRebuild = true;
            event->mergeNodePrimalMergeMicros =
                nonIncreasingElapsedMicros(fallbackStart, fallbackMergeEnd);
            event->mergeNodeEmptyReleaseMicros =
                nonIncreasingElapsedMicros(fallbackReleaseStart, fallbackEnd);
            if (event->mergeNodeLocalFailureSummary.empty() && !message.empty()) {
                event->mergeNodeLocalFailureSummary = message;
            }
            event->applied = true;
        }
    };

    if (dualMaintenance == NonIncreasingDualMaintenance::StrippedThenOracleLocalRecomposition ||
        dualMaintenance == NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition ||
        dualMaintenance == NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition) {
        const bool useRebuildFallback =
            dualMaintenance == NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition;
        const bool usePhaseRebuildFallback =
            dualMaintenance == NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition;
        const bool directCommit = strippedCommitMode == NonIncreasingStrippedCommitMode::Direct;
        NonIncreasingStrippedLocalUpdateOptions effectiveStrippedLocalUpdateOptions =
            strippedLocalUpdateOptions;
        if (useRebuildFallback || usePhaseRebuildFallback) {
            effectiveStrippedLocalUpdateOptions.releaseEmptyNodesAfterPrimalMerge = true;
        }
        if (event != nullptr) {
            event->strippedLocalAttempted = true;
            event->strippedDirectCommit = directCommit;
        }

        auto strippedUpdateAccepted = [&](const NonIncreasingStrippedLocalUpdateResult &strippedUpdate) {
            const bool finalCheckOk =
                !effectiveStrippedLocalUpdateOptions.checkFinalDualPrimalReconstruction ||
                strippedUpdate.dualMatchesPrimalAfterMerge;
            const bool targetCheckOk =
                strippedUpdate.skippedNoImageChange ||
                !effectiveStrippedLocalUpdateOptions.checkMaterializationTargetReconstruction ||
                strippedUpdate.localRecompositionMatchesTarget;
            const bool planAuditOk =
                !effectiveStrippedLocalUpdateOptions.auditPlanOracleDependencies ||
                strippedUpdate.localPlanOracleFree;
            const bool materializationAuditOk =
                !effectiveStrippedLocalUpdateOptions.auditMaterializationOracleDependencies ||
                strippedUpdate.materializationInputOracleFree;
            const bool noEmptyCreatedNodes =
                !(useRebuildFallback || usePhaseRebuildFallback) ||
                strippedUpdate.emptyCreatedNodesAfterCanonicalization.empty();
            const bool canonicalByConstruction =
                strippedUpdate.postCanonicalizationChangesRequired == 0;
            return
                strippedUpdate.applied &&
                finalCheckOk &&
                (strippedUpdate.skippedNoImageChange ||
                 (targetCheckOk &&
                  planAuditOk &&
                  materializationAuditOk &&
                  noEmptyCreatedNodes &&
                  canonicalByConstruction));
        };

        auto recordStrippedUpdate = [&](const NonIncreasingStrippedLocalUpdateResult &strippedUpdate,
                                        bool strippedAccepted) {
            if (event != nullptr) {
                event->strippedLocalApplied = strippedUpdate.applied;
                event->strippedLocalPlanAuditChecked =
                    strippedUpdate.localPlanOracleAuditChecked;
                event->strippedLocalMatchesTarget = strippedUpdate.localRecompositionMatchesTarget;
                event->strippedLocalTargetChecked =
                    strippedUpdate.localRecompositionTargetChecked;
                event->strippedLocalStructureChecked =
                    strippedUpdate.localStructureAgainstRebuildChecked;
                event->strippedLocalStructureMatchesRebuild =
                    strippedUpdate.localStructureMatchesRebuild;
                event->strippedLocalPlanOracleFree = strippedUpdate.localPlanOracleFree;
                event->strippedMaterializationAuditChecked =
                    strippedUpdate.materializationInputOracleAuditChecked;
                event->strippedMaterializationInputOracleFree =
                    strippedUpdate.materializationInputOracleFree;
                event->strippedFinalReconstructionChecked =
                    strippedUpdate.finalReconstructionChecked;
                event->strippedLocalDualMatchesPrimalAfterMerge =
                    strippedUpdate.dualMatchesPrimalAfterMerge;
                event->strippedLocalRecompositionNodes = strippedUpdate.recompositionNodeCount;
                event->strippedLocalCreatedNodes = strippedUpdate.createdNodeCount;
                event->strippedLocalPixelCount = strippedUpdate.localPixelCount;
                event->strippedCreatedConnectorNodes =
                    strippedUpdate.createdConnectorNodeCount;
                event->strippedCreatedHiddenExternalParentConnectorNodes =
                    strippedUpdate.createdHiddenExternalParentConnectorNodeCount;
                event->strippedCreatedOracleBridgeNodes =
                    strippedUpdate.createdOracleBridgeNodeCount;
                event->strippedCreatedStableChildBridgeNodes =
                    strippedUpdate.createdStableChildBridgeNodeCount;
                event->strippedCreatedHiddenSupportConnectorNodes =
                    strippedUpdate.createdHiddenSupportConnectorNodeCount;
                event->strippedCreatedNonOracleHiddenSupportConnectorNodes =
                    strippedUpdate.createdNonOracleHiddenSupportConnectorNodeCount;
                event->strippedCreatedNonOracleParentLinkBridgeNodes =
                    strippedUpdate.createdNonOracleParentLinkBridgeNodeCount;
                event->strippedCreatedNonOracleChildLinkBridgeNodes =
                    strippedUpdate.createdNonOracleChildLinkBridgeNodeCount;
                event->strippedNonOracleHiddenSupportCandidates =
                    strippedUpdate.nonOracleHiddenSupportCandidateCount;
                event->strippedNonOracleHiddenSupportAlreadyRepresented =
                    strippedUpdate.nonOracleHiddenSupportAlreadyRepresentedCount;
                event->strippedNonOracleHiddenSupportMaterializable =
                    strippedUpdate.nonOracleHiddenSupportMaterializableCount;
                event->strippedNonOracleParentLinksTouchingPlan =
                    strippedUpdate.nonOracleParentLinkTouchesPlanCount;
                event->strippedNonOracleParentLinksBothRepresented =
                    strippedUpdate.nonOracleParentLinkBothRepresentedCount;
                event->strippedNonOracleParentLinksBlockedByModifiedSnapshot =
                    strippedUpdate.nonOracleParentLinkBlockedByModifiedSnapshotCount;
                event->strippedNonOracleParentLinksInvalidLevel =
                    strippedUpdate.nonOracleParentLinkInvalidLevelCount;
                event->strippedNonOracleParentLinksAlreadyDirect =
                    strippedUpdate.nonOracleParentLinkAlreadyDirectCount;
                event->strippedNonOracleParentLinksAlreadyIntermediate =
                    strippedUpdate.nonOracleParentLinkAlreadyIntermediateCount;
                event->strippedNonOracleParentLinksMoved =
                    strippedUpdate.nonOracleParentLinkMovedCount;
                event->strippedNonOracleParentLinkBridgeAttempts =
                    strippedUpdate.nonOracleParentLinkBridgeAttemptCount;
                event->strippedNonOracleParentLinkBridgeMaterializable =
                    strippedUpdate.nonOracleParentLinkBridgeMaterializableCount;
                event->strippedNonOracleChildLinkBridgeAttempts =
                    strippedUpdate.nonOracleChildLinkBridgeAttemptCount;
                event->strippedNonOracleChildLinkBridgeInvalidLevel =
                    strippedUpdate.nonOracleChildLinkBridgeInvalidLevelCount;
                event->strippedNonOracleChildLinkBridgeMaterializable =
                    strippedUpdate.nonOracleChildLinkBridgeMaterializableCount;
                event->strippedEmptyCreatedNodesAfterCanonicalization =
                    static_cast<int>(
                        strippedUpdate.emptyCreatedNodesAfterCanonicalization.size());
                event->strippedPostEmptyNodeReleasesRequired =
                    strippedUpdate.postEmptyNodeReleasesRequired;
                event->strippedPostSameLevelParentChildCollapsesRequired =
                    strippedUpdate.postSameLevelParentChildCollapsesRequired;
                event->strippedPostSameLevelSiblingMergesRequired =
                    strippedUpdate.postSameLevelSiblingMergesRequired;
                event->strippedPostComparableSiblingMovesRequired =
                    strippedUpdate.postComparableSiblingMovesRequired;
                event->strippedPostCanonicalizationChangesRequired =
                    strippedUpdate.postCanonicalizationChangesRequired;
                if (collectDetailedEventDiagnostics) {
                    event->strippedCreatedNodeLevels = strippedUpdate.createdNodeLevels;
                    event->strippedCreatedNodePlanIndices =
                        strippedUpdate.createdNodePlanIndices;
                    event->strippedEmptyCreatedNodeDiagnostics =
                        strippedUpdate.emptyCreatedNodesAfterCanonicalization;
                    event->strippedEmptyConnectorNodesAfterCanonicalization =
                        static_cast<int>(std::count_if(
                            event->strippedEmptyCreatedNodeDiagnostics.begin(),
                            event->strippedEmptyCreatedNodeDiagnostics.end(),
                            [](const auto &node) {
                                return node.planIndex < 0;
                            }));
                    event->strippedLocalStructureMismatchSummary =
                        strippedUpdate.localStructureMismatchSummary;
                }
                event->strippedLocalReleasedSourceOwners =
                    strippedUpdate.releasedSourceOwnerCount;
                event->strippedInputExtractionMicros = strippedUpdate.inputExtractionMicros;
                event->strippedPartitionPlanMicros = strippedUpdate.partitionPlanMicros;
                event->strippedRecompositionPlanMicros = strippedUpdate.recompositionPlanMicros;
                event->strippedPlanAuditMicros = strippedUpdate.planAuditMicros;
                event->strippedMaterializationContextMicros =
                    strippedUpdate.materializationContextMicros;
                event->strippedMaterializationAuditMicros =
                    strippedUpdate.materializationAuditMicros;
                event->strippedMaterializationMicros = strippedUpdate.materializationMicros;
                event->strippedPostLocalSnapshotMicros = strippedUpdate.postLocalSnapshotMicros;
                event->strippedPrimalMergeReleaseMicros = strippedUpdate.primalMergeReleaseMicros;
                event->strippedFinalReconstructionCheckMicros =
                    strippedUpdate.finalReconstructionCheckMicros;
                event->strippedMaterializationValidationMicros =
                    strippedUpdate.materializationValidationMicros;
                event->strippedMaterializationCreateNodesMicros =
                    strippedUpdate.materializationCreateNodesMicros;
                event->strippedMaterializationAttachNodesMicros =
                    strippedUpdate.materializationAttachNodesMicros;
                event->strippedMaterializationMoveProperPartsMicros =
                    strippedUpdate.materializationMoveProperPartsMicros;
                event->strippedMaterializationChildBridgeMicros =
                    strippedUpdate.materializationChildBridgeMicros;
                event->strippedMaterializationStructuralRepairMicros =
                    strippedUpdate.materializationStructuralRepairMicros;
                event->strippedMaterializationReleaseMicros =
                    strippedUpdate.materializationReleaseEmptySourceOwnersMicros;
                event->strippedMaterializationRootPromotionMicros =
                    strippedUpdate.materializationRootPromotionMicros;
                event->strippedMaterializationCanonicalizationMicros =
                    strippedUpdate.materializationCanonicalizationMicros;
                event->strippedMaterializationTargetCheckMicros =
                    strippedUpdate.materializationTargetReconstructionCheckMicros;
                event->strippedLocalAccepted = strippedAccepted;
            }
        };

        if (directCommit) {
            const auto preflightCopyStart = std::chrono::steady_clock::now();
            DynamicComponentTree primalTrial = primalTree;
            DynamicComponentTree dualTrial = dualTree;
            const auto preflightCopyEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->strippedLocalCopyMicros =
                    nonIncreasingElapsedMicros(preflightCopyStart, preflightCopyEnd);
            }

            bool preflightAccepted = false;
            const auto preflightStart = std::chrono::steady_clock::now();
            try {
                const auto preflightUpdate = applyStrippedNonOracleLocalMergeUpdate(
                    primalTrial,
                    dualTrial,
                    adjacency,
                    nodeId,
                    effectiveStrippedLocalUpdateOptions);
                preflightAccepted = strippedUpdateAccepted(preflightUpdate);
                recordStrippedUpdate(preflightUpdate, preflightAccepted);
            } catch (const std::exception &e) {
                if (event != nullptr) {
                    event->strippedLocalFailureSummary = e.what();
                }
            }
            const auto preflightEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->strippedLocalUpdateMicros =
                    nonIncreasingElapsedMicros(preflightStart, preflightEnd);
            }
            if (!preflightAccepted) {
                if (event != nullptr && event->strippedLocalFailureSummary.empty()) {
                    if (event->strippedEmptyCreatedNodesAfterCanonicalization != 0) {
                        event->strippedLocalFailureSummary =
                            "direct stripped local update preflight left created nodes without proper parts";
                    } else if (event->strippedPostCanonicalizationChangesRequired != 0) {
                        event->strippedLocalFailureSummary =
                            "direct stripped local update preflight required corrective post-canonicalization";
                    } else {
                        event->strippedLocalFailureSummary =
                            "direct stripped local update preflight did not satisfy the active acceptance checks";
                    }
                }
                if (!useRebuildFallback && !usePhaseRebuildFallback) {
                    throw std::runtime_error(
                        event != nullptr && !event->strippedLocalFailureSummary.empty()
                            ? event->strippedLocalFailureSummary
                            : "direct stripped local update preflight did not satisfy the active acceptance checks");
                }
                if (usePhaseRebuildFallback) {
                    const auto fallbackStart = std::chrono::steady_clock::now();
                    primalTree.mergeNodeIntoParent(nodeId);
                    (void) releaseAllEmptyLeafNodes(primalTree);
                    const auto fallbackEnd = std::chrono::steady_clock::now();
                    if (event != nullptr) {
                        event->strippedLocalFallbackToPhaseRebuild = true;
                        event->strippedPrimalMergeReleaseMicros =
                            nonIncreasingElapsedMicros(fallbackStart, fallbackEnd);
                        event->applied = true;
                    }
                    return;
                }
                const auto fallbackStart = std::chrono::steady_clock::now();
                primalTree.mergeNodeIntoParent(nodeId);
                (void) releaseAllEmptyLeafNodes(primalTree);
                dualTree.build(primalTree.reconstructionImage()->clone(), dualTree.isMaxtree(), adjacency);
                const auto fallbackEnd = std::chrono::steady_clock::now();
                if (event != nullptr) {
                    event->strippedLocalFallbackToRebuild = true;
                    event->strippedRebuildFallbackMicros =
                        nonIncreasingElapsedMicros(fallbackStart, fallbackEnd);
                    event->dualRebuiltAfterLocalMismatch = true;
                    event->strippedLocalDualMatchesPrimalAfterMerge =
                        dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
                    event->applied = true;
                }
                return;
            }

            bool strippedAccepted = false;
            const auto strippedStart = std::chrono::steady_clock::now();
            try {
                const auto strippedUpdate = applyStrippedNonOracleLocalMergeUpdate(
                    primalTree,
                    dualTree,
                    adjacency,
                    nodeId,
                    effectiveStrippedLocalUpdateOptions);
                strippedAccepted = strippedUpdateAccepted(strippedUpdate);
                recordStrippedUpdate(strippedUpdate, strippedAccepted);
            } catch (const std::exception &e) {
                if (event != nullptr) {
                    event->strippedLocalFailureSummary = e.what();
                }
                throw std::runtime_error(std::string("direct stripped local update failed after accepted preflight: ") + e.what());
            }
            const auto strippedEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->strippedLocalUpdateMicros +=
                    nonIncreasingElapsedMicros(strippedStart, strippedEnd);
            }
            if (!strippedAccepted) {
                if (event != nullptr && event->strippedLocalFailureSummary.empty()) {
                    if (event->strippedEmptyCreatedNodesAfterCanonicalization != 0) {
                        event->strippedLocalFailureSummary =
                            "direct stripped local update left created nodes without proper parts after accepted preflight";
                    } else if (event->strippedPostCanonicalizationChangesRequired != 0) {
                        event->strippedLocalFailureSummary =
                            "direct stripped local update required corrective post-canonicalization after accepted preflight";
                    } else {
                        event->strippedLocalFailureSummary =
                            "direct stripped local update did not satisfy the active acceptance checks after accepted preflight";
                    }
                }
                throw std::runtime_error(
                    event != nullptr && !event->strippedLocalFailureSummary.empty()
                        ? event->strippedLocalFailureSummary
                        : "direct stripped local update did not satisfy the active acceptance checks after accepted preflight");
            }
            if (event != nullptr) {
                event->applied = true;
            }
            return;
        }

        const auto copyStart = std::chrono::steady_clock::now();
        DynamicComponentTree primalTrial = primalTree;
        DynamicComponentTree dualTrial = dualTree;
        const auto copyEnd = std::chrono::steady_clock::now();
        if (event != nullptr) {
            event->strippedLocalCopyMicros =
                nonIncreasingElapsedMicros(copyStart, copyEnd);
        }

        bool strippedAccepted = false;
        const auto strippedStart = std::chrono::steady_clock::now();
        try {
            const auto strippedUpdate = applyStrippedNonOracleLocalMergeUpdate(
                primalTrial,
                dualTrial,
                adjacency,
                nodeId,
                effectiveStrippedLocalUpdateOptions);
            strippedAccepted = strippedUpdateAccepted(strippedUpdate);
            recordStrippedUpdate(strippedUpdate, strippedAccepted);
        } catch (const std::exception &e) {
            if (event != nullptr) {
                event->strippedLocalFailureSummary = e.what();
            }
        }
        const auto strippedEnd = std::chrono::steady_clock::now();
        if (event != nullptr) {
            event->strippedLocalUpdateMicros =
                nonIncreasingElapsedMicros(strippedStart, strippedEnd);
        }

        if (strippedAccepted) {
            primalTree = std::move(primalTrial);
            dualTree = std::move(dualTrial);
            if (event != nullptr) {
                event->applied = true;
            }
            return;
        }

        if (event != nullptr && event->strippedLocalFailureSummary.empty()) {
            if (useRebuildFallback &&
                event->strippedEmptyCreatedNodesAfterCanonicalization != 0) {
                event->strippedLocalFailureSummary =
                    "stripped local update left created nodes without proper parts";
            } else if (event->strippedPostCanonicalizationChangesRequired != 0) {
                event->strippedLocalFailureSummary =
                    "stripped local update required corrective post-canonicalization";
            } else {
                event->strippedLocalFailureSummary =
                    "stripped local update did not satisfy the active acceptance checks";
            }
        }

        if (useRebuildFallback || usePhaseRebuildFallback) {
            if (event != nullptr) {
                if (usePhaseRebuildFallback) {
                    event->strippedLocalFallbackToPhaseRebuild = true;
                } else {
                    event->strippedLocalFallbackToRebuild = true;
                }
            }

            const auto fallbackStart = std::chrono::steady_clock::now();
            primalTree.mergeNodeIntoParent(nodeId);
            (void) releaseAllEmptyLeafNodes(primalTree);
            if (usePhaseRebuildFallback) {
                const auto fallbackEnd = std::chrono::steady_clock::now();
                if (event != nullptr) {
                    event->strippedPrimalMergeReleaseMicros =
                        nonIncreasingElapsedMicros(fallbackStart, fallbackEnd);
                    event->applied = true;
                }
                return;
            }
            dualTree.build(primalTree.reconstructionImage()->clone(), dualTree.isMaxtree(), adjacency);
            const auto fallbackEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->strippedRebuildFallbackMicros =
                    nonIncreasingElapsedMicros(fallbackStart, fallbackEnd);
                event->dualRebuiltAfterLocalMismatch = true;
                event->strippedLocalDualMatchesPrimalAfterMerge =
                    dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
                event->applied = true;
            }
            return;
        }

        if (event != nullptr) {
            event->strippedLocalFallbackToOracle = true;
        }

        const auto fallbackStart = std::chrono::steady_clock::now();
        const auto localUpdate = applyOracleGuidedLocalMergeUpdate(
            primalTree,
            dualTree,
            adjacency,
            nodeId);
        const auto fallbackEnd = std::chrono::steady_clock::now();
        if (event != nullptr) {
            event->strippedOracleFallbackMicros =
                nonIncreasingElapsedMicros(fallbackStart, fallbackEnd);
            event->oracleGuidedLocalApplied = localUpdate.localRecompositionApplied;
            event->oracleGuidedLocalMatchesTarget = localUpdate.localRecompositionMatchesTarget;
            event->oracleGuidedLocalStructureMatchesOracle = localUpdate.localStructureMatchesOracle;
            event->oracleGuidedDualMatchesPrimalAfterMerge = localUpdate.dualMatchesPrimalAfterMerge;
            event->oracleGuidedRecompositionNodes = localUpdate.recompositionNodeCount;
            event->oracleGuidedCreatedNodes = localUpdate.createdNodeCount;
            event->oracleGuidedReleasedSourceOwners = localUpdate.releasedSourceOwnerCount;
            event->oracleGuidedStructureMismatchSummary =
                localUpdate.localStructureMismatchSummary;
            event->applied = localUpdate.applied;
        }
        return;
    }

    if (dualMaintenance == NonIncreasingDualMaintenance::PartialF) {
        if (event != nullptr) {
            event->mergeNodeLocalAttempted = true;
            event->partialFAttempted = true;
            event->mergeNodeLocalDualMatchesPrimalBeforeMerge =
                dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
        }

        const auto result = applyNonIncreasingCasfPartialF(
            primalTree,
            dualTree,
            adjust,
            *adjacency,
            nodeId);

        if (event != nullptr) {
            recordNonIncreasingPartialFResult(
                *event,
                result,
                &primalTree,
                &dualTree);
        }
        return;
    }

    if (dualMaintenance == NonIncreasingDualMaintenance::Problem1CloseRecompose) {
        if (collectDetailedEventDiagnostics) {
            for (const auto &component : components) {
                if (component.empty()) {
                    continue;
                }
                const NodeId nodeCa = findNodeCaForExplicitRegion(
                    dualTree,
                    std::span<const PixelId>(component));
                event->componentNodeCa.push_back(nodeCa);
                event->componentAltitudeCa.push_back(
                    nodeCa == InvalidNode ? -1 : dualTree.getAltitude(nodeCa));
            }
        }

        bool problem1RestrictedCandidateAccepted = true;
        std::string problem1RestrictedRejection;
        for (const auto &component : components) {
            if (component.empty()) {
                continue;
            }
            const auto analysis = analyzeNonIncreasingProblem1ExplicitRegionPreflight<AltitudeType>(
                dualTree,
                *adjacency,
                std::span<const PixelId>(component),
                b,
                true);
            if (analysis.acceptsPublishedProblem() && analysis.affectedOwners() > 1) {
                if (event != nullptr) {
                    ++event->mergeNodeLocalProblem1PreflightRejectedComponents;
                    ++event->mergeNodeLocalProblem1PreflightRejectedByPublishedProblem;
                    if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                        event->mergeNodeLocalProblem1FailureMode = "multiple_nodeCa";
                    }
                }
                problem1RestrictedCandidateAccepted = false;
                problem1RestrictedRejection =
                    "problem1-restricted preflight requires a single nodeCa for each component";
                break;
            }
            if (!analysis.acceptsPublishedProblem() || !analysis.readyForLocalUpdater()) {
                if (event != nullptr) {
                    ++event->mergeNodeLocalProblem1PreflightRejectedComponents;
                    if (analysis.acceptsPublishedProblem()) {
                        ++event->mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation;
                    } else {
                        ++event->mergeNodeLocalProblem1PreflightRejectedByPublishedProblem;
                    }
                    if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                        event->mergeNodeLocalProblem1FailureMode =
                            classifyProblem1ExplicitRegionFailureMode(analysis);
                    }
                }

                problem1RestrictedCandidateAccepted = false;
                problem1RestrictedRejection =
                    problem1ExplicitRegionPreflightMessage(analysis.preflight.status);
                break;
            }
        }

        int totalComponentPixels = 0;
        for (const auto &component : components) {
            totalComponentPixels += static_cast<int>(component.size());
        }

        if (event != nullptr) {
            event->mergeNodeLocalAttempted = true;
            event->mergeNodeLocalDualMatchesPrimalBeforeMerge =
                dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
        }

        if (!problem1RestrictedCandidateAccepted) {
            if (event != nullptr) {
                event->mergeNodeLocalRejectedComponents =
                    event->mergeNodeLocalProblem1PreflightRejectedComponents;
                if (event->mergeNodeLocalRejectedComponents <= 0) {
                    event->mergeNodeLocalRejectedComponents = 1;
                }
            }
            fallbackToPrimalMergeWithPhaseRebuild(problem1RestrictedRejection.empty()
                                                      ? "problem1-close-recompose preflight rejected the component set"
                                                      : problem1RestrictedRejection);
            if (event != nullptr && event->mergeNodeLocalFailureSummary.empty()) {
                event->mergeNodeLocalFailureSummary =
                    "problem1-close-recompose preflight rejected the component set";
                if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                    event->mergeNodeLocalProblem1FailureMode = "contract_rejection";
                }
            }
            return;
        }

        try {
            const auto copyStart = std::chrono::steady_clock::now();
            DynamicComponentTree primalTrial = primalTree;
            DynamicComponentTree dualTrial = dualTree;
            const auto copyEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->mergeNodeLocalCopyMicros =
                    nonIncreasingElapsedMicros(copyStart, copyEnd);
            }

            bool localMatchesAfterTrial = false;
            auto trialStart = std::chrono::steady_clock::now();
            const auto localUpdate = applyOracleGuidedLocalMergeUpdate(
                primalTrial,
                dualTrial,
                adjacency,
                nodeId);
            const auto trialEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->mergeNodeLocalApplied = localUpdate.localRecompositionApplied;
                event->mergeNodeLocalUpdateMicros =
                    nonIncreasingElapsedMicros(trialStart, trialEnd);
                event->mergeNodeLocalComponentCount =
                    static_cast<int>(components.size());
                event->mergeNodeLocalPixelCount = totalComponentPixels;
                event->mergeNodeLocalDualMatchesPrimalAfterMerge =
                    localUpdate.dualMatchesPrimalAfterMerge;
                event->oracleGuidedLocalApplied = localUpdate.localRecompositionApplied;
                event->oracleGuidedLocalMatchesTarget = localUpdate.localRecompositionMatchesTarget;
                event->oracleGuidedLocalStructureMatchesOracle = localUpdate.localStructureMatchesOracle;
                event->oracleGuidedDualMatchesPrimalAfterMerge =
                    localUpdate.dualMatchesPrimalAfterMerge;
                event->oracleGuidedRecompositionNodes = localUpdate.recompositionNodeCount;
                event->oracleGuidedCreatedNodes = localUpdate.createdNodeCount;
                event->oracleGuidedReleasedSourceOwners = localUpdate.releasedSourceOwnerCount;
                event->oracleGuidedStructureMismatchSummary =
                    localUpdate.localStructureMismatchSummary;
            }

            localMatchesAfterTrial = localUpdate.localRecompositionApplied &&
                localUpdate.localRecompositionMatchesTarget &&
                localUpdate.dualMatchesPrimalAfterMerge;
            if (localMatchesAfterTrial) {
                primalTree = std::move(primalTrial);
                dualTree = std::move(dualTrial);
                if (event != nullptr) {
                    event->mergeNodeLocalAccepted = true;
                    event->applied = true;
                    event->mergeNodeLocalFallbackToPhaseRebuild = false;
                }
                return;
            }

            if (event != nullptr) {
                event->mergeNodeLocalAccepted = false;
                event->mergeNodeLocalFailureSummary =
                    !localUpdate.localRecompositionMatchesTarget
                        ? "problem1-close-recompose local recomposition did not match target"
                        : "problem1-close-recompose local recomposition did not match primal reconstruction";
            }
        } catch (const std::exception &e) {
            if (event != nullptr) {
                event->mergeNodeLocalFailureSummary =
                    std::string("problem1-close-recompose trial update threw: ") + e.what();
                ++event->mergeNodeLocalRejectedComponents;
                if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                    event->mergeNodeLocalProblem1FailureMode = "trial_exception";
                }
            }
        }

        if (event != nullptr && event->mergeNodeLocalRejectedComponents == 0) {
            ++event->mergeNodeLocalRejectedComponents;
        }

        fallbackToPrimalMergeWithPhaseRebuild(
            event != nullptr && !event->mergeNodeLocalFailureSummary.empty()
                ? event->mergeNodeLocalFailureSummary
                : "problem1-close-recompose trial was rejected; fallback to phase rebuild");
        if (event != nullptr && event->mergeNodeLocalFailureSummary.empty()) {
            event->mergeNodeLocalFailureSummary =
                "problem1-close-recompose trial was rejected; fallback to phase rebuild";
            if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                event->mergeNodeLocalProblem1FailureMode = "trial_rejected_fallback";
            }
        }
        return;
    }

    if (dualMaintenance == NonIncreasingDualMaintenance::OracleGuidedLocalRecomposition) {
        if (useOracleLocalLeafExplicitFastPath && primalTree.getNumChildren(nodeId) == 0) {
            for (const auto &component : components) {
                if (!component.empty()) {
                    adjust.updateTreeWithExplicitRegion(&dualTree, std::span<const PixelId>(component), b);
                }
            }
            primalTree.mergeNodeIntoParent(nodeId);
            (void) releaseAllEmptyLeafNodes(primalTree);
            (void) releaseAllEmptyLeafNodes(dualTree);
            if (!dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage())) {
                dualTree.build(primalTree.reconstructionImage()->clone(), dualTree.isMaxtree(), adjacency);
                if (event != nullptr) {
                    event->dualRebuiltAfterLocalMismatch = true;
                }
            }
            if (event != nullptr) {
                event->oracleGuidedDualMatchesPrimalAfterMerge =
                    dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
                event->applied = true;
            }
            return;
        }

        const auto localUpdate = applyOracleGuidedLocalMergeUpdate(
            primalTree,
            dualTree,
            adjacency,
            nodeId);
        if (event != nullptr) {
            event->oracleGuidedLocalApplied = localUpdate.localRecompositionApplied;
            event->oracleGuidedLocalMatchesTarget = localUpdate.localRecompositionMatchesTarget;
            event->oracleGuidedLocalStructureMatchesOracle = localUpdate.localStructureMatchesOracle;
            event->oracleGuidedDualMatchesPrimalAfterMerge = localUpdate.dualMatchesPrimalAfterMerge;
            event->oracleGuidedRecompositionNodes = localUpdate.recompositionNodeCount;
            event->oracleGuidedCreatedNodes = localUpdate.createdNodeCount;
            event->oracleGuidedReleasedSourceOwners = localUpdate.releasedSourceOwnerCount;
            event->oracleGuidedStructureMismatchSummary =
                localUpdate.localStructureMismatchSummary;
            event->applied = localUpdate.applied;
        }
        return;
    }

    if (dualMaintenance == NonIncreasingDualMaintenance::MergeNodeUpdateThenPhaseRebuild ||
        dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild ||
        dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild ||
        dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild ||
        dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCleanUnchecked) {
        if (event != nullptr) {
            event->mergeNodeLocalAttempted = true;
        }
        NonIncreasingDirectComponentScratch localReleaseScratch;
        auto &releaseScratch =
            componentScratch == nullptr ? localReleaseScratch : *componentScratch;
        std::vector<NodeId> *mergeNodeDualReleaseFocusNodes = nullptr;
        if (dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCleanUnchecked) {
            releaseScratch.releaseFocusMarks.ensureSize(
                static_cast<std::size_t>(dualTree.getNumInternalNodeSlots()));
            releaseScratch.releaseFocusMarks.resetAll();
            releaseScratch.dualReleaseFocusNodes.clear();
            releaseScratch.dualReleaseFocusNodes.reserve(components.size());
            mergeNodeDualReleaseFocusNodes = &releaseScratch.dualReleaseFocusNodes;
        }
        if (collectDetailedEventDiagnostics ||
            dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCleanUnchecked) {
            for (const auto &component : components) {
                if (component.empty()) {
                    continue;
                }
                NodeId nodeCa = InvalidNode;
                if (dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCleanUnchecked) {
                    nodeCa = findNodeCaAndAppendExplicitRegionOwnerFocusNodes(
                        dualTree,
                        std::span<const PixelId>(component),
                        *mergeNodeDualReleaseFocusNodes,
                        &releaseScratch.releaseFocusMarks);
                    addUniqueReleaseFocusNode(
                        *mergeNodeDualReleaseFocusNodes,
                        nodeCa,
                        &releaseScratch.releaseFocusMarks);
                } else {
                    nodeCa = findNodeCaForExplicitRegion(
                        dualTree,
                        std::span<const PixelId>(component));
                }
                if (collectDetailedEventDiagnostics) {
                    event->componentNodeCa.push_back(nodeCa);
                    event->componentAltitudeCa.push_back(
                        nodeCa == InvalidNode ? -1 : dualTree.getAltitude(nodeCa));
                }
            }
        }

        auto mergePrimalOnlyForMergeNodePhaseFallback = [&](const std::string &message) {
            const auto fallbackStart = std::chrono::steady_clock::now();
            if (nodeId != InvalidNode &&
                primalTree.isNode(nodeId) &&
                primalTree.isAlive(nodeId) &&
                nodeId != primalTree.getRoot()) {
                primalTree.mergeNodeIntoParent(nodeId);
            }
            const auto fallbackMergeEnd = std::chrono::steady_clock::now();
            const auto fallbackReleaseStart = std::chrono::steady_clock::now();
            (void) releaseAllEmptyLeafNodes(primalTree);
            const auto fallbackEnd = std::chrono::steady_clock::now();
            if (event != nullptr) {
                event->mergeNodeLocalFallbackToPhaseRebuild = true;
                event->mergeNodePrimalMergeMicros =
                    nonIncreasingElapsedMicros(fallbackStart, fallbackMergeEnd);
                event->mergeNodeEmptyReleaseMicros =
                    nonIncreasingElapsedMicros(fallbackReleaseStart, fallbackEnd);
                event->applied = true;
                if (event->mergeNodeLocalFailureSummary.empty()) {
                    event->mergeNodeLocalFailureSummary = message;
                }
            }
        };

        const auto attemptMergeNodePhaseRebuildFallback = [&](const std::string &message) {
            if (!allowMergeNodePhaseRebuildFallback) {
                if (event != nullptr && event->mergeNodeLocalFailureSummary.empty()) {
                    event->mergeNodeLocalFailureSummary = message;
                }
                if (event != nullptr) {
                    event->mergeNodeLocalSkippedByPhaseRebuild = true;
                }
                return false;
            }
            mergePrimalOnlyForMergeNodePhaseFallback(message);
            return true;
        };

        if (dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCleanUnchecked) {
            try {
                const auto directUpdate = adjust.mergeNodeIntoParentAndUpdate(
                    &primalTree,
                    &dualTree,
                    nodeId,
                    components,
                    event != nullptr
                        ? DualMinMaxTreeIncrementalFilter<AltitudeType>::MergeNodeTimingPolicy::Collect
                        : DualMinMaxTreeIncrementalFilter<AltitudeType>::MergeNodeTimingPolicy::Skip);
                const bool directAccepted =
                    directUpdate.localUpdateApplied &&
                    directUpdate.primalMergeApplied;
                if (directAccepted && attributeDirtyRecorder != nullptr) {
                    appendNonIncreasingDirtyNode(
                        attributeDirtyRecorder->primalDirtySeeds,
                        parentId);
                    appendNonIncreasingDirtyNode(
                        attributeDirtyRecorder->primalReleasedNodes,
                        nodeId);
                    if (attributeDirtyRecorder->dualDirtySeeds != nullptr) {
                        attributeDirtyRecorder->dualDirtySeeds->insert(
                            attributeDirtyRecorder->dualDirtySeeds->end(),
                            directUpdate.dualDirtySeeds.begin(),
                            directUpdate.dualDirtySeeds.end());
                    }
                    if (attributeDirtyRecorder->dualReleasedNodes != nullptr) {
                        attributeDirtyRecorder->dualReleasedNodes->insert(
                            attributeDirtyRecorder->dualReleasedNodes->end(),
                            directUpdate.dualReleasedNodes.begin(),
                            directUpdate.dualReleasedNodes.end());
                    }
                }
                std::chrono::steady_clock::time_point releaseStart;
                if (event != nullptr) {
                    releaseStart = std::chrono::steady_clock::now();
                }
                releaseScratch.primalReleaseFocusNodes.clear();
                releaseScratch.primalReleaseFocusNodes.push_back(parentId);
                (void) releaseFocusedEmptyNonRootNodes(
                    primalTree,
                    std::span<const NodeId>(releaseScratch.primalReleaseFocusNodes),
                    std::span<const NodeId>(),
                    &releaseScratch.releaseCandidateMarks,
                    &releaseScratch.releaseWorklist,
                    attributeDirtyRecorder == nullptr ? nullptr : attributeDirtyRecorder->primalDirtySeeds,
                    attributeDirtyRecorder == nullptr ? nullptr : attributeDirtyRecorder->primalReleasedNodes);
                (void) releaseFocusedEmptyNonRootNodes(
                    dualTree,
                    mergeNodeDualReleaseFocusNodes == nullptr
                        ? std::span<const NodeId>(directUpdate.dualDirtySeeds.data(), directUpdate.dualDirtySeeds.size())
                        : std::span<const NodeId>(*mergeNodeDualReleaseFocusNodes),
                    std::span<const NodeId>(),
                    &releaseScratch.releaseCandidateMarks,
                    &releaseScratch.releaseWorklist,
                    attributeDirtyRecorder == nullptr ? nullptr : attributeDirtyRecorder->dualDirtySeeds,
                    attributeDirtyRecorder == nullptr ? nullptr : attributeDirtyRecorder->dualReleasedNodes);
                long long releaseMicros = 0;
                if (event != nullptr) {
                    releaseMicros =
                        nonIncreasingElapsedMicros(releaseStart, std::chrono::steady_clock::now());
                }
                if (event != nullptr) {
                    event->mergeNodeLocalApplied = directAccepted;
                    event->mergeNodeLocalAccepted = directAccepted;
                    event->mergeNodeLocalDirectCommit = true;
                    event->mergeNodeLocalTargetCreationAllowed = true;
                    event->mergeNodeLocalComponentCount = directUpdate.componentCount;
                    event->mergeNodeLocalPixelCount = directUpdate.regionSize;
                    event->mergeNodeLocalUpdateMicros =
                        directUpdate.collectComponentsUs +
                        directUpdate.localUpdateUs;
                    event->mergeNodePrimalMergeMicros = directUpdate.primalMergeUs;
                    event->mergeNodeEmptyReleaseMicros = releaseMicros;
                    if (!directAccepted) {
                        event->mergeNodeLocalFailureSummary =
                            directUpdate.message.empty()
                                ? "merge-node clean unchecked direct update did not apply"
                                : directUpdate.message;
                    }
                    event->applied = directAccepted;
                }
                return;
            } catch (const std::exception &e) {
                if (event == nullptr) {
                    throw;
                }
                if (event != nullptr) {
                    event->mergeNodeLocalFailureSummary =
                        std::string("merge-node clean unchecked direct update threw: ") +
                        e.what();
                }
                return;
            }

        }

        const bool certifiedCreateBMode =
            dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild;
        const bool certifiedCreateBMcpcMode =
            dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild;
        if (dualMaintenance == NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild ||
            certifiedCreateBMode ||
            certifiedCreateBMcpcMode) {
            const auto certification = certifyMergeNodeExplicitComponents(
                primalTree,
                dualTree,
                *adjacency,
                components,
                b,
                certifiedCreateBMode || certifiedCreateBMcpcMode);
            if (event != nullptr) {
                event->mergeNodeLocalCertificationChecked = certification.checked;
                event->mergeNodeLocalCertificationPassed = certification.passed;
                event->mergeNodeLocalTargetCreationAllowed =
                    certifiedCreateBMode || certifiedCreateBMcpcMode;
                event->mergeNodeLocalTargetCreationRequired =
                    (certifiedCreateBMode || certifiedCreateBMcpcMode) &&
                    certification.missingTargetComponents > 0;
                event->mergeNodeLocalComponentCount = certification.componentCount;
                event->mergeNodeLocalPixelCount = certification.pixelCount;
                event->mergeNodeLocalRejectedComponents = certification.rejectedComponents;
                event->mergeNodeLocalPartialOwnerComponents =
                    certification.partialOwnerComponents;
                event->mergeNodeLocalRawMissingTargetComponents =
                    certification.rawMissingTargetComponents;
                event->mergeNodeLocalMissingTargetComponents =
                    certification.missingTargetComponents;
                event->mergeNodeLocalPartialOwnerAndMissingTargetComponents =
                    certification.partialOwnerAndMissingTargetComponents;
                event->mergeNodeLocalCertificationMicros =
                    certification.certificationMicros;
                event->mergeNodeLocalDualMatchesPrimalBeforeMerge =
                    certification.dualMatchesPrimalBeforeMerge;
                if (!certification.passed) {
                    event->mergeNodeLocalFailureSummary = certification.message;
                }
            }
            if (!certification.passed) {
                if (certifiedCreateBMcpcMode &&
                    certification.partialOwnerComponents > 0) {
                    std::optional<NonIncreasingMcpcDirectPreflightArtifacts>
                        mcpcPreflightArtifacts;
                    NonIncreasingMcpcNoCopyCandidateDiagnostics noCopyCandidate;
                    if (mcpcDirectPreflightDiagnosticsEnabled ||
                        mcpcDirectV1StrictPrecheckEnabled ||
                        mcpcNoCopyDirectCommitEnabled ||
                        mcpcValidatedDirectCommitEnabled) {
                        const auto preflightStart = std::chrono::steady_clock::now();
                        mcpcPreflightArtifacts = buildMcpcDirectPreflightArtifacts(
                            primalTree,
                            dualTree,
                            adjacency,
                            nodeId);
                        const auto preflightEnd = std::chrono::steady_clock::now();
                        const auto &preflight = mcpcPreflightArtifacts->preflight;
                        noCopyCandidate = diagnoseMcpcNoCopyCandidate(preflight);
                        if (event != nullptr) {
                            event->mcpcDirectPreflightAttempted = true;
                            event->mcpcDirectPreflightDirectV1StrictGateEnabled =
                                mcpcDirectV1StrictPrecheckEnabled;
                            event->mcpcNoCopyDirectCommitEnabled =
                                mcpcNoCopyDirectCommitEnabled;
                            event->mcpcValidatedDirectCommitEnabled =
                                mcpcValidatedDirectCommitEnabled;
                            event->mcpcDirectPreflightChecked = preflight.checked;
                            event->mcpcDirectPreflightPassed = preflight.passed;
                            event->mcpcDirectPreflightDirectV1PrecheckEvaluated =
                                preflight.directV1PrecheckEvaluated;
                            event->mcpcDirectPreflightDirectV1PrecheckAccepted =
                                preflight.directV1PrecheckAccepted;
                            event->mcpcDirectPreflightPlanOracleFree =
                                preflight.planOracleFree;
                            event->mcpcDirectPreflightMaterializationOracleFree =
                                preflight.materializationOracleFree;
                            event->mcpcDirectPreflightMaterializationInputsValid =
                                preflight.materializationInputsValid;
                            event->mcpcDirectPreflightSummaryBuilt =
                                preflight.rollbackFreeSummaryBuilt;
                            event->mcpcDirectPreflightAffectedLevels =
                                preflight.summary.affectedLevelCount;
                            event->mcpcDirectPreflightLocalPixels =
                                preflight.summary.localPixelCount;
                            event->mcpcDirectPreflightPlannedNodes =
                                preflight.summary.plannedNodeCount;
                            event->mcpcDirectPreflightPlannedDirectProperParts =
                                preflight.summary.plannedDirectProperPartCount;
                            event->mcpcDirectPreflightPlannedExternalAttachments =
                                preflight.summary.plannedExternalAttachmentCount;
                            event->mcpcDirectPreflightSourceOwnersToBeEmptied =
                                preflight.summary.sourceOwnerToBeEmptiedCount;
                            event->mcpcDirectPreflightStableChildrenToMove =
                                preflight.summary.stableChildToMoveCount;
                            event->mcpcDirectPreflightPreviousParentsPossiblyEmptied =
                                preflight.summary.previousParentPossiblyEmptiedCount;
                            event->mcpcDirectPreflightInvalidStableChildRefs =
                                preflight.summary.invalidStableChildRefCount;
                            event->mcpcDirectPreflightStableChildCycleRisks =
                                preflight.summary.stableChildCycleRiskCount;
                            event->mcpcDirectPreflightPredictedAuxiliaryNodes =
                                preflight.summary.predictedAuxiliaryNodeCount;
                            event->mcpcDirectPreflightPredictedHiddenExternalParentConnectorNodes =
                                preflight.summary.predictedHiddenExternalParentConnectorNodeCount;
                            event->mcpcDirectPreflightPredictedStableChildBridgeNodes =
                                preflight.summary.predictedStableChildBridgeNodeCount;
                            event->mcpcDirectPreflightCandidateNonOracleHiddenSupportConnectors =
                                preflight.summary.candidateNonOracleHiddenSupportConnectorCount;
                            event->mcpcDirectPreflightCandidateNonOracleParentLinkBridges =
                                preflight.summary.candidateNonOracleParentLinkBridgeCount;
                            event->mcpcDirectPreflightCandidateNonOracleChildLinkBridges =
                                preflight.summary.candidateNonOracleChildLinkBridgeCount;
                            event->mcpcDirectPreflightV1PrecheckHiddenSupportCandidates =
                                preflight.directV1Precheck.nonOracleHiddenSupportCandidateCount;
                            event->mcpcDirectPreflightV1PrecheckHiddenSupportAlreadyRepresented =
                                preflight.directV1Precheck.nonOracleHiddenSupportAlreadyRepresentedCount;
                            event->mcpcDirectPreflightV1PrecheckHiddenSupportNotRedundant =
                                preflight.directV1Precheck.nonOracleHiddenSupportNotRedundantCount;
                            event->mcpcDirectPreflightV1PrecheckParentLinksTouchingPlan =
                                preflight.directV1Precheck.nonOracleParentLinkTouchesPlanCount;
                            event->mcpcDirectPreflightV1PrecheckParentLinksBothRepresented =
                                preflight.directV1Precheck.nonOracleParentLinkBothRepresentedCount;
                            event->mcpcDirectPreflightV1PrecheckParentLinksNotBothRepresented =
                                preflight.directV1Precheck.nonOracleParentLinkNotBothRepresentedCount;
                            event->mcpcDirectPreflightV1PrecheckParentLinksBlockedByModifiedSnapshot =
                                preflight.directV1Precheck.nonOracleParentLinkBlockedByModifiedSnapshotCount;
                            event->mcpcDirectPreflightV1PrecheckParentLinksInvalidLevel =
                                preflight.directV1Precheck.nonOracleParentLinkInvalidLevelCount;
                            event->mcpcDirectPreflightPlannedNodeLevels =
                                preflight.summary.plannedNodeLevels;
                            event->mcpcDirectPreflightMicros =
                                nonIncreasingElapsedMicros(preflightStart, preflightEnd);
                            event->mcpcDirectPreflightFailureSummary =
                                preflight.failureSummary;
                            event->mcpcDirectPreflightDirectV1PrecheckRejectionSummary =
                                preflight.directV1PrecheckRejectionSummary;
                            event->mcpcNoCopyCandidateEvaluated =
                                noCopyCandidate.evaluated;
                            event->mcpcNoCopyCandidateAccepted =
                                noCopyCandidate.accepted;
                            event->mcpcNoCopyCandidateRejectionSummary =
                                noCopyCandidate.rejectionSummary;
                        }
                    }

                    if (mcpcDirectV1StrictPrecheckEnabled &&
                        (!mcpcPreflightArtifacts ||
                         !mcpcPreflightArtifacts->preflight.directV1PrecheckAccepted)) {
                        if (event != nullptr) {
                            event->mcpcDirectPreflightDirectV1StrictGateEnabled = true;
                            event->mcpcDirectPreflightDirectV1StrictGateRejected = true;
                            event->mergeNodeLocalSkippedByPhaseRebuild = true;
                            event->mergeNodeLocalFailureSummary =
                                mcpcPreflightArtifacts
                                    ? mcpcPreflightArtifacts->preflight.directV1PrecheckRejectionSummary
                                    : "P1-MCPC direct-v1 strict precheck did not run";
                        }
                        if (!attemptMergeNodePhaseRebuildFallback(
                                mcpcPreflightArtifacts
                                    ? mcpcPreflightArtifacts->preflight.directV1PrecheckRejectionSummary
                                    : "P1-MCPC direct-v1 strict precheck did not run")) {
                            return;
                        }
                        return;
                    }

                    NonIncreasingStrippedLocalUpdateOptions mcpcOptions =
                        strippedLocalUpdateOptions;
                    mcpcOptions.releaseEmptyNodesAfterPrimalMerge = true;

                    auto mcpcUpdateAccepted =
                        [&](const NonIncreasingStrippedLocalUpdateResult &update) {
                            const bool finalCheckOk =
                                !mcpcOptions.checkFinalDualPrimalReconstruction ||
                                update.dualMatchesPrimalAfterMerge;
                            const bool targetCheckOk =
                                update.skippedNoImageChange ||
                                !mcpcOptions.checkMaterializationTargetReconstruction ||
                                update.localRecompositionMatchesTarget;
                            const bool planAuditOk =
                                !mcpcOptions.auditPlanOracleDependencies ||
                                update.localPlanOracleFree;
                            const bool materializationAuditOk =
                                !mcpcOptions.auditMaterializationOracleDependencies ||
                                update.materializationInputOracleFree;
                            const bool noEmptyCreatedNodes =
                                update.emptyCreatedNodesAfterCanonicalization.empty();
                            const bool canonicalByConstruction =
                                update.postCanonicalizationChangesRequired == 0;
                            return update.applied &&
                                   finalCheckOk &&
                                   (update.skippedNoImageChange ||
                                    (targetCheckOk &&
                                     planAuditOk &&
                                     materializationAuditOk &&
                                     noEmptyCreatedNodes &&
                                     canonicalByConstruction));
                        };

                    if (mcpcNoCopyDirectCommitEnabled &&
                        mcpcPreflightArtifacts &&
                        noCopyCandidate.accepted) {
                        bool acceptedByNoCopyDirectCommit = false;
                        bool fallbackToPhaseRebuild = false;
                        const auto directStart = std::chrono::steady_clock::now();
                        try {
                            const auto directUpdate =
                                applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight(
                                    primalTree,
                                    dualTree,
                                    adjacency,
                                    nodeId,
                                    *mcpcPreflightArtifacts,
                                    mcpcOptions);
                            acceptedByNoCopyDirectCommit =
                                mcpcUpdateAccepted(directUpdate);
                            const bool directDualMatchesPrimal =
                                directUpdate.finalReconstructionChecked
                                    ? directUpdate.dualMatchesPrimalAfterMerge
                                    : dualTree.reconstructionImage()->isEqual(
                                          primalTree.reconstructionImage());
                            fallbackToPhaseRebuild = !acceptedByNoCopyDirectCommit;
                            if (event != nullptr) {
                                event->mcpcNoCopyDirectCommitAttempted = true;
                                event->mcpcNoCopyDirectCommitAccepted =
                                    acceptedByNoCopyDirectCommit;
                                event->mcpcNoCopyDirectCommitDualMatchesPrimal =
                                    directDualMatchesPrimal;
                                event->strippedLocalAttempted = true;
                                event->strippedDirectCommit = true;
                                event->strippedLocalApplied = directUpdate.applied;
                                event->strippedLocalAccepted =
                                    acceptedByNoCopyDirectCommit;
                                event->strippedLocalPlanAuditChecked =
                                    directUpdate.localPlanOracleAuditChecked;
                                event->strippedLocalPlanOracleFree =
                                    directUpdate.localPlanOracleFree;
                                event->strippedMaterializationAuditChecked =
                                    directUpdate.materializationInputOracleAuditChecked;
                                event->strippedMaterializationInputOracleFree =
                                    directUpdate.materializationInputOracleFree;
                                event->strippedLocalTargetChecked =
                                    directUpdate.localRecompositionTargetChecked;
                                event->strippedLocalMatchesTarget =
                                    directUpdate.localRecompositionMatchesTarget;
                                event->strippedFinalReconstructionChecked =
                                    directUpdate.finalReconstructionChecked;
                                event->strippedLocalDualMatchesPrimalAfterMerge =
                                    directUpdate.dualMatchesPrimalAfterMerge;
                                event->strippedLocalRecompositionNodes =
                                    directUpdate.recompositionNodeCount;
                                event->strippedLocalCreatedNodes =
                                    directUpdate.createdNodeCount;
                                event->strippedLocalPixelCount =
                                    directUpdate.localPixelCount;
                                event->strippedCreatedConnectorNodes =
                                    directUpdate.createdConnectorNodeCount;
                                event->strippedEmptyCreatedNodesAfterCanonicalization =
                                    static_cast<int>(
                                        directUpdate.emptyCreatedNodesAfterCanonicalization.size());
                                event->strippedPostEmptyNodeReleasesRequired =
                                    directUpdate.postEmptyNodeReleasesRequired;
                                event->strippedPostSameLevelParentChildCollapsesRequired =
                                    directUpdate.postSameLevelParentChildCollapsesRequired;
                                event->strippedPostSameLevelSiblingMergesRequired =
                                    directUpdate.postSameLevelSiblingMergesRequired;
                                event->strippedPostComparableSiblingMovesRequired =
                                    directUpdate.postComparableSiblingMovesRequired;
                                event->strippedPostCanonicalizationChangesRequired =
                                    directUpdate.postCanonicalizationChangesRequired;
                                if (collectDetailedEventDiagnostics) {
                                    event->strippedEmptyCreatedNodeDiagnostics =
                                        directUpdate.emptyCreatedNodesAfterCanonicalization;
                                }
                                event->strippedLocalReleasedSourceOwners =
                                    directUpdate.releasedSourceOwnerCount;
                                event->mergeNodeLocalApplied = directUpdate.applied;
                                event->mergeNodeLocalDualMatchesPrimalAfterMerge =
                                    directDualMatchesPrimal;
                                if (!acceptedByNoCopyDirectCommit) {
                                    if (directUpdate.postCanonicalizationChangesRequired != 0) {
                                        event->mcpcNoCopyDirectCommitFailureSummary =
                                            "P1-MCPC no-copy direct commit required corrective post-canonicalization";
                                    } else {
                                        event->mcpcNoCopyDirectCommitFailureSummary =
                                            "P1-MCPC no-copy direct commit failed acceptance checks";
                                    }
                                }
                            }
                        } catch (const std::exception &e) {
                            fallbackToPhaseRebuild = true;
                            if (event != nullptr) {
                                event->mcpcNoCopyDirectCommitAttempted = true;
                                event->mcpcNoCopyDirectCommitFailureSummary =
                                    std::string("P1-MCPC no-copy direct commit threw: ") +
                                    e.what();
                                event->mergeNodeLocalFailureSummary =
                                    event->mcpcNoCopyDirectCommitFailureSummary;
                            }
                        }
                        const auto directEnd = std::chrono::steady_clock::now();
                        if (event != nullptr) {
                            event->mcpcNoCopyDirectCommitMicros =
                                nonIncreasingElapsedMicros(directStart, directEnd);
                        }

                        if (acceptedByNoCopyDirectCommit) {
                            if (event != nullptr) {
                                event->mcpcNoCopyDirectCommitAccepted = true;
                                event->mergeNodeLocalAccepted = true;
                                event->mergeNodeLocalApplied = true;
                                event->mergeNodeLocalDirectCommit = true;
                                event->strippedDirectCommit = true;
                                event->mergeNodeLocalFailureSummary.clear();
                                event->strippedLocalFailureSummary.clear();
                                event->applied = true;
                            }
                            return;
                        }

                        if (fallbackToPhaseRebuild) {
                            if (event != nullptr) {
                                event->mcpcNoCopyDirectCommitFellBackToPhaseRebuild = true;
                            }
                            if (!attemptMergeNodePhaseRebuildFallback(
                                    event != nullptr &&
                                        !event->mcpcNoCopyDirectCommitFailureSummary.empty()
                                        ? event->mcpcNoCopyDirectCommitFailureSummary
                                        : "P1-MCPC no-copy direct commit fell back to phase rebuild")) {
                                return;
                            }
                            return;
                        }
                    }

                    const auto copyStart = std::chrono::steady_clock::now();
                    DynamicComponentTree primalTrial = primalTree;
                    DynamicComponentTree dualTrial = dualTree;
                    const auto copyEnd = std::chrono::steady_clock::now();
                    if (event != nullptr) {
                        event->strippedLocalAttempted = true;
                        event->strippedDirectCommit = false;
                        event->strippedLocalCopyMicros =
                            nonIncreasingElapsedMicros(copyStart, copyEnd);
                    }

                    bool mcpcAccepted = false;
                    NonIncreasingMcpcDirectSummaryComparison summaryComparison;
                    NonIncreasingMcpcDirectV1AcceptanceDiagnostics directV1;
                    NonIncreasingMcpcDirectCommitCandidateDiagnostics commitCandidate;
                    const auto mcpcStart = std::chrono::steady_clock::now();
                    try {
                        NonIncreasingStrippedLocalUpdateResult strippedUpdate;
                        if (mcpcPreflightArtifacts &&
                            mcpcPreflightArtifacts->preflight.passed) {
                            if (event != nullptr) {
                                event->mcpcDirectPreflightSharedWithStrippedUpdate = true;
                            }
                            strippedUpdate =
                                applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight(
                                    primalTrial,
                                    dualTrial,
                                    adjacency,
                                    nodeId,
                                    *mcpcPreflightArtifacts,
                                    mcpcOptions);
                        } else {
                            strippedUpdate = applyStrippedNonOracleLocalMergeUpdate(
                                primalTrial,
                                dualTrial,
                                adjacency,
                                nodeId,
                                mcpcOptions);
                        }
                        mcpcAccepted = mcpcUpdateAccepted(strippedUpdate);
                        if (mcpcPreflightArtifacts) {
                            summaryComparison =
                                compareMcpcDirectPreflightSummaryWithStrippedUpdate(
                                    mcpcPreflightArtifacts->preflight.summary,
                                    strippedUpdate);
                            directV1 =
                                diagnoseMcpcDirectPreflightV1Acceptance(
                                    mcpcPreflightArtifacts->preflight,
                                    strippedUpdate);
                            commitCandidate =
                                diagnoseMcpcDirectCommitCandidate(
                                    mcpcPreflightArtifacts->preflight,
                                    summaryComparison,
                                    directV1,
                                    strippedUpdate);
                        }
                        if (event != nullptr) {
                            event->strippedLocalApplied = strippedUpdate.applied;
                            event->strippedLocalAccepted = mcpcAccepted;
                            event->strippedLocalPlanAuditChecked =
                                strippedUpdate.localPlanOracleAuditChecked;
                            event->strippedLocalPlanOracleFree =
                                strippedUpdate.localPlanOracleFree;
                            event->strippedMaterializationAuditChecked =
                                strippedUpdate.materializationInputOracleAuditChecked;
                            event->strippedMaterializationInputOracleFree =
                                strippedUpdate.materializationInputOracleFree;
                            event->strippedLocalTargetChecked =
                                strippedUpdate.localRecompositionTargetChecked;
                            event->strippedLocalMatchesTarget =
                                strippedUpdate.localRecompositionMatchesTarget;
                            event->strippedFinalReconstructionChecked =
                                strippedUpdate.finalReconstructionChecked;
                            event->strippedLocalDualMatchesPrimalAfterMerge =
                                strippedUpdate.dualMatchesPrimalAfterMerge;
                            event->strippedLocalRecompositionNodes =
                                strippedUpdate.recompositionNodeCount;
                            event->strippedLocalCreatedNodes =
                                strippedUpdate.createdNodeCount;
                            event->strippedLocalPixelCount =
                                strippedUpdate.localPixelCount;
                            event->strippedCreatedConnectorNodes =
                                strippedUpdate.createdConnectorNodeCount;
                            event->strippedCreatedHiddenExternalParentConnectorNodes =
                                strippedUpdate.createdHiddenExternalParentConnectorNodeCount;
                            event->strippedCreatedOracleBridgeNodes =
                                strippedUpdate.createdOracleBridgeNodeCount;
                            event->strippedCreatedStableChildBridgeNodes =
                                strippedUpdate.createdStableChildBridgeNodeCount;
                            event->strippedCreatedHiddenSupportConnectorNodes =
                                strippedUpdate.createdHiddenSupportConnectorNodeCount;
                            event->strippedCreatedNonOracleHiddenSupportConnectorNodes =
                                strippedUpdate.createdNonOracleHiddenSupportConnectorNodeCount;
                            event->strippedCreatedNonOracleParentLinkBridgeNodes =
                                strippedUpdate.createdNonOracleParentLinkBridgeNodeCount;
                            event->strippedCreatedNonOracleChildLinkBridgeNodes =
                                strippedUpdate.createdNonOracleChildLinkBridgeNodeCount;
                            event->strippedNonOracleHiddenSupportCandidates =
                                strippedUpdate.nonOracleHiddenSupportCandidateCount;
                            event->strippedNonOracleHiddenSupportAlreadyRepresented =
                                strippedUpdate.nonOracleHiddenSupportAlreadyRepresentedCount;
                            event->strippedNonOracleHiddenSupportMaterializable =
                                strippedUpdate.nonOracleHiddenSupportMaterializableCount;
                            event->strippedNonOracleParentLinksTouchingPlan =
                                strippedUpdate.nonOracleParentLinkTouchesPlanCount;
                            event->strippedNonOracleParentLinksBothRepresented =
                                strippedUpdate.nonOracleParentLinkBothRepresentedCount;
                            event->strippedNonOracleParentLinksBlockedByModifiedSnapshot =
                                strippedUpdate.nonOracleParentLinkBlockedByModifiedSnapshotCount;
                            event->strippedNonOracleParentLinksInvalidLevel =
                                strippedUpdate.nonOracleParentLinkInvalidLevelCount;
                            event->strippedNonOracleParentLinksAlreadyDirect =
                                strippedUpdate.nonOracleParentLinkAlreadyDirectCount;
                            event->strippedNonOracleParentLinksAlreadyIntermediate =
                                strippedUpdate.nonOracleParentLinkAlreadyIntermediateCount;
                            event->strippedNonOracleParentLinksMoved =
                                strippedUpdate.nonOracleParentLinkMovedCount;
                            event->strippedNonOracleParentLinkBridgeAttempts =
                                strippedUpdate.nonOracleParentLinkBridgeAttemptCount;
                            event->strippedNonOracleParentLinkBridgeMaterializable =
                                strippedUpdate.nonOracleParentLinkBridgeMaterializableCount;
                            event->strippedNonOracleChildLinkBridgeAttempts =
                                strippedUpdate.nonOracleChildLinkBridgeAttemptCount;
                            event->strippedNonOracleChildLinkBridgeInvalidLevel =
                                strippedUpdate.nonOracleChildLinkBridgeInvalidLevelCount;
                            event->strippedNonOracleChildLinkBridgeMaterializable =
                                strippedUpdate.nonOracleChildLinkBridgeMaterializableCount;
                            event->strippedEmptyCreatedNodesAfterCanonicalization =
                                static_cast<int>(
                                    strippedUpdate.emptyCreatedNodesAfterCanonicalization.size());
                            event->strippedPostEmptyNodeReleasesRequired =
                                strippedUpdate.postEmptyNodeReleasesRequired;
                            event->strippedPostSameLevelParentChildCollapsesRequired =
                                strippedUpdate.postSameLevelParentChildCollapsesRequired;
                            event->strippedPostSameLevelSiblingMergesRequired =
                                strippedUpdate.postSameLevelSiblingMergesRequired;
                            event->strippedPostComparableSiblingMovesRequired =
                                strippedUpdate.postComparableSiblingMovesRequired;
                            event->strippedPostCanonicalizationChangesRequired =
                                strippedUpdate.postCanonicalizationChangesRequired;
                            if (collectDetailedEventDiagnostics) {
                                event->strippedCreatedNodeLevels =
                                    strippedUpdate.createdNodeLevels;
                                event->strippedCreatedNodePlanIndices =
                                    strippedUpdate.createdNodePlanIndices;
                                event->strippedEmptyCreatedNodeDiagnostics =
                                    strippedUpdate.emptyCreatedNodesAfterCanonicalization;
                                event->strippedEmptyConnectorNodesAfterCanonicalization =
                                    static_cast<int>(std::count_if(
                                        event->strippedEmptyCreatedNodeDiagnostics.begin(),
                                        event->strippedEmptyCreatedNodeDiagnostics.end(),
                                        [](const auto &node) {
                                            return node.planIndex < 0;
                                        }));
                            }
                            event->strippedLocalReleasedSourceOwners =
                                strippedUpdate.releasedSourceOwnerCount;
                            event->strippedInputExtractionMicros =
                                strippedUpdate.inputExtractionMicros;
                            event->strippedPartitionPlanMicros =
                                strippedUpdate.partitionPlanMicros;
                            event->strippedRecompositionPlanMicros =
                                strippedUpdate.recompositionPlanMicros;
                            event->strippedPlanAuditMicros =
                                strippedUpdate.planAuditMicros;
                            event->strippedMaterializationContextMicros =
                                strippedUpdate.materializationContextMicros;
                            event->strippedMaterializationAuditMicros =
                                strippedUpdate.materializationAuditMicros;
                            event->strippedMaterializationMicros =
                                strippedUpdate.materializationMicros;
                            event->strippedPostLocalSnapshotMicros =
                                strippedUpdate.postLocalSnapshotMicros;
                            event->strippedPrimalMergeReleaseMicros =
                                strippedUpdate.primalMergeReleaseMicros;
                            event->strippedFinalReconstructionCheckMicros =
                                strippedUpdate.finalReconstructionCheckMicros;
                            if (mcpcPreflightArtifacts) {
                                event->mcpcDirectPreflightSummaryComparedWithStrippedUpdate =
                                    summaryComparison.evaluated;
                                event->mcpcDirectPreflightSummaryMatchesStrippedUpdate =
                                    summaryComparison.matches;
                                event->mcpcDirectPreflightSummaryMismatchCount =
                                    summaryComparison.mismatchCount;
                                event->mcpcDirectPreflightSummaryPlannedNodeCountMismatch =
                                    summaryComparison.plannedNodeCountMismatch;
                                event->mcpcDirectPreflightSummaryDirectProperPartCountMismatch =
                                    summaryComparison.plannedDirectProperPartCountMismatch;
                                event->mcpcDirectPreflightSummaryAuxiliaryNodeCountMismatch =
                                    summaryComparison.predictedAuxiliaryNodeCountMismatch;
                                event->mcpcDirectPreflightSummaryHiddenExternalParentConnectorMismatch =
                                    summaryComparison.predictedHiddenExternalParentConnectorMismatch;
                                event->mcpcDirectPreflightSummaryStableChildBridgeMismatch =
                                    summaryComparison.predictedStableChildBridgeMismatch;
                                event->mcpcDirectPreflightSummaryEmptiedLocalNodeMismatch =
                                    summaryComparison.emptiedLocalNodeCountMismatch;
                                event->mcpcDirectPreflightSummaryMovedStableChildMismatch =
                                    summaryComparison.movedStableChildCountMismatch;
                                event->mcpcDirectPreflightPredictedEmptiedLocalNodes =
                                    summaryComparison.predictedEmptiedLocalNodeCount;
                                event->mcpcDirectPreflightObservedEmptiedLocalNodes =
                                    summaryComparison.observedEmptiedLocalNodeCount;
                                event->mcpcDirectPreflightSummaryMismatchSummary =
                                    summaryComparison.mismatchSummary;

                                event->mcpcDirectPreflightDirectV1Evaluated =
                                    directV1.evaluated;
                                event->mcpcDirectPreflightDirectV1Accepted =
                                    directV1.accepted;
                                event->mcpcDirectPreflightDirectV1RejectionSummary =
                                    directV1.rejectionSummary;

                                event->mcpcDirectCommitCandidateEvaluated =
                                    commitCandidate.evaluated;
                                event->mcpcDirectCommitCandidateAccepted =
                                    commitCandidate.accepted;
                                event->mcpcDirectCommitCandidateRejectionSummary =
                                    commitCandidate.rejectionSummary;
                            }
                        }
                    } catch (const std::exception &e) {
                        if (event != nullptr) {
                            event->strippedLocalFailureSummary = e.what();
                        }
                    }
                    const auto mcpcEnd = std::chrono::steady_clock::now();
                    if (event != nullptr) {
                        event->strippedLocalUpdateMicros =
                            nonIncreasingElapsedMicros(mcpcStart, mcpcEnd);
                    }

                    if (mcpcAccepted) {
                        bool acceptedByValidatedDirectCommit = false;
                        if (mcpcValidatedDirectCommitEnabled &&
                            mcpcPreflightArtifacts &&
                            commitCandidate.accepted) {
                            if (event != nullptr) {
                                event->mcpcValidatedDirectCommitEnabled = true;
                                event->mcpcValidatedDirectCommitAttempted = true;
                            }
                            const auto directStart = std::chrono::steady_clock::now();
                            try {
                                const auto directUpdate =
                                    applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight(
                                        primalTree,
                                        dualTree,
                                        adjacency,
                                        nodeId,
                                        *mcpcPreflightArtifacts,
                                        mcpcOptions);
                                const bool directUpdateAccepted =
                                    mcpcUpdateAccepted(directUpdate);
                                const bool directDualMatchesPrimal =
                                    directUpdate.finalReconstructionChecked
                                        ? directUpdate.dualMatchesPrimalAfterMerge
                                        : dualTree.reconstructionImage()->isEqual(
                                              primalTree.reconstructionImage());
                                const bool directMatchesCopiedOutput =
                                    primalTree.reconstructionImage()->isEqual(
                                        primalTrial.reconstructionImage()) &&
                                    dualTree.reconstructionImage()->isEqual(
                                        dualTrial.reconstructionImage());
                                acceptedByValidatedDirectCommit =
                                    directUpdateAccepted &&
                                    directDualMatchesPrimal &&
                                    directMatchesCopiedOutput;
                                if (event != nullptr) {
                                    event->mcpcValidatedDirectCommitAccepted =
                                        acceptedByValidatedDirectCommit;
                                    event->mcpcValidatedDirectCommitDualMatchesPrimal =
                                        directDualMatchesPrimal;
                                    event->mcpcValidatedDirectCommitMatchesCopiedOutput =
                                        directMatchesCopiedOutput;
                                    if (!acceptedByValidatedDirectCommit) {
                                        event->mcpcValidatedDirectCommitFailureSummary =
                                            "validated P1-MCPC direct commit did not match copied validation";
                                    }
                                }
                            } catch (const std::exception &e) {
                                if (event != nullptr) {
                                    event->mcpcValidatedDirectCommitFailureSummary =
                                        std::string("validated P1-MCPC direct commit threw: ") +
                                        e.what();
                                }
                            }
                            const auto directEnd = std::chrono::steady_clock::now();
                            if (event != nullptr) {
                                event->mcpcValidatedDirectCommitMicros =
                                    nonIncreasingElapsedMicros(directStart, directEnd);
                            }
                        }

                        if (!acceptedByValidatedDirectCommit) {
                            if (event != nullptr &&
                                event->mcpcValidatedDirectCommitAttempted) {
                                event->mcpcValidatedDirectCommitFellBackToCopy = true;
                            }
                            primalTree = std::move(primalTrial);
                            dualTree = std::move(dualTrial);
                        }
                        if (event != nullptr) {
                            event->mergeNodeLocalAccepted = true;
                            event->mergeNodeLocalApplied = true;
                            event->mergeNodeLocalDualMatchesPrimalAfterMerge =
                                acceptedByValidatedDirectCommit
                                    ? event->mcpcValidatedDirectCommitDualMatchesPrimal
                                    : event->strippedLocalDualMatchesPrimalAfterMerge;
                            event->mergeNodeLocalDirectCommit =
                                acceptedByValidatedDirectCommit;
                            event->strippedDirectCommit =
                                acceptedByValidatedDirectCommit;
                            event->mergeNodeLocalFailureSummary.clear();
                            event->strippedLocalFailureSummary.clear();
                            event->applied = true;
                        }
                        return;
                    }

                    if (attemptMergeNodePhaseRebuildFallback(certification.message)) {
                        if (event != nullptr) {
                            event->strippedLocalFallbackToPhaseRebuild = true;
                            if (event->strippedLocalFailureSummary.empty()) {
                                event->strippedLocalFailureSummary =
                                    "merge-node P1-MCPC stripped trial did not satisfy acceptance checks";
                            }
                        }
                        return;
                    }
                }
                return;
            }

            if (certifiedCreateBMode || certifiedCreateBMcpcMode) {
                const auto copyStart = std::chrono::steady_clock::now();
                DynamicComponentTree primalTrial = primalTree;
                DynamicComponentTree dualTrial = dualTree;
                const auto copyEnd = std::chrono::steady_clock::now();
                if (event != nullptr) {
                    event->mergeNodeLocalCopyMicros =
                        nonIncreasingElapsedMicros(copyStart, copyEnd);
                }

                bool copiedUpdateAccepted = false;
                const auto updateStart = std::chrono::steady_clock::now();
                try {
                    DynamicComponentTree *minTrial = dualTrial.isMaxtree() ? &primalTrial : &dualTrial;
                    DynamicComponentTree *maxTrial = dualTrial.isMaxtree() ? &dualTrial : &primalTrial;
                    DualMinMaxTreeIncrementalFilter<AltitudeType> trialAdjust(
                        minTrial,
                        maxTrial,
                        *adjacency);
                    for (const auto &component : components) {
                        if (!component.empty()) {
                            trialAdjust.updateTreeWithExplicitRegion(
                                &dualTrial,
                                std::span<const PixelId>(component),
                                b,
                                true);
                        }
                    }
                    const auto mergeStart = std::chrono::steady_clock::now();
                    primalTrial.mergeNodeIntoParent(nodeId);
                    const auto mergeEnd = std::chrono::steady_clock::now();
                    copiedUpdateAccepted =
                        dualTrial.reconstructionImage()->isEqual(primalTrial.reconstructionImage());
                    if (event != nullptr) {
                        event->mergeNodeLocalApplied = true;
                        event->mergeNodeLocalDualMatchesPrimalAfterMerge =
                            copiedUpdateAccepted;
                        event->mergeNodePrimalMergeMicros =
                            nonIncreasingElapsedMicros(mergeStart, mergeEnd);
                    }
                } catch (const std::exception &e) {
                    if (event != nullptr) {
                        event->mergeNodeLocalFailureSummary = e.what();
                    }
                }
                const auto updateEnd = std::chrono::steady_clock::now();
                if (event != nullptr) {
                    event->mergeNodeLocalUpdateMicros =
                        nonIncreasingElapsedMicros(updateStart, updateEnd) -
                        event->mergeNodePrimalMergeMicros;
                    if (event->mergeNodeLocalUpdateMicros < 0) {
                        event->mergeNodeLocalUpdateMicros = 0;
                    }
                }

                if (copiedUpdateAccepted) {
                    primalTree = std::move(primalTrial);
                    dualTree = std::move(dualTrial);
                    if (event != nullptr) {
                        event->mergeNodeLocalAccepted = true;
                        event->applied = true;
                    }
                    return;
                }

                if (!attemptMergeNodePhaseRebuildFallback(
                        event != nullptr && !event->mergeNodeLocalFailureSummary.empty()
                            ? event->mergeNodeLocalFailureSummary
                            : "certified create-b merge-node copy update did not preserve dual/primal reconstruction equality")) {
                    return;
                }
                return;
            }

            bool directAccepted = false;
            bool directPrimalMerged = false;
            try {
                const auto directUpdate = adjust.mergeNodeIntoParentAndUpdate(
                    &primalTree,
                    &dualTree,
                    nodeId);
                directPrimalMerged = directUpdate.primalMergeApplied;
                const bool dualMatchesPrimalAfterMerge =
                    directUpdate.primalMergeApplied &&
                    directUpdate.localUpdateApplied &&
                    dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
                directAccepted = dualMatchesPrimalAfterMerge;
                if (event != nullptr) {
                    event->mergeNodeLocalApplied =
                        directUpdate.localUpdateApplied && directUpdate.primalMergeApplied;
                    event->mergeNodeLocalDirectCommit = true;
                    event->mergeNodeLocalDualMatchesPrimalAfterMerge =
                        dualMatchesPrimalAfterMerge;
                    event->mergeNodeLocalComponentCount = directUpdate.componentCount;
                    event->mergeNodeLocalPixelCount = directUpdate.regionSize;
                    event->mergeNodeLocalUpdateMicros =
                        directUpdate.collectComponentsUs + directUpdate.localUpdateUs;
                    event->mergeNodePrimalMergeMicros = directUpdate.primalMergeUs;
                    if (!directAccepted) {
                        event->mergeNodeLocalFailureSummary = directUpdate.message.empty()
                            ? "certified merge-node direct update did not preserve dual/primal reconstruction equality"
                            : directUpdate.message;
                    }
                }
            } catch (const std::exception &e) {
                if (!attemptMergeNodePhaseRebuildFallback(
                        std::string("certified merge-node direct update threw: ") + e.what())) {
                    return;
                }
                return;
            }

            if (directAccepted) {
                if (event != nullptr) {
                    event->mergeNodeLocalAccepted = true;
                    event->applied = true;
                }
                return;
            }

            if (!directPrimalMerged) {
                if (!attemptMergeNodePhaseRebuildFallback(
                        event != nullptr && !event->mergeNodeLocalFailureSummary.empty()
                            ? event->mergeNodeLocalFailureSummary
                            : "certified merge-node direct update did not apply the primal merge")) {
                    return;
                }
                return;
            }
            if (event != nullptr) {
                event->mergeNodeLocalFallbackToPhaseRebuild = true;
                event->applied = true;
                if (event->mergeNodeLocalFailureSummary.empty()) {
                    event->mergeNodeLocalFailureSummary =
                        "certified merge-node direct update did not satisfy phase-fallback acceptance checks";
                }
            }
            return;
        }

        const auto copyStart = std::chrono::steady_clock::now();
        DynamicComponentTree primalTrial = primalTree;
        DynamicComponentTree dualTrial = dualTree;
        const auto copyEnd = std::chrono::steady_clock::now();
        if (event != nullptr) {
            event->mergeNodeLocalCopyMicros =
                nonIncreasingElapsedMicros(copyStart, copyEnd);
        }

        bool directAccepted = false;
        try {
            DynamicComponentTree *minTrial = dualTrial.isMaxtree() ? &primalTrial : &dualTrial;
            DynamicComponentTree *maxTrial = dualTrial.isMaxtree() ? &dualTrial : &primalTrial;
            DualMinMaxTreeIncrementalFilter<AltitudeType> trialAdjust(
                minTrial,
                maxTrial,
                *adjacency);
            const auto directUpdate = trialAdjust.mergeNodeIntoParentAndUpdate(
                &primalTrial,
                &dualTrial,
                nodeId);
            const bool dualMatchesPrimalAfterMerge =
                directUpdate.primalMergeApplied &&
                directUpdate.localUpdateApplied &&
                dualTrial.reconstructionImage()->isEqual(primalTrial.reconstructionImage());
            directAccepted = dualMatchesPrimalAfterMerge;
            if (event != nullptr) {
                event->mergeNodeLocalApplied =
                    directUpdate.localUpdateApplied && directUpdate.primalMergeApplied;
                event->mergeNodeLocalDualMatchesPrimalAfterMerge =
                    dualMatchesPrimalAfterMerge;
                event->mergeNodeLocalComponentCount = directUpdate.componentCount;
                event->mergeNodeLocalPixelCount = directUpdate.regionSize;
                event->mergeNodeLocalUpdateMicros =
                    directUpdate.collectComponentsUs + directUpdate.localUpdateUs;
                event->mergeNodePrimalMergeMicros = directUpdate.primalMergeUs;
                if (!directAccepted) {
                    event->mergeNodeLocalFailureSummary = directUpdate.message.empty()
                        ? "merge-node local update did not preserve dual/primal reconstruction equality"
                        : directUpdate.message;
                }
            }
        } catch (const std::exception &e) {
            if (event != nullptr) {
                event->mergeNodeLocalFailureSummary = e.what();
            }
        }

        if (directAccepted) {
            primalTree = std::move(primalTrial);
            dualTree = std::move(dualTrial);
            if (event != nullptr) {
                event->mergeNodeLocalAccepted = true;
                event->applied = true;
            }
            return;
        }

        if (!attemptMergeNodePhaseRebuildFallback(
                "merge-node local update did not satisfy phase-fallback acceptance checks")) {
            return;
        }
        return;
    }

    if (dualMaintenance == NonIncreasingDualMaintenance::ExplicitOnly) {
        if (collectDetailedEventDiagnostics) {
            for (const auto &component : components) {
                if (component.empty()) {
                    continue;
                }
                const NodeId nodeCa = findNodeCaForExplicitRegion(
                    dualTree,
                    std::span<const PixelId>(component));
                event->componentNodeCa.push_back(nodeCa);
                event->componentAltitudeCa.push_back(
                    nodeCa == InvalidNode ? -1 : dualTree.getAltitude(nodeCa));
            }
        }

        bool problem1RestrictedCandidateAccepted = true;
        std::string problem1RestrictedRejection;
        for (const auto &component : components) {
            if (component.empty()) {
                continue;
            }
            const auto analysis = analyzeNonIncreasingProblem1ExplicitRegionPreflight<AltitudeType>(
                dualTree,
                *adjacency,
                std::span<const PixelId>(component),
                b,
                true);
            if (analysis.acceptsPublishedProblem() && analysis.affectedOwners() > 1) {
                if (event != nullptr) {
                    ++event->mergeNodeLocalProblem1PreflightRejectedComponents;
                    ++event->mergeNodeLocalProblem1PreflightRejectedByPublishedProblem;
                    if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                        event->mergeNodeLocalProblem1FailureMode = "multiple_nodeCa";
                    }
                }
                problem1RestrictedCandidateAccepted = false;
                problem1RestrictedRejection =
                    "problem1-restricted requires a single nodeCa for each component";
                break;
            }
            if (!analysis.acceptsPublishedProblem() || !analysis.readyForLocalUpdater()) {
                if (event != nullptr) {
                    ++event->mergeNodeLocalProblem1PreflightRejectedComponents;
                    if (analysis.acceptsPublishedProblem()) {
                        ++event->mergeNodeLocalProblem1PreflightRejectedByUpdaterLimitation;
                    } else {
                        ++event->mergeNodeLocalProblem1PreflightRejectedByPublishedProblem;
                    }
                    if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                        event->mergeNodeLocalProblem1FailureMode =
                            classifyProblem1ExplicitRegionFailureMode(analysis);
                    }
                }

                problem1RestrictedCandidateAccepted = false;
                problem1RestrictedRejection =
                    problem1ExplicitRegionPreflightMessage(analysis.preflight.status);
                break;
            }
        }

        if (event != nullptr) {
            event->mergeNodeLocalAttempted = true;
            event->mergeNodeLocalDualMatchesPrimalBeforeMerge =
                dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
        }

        if (!problem1RestrictedCandidateAccepted) {
            if (event != nullptr) {
                event->mergeNodeLocalRejectedComponents =
                    event->mergeNodeLocalProblem1PreflightRejectedComponents;
                if (event->mergeNodeLocalRejectedComponents <= 0) {
                    event->mergeNodeLocalRejectedComponents = 1;
                }
            }
            fallbackToPrimalMergeWithPhaseRebuild(problem1RestrictedRejection.empty()
                                                      ? "problem1-restricted preflight rejected the component set"
                                                      : problem1RestrictedRejection);
            if (event != nullptr && event->mergeNodeLocalFailureSummary.empty()) {
                event->mergeNodeLocalFailureSummary =
                    "problem1-restricted preflight rejected the component set";
                if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                    event->mergeNodeLocalProblem1FailureMode = "contract_rejection";
                }
            }
            return;
        }

        bool trialAccepted = false;
        try {
            DynamicComponentTree primalTrial = primalTree;
            DynamicComponentTree dualTrial = dualTree;
            DualMinMaxTreeIncrementalFilter<AltitudeType> directAdjust(
                &primalTrial,
                &dualTrial,
                *adjacency);
            const auto directUpdate = directAdjust.mergeNodeIntoParentAndUpdate(
                &primalTrial,
                &dualTrial,
                nodeId);
            const bool dualMatchesAfterTrial = directUpdate.primalMergeApplied &&
                directUpdate.localUpdateApplied &&
                dualTrial.reconstructionImage()->isEqual(primalTrial.reconstructionImage());
            if (event != nullptr) {
                event->mergeNodeLocalApplied =
                    directUpdate.localUpdateApplied && directUpdate.primalMergeApplied;
                event->mergeNodeLocalAccepted = dualMatchesAfterTrial;
                event->mergeNodeLocalDualMatchesPrimalAfterMerge = dualMatchesAfterTrial;
                event->mergeNodeLocalComponentCount = directUpdate.componentCount;
                event->mergeNodeLocalPixelCount = directUpdate.regionSize;
                event->mergeNodeLocalUpdateMicros =
                    directUpdate.collectComponentsUs + directUpdate.localUpdateUs;
                event->mergeNodePrimalMergeMicros = directUpdate.primalMergeUs;
                if (!dualMatchesAfterTrial) {
                    event->mergeNodeLocalFailureSummary =
                        "problem1-restricted trial did not preserve dual/primal reconstruction equality";
                    if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                        event->mergeNodeLocalProblem1FailureMode = "trial_postcheck_failed";
                    }
                }
            }
            trialAccepted = directUpdate.inputValid && directUpdate.primalMergeApplied &&
                directUpdate.localUpdateApplied && dualMatchesAfterTrial;
            if (trialAccepted) {
                primalTree = std::move(primalTrial);
                dualTree = std::move(dualTrial);
                if (event != nullptr) {
                    event->applied = true;
                }
                return;
            }
            if (event != nullptr && event->mergeNodeLocalRejectedComponents == 0) {
                ++event->mergeNodeLocalRejectedComponents;
            }
        } catch (const std::exception &e) {
            if (event != nullptr) {
                event->mergeNodeLocalFailureSummary =
                    std::string("problem1-restricted trial update threw: ") + e.what();
                ++event->mergeNodeLocalRejectedComponents;
                if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                    event->mergeNodeLocalProblem1FailureMode = "trial_exception";
                }
            }
        }

        fallbackToPrimalMergeWithPhaseRebuild(
            event != nullptr && !event->mergeNodeLocalFailureSummary.empty()
                ? event->mergeNodeLocalFailureSummary
                : "problem1-restricted trial was rejected; fallback to phase rebuild");
        if (event != nullptr && event->mergeNodeLocalFailureSummary.empty()) {
            event->mergeNodeLocalFailureSummary =
                "problem1-restricted trial was rejected; fallback to phase rebuild";
            if (event->mergeNodeLocalRejectedComponents <= 0) {
                event->mergeNodeLocalRejectedComponents = 1;
            }
            if (event->mergeNodeLocalProblem1FailureMode.empty()) {
                event->mergeNodeLocalProblem1FailureMode = "trial_rejected_fallback";
            }
        }
        return;
    }

    for (const auto &component : components) {
        if (!component.empty()) {
            if (collectDetailedEventDiagnostics) {
                const NodeId nodeCa = findNodeCaForExplicitRegion(dualTree, std::span<const PixelId>(component));
                event->componentNodeCa.push_back(nodeCa);
                event->componentAltitudeCa.push_back(
                    nodeCa == InvalidNode ? -1 : dualTree.getAltitude(nodeCa));
            }

            bool updatedByMaterializedPath = false;
            if (dualMaintenance == NonIncreasingDualMaintenance::MaterializeThenExplicit) {
                const auto anchored = ::debug::region::buildExplicitAnchoredRegionChange<AltitudeType>(
                    dualTree,
                    std::span<const PixelId>(component),
                    b,
                    true,
                    false);
                if (anchored.change.direction == ::debug::region::RegionUpdateDirection::Lower) {
                    DynamicComponentTree &minTree = dualTree.isMaxtree() ? primalTree : dualTree;
                    DynamicComponentTree &maxTree = dualTree.isMaxtree() ? dualTree : primalTree;
                    (void) ::debug::region::RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingPartialOwners(
                        minTree,
                        maxTree,
                        adjacency,
                        anchored.change);
                    updatedByMaterializedPath = true;
                }
            }

            if (!updatedByMaterializedPath) {
                adjust.updateTreeWithExplicitRegion(&dualTree, std::span<const PixelId>(component), b);
            }
        }
    }

    primalTree.mergeNodeIntoParent(nodeId);
    if (dualMaintenance == NonIncreasingDualMaintenance::MaterializeThenExplicit &&
        repairDualMismatchAfterMerge &&
        !dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage())) {
        dualTree.build(primalTree.reconstructionImage()->clone(), dualTree.isMaxtree(), adjacency);
        if (event != nullptr) {
            event->dualRebuiltAfterLocalMismatch = true;
        }
    }

    if (event != nullptr) {
        event->applied = true;
    }
}

inline ImageUInt8Ptr applyNaiveAttributeMergeThreshold(ImageUInt8Ptr input,
                                                       AdjacencyRelationPtr adjacency,
                                                       float threshold,
                                                       NonIncreasingAttributeKind attributeKind) {
    DynamicComponentTree maxTree(input, true, adjacency);
    auto maxAttribute = computeNonIncreasingAttribute(maxTree, attributeKind);
    for (NodeId nodeId : collectNonIncreasingMergeCandidates(maxTree, maxAttribute, threshold)) {
        if (nodeId != maxTree.getRoot() && maxTree.isNode(nodeId) && maxTree.isAlive(nodeId)) {
            maxTree.mergeNodeIntoParent(nodeId);
        }
    }

    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adjacency);
    auto minAttribute = computeNonIncreasingAttribute(minTree, attributeKind);
    for (NodeId nodeId : collectNonIncreasingMergeCandidates(minTree, minAttribute, threshold)) {
        if (nodeId != minTree.getRoot() && minTree.isNode(nodeId) && minTree.isAlive(nodeId)) {
            minTree.mergeNodeIntoParent(nodeId);
        }
    }

    return minTree.reconstructionImage();
}

inline ImageUInt8Ptr applyNaiveGrayVarianceMergeThreshold(ImageUInt8Ptr input,
                                                          AdjacencyRelationPtr adjacency,
                                                          float threshold) {
    return applyNaiveAttributeMergeThreshold(
        input,
        adjacency,
        threshold,
        NonIncreasingAttributeKind::GrayVariance);
}

inline ImageUInt8Ptr applyNaiveAttributeMergeThresholdPrefix(ImageUInt8Ptr input,
                                                             AdjacencyRelationPtr adjacency,
                                                             float threshold,
                                                             int maxAppliedMergeEvents,
                                                             NonIncreasingAttributeKind attributeKind) {
    if (maxAppliedMergeEvents <= 0) {
        return applyNaiveAttributeMergeThreshold(input, adjacency, threshold, attributeKind);
    }

    int appliedMergeEvents = 0;
    DynamicComponentTree maxTree(input, true, adjacency);
    auto maxAttribute = computeNonIncreasingAttribute(maxTree, attributeKind);
    for (NodeId nodeId : collectNonIncreasingMergeCandidates(maxTree, maxAttribute, threshold)) {
        if (nodeId != maxTree.getRoot() && maxTree.isNode(nodeId) && maxTree.isAlive(nodeId)) {
            maxTree.mergeNodeIntoParent(nodeId);
            ++appliedMergeEvents;
            if (appliedMergeEvents >= maxAppliedMergeEvents) {
                return maxTree.reconstructionImage();
            }
        }
    }

    auto current = maxTree.reconstructionImage();
    DynamicComponentTree minTree(current, false, adjacency);
    auto minAttribute = computeNonIncreasingAttribute(minTree, attributeKind);
    for (NodeId nodeId : collectNonIncreasingMergeCandidates(minTree, minAttribute, threshold)) {
        if (nodeId != minTree.getRoot() && minTree.isNode(nodeId) && minTree.isAlive(nodeId)) {
            minTree.mergeNodeIntoParent(nodeId);
            ++appliedMergeEvents;
            if (appliedMergeEvents >= maxAppliedMergeEvents) {
                return minTree.reconstructionImage();
            }
        }
    }

    return minTree.reconstructionImage();
}

inline ImageUInt8Ptr applyNaiveGrayVarianceMergeThresholdPrefix(ImageUInt8Ptr input,
                                                                AdjacencyRelationPtr adjacency,
                                                                float threshold,
                                                                int maxAppliedMergeEvents) {
    return applyNaiveAttributeMergeThresholdPrefix(
        input,
        adjacency,
        threshold,
        maxAppliedMergeEvents,
        NonIncreasingAttributeKind::GrayVariance);
}

inline ImageUInt8Ptr applyNaiveAttributeMergeSequence(
    ImageUInt8Ptr input,
    AdjacencyRelationPtr adjacency,
    const std::vector<float> &thresholds,
    NonIncreasingAttributeKind attributeKind) {
    auto current = input->clone();
    for (float threshold : thresholds) {
        current = applyNaiveAttributeMergeThreshold(
            current,
            adjacency,
            threshold,
            attributeKind);
    }
    return current;
}

inline ImageUInt8Ptr applyNaiveGrayVarianceMergeSequence(ImageUInt8Ptr input,
                                                         AdjacencyRelationPtr adjacency,
                                                         const std::vector<float> &thresholds) {
    return applyNaiveAttributeMergeSequence(
        input,
        adjacency,
        thresholds,
        NonIncreasingAttributeKind::GrayVariance);
}

class NonIncreasingCasfKernel {
private:
    ImageUInt8Ptr inputImage_;
    AdjacencyRelationPtr adjacency_;
    DynamicComponentTree maxTree_;
    DynamicComponentTree minTree_;
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<AltitudeType>> adjust_;
    std::unique_ptr<DynamicAttributeComputer> maxAttributeComputer_;
    std::unique_ptr<DynamicAttributeComputer> minAttributeComputer_;
    std::vector<float> maxAttribute_;
    std::vector<float> minAttribute_;
    NonIncreasingDualMaintenance dualMaintenance_ = NonIncreasingDualMaintenance::PartialF;
    NonIncreasingAttributeKind attributeKind_ = NonIncreasingAttributeKind::GrayVariance;
    NonIncreasingCasfDebugOptionsState debugOptions_;
    bool invariantMismatchFound_ = false;
    bool firstInvariantMismatchOracleFound_ = false;
    bool localOracleCaptured_ = false;
    bool mergeEventLimitReached_ = false;
    bool timeBudgetReached_ = false;
    bool phaseFallbackStopReached_ = false;
    bool phaseReleaseRejectionStopReached_ = false;
    bool localUpdateRejectionStopReached_ = false;
    int phaseDualRebuilds_ = 0;
    int phaseCandidateGuardRebuilds_ = 0;
    int phaseCandidateGuardCanonicalizationChanges_ = 0;
    std::size_t accumulatedNodeStorageCapacityGrowthEvents_ = 0;
    std::size_t accumulatedNodeStorageVectorReallocationEvents_ = 0;
    std::size_t accumulatedNodeStorageReallocatedElements_ = 0;
    std::size_t accumulatedNodeStorageReallocatedBytes_ = 0;
    NonIncreasingCasfDebugState debugState_;
    NonIncreasingDirectComponentScratch componentScratch_;
    NonIncreasingCandidateCollectionScratch candidateScratch_;
    std::vector<NodeId> maxAttributeDirtySeeds_;
    std::vector<NodeId> minAttributeDirtySeeds_;
    std::vector<NodeId> maxAttributeReleasedNodes_;
    std::vector<NodeId> minAttributeReleasedNodes_;

    std::unique_ptr<DynamicAttributeComputer> makeAttributeComputer(
        DynamicComponentTree *tree) const {
        if (attributeKind_ == NonIncreasingAttributeKind::BoxAreaPerimeterRatio) {
            return std::make_unique<DynamicBoundingBoxComputer>(
                tree,
                Attribute::AREA_PER_BOX_PERIMETER);
        }
        return std::make_unique<DynamicGrayLevelVarianceComputer>(tree);
    }

    void computeCurrentAttribute(DynamicAttributeComputer &computer,
                                 std::span<float> buffer) const {
        if (attributeKind_ == NonIncreasingAttributeKind::BoxAreaPerimeterRatio) {
            static_cast<DynamicBoundingBoxComputer &>(computer).compute(buffer);
            return;
        }
        static_cast<DynamicGrayLevelVarianceComputer &>(computer).compute(buffer);
    }

    static std::size_t attributeBufferSizeForTree(const DynamicComponentTree &tree) {
        return static_cast<std::size_t>(tree.getNumInternalNodeSlots());
    }

    void ensureAttributeRuntime() {
        if (!adjust_) {
            adjust_ = std::make_unique<DualMinMaxTreeIncrementalFilter<AltitudeType>>(
                &minTree_,
                &maxTree_,
                *adjacency_);
        }
        if (!maxAttributeComputer_) {
            maxAttributeComputer_ = makeAttributeComputer(&maxTree_);
        }
        if (!minAttributeComputer_) {
            minAttributeComputer_ = makeAttributeComputer(&minTree_);
        }
    }

    void ensureAttributeBuffersAndBind() {
        maxAttribute_.resize(
            attributeBufferSizeForTree(maxTree_),
            0.0f);
        minAttribute_.resize(
            attributeBufferSizeForTree(minTree_),
            0.0f);
    }

    void clearAttributeDirtyState() {
        maxAttributeDirtySeeds_.clear();
        minAttributeDirtySeeds_.clear();
        maxAttributeReleasedNodes_.clear();
        minAttributeReleasedNodes_.clear();
    }

    void refreshAttributeState() {
        ensureAttributeRuntime();
        ensureAttributeBuffersAndBind();
        computeCurrentAttribute(*maxAttributeComputer_, std::span<float>(maxAttribute_));
        computeCurrentAttribute(*minAttributeComputer_, std::span<float>(minAttribute_));
        clearAttributeDirtyState();
    }

    bool canRefreshAttributesLocallyAfterMergeNodePhase() const {
        return attributeKind_ == NonIncreasingAttributeKind::BoxAreaPerimeterRatio &&
               dualMaintenance_ == NonIncreasingDualMaintenance::MergeNodeCleanUnchecked;
    }

    void refreshTreeAttributeDirtyClosure(DynamicComponentTree &tree,
                                          DynamicAttributeComputer &computer,
                                          std::vector<NodeId> &dirtySeeds,
                                          std::vector<NodeId> &releasedNodes,
                                          std::vector<float> &buffer) {
        sortUniqueNonIncreasingNodes(dirtySeeds);
        sortUniqueNonIncreasingNodes(releasedNodes);
        if (dirtySeeds.empty() && releasedNodes.empty()) {
            return;
        }
        (void) ::debug::region::AttributeDirtyRecomputer::recomputeAncestorClosure(
            tree,
            computer,
            std::span<const NodeId>(dirtySeeds.data(), dirtySeeds.size()),
            std::span<const NodeId>(releasedNodes.data(), releasedNodes.size()),
            buffer,
            true);
    }

    void refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase completedPhase) {
        if (!canRefreshAttributesLocallyAfterMergeNodePhase() ||
            !adjust_ ||
            !maxAttributeComputer_ ||
            !minAttributeComputer_) {
            refreshAttributeState();
            return;
        }

        ensureAttributeBuffersAndBind();
        (void) completedPhase;
        refreshTreeAttributeDirtyClosure(maxTree_,
                                         *maxAttributeComputer_,
                                         maxAttributeDirtySeeds_,
                                         maxAttributeReleasedNodes_,
                                         maxAttribute_);
        refreshTreeAttributeDirtyClosure(minTree_,
                                         *minAttributeComputer_,
                                         minAttributeDirtySeeds_,
                                         minAttributeReleasedNodes_,
                                         minAttribute_);
        clearAttributeDirtyState();
    }

    NonIncreasingAttributeDirtyRecorder makeMaxToMinAttributeDirtyRecorder() {
        return NonIncreasingAttributeDirtyRecorder{
            &maxAttributeDirtySeeds_,
            &maxAttributeReleasedNodes_,
            &minAttributeDirtySeeds_,
            &minAttributeReleasedNodes_};
    }

    NonIncreasingAttributeDirtyRecorder makeMinToMaxAttributeDirtyRecorder() {
        return NonIncreasingAttributeDirtyRecorder{
            &minAttributeDirtySeeds_,
            &minAttributeReleasedNodes_,
            &maxAttributeDirtySeeds_,
            &maxAttributeReleasedNodes_};
    }

    void accumulateMinTreeNodeStorageReallocationCounters() {
        accumulatedNodeStorageCapacityGrowthEvents_ +=
            minTree_.getNodeStorageCapacityGrowthEvents();
        accumulatedNodeStorageVectorReallocationEvents_ +=
            minTree_.getNodeStorageVectorReallocationEvents();
        accumulatedNodeStorageReallocatedElements_ +=
            minTree_.getNodeStorageReallocatedElements();
        accumulatedNodeStorageReallocatedBytes_ +=
            minTree_.getNodeStorageReallocatedBytes();
    }

    void accumulateMaxTreeNodeStorageReallocationCounters() {
        accumulatedNodeStorageCapacityGrowthEvents_ +=
            maxTree_.getNodeStorageCapacityGrowthEvents();
        accumulatedNodeStorageVectorReallocationEvents_ +=
            maxTree_.getNodeStorageVectorReallocationEvents();
        accumulatedNodeStorageReallocatedElements_ +=
            maxTree_.getNodeStorageReallocatedElements();
        accumulatedNodeStorageReallocatedBytes_ +=
            maxTree_.getNodeStorageReallocatedBytes();
    }

    void rebuildMinTreeFromMaxTreeImage() {
        accumulateMinTreeNodeStorageReallocationCounters();
        minTree_ = DynamicComponentTree(maxTree_.reconstructionImage(), false, adjacency_);
        minAttributeComputer_.reset();
        adjust_.reset();
    }

    void rebuildMaxTreeFromMinTreeImage() {
        accumulateMaxTreeNodeStorageReallocationCounters();
        maxTree_ = DynamicComponentTree(minTree_.reconstructionImage(), true, adjacency_);
        maxAttributeComputer_.reset();
        adjust_.reset();
    }

    bool reconstructionsMatch() const {
        return maxTree_.reconstructionImage()->isEqual(minTree_.reconstructionImage());
    }

    static bool isRejectedLocalUpdateEvent(const NonIncreasingMergeEvent &event) {
        if (event.applied) {
            return false;
        }
        return (event.partialFAttempted && event.partialFRejected) ||
               event.mergeNodeLocalAttempted ||
               event.strippedLocalAttempted;
    }

    bool shouldCollectDetailedMergeEventDiagnostics() const {
        return debugOptions_.collectMergeEvents ||
               debugOptions_.stopOnFirstInvariantMismatch ||
               debugOptions_.stopOnFirstPhaseRebuildFallback ||
               debugOptions_.collectFirstInvariantMismatchOracle ||
               debugOptions_.localOracleCaptureTarget.enabled ||
               debugOptions_.maxAppliedMergeEventsPerThreshold > 0 ||
               debugOptions_.maxThresholdMicros > 0;
    }

    bool shouldUseOperationalMergeEvent() const {
        if (dualMaintenance_ != NonIncreasingDualMaintenance::MergeNodeCleanUnchecked) {
            return true;
        }
        return debugOptions_.operationalMergeEventInstrumentationEnabled ||
               shouldCollectDetailedMergeEventDiagnostics();
    }

    void rebuildMinTreeIfMaterializedDualMismatch() {
        if (dualMaintenance_ == NonIncreasingDualMaintenance::MaterializeThenExplicit &&
            !reconstructionsMatch()) {
            rebuildMinTreeFromMaxTreeImage();
            ++phaseDualRebuilds_;
        }
    }

    void rebuildMaxTreeIfMaterializedDualMismatch() {
        if (dualMaintenance_ == NonIncreasingDualMaintenance::MaterializeThenExplicit &&
            !reconstructionsMatch()) {
            rebuildMaxTreeFromMinTreeImage();
            ++phaseDualRebuilds_;
        }
    }

    bool recordAndMaybeStopAfterMerge(const NonIncreasingMergeEvent &event,
                                      const LocalMergeOracle *oracleBeforeMerge = nullptr,
                                      bool skipInvariantCheck = false) {
        if (debugOptions_.collectMergeEvents &&
            (event.applied ||
             event.partialFAttempted ||
             event.mergeNodeLocalAttempted ||
             event.strippedLocalAttempted)) {
            debugState_.mergeEvents.push_back(event);
        }

        if (debugOptions_.stopOnFirstPhaseRebuildFallback &&
            !phaseFallbackStopReached_ &&
            event.applied &&
            ((event.strippedLocalFallbackToPhaseRebuild &&
              !event.strippedLocalSkippedByPhaseRebuild) ||
             (event.mergeNodeLocalFallbackToPhaseRebuild &&
              !event.mergeNodeLocalSkippedByPhaseRebuild))) {
            phaseFallbackStopReached_ = true;
            return true;
        }

        if (skipInvariantCheck) {
            return false;
        }

        if (!debugOptions_.stopOnFirstInvariantMismatch || invariantMismatchFound_ || !event.applied) {
            return false;
        }

        if (reconstructionsMatch()) {
            return false;
        }

        debugState_.firstInvariantMismatchEvent = event;
        if (debugOptions_.collectFirstInvariantMismatchOracle && oracleBeforeMerge != nullptr) {
            debugState_.firstInvariantMismatchOracle = *oracleBeforeMerge;
            firstInvariantMismatchOracleFound_ = true;
        }
        invariantMismatchFound_ = true;
        return true;
    }

    bool shouldCaptureLocalOracle(float threshold,
                                  NonIncreasingPhase phase,
                                  NodeId nodeId) const {
        return debugOptions_.localOracleCaptureTarget.enabled &&
               !localOracleCaptured_ &&
               std::abs(threshold - debugOptions_.localOracleCaptureTarget.threshold) < 1.0e-6f &&
               phase == debugOptions_.localOracleCaptureTarget.phase &&
               nodeId == debugOptions_.localOracleCaptureTarget.nodeId;
    }

    bool recordPhaseReleaseRejectionStop(
        const NonIncreasingPhaseOracleComparison &comparison) {
        if (!debugOptions_.stopOnFirstPhaseReleaseRejection ||
            phaseReleaseRejectionStopReached_ ||
            !comparison.releasePredicateEvaluated ||
            !comparison.releaseFocusAvailable ||
            comparison.releasePredicateSatisfied) {
            return false;
        }
        phaseReleaseRejectionStopReached_ = true;
        return true;
    }

    bool checkPhaseCandidateEquivalence(float threshold,
                                        NonIncreasingPhase completedPhase,
                                        bool beforeMaxToMin = false) {
        if (!debugOptions_.collectPhaseOracleComparisons &&
            !debugOptions_.phaseCandidateEquivalenceGuardEnabled &&
            !debugOptions_.stopOnFirstPhaseReleaseRejection) {
            return false;
        }

        auto requiresRepair = [](const NonIncreasingPhaseOracleComparison &comparison) {
            return !comparison.localReconstructionMatchesRebuild ||
                   !comparison.candidateSetMatchesRebuild ||
                   !comparison.candidateScanOrderMatchesRebuild;
        };

        if (completedPhase == NonIncreasingPhase::MaxToMin) {
            DynamicComponentTree rebuiltMinTree(maxTree_.reconstructionImage(), false, adjacency_);
            auto comparison = compareNonIncreasingPhaseAgainstRebuild(
                minTree_,
                rebuiltMinTree,
                threshold,
                completedPhase,
                beforeMaxToMin,
                adjacency_.get(),
                attributeKind_);
            const bool requiresRebuild = requiresRepair(comparison);
            if (debugOptions_.collectPhaseOracleComparisons) {
                debugState_.phaseOracleComparisons.push_back(comparison);
            }
            if (recordPhaseReleaseRejectionStop(comparison)) {
                if (!debugOptions_.collectPhaseOracleComparisons) {
                    debugState_.phaseOracleComparisons.push_back(std::move(comparison));
                }
                return false;
            }
            if (debugOptions_.phaseCandidateEquivalenceGuardEnabled && requiresRebuild) {
                auto [focusNodes, localPixels] =
                    buildPhaseCandidateGuardFocus(minTree_, comparison);
                const int canonicalizationChanges = focusNodes.empty()
                                                        ? 0
                                                        : canonicalizeAfterLocalRecomposition(
                                                              minTree_,
                                                              adjacency_.get(),
                                                              std::span<const NodeId>(focusNodes),
                                                              std::span<const NodeId>{},
                                                              std::span<const PixelId>(localPixels));
                if (canonicalizationChanges > 0) {
                    phaseCandidateGuardCanonicalizationChanges_ += canonicalizationChanges;
                    minAttributeComputer_.reset();
                    adjust_.reset();
                    const auto repairedComparison = compareNonIncreasingPhaseAgainstRebuild(
                        minTree_,
                        rebuiltMinTree,
                        threshold,
                        completedPhase,
                        beforeMaxToMin,
                        adjacency_.get(),
                        attributeKind_);
                    if (!requiresRepair(repairedComparison)) {
                        return true;
                    }
                }
                accumulateMinTreeNodeStorageReallocationCounters();
                minTree_ = std::move(rebuiltMinTree);
                minAttributeComputer_.reset();
                adjust_.reset();
                ++phaseDualRebuilds_;
                ++phaseCandidateGuardRebuilds_;
                return true;
            }
        } else {
            DynamicComponentTree rebuiltMaxTree(minTree_.reconstructionImage(), true, adjacency_);
            auto comparison = compareNonIncreasingPhaseAgainstRebuild(
                maxTree_,
                rebuiltMaxTree,
                threshold,
                completedPhase,
                beforeMaxToMin,
                adjacency_.get(),
                attributeKind_);
            const bool requiresRebuild = requiresRepair(comparison);
            if (debugOptions_.collectPhaseOracleComparisons) {
                debugState_.phaseOracleComparisons.push_back(comparison);
            }
            if (recordPhaseReleaseRejectionStop(comparison)) {
                if (!debugOptions_.collectPhaseOracleComparisons) {
                    debugState_.phaseOracleComparisons.push_back(std::move(comparison));
                }
                return false;
            }
            if (debugOptions_.phaseCandidateEquivalenceGuardEnabled && requiresRebuild) {
                auto [focusNodes, localPixels] =
                    buildPhaseCandidateGuardFocus(maxTree_, comparison);
                const int canonicalizationChanges = focusNodes.empty()
                                                        ? 0
                                                        : canonicalizeAfterLocalRecomposition(
                                                              maxTree_,
                                                              adjacency_.get(),
                                                              std::span<const NodeId>(focusNodes),
                                                              std::span<const NodeId>{},
                                                              std::span<const PixelId>(localPixels));
                if (canonicalizationChanges > 0) {
                    phaseCandidateGuardCanonicalizationChanges_ += canonicalizationChanges;
                    maxAttributeComputer_.reset();
                    adjust_.reset();
                    const auto repairedComparison = compareNonIncreasingPhaseAgainstRebuild(
                        maxTree_,
                        rebuiltMaxTree,
                        threshold,
                        completedPhase,
                        beforeMaxToMin,
                        adjacency_.get(),
                        attributeKind_);
                    if (!requiresRepair(repairedComparison)) {
                        return true;
                    }
                }
                accumulateMaxTreeNodeStorageReallocationCounters();
                maxTree_ = std::move(rebuiltMaxTree);
                maxAttributeComputer_.reset();
                adjust_.reset();
                ++phaseDualRebuilds_;
                ++phaseCandidateGuardRebuilds_;
                return true;
            }
        }
        return false;
    }

public:
    explicit NonIncreasingCasfKernel(ImageUInt8Ptr image,
                                     double radiusAdj,
                                     NonIncreasingDualMaintenance dualMaintenance = NonIncreasingDualMaintenance::PartialF,
                                     NonIncreasingAttributeKind attributeKind = NonIncreasingAttributeKind::GrayVariance)
        : inputImage_(image ? image->clone() : nullptr),
          adjacency_(image ? std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radiusAdj) : nullptr),
          maxTree_(inputImage_, true, adjacency_),
          minTree_(inputImage_, false, adjacency_),
          dualMaintenance_(dualMaintenance),
          attributeKind_(attributeKind) {
        if (image == nullptr) {
            throw std::runtime_error("NonIncreasingCasfKernel requires a valid image");
        }
        refreshAttributeState();
    }

    void applyThreshold(float threshold) {
        debugState_.mergeEvents.clear();
        debugState_.phaseOracleComparisons.clear();
        invariantMismatchFound_ = false;
        firstInvariantMismatchOracleFound_ = false;
        phaseDualRebuilds_ = 0;
        phaseCandidateGuardRebuilds_ = 0;
        phaseCandidateGuardCanonicalizationChanges_ = 0;
        mergeEventLimitReached_ = false;
        timeBudgetReached_ = false;
        phaseFallbackStopReached_ = false;
        phaseReleaseRejectionStopReached_ = false;
        localUpdateRejectionStopReached_ = false;
        const bool usePhaseRebuildFallback =
            usesNonIncreasingPhaseRebuildFallback(dualMaintenance_);
        bool maxPhaseRequiresDualRebuild = false;
        bool minPhaseRequiresDualRebuild = false;
        int appliedMergeEvents = 0;
        const auto thresholdStart = std::chrono::steady_clock::now();
        auto reachedMergeEventLimit = [&]() {
            if (debugOptions_.maxAppliedMergeEventsPerThreshold <= 0) {
                return false;
            }
            ++appliedMergeEvents;
            if (appliedMergeEvents >= debugOptions_.maxAppliedMergeEventsPerThreshold) {
                mergeEventLimitReached_ = true;
                return true;
            }
            return false;
        };
        auto reachedTimeBudget = [&]() {
            if (debugOptions_.maxThresholdMicros <= 0) {
                return false;
            }
            const auto elapsed = nonIncreasingElapsedMicros(
                thresholdStart,
                std::chrono::steady_clock::now());
            if (elapsed >= debugOptions_.maxThresholdMicros) {
                timeBudgetReached_ = true;
                return true;
            }
            return false;
        };

        if (checkPhaseCandidateEquivalence(threshold, NonIncreasingPhase::MinToMax, true)) {
            refreshAttributeState();
        }
        if (phaseReleaseRejectionStopReached_) {
            refreshAttributeState();
            return;
        }
        if (reachedTimeBudget()) {
            return;
        }
        auto maxCandidates = collectNonIncreasingMergeCandidates(
            maxTree_,
            maxAttribute_,
            threshold,
            &candidateScratch_);
        auto rebuildMinAfterDirtyMaxPhase = [&]() {
            if (maxPhaseRequiresDualRebuild) {
                rebuildMinTreeFromMaxTreeImage();
                ++phaseDualRebuilds_;
                maxPhaseRequiresDualRebuild = false;
            }
        };
        for (NodeId nodeId : maxCandidates) {
            std::optional<NonIncreasingMergeEvent> eventStorage;
            NonIncreasingMergeEvent *eventPtr = nullptr;
            std::optional<LocalMergeOracle> oracleBeforeMerge;
            const bool validCandidate =
                nodeId != InvalidNode &&
                maxTree_.isNode(nodeId) &&
                maxTree_.isAlive(nodeId) &&
                static_cast<std::size_t>(nodeId) < maxAttribute_.size();
            if (validCandidate) {
                const bool captureTarget = shouldCaptureLocalOracle(threshold, NonIncreasingPhase::MaxToMin, nodeId);
                if (!maxPhaseRequiresDualRebuild &&
                    (captureTarget ||
                     (debugOptions_.collectFirstInvariantMismatchOracle && debugOptions_.stopOnFirstInvariantMismatch))) {
                    oracleBeforeMerge = buildLocalMergeOracle(maxTree_, minTree_, adjacency_, nodeId);
                }
                if (captureTarget && oracleBeforeMerge) {
                    debugState_.capturedLocalOracle = *oracleBeforeMerge;
                    localOracleCaptured_ = true;
                }
            }
            if (validCandidate && shouldUseOperationalMergeEvent()) {
                eventStorage.emplace();
                NonIncreasingMergeEvent &event = *eventStorage;
                event.collectDetailedDiagnostics =
                    shouldCollectDetailedMergeEventDiagnostics();
                event.threshold = threshold;
                event.phase = NonIncreasingPhase::MaxToMin;
                event.attributeValue = maxAttribute_[static_cast<std::size_t>(nodeId)];
                eventPtr = &event;
            }

            auto maxToMinAttributeDirtyRecorder = makeMaxToMinAttributeDirtyRecorder();
            mergeNodeIntoParentWithExplicitDualUpdate(
                maxTree_,
                minTree_,
                *adjust_,
                adjacency_,
                dualMaintenance_,
                nodeId,
                eventPtr,
                debugOptions_.stopOnFirstInvariantMismatch,
                debugOptions_.leafExplicitFastPathEnabled,
                debugOptions_.strippedLocalUpdateOptions,
                debugOptions_.strippedCommitMode,
                usePhaseRebuildFallback && maxPhaseRequiresDualRebuild,
                debugOptions_.mcpcDirectPreflightDiagnosticsEnabled,
                debugOptions_.mcpcDirectV1StrictPrecheckEnabled,
                debugOptions_.mcpcNoCopyDirectCommitEnabled,
                debugOptions_.mcpcValidatedDirectCommitEnabled,
                false,
                &componentScratch_,
                canRefreshAttributesLocallyAfterMergeNodePhase()
                    ? &maxToMinAttributeDirtyRecorder
                    : nullptr);
            if (eventPtr != nullptr &&
                (eventPtr->strippedLocalFallbackToPhaseRebuild ||
                 eventPtr->mergeNodeLocalFallbackToPhaseRebuild)) {
                maxPhaseRequiresDualRebuild = true;
            }
            const bool rejectedLocalUpdate =
                eventPtr != nullptr && isRejectedLocalUpdateEvent(*eventPtr);
            if (eventPtr != nullptr &&
                recordAndMaybeStopAfterMerge(
                    *eventPtr,
                    oracleBeforeMerge ? &*oracleBeforeMerge : nullptr,
                    maxPhaseRequiresDualRebuild)) {
                if (phaseFallbackStopReached_) {
                    rebuildMinAfterDirtyMaxPhase();
                    refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MaxToMin);
                }
                return;
            }
            if (rejectedLocalUpdate) {
                localUpdateRejectionStopReached_ = true;
                refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MaxToMin);
                return;
            }
            if (eventPtr != nullptr && eventPtr->applied && reachedMergeEventLimit()) {
                rebuildMinAfterDirtyMaxPhase();
                refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MaxToMin);
                return;
            }
            if (eventPtr != nullptr && eventPtr->applied && reachedTimeBudget()) {
                rebuildMinAfterDirtyMaxPhase();
                refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MaxToMin);
                return;
            }
        }

        if (dualMaintenance_ == NonIncreasingDualMaintenance::RebuildAfterPhase) {
            rebuildMinTreeFromMaxTreeImage();
            ++phaseDualRebuilds_;
        } else if (maxPhaseRequiresDualRebuild) {
            rebuildMinAfterDirtyMaxPhase();
        } else {
            rebuildMinTreeIfMaterializedDualMismatch();
        }
        if (checkPhaseCandidateEquivalence(threshold, NonIncreasingPhase::MaxToMin)) {
            // The local dual can be image-equivalent but expose different
            // non-increasing candidates. Rebuild before collecting the next
            // phase candidates.
        }
        refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MaxToMin);
        if (phaseReleaseRejectionStopReached_) {
            return;
        }
        if (reachedTimeBudget()) {
            return;
        }

        auto minCandidates = collectNonIncreasingMergeCandidates(
            minTree_,
            minAttribute_,
            threshold,
            &candidateScratch_);
        auto rebuildMaxAfterDirtyMinPhase = [&]() {
            if (minPhaseRequiresDualRebuild) {
                rebuildMaxTreeFromMinTreeImage();
                ++phaseDualRebuilds_;
                minPhaseRequiresDualRebuild = false;
            }
        };
        for (NodeId nodeId : minCandidates) {
            std::optional<NonIncreasingMergeEvent> eventStorage;
            NonIncreasingMergeEvent *eventPtr = nullptr;
            std::optional<LocalMergeOracle> oracleBeforeMerge;
            const bool validCandidate =
                nodeId != InvalidNode &&
                minTree_.isNode(nodeId) &&
                minTree_.isAlive(nodeId) &&
                static_cast<std::size_t>(nodeId) < minAttribute_.size();
            if (validCandidate) {
                const bool captureTarget = shouldCaptureLocalOracle(threshold, NonIncreasingPhase::MinToMax, nodeId);
                if (!minPhaseRequiresDualRebuild &&
                    (captureTarget ||
                     (debugOptions_.collectFirstInvariantMismatchOracle && debugOptions_.stopOnFirstInvariantMismatch))) {
                    oracleBeforeMerge = buildLocalMergeOracle(minTree_, maxTree_, adjacency_, nodeId);
                }
                if (captureTarget && oracleBeforeMerge) {
                    debugState_.capturedLocalOracle = *oracleBeforeMerge;
                    localOracleCaptured_ = true;
                }
            }
            if (validCandidate && shouldUseOperationalMergeEvent()) {
                eventStorage.emplace();
                NonIncreasingMergeEvent &event = *eventStorage;
                event.collectDetailedDiagnostics =
                    shouldCollectDetailedMergeEventDiagnostics();
                event.threshold = threshold;
                event.phase = NonIncreasingPhase::MinToMax;
                event.attributeValue = minAttribute_[static_cast<std::size_t>(nodeId)];
                eventPtr = &event;
            }

            auto minToMaxAttributeDirtyRecorder = makeMinToMaxAttributeDirtyRecorder();
            mergeNodeIntoParentWithExplicitDualUpdate(
                minTree_,
                maxTree_,
                *adjust_,
                adjacency_,
                dualMaintenance_,
                nodeId,
                eventPtr,
                debugOptions_.stopOnFirstInvariantMismatch,
                debugOptions_.leafExplicitFastPathEnabled,
                debugOptions_.strippedLocalUpdateOptions,
                debugOptions_.strippedCommitMode,
                usePhaseRebuildFallback && minPhaseRequiresDualRebuild,
                debugOptions_.mcpcDirectPreflightDiagnosticsEnabled,
                debugOptions_.mcpcDirectV1StrictPrecheckEnabled,
                debugOptions_.mcpcNoCopyDirectCommitEnabled,
                debugOptions_.mcpcValidatedDirectCommitEnabled,
                false,
                &componentScratch_,
                canRefreshAttributesLocallyAfterMergeNodePhase()
                    ? &minToMaxAttributeDirtyRecorder
                    : nullptr);
            if (eventPtr != nullptr &&
                (eventPtr->strippedLocalFallbackToPhaseRebuild ||
                 eventPtr->mergeNodeLocalFallbackToPhaseRebuild)) {
                minPhaseRequiresDualRebuild = true;
            }
            const bool rejectedLocalUpdate =
                eventPtr != nullptr && isRejectedLocalUpdateEvent(*eventPtr);
            if (eventPtr != nullptr &&
                recordAndMaybeStopAfterMerge(
                    *eventPtr,
                    oracleBeforeMerge ? &*oracleBeforeMerge : nullptr,
                    minPhaseRequiresDualRebuild)) {
                if (phaseFallbackStopReached_) {
                    rebuildMaxAfterDirtyMinPhase();
                    refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MinToMax);
                }
                return;
            }
            if (rejectedLocalUpdate) {
                localUpdateRejectionStopReached_ = true;
                refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MinToMax);
                return;
            }
            if (eventPtr != nullptr && eventPtr->applied && reachedMergeEventLimit()) {
                rebuildMaxAfterDirtyMinPhase();
                refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MinToMax);
                return;
            }
            if (eventPtr != nullptr && eventPtr->applied && reachedTimeBudget()) {
                rebuildMaxAfterDirtyMinPhase();
                refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MinToMax);
                return;
            }
        }

        if (dualMaintenance_ == NonIncreasingDualMaintenance::RebuildAfterPhase) {
            rebuildMaxTreeFromMinTreeImage();
            ++phaseDualRebuilds_;
        } else if (minPhaseRequiresDualRebuild) {
            rebuildMaxAfterDirtyMinPhase();
        } else {
            rebuildMaxTreeIfMaterializedDualMismatch();
        }
        if (checkPhaseCandidateEquivalence(threshold, NonIncreasingPhase::MinToMax)) {
            // Keep the max-tree canonical for the next threshold.
        }
        refreshAttributeStateAfterLocalMergeNodePhase(NonIncreasingPhase::MinToMax);
    }

    ImageUInt8Ptr filter(const std::vector<float> &thresholds) {
        for (float threshold : thresholds) {
            applyThreshold(threshold);
        }
        return minTree_.reconstructionImage();
    }

    ImageUInt8Ptr reconstructImage() const {
        return minTree_.reconstructionImage();
    }

    void setMergeEventCollectionEnabled(bool enabled) {
        debugOptions_.collectMergeEvents = enabled;
        if (!enabled) {
            debugState_.mergeEvents.clear();
        }
    }

    void setPhaseOracleComparisonEnabled(bool enabled) {
        debugOptions_.collectPhaseOracleComparisons = enabled;
        if (!enabled) {
            debugState_.phaseOracleComparisons.clear();
        }
    }

    void setPhaseCandidateEquivalenceGuardEnabled(bool enabled) {
        debugOptions_.phaseCandidateEquivalenceGuardEnabled = enabled;
        phaseCandidateGuardRebuilds_ = 0;
        phaseCandidateGuardCanonicalizationChanges_ = 0;
    }

    void setStopOnFirstInvariantMismatchEnabled(bool enabled) {
        debugOptions_.stopOnFirstInvariantMismatch = enabled;
        invariantMismatchFound_ = false;
        firstInvariantMismatchOracleFound_ = false;
    }

    void setStopOnFirstPhaseRebuildFallbackEnabled(bool enabled) {
        debugOptions_.stopOnFirstPhaseRebuildFallback = enabled;
        phaseFallbackStopReached_ = false;
    }

    void setStopOnFirstPhaseReleaseRejectionEnabled(bool enabled) {
        debugOptions_.stopOnFirstPhaseReleaseRejection = enabled;
        phaseReleaseRejectionStopReached_ = false;
    }

    void setFirstInvariantMismatchOracleCollectionEnabled(bool enabled) {
        debugOptions_.collectFirstInvariantMismatchOracle = enabled;
        if (!enabled) {
            firstInvariantMismatchOracleFound_ = false;
        }
    }

    void setLocalOracleCaptureTarget(float threshold,
                                     NonIncreasingPhase phase,
                                     NodeId nodeId) {
        debugOptions_.localOracleCaptureTarget.threshold = threshold;
        debugOptions_.localOracleCaptureTarget.phase = phase;
        debugOptions_.localOracleCaptureTarget.nodeId = nodeId;
        debugOptions_.localOracleCaptureTarget.enabled = true;
        localOracleCaptured_ = false;
    }

    void clearLocalOracleCaptureTarget() {
        debugOptions_.localOracleCaptureTarget = NonIncreasingLocalOracleCaptureTarget{};
        localOracleCaptured_ = false;
    }

    void setMaxAppliedMergeEventsPerThreshold(int maxEvents) {
        if (maxEvents < 0) {
            throw std::runtime_error("max applied merge events per threshold must be non-negative");
        }
        debugOptions_.maxAppliedMergeEventsPerThreshold = maxEvents;
    }

    void setMaxThresholdMicros(long long maxMicros) {
        if (maxMicros < 0) {
            throw std::runtime_error("max threshold time budget must be non-negative");
        }
        debugOptions_.maxThresholdMicros = maxMicros;
    }

    void setLeafExplicitFastPathEnabled(bool enabled) {
        debugOptions_.leafExplicitFastPathEnabled = enabled;
    }

    void setMcpcDirectPreflightDiagnosticsEnabled(bool enabled) {
        debugOptions_.mcpcDirectPreflightDiagnosticsEnabled = enabled;
    }

    void setMcpcDirectV1StrictPrecheckEnabled(bool enabled) {
        debugOptions_.mcpcDirectV1StrictPrecheckEnabled = enabled;
    }

    void setMcpcNoCopyDirectCommitEnabled(bool enabled) {
        debugOptions_.mcpcNoCopyDirectCommitEnabled = enabled;
    }

    void setMcpcValidatedDirectCommitEnabled(bool enabled) {
        debugOptions_.mcpcValidatedDirectCommitEnabled = enabled;
    }

    void setOperationalMergeEventInstrumentationEnabled(bool enabled) {
        debugOptions_.operationalMergeEventInstrumentationEnabled = enabled;
    }

    void setStrippedLocalUpdateOptions(const NonIncreasingStrippedLocalUpdateOptions &options) {
        debugOptions_.strippedLocalUpdateOptions = options;
    }

    void setStrippedCommitMode(NonIncreasingStrippedCommitMode mode) {
        debugOptions_.strippedCommitMode = mode;
    }

    bool lastMergeEventLimitReached() const {
        return mergeEventLimitReached_;
    }

    bool lastTimeBudgetReached() const {
        return timeBudgetReached_;
    }

    bool lastPhaseFallbackStopReached() const {
        return phaseFallbackStopReached_;
    }

    bool lastPhaseReleaseRejectionStopReached() const {
        return phaseReleaseRejectionStopReached_;
    }

    bool lastLocalUpdateRejectionStopReached() const {
        return localUpdateRejectionStopReached_;
    }

    const std::vector<NonIncreasingMergeEvent> &lastMergeEvents() const {
        return debugState_.mergeEvents;
    }

    const std::vector<NonIncreasingPhaseOracleComparison> &lastPhaseOracleComparisons() const {
        return debugState_.phaseOracleComparisons;
    }

    int lastPhaseDualRebuilds() const {
        return phaseDualRebuilds_;
    }

    int lastPhaseCandidateGuardRebuilds() const {
        return phaseCandidateGuardRebuilds_;
    }

    int lastPhaseCandidateGuardCanonicalizationChanges() const {
        return phaseCandidateGuardCanonicalizationChanges_;
    }

    bool hasFirstInvariantMismatch() const {
        return invariantMismatchFound_;
    }

    const NonIncreasingMergeEvent &firstInvariantMismatchEvent() const {
        return debugState_.firstInvariantMismatchEvent;
    }

    bool hasFirstInvariantMismatchOracle() const {
        return firstInvariantMismatchOracleFound_;
    }

    const LocalMergeOracle &firstInvariantMismatchOracle() const {
        return debugState_.firstInvariantMismatchOracle;
    }

    bool hasCapturedLocalOracle() const {
        return localOracleCaptured_;
    }

    const LocalMergeOracle &capturedLocalOracle() const {
        return debugState_.capturedLocalOracle;
    }

    std::size_t nodeStorageCapacityGrowthEvents() const {
        return accumulatedNodeStorageCapacityGrowthEvents_ +
               maxTree_.getNodeStorageCapacityGrowthEvents() +
               minTree_.getNodeStorageCapacityGrowthEvents();
    }

    std::size_t nodeStorageVectorReallocationEvents() const {
        return accumulatedNodeStorageVectorReallocationEvents_ +
               maxTree_.getNodeStorageVectorReallocationEvents() +
               minTree_.getNodeStorageVectorReallocationEvents();
    }

    std::size_t nodeStorageReallocatedElements() const {
        return accumulatedNodeStorageReallocatedElements_ +
               maxTree_.getNodeStorageReallocatedElements() +
               minTree_.getNodeStorageReallocatedElements();
    }

    std::size_t nodeStorageReallocatedBytes() const {
        return accumulatedNodeStorageReallocatedBytes_ +
               maxTree_.getNodeStorageReallocatedBytes() +
               minTree_.getNodeStorageReallocatedBytes();
    }

    DynamicComponentTree &maxTree() { return maxTree_; }
    DynamicComponentTree &minTree() { return minTree_; }
    const DynamicComponentTree &maxTree() const { return maxTree_; }
    const DynamicComponentTree &minTree() const { return minTree_; }
    const std::vector<float> &maxAttribute() const { return maxAttribute_; }
    const std::vector<float> &minAttribute() const { return minAttribute_; }
};

} // namespace morphotree::details::non_increasing
