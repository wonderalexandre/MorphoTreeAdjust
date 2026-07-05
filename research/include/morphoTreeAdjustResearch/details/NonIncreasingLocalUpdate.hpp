#pragma once

#include <algorithm>
#include <chrono>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "../AdjacencyRelation.hpp"
#include "../Common.hpp"
#include "../DynamicComponentTree.hpp"
#include "NonIncreasingLocalOracle.hpp"
#include "NonIncreasingLocalPartition.hpp"
#include "NonIncreasingLocalRecomposer.hpp"
#include "NonIncreasingStrippedLocalRecomposition.hpp"

namespace morphotree::details::non_increasing {

inline long long localUpdateElapsedMicros(std::chrono::steady_clock::time_point start,
                                          std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

struct NonIncreasingOracleGuidedLocalUpdateResult {
    bool applied = false;
    bool skippedNoImageChange = false;
    bool oracleBuilt = false;
    bool preDualMatchesPrimal = false;
    bool oracleDualMatchesUpdatedImage = false;
    bool localRecompositionApplied = false;
    bool localRecompositionMatchesTarget = false;
    bool localStructureMatchesOracle = false;
    bool dualMatchesPrimalAfterMerge = false;
    int affectedLevelCount = 0;
    int recompositionNodeCount = 0;
    int createdNodeCount = 0;
    int releasedSourceOwnerCount = 0;
    int postEmptyNodeReleasesRequired = 0;
    int postSameLevelParentChildCollapsesRequired = 0;
    int postSameLevelSiblingMergesRequired = 0;
    int postComparableSiblingMovesRequired = 0;
    int postCanonicalizationChangesRequired = 0;
    int oracleLocalNodeCount = 0;
    int postLocalNodeCount = 0;
    std::string localStructureMismatchSummary;
};

struct NonIncreasingStrippedLocalUpdateResult {
    bool applied = false;
    bool skippedNoImageChange = false;
    bool inputBuilt = false;
    bool localPlanOracleAuditChecked = false;
    bool localPlanOracleFree = false;
    bool materializationInputOracleAuditChecked = false;
    bool materializationInputOracleFree = false;
    bool localRecompositionApplied = false;
    bool localRecompositionTargetChecked = false;
    bool localRecompositionMatchesTarget = false;
    bool localStructureAgainstRebuildChecked = false;
    bool localStructureMatchesRebuild = false;
    bool finalReconstructionChecked = false;
    bool dualMatchesPrimalAfterMerge = false;
    int affectedLevelCount = 0;
    int localPixelCount = 0;
    int preLocalNodeCount = 0;
    int recompositionNodeCount = 0;
    int createdNodeCount = 0;
    int createdPlannedNodeCount = 0;
    std::vector<int> createdNodeLevels;
    std::vector<int> createdNodePlanIndices;
    int createdConnectorNodeCount = 0;
    int createdHiddenExternalParentConnectorNodeCount = 0;
    int createdOracleBridgeNodeCount = 0;
    int createdStableChildBridgeNodeCount = 0;
    int createdHiddenSupportConnectorNodeCount = 0;
    int createdNonOracleHiddenSupportConnectorNodeCount = 0;
    int createdNonOracleParentLinkBridgeNodeCount = 0;
    int createdNonOracleChildLinkBridgeNodeCount = 0;
    int nonOracleHiddenSupportCandidateCount = 0;
    int nonOracleHiddenSupportAlreadyRepresentedCount = 0;
    int nonOracleHiddenSupportMaterializableCount = 0;
    int nonOracleParentLinkTouchesPlanCount = 0;
    int nonOracleParentLinkBothRepresentedCount = 0;
    int nonOracleParentLinkBlockedByModifiedSnapshotCount = 0;
    int nonOracleParentLinkInvalidLevelCount = 0;
    int nonOracleParentLinkAlreadyDirectCount = 0;
    int nonOracleParentLinkAlreadyIntermediateCount = 0;
    int nonOracleParentLinkMovedCount = 0;
    int nonOracleParentLinkBridgeAttemptCount = 0;
    int nonOracleParentLinkBridgeMaterializableCount = 0;
    int nonOracleChildLinkBridgeAttemptCount = 0;
    int nonOracleChildLinkBridgeInvalidLevelCount = 0;
    int nonOracleChildLinkBridgeMaterializableCount = 0;
    int movedDirectProperPartCount = 0;
    int movedStableChildNodeCount = 0;
    int emptiedSourceOwnerCount = 0;
    int releasedSourceOwnerCount = 0;
    int postEmptyNodeReleasesRequired = 0;
    int postSameLevelParentChildCollapsesRequired = 0;
    int postSameLevelSiblingMergesRequired = 0;
    int postComparableSiblingMovesRequired = 0;
    int postCanonicalizationChangesRequired = 0;
    int postLocalNodeCount = 0;
    std::vector<LocalEmptyCreatedNodeDiagnostic> emptyCreatedNodesAfterCanonicalization;
    std::string localStructureMismatchSummary;
    long long inputExtractionMicros = 0;
    long long partitionPlanMicros = 0;
    long long recompositionPlanMicros = 0;
    long long planAuditMicros = 0;
    long long materializationContextMicros = 0;
    long long materializationAuditMicros = 0;
    long long materializationMicros = 0;
    long long postLocalSnapshotMicros = 0;
    long long primalMergeReleaseMicros = 0;
    long long finalReconstructionCheckMicros = 0;
    long long materializationValidationMicros = 0;
    long long materializationCreateNodesMicros = 0;
    long long materializationAttachNodesMicros = 0;
    long long materializationMoveProperPartsMicros = 0;
    long long materializationChildBridgeMicros = 0;
    long long materializationStructuralRepairMicros = 0;
    long long materializationReleaseEmptySourceOwnersMicros = 0;
    long long materializationRootPromotionMicros = 0;
    long long materializationCanonicalizationMicros = 0;
    long long materializationTargetReconstructionCheckMicros = 0;
};

struct NonIncreasingStrippedLocalUpdateOptions {
    bool checkInputPreReconstruction = true;
    bool buildFullUpdatedImage = true;
    bool auditPlanOracleDependencies = true;
    bool auditMaterializationOracleDependencies = true;
    bool checkPreMaterializationReconstruction = true;
    bool checkMaterializationTargetReconstruction = true;
    bool compareLocalStructureAgainstRebuild = false;
    bool checkNoEmptyLiveNodesAfterCanonicalization = false;
    bool checkFinalDualPrimalReconstruction = true;
    bool releaseEmptyNodesAfterPrimalMerge = true;
    bool releaseEmptySourceOwnersWithResidualChildren = false;
    bool testThrowAfterMaterializationBeforePrimalMerge = false;
};

inline void recordStrippedCreatedRecompositionNodes(
    NonIncreasingStrippedLocalUpdateResult &result,
    const LocalRecompositionMaterializationResult &materialization) {
    result.createdNodeCount = materialization.numCreatedNodes();
    result.createdPlannedNodeCount = 0;
    result.movedDirectProperPartCount = 0;
    result.movedStableChildNodeCount = 0;
    result.createdNodeLevels.clear();
    result.createdNodePlanIndices.clear();
    result.createdConnectorNodeCount = 0;
    result.createdHiddenExternalParentConnectorNodeCount = 0;
    result.createdOracleBridgeNodeCount = 0;
    result.createdStableChildBridgeNodeCount = 0;
    result.createdHiddenSupportConnectorNodeCount = 0;
    result.createdNonOracleHiddenSupportConnectorNodeCount = 0;
    result.createdNonOracleParentLinkBridgeNodeCount = 0;
    result.createdNonOracleChildLinkBridgeNodeCount = 0;
    result.createdNodeLevels.reserve(materialization.createdNodes.size());
    result.createdNodePlanIndices.reserve(materialization.createdNodes.size());
    for (const auto &created : materialization.createdNodes) {
        result.createdNodeLevels.push_back(static_cast<int>(created.level));
        result.createdNodePlanIndices.push_back(created.planIndex);
        if (created.planIndex >= 0) {
            ++result.createdPlannedNodeCount;
            result.movedDirectProperPartCount +=
                static_cast<int>(created.movedDirectProperParts.size());
            result.movedStableChildNodeCount +=
                static_cast<int>(created.movedChildNodes.size());
            continue;
        }

        ++result.createdConnectorNodeCount;
        switch (created.kind) {
        case LocalCreatedRecompositionNodeKind::HiddenExternalParentConnector:
            ++result.createdHiddenExternalParentConnectorNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::OracleBridge:
            ++result.createdOracleBridgeNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::StableChildBridge:
            ++result.createdStableChildBridgeNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::HiddenSupportConnector:
            ++result.createdHiddenSupportConnectorNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::NonOracleHiddenSupportConnector:
            ++result.createdNonOracleHiddenSupportConnectorNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::NonOracleParentLinkBridge:
            ++result.createdNonOracleParentLinkBridgeNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::NonOracleChildLinkBridge:
            ++result.createdNonOracleChildLinkBridgeNodeCount;
            break;
        case LocalCreatedRecompositionNodeKind::Planned:
            break;
        }
    }
}

inline void recordStrippedMaterializationGuardDiagnostics(
    NonIncreasingStrippedLocalUpdateResult &result,
    const LocalRecompositionMaterializationResult &materialization) {
    result.nonOracleHiddenSupportCandidateCount =
        materialization.nonOracleHiddenSupportCandidateCount;
    result.nonOracleHiddenSupportAlreadyRepresentedCount =
        materialization.nonOracleHiddenSupportAlreadyRepresentedCount;
    result.nonOracleHiddenSupportMaterializableCount =
        materialization.nonOracleHiddenSupportMaterializableCount;
    result.nonOracleParentLinkTouchesPlanCount =
        materialization.nonOracleParentLinkTouchesPlanCount;
    result.nonOracleParentLinkBothRepresentedCount =
        materialization.nonOracleParentLinkBothRepresentedCount;
    result.nonOracleParentLinkBlockedByModifiedSnapshotCount =
        materialization.nonOracleParentLinkBlockedByModifiedSnapshotCount;
    result.nonOracleParentLinkInvalidLevelCount =
        materialization.nonOracleParentLinkInvalidLevelCount;
    result.nonOracleParentLinkAlreadyDirectCount =
        materialization.nonOracleParentLinkAlreadyDirectCount;
    result.nonOracleParentLinkAlreadyIntermediateCount =
        materialization.nonOracleParentLinkAlreadyIntermediateCount;
    result.nonOracleParentLinkMovedCount =
        materialization.nonOracleParentLinkMovedCount;
    result.nonOracleParentLinkBridgeAttemptCount =
        materialization.nonOracleParentLinkBridgeAttemptCount;
    result.nonOracleParentLinkBridgeMaterializableCount =
        materialization.nonOracleParentLinkBridgeMaterializableCount;
    result.nonOracleChildLinkBridgeAttemptCount =
        materialization.nonOracleChildLinkBridgeAttemptCount;
    result.nonOracleChildLinkBridgeInvalidLevelCount =
        materialization.nonOracleChildLinkBridgeInvalidLevelCount;
    result.nonOracleChildLinkBridgeMaterializableCount =
        materialization.nonOracleChildLinkBridgeMaterializableCount;
}

struct LocalDualNodeShape {
    int altitude = 0;
    std::vector<PixelId> directLocalProperParts;
    std::vector<PixelId> supportLocalPixels;
    bool parentInSnapshot = false;
    int parentAltitude = 0;
    std::vector<PixelId> parentSupportLocalPixels;
    std::vector<std::pair<int, std::vector<PixelId>>> children;
};

inline bool operator<(const LocalDualNodeShape &lhs, const LocalDualNodeShape &rhs) {
    return std::tie(lhs.altitude,
                    lhs.directLocalProperParts,
                    lhs.supportLocalPixels,
                    lhs.parentInSnapshot,
                    lhs.parentAltitude,
                    lhs.parentSupportLocalPixels,
                    lhs.children) <
           std::tie(rhs.altitude,
                    rhs.directLocalProperParts,
                    rhs.supportLocalPixels,
                    rhs.parentInSnapshot,
                    rhs.parentAltitude,
                    rhs.parentSupportLocalPixels,
                    rhs.children);
}

inline bool operator==(const LocalDualNodeShape &lhs, const LocalDualNodeShape &rhs) {
    return std::tie(lhs.altitude,
                    lhs.directLocalProperParts,
                    lhs.supportLocalPixels,
                    lhs.parentInSnapshot,
                    lhs.parentAltitude,
                    lhs.parentSupportLocalPixels,
                    lhs.children) ==
           std::tie(rhs.altitude,
                    rhs.directLocalProperParts,
                    rhs.supportLocalPixels,
                    rhs.parentInSnapshot,
                    rhs.parentAltitude,
                    rhs.parentSupportLocalPixels,
                    rhs.children);
}

inline const LocalDualNodeSnapshot *findSnapshotByNodeId(std::span<const LocalDualNodeSnapshot> snapshots,
                                                         NodeId nodeId) {
    for (const auto &snapshot : snapshots) {
        if (snapshot.nodeId == nodeId) {
            return &snapshot;
        }
    }
    return nullptr;
}

inline std::vector<LocalDualNodeShape> buildLocalDualNodeShapes(
    std::span<const LocalDualNodeSnapshot> snapshots) {
    std::vector<LocalDualNodeShape> shapes;
    shapes.reserve(snapshots.size());

    for (const auto &snapshot : snapshots) {
        LocalDualNodeShape shape;
        shape.altitude = snapshot.altitude;
        shape.directLocalProperParts = snapshot.directLocalProperParts;
        shape.supportLocalPixels = snapshot.supportLocalPixels;
        shape.parentInSnapshot = snapshot.parentInSnapshot;

        if (snapshot.parentInSnapshot) {
            const auto *parent = findSnapshotByNodeId(snapshots, snapshot.parentId);
            if (parent != nullptr) {
                shape.parentAltitude = parent->altitude;
                shape.parentSupportLocalPixels = parent->supportLocalPixels;
            }
        }

        for (NodeId childId : snapshot.childrenInSnapshot) {
            const auto *child = findSnapshotByNodeId(snapshots, childId);
            if (child != nullptr) {
                shape.children.push_back({child->altitude, child->supportLocalPixels});
            }
        }
        std::sort(shape.children.begin(), shape.children.end());
        shapes.push_back(std::move(shape));
    }

    std::sort(shapes.begin(), shapes.end());
    return shapes;
}

inline bool localDualNodeShapesMatch(std::span<const LocalDualNodeSnapshot> lhs,
                                     std::span<const LocalDualNodeSnapshot> rhs) {
    const auto lhsShapes = buildLocalDualNodeShapes(lhs);
    const auto rhsShapes = buildLocalDualNodeShapes(rhs);
    return lhsShapes == rhsShapes;
}

inline void appendLocalDualNodeShape(std::ostringstream &out,
                                     const LocalDualNodeShape &shape) {
    out << "level=" << shape.altitude
        << " direct=";
    appendPixelVector(out, std::span<const PixelId>(shape.directLocalProperParts));
    out << " support=";
    appendPixelVector(out, std::span<const PixelId>(shape.supportLocalPixels));
    out << " parent_in_snapshot=" << (shape.parentInSnapshot ? "true" : "false");
    if (shape.parentInSnapshot) {
        out << " parent_level=" << shape.parentAltitude
            << " parent_support=";
        appendPixelVector(out, std::span<const PixelId>(shape.parentSupportLocalPixels));
    }
    out << " children=[";
    for (std::size_t i = 0; i < shape.children.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        out << "{level=" << shape.children[i].first
            << " support=";
        appendPixelVector(out, std::span<const PixelId>(shape.children[i].second));
        out << "}";
    }
    out << "]";
}

inline std::string localDualNodeShapeMismatchSummary(
    std::span<const LocalDualNodeSnapshot> postLocal,
    std::span<const LocalDualNodeSnapshot> oracleLocal) {
    const auto postShapes = buildLocalDualNodeShapes(postLocal);
    const auto oracleShapes = buildLocalDualNodeShapes(oracleLocal);
    std::ostringstream out;
    out << "post_shapes=" << postShapes.size()
        << " oracle_shapes=" << oracleShapes.size();

    auto appendFirstMissing = [&](const char *label,
                                  const std::vector<LocalDualNodeShape> &lhs,
                                  const std::vector<LocalDualNodeShape> &rhs) {
        for (const auto &shape : lhs) {
            if (std::find(rhs.begin(), rhs.end(), shape) == rhs.end()) {
                out << "\n    " << label << ": ";
                appendLocalDualNodeShape(out, shape);
                return;
            }
        }
        out << "\n    " << label << ": none";
    };

    appendFirstMissing("first_post_only_shape", postShapes, oracleShapes);
    appendFirstMissing("first_oracle_only_shape", oracleShapes, postShapes);
    auto appendShapeList = [&](const char *label,
                               const std::vector<LocalDualNodeShape> &shapes) {
        out << "\n    " << label << "=[";
        const std::size_t limit = std::min<std::size_t>(shapes.size(), 8);
        for (std::size_t i = 0; i < limit; ++i) {
            if (i != 0) {
                out << "; ";
            }
            appendLocalDualNodeShape(out, shapes[i]);
        }
        if (shapes.size() > limit) {
            out << "; ...";
        }
        out << "]";
    };
    appendShapeList("post_shape_list", postShapes);
    appendShapeList("oracle_shape_list", oracleShapes);
    return out.str();
}

inline std::vector<LocalDualNodeSnapshot> buildPostUpdateLocalDualSnapshot(
    const DynamicComponentTree &dualTree,
    const LocalMergeOracle &oracle) {
    const int domainSize = dualTree.getNumTotalProperParts();
    std::vector<PixelId> modifiedPixels;
    std::vector<PixelId> preservedPixels;
    std::vector<PixelId> boundaryPixels;

    for (const auto &owner : oracle.pixelOwners) {
        switch (owner.role) {
        case LocalMergePixelRole::Modified:
            modifiedPixels.push_back(owner.pixelId);
            break;
        case LocalMergePixelRole::Preserved:
            preservedPixels.push_back(owner.pixelId);
            break;
        case LocalMergePixelRole::Boundary:
            boundaryPixels.push_back(owner.pixelId);
            break;
        }
    }

    modifiedPixels = sortedUniquePixels(std::move(modifiedPixels));
    preservedPixels = sortedUniquePixels(std::move(preservedPixels));
    boundaryPixels = sortedUniquePixels(std::move(boundaryPixels));

    const auto modifiedMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(modifiedPixels),
        "post-update local snapshot modified pixels");
    const auto preservedMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(preservedPixels),
        "post-update local snapshot preserved pixels");
    const auto boundaryMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(boundaryPixels),
        "post-update local snapshot boundary pixels");

    return buildLocalDualNodeSnapshot(
        dualTree,
        std::span<const PixelId>(oracle.localPixels),
        modifiedMask,
        preservedMask,
        boundaryMask);
}

inline std::vector<LocalDualNodeSnapshot> buildPostUpdateLocalDualSnapshot(
    const DynamicComponentTree &dualTree,
    const LocalPartitionInput &input) {
    const int domainSize = dualTree.getNumTotalProperParts();
    std::vector<PixelId> modifiedPixels;
    std::vector<PixelId> preservedPixels;
    std::vector<PixelId> boundaryPixels;

    for (const auto &metadata : input.pixelMetadata) {
        switch (metadata.role) {
        case LocalMergePixelRole::Modified:
            modifiedPixels.push_back(metadata.pixelId);
            break;
        case LocalMergePixelRole::Preserved:
            preservedPixels.push_back(metadata.pixelId);
            break;
        case LocalMergePixelRole::Boundary:
            boundaryPixels.push_back(metadata.pixelId);
            break;
        }
    }

    modifiedPixels = sortedUniquePixels(std::move(modifiedPixels));
    preservedPixels = sortedUniquePixels(std::move(preservedPixels));
    boundaryPixels = sortedUniquePixels(std::move(boundaryPixels));

    const auto modifiedMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(modifiedPixels),
        "post-update local input snapshot modified pixels");
    const auto preservedMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(preservedPixels),
        "post-update local input snapshot preserved pixels");
    const auto boundaryMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(boundaryPixels),
        "post-update local input snapshot boundary pixels");

    return buildLocalDualNodeSnapshot(
        dualTree,
        std::span<const PixelId>(input.localPixels),
        modifiedMask,
        preservedMask,
        boundaryMask);
}

inline NonIncreasingStrippedLocalUpdateResult applyStrippedNonOracleLocalMergeUpdate(
    DynamicComponentTree &primalTree,
    DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId nodeId,
    const NonIncreasingStrippedLocalUpdateOptions &options = {}) {
    if (!adjacency) {
        throw std::runtime_error("applyStrippedNonOracleLocalMergeUpdate requires a valid adjacency relation");
    }
    if (nodeId == InvalidNode ||
        nodeId == primalTree.getRoot() ||
        !primalTree.isNode(nodeId) ||
        !primalTree.isAlive(nodeId)) {
        return {};
    }

    const NodeId parentId = primalTree.getNodeParent(nodeId);
    if (parentId == InvalidNode || parentId == nodeId || !primalTree.isAlive(parentId)) {
        return {};
    }

    NonIncreasingStrippedLocalUpdateResult result;

    const auto a = static_cast<AltitudeType>(primalTree.getAltitude(nodeId));
    const auto b = static_cast<AltitudeType>(primalTree.getAltitude(parentId));
    if (a == b || primalTree.getNumProperParts(nodeId) == 0) {
        const auto mergeReleaseStart = std::chrono::steady_clock::now();
        primalTree.mergeNodeIntoParent(nodeId);
        if (options.releaseEmptyNodesAfterPrimalMerge) {
            result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(primalTree);
            result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(dualTree);
        }
        result.primalMergeReleaseMicros =
            localUpdateElapsedMicros(mergeReleaseStart, std::chrono::steady_clock::now());
        result.applied = true;
        result.skippedNoImageChange = true;
        result.localRecompositionMatchesTarget = true;
        const auto finalCheckStart = std::chrono::steady_clock::now();
        if (options.checkFinalDualPrimalReconstruction) {
            result.dualMatchesPrimalAfterMerge =
                dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
            result.finalReconstructionChecked = true;
        }
        result.finalReconstructionCheckMicros =
            localUpdateElapsedMicros(finalCheckStart, std::chrono::steady_clock::now());
        return result;
    }

    std::optional<LocalMergeOracle> rebuildOracle;
    if (options.compareLocalStructureAgainstRebuild) {
        rebuildOracle = buildLocalMergeOracle(primalTree, dualTree, adjacency, nodeId);
        if (!rebuildOracle->preDualMatchesPrimal) {
            throw std::runtime_error(
                "stripped structure comparison requires matching primal/dual reconstructions before merge");
        }
        if (!rebuildOracle->oracleDualMatchesUpdatedImage) {
            throw std::runtime_error(
                "stripped structure comparison could not build a valid rebuilt dual oracle");
        }
    }

    if (options.checkMaterializationTargetReconstruction &&
        !options.buildFullUpdatedImage) {
        throw std::runtime_error(
            "stripped target reconstruction checks require a full updated image");
    }

    auto stageStart = std::chrono::steady_clock::now();
    LocalPartitionInputOptions inputOptions;
    inputOptions.checkPreReconstruction = options.checkInputPreReconstruction;
    inputOptions.buildFullUpdatedImage = options.buildFullUpdatedImage;
    const auto input = makeLocalPartitionInputFromLiveTrees(
        primalTree,
        dualTree,
        adjacency,
        nodeId,
        inputOptions);
    result.inputExtractionMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
    result.inputBuilt = true;
    result.affectedLevelCount = static_cast<int>(affectedLevelsForLocalMerge(input.contraction).size());
    result.localPixelCount = static_cast<int>(input.localPixels.size());
    result.preLocalNodeCount = static_cast<int>(input.preDualNodes.size());

    stageStart = std::chrono::steady_clock::now();
    auto partitionPlan = buildLocalPartitionPlanWithoutOracle(
        input,
        *adjacency);
    result.partitionPlanMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());

    stageStart = std::chrono::steady_clock::now();
    auto recompositionPlan = buildLocalRecompositionPlanWithoutOracle(
        input,
        partitionPlan);
    result.recompositionPlanMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
    result.recompositionNodeCount = static_cast<int>(recompositionPlan.nodes.size());

    stageStart = std::chrono::steady_clock::now();
    if (options.auditPlanOracleDependencies) {
        const auto dependencyReport = auditLocalPlanOracleDependencies(
            &partitionPlan,
            &recompositionPlan);
        result.localPlanOracleFree = dependencyReport.isOracleFree();
        result.localPlanOracleAuditChecked = true;
    }
    result.planAuditMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
    if (options.auditPlanOracleDependencies && !result.localPlanOracleFree) {
        throw std::runtime_error("stripped non-oracle recomposition build exposed oracle-guided fields");
    }

    stageStart = std::chrono::steady_clock::now();
    const auto materializationContext = makeStrippedNonOracleMaterializationContext(
        dualTree,
        input,
        options.checkPreMaterializationReconstruction);
    result.materializationContextMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());

    stageStart = std::chrono::steady_clock::now();
    if (options.auditMaterializationOracleDependencies) {
        const auto materializationDependencyReport = auditLocalMaterializationOracleDependencies(
            materializationContext,
            recompositionPlan);
        result.materializationInputOracleFree =
            materializationDependencyReport.isOracleFreeMaterializationInput();
        result.materializationInputOracleAuditChecked = true;
    }
    result.materializationAuditMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
    if (options.auditMaterializationOracleDependencies &&
        !result.materializationInputOracleFree) {
        throw std::runtime_error("stripped non-oracle materialization exposed oracle-guided fields");
    }

    stageStart = std::chrono::steady_clock::now();
    LocalRecompositionMaterializationOptions materializationOptions;
    materializationOptions.checkPreReconstruction =
        options.checkPreMaterializationReconstruction;
    materializationOptions.checkTargetReconstruction =
        options.checkMaterializationTargetReconstruction;
    materializationOptions.checkNoEmptyLiveNodesAfterCanonicalization =
        options.checkNoEmptyLiveNodesAfterCanonicalization;
    const auto materialization = materializeLocalRecompositionNodes(
        dualTree,
        materializationContext,
        recompositionPlan,
        adjacency.get(),
        materializationOptions);
    result.materializationMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
    result.materializationValidationMicros = materialization.validationMicros;
    result.materializationCreateNodesMicros = materialization.createNodesMicros;
    result.materializationAttachNodesMicros = materialization.attachNodesMicros;
    result.materializationMoveProperPartsMicros = materialization.moveProperPartsMicros;
    result.materializationChildBridgeMicros = materialization.childBridgeMicros;
    result.materializationStructuralRepairMicros = materialization.structuralRepairMicros;
    result.materializationReleaseEmptySourceOwnersMicros =
        materialization.releaseEmptySourceOwnersMicros;
    result.materializationRootPromotionMicros = materialization.rootPromotionMicros;
    result.materializationCanonicalizationMicros = materialization.canonicalizationMicros;
    result.materializationTargetReconstructionCheckMicros =
        materialization.targetReconstructionCheckMicros;
    result.localRecompositionApplied = true;
    result.localRecompositionTargetChecked = materialization.targetReconstructionChecked;
    result.localRecompositionMatchesTarget = materialization.reconstructionMatchesTarget;
    result.postEmptyNodeReleasesRequired =
        materialization.canonicalizationBreakdown.emptyNodeReleases;
    result.postSameLevelParentChildCollapsesRequired =
        materialization.canonicalizationBreakdown.sameLevelParentChildCollapses;
    result.postSameLevelSiblingMergesRequired =
        materialization.canonicalizationBreakdown.sameLevelSiblingMerges;
    result.postComparableSiblingMovesRequired =
        materialization.canonicalizationBreakdown.comparableSiblingMoves;
    result.postCanonicalizationChangesRequired =
        materialization.canonicalizationBreakdown.total();

    if (options.testThrowAfterMaterializationBeforePrimalMerge) {
        throw std::runtime_error(
            "forced test failure after materialization before primal merge");
    }

    recordStrippedCreatedRecompositionNodes(result, materialization);
    recordStrippedMaterializationGuardDiagnostics(result, materialization);
    result.emptyCreatedNodesAfterCanonicalization =
        materialization.emptyCreatedNodesAfterCanonicalization;
    result.emptiedSourceOwnerCount =
        static_cast<int>(materialization.emptiedSourceOwners.size());
    result.releasedSourceOwnerCount =
        static_cast<int>(materialization.releasedSourceOwners.size());
    if (rebuildOracle) {
        stageStart = std::chrono::steady_clock::now();
        const auto postLocalOracleWindow =
            buildPostUpdateLocalDualSnapshot(dualTree, *rebuildOracle);
        result.postLocalSnapshotMicros =
            localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
        result.postLocalNodeCount = static_cast<int>(postLocalOracleWindow.size());
        result.localStructureAgainstRebuildChecked = true;
        result.localStructureMatchesRebuild = localDualNodeShapesMatch(
            std::span<const LocalDualNodeSnapshot>(postLocalOracleWindow),
            std::span<const LocalDualNodeSnapshot>(rebuildOracle->oracleDualNodes));
        if (!result.localStructureMatchesRebuild) {
            result.localStructureMismatchSummary = localDualNodeShapeMismatchSummary(
                std::span<const LocalDualNodeSnapshot>(postLocalOracleWindow),
                std::span<const LocalDualNodeSnapshot>(rebuildOracle->oracleDualNodes));
        }
    }

    if (options.checkMaterializationTargetReconstruction &&
        !materialization.reconstructionMatchesTarget) {
        throw std::runtime_error("stripped non-oracle local recomposition did not reconstruct the local target image");
    }

    const auto mergeReleaseStart = std::chrono::steady_clock::now();
    primalTree.mergeNodeIntoParent(nodeId);
    if (options.releaseEmptyNodesAfterPrimalMerge) {
        result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(primalTree);
        result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(dualTree);
    }
    result.primalMergeReleaseMicros =
        localUpdateElapsedMicros(mergeReleaseStart, std::chrono::steady_clock::now());
    result.applied = true;
    const auto finalCheckStart = std::chrono::steady_clock::now();
    if (options.checkFinalDualPrimalReconstruction) {
        result.dualMatchesPrimalAfterMerge =
            dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
        result.finalReconstructionChecked = true;
    }
    result.finalReconstructionCheckMicros =
        localUpdateElapsedMicros(finalCheckStart, std::chrono::steady_clock::now());
    return result;
}

inline NonIncreasingOracleGuidedLocalUpdateResult applyOracleGuidedLocalMergeUpdate(
    DynamicComponentTree &primalTree,
    DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId nodeId) {
    if (!adjacency) {
        throw std::runtime_error("applyOracleGuidedLocalMergeUpdate requires a valid adjacency relation");
    }
    if (nodeId == InvalidNode ||
        nodeId == primalTree.getRoot() ||
        !primalTree.isNode(nodeId) ||
        !primalTree.isAlive(nodeId)) {
        return {};
    }

    const NodeId parentId = primalTree.getNodeParent(nodeId);
    if (parentId == InvalidNode || parentId == nodeId || !primalTree.isAlive(parentId)) {
        return {};
    }

    NonIncreasingOracleGuidedLocalUpdateResult result;

    const auto a = static_cast<AltitudeType>(primalTree.getAltitude(nodeId));
    const auto b = static_cast<AltitudeType>(primalTree.getAltitude(parentId));
    if (a == b || primalTree.getNumProperParts(nodeId) == 0) {
        primalTree.mergeNodeIntoParent(nodeId);
        result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(primalTree);
        result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(dualTree);
        result.applied = true;
        result.skippedNoImageChange = true;
        result.localRecompositionMatchesTarget = true;
        result.dualMatchesPrimalAfterMerge =
            dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
        return result;
    }

    const auto oracle = buildLocalMergeOracle(primalTree, dualTree, adjacency, nodeId);
    result.oracleBuilt = true;
    result.preDualMatchesPrimal = oracle.preDualMatchesPrimal;
    result.oracleDualMatchesUpdatedImage = oracle.oracleDualMatchesUpdatedImage;
    result.oracleLocalNodeCount = static_cast<int>(oracle.oracleDualNodes.size());

    if (!oracle.preDualMatchesPrimal) {
        throw std::runtime_error("oracle-guided local update requires matching primal/dual reconstructions before merge");
    }
    if (!oracle.oracleDualMatchesUpdatedImage) {
        throw std::runtime_error("oracle-guided local update could not build a valid rebuilt dual oracle");
    }

    const auto partitionPlan = buildLocalPartitionPlan(oracle, *adjacency);
    const auto recompositionPlan = buildLocalRecompositionPlan(oracle, partitionPlan);
    result.affectedLevelCount = static_cast<int>(partitionPlan.affectedLevels.size());
    result.recompositionNodeCount = static_cast<int>(recompositionPlan.nodes.size());

    const auto materialization = materializeLocalRecompositionNodes(
        dualTree,
        oracle,
        recompositionPlan,
        adjacency.get());
    result.localRecompositionApplied = true;
    result.localRecompositionMatchesTarget = materialization.reconstructionMatchesTarget;
    result.createdNodeCount = materialization.numCreatedNodes();
    result.releasedSourceOwnerCount = static_cast<int>(materialization.releasedSourceOwners.size());
    result.postEmptyNodeReleasesRequired =
        materialization.canonicalizationBreakdown.emptyNodeReleases;
    result.postSameLevelParentChildCollapsesRequired =
        materialization.canonicalizationBreakdown.sameLevelParentChildCollapses;
    result.postSameLevelSiblingMergesRequired =
        materialization.canonicalizationBreakdown.sameLevelSiblingMerges;
    result.postComparableSiblingMovesRequired =
        materialization.canonicalizationBreakdown.comparableSiblingMoves;
    result.postCanonicalizationChangesRequired =
        materialization.canonicalizationBreakdown.total();
    const auto postLocalSnapshot = buildPostUpdateLocalDualSnapshot(dualTree, oracle);
    result.postLocalNodeCount = static_cast<int>(postLocalSnapshot.size());
    result.localStructureMatchesOracle = localDualNodeShapesMatch(
        std::span<const LocalDualNodeSnapshot>(postLocalSnapshot),
        std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes));
    if (!result.localStructureMatchesOracle) {
        result.localStructureMismatchSummary = localDualNodeShapeMismatchSummary(
            std::span<const LocalDualNodeSnapshot>(postLocalSnapshot),
            std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes));
    }

    if (!materialization.reconstructionMatchesTarget) {
        throw std::runtime_error("oracle-guided local recomposition did not reconstruct the local target image");
    }

    primalTree.mergeNodeIntoParent(nodeId);
    result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(primalTree);
    result.releasedSourceOwnerCount += releaseAllEmptyLeafNodes(dualTree);
    result.applied = true;
    result.dualMatchesPrimalAfterMerge =
        dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage());
    return result;
}

} // namespace morphotree::details::non_increasing
