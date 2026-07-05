#pragma once

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../AdjacencyRelation.hpp"
#include "../Common.hpp"
#include "../DynamicComponentTree.hpp"
#include "NonIncreasingLocalOracle.hpp"
#include "NonIncreasingLocalPartition.hpp"

namespace morphotree::details::non_increasing {

inline constexpr std::size_t LocalCanonicalizationPixelFastPathLimit = 8192;

inline long long localRecomposerElapsedMicros(std::chrono::steady_clock::time_point start,
                                              std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

enum class LocalCreatedRecompositionNodeKind {
    Planned,
    HiddenExternalParentConnector,
    OracleBridge,
    StableChildBridge,
    HiddenSupportConnector,
    NonOracleHiddenSupportConnector,
    NonOracleParentLinkBridge,
    NonOracleChildLinkBridge
};

struct LocalCreatedRecompositionNode {
    int planIndex = -1;
    NodeId nodeId = InvalidNode;
    NodeId parentNodeId = InvalidNode;
    AltitudeType level = 0;
    LocalCreatedRecompositionNodeKind kind = LocalCreatedRecompositionNodeKind::Planned;
    std::vector<PixelId> movedDirectProperParts;
    std::vector<NodeId> sourceOwnersBefore;
    std::vector<NodeId> movedChildNodes;
};

struct LocalEmptyCreatedNodeDiagnostic {
    int planIndex = -1;
    NodeId nodeId = InvalidNode;
    NodeId parentNodeId = InvalidNode;
    int level = 0;
    int childCount = 0;
    int movedDirectProperPartCount = 0;
    int movedChildNodeCount = 0;
};

struct LocalCanonicalizationBreakdown {
    int emptyNodeReleases = 0;
    int sameLevelParentChildCollapses = 0;
    int sameLevelSiblingMerges = 0;
    int comparableSiblingMoves = 0;

    int total() const {
        return emptyNodeReleases +
               sameLevelParentChildCollapses +
               sameLevelSiblingMerges +
               comparableSiblingMoves;
    }
};

struct LocalRecompositionMaterializationResult {
    std::vector<LocalCreatedRecompositionNode> createdNodes;
    std::vector<NodeId> emptiedSourceOwners;
    std::vector<NodeId> releasedSourceOwners;
    std::vector<LocalEmptyCreatedNodeDiagnostic> emptyCreatedNodesAfterCanonicalization;
    LocalCanonicalizationBreakdown canonicalizationBreakdown;
    bool reconstructionMatchesTarget = false;
    bool targetReconstructionChecked = false;
    int emptyLiveNodesAfterCanonicalization = 0;
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
    long long validationMicros = 0;
    long long createNodesMicros = 0;
    long long attachNodesMicros = 0;
    long long moveProperPartsMicros = 0;
    long long childBridgeMicros = 0;
    long long structuralRepairMicros = 0;
    long long releaseEmptySourceOwnersMicros = 0;
    long long rootPromotionMicros = 0;
    long long canonicalizationMicros = 0;
    long long targetReconstructionCheckMicros = 0;

    int numCreatedNodes() const {
        return static_cast<int>(createdNodes.size());
    }

    int numEmptyConnectorNodesAfterCanonicalization() const {
        return static_cast<int>(std::count_if(
            emptyCreatedNodesAfterCanonicalization.begin(),
            emptyCreatedNodesAfterCanonicalization.end(),
            [](const auto &node) {
                return node.planIndex < 0;
            }));
    }
};

struct LocalRecompositionMaterializationOptions {
    bool checkPreReconstruction = true;
    bool checkTargetReconstruction = true;
    bool checkNoEmptyLiveNodesAfterCanonicalization = false;
    bool releaseEmptySourceOwnersWithResidualChildren = false;
};

struct StableOracleChildAttachment {
    NodeId oracleChildId = InvalidNode;
    NodeId childNodeId = InvalidNode;
};

inline std::vector<NodeId> sortedUniqueValidNodes(std::vector<NodeId> nodes) {
    nodes.erase(std::remove(nodes.begin(), nodes.end(), InvalidNode), nodes.end());
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return nodes;
}

inline std::vector<NodeId> stableUniqueValidNodes(std::vector<NodeId> nodes) {
    std::vector<NodeId> uniqueNodes;
    uniqueNodes.reserve(nodes.size());
    std::unordered_set<NodeId> seen;
    seen.reserve(nodes.size());
    for (NodeId nodeId : nodes) {
        if (nodeId == InvalidNode) {
            continue;
        }
        if (seen.insert(nodeId).second) {
            uniqueNodes.push_back(nodeId);
        }
    }
    return uniqueNodes;
}

inline int countAliveNodesWithoutProperParts(const DynamicComponentTree &tree) {
    int count = 0;
    const int nodeSlots = tree.getNumInternalNodeSlots();
    for (NodeId nodeId = 0; nodeId < nodeSlots; ++nodeId) {
        if (!tree.isNode(nodeId) ||
            !tree.isAlive(nodeId) ||
            tree.getNumProperParts(nodeId) != 0) {
            continue;
        }
        ++count;
    }
    return count;
}

inline std::string describeAliveNodesWithoutProperParts(const DynamicComponentTree &tree,
                                                        int maxExamples = 8) {
    std::ostringstream out;
    int count = 0;
    const int nodeSlots = tree.getNumInternalNodeSlots();
    for (NodeId nodeId = 0; nodeId < nodeSlots; ++nodeId) {
        if (!tree.isNode(nodeId) ||
            !tree.isAlive(nodeId) ||
            tree.getNumProperParts(nodeId) != 0) {
            continue;
        }
        if (count < maxExamples) {
            out << " id=" << nodeId
                << " parent=" << tree.getNodeParent(nodeId)
                << " altitude=" << tree.getAltitude(nodeId)
                << " children=" << tree.getNumChildren(nodeId);
        }
        ++count;
    }
    out << " count=" << count;
    return out.str();
}

inline std::vector<LocalEmptyCreatedNodeDiagnostic> collectEmptyCreatedNodeDiagnostics(
    const DynamicComponentTree &tree,
    std::span<const LocalCreatedRecompositionNode> createdNodes) {
    std::vector<LocalEmptyCreatedNodeDiagnostic> diagnostics;
    for (const auto &created : createdNodes) {
        if (created.nodeId == InvalidNode ||
            !tree.isNode(created.nodeId) ||
            !tree.isAlive(created.nodeId) ||
            tree.getNumProperParts(created.nodeId) != 0) {
            continue;
        }
        LocalEmptyCreatedNodeDiagnostic diagnostic;
        diagnostic.planIndex = created.planIndex;
        diagnostic.nodeId = created.nodeId;
        diagnostic.parentNodeId = tree.getNodeParent(created.nodeId);
        diagnostic.level = tree.getAltitude(created.nodeId);
        diagnostic.childCount = tree.getNumChildren(created.nodeId);
        diagnostic.movedDirectProperPartCount =
            static_cast<int>(created.movedDirectProperParts.size());
        diagnostic.movedChildNodeCount =
            static_cast<int>(created.movedChildNodes.size());
        diagnostics.push_back(diagnostic);
    }
    return diagnostics;
}

inline NodeId chooseClosestValidExternalParent(const DynamicComponentTree &tree,
                                               const LocalRecompositionNodePlan &nodePlan,
                                               std::vector<NodeId> candidates) {
    candidates = sortedUniqueValidNodes(std::move(candidates));
    NodeId best = InvalidNode;
    int bestLevelDistance = 0;
    for (NodeId candidate : candidates) {
        if (!tree.isAlive(candidate)) {
            continue;
        }
        const auto candidateLevel = static_cast<AltitudeType>(tree.getAltitude(candidate));
        if (!validLocalParentLevel(nodePlan.level, candidateLevel, tree.isMaxtree())) {
            continue;
        }
        const int levelDistance = std::abs(static_cast<int>(candidateLevel) - static_cast<int>(nodePlan.level));
        if (best == InvalidNode || levelDistance < bestLevelDistance) {
            best = candidate;
            bestLevelDistance = levelDistance;
        }
    }
    return best;
}

inline std::vector<NodeId> collectPreOwnersForOracleOwner(const LocalMergeOracle &oracle,
                                                          NodeId oracleOwnerId) {
    std::vector<NodeId> owners;
    if (oracleOwnerId == InvalidNode) {
        return owners;
    }
    for (const auto &pixelOwner : oracle.pixelOwners) {
        if (pixelOwner.oracleDualOwner == oracleOwnerId) {
            owners.push_back(pixelOwner.preDualOwner);
        }
    }
    return sortedUniqueValidNodes(std::move(owners));
}

inline const LocalDualNodeSnapshot *findLocalSnapshotByNodeId(
    std::span<const LocalDualNodeSnapshot> snapshots,
    NodeId nodeId) {
    for (const auto &snapshot : snapshots) {
        if (snapshot.nodeId == nodeId) {
            return &snapshot;
        }
    }
    return nullptr;
}

inline NodeId findStablePreNodeForOracleNode(const LocalMergeOracle &oracle,
                                             NodeId oracleNodeId) {
    const auto *oracleSnapshot = findLocalSnapshotByNodeId(
        std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
        oracleNodeId);
    if (oracleSnapshot == nullptr) {
        return InvalidNode;
    }

    for (const auto &preSnapshot : oracle.preDualNodes) {
        if (preSnapshot.altitude == oracleSnapshot->altitude &&
            preSnapshot.supportLocalPixels == oracleSnapshot->supportLocalPixels) {
            return preSnapshot.nodeId;
        }
    }
    return InvalidNode;
}

inline void recordEmptyPreviousParent(DynamicComponentTree &tree,
                                      NodeId previousParent,
                                      NodeId newParent,
                                      std::vector<NodeId> &emptiedSourceOwners) {
    if (previousParent != InvalidNode &&
        previousParent != newParent &&
        tree.isAlive(previousParent) &&
        !tree.isRoot(previousParent) &&
        tree.getNumProperParts(previousParent) == 0 &&
        tree.getNumChildren(previousParent) == 0) {
        emptiedSourceOwners.push_back(previousParent);
    }
}

inline std::vector<NodeId> collectExistingChildNodesForOracleNode(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    NodeId oracleNodeId,
    NodeId targetParentNodeId,
    NodeId excludedNodeId = InvalidNode) {
    std::vector<NodeId> candidates;
    const NodeId stableNode = findStablePreNodeForOracleNode(oracle, oracleNodeId);
    if (stableNode != InvalidNode) {
        candidates.push_back(stableNode);
    }
    auto ownerCandidates = collectPreOwnersForOracleOwner(oracle, oracleNodeId);
    candidates.insert(candidates.end(), ownerCandidates.begin(), ownerCandidates.end());
    candidates = sortedUniqueValidNodes(std::move(candidates));

    std::vector<NodeId> validCandidates;
    if (targetParentNodeId == InvalidNode || !tree.isAlive(targetParentNodeId)) {
        return validCandidates;
    }
    const auto parentLevel = static_cast<AltitudeType>(tree.getAltitude(targetParentNodeId));
    for (NodeId candidate : candidates) {
        if (candidate == InvalidNode ||
            candidate == excludedNodeId ||
            candidate == targetParentNodeId ||
            !tree.isAlive(candidate) ||
            tree.isRoot(candidate)) {
            continue;
        }
        if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(candidate)),
                                   parentLevel,
                                   tree.isMaxtree())) {
            continue;
        }
        validCandidates.push_back(candidate);
    }
    return sortedUniqueValidNodes(std::move(validCandidates));
}

inline NodeId externalAttachmentForRecompositionNode(const DynamicComponentTree &tree,
                                                     const LocalMergeOracle &oracle,
                                                     const LocalRecompositionNodePlan &nodePlan) {
    const NodeId stableSnapshotParent = findStablePreNodeForOracleNode(
        oracle,
        nodePlan.oracleParentId);
    if (stableSnapshotParent != InvalidNode &&
        tree.isAlive(stableSnapshotParent) &&
        validLocalParentLevel(nodePlan.level,
                              static_cast<AltitudeType>(tree.getAltitude(stableSnapshotParent)),
                              tree.isMaxtree())) {
        return stableSnapshotParent;
    }

    const NodeId mappedParent = chooseClosestValidExternalParent(
        tree,
        nodePlan,
        collectPreOwnersForOracleOwner(oracle, nodePlan.oracleParentId));
    if (mappedParent != InvalidNode) {
        return mappedParent;
    }

    const NodeId bridgeParent = chooseClosestValidExternalParent(
        tree,
        nodePlan,
        nodePlan.preDualOwners);
    if (bridgeParent != InvalidNode) {
        return bridgeParent;
    }

    const NodeId root = tree.getRoot();
    if (root == InvalidNode || !tree.isAlive(root)) {
        throw std::runtime_error("local recomposition requires a live root attachment point");
    }
    if (!validLocalParentLevel(nodePlan.level,
                               static_cast<AltitudeType>(tree.getAltitude(root)),
                               tree.isMaxtree())) {
        throw std::runtime_error("local recomposition root is not a valid external parent for the node level");
    }
    return root;
}

inline std::vector<StableOracleChildAttachment> collectStableChildAttachmentsForRecompositionNode(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionNodePlan &nodePlan,
    NodeId createdNodeId) {
    std::vector<StableOracleChildAttachment> attachments;
    for (NodeId oracleChildId : nodePlan.oracleChildrenInSnapshot) {
        const auto candidates = collectExistingChildNodesForOracleNode(
            tree,
            oracle,
            oracleChildId,
            createdNodeId,
            createdNodeId);
        for (NodeId candidate : candidates) {
            attachments.push_back({oracleChildId, candidate});
        }
    }

    std::sort(attachments.begin(), attachments.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.oracleChildId != rhs.oracleChildId) {
            return lhs.oracleChildId < rhs.oracleChildId;
        }
        return lhs.childNodeId < rhs.childNodeId;
    });
    attachments.erase(
        std::unique(attachments.begin(),
                    attachments.end(),
                    [](const auto &lhs, const auto &rhs) {
                        return lhs.oracleChildId == rhs.oracleChildId &&
                               lhs.childNodeId == rhs.childNodeId;
                    }),
        attachments.end());
    return attachments;
}

inline std::vector<NodeId> moveStablePreChildNodesForRecompositionNode(
    DynamicComponentTree &tree,
    const LocalRecompositionNodePlan &nodePlan,
    NodeId createdNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners) {
    std::vector<NodeId> movedChildren;
    if (createdNodeId == InvalidNode || !tree.isAlive(createdNodeId)) {
        return movedChildren;
    }

    const auto createdLevel = static_cast<AltitudeType>(tree.getAltitude(createdNodeId));
    for (NodeId stableChildId : sortedUniqueValidNodes(nodePlan.stableChildPreNodes)) {
        if (stableChildId == createdNodeId ||
            !tree.isAlive(stableChildId) ||
            tree.isRoot(stableChildId)) {
            continue;
        }
        if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(stableChildId)),
                                   createdLevel,
                                   tree.isMaxtree())) {
            continue;
        }

        if (tree.getNodeParent(stableChildId) != createdNodeId) {
            const NodeId previousParent = tree.getNodeParent(stableChildId);
            tree.moveNode(stableChildId, createdNodeId);
            recordEmptyPreviousParent(
                tree,
                previousParent,
                createdNodeId,
                emptiedSourceOwners);
        }
        canonicalFocusNodeIds.push_back(stableChildId);
        movedChildren.push_back(stableChildId);
    }

    return sortedUniqueValidNodes(std::move(movedChildren));
}

inline bool oracleNodeIsMaterializedByPlan(const LocalRecompositionPlan &plan,
                                           NodeId oracleNodeId) {
    if (oracleNodeId == InvalidNode) {
        return false;
    }
    return std::any_of(plan.nodes.begin(), plan.nodes.end(), [&](const auto &nodePlan) {
        return nodePlan.oracleNodeId == oracleNodeId;
    });
}

inline NodeId findMaterializedBridgeNode(
    const std::vector<std::pair<NodeId, NodeId>> &bridgeNodes,
    NodeId oracleNodeId) {
    for (const auto &[oracleId, bridgeId] : bridgeNodes) {
        if (oracleId == oracleNodeId) {
            return bridgeId;
        }
    }
    return InvalidNode;
}

inline void materializeMissingOracleDescendantBridges(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    NodeId oracleNodeId,
    NodeId materializedNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &bridgeNodes);

inline NodeId materializeOracleBridgeSubtree(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    NodeId oracleNodeId,
    NodeId parentNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &bridgeNodes) {
    if (oracleNodeId == InvalidNode ||
        parentNodeId == InvalidNode ||
        oracleNodeIsMaterializedByPlan(plan, oracleNodeId)) {
        return InvalidNode;
    }

    const NodeId existingBridge = findMaterializedBridgeNode(bridgeNodes, oracleNodeId);
    if (existingBridge != InvalidNode) {
        if (tree.isAlive(existingBridge) && tree.getNodeParent(existingBridge) != parentNodeId) {
            tree.moveNode(existingBridge, parentNodeId);
        }
        return existingBridge;
    }

    const auto *oracleSnapshot = findLocalSnapshotByNodeId(
        std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
        oracleNodeId);
    if (oracleSnapshot == nullptr) {
        return InvalidNode;
    }

    const auto bridgeLevel = static_cast<AltitudeType>(oracleSnapshot->altitude);
    if (!tree.isAlive(parentNodeId) ||
        !validLocalParentLevel(bridgeLevel,
                               static_cast<AltitudeType>(tree.getAltitude(parentNodeId)),
                               tree.isMaxtree())) {
        return InvalidNode;
    }

    const NodeId bridgeNodeId = tree.createDetachedNode(static_cast<int>(bridgeLevel));
    if (bridgeNodeId == InvalidNode) {
        throw std::runtime_error("local recomposition materialization could not create an oracle bridge node");
    }
    tree.attachNode(parentNodeId, bridgeNodeId);
    canonicalFocusNodeIds.push_back(bridgeNodeId);
    bridgeNodes.push_back({oracleNodeId, bridgeNodeId});

    LocalCreatedRecompositionNode created;
    created.planIndex = -1;
    created.nodeId = bridgeNodeId;
    created.parentNodeId = parentNodeId;
    created.level = bridgeLevel;
    created.kind = LocalCreatedRecompositionNodeKind::OracleBridge;

    for (PixelId pixelId : oracleSnapshot->directLocalProperParts) {
        const NodeId sourceOwner = tree.getSmallestComponent(pixelId);
        if (sourceOwner == InvalidNode || !tree.isAlive(sourceOwner)) {
            throw std::runtime_error("local recomposition bridge cannot locate a live source owner");
        }
        created.sourceOwnersBefore.push_back(sourceOwner);
        if (sourceOwner != bridgeNodeId) {
            tree.moveProperPart(bridgeNodeId, sourceOwner, pixelId);
            created.movedDirectProperParts.push_back(pixelId);
            if (tree.getNumProperParts(sourceOwner) == 0) {
                emptiedSourceOwners.push_back(sourceOwner);
            }
        }
    }
    created.movedDirectProperParts =
        sortedUniquePixels(std::move(created.movedDirectProperParts));
    created.sourceOwnersBefore =
        sortedUniqueValidNodes(std::move(created.sourceOwnersBefore));

    for (NodeId oracleChildId : oracleSnapshot->childrenInSnapshot) {
        if (oracleNodeIsMaterializedByPlan(plan, oracleChildId)) {
            continue;
        }

        auto candidates = collectExistingChildNodesForOracleNode(
            tree,
            oracle,
            oracleChildId,
            bridgeNodeId,
            bridgeNodeId);
        bool movedExistingChild = false;
        for (NodeId candidate : candidates) {
            if (tree.getNodeParent(candidate) != bridgeNodeId) {
                const NodeId previousParent = tree.getNodeParent(candidate);
                tree.moveNode(candidate, bridgeNodeId);
                recordEmptyPreviousParent(
                    tree,
                    previousParent,
                    bridgeNodeId,
                    emptiedSourceOwners);
            }
            canonicalFocusNodeIds.push_back(candidate);
            created.movedChildNodes.push_back(candidate);
            materializeMissingOracleDescendantBridges(
                tree,
                oracle,
                plan,
                oracleChildId,
                candidate,
                canonicalFocusNodeIds,
                emptiedSourceOwners,
                result,
                bridgeNodes);
            movedExistingChild = true;
        }

        if (!movedExistingChild) {
            const NodeId childBridge = materializeOracleBridgeSubtree(
                tree,
                oracle,
                plan,
                oracleChildId,
                bridgeNodeId,
                canonicalFocusNodeIds,
                emptiedSourceOwners,
                result,
                bridgeNodes);
            if (childBridge != InvalidNode) {
                created.movedChildNodes.push_back(childBridge);
            }
        }
    }

    created.movedChildNodes =
        sortedUniqueValidNodes(std::move(created.movedChildNodes));
    result.createdNodes.push_back(std::move(created));
    return bridgeNodeId;
}

inline void materializeMissingOracleDescendantBridges(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    NodeId oracleNodeId,
    NodeId materializedNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &bridgeNodes) {
    if (oracleNodeId == InvalidNode ||
        materializedNodeId == InvalidNode ||
        !tree.isAlive(materializedNodeId)) {
        return;
    }

    const auto *oracleSnapshot = findLocalSnapshotByNodeId(
        std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
        oracleNodeId);
    if (oracleSnapshot == nullptr) {
        return;
    }

    for (NodeId oracleChildId : oracleSnapshot->childrenInSnapshot) {
        if (oracleNodeIsMaterializedByPlan(plan, oracleChildId)) {
            continue;
        }

        const auto existingChildren = collectExistingChildNodesForOracleNode(
            tree,
            oracle,
            oracleChildId,
            materializedNodeId,
            materializedNodeId);
        if (!existingChildren.empty()) {
            for (NodeId childNodeId : existingChildren) {
                if (tree.getNodeParent(childNodeId) != materializedNodeId) {
                    const NodeId previousParent = tree.getNodeParent(childNodeId);
                    tree.moveNode(childNodeId, materializedNodeId);
                    recordEmptyPreviousParent(
                        tree,
                        previousParent,
                        materializedNodeId,
                        emptiedSourceOwners);
                }
                canonicalFocusNodeIds.push_back(childNodeId);
                materializeMissingOracleDescendantBridges(
                    tree,
                    oracle,
                    plan,
                    oracleChildId,
                    childNodeId,
                    canonicalFocusNodeIds,
                    emptiedSourceOwners,
                    result,
                    bridgeNodes);
            }
            continue;
        }

        (void) materializeOracleBridgeSubtree(
            tree,
            oracle,
            plan,
            oracleChildId,
            materializedNodeId,
            canonicalFocusNodeIds,
            emptiedSourceOwners,
            result,
            bridgeNodes);
    }
}

inline std::vector<NodeId> materializeMissingOracleChildBridges(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    const LocalRecompositionNodePlan &nodePlan,
    NodeId createdNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &bridgeNodes) {
    std::vector<NodeId> bridges;
    for (NodeId oracleChildId : nodePlan.oracleChildrenInSnapshot) {
        if (oracleNodeIsMaterializedByPlan(plan, oracleChildId)) {
            continue;
        }
        const auto existingCandidates = collectPreOwnersForOracleOwner(oracle, oracleChildId);
        if (!existingCandidates.empty()) {
            continue;
        }
        const NodeId bridgeNodeId = materializeOracleBridgeSubtree(
            tree,
            oracle,
            plan,
            oracleChildId,
            createdNodeId,
            canonicalFocusNodeIds,
            emptiedSourceOwners,
            result,
            bridgeNodes);
        if (bridgeNodeId != InvalidNode) {
            bridges.push_back(bridgeNodeId);
        }
    }
    return sortedUniqueValidNodes(std::move(bridges));
}

inline NodeId findPlanMaterializedNodeForOracleNode(const LocalRecompositionPlan &plan,
                                                    std::span<const NodeId> createdNodeIds,
                                                    NodeId oracleNodeId) {
    if (oracleNodeId == InvalidNode) {
        return InvalidNode;
    }
    for (std::size_t index = 0; index < plan.nodes.size() && index < createdNodeIds.size(); ++index) {
        if (plan.nodes[index].oracleNodeId == oracleNodeId) {
            return createdNodeIds[index];
        }
    }
    return InvalidNode;
}

inline NodeId findMaterializedOrExistingOracleNode(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    const std::vector<std::pair<NodeId, NodeId>> &bridgeNodes,
    NodeId oracleNodeId) {
    const NodeId planNode = findPlanMaterializedNodeForOracleNode(
        plan,
        createdNodeIds,
        oracleNodeId);
    if (planNode != InvalidNode && tree.isAlive(planNode)) {
        return planNode;
    }

    const NodeId bridgeNode = findMaterializedBridgeNode(bridgeNodes, oracleNodeId);
    if (bridgeNode != InvalidNode && tree.isAlive(bridgeNode)) {
        return bridgeNode;
    }

    const auto *oracleSnapshot = findLocalSnapshotByNodeId(
        std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
        oracleNodeId);
    if (oracleSnapshot == nullptr) {
        return InvalidNode;
    }

    auto candidateMatchesOracleLevel = [&](NodeId candidate) {
        return candidate != InvalidNode &&
               tree.isNode(candidate) &&
               tree.isAlive(candidate) &&
               tree.getAltitude(candidate) == oracleSnapshot->altitude;
    };

    const NodeId stableNode = findStablePreNodeForOracleNode(oracle, oracleNodeId);
    if (candidateMatchesOracleLevel(stableNode)) {
        return stableNode;
    }

    const auto candidates = collectPreOwnersForOracleOwner(oracle, oracleNodeId);
    NodeId selected = InvalidNode;
    bool ambiguous = false;
    for (NodeId candidate : candidates) {
        if (!candidateMatchesOracleLevel(candidate)) {
            continue;
        }
        if (selected != InvalidNode && selected != candidate) {
            ambiguous = true;
            break;
        }
        selected = candidate;
    }
    if (!ambiguous && selected != InvalidNode) {
        return selected;
    }

    for (NodeId oracleChildId : oracleSnapshot->childrenInSnapshot) {
        const NodeId childNodeId = findMaterializedOrExistingOracleNode(
            tree,
            oracle,
            plan,
            createdNodeIds,
            bridgeNodes,
            oracleChildId);
        if (childNodeId == InvalidNode ||
            !tree.isAlive(childNodeId) ||
            tree.isRoot(childNodeId)) {
            continue;
        }
        const NodeId currentParent = tree.getNodeParent(childNodeId);
        if (currentParent != InvalidNode &&
            currentParent != childNodeId &&
            tree.isAlive(currentParent) &&
            tree.getAltitude(currentParent) == oracleSnapshot->altitude) {
            return currentParent;
        }
    }
    return InvalidNode;
}

inline NodeId materializeMissingOracleParentBridgeAboveChild(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    NodeId parentOracleNodeId,
    NodeId childNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &bridgeNodes) {
    if (parentOracleNodeId == InvalidNode ||
        childNodeId == InvalidNode ||
        !tree.isAlive(childNodeId) ||
        tree.isRoot(childNodeId)) {
        return InvalidNode;
    }

    const auto *parentSnapshot = findLocalSnapshotByNodeId(
        std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
        parentOracleNodeId);
    if (parentSnapshot == nullptr) {
        return InvalidNode;
    }

    const NodeId currentParent = tree.getNodeParent(childNodeId);
    if (currentParent == InvalidNode ||
        currentParent == childNodeId ||
        !tree.isAlive(currentParent)) {
        return InvalidNode;
    }

    const auto bridgeLevel = static_cast<AltitudeType>(parentSnapshot->altitude);
    if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(childNodeId)),
                               bridgeLevel,
                               tree.isMaxtree()) ||
        !validLocalParentLevel(bridgeLevel,
                               static_cast<AltitudeType>(tree.getAltitude(currentParent)),
                               tree.isMaxtree())) {
        return InvalidNode;
    }

    return materializeOracleBridgeSubtree(
        tree,
        oracle,
        plan,
        parentOracleNodeId,
        currentParent,
        canonicalFocusNodeIds,
        emptiedSourceOwners,
        result,
        bridgeNodes);
}

inline bool oracleParentIsOutsideLocalSnapshot(const LocalMergeOracle &oracle,
                                               NodeId oracleParentId) {
    return oracleParentId != InvalidNode &&
           findLocalSnapshotByNodeId(
               std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
               oracleParentId) == nullptr;
}

inline int hiddenConnectorLevelBetween(int parentLevel,
                                       int childLevel,
                                       bool isMaxtree);

inline int hiddenExternalParentConnectorLevel(int childLevel,
                                              int externalParentLevel,
                                              bool isMaxtree) {
    if (isMaxtree) {
        if (childLevel - 1 <= externalParentLevel) {
            return -1;
        }
        return childLevel - 1;
    }
    if (childLevel + 1 >= externalParentLevel) {
        return -1;
    }
    return childLevel + 1;
}

inline NodeId materializeHiddenExternalParentConnector(
    DynamicComponentTree &tree,
    const LocalRecompositionNodePlan &nodePlan,
    NodeId externalParentId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<LocalCreatedRecompositionNode> &hiddenCreatedNodes) {
    if (externalParentId == InvalidNode ||
        !tree.isAlive(externalParentId)) {
        return InvalidNode;
    }

    const int connectorLevel = hiddenExternalParentConnectorLevel(
        static_cast<int>(nodePlan.level),
        tree.getAltitude(externalParentId),
        tree.isMaxtree());
    if (connectorLevel < 0) {
        return InvalidNode;
    }

    const NodeId connectorId = tree.createDetachedNode(connectorLevel);
    if (connectorId == InvalidNode) {
        throw std::runtime_error("local recomposition could not create a hidden external parent connector");
    }
    tree.attachNode(externalParentId, connectorId);
    canonicalFocusNodeIds.push_back(connectorId);

    LocalCreatedRecompositionNode created;
    created.planIndex = -1;
    created.nodeId = connectorId;
    created.parentNodeId = externalParentId;
    created.level = static_cast<AltitudeType>(connectorLevel);
    created.kind = LocalCreatedRecompositionNodeKind::HiddenExternalParentConnector;
    hiddenCreatedNodes.push_back(std::move(created));
    return connectorId;
}

inline bool shouldMaterializeNonOracleHiddenExternalParentConnector(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionNodePlan &nodePlan,
    NodeId externalParentId) {
    if (nodePlan.oracleNodeId != InvalidNode ||
        externalParentId == InvalidNode ||
        !tree.isAlive(externalParentId) ||
        tree.isRoot(externalParentId)) {
        return false;
    }

    const auto *externalParentSnapshot = findLocalSnapshotByNodeId(
        std::span<const LocalDualNodeSnapshot>(oracle.preDualNodes),
        externalParentId);
    if (externalParentSnapshot == nullptr ||
        externalParentSnapshot->supportLocalPixels.empty()) {
        return false;
    }

    auto parentSupport = sortedUniquePixels(externalParentSnapshot->supportLocalPixels);
    auto nodeSupport = sortedUniquePixels(nodePlan.supportPixels);
    if (nodeSupport.empty() ||
        equalPixelSets(std::span<const PixelId>(parentSupport),
                       std::span<const PixelId>(nodeSupport)) ||
        !containsAllPixels(std::span<const PixelId>(parentSupport),
                           std::span<const PixelId>(nodeSupport))) {
        return false;
    }

    return hiddenExternalParentConnectorLevel(
               static_cast<int>(nodePlan.level),
               tree.getAltitude(externalParentId),
               tree.isMaxtree()) >= 0;
}

inline bool supportContainsPixel(const DynamicComponentTree &tree,
                                 NodeId nodeId,
                                 PixelId pixelId) {
    if (nodeId == InvalidNode ||
        pixelId < 0 ||
        pixelId >= tree.getNumTotalProperParts() ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId)) {
        return false;
    }

    NodeId current = tree.getSmallestComponent(pixelId);
    while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
        if (current == nodeId) {
            return true;
        }
        const NodeId parentId = tree.getNodeParent(current);
        if (parentId == current) {
            break;
        }
        current = parentId;
    }
    return false;
}

inline NodeId closestDetachedOracleSupportAncestor(const LocalMergeOracle &oracle,
                                                   PixelId pixelId,
                                                   NodeId targetOracleNodeId,
                                                   int targetLevel,
                                                   bool isMaxtree) {
    NodeId best = InvalidNode;
    int bestDistance = 0;
    for (const auto &snapshot : oracle.oracleDualNodes) {
        if (snapshot.nodeId == targetOracleNodeId ||
            snapshot.parentInSnapshot ||
            !std::binary_search(snapshot.supportLocalPixels.begin(),
                                snapshot.supportLocalPixels.end(),
                                pixelId)) {
            continue;
        }
        if (!validLocalParentLevel(static_cast<AltitudeType>(snapshot.altitude),
                                   static_cast<AltitudeType>(targetLevel),
                                   isMaxtree)) {
            continue;
        }

        const int distance = std::abs(targetLevel - snapshot.altitude);
        if (best == InvalidNode || distance < bestDistance) {
            best = snapshot.nodeId;
            bestDistance = distance;
        }
    }
    return best;
}

inline int hiddenConnectorLevelBetween(int parentLevel,
                                       int childLevel,
                                       bool isMaxtree) {
    if (isMaxtree) {
        if (parentLevel + 1 >= childLevel) {
            return -1;
        }
        return parentLevel + 1;
    }
    if (parentLevel - 1 <= childLevel) {
        return -1;
    }
    return parentLevel - 1;
}

inline int intermediateLocalBridgeLevelBetween(int parentLevel,
                                               int childLevel,
                                               bool isMaxtree) {
    if (!validLocalParentLevel(static_cast<AltitudeType>(childLevel),
                               static_cast<AltitudeType>(parentLevel),
                               isMaxtree)) {
        return -1;
    }

    const int low = std::min(parentLevel, childLevel);
    const int high = std::max(parentLevel, childLevel);
    if (high - low <= 1) {
        return -1;
    }
    return low + ((high - low) / 2);
}

inline bool snapshotSupportStrictlyContains(const LocalDualNodeSnapshot &parentSnapshot,
                                            const LocalDualNodeSnapshot &childSnapshot) {
    auto parentSupport = sortedUniquePixels(parentSnapshot.supportLocalPixels);
    auto childSupport = sortedUniquePixels(childSnapshot.supportLocalPixels);
    return !equalPixelSets(std::span<const PixelId>(parentSupport),
                           std::span<const PixelId>(childSupport)) &&
           containsAllPixels(std::span<const PixelId>(parentSupport),
                             std::span<const PixelId>(childSupport));
}

inline void materializeStablePreDescendantBridges(
    DynamicComponentTree &tree,
    std::span<const LocalDualNodeSnapshot> preSnapshots,
    NodeId parentNodeId,
    const LocalDualNodeSnapshot &parentSnapshot,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result) {
    if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
        return;
    }

    for (NodeId childNodeId : parentSnapshot.childrenInSnapshot) {
        const auto *childSnapshot = findLocalSnapshotByNodeId(preSnapshots, childNodeId);
        if (childSnapshot == nullptr ||
            childNodeId == InvalidNode ||
            !tree.isAlive(childNodeId) ||
            tree.isRoot(childNodeId) ||
            tree.getNodeParent(childNodeId) != parentNodeId ||
            !snapshotSupportStrictlyContains(parentSnapshot, *childSnapshot)) {
            continue;
        }

        const int bridgeLevel = intermediateLocalBridgeLevelBetween(
            tree.getAltitude(parentNodeId),
            tree.getAltitude(childNodeId),
            tree.isMaxtree());
        if (bridgeLevel >= 0 &&
            validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(childNodeId)),
                                  static_cast<AltitudeType>(bridgeLevel),
                                  tree.isMaxtree()) &&
            validLocalParentLevel(static_cast<AltitudeType>(bridgeLevel),
                                  static_cast<AltitudeType>(tree.getAltitude(parentNodeId)),
                                  tree.isMaxtree())) {
            canonicalFocusNodeIds.push_back(parentNodeId);
            canonicalFocusNodeIds.push_back(childNodeId);
        }

        materializeStablePreDescendantBridges(
            tree,
            preSnapshots,
            childNodeId,
            *childSnapshot,
            canonicalFocusNodeIds,
            emptiedSourceOwners,
            result);
    }
}

inline void materializeStablePreChildDescendantBridges(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    std::span<const NodeId> stableChildPreNodes,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result) {
    if (oracle.preDualNodes.empty()) {
        return;
    }

    const std::span<const LocalDualNodeSnapshot> preSnapshots(oracle.preDualNodes);
    for (NodeId stableChildId : stableChildPreNodes) {
        const auto *stableSnapshot = findLocalSnapshotByNodeId(preSnapshots, stableChildId);
        if (stableSnapshot == nullptr) {
            continue;
        }
        materializeStablePreDescendantBridges(
            tree,
            preSnapshots,
            stableChildId,
            *stableSnapshot,
            canonicalFocusNodeIds,
            emptiedSourceOwners,
            result);
    }
}

// Hidden connectors make detached oracle branches contribute to an ancestor
// support without turning the detached branch root into a visible snapshot child.
inline void materializeHiddenSupportConnectors(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    const std::vector<std::pair<NodeId, NodeId>> &bridgeNodes) {
    std::vector<NodeId> movedBranchRoots;
    std::vector<const LocalDualNodeSnapshot *> targetSnapshots;
    targetSnapshots.reserve(oracle.oracleDualNodes.size());
    for (const auto &snapshot : oracle.oracleDualNodes) {
        targetSnapshots.push_back(&snapshot);
    }
    std::sort(targetSnapshots.begin(), targetSnapshots.end(), [&](const auto *lhs, const auto *rhs) {
        if (lhs->altitude != rhs->altitude) {
            return tree.isMaxtree() ? lhs->altitude > rhs->altitude
                                    : lhs->altitude < rhs->altitude;
        }
        return lhs->nodeId < rhs->nodeId;
    });

    for (const auto *targetSnapshot : targetSnapshots) {
        const NodeId targetNodeId = findMaterializedOrExistingOracleNode(
            tree,
            oracle,
            plan,
            createdNodeIds,
            bridgeNodes,
            targetSnapshot->nodeId);
        if (targetNodeId == InvalidNode || !tree.isAlive(targetNodeId)) {
            continue;
        }

        for (PixelId pixelId : targetSnapshot->supportLocalPixels) {
            const NodeId detachedAncestorId = closestDetachedOracleSupportAncestor(
                oracle,
                pixelId,
                targetSnapshot->nodeId,
                targetSnapshot->altitude,
                tree.isMaxtree());
            if (detachedAncestorId == InvalidNode ||
                detachedAncestorId == targetSnapshot->nodeId) {
                continue;
            }

            const NodeId branchRootId = findMaterializedOrExistingOracleNode(
                tree,
                oracle,
                plan,
                createdNodeIds,
                bridgeNodes,
                detachedAncestorId);
            if (branchRootId == InvalidNode ||
                branchRootId == targetNodeId ||
                !tree.isAlive(branchRootId) ||
                tree.isRoot(branchRootId) ||
                std::find(movedBranchRoots.begin(),
                          movedBranchRoots.end(),
                          branchRootId) != movedBranchRoots.end()) {
                continue;
            }

            const bool targetAlreadyContainsPixel = supportContainsPixel(tree, targetNodeId, pixelId);
            const bool branchIsDirectVisibleChild = tree.getNodeParent(branchRootId) == targetNodeId;
            if (targetAlreadyContainsPixel && !branchIsDirectVisibleChild) {
                continue;
            }

            const int connectorLevel = hiddenConnectorLevelBetween(
                tree.getAltitude(targetNodeId),
                tree.getAltitude(branchRootId),
                tree.isMaxtree());
            if (connectorLevel < 0) {
                continue;
            }

            const NodeId connectorId = tree.createDetachedNode(connectorLevel);
            if (connectorId == InvalidNode) {
                throw std::runtime_error("local recomposition could not create a hidden support connector");
            }
            tree.attachNode(targetNodeId, connectorId);
            const NodeId previousParent = tree.getNodeParent(branchRootId);
            tree.moveNode(branchRootId, connectorId);
            recordEmptyPreviousParent(
                tree,
                previousParent,
                connectorId,
                emptiedSourceOwners);
            canonicalFocusNodeIds.push_back(connectorId);
            canonicalFocusNodeIds.push_back(branchRootId);
            movedBranchRoots.push_back(branchRootId);

            LocalCreatedRecompositionNode created;
            created.planIndex = -1;
            created.nodeId = connectorId;
            created.parentNodeId = targetNodeId;
            created.level = static_cast<AltitudeType>(connectorLevel);
            created.kind = LocalCreatedRecompositionNodeKind::HiddenSupportConnector;
            created.movedChildNodes = {branchRootId};
            result.createdNodes.push_back(std::move(created));
        }
    }
}

inline bool preSnapshotIsNonOracleHiddenSupportCandidate(
    const DynamicComponentTree &tree,
    const LocalDualNodeSnapshot &snapshot,
    const LocalRecompositionNodePlan &nodePlan,
    NodeId targetNodeId) {
    if (snapshot.nodeId == InvalidNode ||
        snapshot.nodeId == targetNodeId ||
        snapshot.parentInSnapshot ||
        snapshot.supportLocalPixels.empty() ||
        snapshot.supportTouchesModified ||
        !snapshot.supportTouchesBoundary ||
        !tree.isAlive(snapshot.nodeId) ||
        tree.isRoot(snapshot.nodeId)) {
        return false;
    }

    auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
    auto targetSupport = sortedUniquePixels(nodePlan.supportPixels);
    if (targetSupport.empty() ||
        equalPixelSets(std::span<const PixelId>(snapshotSupport),
                       std::span<const PixelId>(targetSupport)) ||
        !containsAllPixels(std::span<const PixelId>(targetSupport),
                           std::span<const PixelId>(snapshotSupport))) {
        return false;
    }

    return hiddenConnectorLevelBetween(
               tree.getAltitude(targetNodeId),
               tree.getAltitude(snapshot.nodeId),
               tree.isMaxtree()) >= 0;
}

inline void materializeNonOracleHiddenSupportConnectors(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result) {
    std::vector<NodeId> movedBranchRoots;
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size() && planIndex < createdNodeIds.size(); ++planIndex) {
        const auto &nodePlan = plan.nodes[planIndex];
        if (nodePlan.oracleNodeId != InvalidNode) {
            continue;
        }
        const NodeId targetNodeId = createdNodeIds[planIndex];
        if (targetNodeId == InvalidNode || !tree.isAlive(targetNodeId)) {
            continue;
        }

        for (const auto &snapshot : oracle.preDualNodes) {
            const NodeId branchRootId = snapshot.nodeId;
            if (!preSnapshotIsNonOracleHiddenSupportCandidate(
                    tree,
                    snapshot,
                    nodePlan,
                    targetNodeId) ||
                std::find(movedBranchRoots.begin(),
                          movedBranchRoots.end(),
                          branchRootId) != movedBranchRoots.end()) {
                continue;
            }
            ++result.nonOracleHiddenSupportCandidateCount;

            bool targetAlreadyContainsBranchSupport = true;
            for (PixelId pixelId : snapshot.supportLocalPixels) {
                if (!supportContainsPixel(tree, targetNodeId, pixelId)) {
                    targetAlreadyContainsBranchSupport = false;
                    break;
                }
            }
            const bool branchIsDirectVisibleChild = tree.getNodeParent(branchRootId) == targetNodeId;
            if (targetAlreadyContainsBranchSupport && !branchIsDirectVisibleChild) {
                ++result.nonOracleHiddenSupportAlreadyRepresentedCount;
                continue;
            }

            const int connectorLevel = hiddenConnectorLevelBetween(
                tree.getAltitude(targetNodeId),
                tree.getAltitude(branchRootId),
                tree.isMaxtree());
            if (connectorLevel < 0) {
                continue;
            }

            ++result.nonOracleHiddenSupportMaterializableCount;
            const NodeId connectorId = tree.createDetachedNode(connectorLevel);
            if (connectorId == InvalidNode) {
                throw std::runtime_error("local recomposition could not create a non-oracle hidden support connector");
            }
            tree.attachNode(targetNodeId, connectorId);
            const NodeId previousParent = tree.getNodeParent(branchRootId);
            tree.moveNode(branchRootId, connectorId);
            recordEmptyPreviousParent(
                tree,
                previousParent,
                connectorId,
                emptiedSourceOwners);
            canonicalFocusNodeIds.push_back(connectorId);
            canonicalFocusNodeIds.push_back(branchRootId);
            movedBranchRoots.push_back(branchRootId);

            LocalCreatedRecompositionNode created;
            created.planIndex = -1;
            created.nodeId = connectorId;
            created.parentNodeId = targetNodeId;
            created.level = static_cast<AltitudeType>(connectorLevel);
            created.kind = LocalCreatedRecompositionNodeKind::NonOracleHiddenSupportConnector;
            created.movedChildNodes = {branchRootId};
            result.createdNodes.push_back(std::move(created));
        }
    }
}

inline NodeId findNonOracleMaterializedNodeForPreSnapshot(
    const DynamicComponentTree &tree,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    const LocalDualNodeSnapshot &snapshot,
    const std::vector<std::pair<NodeId, NodeId>> &preBridgeNodes = {}) {
    const NodeId bridgeNodeId = findMaterializedBridgeNode(preBridgeNodes, snapshot.nodeId);
    if (bridgeNodeId != InvalidNode && tree.isAlive(bridgeNodeId)) {
        return bridgeNodeId;
    }

    auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size() && planIndex < createdNodeIds.size(); ++planIndex) {
        const auto &nodePlan = plan.nodes[planIndex];
        if (nodePlan.oracleNodeId != InvalidNode ||
            nodePlan.level != static_cast<AltitudeType>(snapshot.altitude)) {
            continue;
        }
        auto nodeSupport = sortedUniquePixels(nodePlan.supportPixels);
        if (equalPixelSets(std::span<const PixelId>(nodeSupport),
                           std::span<const PixelId>(snapshotSupport)) &&
            createdNodeIds[planIndex] != InvalidNode &&
            tree.isAlive(createdNodeIds[planIndex])) {
            return createdNodeIds[planIndex];
        }
    }

    if (snapshot.nodeId != InvalidNode &&
        tree.isNode(snapshot.nodeId) &&
        tree.isAlive(snapshot.nodeId) &&
        tree.getAltitude(snapshot.nodeId) == snapshot.altitude) {
        return snapshot.nodeId;
    }
    return InvalidNode;
}

inline bool preSnapshotParentLinkTouchesNonOraclePlan(
    const LocalRecompositionPlan &plan,
    const LocalDualNodeSnapshot &childSnapshot,
    const LocalDualNodeSnapshot &parentSnapshot) {
    auto childSupport = sortedUniquePixels(childSnapshot.supportLocalPixels);
    auto parentSupport = sortedUniquePixels(parentSnapshot.supportLocalPixels);
    for (const auto &nodePlan : plan.nodes) {
        if (nodePlan.oracleNodeId != InvalidNode) {
            continue;
        }
        if (containsAllPixels(std::span<const PixelId>(nodePlan.supportPixels),
                              std::span<const PixelId>(childSupport)) &&
            containsAllPixels(std::span<const PixelId>(nodePlan.supportPixels),
                              std::span<const PixelId>(parentSupport))) {
            return true;
        }
    }
    return false;
}

inline bool currentParentIsLocalIntermediateBetween(const DynamicComponentTree &tree,
                                                    NodeId currentParentId,
                                                    NodeId childNodeId,
                                                    NodeId expectedParentId) {
    if (currentParentId == InvalidNode ||
        currentParentId == expectedParentId ||
        childNodeId == InvalidNode ||
        expectedParentId == InvalidNode ||
        !tree.isAlive(currentParentId) ||
        !tree.isAlive(childNodeId) ||
        !tree.isAlive(expectedParentId) ||
        tree.getNumProperParts(currentParentId) != 0 ||
        tree.getNodeParent(currentParentId) != expectedParentId) {
        return false;
    }
    return validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(childNodeId)),
                                 static_cast<AltitudeType>(tree.getAltitude(currentParentId)),
                                 tree.isMaxtree()) &&
           validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(currentParentId)),
                                 static_cast<AltitudeType>(tree.getAltitude(expectedParentId)),
                                 tree.isMaxtree());
}

inline NodeId materializeMissingNonOraclePreParentBridgeAboveChild(
    DynamicComponentTree &tree,
    const LocalDualNodeSnapshot &parentSnapshot,
    const LocalDualNodeSnapshot &childSnapshot,
    NodeId childNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &preBridgeNodes) {
    ++result.nonOracleParentLinkBridgeAttemptCount;
    if (childNodeId == InvalidNode ||
        !tree.isAlive(childNodeId) ||
        tree.isRoot(childNodeId) ||
        !parentSnapshot.directLocalProperParts.empty()) {
        return InvalidNode;
    }

    auto parentSupport = sortedUniquePixels(parentSnapshot.supportLocalPixels);
    auto childSupport = sortedUniquePixels(childSnapshot.supportLocalPixels);
    if (parentSupport.empty() ||
        childSupport.empty() ||
        !containsAllPixels(std::span<const PixelId>(parentSupport),
                           std::span<const PixelId>(childSupport))) {
        return InvalidNode;
    }

    const NodeId currentParentId = tree.getNodeParent(childNodeId);
    if (currentParentId == InvalidNode ||
        currentParentId == childNodeId ||
        !tree.isAlive(currentParentId)) {
        return InvalidNode;
    }

    const auto bridgeLevel = static_cast<AltitudeType>(parentSnapshot.altitude);
    if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(childNodeId)),
                               bridgeLevel,
                               tree.isMaxtree()) ||
        !validLocalParentLevel(bridgeLevel,
                               static_cast<AltitudeType>(tree.getAltitude(currentParentId)),
                               tree.isMaxtree())) {
        ++result.nonOracleParentLinkInvalidLevelCount;
        return InvalidNode;
    }

    ++result.nonOracleParentLinkBridgeMaterializableCount;
    const NodeId bridgeNodeId = tree.createDetachedNode(static_cast<int>(bridgeLevel));
    if (bridgeNodeId == InvalidNode) {
        throw std::runtime_error("local recomposition could not create a non-oracle parent-link bridge");
    }
    tree.attachNode(currentParentId, bridgeNodeId);
    tree.moveNode(childNodeId, bridgeNodeId);
    recordEmptyPreviousParent(
        tree,
        currentParentId,
        bridgeNodeId,
        emptiedSourceOwners);
    canonicalFocusNodeIds.push_back(bridgeNodeId);
    canonicalFocusNodeIds.push_back(childNodeId);

    LocalCreatedRecompositionNode created;
    created.planIndex = -1;
    created.nodeId = bridgeNodeId;
    created.parentNodeId = currentParentId;
    created.level = bridgeLevel;
    created.kind = LocalCreatedRecompositionNodeKind::NonOracleParentLinkBridge;
    created.movedChildNodes = {childNodeId};
    result.createdNodes.push_back(std::move(created));
    preBridgeNodes.push_back({parentSnapshot.nodeId, bridgeNodeId});
    return bridgeNodeId;
}

inline NodeId materializeMissingNonOraclePreChildBridgeBelowParent(
    DynamicComponentTree &tree,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    std::span<const LocalDualNodeSnapshot> preSnapshots,
    const LocalDualNodeSnapshot &childSnapshot,
    NodeId parentNodeId,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &preBridgeNodes) {
    ++result.nonOracleChildLinkBridgeAttemptCount;
    if (parentNodeId == InvalidNode ||
        !tree.isAlive(parentNodeId) ||
        !childSnapshot.directLocalProperParts.empty() ||
        childSnapshot.supportTouchesModified ||
        childSnapshot.childrenInSnapshot.empty()) {
        return InvalidNode;
    }

    auto childSupport = sortedUniquePixels(childSnapshot.supportLocalPixels);
    if (childSupport.empty()) {
        return InvalidNode;
    }

    const auto bridgeLevel = static_cast<AltitudeType>(childSnapshot.altitude);
    if (!validLocalParentLevel(bridgeLevel,
                               static_cast<AltitudeType>(tree.getAltitude(parentNodeId)),
                               tree.isMaxtree())) {
        ++result.nonOracleChildLinkBridgeInvalidLevelCount;
        return InvalidNode;
    }

    std::vector<NodeId> liveChildrenToMove;
    for (NodeId descendantSnapshotId : childSnapshot.childrenInSnapshot) {
        const auto *descendantSnapshot = findLocalSnapshotByNodeId(
            preSnapshots,
            descendantSnapshotId);
        if (descendantSnapshot == nullptr) {
            continue;
        }

        auto descendantSupport = sortedUniquePixels(descendantSnapshot->supportLocalPixels);
        if (descendantSupport.empty() ||
            !containsAllPixels(std::span<const PixelId>(childSupport),
                               std::span<const PixelId>(descendantSupport))) {
            continue;
        }

        const NodeId descendantNodeId = findNonOracleMaterializedNodeForPreSnapshot(
            tree,
            plan,
            createdNodeIds,
            *descendantSnapshot,
            preBridgeNodes);
        if (descendantNodeId == InvalidNode ||
            descendantNodeId == parentNodeId ||
            !tree.isAlive(descendantNodeId) ||
            tree.isRoot(descendantNodeId) ||
            !validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(descendantNodeId)),
                                   bridgeLevel,
                                   tree.isMaxtree())) {
            continue;
        }
        liveChildrenToMove.push_back(descendantNodeId);
    }
    liveChildrenToMove = sortedUniqueValidNodes(std::move(liveChildrenToMove));
    if (liveChildrenToMove.empty()) {
        return InvalidNode;
    }

    ++result.nonOracleChildLinkBridgeMaterializableCount;
    const NodeId bridgeNodeId = tree.createDetachedNode(static_cast<int>(bridgeLevel));
    if (bridgeNodeId == InvalidNode) {
        throw std::runtime_error("local recomposition could not create a non-oracle child-link bridge");
    }
    tree.attachNode(parentNodeId, bridgeNodeId);
    canonicalFocusNodeIds.push_back(bridgeNodeId);

    LocalCreatedRecompositionNode created;
    created.planIndex = -1;
    created.nodeId = bridgeNodeId;
    created.parentNodeId = parentNodeId;
    created.level = bridgeLevel;
    created.kind = LocalCreatedRecompositionNodeKind::NonOracleChildLinkBridge;

    for (NodeId liveChildNodeId : liveChildrenToMove) {
        const NodeId previousParentId = tree.getNodeParent(liveChildNodeId);
        if (previousParentId != bridgeNodeId) {
            tree.moveNode(liveChildNodeId, bridgeNodeId);
            recordEmptyPreviousParent(
                tree,
                previousParentId,
                bridgeNodeId,
                emptiedSourceOwners);
        }
        canonicalFocusNodeIds.push_back(liveChildNodeId);
        created.movedChildNodes.push_back(liveChildNodeId);
    }

    created.movedChildNodes =
        sortedUniqueValidNodes(std::move(created.movedChildNodes));
    result.createdNodes.push_back(std::move(created));
    preBridgeNodes.push_back({childSnapshot.nodeId, bridgeNodeId});
    return bridgeNodeId;
}

inline void enforceNonOraclePreSnapshotParentLinks(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result) {
    std::vector<const LocalDualNodeSnapshot *> snapshots;
    snapshots.reserve(oracle.preDualNodes.size());
    const std::span<const LocalDualNodeSnapshot> preSnapshots(oracle.preDualNodes);
    for (const auto &snapshot : oracle.preDualNodes) {
        if (snapshot.parentInSnapshot &&
            snapshot.parentId != InvalidNode &&
            snapshot.parentId != snapshot.nodeId) {
            snapshots.push_back(&snapshot);
        }
    }

    std::sort(snapshots.begin(), snapshots.end(), [&](const auto *lhs, const auto *rhs) {
        if (lhs->altitude != rhs->altitude) {
            return tree.isMaxtree() ? lhs->altitude < rhs->altitude
                                    : lhs->altitude > rhs->altitude;
        }
        return lhs->nodeId < rhs->nodeId;
    });

    std::vector<std::pair<NodeId, NodeId>> preBridgeNodes;
    for (const auto *childSnapshot : snapshots) {
        const auto *parentSnapshot = findLocalSnapshotByNodeId(
            preSnapshots,
            childSnapshot->parentId);
        if (parentSnapshot == nullptr ||
            !preSnapshotParentLinkTouchesNonOraclePlan(plan, *childSnapshot, *parentSnapshot)) {
            continue;
        }
        ++result.nonOracleParentLinkTouchesPlanCount;

        const NodeId parentNodeId = findNonOracleMaterializedNodeForPreSnapshot(
            tree,
            plan,
            createdNodeIds,
            *parentSnapshot,
            preBridgeNodes);
        const NodeId childNodeId = findNonOracleMaterializedNodeForPreSnapshot(
            tree,
            plan,
            createdNodeIds,
            *childSnapshot,
            preBridgeNodes);
        if (parentNodeId == InvalidNode) {
            if (childNodeId != InvalidNode &&
                !parentSnapshot->supportTouchesModified &&
                !(childSnapshot->supportTouchesModified && childNodeId == childSnapshot->nodeId)) {
                (void) materializeMissingNonOraclePreParentBridgeAboveChild(
                    tree,
                    *parentSnapshot,
                    *childSnapshot,
                    childNodeId,
                    canonicalFocusNodeIds,
                    emptiedSourceOwners,
                    result,
                    preBridgeNodes);
            }
            continue;
        }

        if (childNodeId == InvalidNode) {
            (void) materializeMissingNonOraclePreChildBridgeBelowParent(
                tree,
                plan,
                createdNodeIds,
                preSnapshots,
                *childSnapshot,
                parentNodeId,
                canonicalFocusNodeIds,
                emptiedSourceOwners,
                result,
                preBridgeNodes);
            continue;
        }

        if (parentNodeId == childNodeId ||
            !tree.isAlive(parentNodeId) ||
            !tree.isAlive(childNodeId) ||
            tree.isRoot(childNodeId)) {
            continue;
        }
        ++result.nonOracleParentLinkBothRepresentedCount;
        if ((parentSnapshot->supportTouchesModified && parentNodeId == parentSnapshot->nodeId) ||
            (childSnapshot->supportTouchesModified && childNodeId == childSnapshot->nodeId)) {
            ++result.nonOracleParentLinkBlockedByModifiedSnapshotCount;
            continue;
        }

        if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(childNodeId)),
                                   static_cast<AltitudeType>(tree.getAltitude(parentNodeId)),
                                   tree.isMaxtree())) {
            ++result.nonOracleParentLinkInvalidLevelCount;
            continue;
        }

        const NodeId currentParentId = tree.getNodeParent(childNodeId);
        if (currentParentId == parentNodeId) {
            ++result.nonOracleParentLinkAlreadyDirectCount;
            canonicalFocusNodeIds.push_back(childNodeId);
            continue;
        }
        if (currentParentIsLocalIntermediateBetween(
                tree,
                currentParentId,
                childNodeId,
                parentNodeId)) {
            ++result.nonOracleParentLinkAlreadyIntermediateCount;
            canonicalFocusNodeIds.push_back(currentParentId);
            canonicalFocusNodeIds.push_back(childNodeId);
            continue;
        }

        ++result.nonOracleParentLinkMovedCount;
        tree.moveNode(childNodeId, parentNodeId);
        canonicalFocusNodeIds.push_back(childNodeId);
        recordEmptyPreviousParent(
            tree,
            currentParentId,
            parentNodeId,
            emptiedSourceOwners);
    }
}

inline void enforceOracleSnapshotParentLinks(
    DynamicComponentTree &tree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    std::span<const NodeId> createdNodeIds,
    std::vector<NodeId> &canonicalFocusNodeIds,
    std::vector<NodeId> &emptiedSourceOwners,
    LocalRecompositionMaterializationResult &result,
    std::vector<std::pair<NodeId, NodeId>> &bridgeNodes) {
    std::vector<const LocalDualNodeSnapshot *> snapshots;
    snapshots.reserve(oracle.oracleDualNodes.size());
    for (const auto &snapshot : oracle.oracleDualNodes) {
        if (snapshot.parentInSnapshot &&
            snapshot.parentId != InvalidNode &&
            snapshot.parentId != snapshot.nodeId) {
            snapshots.push_back(&snapshot);
        }
    }

    std::sort(snapshots.begin(), snapshots.end(), [&](const auto *lhs, const auto *rhs) {
        if (lhs->altitude != rhs->altitude) {
            return tree.isMaxtree() ? lhs->altitude < rhs->altitude
                                    : lhs->altitude > rhs->altitude;
        }
        return lhs->nodeId < rhs->nodeId;
    });

    bool changed = true;
    for (std::size_t pass = 0; changed && pass <= snapshots.size(); ++pass) {
        changed = false;
        for (const auto *snapshot : snapshots) {
            const NodeId parentNodeId = findMaterializedOrExistingOracleNode(
                tree,
                oracle,
                plan,
                createdNodeIds,
                bridgeNodes,
                snapshot->parentId);
            if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
                const NodeId childNodeId = findMaterializedOrExistingOracleNode(
                    tree,
                    oracle,
                    plan,
                    createdNodeIds,
                    bridgeNodes,
                    snapshot->nodeId);
                if (childNodeId != InvalidNode &&
                    materializeMissingOracleParentBridgeAboveChild(
                        tree,
                        oracle,
                        plan,
                        snapshot->parentId,
                        childNodeId,
                        canonicalFocusNodeIds,
                        emptiedSourceOwners,
                        result,
                        bridgeNodes) != InvalidNode) {
                    changed = true;
                }
                continue;
            }

            NodeId childNodeId = findMaterializedOrExistingOracleNode(
                tree,
                oracle,
                plan,
                createdNodeIds,
                bridgeNodes,
                snapshot->nodeId);
            if (childNodeId == InvalidNode) {
                childNodeId = materializeOracleBridgeSubtree(
                    tree,
                    oracle,
                    plan,
                    snapshot->nodeId,
                    parentNodeId,
                    canonicalFocusNodeIds,
                    emptiedSourceOwners,
                    result,
                    bridgeNodes);
                if (childNodeId != InvalidNode) {
                    changed = true;
                }
                continue;
            }

            if (childNodeId == parentNodeId ||
                !tree.isAlive(childNodeId) ||
                tree.isRoot(childNodeId)) {
                continue;
            }
            if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(childNodeId)),
                                       static_cast<AltitudeType>(tree.getAltitude(parentNodeId)),
                                       tree.isMaxtree())) {
                continue;
            }
            if (tree.getNodeParent(childNodeId) == parentNodeId) {
                canonicalFocusNodeIds.push_back(childNodeId);
                continue;
            }

            const NodeId previousParent = tree.getNodeParent(childNodeId);
            tree.moveNode(childNodeId, parentNodeId);
            canonicalFocusNodeIds.push_back(childNodeId);
            recordEmptyPreviousParent(
                tree,
                previousParent,
                parentNodeId,
                emptiedSourceOwners);
            changed = true;
        }
    }
}

inline void releaseEmptySourceOwners(DynamicComponentTree &tree,
                                     std::vector<NodeId> emptiedSourceOwners,
                                     std::vector<NodeId> &releasedSourceOwners,
                                     bool releaseOwnersWithResidualChildren = false,
                                     std::span<const NodeId> autoPromotableResidualChildren = {}) {
    emptiedSourceOwners = sortedUniqueValidNodes(std::move(emptiedSourceOwners));
    bool releasedAny = true;
    while (releasedAny) {
        releasedAny = false;
        for (NodeId sourceOwner : emptiedSourceOwners) {
            if (sourceOwner == InvalidNode ||
                !tree.isAlive(sourceOwner) ||
                tree.isRoot(sourceOwner) ||
                tree.getNumProperParts(sourceOwner) != 0) {
                continue;
            }
            const NodeId parentId = tree.getNodeParent(sourceOwner);
            if (parentId == InvalidNode || parentId == sourceOwner || !tree.isAlive(parentId)) {
                continue;
            }
            if (tree.getNumChildren(sourceOwner) != 0) {
                if (!releaseOwnersWithResidualChildren) {
                    bool allChildrenAutoPromotable = true;
                    for (NodeId childId : tree.getChildren(sourceOwner)) {
                        if (std::find(autoPromotableResidualChildren.begin(),
                                      autoPromotableResidualChildren.end(),
                                      childId) == autoPromotableResidualChildren.end()) {
                            allChildrenAutoPromotable = false;
                            break;
                        }
                    }
                    if (!allChildrenAutoPromotable) {
                        continue;
                    }
                }
                tree.moveChildren(parentId, sourceOwner);
            }
            tree.removeChild(parentId, sourceOwner, true);
            releasedSourceOwners.push_back(sourceOwner);
            releasedAny = true;
        }
    }
    releasedSourceOwners = sortedUniqueValidNodes(std::move(releasedSourceOwners));
}

inline bool containsNode(std::span<const NodeId> nodes, NodeId nodeId) {
    return std::find(nodes.begin(), nodes.end(), nodeId) != nodes.end();
}

inline int releaseAllEmptyNonRootNodes(DynamicComponentTree &tree,
                                       std::span<const NodeId> protectedNodes) {
    int released = 0;
    bool releasedAny = true;
    while (releasedAny) {
        releasedAny = false;
        const int nodeSlots = tree.getNumInternalNodeSlots();
        for (NodeId nodeId = 0; nodeId < nodeSlots; ++nodeId) {
            if (nodeId == InvalidNode ||
                !tree.isNode(nodeId) ||
                !tree.isAlive(nodeId) ||
                tree.isRoot(nodeId) ||
                containsNode(protectedNodes, nodeId) ||
                tree.getNumProperParts(nodeId) != 0) {
                continue;
            }
            const NodeId parentId = tree.getNodeParent(nodeId);
            if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
                continue;
            }
            tree.moveChildren(parentId, nodeId);
            tree.removeChild(parentId, nodeId, true);
            ++released;
            releasedAny = true;
        }
    }
    return released;
}

inline void addLocalCanonicalizationCandidate(DynamicComponentTree &tree,
                                              std::vector<NodeId> &candidates,
                                              NodeId nodeId) {
    if (nodeId == InvalidNode ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId)) {
        return;
    }
    candidates.push_back(nodeId);
    if (!tree.isRoot(nodeId)) {
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId != InvalidNode && parentId != nodeId && tree.isAlive(parentId)) {
            candidates.push_back(parentId);
            if (!tree.isRoot(parentId)) {
                const NodeId grandParentId = tree.getNodeParent(parentId);
                if (grandParentId != InvalidNode &&
                    grandParentId != parentId &&
                    tree.isAlive(grandParentId)) {
                    candidates.push_back(grandParentId);
                }
            }
        }
    }
    for (NodeId childId : tree.getChildren(nodeId)) {
        if (childId != InvalidNode && tree.isAlive(childId)) {
            candidates.push_back(childId);
        }
    }
}

inline std::vector<NodeId> collectLocalCanonicalizationCandidates(
    DynamicComponentTree &tree,
    std::span<const NodeId> focusNodes) {
    std::vector<NodeId> candidates;
    for (NodeId focusNodeId : focusNodes) {
        addLocalCanonicalizationCandidate(tree, candidates, focusNodeId);
    }
    return sortedUniqueValidNodes(std::move(candidates));
}

inline void enqueueFocusedEmptyReleaseCandidate(DynamicComponentTree &tree,
                                                NodeId nodeId,
                                                std::vector<NodeId> &worklist,
                                                GenerationStampSet &queuedMarks) {
    if (nodeId == InvalidNode ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId)) {
        return;
    }
    const auto index = static_cast<std::size_t>(nodeId);
    if (queuedMarks.isMarked(index)) {
        return;
    }
    queuedMarks.mark(index);
    worklist.push_back(nodeId);
}

inline void enqueueFocusedEmptyReleaseNeighborhood(DynamicComponentTree &tree,
                                                   NodeId nodeId,
                                                   std::vector<NodeId> &worklist,
                                                   GenerationStampSet &queuedMarks) {
    if (nodeId == InvalidNode ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId)) {
        return;
    }

    enqueueFocusedEmptyReleaseCandidate(tree, nodeId, worklist, queuedMarks);
    if (!tree.isRoot(nodeId)) {
        const NodeId parentId = tree.getNodeParent(nodeId);
        enqueueFocusedEmptyReleaseCandidate(tree, parentId, worklist, queuedMarks);
        if (parentId != InvalidNode &&
            parentId != nodeId &&
            tree.isNode(parentId) &&
            tree.isAlive(parentId) &&
            !tree.isRoot(parentId)) {
            enqueueFocusedEmptyReleaseCandidate(
                tree,
                tree.getNodeParent(parentId),
                worklist,
                queuedMarks);
        }
    }
    for (NodeId childId : tree.getChildren(nodeId)) {
        enqueueFocusedEmptyReleaseCandidate(tree, childId, worklist, queuedMarks);
    }
}

inline int releaseFocusedEmptyNonRootNodes(DynamicComponentTree &tree,
                                           std::span<const NodeId> focusNodes,
                                           std::span<const NodeId> protectedNodes,
                                           GenerationStampSet *queuedMarksScratch = nullptr,
                                           std::vector<NodeId> *worklistScratch = nullptr,
                                           std::vector<NodeId> *dirtySeeds = nullptr,
                                           std::vector<NodeId> *releasedNodes = nullptr) {
    int released = 0;
    GenerationStampSet localQueuedMarks;
    std::vector<NodeId> localWorklist;
    auto &queuedMarks =
        queuedMarksScratch == nullptr ? localQueuedMarks : *queuedMarksScratch;
    auto &worklist =
        worklistScratch == nullptr ? localWorklist : *worklistScratch;

    queuedMarks.ensureSize(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
    queuedMarks.resetAll();
    worklist.clear();

    for (NodeId focusNodeId : focusNodes) {
        enqueueFocusedEmptyReleaseNeighborhood(tree, focusNodeId, worklist, queuedMarks);
    }

    std::size_t head = 0;
    while (head < worklist.size()) {
        const NodeId nodeId = worklist[head++];
        if (nodeId == InvalidNode ||
            !tree.isNode(nodeId) ||
            !tree.isAlive(nodeId) ||
            tree.isRoot(nodeId) ||
            containsNode(protectedNodes, nodeId) ||
            tree.getNumProperParts(nodeId) != 0) {
            continue;
        }

        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode ||
            parentId == nodeId ||
            !tree.isNode(parentId) ||
            !tree.isAlive(parentId)) {
            continue;
        }

        if (!tree.isRoot(parentId)) {
            enqueueFocusedEmptyReleaseCandidate(
                tree,
                tree.getNodeParent(parentId),
                worklist,
                queuedMarks);
        }
        for (NodeId childId : tree.getChildren(nodeId)) {
            enqueueFocusedEmptyReleaseCandidate(tree, childId, worklist, queuedMarks);
        }
        enqueueFocusedEmptyReleaseCandidate(tree, parentId, worklist, queuedMarks);

        if (dirtySeeds != nullptr) {
            dirtySeeds->push_back(parentId);
        }
        if (releasedNodes != nullptr) {
            releasedNodes->push_back(nodeId);
        }
        tree.moveChildren(parentId, nodeId);
        tree.removeChild(parentId, nodeId, true);
        ++released;

        enqueueFocusedEmptyReleaseNeighborhood(tree, parentId, worklist, queuedMarks);
    }
    return released;
}

inline int releaseAllEmptyNonRootNodes(DynamicComponentTree &tree) {
    return releaseAllEmptyNonRootNodes(tree, {});
}

inline int releaseAllEmptyLeafNodes(DynamicComponentTree &tree) {
    return releaseAllEmptyNonRootNodes(tree);
}

inline bool promoteEmptyRootWithSingleChild(DynamicComponentTree &tree) {
    const NodeId oldRoot = tree.getRoot();
    if (oldRoot == InvalidNode ||
        !tree.isNode(oldRoot) ||
        !tree.isAlive(oldRoot) ||
        tree.getNumProperParts(oldRoot) != 0 ||
        tree.getNumChildren(oldRoot) != 1) {
        return false;
    }
    NodeId newRoot = InvalidNode;
    for (NodeId childId : tree.getChildren(oldRoot)) {
        newRoot = childId;
        break;
    }
    if (newRoot == InvalidNode || !tree.isNode(newRoot) || !tree.isAlive(newRoot)) {
        return false;
    }
    tree.setRoot(newRoot);
    if (tree.isNode(oldRoot) &&
        tree.isAlive(oldRoot) &&
        tree.getNumProperParts(oldRoot) == 0 &&
        tree.getNumChildren(oldRoot) == 0) {
        tree.releaseNode(oldRoot);
    }
    return true;
}

inline void collapseNodeIntoParent(DynamicComponentTree &tree, NodeId nodeId) {
    if (nodeId == InvalidNode ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId) ||
        tree.isRoot(nodeId)) {
        return;
    }
    const NodeId parentId = tree.getNodeParent(nodeId);
    if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
        return;
    }
    tree.moveChildren(parentId, nodeId);
    tree.moveProperParts(parentId, nodeId);
    tree.removeChild(parentId, nodeId, true);
}

inline bool collapseSameLevelParentChildNodes(DynamicComponentTree &tree) {
    const int nodeSlots = tree.getNumInternalNodeSlots();
    for (NodeId nodeId = 0; nodeId < nodeSlots; ++nodeId) {
        if (nodeId == InvalidNode ||
            !tree.isNode(nodeId) ||
            !tree.isAlive(nodeId) ||
            tree.isRoot(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
            continue;
        }
        if (tree.getAltitude(nodeId) == tree.getAltitude(parentId)) {
            collapseNodeIntoParent(tree, nodeId);
            return true;
        }
    }
    return false;
}

inline bool collapseFocusedSameLevelParentChildNodes(DynamicComponentTree &tree,
                                                     std::span<const NodeId> focusNodes) {
    auto candidates = collectLocalCanonicalizationCandidates(tree, focusNodes);
    for (NodeId nodeId : candidates) {
        if (nodeId == InvalidNode ||
            !tree.isNode(nodeId) ||
            !tree.isAlive(nodeId) ||
            tree.isRoot(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
            continue;
        }
        if (tree.getAltitude(nodeId) == tree.getAltitude(parentId)) {
            collapseNodeIntoParent(tree, nodeId);
            return true;
        }
    }
    return false;
}

inline bool supportsAreAdjacent(const DynamicComponentTree &tree,
                                AdjacencyRelation &adjacency,
                                NodeId lhs,
                                NodeId rhs) {
    const int domainSize = tree.getNumTotalProperParts();
    const auto lhsPixels = collectNodeSupportPixels(tree, lhs);
    if (lhsPixels.empty()) {
        return false;
    }

    auto nodeSupportContainsPixel = [&](NodeId nodeId, PixelId pixelId) {
        if (pixelId < 0 || pixelId >= domainSize) {
            return false;
        }
        NodeId current = tree.getSmallestComponent(pixelId);
        while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
            if (current == nodeId) {
                return true;
            }
            const NodeId parentId = tree.getNodeParent(current);
            if (parentId == current) {
                break;
            }
            current = parentId;
        }
        return false;
    };

    for (PixelId pixelId : lhsPixels) {
        if (nodeSupportContainsPixel(rhs, pixelId)) {
            return true;
        }
        for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
            if (nodeSupportContainsPixel(rhs, neighbor)) {
                return true;
            }
        }
    }
    return false;
}

inline bool nodeSupportContainsPixel(const DynamicComponentTree &tree,
                                     NodeId nodeId,
                                     PixelId pixelId) {
    if (nodeId == InvalidNode ||
        pixelId < 0 ||
        pixelId >= tree.getNumTotalProperParts() ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId)) {
        return false;
    }
    NodeId current = tree.getSmallestComponent(pixelId);
    while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
        if (current == nodeId) {
            return true;
        }
        const NodeId parentId = tree.getNodeParent(current);
        if (parentId == current) {
            break;
        }
        current = parentId;
    }
    return false;
}

inline bool validFocusedComparableSiblingMove(const DynamicComponentTree &tree,
                                              NodeId focusNodeId,
                                              NodeId siblingId) {
    if (focusNodeId == InvalidNode ||
        siblingId == InvalidNode ||
        !tree.isNode(focusNodeId) ||
        !tree.isNode(siblingId) ||
        !tree.isAlive(focusNodeId) ||
        !tree.isAlive(siblingId) ||
        tree.getAltitude(focusNodeId) == tree.getAltitude(siblingId)) {
        return false;
    }

    NodeId container = focusNodeId;
    NodeId nested = siblingId;
    if (tree.isMaxtree()) {
        if (tree.getAltitude(container) > tree.getAltitude(nested)) {
            std::swap(container, nested);
        }
    } else if (tree.getAltitude(container) < tree.getAltitude(nested)) {
        std::swap(container, nested);
    }

    return validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(nested)),
                                 static_cast<AltitudeType>(tree.getAltitude(container)),
                                 tree.isMaxtree());
}

inline std::vector<NodeId> collectActionableCanonicalSiblings(const DynamicComponentTree &tree,
                                                              NodeId focusNodeId,
                                                              NodeId parentId) {
    std::vector<NodeId> siblings;
    if (focusNodeId == InvalidNode ||
        parentId == InvalidNode ||
        !tree.isNode(focusNodeId) ||
        !tree.isAlive(focusNodeId) ||
        !tree.isAlive(parentId)) {
        return siblings;
    }

    for (NodeId siblingId : tree.getChildren(parentId)) {
        if (siblingId == InvalidNode ||
            siblingId == focusNodeId ||
            !tree.isAlive(siblingId) ||
            tree.getNodeParent(siblingId) != parentId) {
            continue;
        }
        if (tree.getAltitude(focusNodeId) == tree.getAltitude(siblingId) ||
            validFocusedComparableSiblingMove(tree, focusNodeId, siblingId)) {
            siblings.push_back(siblingId);
        }
    }
    return sortedUniqueValidNodes(std::move(siblings));
}

inline int countNodeSupportProperPartsUpTo(const DynamicComponentTree &tree,
                                           NodeId nodeId,
                                           int limit,
                                           std::vector<NodeId> *supportNodeIds = nullptr) {
    if (nodeId == InvalidNode ||
        !tree.isNode(nodeId) ||
        !tree.isAlive(nodeId)) {
        return 0;
    }

    int count = 0;
    for (NodeId supportNodeId : tree.getNodeSubtree(nodeId)) {
        if (supportNodeId == InvalidNode || !tree.isAlive(supportNodeId)) {
            continue;
        }
        if (supportNodeIds != nullptr) {
            supportNodeIds->push_back(supportNodeId);
        }
        count += tree.getNumProperParts(supportNodeId);
        if (limit >= 0 && count >= limit) {
            return count;
        }
    }
    return count;
}

struct LocalCanonicalizationScratch {
    std::vector<unsigned int> directChildOwnerEpoch;
    std::vector<NodeId> directChildByOwner;
    std::vector<unsigned int> focusSupportEpoch;
    std::vector<unsigned int> actionableSiblingEpoch;
    std::vector<int> actionableSiblingIndex;
    unsigned int currentDirectChildOwnerEpoch = 1;
    unsigned int currentFocusSupportEpoch = 1;
    unsigned int currentActionableSiblingEpoch = 1;

    void ensureNodeCapacity(const DynamicComponentTree &tree) {
        const auto nodeSlots = static_cast<std::size_t>(tree.getNumInternalNodeSlots());
        if (directChildOwnerEpoch.size() < nodeSlots) {
            directChildOwnerEpoch.resize(nodeSlots, 0);
            directChildByOwner.resize(nodeSlots, InvalidNode);
            focusSupportEpoch.resize(nodeSlots, 0);
            actionableSiblingEpoch.resize(nodeSlots, 0);
            actionableSiblingIndex.resize(nodeSlots, -1);
        }
    }

    void beginDirectChildOwnerCache(const DynamicComponentTree &tree) {
        ensureNodeCapacity(tree);
        ++currentDirectChildOwnerEpoch;
        if (currentDirectChildOwnerEpoch == 0) {
            std::fill(directChildOwnerEpoch.begin(), directChildOwnerEpoch.end(), 0);
            currentDirectChildOwnerEpoch = 1;
        }
    }

    bool cachedDirectChildForOwner(NodeId ownerId, NodeId &directChildId) const {
        if (ownerId == InvalidNode ||
            ownerId < 0 ||
            static_cast<std::size_t>(ownerId) >= directChildOwnerEpoch.size() ||
            directChildOwnerEpoch[static_cast<std::size_t>(ownerId)] != currentDirectChildOwnerEpoch) {
            return false;
        }
        directChildId = directChildByOwner[static_cast<std::size_t>(ownerId)];
        return true;
    }

    void cacheDirectChildForOwner(NodeId ownerId, NodeId directChildId) {
        if (ownerId == InvalidNode ||
            ownerId < 0 ||
            static_cast<std::size_t>(ownerId) >= directChildOwnerEpoch.size()) {
            return;
        }
        const auto index = static_cast<std::size_t>(ownerId);
        directChildOwnerEpoch[index] = currentDirectChildOwnerEpoch;
        directChildByOwner[index] = directChildId;
    }

    void markFocusSupportNodes(const DynamicComponentTree &tree,
                               std::span<const NodeId> nodeIds) {
        ensureNodeCapacity(tree);
        ++currentFocusSupportEpoch;
        if (currentFocusSupportEpoch == 0) {
            std::fill(focusSupportEpoch.begin(), focusSupportEpoch.end(), 0);
            currentFocusSupportEpoch = 1;
        }
        for (NodeId nodeId : nodeIds) {
            if (nodeId != InvalidNode &&
                nodeId >= 0 &&
                static_cast<std::size_t>(nodeId) < focusSupportEpoch.size()) {
                focusSupportEpoch[static_cast<std::size_t>(nodeId)] =
                    currentFocusSupportEpoch;
            }
        }
    }

    bool containsFocusSupportOwner(NodeId ownerId) const {
        return ownerId != InvalidNode &&
               ownerId >= 0 &&
               static_cast<std::size_t>(ownerId) < focusSupportEpoch.size() &&
               focusSupportEpoch[static_cast<std::size_t>(ownerId)] == currentFocusSupportEpoch;
    }

    void markActionableSiblings(const DynamicComponentTree &tree,
                                std::span<const NodeId> siblings) {
        ensureNodeCapacity(tree);
        ++currentActionableSiblingEpoch;
        if (currentActionableSiblingEpoch == 0) {
            std::fill(actionableSiblingEpoch.begin(), actionableSiblingEpoch.end(), 0);
            currentActionableSiblingEpoch = 1;
        }
        for (std::size_t index = 0; index < siblings.size(); ++index) {
            const NodeId nodeId = siblings[index];
            if (nodeId == InvalidNode ||
                nodeId < 0 ||
                static_cast<std::size_t>(nodeId) >= actionableSiblingEpoch.size()) {
                continue;
            }
            const auto nodeIndex = static_cast<std::size_t>(nodeId);
            actionableSiblingEpoch[nodeIndex] = currentActionableSiblingEpoch;
            actionableSiblingIndex[nodeIndex] = static_cast<int>(index);
        }
    }

    bool actionableSiblingLookup(NodeId nodeId, std::size_t &index) const {
        if (nodeId == InvalidNode ||
            nodeId < 0 ||
            static_cast<std::size_t>(nodeId) >= actionableSiblingEpoch.size() ||
            actionableSiblingEpoch[static_cast<std::size_t>(nodeId)] != currentActionableSiblingEpoch) {
            return false;
        }
        const int storedIndex = actionableSiblingIndex[static_cast<std::size_t>(nodeId)];
        if (storedIndex < 0) {
            return false;
        }
        index = static_cast<std::size_t>(storedIndex);
        return true;
    }
};

inline std::vector<NodeId> collectSupportBoundarySiblings(
    const DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    NodeId focusNodeId,
    NodeId parentId,
    LocalCanonicalizationScratch *scratch = nullptr,
    std::span<const PixelId> localPixels = {}) {
    std::vector<NodeId> siblings;
    if (focusNodeId == InvalidNode ||
        parentId == InvalidNode ||
        !tree.isNode(focusNodeId) ||
        !tree.isAlive(focusNodeId) ||
        tree.getNumProperParts(focusNodeId) == 0) {
        return siblings;
    }
    const auto actionableSiblings = collectActionableCanonicalSiblings(
        tree,
        focusNodeId,
        parentId);
    if (actionableSiblings.empty()) {
        return siblings;
    }

    if (!localPixels.empty() &&
        localPixels.size() <= LocalCanonicalizationPixelFastPathLimit) {
        if (scratch != nullptr) {
            scratch->markActionableSiblings(
                tree,
                std::span<const NodeId>(actionableSiblings));
        }
        auto isActionableSibling = [&](NodeId nodeId, std::size_t &index) {
            if (scratch != nullptr) {
                return scratch->actionableSiblingLookup(nodeId, index);
            }
            const auto found = std::lower_bound(
                actionableSiblings.begin(),
                actionableSiblings.end(),
                nodeId);
            if (found == actionableSiblings.end() || *found != nodeId) {
                return false;
            }
            index = static_cast<std::size_t>(std::distance(actionableSiblings.begin(), found));
            return true;
        };
        std::vector<unsigned char> foundActionableSibling(actionableSiblings.size(), 0);
        std::size_t foundActionableCount = 0;
        auto addSibling = [&](NodeId siblingId) {
            std::size_t foundIndex = 0;
            if (!isActionableSibling(siblingId, foundIndex)) {
                return false;
            }
            if (foundActionableSibling[foundIndex] == 0) {
                foundActionableSibling[foundIndex] = 1;
                ++foundActionableCount;
            }
            siblings.push_back(siblingId);
            return foundActionableCount == actionableSiblings.size();
        };

        struct LocalPixelRelation {
            bool known = false;
            bool inFocus = false;
            NodeId siblingId = InvalidNode;
        };
        std::unordered_map<NodeId, LocalPixelRelation> relationByOwner;
        relationByOwner.reserve(std::max<std::size_t>(16, localPixels.size() * 2));
        auto classifyPixel = [&](PixelId pixelId) {
            LocalPixelRelation relation;
            if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
                return relation;
            }
            const NodeId ownerId = tree.getSmallestComponent(pixelId);
            if (ownerId == InvalidNode || !tree.isNode(ownerId) || !tree.isAlive(ownerId)) {
                return relation;
            }
            if (const auto found = relationByOwner.find(ownerId);
                found != relationByOwner.end()) {
                return found->second;
            }

            relation.known = true;
            NodeId current = ownerId;
            while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
                if (current == focusNodeId) {
                    relation.inFocus = true;
                    break;
                }
                const NodeId currentParent = tree.getNodeParent(current);
                if (currentParent == parentId) {
                    std::size_t localSiblingIndex = 0;
                    if (isActionableSibling(current, localSiblingIndex)) {
                        relation.siblingId = current;
                    }
                    break;
                }
                if (currentParent == current) {
                    break;
                }
                current = currentParent;
            }
            relationByOwner.emplace(ownerId, relation);
            return relation;
        };

        for (PixelId pixelId : localPixels) {
            const auto pixelRelation = classifyPixel(pixelId);
            if (!pixelRelation.known ||
                (!pixelRelation.inFocus && pixelRelation.siblingId == InvalidNode)) {
                continue;
            }
            for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
                const auto neighborRelation = classifyPixel(neighbor);
                if (!neighborRelation.known) {
                    continue;
                }
                if (pixelRelation.inFocus) {
                    if (neighborRelation.siblingId != InvalidNode &&
                        addSibling(neighborRelation.siblingId)) {
                        return sortedUniqueValidNodes(std::move(siblings));
                    }
                } else if (neighborRelation.inFocus &&
                           addSibling(pixelRelation.siblingId)) {
                    return sortedUniqueValidNodes(std::move(siblings));
                }
            }
        }
        return sortedUniqueValidNodes(std::move(siblings));
    }

    std::vector<NodeId> focusSupportNodeIds;
    const int focusSupportSize = countNodeSupportProperPartsUpTo(
        tree,
        focusNodeId,
        -1,
        &focusSupportNodeIds);
    if (focusSupportSize <= 0) {
        return siblings;
    }

    int actionableSupportSize = 0;
    for (NodeId siblingId : actionableSiblings) {
        actionableSupportSize += countNodeSupportProperPartsUpTo(
            tree,
            siblingId,
            focusSupportSize - actionableSupportSize);
        if (actionableSupportSize >= focusSupportSize) {
            break;
        }
    }
    if (actionableSupportSize > 0 && actionableSupportSize < focusSupportSize) {
        if (scratch != nullptr) {
            scratch->markFocusSupportNodes(
                tree,
                std::span<const NodeId>(focusSupportNodeIds));
        } else {
            focusSupportNodeIds = sortedUniqueValidNodes(std::move(focusSupportNodeIds));
        }
        auto focusSupportContainsOwnerPixel = [&](PixelId pixelId) {
            if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
                return false;
            }
            const NodeId ownerId = tree.getSmallestComponent(pixelId);
            if (ownerId == InvalidNode || !tree.isNode(ownerId) || !tree.isAlive(ownerId)) {
                return false;
            }
            if (scratch != nullptr) {
                return scratch->containsFocusSupportOwner(ownerId);
            }
            if (focusSupportNodeIds.size() == 1) {
                return ownerId == focusSupportNodeIds.front();
            }
            return std::binary_search(
                focusSupportNodeIds.begin(),
                focusSupportNodeIds.end(),
                ownerId);
        };
        for (NodeId siblingId : actionableSiblings) {
            bool found = false;
            for (NodeId supportNodeId : tree.getNodeSubtree(siblingId)) {
                if (supportNodeId == InvalidNode || !tree.isAlive(supportNodeId)) {
                    continue;
                }
                for (PixelId pixelId : tree.getProperParts(supportNodeId)) {
                    for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
                        if (focusSupportContainsOwnerPixel(neighbor)) {
                            siblings.push_back(siblingId);
                            found = true;
                            break;
                        }
                    }
                    if (found) {
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
        }
        return sortedUniqueValidNodes(std::move(siblings));
    }

    std::unordered_map<NodeId, NodeId> directChildByOwnerCache;
    if (scratch != nullptr) {
        scratch->beginDirectChildOwnerCache(tree);
    } else {
        directChildByOwnerCache.reserve(std::min<std::size_t>(
            static_cast<std::size_t>(tree.getNumTotalProperParts()),
            std::max<std::size_t>(16, static_cast<std::size_t>(focusSupportSize) * 4)));
    }
    auto directChildOfParentContainingPixel = [&](PixelId pixelId) {
        if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
            return InvalidNode;
        }
        NodeId current = tree.getSmallestComponent(pixelId);
        if (current == InvalidNode || !tree.isNode(current) || !tree.isAlive(current)) {
            return InvalidNode;
        }
        if (scratch != nullptr) {
            NodeId cachedDirectChildId = InvalidNode;
            if (scratch->cachedDirectChildForOwner(current, cachedDirectChildId)) {
                return cachedDirectChildId;
            }
        } else {
            if (const auto found = directChildByOwnerCache.find(current);
                found != directChildByOwnerCache.end()) {
                return found->second;
            }
        }
        const NodeId ownerId = current;
        NodeId directChildId = InvalidNode;
        while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
            if (current == focusNodeId) {
                break;
            }
            const NodeId currentParent = tree.getNodeParent(current);
            if (currentParent == parentId) {
                directChildId = current;
                break;
            }
            if (currentParent == current) {
                break;
            }
            current = currentParent;
        }
        if (directChildId != InvalidNode &&
            !std::binary_search(actionableSiblings.begin(), actionableSiblings.end(), directChildId)) {
            directChildId = InvalidNode;
        }
        if (scratch != nullptr) {
            scratch->cacheDirectChildForOwner(ownerId, directChildId);
        } else {
            directChildByOwnerCache.emplace(ownerId, directChildId);
        }
        return directChildId;
    };

    if (scratch != nullptr) {
        scratch->markActionableSiblings(
            tree,
            std::span<const NodeId>(actionableSiblings));
    }
    std::vector<unsigned char> foundActionableSibling(actionableSiblings.size(), 0);
    std::size_t foundActionableCount = 0;
    for (NodeId supportNodeId : tree.getNodeSubtree(focusNodeId)) {
        if (supportNodeId == InvalidNode || !tree.isAlive(supportNodeId)) {
            continue;
        }
        for (PixelId pixelId : tree.getProperParts(supportNodeId)) {
            for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
                const NodeId siblingId = directChildOfParentContainingPixel(neighbor);
                if (siblingId != InvalidNode && siblingId != focusNodeId) {
                    std::size_t foundIndex = 0;
                    if (scratch != nullptr) {
                        if (!scratch->actionableSiblingLookup(siblingId, foundIndex)) {
                            continue;
                        }
                    } else {
                        const auto found = std::lower_bound(
                            actionableSiblings.begin(),
                            actionableSiblings.end(),
                            siblingId);
                        if (found == actionableSiblings.end() || *found != siblingId) {
                            continue;
                        }
                        foundIndex = static_cast<std::size_t>(
                            std::distance(actionableSiblings.begin(), found));
                    }
                    if (foundActionableSibling[foundIndex] == 0) {
                        foundActionableSibling[foundIndex] = 1;
                        ++foundActionableCount;
                    }
                    siblings.push_back(siblingId);
                    if (foundActionableCount == actionableSiblings.size()) {
                        return sortedUniqueValidNodes(std::move(siblings));
                    }
                }
            }
        }
    }
    return sortedUniqueValidNodes(std::move(siblings));
}

inline bool mergeFirstAdjacentSameLevelSiblingPair(DynamicComponentTree &tree,
                                                   AdjacencyRelation &adjacency) {
    const int nodeSlots = tree.getNumInternalNodeSlots();
    for (NodeId parentId = 0; parentId < nodeSlots; ++parentId) {
        if (parentId == InvalidNode ||
            !tree.isNode(parentId) ||
            !tree.isAlive(parentId)) {
            continue;
        }

        std::vector<NodeId> children;
        for (NodeId childId : tree.getChildren(parentId)) {
            if (childId != InvalidNode && tree.isAlive(childId)) {
                children.push_back(childId);
            }
        }

        for (std::size_t i = 0; i < children.size(); ++i) {
            const NodeId keepId = children[i];
            if (!tree.isAlive(keepId)) {
                continue;
            }
            for (std::size_t j = i + 1; j < children.size(); ++j) {
                const NodeId mergeId = children[j];
                if (!tree.isAlive(mergeId) ||
                    tree.getNodeParent(keepId) != parentId ||
                    tree.getNodeParent(mergeId) != parentId ||
                    tree.getAltitude(keepId) != tree.getAltitude(mergeId)) {
                    continue;
                }
                if (!supportsAreAdjacent(tree, adjacency, keepId, mergeId)) {
                    continue;
                }
                tree.moveChildren(keepId, mergeId);
                tree.moveProperParts(keepId, mergeId);
                tree.removeChild(parentId, mergeId, true);
                return true;
            }
        }
    }
    return false;
}

inline bool mergeFirstAdjacentSameLevelSiblingPair(DynamicComponentTree &tree,
                                                   AdjacencyRelation &adjacency,
                                                   std::span<const NodeId> focusNodes) {
    if (focusNodes.empty()) {
        return mergeFirstAdjacentSameLevelSiblingPair(tree, adjacency);
    }

    for (NodeId keepId : focusNodes) {
        if (keepId == InvalidNode ||
            !tree.isNode(keepId) ||
            !tree.isAlive(keepId) ||
            tree.isRoot(keepId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(keepId);
        if (parentId == InvalidNode || parentId == keepId || !tree.isAlive(parentId)) {
            continue;
        }

        const auto boundarySiblings = collectSupportBoundarySiblings(
            tree,
            adjacency,
            keepId,
            parentId);
        if (!boundarySiblings.empty()) {
            for (NodeId mergeId : tree.getChildren(parentId)) {
                if (mergeId == InvalidNode ||
                    !std::binary_search(boundarySiblings.begin(), boundarySiblings.end(), mergeId) ||
                    !tree.isAlive(mergeId) ||
                    tree.getNodeParent(keepId) != parentId ||
                    tree.getNodeParent(mergeId) != parentId ||
                    tree.getAltitude(keepId) != tree.getAltitude(mergeId)) {
                    continue;
                }
                tree.moveChildren(keepId, mergeId);
                tree.moveProperParts(keepId, mergeId);
                tree.removeChild(parentId, mergeId, true);
                return true;
            }
        }
    }
    return false;
}

inline bool moveFirstAdjacentComparableSiblingPair(DynamicComponentTree &tree,
                                                   AdjacencyRelation &adjacency) {
    const int nodeSlots = tree.getNumInternalNodeSlots();
    for (NodeId parentId = 0; parentId < nodeSlots; ++parentId) {
        if (parentId == InvalidNode ||
            !tree.isNode(parentId) ||
            !tree.isAlive(parentId)) {
            continue;
        }

        std::vector<NodeId> children;
        for (NodeId childId : tree.getChildren(parentId)) {
            if (childId != InvalidNode && tree.isAlive(childId)) {
                children.push_back(childId);
            }
        }

        for (std::size_t i = 0; i < children.size(); ++i) {
            const NodeId lhs = children[i];
            if (!tree.isAlive(lhs)) {
                continue;
            }
            for (std::size_t j = i + 1; j < children.size(); ++j) {
                const NodeId rhs = children[j];
                if (!tree.isAlive(rhs) ||
                    tree.getNodeParent(lhs) != parentId ||
                    tree.getNodeParent(rhs) != parentId ||
                    tree.getAltitude(lhs) == tree.getAltitude(rhs)) {
                    continue;
                }

                NodeId container = lhs;
                NodeId nested = rhs;
                if (tree.isMaxtree()) {
                    if (tree.getAltitude(container) > tree.getAltitude(nested)) {
                        std::swap(container, nested);
                    }
                } else if (tree.getAltitude(container) < tree.getAltitude(nested)) {
                    std::swap(container, nested);
                }

                if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(nested)),
                                           static_cast<AltitudeType>(tree.getAltitude(container)),
                                           tree.isMaxtree())) {
                    continue;
                }
                if (!supportsAreAdjacent(tree, adjacency, container, nested)) {
                    continue;
                }

                tree.moveNode(nested, container);
                return true;
            }
        }
    }
    return false;
}

inline bool moveFirstAdjacentComparableSiblingPair(DynamicComponentTree &tree,
                                                   AdjacencyRelation &adjacency,
                                                   std::span<const NodeId> focusNodes) {
    if (focusNodes.empty()) {
        return moveFirstAdjacentComparableSiblingPair(tree, adjacency);
    }

    for (NodeId focusNodeId : focusNodes) {
        if (focusNodeId == InvalidNode ||
            !tree.isNode(focusNodeId) ||
            !tree.isAlive(focusNodeId) ||
            tree.isRoot(focusNodeId)) {
            continue;
        }

        const NodeId parentId = tree.getNodeParent(focusNodeId);
        if (parentId == InvalidNode || parentId == focusNodeId || !tree.isAlive(parentId)) {
            continue;
        }

        const auto boundarySiblings = collectSupportBoundarySiblings(
            tree,
            adjacency,
            focusNodeId,
            parentId);
        if (!boundarySiblings.empty()) {
            for (NodeId siblingId : tree.getChildren(parentId)) {
                if (siblingId == InvalidNode ||
                    !std::binary_search(boundarySiblings.begin(), boundarySiblings.end(), siblingId) ||
                    !tree.isAlive(focusNodeId) ||
                    !tree.isAlive(siblingId) ||
                    tree.getNodeParent(focusNodeId) != parentId ||
                    tree.getNodeParent(siblingId) != parentId ||
                    tree.getAltitude(focusNodeId) == tree.getAltitude(siblingId)) {
                    continue;
                }

                NodeId container = focusNodeId;
                NodeId nested = siblingId;
                if (tree.isMaxtree()) {
                    if (tree.getAltitude(container) > tree.getAltitude(nested)) {
                        std::swap(container, nested);
                    }
                } else if (tree.getAltitude(container) < tree.getAltitude(nested)) {
                    std::swap(container, nested);
                }

                if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(nested)),
                                           static_cast<AltitudeType>(tree.getAltitude(container)),
                                           tree.isMaxtree())) {
                    continue;
                }
                tree.moveNode(nested, container);
                return true;
            }
        }
    }
    return false;
}

struct LocalFocusBoundarySiblings {
    NodeId focusNodeId = InvalidNode;
    NodeId parentId = InvalidNode;
    std::vector<NodeId> siblingIds;
};

inline LocalFocusBoundarySiblings collectFocusBoundarySiblings(
    const DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    NodeId focusNodeId,
    LocalCanonicalizationScratch *scratch = nullptr,
    std::span<const PixelId> localPixels = {}) {
    if (focusNodeId == InvalidNode ||
        !tree.isNode(focusNodeId) ||
        !tree.isAlive(focusNodeId) ||
        tree.isRoot(focusNodeId)) {
        return {};
    }

    const NodeId parentId = tree.getNodeParent(focusNodeId);
    if (parentId == InvalidNode || parentId == focusNodeId || !tree.isAlive(parentId)) {
        return {};
    }

    auto siblingIds = collectSupportBoundarySiblings(
        tree,
        adjacency,
        focusNodeId,
        parentId,
        scratch,
        localPixels);
    if (siblingIds.empty()) {
        return {};
    }
    return {focusNodeId, parentId, std::move(siblingIds)};
}

inline bool mergeFocusedAdjacentSameLevelSiblingPair(
    DynamicComponentTree &tree,
    const LocalFocusBoundarySiblings &record) {
    const NodeId keepId = record.focusNodeId;
    const NodeId parentId = record.parentId;
    if (keepId == InvalidNode ||
        parentId == InvalidNode ||
        !tree.isNode(keepId) ||
        !tree.isAlive(keepId) ||
        !tree.isAlive(parentId) ||
        tree.getNodeParent(keepId) != parentId) {
        return false;
    }

    for (NodeId mergeId : tree.getChildren(parentId)) {
        if (mergeId == InvalidNode ||
            !std::binary_search(record.siblingIds.begin(), record.siblingIds.end(), mergeId) ||
            !tree.isAlive(mergeId) ||
            tree.getNodeParent(keepId) != parentId ||
            tree.getNodeParent(mergeId) != parentId ||
            tree.getAltitude(keepId) != tree.getAltitude(mergeId)) {
            continue;
        }
        tree.moveChildren(keepId, mergeId);
        tree.moveProperParts(keepId, mergeId);
        tree.removeChild(parentId, mergeId, true);
        return true;
    }
    return false;
}

inline bool mergeFirstFocusedAdjacentSameLevelSiblingPair(
    DynamicComponentTree &tree,
    std::span<const LocalFocusBoundarySiblings> records) {
    for (const auto &record : records) {
        if (mergeFocusedAdjacentSameLevelSiblingPair(tree, record)) {
            return true;
        }
    }
    return false;
}

inline bool moveFirstFocusedAdjacentComparableSiblingPair(
    DynamicComponentTree &tree,
    std::span<const LocalFocusBoundarySiblings> records) {
    for (const auto &record : records) {
        const NodeId focusNodeId = record.focusNodeId;
        const NodeId parentId = record.parentId;
        if (focusNodeId == InvalidNode ||
            parentId == InvalidNode ||
            !tree.isNode(focusNodeId) ||
            !tree.isAlive(focusNodeId) ||
            !tree.isAlive(parentId) ||
            tree.getNodeParent(focusNodeId) != parentId) {
            continue;
        }

        for (NodeId siblingId : tree.getChildren(parentId)) {
            if (siblingId == InvalidNode ||
                !std::binary_search(record.siblingIds.begin(), record.siblingIds.end(), siblingId) ||
                !tree.isAlive(focusNodeId) ||
                !tree.isAlive(siblingId) ||
                tree.getNodeParent(focusNodeId) != parentId ||
                tree.getNodeParent(siblingId) != parentId ||
                tree.getAltitude(focusNodeId) == tree.getAltitude(siblingId)) {
                continue;
            }

            NodeId container = focusNodeId;
            NodeId nested = siblingId;
            if (tree.isMaxtree()) {
                if (tree.getAltitude(container) > tree.getAltitude(nested)) {
                    std::swap(container, nested);
                }
            } else if (tree.getAltitude(container) < tree.getAltitude(nested)) {
                std::swap(container, nested);
            }

            if (!validLocalParentLevel(static_cast<AltitudeType>(tree.getAltitude(nested)),
                                       static_cast<AltitudeType>(tree.getAltitude(container)),
                                       tree.isMaxtree())) {
                continue;
            }
            tree.moveNode(nested, container);
            return true;
        }
    }
    return false;
}

enum class LocalCanonicalSiblingAction {
    None,
    MergeSameLevelSibling,
    MoveComparableSibling,
};

inline LocalCanonicalSiblingAction canonicalizeFirstFocusedAdjacentSiblingPairAction(
    DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    std::span<const NodeId> focusNodes,
    LocalCanonicalizationScratch *scratch = nullptr,
    std::span<const PixelId> localPixels = {}) {
    std::vector<LocalFocusBoundarySiblings> recordsWithoutMerge;
    recordsWithoutMerge.reserve(focusNodes.size());

    for (NodeId focusNodeId : focusNodes) {
        auto record = collectFocusBoundarySiblings(
            tree,
            adjacency,
            focusNodeId,
            scratch,
            localPixels);
        if (record.siblingIds.empty()) {
            continue;
        }
        if (mergeFocusedAdjacentSameLevelSiblingPair(tree, record)) {
            return LocalCanonicalSiblingAction::MergeSameLevelSibling;
        }
        recordsWithoutMerge.push_back(std::move(record));
    }
    if (moveFirstFocusedAdjacentComparableSiblingPair(
            tree,
            std::span<const LocalFocusBoundarySiblings>(recordsWithoutMerge))) {
        return LocalCanonicalSiblingAction::MoveComparableSibling;
    }
    return LocalCanonicalSiblingAction::None;
}

inline bool canonicalizeFirstFocusedAdjacentSiblingPair(DynamicComponentTree &tree,
                                                        AdjacencyRelation &adjacency,
                                                        std::span<const NodeId> focusNodes,
                                                        LocalCanonicalizationScratch *scratch = nullptr,
                                                        std::span<const PixelId> localPixels = {}) {
    return canonicalizeFirstFocusedAdjacentSiblingPairAction(
               tree,
               adjacency,
               focusNodes,
               scratch,
               localPixels) != LocalCanonicalSiblingAction::None;
}

inline bool isOracleRootPlan(const LocalRecompositionNodePlan &nodePlan) {
    return nodePlan.oracleNodeId != InvalidNode &&
           nodePlan.oracleParentId == nodePlan.oracleNodeId;
}

inline int promoteOracleRootPlan(DynamicComponentTree &tree,
                                 const LocalRecompositionPlan &plan,
                                 std::span<const NodeId> createdNodeIds) {
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size(); ++planIndex) {
        if (!isOracleRootPlan(plan.nodes[planIndex])) {
            continue;
        }
        const NodeId newRoot = createdNodeIds[planIndex];
        if (newRoot == InvalidNode || !tree.isAlive(newRoot) || tree.isRoot(newRoot)) {
            return 0;
        }
        const NodeId oldRoot = tree.getRoot();
        tree.setRoot(newRoot);
        if (oldRoot != InvalidNode &&
            oldRoot != newRoot &&
            tree.isNode(oldRoot) &&
            tree.isAlive(oldRoot)) {
            tree.moveChildren(newRoot, oldRoot);
            if (tree.getNumProperParts(oldRoot) > 0 &&
                tree.getAltitude(oldRoot) == tree.getAltitude(newRoot)) {
                tree.moveProperParts(newRoot, oldRoot);
            }
            if (tree.getNumProperParts(oldRoot) == 0 &&
                tree.getNumChildren(oldRoot) == 0) {
                tree.releaseNode(oldRoot);
            }
        }
        return 1;
    }
    return 0;
}

inline int canonicalizeAfterLocalRecomposition(DynamicComponentTree &tree,
                                               AdjacencyRelation *adjacency,
                                               std::span<const NodeId> focusNodesInput = {},
                                               std::span<const NodeId> protectedNodesInput = {},
                                               std::span<const PixelId> localPixels = {},
                                               LocalCanonicalizationBreakdown *breakdownOut = nullptr) {
    std::vector<NodeId> compactProtectedNodes;
    std::vector<NodeId> compactFocusNodes;
    std::span<const NodeId> focusNodes = focusNodesInput;
    std::span<const NodeId> protectedNodes = protectedNodesInput;
    if (!focusNodesInput.empty()) {
        compactFocusNodes.assign(focusNodesInput.begin(), focusNodesInput.end());
        compactFocusNodes = stableUniqueValidNodes(std::move(compactFocusNodes));
        focusNodes = std::span<const NodeId>(compactFocusNodes);
    }
    if (!protectedNodesInput.empty()) {
        compactProtectedNodes.assign(protectedNodesInput.begin(), protectedNodesInput.end());
        compactProtectedNodes = stableUniqueValidNodes(std::move(compactProtectedNodes));
        protectedNodes = std::span<const NodeId>(compactProtectedNodes);
    }

    LocalCanonicalizationBreakdown breakdown;
    bool changed = true;
    LocalCanonicalizationScratch scratch;
    while (changed) {
        changed = false;
        const int released = focusNodes.empty()
                                 ? releaseAllEmptyNonRootNodes(tree, protectedNodes)
                                 : releaseFocusedEmptyNonRootNodes(tree, focusNodes, protectedNodes);
        if (released > 0) {
            breakdown.emptyNodeReleases += released;
            changed = true;
            continue;
        }
        if (promoteEmptyRootWithSingleChild(tree)) {
            ++breakdown.emptyNodeReleases;
            changed = true;
            continue;
        }
        if (focusNodes.empty()
                ? collapseSameLevelParentChildNodes(tree)
                : collapseFocusedSameLevelParentChildNodes(tree, focusNodes)) {
            ++breakdown.sameLevelParentChildCollapses;
            changed = true;
            continue;
        }
        if (adjacency != nullptr) {
            if (focusNodes.empty()) {
                if (mergeFirstAdjacentSameLevelSiblingPair(tree, *adjacency)) {
                    ++breakdown.sameLevelSiblingMerges;
                    changed = true;
                    continue;
                }
                if (moveFirstAdjacentComparableSiblingPair(tree, *adjacency)) {
                    ++breakdown.comparableSiblingMoves;
                    changed = true;
                }
            } else {
                const LocalCanonicalSiblingAction siblingAction =
                    canonicalizeFirstFocusedAdjacentSiblingPairAction(
                        tree,
                        *adjacency,
                        focusNodes,
                        &scratch,
                        localPixels);
                if (siblingAction == LocalCanonicalSiblingAction::MergeSameLevelSibling) {
                    ++breakdown.sameLevelSiblingMerges;
                    changed = true;
                    continue;
                }
                if (siblingAction == LocalCanonicalSiblingAction::MoveComparableSibling) {
                    ++breakdown.comparableSiblingMoves;
                    changed = true;
                }
            }
        }
    }
    if (breakdownOut != nullptr) {
        *breakdownOut = breakdown;
    }
    return breakdown.total();
}

inline void validateLocalRecompositionMaterializationInputs(
    const DynamicComponentTree &dualTree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    const LocalRecompositionMaterializationOptions &options = {}) {
    if (!localMergeOracleHasUpdatedTarget(oracle)) {
        throw std::runtime_error("local recomposition materialization requires target values");
    }
    if (options.checkPreReconstruction && !oracle.preImage) {
        throw std::runtime_error("local recomposition materialization requires a pre-update image");
    }
    if (dualTree.isMaxtree() != oracle.contraction.dualIsMaxtree) {
        throw std::runtime_error("local recomposition materialization received the wrong dual-tree polarity");
    }
    if (localMergeOracleDomainSize(oracle) != dualTree.getNumTotalProperParts() ||
        (oracle.preImage &&
         oracle.preImage->getSize() != dualTree.getNumTotalProperParts())) {
        throw std::runtime_error("local recomposition materialization received images with incompatible size");
    }
    if (options.checkTargetReconstruction && !oracle.updatedImage) {
        throw std::runtime_error(
            "local recomposition target reconstruction check requires a full target image");
    }
    if (options.checkPreReconstruction &&
        !dualTree.reconstructionImage()->isEqual(oracle.preImage)) {
        throw std::runtime_error("local recomposition materialization requires the pre-update dual tree");
    }

    std::size_t claimedPixelCapacity = 0;
    for (const auto &nodePlan : plan.nodes) {
        claimedPixelCapacity += nodePlan.directProperParts.size();
    }
    std::vector<PixelId> claimedPixels;
    claimedPixels.reserve(claimedPixelCapacity);
    for (const auto &nodePlan : plan.nodes) {
        if (nodePlan.directProperParts.empty() &&
            nodePlan.childPlanIndices.empty() &&
            nodePlan.oracleChildrenInSnapshot.empty() &&
            nodePlan.stableChildPreNodes.empty()) {
            throw std::runtime_error("local recomposition materialization rejects empty isolated recomposition nodes");
        }
        for (PixelId pixelId : nodePlan.directProperParts) {
            if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
                throw std::runtime_error("local recomposition materialization received an out-of-bounds pixel");
            }
            if (localMergeOracleUpdatedPixelValue(oracle, pixelId) != nodePlan.level) {
                throw std::runtime_error("local recomposition direct proper part does not match the planned node level");
            }
            claimedPixels.push_back(pixelId);
        }
    }
    std::sort(claimedPixels.begin(), claimedPixels.end());
    if (std::adjacent_find(claimedPixels.begin(), claimedPixels.end()) != claimedPixels.end()) {
        throw std::runtime_error("local recomposition materialization received duplicated direct proper parts");
    }
}

inline LocalRecompositionMaterializationResult materializeLocalRecompositionNodes(
    DynamicComponentTree &dualTree,
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &plan,
    AdjacencyRelation *adjacency = nullptr,
    const LocalRecompositionMaterializationOptions &options = {}) {
    auto stageStart = std::chrono::steady_clock::now();
    validateLocalRecompositionMaterializationInputs(dualTree, oracle, plan, options);

    LocalRecompositionMaterializationResult result;
    result.validationMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());
    std::vector<NodeId> createdNodeIds(plan.nodes.size(), InvalidNode);
    std::vector<NodeId> canonicalFocusNodeIds;
    std::vector<LocalCreatedRecompositionNode> hiddenCreatedNodes;
    canonicalFocusNodeIds.reserve(plan.nodes.size());

    stageStart = std::chrono::steady_clock::now();
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size(); ++planIndex) {
        const auto &nodePlan = plan.nodes[planIndex];
        const NodeId createdNodeId = dualTree.createDetachedNode(static_cast<int>(nodePlan.level));
        if (createdNodeId == InvalidNode) {
            throw std::runtime_error("local recomposition materialization could not create a node");
        }
        createdNodeIds[planIndex] = createdNodeId;
        canonicalFocusNodeIds.push_back(createdNodeId);
    }
    result.createNodesMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());

    stageStart = std::chrono::steady_clock::now();
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size(); ++planIndex) {
        const auto &nodePlan = plan.nodes[planIndex];
        const NodeId createdNodeId = createdNodeIds[planIndex];
        NodeId parentNodeId = InvalidNode;
        if (nodePlan.parentPlanIndex >= 0) {
            const auto parentIndex = static_cast<std::size_t>(nodePlan.parentPlanIndex);
            if (parentIndex >= createdNodeIds.size() || createdNodeIds[parentIndex] == InvalidNode) {
                throw std::runtime_error("local recomposition materialization has an invalid local parent index");
            }
            if (!validLocalParentLevel(nodePlan.level,
                                       plan.nodes[parentIndex].level,
                                       dualTree.isMaxtree())) {
                throw std::runtime_error("local recomposition materialization has an invalid local parent level");
            }
            parentNodeId = createdNodeIds[parentIndex];
        } else {
            parentNodeId = externalAttachmentForRecompositionNode(dualTree, oracle, nodePlan);
            if (nodePlan.oracleNodeId != InvalidNode &&
                oracleParentIsOutsideLocalSnapshot(oracle, nodePlan.oracleParentId)) {
                const NodeId hiddenParentId = materializeHiddenExternalParentConnector(
                    dualTree,
                    nodePlan,
                    parentNodeId,
                    canonicalFocusNodeIds,
                    hiddenCreatedNodes);
                if (hiddenParentId != InvalidNode) {
                    parentNodeId = hiddenParentId;
                }
            } else if (shouldMaterializeNonOracleHiddenExternalParentConnector(
                           dualTree,
                           oracle,
                           nodePlan,
                           parentNodeId)) {
                const NodeId hiddenParentId = materializeHiddenExternalParentConnector(
                    dualTree,
                    nodePlan,
                    parentNodeId,
                    canonicalFocusNodeIds,
                    hiddenCreatedNodes);
                if (hiddenParentId != InvalidNode) {
                    parentNodeId = hiddenParentId;
                }
            }
        }
        dualTree.attachNode(parentNodeId, createdNodeId);
    }
    result.attachNodesMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());

    std::vector<NodeId> emptiedSourceOwners;
    stageStart = std::chrono::steady_clock::now();
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size(); ++planIndex) {
        const auto &nodePlan = plan.nodes[planIndex];
        const NodeId createdNodeId = createdNodeIds[planIndex];

        LocalCreatedRecompositionNode created;
        created.planIndex = static_cast<int>(planIndex);
        created.nodeId = createdNodeId;
        created.parentNodeId = dualTree.getNodeParent(createdNodeId);
        created.level = nodePlan.level;
        created.movedDirectProperParts = nodePlan.directProperParts;

        for (PixelId pixelId : nodePlan.directProperParts) {
            const NodeId sourceOwner = dualTree.getSmallestComponent(pixelId);
            if (sourceOwner == InvalidNode || !dualTree.isAlive(sourceOwner)) {
                throw std::runtime_error("local recomposition materialization cannot locate a live source owner");
            }
            created.sourceOwnersBefore.push_back(sourceOwner);
            if (sourceOwner != createdNodeId) {
                dualTree.moveProperPart(createdNodeId, sourceOwner, pixelId);
                if (dualTree.getNumProperParts(sourceOwner) == 0) {
                    emptiedSourceOwners.push_back(sourceOwner);
                }
            }
        }
        created.sourceOwnersBefore = sortedUniqueValidNodes(std::move(created.sourceOwnersBefore));
        result.createdNodes.push_back(std::move(created));
    }
    result.createdNodes.insert(
        result.createdNodes.end(),
        hiddenCreatedNodes.begin(),
        hiddenCreatedNodes.end());
    result.moveProperPartsMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());

    std::vector<std::pair<NodeId, NodeId>> bridgeNodes;
    stageStart = std::chrono::steady_clock::now();
    for (std::size_t planIndex = 0; planIndex < plan.nodes.size(); ++planIndex) {
        const auto &nodePlan = plan.nodes[planIndex];
        const NodeId createdNodeId = createdNodeIds[planIndex];
        auto stableChildAttachments = collectStableChildAttachmentsForRecompositionNode(
            dualTree,
            oracle,
            nodePlan,
            createdNodeId);
        std::vector<NodeId> stablePreMovedChildren = moveStablePreChildNodesForRecompositionNode(
            dualTree,
            nodePlan,
            createdNodeId,
            canonicalFocusNodeIds,
            emptiedSourceOwners);
        std::vector<NodeId> movedChildren = stablePreMovedChildren;
        materializeStablePreChildDescendantBridges(
            dualTree,
            oracle,
            std::span<const NodeId>(stablePreMovedChildren),
            canonicalFocusNodeIds,
            emptiedSourceOwners,
            result);
        for (const auto &attachment : stableChildAttachments) {
            const NodeId childNodeId = attachment.childNodeId;
            if (dualTree.isAlive(childNodeId) &&
                dualTree.getNodeParent(childNodeId) != createdNodeId) {
                const NodeId previousParent = dualTree.getNodeParent(childNodeId);
                dualTree.moveNode(childNodeId, createdNodeId);
                canonicalFocusNodeIds.push_back(childNodeId);
                recordEmptyPreviousParent(
                    dualTree,
                    previousParent,
                    createdNodeId,
                    emptiedSourceOwners);
            } else if (dualTree.isAlive(childNodeId)) {
                canonicalFocusNodeIds.push_back(childNodeId);
            }
            movedChildren.push_back(childNodeId);
            materializeMissingOracleDescendantBridges(
                dualTree,
                oracle,
                plan,
                attachment.oracleChildId,
                childNodeId,
                canonicalFocusNodeIds,
                emptiedSourceOwners,
                result,
                bridgeNodes);
        }
        auto bridgeChildren = materializeMissingOracleChildBridges(
            dualTree,
            oracle,
            plan,
            nodePlan,
            createdNodeId,
            canonicalFocusNodeIds,
            emptiedSourceOwners,
            result,
            bridgeNodes);
        movedChildren.insert(movedChildren.end(), bridgeChildren.begin(), bridgeChildren.end());
        result.createdNodes[planIndex].movedChildNodes =
            sortedUniqueValidNodes(std::move(movedChildren));
    }
    result.childBridgeMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());

    stageStart = std::chrono::steady_clock::now();
    enforceOracleSnapshotParentLinks(
        dualTree,
        oracle,
        plan,
        std::span<const NodeId>(createdNodeIds),
        canonicalFocusNodeIds,
        emptiedSourceOwners,
        result,
        bridgeNodes);

    materializeHiddenSupportConnectors(
        dualTree,
        oracle,
        plan,
        std::span<const NodeId>(createdNodeIds),
        canonicalFocusNodeIds,
        emptiedSourceOwners,
        result,
        bridgeNodes);

    materializeNonOracleHiddenSupportConnectors(
        dualTree,
        oracle,
        plan,
        std::span<const NodeId>(createdNodeIds),
        canonicalFocusNodeIds,
        emptiedSourceOwners,
        result);

    enforceNonOraclePreSnapshotParentLinks(
        dualTree,
        oracle,
        plan,
        std::span<const NodeId>(createdNodeIds),
        canonicalFocusNodeIds,
        emptiedSourceOwners,
        result);
    result.structuralRepairMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());

    result.emptiedSourceOwners = sortedUniqueValidNodes(std::move(emptiedSourceOwners));
    std::vector<NodeId> autoPromotableResidualChildren = canonicalFocusNodeIds;
    autoPromotableResidualChildren.insert(
        autoPromotableResidualChildren.end(),
        createdNodeIds.begin(),
        createdNodeIds.end());
    for (const auto &created : result.createdNodes) {
        autoPromotableResidualChildren.push_back(created.nodeId);
        autoPromotableResidualChildren.insert(
            autoPromotableResidualChildren.end(),
            created.movedChildNodes.begin(),
            created.movedChildNodes.end());
    }
    autoPromotableResidualChildren =
        stableUniqueValidNodes(std::move(autoPromotableResidualChildren));
    stageStart = std::chrono::steady_clock::now();
    releaseEmptySourceOwners(
        dualTree,
        result.emptiedSourceOwners,
        result.releasedSourceOwners,
        options.releaseEmptySourceOwnersWithResidualChildren,
        std::span<const NodeId>(autoPromotableResidualChildren));
    result.releaseEmptySourceOwnersMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());
    stageStart = std::chrono::steady_clock::now();
    (void) promoteOracleRootPlan(
        dualTree,
        plan,
        std::span<const NodeId>(createdNodeIds));
    result.rootPromotionMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());
    stageStart = std::chrono::steady_clock::now();
    std::vector<NodeId> canonicalizationFocusNodeIds = canonicalFocusNodeIds;
    for (PixelId pixelId : oracle.localPixels) {
        if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
            continue;
        }
        NodeId current = dualTree.getSmallestComponent(pixelId);
        while (current != InvalidNode &&
               dualTree.isNode(current) &&
               dualTree.isAlive(current)) {
            canonicalizationFocusNodeIds.push_back(current);
            const NodeId parentId = dualTree.getNodeParent(current);
            if (parentId == current || parentId == InvalidNode || !dualTree.isAlive(parentId)) {
                break;
            }
            current = parentId;
        }
    }
    canonicalizationFocusNodeIds =
        stableUniqueValidNodes(std::move(canonicalizationFocusNodeIds));
    std::span<const NodeId> protectedCanonicalizationNodes;
    std::vector<NodeId> protectedResidualSourceOwners;
    if (!options.releaseEmptySourceOwnersWithResidualChildren) {
        auto parentIsPlannedPreOwner = [&](NodeId parentId) {
            if (parentId == InvalidNode) {
                return false;
            }
            for (const auto &nodePlan : plan.nodes) {
                if (std::find(nodePlan.preDualOwners.begin(),
                              nodePlan.preDualOwners.end(),
                              parentId) != nodePlan.preDualOwners.end()) {
                    return true;
                }
            }
            return false;
        };
        for (NodeId sourceOwner : result.emptiedSourceOwners) {
            if (sourceOwner != InvalidNode &&
                dualTree.isNode(sourceOwner) &&
                dualTree.isAlive(sourceOwner) &&
                dualTree.getNumChildren(sourceOwner) != 0 &&
                parentIsPlannedPreOwner(dualTree.getNodeParent(sourceOwner))) {
                protectedResidualSourceOwners.push_back(sourceOwner);
            }
        }
        protectedResidualSourceOwners =
            stableUniqueValidNodes(std::move(protectedResidualSourceOwners));
        protectedCanonicalizationNodes =
            std::span<const NodeId>(protectedResidualSourceOwners);
    }
    LocalCanonicalizationBreakdown canonicalizationBreakdown;
    (void) canonicalizeAfterLocalRecomposition(
        dualTree,
        adjacency,
        std::span<const NodeId>(canonicalizationFocusNodeIds),
        protectedCanonicalizationNodes,
        std::span<const PixelId>(oracle.localPixels),
        &canonicalizationBreakdown);
    result.canonicalizationBreakdown = canonicalizationBreakdown;
    result.canonicalizationMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());
    result.emptyCreatedNodesAfterCanonicalization =
        collectEmptyCreatedNodeDiagnostics(
            dualTree,
            std::span<const LocalCreatedRecompositionNode>(result.createdNodes));
    if (options.checkNoEmptyLiveNodesAfterCanonicalization) {
        result.emptyLiveNodesAfterCanonicalization = countAliveNodesWithoutProperParts(dualTree);
        if (result.emptyLiveNodesAfterCanonicalization != 0) {
            throw std::runtime_error(
                "local recomposition materialization left alive nodes without proper parts:" +
                describeAliveNodesWithoutProperParts(dualTree));
        }
    }
    stageStart = std::chrono::steady_clock::now();
    if (options.checkTargetReconstruction) {
        result.reconstructionMatchesTarget = dualTree.reconstructionImage()->isEqual(oracle.updatedImage);
        result.targetReconstructionChecked = true;
    }
    result.targetReconstructionCheckMicros =
        localRecomposerElapsedMicros(stageStart, std::chrono::steady_clock::now());
    return result;
}

inline std::string debugString(const LocalRecompositionMaterializationResult &result) {
    std::ostringstream out;
    out << "LocalRecompositionMaterializationResult\n";
    out << "created_nodes=" << result.createdNodes.size()
        << " reconstruction_matches_target=" << (result.reconstructionMatchesTarget ? "true" : "false")
        << " empty_live_nodes_after_canonicalization=" << result.emptyLiveNodesAfterCanonicalization
        << " canonicalization_empty_node_releases=" << result.canonicalizationBreakdown.emptyNodeReleases
        << " canonicalization_same_level_parent_child_collapses="
        << result.canonicalizationBreakdown.sameLevelParentChildCollapses
        << " canonicalization_same_level_sibling_merges="
        << result.canonicalizationBreakdown.sameLevelSiblingMerges
        << " canonicalization_comparable_sibling_moves="
        << result.canonicalizationBreakdown.comparableSiblingMoves
        << "\n";
    for (const auto &created : result.createdNodes) {
        out << "  plan_index=" << created.planIndex
            << " node=" << created.nodeId
            << " parent=" << created.parentNodeId
            << " level=" << static_cast<int>(created.level)
            << " moved=";
        appendPixelVector(out, std::span<const PixelId>(created.movedDirectProperParts));
        out << " source_owners=";
        appendNodeVector(out, std::span<const NodeId>(created.sourceOwnersBefore));
        out << " moved_children=";
        appendNodeVector(out, std::span<const NodeId>(created.movedChildNodes));
        out << "\n";
    }
    out << "emptied_source_owners=";
    appendNodeVector(out, std::span<const NodeId>(result.emptiedSourceOwners));
    out << "\n";
    out << "released_source_owners=";
    appendNodeVector(out, std::span<const NodeId>(result.releasedSourceOwners));
    out << "\n";
    return out.str();
}

} // namespace morphotree::details::non_increasing
