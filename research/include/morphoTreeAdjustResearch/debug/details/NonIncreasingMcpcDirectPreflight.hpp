#pragma once

#include <algorithm>
#include <map>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"
#include "../../details/NonIncreasingLocalPartition.hpp"
#include "../../details/NonIncreasingLocalPartitionPlanner.hpp"
#include "../../details/NonIncreasingLocalRecomposer.hpp"
#include "../../details/NonIncreasingLocalUpdate.hpp"
#include "../../details/NonIncreasingStrippedLocalRecomposition.hpp"

namespace debug::non_increasing {

using namespace ::morphotree::details::non_increasing;

struct NonIncreasingMcpcDirectMutationSummary {
    int affectedLevelCount = 0;
    int localPixelCount = 0;
    int preLocalNodeCount = 0;
    int partitionRecomposedComponentCount = 0;
    int plannedNodeCount = 0;
    int plannedDirectProperPartCount = 0;
    int plannedStableChildRefCount = 0;
    int plannedExternalAttachmentCount = 0;
    int sourceOwnerTouchedCount = 0;
    int sourceOwnerToBeEmptiedCount = 0;
    int stableChildToMoveCount = 0;
    int previousParentTouchedCount = 0;
    int previousParentPossiblyEmptiedCount = 0;
    int invalidStableChildRefCount = 0;
    int stableChildCycleRiskCount = 0;
    int predictedAuxiliaryNodeCount = 0;
    int predictedHiddenExternalParentConnectorNodeCount = 0;
    int predictedStableChildBridgeNodeCount = 0;
    int candidateNonOracleHiddenSupportConnectorCount = 0;
    int candidateNonOracleParentLinkBridgeCount = 0;
    int candidateNonOracleChildLinkBridgeCount = 0;
    std::vector<int> plannedNodeLevels;
};

struct NonIncreasingMcpcDirectV1PrecheckSummary {
    int predictedAuxiliaryNodeCount = 0;
    int sameLevelPlannedParentChildCount = 0;
    int sameLevelExternalAttachmentCount = 0;
    int nonOracleHiddenSupportCandidateCount = 0;
    int nonOracleHiddenSupportAlreadyRepresentedCount = 0;
    int nonOracleHiddenSupportNotRedundantCount = 0;
    int nonOracleParentLinkTouchesPlanCount = 0;
    int nonOracleParentLinkBothRepresentedCount = 0;
    int nonOracleParentLinkNotBothRepresentedCount = 0;
    int nonOracleParentLinkBlockedByModifiedSnapshotCount = 0;
    int nonOracleParentLinkInvalidLevelCount = 0;
};

struct NonIncreasingMcpcDirectPreflightResult {
    bool checked = false;
    bool inputBuilt = false;
    bool partitionPlanBuilt = false;
    bool recompositionPlanBuilt = false;
    bool planOracleAuditChecked = false;
    bool planOracleFree = false;
    bool materializationContextBuilt = false;
    bool materializationOracleAuditChecked = false;
    bool materializationOracleFree = false;
    bool materializationInputsValid = false;
    bool rollbackFreeSummaryBuilt = false;
    bool passed = false;
    bool directV1PrecheckEvaluated = false;
    bool directV1PrecheckAccepted = false;
    NonIncreasingMcpcDirectMutationSummary summary;
    NonIncreasingMcpcDirectV1PrecheckSummary directV1Precheck;
    std::string failureSummary;
    std::string directV1PrecheckRejectionSummary;
};

struct NonIncreasingMcpcDirectV1AcceptanceDiagnostics {
    bool evaluated = false;
    bool accepted = false;
    std::string rejectionSummary;
};

struct NonIncreasingMcpcDirectSummaryComparison {
    bool evaluated = false;
    bool matches = false;
    int mismatchCount = 0;
    int plannedNodeCountMismatch = 0;
    int plannedDirectProperPartCountMismatch = 0;
    int predictedAuxiliaryNodeCountMismatch = 0;
    int predictedHiddenExternalParentConnectorMismatch = 0;
    int predictedStableChildBridgeMismatch = 0;
    int emptiedLocalNodeCountMismatch = 0;
    int movedStableChildCountMismatch = 0;
    int predictedEmptiedLocalNodeCount = 0;
    int observedEmptiedLocalNodeCount = 0;
    std::string mismatchSummary;
};

struct NonIncreasingMcpcNoCopyCandidateDiagnostics {
    bool evaluated = false;
    bool accepted = false;
    std::string rejectionSummary;
};

struct NonIncreasingMcpcDirectCommitCandidateDiagnostics {
    bool evaluated = false;
    bool accepted = false;
    std::string rejectionSummary;
};

struct NonIncreasingMcpcD12DomainDiagnostics {
    bool evaluated = false;
    bool accepted = false;
    std::string rejectReasonCode;
    std::string rejectionSummary;
};

struct NonIncreasingMcpcLiveMutationDryRun {
    bool evaluated = false;
    bool accepted = false;
    std::string rejectReasonCode;
    std::string rejectionSummary;
    int plannedNodeCount = 0;
    int plannedDirectProperPartCount = 0;
    int plannedExternalAttachmentCount = 0;
    int plannedStableChildRefCount = 0;
    int sourceOwnerToBeEmptiedCount = 0;
    int previousParentPossiblyEmptiedCount = 0;
    int stableChildToMoveCount = 0;
    int predictedAuxiliaryNodeCount = 0;
    int predictedHiddenExternalParentConnectorNodeCount = 0;
    int predictedStableChildBridgeNodeCount = 0;
    int candidateNonOracleHiddenSupportConnectorCount = 0;
    int candidateNonOracleParentLinkBridgeCount = 0;
    int candidateNonOracleChildLinkBridgeCount = 0;
};

struct NonIncreasingMcpcDirectPreflightArtifacts {
    NonIncreasingMcpcDirectPreflightResult preflight;
    LocalPartitionInput input;
    LocalPartitionPlan partitionPlan;
    LocalRecompositionPlan recompositionPlan;
    LocalMergeOracle materializationContext;
};

inline std::vector<NodeId> mcpcPreflightSortedUniqueNodes(std::vector<NodeId> nodes) {
    nodes.erase(std::remove(nodes.begin(), nodes.end(), InvalidNode), nodes.end());
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return nodes;
}

inline bool mcpcPreflightPixelSetsIntersect(std::span<const PixelId> lhs,
                                            std::span<const PixelId> rhs) {
    auto left = sortedUniquePixels(std::vector<PixelId>(lhs.begin(), lhs.end()));
    auto right = sortedUniquePixels(std::vector<PixelId>(rhs.begin(), rhs.end()));
    std::size_t leftIndex = 0;
    std::size_t rightIndex = 0;
    while (leftIndex < left.size() && rightIndex < right.size()) {
        if (left[leftIndex] == right[rightIndex]) {
            return true;
        }
        if (left[leftIndex] < right[rightIndex]) {
            ++leftIndex;
        } else {
            ++rightIndex;
        }
    }
    return false;
}

inline bool mcpcPreflightNodeIsAncestorOf(const DynamicComponentTree &tree,
                                          NodeId ancestorId,
                                          NodeId nodeId) {
    if (ancestorId == InvalidNode || nodeId == InvalidNode ||
        !tree.isAlive(ancestorId) || !tree.isAlive(nodeId)) {
        return false;
    }
    NodeId current = nodeId;
    while (current != InvalidNode && tree.isAlive(current)) {
        if (current == ancestorId) {
            return true;
        }
        const NodeId parent = tree.getNodeParent(current);
        if (parent == current) {
            break;
        }
        current = parent;
    }
    return false;
}

inline int countMcpcPreflightStableChildBridgeDescendants(
    const DynamicComponentTree &tree,
    std::span<const LocalDualNodeSnapshot> preSnapshots,
    NodeId parentNodeId,
    const LocalDualNodeSnapshot &parentSnapshot) {
    if (parentNodeId == InvalidNode || !tree.isAlive(parentNodeId)) {
        return 0;
    }

    int count = 0;
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
            ++count;
        }

        count += countMcpcPreflightStableChildBridgeDescendants(
            tree,
            preSnapshots,
            childNodeId,
            *childSnapshot);
    }
    return count;
}

inline int countMcpcPreflightStableChildBridgeNodes(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &materializationContext,
    const LocalRecompositionPlan &recompositionPlan) {
    (void) tree;
    (void) materializationContext;
    (void) recompositionPlan;
    // Stable-child gaps are canonical by construction: the local materializer
    // keeps the existing parent-child edge and does not create an empty bridge.
    return 0;
}

inline int countMcpcPreflightHiddenExternalParentConnectors(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &materializationContext,
    const LocalRecompositionPlan &recompositionPlan) {
    int count = 0;
    for (const auto &nodePlan : recompositionPlan.nodes) {
        if (nodePlan.parentPlanIndex >= 0) {
            continue;
        }

        const NodeId externalParentId = externalAttachmentForRecompositionNode(
            tree,
            materializationContext,
            nodePlan);
        bool predicted = false;
        if (nodePlan.oracleNodeId != InvalidNode &&
            oracleParentIsOutsideLocalSnapshot(
                materializationContext,
                nodePlan.oracleParentId)) {
            predicted =
                externalParentId != InvalidNode &&
                tree.isAlive(externalParentId) &&
                hiddenExternalParentConnectorLevel(
                    static_cast<int>(nodePlan.level),
                    tree.getAltitude(externalParentId),
                    tree.isMaxtree()) >= 0;
        } else {
            predicted = shouldMaterializeNonOracleHiddenExternalParentConnector(
                tree,
                materializationContext,
                nodePlan,
                externalParentId);
        }
        if (predicted) {
            ++count;
        }
    }
    return count;
}

inline bool mcpcPreflightPlanNodeMatchesPreSnapshot(
    const LocalRecompositionNodePlan &nodePlan,
    const LocalDualNodeSnapshot &snapshot) {
    if (nodePlan.oracleNodeId != InvalidNode ||
        nodePlan.level != static_cast<AltitudeType>(snapshot.altitude)) {
        return false;
    }
    auto nodeSupport = sortedUniquePixels(nodePlan.supportPixels);
    auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
    return equalPixelSets(std::span<const PixelId>(nodeSupport),
                          std::span<const PixelId>(snapshotSupport));
}

inline bool mcpcPreflightPlanHasNodeForPreSnapshot(
    const LocalRecompositionPlan &recompositionPlan,
    const LocalDualNodeSnapshot &snapshot) {
    for (const auto &nodePlan : recompositionPlan.nodes) {
        if (mcpcPreflightPlanNodeMatchesPreSnapshot(nodePlan, snapshot)) {
            return true;
        }
    }
    return false;
}

inline bool mcpcPreflightPlanHasStableChild(const LocalRecompositionNodePlan &nodePlan,
                                            NodeId nodeId) {
    return std::find(nodePlan.stableChildPreNodes.begin(),
                     nodePlan.stableChildPreNodes.end(),
                     nodeId) != nodePlan.stableChildPreNodes.end();
}

inline bool mcpcPreflightPredictsPreSnapshotRepresentation(
    const DynamicComponentTree &tree,
    const LocalRecompositionPlan &recompositionPlan,
    const LocalDualNodeSnapshot &snapshot,
    bool *representedByLiveSnapshot = nullptr) {
    if (representedByLiveSnapshot != nullptr) {
        *representedByLiveSnapshot = false;
    }
    if (mcpcPreflightPlanHasNodeForPreSnapshot(recompositionPlan, snapshot)) {
        return true;
    }
    if (snapshot.nodeId != InvalidNode &&
        tree.isNode(snapshot.nodeId) &&
        tree.isAlive(snapshot.nodeId) &&
        tree.getAltitude(snapshot.nodeId) == snapshot.altitude) {
        if (representedByLiveSnapshot != nullptr) {
            *representedByLiveSnapshot = true;
        }
        return true;
    }
    return false;
}

inline int countMcpcPreflightNonOracleHiddenSupportCandidates(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &materializationContext,
    const LocalRecompositionPlan &recompositionPlan) {
    std::vector<NodeId> countedBranchRoots;
    int count = 0;
    for (const auto &nodePlan : recompositionPlan.nodes) {
        if (nodePlan.oracleNodeId != InvalidNode ||
            nodePlan.supportPixels.empty()) {
            continue;
        }

        auto targetSupport = sortedUniquePixels(nodePlan.supportPixels);
        for (const auto &snapshot : materializationContext.preDualNodes) {
            const NodeId branchRootId = snapshot.nodeId;
            if (std::find(countedBranchRoots.begin(),
                          countedBranchRoots.end(),
                          branchRootId) != countedBranchRoots.end() ||
                branchRootId == InvalidNode ||
                snapshot.parentInSnapshot ||
                snapshot.supportLocalPixels.empty() ||
                snapshot.supportTouchesModified ||
                !snapshot.supportTouchesBoundary ||
                !tree.isAlive(branchRootId) ||
                tree.isRoot(branchRootId)) {
                continue;
            }

            auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
            if (equalPixelSets(std::span<const PixelId>(snapshotSupport),
                               std::span<const PixelId>(targetSupport)) ||
                !containsAllPixels(std::span<const PixelId>(targetSupport),
                                   std::span<const PixelId>(snapshotSupport))) {
                continue;
            }

            if (hiddenConnectorLevelBetween(
                    static_cast<int>(nodePlan.level),
                    snapshot.altitude,
                    tree.isMaxtree()) < 0) {
                continue;
            }

            countedBranchRoots.push_back(branchRootId);
            ++count;
        }
    }
    return count;
}

inline NonIncreasingMcpcDirectV1PrecheckSummary summarizeMcpcDirectV1Precheck(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &materializationContext,
    const LocalRecompositionPlan &recompositionPlan) {
    NonIncreasingMcpcDirectV1PrecheckSummary summary;

    for (const auto &nodePlan : recompositionPlan.nodes) {
        if (nodePlan.parentPlanIndex >= 0) {
            const auto parentIndex = static_cast<std::size_t>(nodePlan.parentPlanIndex);
            if (parentIndex < recompositionPlan.nodes.size() &&
                recompositionPlan.nodes[parentIndex].level == nodePlan.level) {
                ++summary.sameLevelPlannedParentChildCount;
            }
            continue;
        }

        const NodeId externalParent =
            externalAttachmentForRecompositionNode(tree, materializationContext, nodePlan);
        if (externalParent != InvalidNode &&
            tree.isAlive(externalParent) &&
            tree.getAltitude(externalParent) == static_cast<int>(nodePlan.level)) {
            ++summary.sameLevelExternalAttachmentCount;
        }
    }

    std::vector<NodeId> countedHiddenSupportBranchRoots;
    for (const auto &nodePlan : recompositionPlan.nodes) {
        if (nodePlan.oracleNodeId != InvalidNode ||
            nodePlan.supportPixels.empty()) {
            continue;
        }

        auto targetSupport = sortedUniquePixels(nodePlan.supportPixels);
        for (const auto &snapshot : materializationContext.preDualNodes) {
            const NodeId branchRootId = snapshot.nodeId;
            if (std::find(countedHiddenSupportBranchRoots.begin(),
                          countedHiddenSupportBranchRoots.end(),
                          branchRootId) != countedHiddenSupportBranchRoots.end() ||
                branchRootId == InvalidNode ||
                snapshot.parentInSnapshot ||
                snapshot.supportLocalPixels.empty() ||
                snapshot.supportTouchesModified ||
                !snapshot.supportTouchesBoundary ||
                !tree.isAlive(branchRootId) ||
                tree.isRoot(branchRootId)) {
                continue;
            }

            auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
            if (equalPixelSets(std::span<const PixelId>(snapshotSupport),
                               std::span<const PixelId>(targetSupport)) ||
                !containsAllPixels(std::span<const PixelId>(targetSupport),
                                   std::span<const PixelId>(snapshotSupport)) ||
                hiddenConnectorLevelBetween(
                    static_cast<int>(nodePlan.level),
                    snapshot.altitude,
                    tree.isMaxtree()) < 0) {
                continue;
            }

            countedHiddenSupportBranchRoots.push_back(branchRootId);
            ++summary.nonOracleHiddenSupportCandidateCount;
            if (mcpcPreflightPlanHasStableChild(nodePlan, branchRootId)) {
                ++summary.nonOracleHiddenSupportNotRedundantCount;
            } else {
                ++summary.nonOracleHiddenSupportAlreadyRepresentedCount;
            }
        }
    }

    const std::span<const LocalDualNodeSnapshot> preSnapshots(
        materializationContext.preDualNodes);
    for (const auto &childSnapshot : materializationContext.preDualNodes) {
        if (!childSnapshot.parentInSnapshot ||
            childSnapshot.parentId == InvalidNode ||
            childSnapshot.parentId == childSnapshot.nodeId) {
            continue;
        }

        const auto *parentSnapshot = findLocalSnapshotByNodeId(
            preSnapshots,
            childSnapshot.parentId);
        if (parentSnapshot == nullptr ||
            !preSnapshotParentLinkTouchesNonOraclePlan(
                recompositionPlan,
                childSnapshot,
                *parentSnapshot)) {
            continue;
        }
        ++summary.nonOracleParentLinkTouchesPlanCount;

        bool parentRepresentedByLiveSnapshot = false;
        bool childRepresentedByLiveSnapshot = false;
        const bool parentRepresented =
            mcpcPreflightPredictsPreSnapshotRepresentation(
                tree,
                recompositionPlan,
                *parentSnapshot,
                &parentRepresentedByLiveSnapshot);
        const bool childRepresented =
            mcpcPreflightPredictsPreSnapshotRepresentation(
                tree,
                recompositionPlan,
                childSnapshot,
                &childRepresentedByLiveSnapshot);
        if (!parentRepresented || !childRepresented) {
            ++summary.nonOracleParentLinkNotBothRepresentedCount;
            continue;
        }

        ++summary.nonOracleParentLinkBothRepresentedCount;
        if ((parentSnapshot->supportTouchesModified && parentRepresentedByLiveSnapshot) ||
            (childSnapshot.supportTouchesModified && childRepresentedByLiveSnapshot)) {
            ++summary.nonOracleParentLinkBlockedByModifiedSnapshotCount;
            continue;
        }

        if (!validLocalParentLevel(static_cast<AltitudeType>(childSnapshot.altitude),
                                   static_cast<AltitudeType>(parentSnapshot->altitude),
                                   tree.isMaxtree())) {
            ++summary.nonOracleParentLinkInvalidLevelCount;
        }
    }

    return summary;
}

inline bool evaluateMcpcDirectV1PrecheckAcceptance(
    const NonIncreasingMcpcDirectPreflightResult &result,
    std::string *rejectionSummary = nullptr) {
    std::vector<std::string> reasons;
    if (!result.passed) {
        reasons.push_back("base preflight did not pass");
    }
    if (result.directV1Precheck.predictedAuxiliaryNodeCount != 0) {
        reasons.push_back("predicted auxiliary nodes require canonical construction");
    }
    if (result.directV1Precheck.sameLevelPlannedParentChildCount != 0) {
        reasons.push_back("same-level planned parent-child relation requires canonical construction");
    }
    if (result.directV1Precheck.sameLevelExternalAttachmentCount != 0) {
        reasons.push_back("same-level external attachment requires canonical construction");
    }
    if (result.directV1Precheck.nonOracleHiddenSupportNotRedundantCount != 0) {
        reasons.push_back("non-oracle hidden-support candidate was not redundant");
    }
    if (result.directV1Precheck.nonOracleParentLinkNotBothRepresentedCount != 0) {
        reasons.push_back("non-oracle parent link touched plan without both endpoints represented");
    }
    if (result.directV1Precheck.nonOracleParentLinkInvalidLevelCount != 0) {
        reasons.push_back("invalid non-oracle parent-link level");
    }

    std::ostringstream stream;
    for (std::size_t i = 0; i < reasons.size(); ++i) {
        if (i != 0) {
            stream << "; ";
        }
        stream << reasons[i];
    }
    if (rejectionSummary != nullptr) {
        *rejectionSummary = stream.str();
    }
    return reasons.empty();
}

inline int countMcpcPreflightNonOracleParentLinkBridgeCandidates(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &materializationContext,
    const LocalRecompositionPlan &recompositionPlan) {
    const std::span<const LocalDualNodeSnapshot> preSnapshots(
        materializationContext.preDualNodes);
    int count = 0;
    for (const auto &childSnapshot : materializationContext.preDualNodes) {
        if (!childSnapshot.parentInSnapshot ||
            childSnapshot.parentId == InvalidNode ||
            childSnapshot.parentId == childSnapshot.nodeId) {
            continue;
        }

        const auto *parentSnapshot = findLocalSnapshotByNodeId(
            preSnapshots,
            childSnapshot.parentId);
        if (parentSnapshot == nullptr ||
            parentSnapshot->supportTouchesModified ||
            !parentSnapshot->directLocalProperParts.empty() ||
            !preSnapshotParentLinkTouchesNonOraclePlan(
                recompositionPlan,
                childSnapshot,
                *parentSnapshot)) {
            continue;
        }

        auto parentSupport = sortedUniquePixels(parentSnapshot->supportLocalPixels);
        auto childSupport = sortedUniquePixels(childSnapshot.supportLocalPixels);
        if (parentSupport.empty() ||
            childSupport.empty() ||
            !containsAllPixels(std::span<const PixelId>(parentSupport),
                               std::span<const PixelId>(childSupport))) {
            continue;
        }

        const bool parentPlanned = mcpcPreflightPlanHasNodeForPreSnapshot(
            recompositionPlan,
            *parentSnapshot);
        const bool childPlanned = mcpcPreflightPlanHasNodeForPreSnapshot(
            recompositionPlan,
            childSnapshot);
        if (parentPlanned ||
            (childSnapshot.supportTouchesModified && !childPlanned) ||
            !validLocalParentLevel(
                static_cast<AltitudeType>(childSnapshot.altitude),
                static_cast<AltitudeType>(parentSnapshot->altitude),
                tree.isMaxtree())) {
            continue;
        }

        ++count;
    }
    return count;
}

inline int countMcpcPreflightNonOracleChildLinkBridgeCandidates(
    const DynamicComponentTree &tree,
    const LocalMergeOracle &materializationContext,
    const LocalRecompositionPlan &recompositionPlan) {
    const std::span<const LocalDualNodeSnapshot> preSnapshots(
        materializationContext.preDualNodes);
    int count = 0;
    for (const auto &childSnapshot : materializationContext.preDualNodes) {
        if (!childSnapshot.parentInSnapshot ||
            childSnapshot.parentId == InvalidNode ||
            childSnapshot.parentId == childSnapshot.nodeId ||
            !childSnapshot.directLocalProperParts.empty() ||
            childSnapshot.supportTouchesModified ||
            childSnapshot.childrenInSnapshot.empty() ||
            childSnapshot.supportLocalPixels.empty() ||
            mcpcPreflightPlanHasNodeForPreSnapshot(
                recompositionPlan,
                childSnapshot)) {
            continue;
        }

        const auto *parentSnapshot = findLocalSnapshotByNodeId(
            preSnapshots,
            childSnapshot.parentId);
        if (parentSnapshot == nullptr ||
            !preSnapshotParentLinkTouchesNonOraclePlan(
                recompositionPlan,
                childSnapshot,
                *parentSnapshot) ||
            !validLocalParentLevel(
                static_cast<AltitudeType>(childSnapshot.altitude),
                static_cast<AltitudeType>(parentSnapshot->altitude),
                tree.isMaxtree())) {
            continue;
        }

        auto childSupport = sortedUniquePixels(childSnapshot.supportLocalPixels);
        bool hasPlannedDescendant = false;
        for (NodeId descendantSnapshotId : childSnapshot.childrenInSnapshot) {
            const auto *descendantSnapshot = findLocalSnapshotByNodeId(
                preSnapshots,
                descendantSnapshotId);
            if (descendantSnapshot == nullptr) {
                continue;
            }

            auto descendantSupport =
                sortedUniquePixels(descendantSnapshot->supportLocalPixels);
            if (!descendantSupport.empty() &&
                containsAllPixels(std::span<const PixelId>(childSupport),
                                  std::span<const PixelId>(descendantSupport)) &&
                mcpcPreflightPlanHasNodeForPreSnapshot(
                    recompositionPlan,
                    *descendantSnapshot) &&
                validLocalParentLevel(
                    static_cast<AltitudeType>(descendantSnapshot->altitude),
                    static_cast<AltitudeType>(childSnapshot.altitude),
                    tree.isMaxtree())) {
                hasPlannedDescendant = true;
                break;
            }
        }
        if (!hasPlannedDescendant) {
            continue;
        }

        ++count;
    }
    return count;
}

inline NonIncreasingMcpcDirectMutationSummary summarizeMcpcDirectMutationPlan(
    const DynamicComponentTree &dualTree,
    const LocalPartitionInput &input,
    const LocalPartitionPlan &partitionPlan,
    const LocalRecompositionPlan &recompositionPlan,
    const LocalMergeOracle &materializationContext) {
    NonIncreasingMcpcDirectMutationSummary summary;
    summary.affectedLevelCount = static_cast<int>(recompositionPlan.affectedLevels.size());
    summary.localPixelCount = static_cast<int>(input.localPixels.size());
    summary.preLocalNodeCount = static_cast<int>(input.preDualNodes.size());
    summary.partitionRecomposedComponentCount =
        static_cast<int>(partitionPlan.recomposedComponents.size());
    summary.plannedNodeCount = static_cast<int>(recompositionPlan.nodes.size());
    summary.plannedNodeLevels.reserve(recompositionPlan.nodes.size());

    std::map<NodeId, int> movedProperPartsBySourceOwner;
    std::vector<NodeId> stableChildRefs;
    std::map<NodeId, int> stableChildrenByPreviousParent;
    std::vector<NodeId> previousParents;

    for (const auto &nodePlan : recompositionPlan.nodes) {
        summary.plannedNodeLevels.push_back(static_cast<int>(nodePlan.level));
        summary.plannedDirectProperPartCount +=
            static_cast<int>(nodePlan.directProperParts.size());
        summary.plannedStableChildRefCount +=
            static_cast<int>(nodePlan.stableChildPreNodes.size());
        if (nodePlan.parentPlanIndex < 0) {
            ++summary.plannedExternalAttachmentCount;
        }

        for (PixelId pixelId : nodePlan.directProperParts) {
            if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
                continue;
            }
            const NodeId sourceOwner = dualTree.getSmallestComponent(pixelId);
            if (sourceOwner != InvalidNode && dualTree.isAlive(sourceOwner)) {
                ++movedProperPartsBySourceOwner[sourceOwner];
            }
        }

        for (NodeId stableChildId : nodePlan.stableChildPreNodes) {
            if (stableChildId == InvalidNode || !dualTree.isAlive(stableChildId) ||
                dualTree.isRoot(stableChildId)) {
                ++summary.invalidStableChildRefCount;
                continue;
            }
            stableChildRefs.push_back(stableChildId);
            const NodeId previousParent = dualTree.getNodeParent(stableChildId);
            if (previousParent != InvalidNode && dualTree.isAlive(previousParent)) {
                previousParents.push_back(previousParent);
                ++stableChildrenByPreviousParent[previousParent];
            }
        }
    }

    for (const auto &[sourceOwner, movedCount] : movedProperPartsBySourceOwner) {
        ++summary.sourceOwnerTouchedCount;
        if (dualTree.getNumProperParts(sourceOwner) == movedCount) {
            ++summary.sourceOwnerToBeEmptiedCount;
        }
    }

    stableChildRefs = mcpcPreflightSortedUniqueNodes(std::move(stableChildRefs));
    previousParents = mcpcPreflightSortedUniqueNodes(std::move(previousParents));
    summary.stableChildToMoveCount = static_cast<int>(stableChildRefs.size());
    summary.previousParentTouchedCount = static_cast<int>(previousParents.size());
    for (NodeId previousParent : previousParents) {
        const int movedChildren = stableChildrenByPreviousParent[previousParent];
        if (dualTree.getNumProperParts(previousParent) == 0 &&
            dualTree.getNumChildren(previousParent) == movedChildren) {
            ++summary.previousParentPossiblyEmptiedCount;
        }
    }

    for (NodeId stableChildId : stableChildRefs) {
        for (NodeId previousParent : previousParents) {
            if (stableChildId != previousParent &&
                mcpcPreflightNodeIsAncestorOf(dualTree, stableChildId, previousParent)) {
                ++summary.stableChildCycleRiskCount;
            }
        }
    }

    std::sort(summary.plannedNodeLevels.begin(), summary.plannedNodeLevels.end());
    summary.predictedHiddenExternalParentConnectorNodeCount =
        countMcpcPreflightHiddenExternalParentConnectors(
            dualTree,
            materializationContext,
            recompositionPlan);
    summary.predictedStableChildBridgeNodeCount =
        countMcpcPreflightStableChildBridgeNodes(
            dualTree,
            materializationContext,
            recompositionPlan);
    summary.predictedAuxiliaryNodeCount =
        summary.predictedHiddenExternalParentConnectorNodeCount +
        summary.predictedStableChildBridgeNodeCount;
    summary.candidateNonOracleHiddenSupportConnectorCount =
        countMcpcPreflightNonOracleHiddenSupportCandidates(
            dualTree,
            materializationContext,
            recompositionPlan);
    summary.candidateNonOracleParentLinkBridgeCount =
        countMcpcPreflightNonOracleParentLinkBridgeCandidates(
            dualTree,
            materializationContext,
            recompositionPlan);
    summary.candidateNonOracleChildLinkBridgeCount =
        countMcpcPreflightNonOracleChildLinkBridgeCandidates(
            dualTree,
            materializationContext,
            recompositionPlan);
    return summary;
}

inline NonIncreasingMcpcDirectV1AcceptanceDiagnostics diagnoseMcpcDirectPreflightV1Acceptance(
    const NonIncreasingMcpcDirectPreflightResult &result,
    const NonIncreasingStrippedLocalUpdateResult &strippedUpdate) {
    NonIncreasingMcpcDirectV1AcceptanceDiagnostics diagnostics;
    diagnostics.evaluated = true;
    std::vector<std::string> reasons;
    if (!result.passed) {
        reasons.push_back("base preflight did not pass");
    }
    if (!strippedUpdate.applied) {
        reasons.push_back("stripped local update did not apply");
    }
    if (strippedUpdate.localRecompositionTargetChecked &&
        !strippedUpdate.localRecompositionMatchesTarget) {
        reasons.push_back("stripped local update did not match target image");
    }
    if (strippedUpdate.finalReconstructionChecked &&
        !strippedUpdate.dualMatchesPrimalAfterMerge) {
        reasons.push_back("stripped local update did not preserve dual/primal reconstruction");
    }
    if (strippedUpdate.createdConnectorNodeCount != 0) {
        reasons.push_back("stripped local update created auxiliary connectors");
    }
    if (strippedUpdate.postCanonicalizationChangesRequired != 0) {
        reasons.push_back("stripped local update required corrective post-canonicalization");
    }
    if (strippedUpdate.nonOracleHiddenSupportMaterializableCount != 0) {
        reasons.push_back("materializable non-oracle hidden-support connector");
    }
    if (strippedUpdate.nonOracleParentLinkBridgeMaterializableCount != 0) {
        reasons.push_back("materializable non-oracle parent-link bridge");
    }
    if (strippedUpdate.nonOracleChildLinkBridgeMaterializableCount != 0) {
        reasons.push_back("materializable non-oracle child-link bridge");
    }
    if (strippedUpdate.nonOracleParentLinkInvalidLevelCount != 0) {
        reasons.push_back("invalid non-oracle parent-link level");
    }
    if (strippedUpdate.nonOracleChildLinkBridgeInvalidLevelCount != 0) {
        reasons.push_back("invalid non-oracle child-link level");
    }
    if (strippedUpdate.nonOracleParentLinkTouchesPlanCount !=
        strippedUpdate.nonOracleParentLinkBothRepresentedCount) {
        reasons.push_back("non-oracle parent link touched plan without both endpoints represented");
    }
    if (strippedUpdate.nonOracleHiddenSupportCandidateCount !=
        strippedUpdate.nonOracleHiddenSupportAlreadyRepresentedCount) {
        reasons.push_back("non-oracle hidden-support candidate was not redundant");
    }

    std::ostringstream stream;
    for (std::size_t i = 0; i < reasons.size(); ++i) {
        if (i != 0) {
            stream << "; ";
        }
        stream << reasons[i];
    }
    diagnostics.rejectionSummary = stream.str();
    diagnostics.accepted = reasons.empty();
    return diagnostics;
}

inline NonIncreasingMcpcDirectSummaryComparison compareMcpcDirectPreflightSummaryWithStrippedUpdate(
    const NonIncreasingMcpcDirectMutationSummary &summary,
    const NonIncreasingStrippedLocalUpdateResult &strippedUpdate) {
    NonIncreasingMcpcDirectSummaryComparison comparison;
    comparison.evaluated = true;
    comparison.predictedEmptiedLocalNodeCount =
        summary.sourceOwnerToBeEmptiedCount +
        summary.previousParentPossiblyEmptiedCount;
    comparison.observedEmptiedLocalNodeCount =
        strippedUpdate.emptiedSourceOwnerCount;

    std::vector<std::string> mismatches;
    const auto checkCount = [&](const char *label,
                                int expected,
                                int observed,
                                int &flag) {
        if (expected == observed) {
            return;
        }
        flag = 1;
        ++comparison.mismatchCount;
        std::ostringstream out;
        out << label << " expected=" << expected << " observed=" << observed;
        mismatches.push_back(out.str());
    };

    checkCount("planned_nodes",
               summary.plannedNodeCount,
               strippedUpdate.createdPlannedNodeCount,
               comparison.plannedNodeCountMismatch);
    checkCount("planned_direct_proper_parts",
               summary.plannedDirectProperPartCount,
               strippedUpdate.movedDirectProperPartCount,
               comparison.plannedDirectProperPartCountMismatch);
    checkCount("predicted_auxiliary_nodes",
               summary.predictedAuxiliaryNodeCount,
               strippedUpdate.createdConnectorNodeCount,
               comparison.predictedAuxiliaryNodeCountMismatch);
    checkCount("predicted_hidden_external_parent_connectors",
               summary.predictedHiddenExternalParentConnectorNodeCount,
               strippedUpdate.createdHiddenExternalParentConnectorNodeCount,
               comparison.predictedHiddenExternalParentConnectorMismatch);
    checkCount("predicted_stable_child_bridges",
               summary.predictedStableChildBridgeNodeCount,
               strippedUpdate.createdStableChildBridgeNodeCount,
               comparison.predictedStableChildBridgeMismatch);
    checkCount("emptied_local_nodes",
               comparison.predictedEmptiedLocalNodeCount,
               comparison.observedEmptiedLocalNodeCount,
               comparison.emptiedLocalNodeCountMismatch);
    checkCount("stable_children_to_move",
               summary.stableChildToMoveCount,
               strippedUpdate.movedStableChildNodeCount,
               comparison.movedStableChildCountMismatch);

    comparison.matches = comparison.mismatchCount == 0;
    for (std::size_t i = 0; i < mismatches.size(); ++i) {
        if (i != 0) {
            comparison.mismatchSummary += "; ";
        }
        comparison.mismatchSummary += mismatches[i];
    }
    return comparison;
}

inline NonIncreasingMcpcDirectCommitCandidateDiagnostics diagnoseMcpcDirectCommitCandidate(
    const NonIncreasingMcpcDirectPreflightResult &preflight,
    const NonIncreasingMcpcDirectSummaryComparison &summaryComparison,
    const NonIncreasingMcpcDirectV1AcceptanceDiagnostics &directV1,
    const NonIncreasingStrippedLocalUpdateResult &strippedUpdate) {
    NonIncreasingMcpcDirectCommitCandidateDiagnostics diagnostics;
    diagnostics.evaluated = true;

    std::vector<std::string> reasons;
    if (!preflight.checked || !preflight.passed) {
        reasons.push_back("preflight did not pass");
    }
    if (!preflight.planOracleAuditChecked || !preflight.planOracleFree) {
        reasons.push_back("plan is not certified oracle-free");
    }
    if (!preflight.materializationOracleAuditChecked ||
        !preflight.materializationOracleFree) {
        reasons.push_back("materialization input is not certified oracle-free");
    }
    if (!preflight.materializationInputsValid) {
        reasons.push_back("materialization inputs were not validated");
    }
    if (!preflight.rollbackFreeSummaryBuilt) {
        reasons.push_back("rollback-free summary was not built");
    }
    if (!preflight.directV1PrecheckEvaluated ||
        !preflight.directV1PrecheckAccepted) {
        reasons.push_back("direct-v1 precheck did not accept");
    }
    if (!directV1.evaluated || !directV1.accepted) {
        reasons.push_back("direct-v1 copied diagnostics did not accept");
    }
    if (!summaryComparison.evaluated || !summaryComparison.matches) {
        reasons.push_back("direct summary does not match copied materialization");
    }
    if (!strippedUpdate.applied || !strippedUpdate.localRecompositionApplied) {
        reasons.push_back("copied local recomposition was not applied");
    }
    if (!strippedUpdate.localRecompositionTargetChecked ||
        !strippedUpdate.localRecompositionMatchesTarget) {
        reasons.push_back("copied local recomposition did not verify target reconstruction");
    }
    if (!strippedUpdate.finalReconstructionChecked ||
        !strippedUpdate.dualMatchesPrimalAfterMerge) {
        reasons.push_back("copied trial did not verify final dual/primal reconstruction");
    }
    if (!strippedUpdate.localPlanOracleAuditChecked ||
        !strippedUpdate.localPlanOracleFree) {
        reasons.push_back("copied plan audit did not certify oracle-free plan");
    }
    if (!strippedUpdate.materializationInputOracleAuditChecked ||
        !strippedUpdate.materializationInputOracleFree) {
        reasons.push_back("copied materialization audit did not certify oracle-free input");
    }
    if (!strippedUpdate.emptyCreatedNodesAfterCanonicalization.empty()) {
        reasons.push_back("copied materialization left created nodes empty");
    }
    if (strippedUpdate.postCanonicalizationChangesRequired != 0) {
        reasons.push_back("copied materialization required corrective post-canonicalization");
    }
    if (preflight.summary.invalidStableChildRefCount != 0) {
        reasons.push_back("preflight reported invalid stable child refs");
    }
    if (preflight.summary.stableChildCycleRiskCount != 0) {
        reasons.push_back("preflight reported stable child cycle risks");
    }

    std::ostringstream stream;
    for (std::size_t i = 0; i < reasons.size(); ++i) {
        if (i != 0) {
            stream << "; ";
        }
        stream << reasons[i];
    }
    diagnostics.rejectionSummary = stream.str();
    diagnostics.accepted = reasons.empty();
    return diagnostics;
}

inline NonIncreasingMcpcNoCopyCandidateDiagnostics diagnoseMcpcNoCopyCandidate(
    const NonIncreasingMcpcDirectPreflightResult &preflight) {
    NonIncreasingMcpcNoCopyCandidateDiagnostics diagnostics;
    diagnostics.evaluated = true;

    std::vector<std::string> reasons;
    if (!preflight.checked || !preflight.passed) {
        reasons.push_back("preflight did not pass");
    }
    if (!preflight.inputBuilt ||
        !preflight.partitionPlanBuilt ||
        !preflight.recompositionPlanBuilt) {
        reasons.push_back("local recomposition plan was not built");
    }
    if (!preflight.planOracleAuditChecked || !preflight.planOracleFree) {
        reasons.push_back("plan is not certified oracle-free");
    }
    if (!preflight.materializationContextBuilt) {
        reasons.push_back("materialization context was not built");
    }
    if (!preflight.materializationOracleAuditChecked ||
        !preflight.materializationOracleFree) {
        reasons.push_back("materialization input is not certified oracle-free");
    }
    if (!preflight.materializationInputsValid) {
        reasons.push_back("materialization inputs were not validated");
    }
    if (!preflight.rollbackFreeSummaryBuilt) {
        reasons.push_back("rollback-free summary was not built");
    }
    if (!preflight.directV1PrecheckEvaluated ||
        !preflight.directV1PrecheckAccepted) {
        reasons.push_back("direct-v1 precheck did not accept");
    }
    if (preflight.summary.invalidStableChildRefCount != 0) {
        reasons.push_back("preflight reported invalid stable child refs");
    }
    if (preflight.summary.stableChildCycleRiskCount != 0) {
        reasons.push_back("preflight reported stable child cycle risks");
    }

    std::ostringstream stream;
    for (std::size_t i = 0; i < reasons.size(); ++i) {
        if (i != 0) {
            stream << "; ";
        }
        stream << reasons[i];
    }
    diagnostics.rejectionSummary = stream.str();
    diagnostics.accepted = reasons.empty();
    return diagnostics;
}

inline NonIncreasingMcpcD12DomainDiagnostics diagnoseMcpcD12Domain(
    const NonIncreasingMcpcDirectPreflightResult &preflight) {
    NonIncreasingMcpcD12DomainDiagnostics diagnostics;
    diagnostics.evaluated = true;

    const auto reject = [&](const std::string &code,
                            const std::string &summary) {
        diagnostics.accepted = false;
        diagnostics.rejectReasonCode = code;
        diagnostics.rejectionSummary = summary;
        return diagnostics;
    };

    if (!preflight.checked) {
        return reject("d12_base_preflight_failed",
                      "D12 requires a checked P1-MCPC direct preflight");
    }
    if (!preflight.inputBuilt ||
        !preflight.partitionPlanBuilt ||
        !preflight.recompositionPlanBuilt) {
        return reject("d12_plan_incomplete",
                      "D12 requires input, partition plan, and recomposition plan");
    }
    if (!preflight.planOracleAuditChecked ||
        !preflight.planOracleFree) {
        return reject("d12_plan_oracle_not_free",
                      "D12 requires an oracle-free recomposition plan");
    }
    if (!preflight.materializationContextBuilt) {
        return reject("d12_materialization_context_missing",
                      "D12 requires a read-only materialization context");
    }
    if (!preflight.materializationOracleAuditChecked ||
        !preflight.materializationOracleFree) {
        return reject("d12_materialization_oracle_not_free",
                      "D12 requires oracle-free materialization input");
    }
    if (!preflight.materializationInputsValid) {
        return reject("d12_materialization_inputs_invalid",
                      "D12 requires valid materialization inputs");
    }
    if (!preflight.rollbackFreeSummaryBuilt) {
        return reject("d12_rollback_summary_missing",
                      "D12 requires a rollback-free mutation summary");
    }
    if (preflight.summary.invalidStableChildRefCount != 0) {
        return reject("d12_invalid_stable_child_ref",
                      "D12 rejects invalid stable-child references");
    }
    if (preflight.summary.stableChildCycleRiskCount != 0) {
        return reject("d12_stable_child_cycle_risk",
                      "D12 rejects stable-child cycle risks");
    }
    if (!preflight.passed) {
        return reject("d12_base_preflight_failed",
                      preflight.failureSummary.empty()
                          ? "D12 base preflight did not pass"
                          : preflight.failureSummary);
    }
    if (!preflight.directV1PrecheckEvaluated) {
        return reject("d12_direct_v1_precheck_missing",
                      "D12 requires direct-v1 precheck diagnostics");
    }
    if (preflight.directV1Precheck.predictedAuxiliaryNodeCount != 0) {
        return reject("d12_auxiliary_nodes_required",
                      "D12 rejects events that require auxiliary connector nodes");
    }
    if (preflight.directV1Precheck.sameLevelPlannedParentChildCount != 0) {
        return reject("d12_same_level_planned_parent_child",
                      "D12 rejects planned same-level parent-child relations");
    }
    if (preflight.directV1Precheck.sameLevelExternalAttachmentCount != 0) {
        return reject("d12_same_level_external_attachment",
                      "D12 rejects same-level external attachments");
    }
    if (preflight.directV1Precheck.nonOracleHiddenSupportNotRedundantCount != 0) {
        return reject("d12_hidden_support_not_redundant",
                      "D12 rejects non-oracle hidden support that is not redundant");
    }
    if (preflight.directV1Precheck.nonOracleParentLinkNotBothRepresentedCount != 0) {
        return reject("d12_parent_link_not_represented",
                      "D12 rejects parent links without both endpoints represented");
    }
    if (preflight.directV1Precheck.nonOracleParentLinkInvalidLevelCount != 0) {
        return reject("d12_parent_link_invalid_level",
                      "D12 rejects invalid non-oracle parent-link levels");
    }
    if (!preflight.directV1PrecheckAccepted) {
        return reject("d12_direct_v1_precheck_rejected",
                      preflight.directV1PrecheckRejectionSummary.empty()
                          ? "D12 direct-v1 precheck rejected the event"
                          : preflight.directV1PrecheckRejectionSummary);
    }

    return reject("d12_exception_free_unproven",
                  "D12 exception-free live mutation lemma is not proved for P1-MCPC");
}

inline NodeId mcpcDryRunExternalAttachmentForPlanIndex(
    const DynamicComponentTree &dualTree,
    const NonIncreasingMcpcDirectPreflightArtifacts &artifacts,
    std::size_t planIndex) {
    if (planIndex >= artifacts.recompositionPlan.nodes.size()) {
        return InvalidNode;
    }

    std::vector<unsigned char> seen(artifacts.recompositionPlan.nodes.size(), 0);
    std::size_t currentIndex = planIndex;
    while (currentIndex < artifacts.recompositionPlan.nodes.size()) {
        if (seen[currentIndex] != 0) {
            return InvalidNode;
        }
        seen[currentIndex] = 1;

        const auto &nodePlan = artifacts.recompositionPlan.nodes[currentIndex];
        if (nodePlan.parentPlanIndex < 0) {
            return externalAttachmentForRecompositionNode(
                dualTree,
                artifacts.materializationContext,
                nodePlan);
        }
        const auto parentIndex = static_cast<std::size_t>(nodePlan.parentPlanIndex);
        if (parentIndex >= artifacts.recompositionPlan.nodes.size() ||
            parentIndex == currentIndex) {
            return InvalidNode;
        }
        currentIndex = parentIndex;
    }
    return InvalidNode;
}

inline bool mcpcDryRunEmptyOwnerReleaseIsSimple(
    const DynamicComponentTree &dualTree,
    const NonIncreasingMcpcDirectPreflightArtifacts &artifacts,
    std::string *rejectReasonCode = nullptr,
    std::string *rejectionSummary = nullptr) {
    const auto reject = [&](const std::string &code,
                            const std::string &summary) {
        if (rejectReasonCode != nullptr) {
            *rejectReasonCode = code;
        }
        if (rejectionSummary != nullptr) {
            *rejectionSummary = summary;
        }
        return false;
    };

    const auto &summary = artifacts.preflight.summary;
    if (summary.sourceOwnerToBeEmptiedCount == 0 &&
        summary.previousParentPossiblyEmptiedCount == 0) {
        return true;
    }

    std::map<NodeId, int> movedProperPartsBySourceOwner;
    std::map<NodeId, int> stableChildrenByPreviousParent;
    std::map<NodeId, std::vector<NodeId>> stableChildNodesByPreviousParent;
    std::vector<NodeId> previousParents;
    for (const auto &nodePlan : artifacts.recompositionPlan.nodes) {
        for (PixelId pixelId : nodePlan.directProperParts) {
            if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
                return reject(
                    "d12_dry_run_empty_owner_invalid_direct_proper_part",
                    "D12 empty-owner release found an invalid direct proper part");
            }
            const NodeId sourceOwner = dualTree.getSmallestComponent(pixelId);
            if (sourceOwner == InvalidNode || !dualTree.isAlive(sourceOwner)) {
                return reject(
                    "d12_dry_run_empty_owner_non_live_source_owner",
                    "D12 empty-owner release found a non-live source owner");
            }
            ++movedProperPartsBySourceOwner[sourceOwner];
        }

        for (NodeId stableChildId : nodePlan.stableChildPreNodes) {
            if (stableChildId == InvalidNode ||
                !dualTree.isAlive(stableChildId) ||
                dualTree.isRoot(stableChildId)) {
                return reject(
                    "d12_dry_run_empty_owner_invalid_stable_child",
                    "D12 empty-owner release found an invalid stable child");
            }
            const NodeId previousParent = dualTree.getNodeParent(stableChildId);
            if (previousParent != InvalidNode && dualTree.isAlive(previousParent)) {
                previousParents.push_back(previousParent);
                ++stableChildrenByPreviousParent[previousParent];
                stableChildNodesByPreviousParent[previousParent].push_back(stableChildId);
            }
        }
    }

    const auto childIsMovedFromSourceOwner = [&](NodeId sourceOwner, NodeId childId) {
        auto it = stableChildNodesByPreviousParent.find(sourceOwner);
        if (it == stableChildNodesByPreviousParent.end()) {
            return false;
        }
        auto children = sortedUniqueValidNodes(it->second);
        return std::find(children.begin(), children.end(), childId) != children.end();
    };

    const auto residualChildCanBePromotedToSourceParent =
        [&](NodeId sourceOwner, NodeId childId, NodeId promotedParent) {
            const auto *childSnapshot = findLocalSnapshotByNodeId(
                std::span<const LocalDualNodeSnapshot>(
                    artifacts.materializationContext.preDualNodes),
                childId);
            if (childSnapshot == nullptr ||
                childSnapshot->supportLocalPixels.empty()) {
                return reject(
                    "d12_dry_run_empty_owner_with_children_residual_child_snapshot_missing",
                    "D12 residual child promotion requires a non-empty local pre-snapshot");
            }
            if (childSnapshot->supportTouchesModified) {
                return reject(
                    "d12_dry_run_empty_owner_with_children_residual_child_touches_modified",
                    "D12 residual child promotion rejects residual children that touch modified pixels");
            }
            if (!validLocalParentLevel(
                    static_cast<AltitudeType>(dualTree.getAltitude(childId)),
                    static_cast<AltitudeType>(dualTree.getAltitude(promotedParent)),
                    dualTree.isMaxtree())) {
                return reject(
                    "d12_dry_run_empty_owner_with_children_invalid_child_level",
                    "D12 residual child promotion found an invalid promoted parent level");
            }
            if (mcpcPreflightNodeIsAncestorOf(dualTree, childId, promotedParent)) {
                return reject(
                    "d12_dry_run_empty_owner_with_children_cycle_risk",
                    "D12 residual child promotion would create a cycle");
            }

            auto childSupport =
                sortedUniquePixels(childSnapshot->supportLocalPixels);
            for (const auto &nodePlan : artifacts.recompositionPlan.nodes) {
                if (nodePlan.supportPixels.empty()) {
                    continue;
                }
                if (mcpcPreflightPixelSetsIntersect(
                        std::span<const PixelId>(nodePlan.supportPixels),
                        std::span<const PixelId>(childSupport))) {
                    return reject(
                        "d12_dry_run_empty_owner_with_children_residual_child_support_overlaps_plan",
                        "D12 residual child promotion rejects residual children whose support intersects a planned node");
                }
            }

            (void) sourceOwner;
            return true;
        };

    const auto sourceOwnerChildrenCanBeReleased = [&](NodeId sourceOwner) {
        const NodeId promotedParent = dualTree.getNodeParent(sourceOwner);
        if (promotedParent == InvalidNode ||
            promotedParent == sourceOwner ||
            !dualTree.isAlive(promotedParent)) {
            return reject(
                "d12_dry_run_empty_owner_with_children_invalid_parent",
                "D12 empty-owner with children requires a live promoted parent");
        }

        std::vector<NodeId> seenChildren;
        for (NodeId childId : dualTree.getChildren(sourceOwner)) {
            if (childId == InvalidNode ||
                !dualTree.isNode(childId) ||
                !dualTree.isAlive(childId) ||
                childId == sourceOwner ||
                childId == promotedParent ||
                dualTree.getNodeParent(childId) != sourceOwner) {
                return reject(
                    "d12_dry_run_empty_owner_with_children_invalid_child",
                    "D12 empty-owner with children found an invalid residual child");
            }
            seenChildren.push_back(childId);
            if (!childIsMovedFromSourceOwner(sourceOwner, childId) &&
                !residualChildCanBePromotedToSourceParent(
                    sourceOwner,
                    childId,
                    promotedParent)) {
                return false;
            }
        }

        seenChildren = sortedUniqueValidNodes(std::move(seenChildren));
        if (static_cast<int>(seenChildren.size()) != dualTree.getNumChildren(sourceOwner)) {
            return reject(
                "d12_dry_run_empty_owner_with_children_duplicate_child",
                "D12 empty-owner with children found duplicated residual children");
        }

        for (std::size_t planIndex = 0;
             planIndex < artifacts.recompositionPlan.nodes.size();
             ++planIndex) {
            const auto &nodePlan = artifacts.recompositionPlan.nodes[planIndex];
            if (nodePlan.parentPlanIndex >= 0) {
                continue;
            }
            NodeId externalAttachment = InvalidNode;
            try {
                externalAttachment =
                    mcpcDryRunExternalAttachmentForPlanIndex(dualTree, artifacts, planIndex);
            } catch (const std::exception &e) {
                return reject(
                    "d12_dry_run_external_attachment_invalid",
                    std::string("D12 empty-owner with children could not resolve an external attachment: ") +
                        e.what());
            }
            if (externalAttachment != sourceOwner) {
                continue;
            }
            return reject(
                "d12_dry_run_empty_owner_with_children_external_dependency",
                "D12 empty-owner with children rejects planned nodes attached through the released owner");
        }
        return true;
    };

    int sourceOwnersToRelease = 0;
    for (const auto &[sourceOwner, movedCount] : movedProperPartsBySourceOwner) {
        if (dualTree.getNumProperParts(sourceOwner) != movedCount) {
            continue;
        }
        ++sourceOwnersToRelease;
        if (dualTree.isRoot(sourceOwner)) {
            return reject(
                "d12_dry_run_empty_owner_root_source_owner",
                "D12 empty-owner release cannot release the root source owner");
        }
        if (dualTree.getNumChildren(sourceOwner) != 0) {
            if (!sourceOwnerChildrenCanBeReleased(sourceOwner)) {
                return false;
            }
        }
    }
    if (sourceOwnersToRelease != summary.sourceOwnerToBeEmptiedCount) {
        return reject(
            "d12_dry_run_empty_owner_source_count_mismatch",
            "D12 empty-owner release source-owner count does not match the summary");
    }

    previousParents = mcpcPreflightSortedUniqueNodes(std::move(previousParents));
    int previousParentsToRelease = 0;
    for (NodeId previousParent : previousParents) {
        const int movedChildren = stableChildrenByPreviousParent[previousParent];
        if (dualTree.getNumProperParts(previousParent) != 0 ||
            dualTree.getNumChildren(previousParent) != movedChildren) {
            continue;
        }
        ++previousParentsToRelease;
        if (dualTree.isRoot(previousParent)) {
            return reject(
                "d12_dry_run_empty_owner_root_previous_parent",
                "D12 empty-owner release cannot release the root previous parent");
        }
    }
    if (previousParentsToRelease != summary.previousParentPossiblyEmptiedCount) {
        return reject(
            "d12_dry_run_empty_owner_previous_parent_count_mismatch",
            "D12 empty-owner release previous-parent count does not match the summary");
    }

    return true;
}

inline NonIncreasingMcpcLiveMutationDryRun diagnoseMcpcLiveMutationDryRun(
    const DynamicComponentTree &dualTree,
    const NonIncreasingMcpcDirectPreflightArtifacts &artifacts) {
    NonIncreasingMcpcLiveMutationDryRun dryRun;
    dryRun.evaluated = true;
    const auto &preflight = artifacts.preflight;
    dryRun.plannedNodeCount = preflight.summary.plannedNodeCount;
    dryRun.plannedDirectProperPartCount =
        preflight.summary.plannedDirectProperPartCount;
    dryRun.plannedExternalAttachmentCount =
        preflight.summary.plannedExternalAttachmentCount;
    dryRun.plannedStableChildRefCount =
        preflight.summary.plannedStableChildRefCount;
    dryRun.sourceOwnerToBeEmptiedCount =
        preflight.summary.sourceOwnerToBeEmptiedCount;
    dryRun.previousParentPossiblyEmptiedCount =
        preflight.summary.previousParentPossiblyEmptiedCount;
    dryRun.stableChildToMoveCount =
        preflight.summary.stableChildToMoveCount;
    dryRun.predictedAuxiliaryNodeCount =
        preflight.summary.predictedAuxiliaryNodeCount;
    dryRun.predictedHiddenExternalParentConnectorNodeCount =
        preflight.summary.predictedHiddenExternalParentConnectorNodeCount;
    dryRun.predictedStableChildBridgeNodeCount =
        preflight.summary.predictedStableChildBridgeNodeCount;
    dryRun.candidateNonOracleHiddenSupportConnectorCount =
        preflight.summary.candidateNonOracleHiddenSupportConnectorCount;
    dryRun.candidateNonOracleParentLinkBridgeCount =
        preflight.summary.candidateNonOracleParentLinkBridgeCount;
    dryRun.candidateNonOracleChildLinkBridgeCount =
        preflight.summary.candidateNonOracleChildLinkBridgeCount;

    const auto reject = [&](const std::string &code,
                            const std::string &summary) {
        dryRun.accepted = false;
        dryRun.rejectReasonCode = code;
        dryRun.rejectionSummary = summary;
        return dryRun;
    };

    const auto d12 = diagnoseMcpcD12Domain(preflight);
    if (!d12.evaluated ||
        d12.rejectReasonCode != "d12_exception_free_unproven") {
        return reject(
            d12.rejectReasonCode.empty() ? "d12_read_only_predicate_failed"
                                         : d12.rejectReasonCode,
            d12.rejectionSummary.empty() ? "D12 read-only predicates did not pass"
                                         : d12.rejectionSummary);
    }
    if (artifacts.recompositionPlan.nodes.empty()) {
        return reject("d12_dry_run_empty_plan",
                      "D12 live dry-run requires a non-empty recomposition plan");
    }
    if (preflight.summary.candidateNonOracleHiddenSupportConnectorCount !=
            preflight.directV1Precheck.nonOracleHiddenSupportCandidateCount ||
        preflight.directV1Precheck.nonOracleHiddenSupportCandidateCount !=
            preflight.directV1Precheck.nonOracleHiddenSupportAlreadyRepresentedCount) {
        return reject("d12_dry_run_non_oracle_hidden_support_connector_present",
                      "D12 hidden-support live subdomain requires all non-oracle hidden-support candidates to be redundant");
    }
    if (preflight.summary.predictedAuxiliaryNodeCount !=
        preflight.summary.predictedHiddenExternalParentConnectorNodeCount +
            preflight.summary.predictedStableChildBridgeNodeCount) {
        return reject("d12_dry_run_auxiliary_nodes_present",
                      "D12 live dry-run found unclassified auxiliary nodes");
    }
    if (preflight.summary.plannedStableChildRefCount !=
        preflight.summary.stableChildToMoveCount) {
        return reject("d12_dry_run_stable_child_move_present",
                      "D12 stable-child live subdomain requires unique stable-child refs");
    }

    std::vector<unsigned char> claimed(
        static_cast<std::size_t>(dualTree.getNumTotalProperParts()),
        0);
    std::vector<unsigned char> claimedStableChild(
        static_cast<std::size_t>(dualTree.getNumInternalNodeSlots()),
        0);
    for (std::size_t i = 0; i < artifacts.recompositionPlan.nodes.size(); ++i) {
        const auto &nodePlan = artifacts.recompositionPlan.nodes[i];
        if (nodePlan.parentPlanIndex < 0) {
            NodeId externalParentId = InvalidNode;
            try {
                externalParentId = externalAttachmentForRecompositionNode(
                    dualTree,
                    artifacts.materializationContext,
                    nodePlan);
            } catch (const std::exception &e) {
                return reject(
                    "d12_dry_run_external_attachment_invalid",
                    std::string("D12 dry-run could not resolve an external attachment: ") +
                        e.what());
            }
            if (externalParentId == InvalidNode ||
                !dualTree.isAlive(externalParentId) ||
                !validLocalParentLevel(
                    nodePlan.level,
                    static_cast<AltitudeType>(dualTree.getAltitude(externalParentId)),
                    dualTree.isMaxtree())) {
                return reject("d12_dry_run_external_attachment_invalid",
                              "D12 dry-run found an invalid external attachment");
            }
            if (dualTree.getNumProperParts(externalParentId) == 0 &&
                dualTree.getNumChildren(externalParentId) == 0) {
                return reject("d12_dry_run_external_attachment_empty",
                              "D12 dry-run rejects empty external attachment parents");
            }
        } else {
            const auto parentIndex = static_cast<std::size_t>(nodePlan.parentPlanIndex);
            if (parentIndex >= artifacts.recompositionPlan.nodes.size() ||
                parentIndex == i) {
                return reject("d12_dry_run_invalid_parent_plan",
                              "D12 dry-run found an invalid local parent index");
            }
            const auto &parentPlan = artifacts.recompositionPlan.nodes[parentIndex];
            if (!validLocalParentLevel(nodePlan.level,
                                       parentPlan.level,
                                       dualTree.isMaxtree())) {
                return reject("d12_dry_run_invalid_parent_level",
                              "D12 dry-run found an invalid local parent level");
            }
        }
        if (nodePlan.directProperParts.empty()) {
            return reject("d12_dry_run_empty_planned_node",
                          "D12 first live subdomain requires each planned node to move direct proper parts");
        }
        if (!nodePlan.oracleChildrenInSnapshot.empty()) {
            return reject("d12_dry_run_child_move_present",
                          "D12 first live subdomain excludes child moves and child bridge materialization");
        }
        for (NodeId stableChildId : nodePlan.stableChildPreNodes) {
            if (stableChildId == InvalidNode ||
                !dualTree.isAlive(stableChildId) ||
                dualTree.isRoot(stableChildId)) {
                return reject("d12_invalid_stable_child_ref",
                              "D12 stable-child live subdomain found an invalid stable-child ref");
            }
            if (stableChildId < 0 ||
                stableChildId >= dualTree.getNumInternalNodeSlots()) {
                return reject("d12_invalid_stable_child_ref",
                              "D12 stable-child live subdomain found an out-of-range stable-child ref");
            }
            const auto stableChildIndex = static_cast<std::size_t>(stableChildId);
            if (claimedStableChild[stableChildIndex] != 0) {
                return reject("d12_dry_run_stable_child_move_present",
                              "D12 stable-child live subdomain rejects duplicated stable-child refs");
            }
            claimedStableChild[stableChildIndex] = 1;
            if (!validLocalParentLevel(
                    static_cast<AltitudeType>(dualTree.getAltitude(stableChildId)),
                    nodePlan.level,
                    dualTree.isMaxtree())) {
                return reject("d12_dry_run_stable_child_move_present",
                              "D12 stable-child live subdomain found an invalid stable-child target level");
            }

            NodeId externalAttachment = InvalidNode;
            try {
                externalAttachment =
                    mcpcDryRunExternalAttachmentForPlanIndex(dualTree, artifacts, i);
            } catch (const std::exception &e) {
                return reject(
                    "d12_dry_run_external_attachment_invalid",
                    std::string("D12 dry-run could not resolve a stable-child attachment: ") +
                        e.what());
            }
            if (externalAttachment != InvalidNode &&
                mcpcPreflightNodeIsAncestorOf(
                    dualTree,
                    stableChildId,
                    externalAttachment)) {
                return reject("d12_stable_child_cycle_risk",
                              "D12 stable-child live subdomain would attach through a descendant");
            }
        }
        for (PixelId pixelId : nodePlan.directProperParts) {
            if (pixelId < 0 || pixelId >= dualTree.getNumTotalProperParts()) {
                return reject("d12_dry_run_invalid_pixel",
                              "D12 dry-run found an out-of-bounds proper part");
            }
            const auto pixelIndex = static_cast<std::size_t>(pixelId);
            if (claimed[pixelIndex] != 0) {
                return reject("d12_dry_run_duplicated_pixel",
                              "D12 dry-run found duplicated direct proper parts");
            }
            claimed[pixelIndex] = 1;
        }
    }

    std::string emptyOwnerRejectCode;
    std::string emptyOwnerRejection;
    if (!mcpcDryRunEmptyOwnerReleaseIsSimple(
            dualTree,
            artifacts,
            &emptyOwnerRejectCode,
            &emptyOwnerRejection)) {
        return reject(emptyOwnerRejectCode.empty()
                          ? "d12_dry_run_empty_owner_present"
                          : emptyOwnerRejectCode,
                      emptyOwnerRejection.empty()
                          ? "D12 empty-owner release predicate rejected the event"
                          : emptyOwnerRejection);
    }

    dryRun.accepted = true;
    dryRun.rejectReasonCode.clear();
    dryRun.rejectionSummary =
        "D12 first live subdomain accepted by structural dry-run";
    return dryRun;
}

inline NonIncreasingMcpcDirectPreflightArtifacts buildMcpcDirectPreflightArtifacts(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId nodeId) {
    NonIncreasingMcpcDirectPreflightArtifacts artifacts;
    auto &result = artifacts.preflight;
    result.checked = true;

    try {
        if (!adjacency) {
            result.failureSummary = "P1-MCPC direct preflight requires a valid adjacency relation";
            return artifacts;
        }
        if (nodeId == InvalidNode ||
            nodeId == primalTree.getRoot() ||
            !primalTree.isNode(nodeId) ||
            !primalTree.isAlive(nodeId)) {
            result.failureSummary = "P1-MCPC direct preflight requires a live non-root primal node";
            return artifacts;
        }
        const NodeId parentId = primalTree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !primalTree.isAlive(parentId)) {
            result.failureSummary = "P1-MCPC direct preflight requires a live primal parent";
            return artifacts;
        }
        if (primalTree.getAltitude(nodeId) == primalTree.getAltitude(parentId) ||
            primalTree.getNumProperParts(nodeId) == 0) {
            result.failureSummary = "P1-MCPC direct preflight requires a non-empty level-changing merge";
            return artifacts;
        }

        LocalPartitionInputOptions inputOptions;
        inputOptions.checkPreReconstruction = true;
        inputOptions.buildFullUpdatedImage = true;
        artifacts.input = makeLocalPartitionInputFromLiveTrees(
            primalTree,
            dualTree,
            adjacency,
            nodeId,
            inputOptions);
        result.inputBuilt = true;

        artifacts.partitionPlan = buildLocalPartitionPlanWithoutOracle(
            artifacts.input,
            *adjacency);
        result.partitionPlanBuilt = true;
        artifacts.recompositionPlan = buildLocalRecompositionPlanWithoutOracle(
            artifacts.input,
            artifacts.partitionPlan);
        result.recompositionPlanBuilt = true;

        const auto planDependencyReport = auditLocalPlanOracleDependencies(
            &artifacts.partitionPlan,
            &artifacts.recompositionPlan);
        result.planOracleAuditChecked = true;
        result.planOracleFree = planDependencyReport.isOracleFree();

        artifacts.materializationContext =
            makeStrippedNonOracleMaterializationContext(dualTree, artifacts.input);
        result.materializationContextBuilt = true;

        const auto materializationDependencyReport =
            auditLocalMaterializationOracleDependencies(
                artifacts.materializationContext,
                artifacts.recompositionPlan);
        result.materializationOracleAuditChecked = true;
        result.materializationOracleFree =
            materializationDependencyReport.isOracleFreeMaterializationInput();

        LocalRecompositionMaterializationOptions materializationOptions;
        materializationOptions.checkPreReconstruction = true;
        materializationOptions.checkTargetReconstruction = true;
        materializationOptions.checkNoEmptyLiveNodesAfterCanonicalization = false;
        validateLocalRecompositionMaterializationInputs(
            dualTree,
            artifacts.materializationContext,
            artifacts.recompositionPlan,
            materializationOptions);
        result.materializationInputsValid = true;

        result.summary = summarizeMcpcDirectMutationPlan(
            dualTree,
            artifacts.input,
            artifacts.partitionPlan,
            artifacts.recompositionPlan,
            artifacts.materializationContext);
        result.rollbackFreeSummaryBuilt = true;
        result.passed =
            result.planOracleFree &&
            result.materializationOracleFree &&
            result.materializationInputsValid &&
            result.summary.invalidStableChildRefCount == 0 &&
            result.summary.stableChildCycleRiskCount == 0;
        if (!result.passed && result.failureSummary.empty()) {
            result.failureSummary = "P1-MCPC direct preflight did not satisfy rollback-free checks";
        }
        result.directV1Precheck = summarizeMcpcDirectV1Precheck(
            dualTree,
            artifacts.materializationContext,
            artifacts.recompositionPlan);
        result.directV1Precheck.predictedAuxiliaryNodeCount =
            result.summary.predictedAuxiliaryNodeCount;
        result.directV1PrecheckEvaluated = true;
        result.directV1PrecheckAccepted =
            evaluateMcpcDirectV1PrecheckAcceptance(
                result,
                &result.directV1PrecheckRejectionSummary);
    } catch (const std::exception &e) {
        result.failureSummary = e.what();
        if (!result.directV1PrecheckEvaluated) {
            result.directV1PrecheckAccepted =
                evaluateMcpcDirectV1PrecheckAcceptance(
                    result,
                    &result.directV1PrecheckRejectionSummary);
        }
    }

    return artifacts;
}

inline NonIncreasingMcpcDirectPreflightResult preflightMcpcDirectMergeNodeUpdate(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId nodeId) {
    return buildMcpcDirectPreflightArtifacts(
        primalTree,
        dualTree,
        adjacency,
        nodeId).preflight;
}

inline NonIncreasingMcpcD12DomainDiagnostics diagnoseMcpcD12Domain(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId nodeId) {
    return diagnoseMcpcD12Domain(
        buildMcpcDirectPreflightArtifacts(
            primalTree,
            dualTree,
            adjacency,
            nodeId).preflight);
}

inline NonIncreasingStrippedLocalUpdateResult applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight(
    DynamicComponentTree &primalTree,
    DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId nodeId,
    const NonIncreasingMcpcDirectPreflightArtifacts &artifacts,
    const NonIncreasingStrippedLocalUpdateOptions &options = {}) {
    if (!adjacency) {
        throw std::runtime_error(
            "applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight requires a valid adjacency relation");
    }
    if (!artifacts.preflight.passed) {
        throw std::runtime_error(
            "applyStrippedNonOracleLocalMergeUpdateFromMcpcPreflight requires a passing preflight");
    }

    NonIncreasingStrippedLocalUpdateResult result;
    result.inputBuilt = artifacts.preflight.inputBuilt;
    result.affectedLevelCount =
        static_cast<int>(artifacts.recompositionPlan.affectedLevels.size());
    result.localPixelCount = static_cast<int>(artifacts.input.localPixels.size());
    result.preLocalNodeCount = static_cast<int>(artifacts.input.preDualNodes.size());
    result.recompositionNodeCount =
        static_cast<int>(artifacts.recompositionPlan.nodes.size());
    result.localPlanOracleAuditChecked = artifacts.preflight.planOracleAuditChecked;
    result.localPlanOracleFree = artifacts.preflight.planOracleFree;
    result.materializationInputOracleAuditChecked =
        artifacts.preflight.materializationOracleAuditChecked;
    result.materializationInputOracleFree =
        artifacts.preflight.materializationOracleFree;

    if (options.auditPlanOracleDependencies && !result.localPlanOracleFree) {
        throw std::runtime_error(
            "prepared P1-MCPC recomposition plan exposed oracle-guided fields");
    }
    if (options.auditMaterializationOracleDependencies &&
        !result.materializationInputOracleFree) {
        throw std::runtime_error(
            "prepared P1-MCPC materialization input exposed oracle-guided fields");
    }

    std::optional<LocalMergeOracle> rebuildOracle;
    if (options.compareLocalStructureAgainstRebuild) {
        rebuildOracle = buildLocalMergeOracle(primalTree, dualTree, adjacency, nodeId);
        if (!rebuildOracle->preDualMatchesPrimal) {
            throw std::runtime_error(
                "prepared stripped structure comparison requires matching primal/dual reconstructions before merge");
        }
        if (!rebuildOracle->oracleDualMatchesUpdatedImage) {
            throw std::runtime_error(
                "prepared stripped structure comparison could not build a valid rebuilt dual oracle");
        }
    }

    auto stageStart = std::chrono::steady_clock::now();
    LocalRecompositionMaterializationOptions materializationOptions;
    materializationOptions.checkPreReconstruction =
        options.checkPreMaterializationReconstruction;
    materializationOptions.checkTargetReconstruction =
        options.checkMaterializationTargetReconstruction;
    materializationOptions.checkNoEmptyLiveNodesAfterCanonicalization =
        options.checkNoEmptyLiveNodesAfterCanonicalization;
    materializationOptions.releaseEmptySourceOwnersWithResidualChildren =
        options.releaseEmptySourceOwnersWithResidualChildren;
    const auto materialization = materializeLocalRecompositionNodes(
        dualTree,
        artifacts.materializationContext,
        artifacts.recompositionPlan,
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
    recordStrippedCreatedRecompositionNodes(result, materialization);
    recordStrippedMaterializationGuardDiagnostics(result, materialization);
    result.emptyCreatedNodesAfterCanonicalization =
        materialization.emptyCreatedNodesAfterCanonicalization;
    result.emptiedSourceOwnerCount =
        static_cast<int>(materialization.emptiedSourceOwners.size());
    result.releasedSourceOwnerCount =
        static_cast<int>(materialization.releasedSourceOwners.size());

    stageStart = std::chrono::steady_clock::now();
    const auto postLocalSnapshot = buildPostUpdateLocalDualSnapshot(
        dualTree,
        artifacts.input);
    result.postLocalSnapshotMicros =
        localUpdateElapsedMicros(stageStart, std::chrono::steady_clock::now());
    result.postLocalNodeCount = static_cast<int>(postLocalSnapshot.size());

    if (rebuildOracle) {
        const auto postLocalOracleWindow =
            buildPostUpdateLocalDualSnapshot(dualTree, *rebuildOracle);
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
        throw std::runtime_error(
            "prepared P1-MCPC local recomposition did not reconstruct the local target image");
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

} // namespace debug::non_increasing
