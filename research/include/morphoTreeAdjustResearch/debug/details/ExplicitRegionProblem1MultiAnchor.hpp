#pragma once

#include <algorithm>
#include <map>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "../Region.hpp"
#include "ExplicitRegionProblem1Contract.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem1 {

using ::debug::region::collectLocalConnectedComponents;

/**
 * @brief Diagnostic geometry for extending Problem 1 beyond a single direct
 * `nodeCa`.
 *
 * This type intentionally distinguishes three notions of locality:
 *
 * - owners of the pixels in `C`;
 * - the sparse merge workspace emitted by `RegionMergeWorkspace`;
 * - the general workspace obtained by adding the owner-to-interval ancestor
 *   paths to the sparse merge workspace.
 *
 * The formal H6' contract used by the alignment plan is attached to the third
 * notion. For non-empty, duplicate-free, connected `C` accepted by the
 * Problem 1 contract, after local materialization of partial owners, the
 * general workspace should have a single anchor. That anchor is `closureRoot`.
 *
 * Multiple owners, multiple extremal owners, disconnected owner partitions, or
 * multiple sparse workspace anchors are useful diagnostics, but they are not
 * by themselves evidence that a connected Problem 1 event needs multiple
 * general workspaces.
 */
template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1MultiAnchorAnalysis {
    bool analyzed = false;
    bool acceptedByPublishedProblem = false;
    bool materializationAttempted = false;
    bool materializationApplied = false;
    bool hasPostMaterializationContract = false;
    ExplicitRegionProblem1Contract<PixelType> contract;
    ExplicitRegionProblem1Contract<PixelType> postMaterializationContract;
    std::vector<NodeId> owners;
    std::vector<NodeId> extremalOwners;
    std::vector<NodeId> sparseWorkspaceAnchors;
    std::vector<NodeId> ownerIntervalWorkspaceAnchors;
    // H6' is checked through this vector: connected/materialized Problem 1
    // events should expose exactly one general anchor, named closureRoot.
    std::vector<NodeId> generalWorkspaceAnchors;
    NodeId closureRoot = InvalidNode;
    int nodeCaToClosureRootPathLength = -1;
    // H7 diagnostic: closureRoot may differ from nodeCa, but only through a
    // chain of live nodes at the same level a.
    bool nodeCaToClosureRootSameLevelPath = false;
    int ownerCount = 0;
    int extremalOwnerCount = 0;
    int ownerPartitionComponents = 0;
    int disconnectedOwnerPartitions = 0;
    int ownerIncidenceEdges = 0;
    int sparseWorkspaceNodeCount = 0;
    int sparseWorkspaceAnchorCount = 0;
    int ownerIntervalWorkspaceNodeCount = 0;
    int ownerIntervalWorkspaceAnchorCount = 0;
    int generalWorkspaceNodeCount = 0;
    int generalWorkspaceAnchorCount = 0;
    int ownerIntervalBucketMissingNodeCount = 0;
    int ownerIntervalBucketClosurePathMissingNodeCount = 0;
    int ownerIntervalBucketCOnlyProperPartMissingNodeCount = 0;
    int ownerIntervalBucketUnexplainedMissingNodeCount = 0;
    int frontierNodeCount = 0;
    std::string message;

    bool hasMultipleExtremalOwners() const {
        return extremalOwnerCount > 1;
    }

    bool hasMultipleOwners() const {
        return ownerCount > 1;
    }

    bool hasDisconnectedOwnerPartition() const {
        return disconnectedOwnerPartitions > 0;
    }

    bool hasMultipleSparseWorkspaceAnchors() const {
        return sparseWorkspaceAnchorCount > 1;
    }

    bool hasMultipleGeneralWorkspaceAnchors() const {
        return generalWorkspaceAnchorCount > 1;
    }

    bool hasClosureRoot() const {
        return closureRoot != InvalidNode;
    }

    bool closureRootMatchesActiveNodeCa() const {
        return hasClosureRoot() && closureRoot == activeContract().nodeCa;
    }

    bool closureRootDiffersFromActiveNodeCa() const {
        return hasClosureRoot() &&
               activeContract().nodeCa != InvalidNode &&
               closureRoot != activeContract().nodeCa;
    }

    bool nodeCaReachesClosureRoot() const {
        return nodeCaToClosureRootPathLength >= 0;
    }

    bool hasUnexplainedOwnerIntervalBucketGap() const {
        return ownerIntervalBucketUnexplainedMissingNodeCount > 0;
    }

    bool isMultiAnchorCandidate() const {
        return hasMultipleExtremalOwners() ||
               hasDisconnectedOwnerPartition() ||
               hasMultipleSparseWorkspaceAnchors() ||
               hasMultipleGeneralWorkspaceAnchors();
    }

    bool isMultiOwnerSingleWorkspaceCandidate() const {
        return hasMultipleOwners() && !hasMultipleGeneralWorkspaceAnchors();
    }

    const ExplicitRegionProblem1Contract<PixelType> &activeContract() const {
        return hasPostMaterializationContract ? postMaterializationContract : contract;
    }
};

namespace detail {

template<typename PixelType = AltitudeType>
RegionLevelChange<PixelType> makeProblem1ChangeFromContract(
    std::span<const PixelId> pixelsC,
    const ExplicitRegionProblem1Contract<PixelType> &contract) {
    RegionLevelChange<PixelType> change;
    change.pixelsC.assign(pixelsC.begin(), pixelsC.end());
    change.a = contract.a;
    change.b = contract.b;
    change.direction = contract.direction;
    change.assumeConnected = true;
    change.assumeMaterializedAsProperParts =
        !contract.requiresLocalMaterialization();
    return change;
}

inline void addUniqueNode(std::vector<NodeId> &nodes,
                          std::vector<unsigned char> &marks,
                          NodeId nodeId) {
    if (nodeId == InvalidNode ||
        nodeId < 0 ||
        nodeId >= static_cast<NodeId>(marks.size()) ||
        marks[static_cast<size_t>(nodeId)] != 0) {
        return;
    }

    marks[static_cast<size_t>(nodeId)] = 1;
    nodes.push_back(nodeId);
}

template<typename PixelType = AltitudeType>
bool isLevelInsideProblem1Interval(const DynamicComponentTree &tree,
                                   NodeId nodeId,
                                   PixelType a,
                                   PixelType b) {
    const PixelType level = static_cast<PixelType>(tree.getAltitude(nodeId));
    return tree.isMaxtree()
        ? (level >= a && level <= b)
        : (level <= a && level >= b);
}

template<typename PixelType = AltitudeType>
std::vector<NodeId> collectIntervalAncestorWorkspace(
    const DynamicComponentTree &tree,
    const std::vector<NodeId> &seeds,
    PixelType a,
    PixelType b) {
    std::vector<NodeId> nodes;
    std::vector<unsigned char> marks(
        static_cast<size_t>(tree.getNumInternalNodeSlots()),
        0);

    for (NodeId seed : seeds) {
        NodeId current = seed;
        while (current != InvalidNode &&
               tree.isAlive(current) &&
               isLevelInsideProblem1Interval<PixelType>(tree, current, a, b)) {
            addUniqueNode(nodes, marks, current);

            const NodeId parent = tree.getNodeParent(current);
            if (parent == current) {
                break;
            }
            current = parent;
        }
    }

    return nodes;
}

inline std::vector<NodeId> collectWorkspaceAnchors(const DynamicComponentTree &tree,
                                                   const std::vector<NodeId> &nodes) {
    std::vector<NodeId> anchors;
    std::vector<unsigned char> marks(
        static_cast<size_t>(tree.getNumInternalNodeSlots()),
        0);

    for (NodeId nodeId : nodes) {
        addUniqueNode(anchors, marks, nodeId);
    }

    anchors.clear();
    for (NodeId nodeId : nodes) {
        const NodeId parent = tree.getNodeParent(nodeId);
        if (parent == nodeId ||
            parent == InvalidNode ||
            parent < 0 ||
            parent >= static_cast<NodeId>(marks.size()) ||
            marks[static_cast<size_t>(parent)] == 0) {
            anchors.push_back(nodeId);
        }
    }

    return anchors;
}

inline int pathLengthToAncestor(const DynamicComponentTree &tree,
                                NodeId nodeId,
                                NodeId ancestorId) {
    if (nodeId == InvalidNode || ancestorId == InvalidNode) {
        return -1;
    }

    int length = 0;
    NodeId current = nodeId;
    while (current != InvalidNode && tree.isAlive(current)) {
        if (current == ancestorId) {
            return length;
        }
        const NodeId parent = tree.getNodeParent(current);
        if (parent == current) {
            break;
        }
        current = parent;
        ++length;
    }

    return -1;
}

inline bool pathToAncestorHasSingleLevel(const DynamicComponentTree &tree,
                                         NodeId nodeId,
                                         NodeId ancestorId) {
    if (nodeId == InvalidNode || ancestorId == InvalidNode) {
        return false;
    }

    const int level = tree.getAltitude(nodeId);
    NodeId current = nodeId;
    while (current != InvalidNode && tree.isAlive(current)) {
        if (tree.getAltitude(current) != level) {
            return false;
        }
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

inline bool isNodeOnPathToAncestor(const DynamicComponentTree &tree,
                                   NodeId startNodeId,
                                   NodeId ancestorId,
                                   NodeId candidateId) {
    if (startNodeId == InvalidNode ||
        ancestorId == InvalidNode ||
        candidateId == InvalidNode) {
        return false;
    }

    NodeId current = startNodeId;
    while (current != InvalidNode && tree.isAlive(current)) {
        if (current == candidateId) {
            return true;
        }
        if (current == ancestorId) {
            return candidateId == ancestorId;
        }
        const NodeId parent = tree.getNodeParent(current);
        if (parent == current) {
            break;
        }
        current = parent;
    }

    return false;
}

inline bool directProperPartsAreInsideRegion(const DynamicComponentTree &tree,
                                             NodeId nodeId,
                                             const std::vector<unsigned char> &inRegion) {
    if (nodeId == InvalidNode || !tree.isAlive(nodeId)) {
        return false;
    }

    bool sawProperPart = false;
    for (PixelId pixelId : tree.getProperParts(nodeId)) {
        sawProperPart = true;
        if (pixelId < 0 ||
            pixelId >= static_cast<PixelId>(inRegion.size()) ||
            inRegion[static_cast<size_t>(pixelId)] == 0) {
            return false;
        }
    }

    return sawProperPart;
}

} // namespace detail

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1MultiAnchorAnalysis<PixelType>
analyzeExplicitRegionProblem1MultiAnchor(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false,
    bool allowLocalMaterialization = true) {
    ExplicitRegionProblem1MultiAnchorAnalysis<PixelType> out;

    if (!input) {
        out.message = "multi-anchor analysis requires a valid input image";
        return out;
    }
    if (!adjacency) {
        out.message = "multi-anchor analysis requires a valid adjacency relation";
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
        out.message = out.contract.rejectionReason;
        return out;
    }

    out.acceptedByPublishedProblem = true;
    const ExplicitRegionProblem1Contract<PixelType> *activeContract = &out.contract;

    if (allowLocalMaterialization && out.contract.requiresLocalMaterialization()) {
        out.materializationAttempted = true;
        const auto materializationPlan = analyzeRegionMaterialization(
            dualTree,
            *adjacency,
            pixelsC);
        if (materializationPlan.requiresLocalMaterialization()) {
            materializePartialOwners(dualTree, materializationPlan);
            out.materializationApplied = true;
            out.postMaterializationContract = checkExplicitRegionProblem1Contract<PixelType>(
                dualTree,
                *adjacency,
                pixelsC,
                b,
                true);
            out.hasPostMaterializationContract = true;
            activeContract = &out.postMaterializationContract;
        }
    }

    if (!activeContract->acceptsPublishedProblem()) {
        out.message = activeContract->rejectionReason;
        return out;
    }

    const int domainSize = dualTree.getNumTotalProperParts();
    std::vector<unsigned char> inRegion(static_cast<size_t>(domainSize), 0);
    std::map<NodeId, std::vector<PixelId>> pixelsByOwner;

    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            out.message = "multi-anchor analysis received an out-of-bounds pixel";
            return out;
        }
        inRegion[static_cast<size_t>(pixelId)] = 1;

        const NodeId owner = dualTree.getSmallestComponent(pixelId);
        if (owner == InvalidNode || !dualTree.isAlive(owner)) {
            out.message = "multi-anchor analysis cannot locate a live owner";
            return out;
        }
        pixelsByOwner[owner].push_back(pixelId);
    }

    out.owners.reserve(pixelsByOwner.size());
    for (const auto &[owner, ownerPixels] : pixelsByOwner) {
        out.owners.push_back(owner);
        if (dualTree.getAltitude(owner) == static_cast<int>(activeContract->a)) {
            out.extremalOwners.push_back(owner);
        }

        const auto ownerComponents = collectLocalConnectedComponents(
            *adjacency,
            std::span<const PixelId>(ownerPixels),
            domainSize);
        out.ownerPartitionComponents += static_cast<int>(ownerComponents.size());
        if (ownerComponents.size() > 1) {
            ++out.disconnectedOwnerPartitions;
        }
    }
    out.ownerCount = static_cast<int>(out.owners.size());
    out.extremalOwnerCount = static_cast<int>(out.extremalOwners.size());

    std::set<std::pair<NodeId, NodeId>> incidenceEdges;
    for (PixelId pixelId : pixelsC) {
        const NodeId owner = dualTree.getSmallestComponent(pixelId);
        for (PixelId neighbor : adjacency->getNeighborPixels(pixelId)) {
            if (neighbor < 0 ||
                neighbor >= domainSize ||
                inRegion[static_cast<size_t>(neighbor)] == 0) {
                continue;
            }
            const NodeId neighborOwner = dualTree.getSmallestComponent(neighbor);
            if (owner == neighborOwner ||
                owner == InvalidNode ||
                neighborOwner == InvalidNode) {
                continue;
            }

            incidenceEdges.insert(std::minmax(owner, neighborOwner));
        }
    }
    out.ownerIncidenceEdges = static_cast<int>(incidenceEdges.size());

    const RegionLevelChange<PixelType> change =
        detail::makeProblem1ChangeFromContract<PixelType>(pixelsC, *activeContract);
    const RegionMergeWorkspace<PixelType> workspace =
        buildRegionMergeWorkspace(dualTree, *adjacency, change);

    const std::vector<NodeId> ownerIntervalWorkspaceNodes =
        detail::collectIntervalAncestorWorkspace<PixelType>(
            dualTree,
            out.owners,
            activeContract->a,
            activeContract->b);
    out.ownerIntervalWorkspaceNodeCount =
        static_cast<int>(ownerIntervalWorkspaceNodes.size());
    out.ownerIntervalWorkspaceAnchors =
        detail::collectWorkspaceAnchors(dualTree, ownerIntervalWorkspaceNodes);
    out.ownerIntervalWorkspaceAnchorCount =
        static_cast<int>(out.ownerIntervalWorkspaceAnchors.size());

    std::vector<NodeId> workspaceNodes;
    std::vector<unsigned char> workspaceMarks(
        static_cast<size_t>(dualTree.getNumInternalNodeSlots()),
        0);
    std::vector<unsigned char> mergeBucketMarks(
        static_cast<size_t>(dualTree.getNumInternalNodeSlots()),
        0);
    detail::addUniqueNode(
        workspaceNodes,
        workspaceMarks,
        activeContract->nodeCa);
    for (NodeId nodeId : workspace.mergeNodesIntervalMinus) {
        detail::addUniqueNode(workspaceNodes, workspaceMarks, nodeId);
        if (nodeId >= 0 && nodeId < static_cast<NodeId>(mergeBucketMarks.size())) {
            mergeBucketMarks[static_cast<size_t>(nodeId)] = 1;
        }
    }

    out.sparseWorkspaceNodeCount = static_cast<int>(workspaceNodes.size());
    out.frontierNodeCount = static_cast<int>(workspace.frontierNodesAboveB.size());
    for (NodeId nodeId : workspaceNodes) {
        NodeId parent = dualTree.getNodeParent(nodeId);
        if (parent == nodeId ||
            parent == InvalidNode ||
            parent < 0 ||
            parent >= static_cast<NodeId>(workspaceMarks.size()) ||
            workspaceMarks[static_cast<size_t>(parent)] == 0) {
            out.sparseWorkspaceAnchors.push_back(nodeId);
        }
    }
    out.sparseWorkspaceAnchorCount =
        static_cast<int>(out.sparseWorkspaceAnchors.size());

    // General workspace used by H6': union the owner ancestor paths in the
    // oriented [a,b] interval with all concrete merge buckets. For connected C,
    // all those paths represent the same threshold component at level a, so the
    // resulting tree-induced subgraph should have exactly one anchor.
    std::vector<NodeId> generalWorkspaceNodes = ownerIntervalWorkspaceNodes;
    std::vector<unsigned char> generalWorkspaceMarks(
        static_cast<size_t>(dualTree.getNumInternalNodeSlots()),
        0);
    for (NodeId nodeId : generalWorkspaceNodes) {
        if (nodeId >= 0 && nodeId < static_cast<NodeId>(generalWorkspaceMarks.size())) {
            generalWorkspaceMarks[static_cast<size_t>(nodeId)] = 1;
        }
    }
    for (NodeId nodeId : workspace.mergeNodesIntervalMinus) {
        detail::addUniqueNode(
            generalWorkspaceNodes,
            generalWorkspaceMarks,
            nodeId);
    }
    out.generalWorkspaceNodeCount =
        static_cast<int>(generalWorkspaceNodes.size());
    out.generalWorkspaceAnchors =
        detail::collectWorkspaceAnchors(dualTree, generalWorkspaceNodes);
    out.generalWorkspaceAnchorCount =
        static_cast<int>(out.generalWorkspaceAnchors.size());
    if (out.generalWorkspaceAnchorCount == 1) {
        out.closureRoot = out.generalWorkspaceAnchors.front();
        out.nodeCaToClosureRootPathLength =
            detail::pathLengthToAncestor(
                dualTree,
                activeContract->nodeCa,
                out.closureRoot);
        out.nodeCaToClosureRootSameLevelPath =
            detail::pathToAncestorHasSingleLevel(
                dualTree,
                activeContract->nodeCa,
                out.closureRoot);
    }

    // H8 diagnostic: owner-interval nodes missing from concrete merge buckets
    // are safe only when they are accounted for by the same-level closure path
    // or by direct proper parts entirely inside C. Any remaining gap means the
    // current workspace has no concrete merge representative for an affected
    // owner-interval node and must not be released as a local Problem 1 update.
    for (NodeId nodeId : ownerIntervalWorkspaceNodes) {
        const bool inMergeBucket =
            nodeId >= 0 &&
            nodeId < static_cast<NodeId>(mergeBucketMarks.size()) &&
            mergeBucketMarks[static_cast<size_t>(nodeId)] != 0;
        if (inMergeBucket) {
            continue;
        }

        ++out.ownerIntervalBucketMissingNodeCount;
        if (detail::isNodeOnPathToAncestor(
                dualTree,
                activeContract->nodeCa,
                out.closureRoot,
                nodeId)) {
            ++out.ownerIntervalBucketClosurePathMissingNodeCount;
        } else if (detail::directProperPartsAreInsideRegion(
                       dualTree,
                       nodeId,
                       inRegion)) {
            ++out.ownerIntervalBucketCOnlyProperPartMissingNodeCount;
        } else {
            ++out.ownerIntervalBucketUnexplainedMissingNodeCount;
        }
    }

    out.analyzed = true;
    out.message = "multi-anchor analysis completed";
    return out;
}

} // namespace debug::problem1
