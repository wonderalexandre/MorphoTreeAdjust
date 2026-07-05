#pragma once

#include <algorithm>
#include <functional>
#include <map>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

#include "Algorithm2FlatRaiseEligibility.hpp"
#include "RegionLevelChange.hpp"
#include "RegionMaterialization.hpp"
#include "RegionMergeWorkspace.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem2 {

using ::debug::region::LocalMaterializationResult;
using ::debug::region::RegionLevelChange;
using ::debug::region::RegionMaterializationResult;
using ::debug::region::RegionMergeWorkspace;
using ::debug::region::buildRegionMergeWorkspace;
using ::debug::region::directionFromLevels;
using ::debug::region::isConnectedRegion;
using ::debug::region::materializePartialOwners;
using ::debug::region::materializePartialOwnersInsideOnly;

struct Algorithm2CreatedNodeAtLevel {
    AltitudeType level = 0;
    NodeId nodeId = InvalidNode;
    bool isQb = false;
};

struct Algorithm2XiBulletLevel {
    AltitudeType level = 0;
    std::vector<NodeId> nodes;
    NodeId createdNodeQh = InvalidNode;
};

struct Algorithm2FlatRaiseApplyResult {
    ImageUInt8Ptr reconstruction;
    RegionMergeWorkspace<AltitudeType> workspace;
    NodeId sourceNodeQa = InvalidNode;
    NodeId nodeQb = InvalidNode;
    NodeId recompositionRoot = InvalidNode;
    std::vector<Algorithm2XiBulletLevel> xiBulletByLevel;
    std::vector<NodeId> xiBulletNodes;
    std::vector<NodeId> xiCircleNodes;
    std::vector<AltitudeType> affectedLevels;
    std::vector<NodeId> createdNodes;
    std::vector<Algorithm2CreatedNodeAtLevel> createdNodesByLevel;
    std::vector<NodeId> absorbedNodes;
    std::vector<NodeId> releasedNodes;
    std::vector<NodeId> attributeDirtySeeds;
    bool materializedBeforeUpdate = false;
    bool rhoQaEqualsPBeforeUpdate = false;
    bool releasedSourceNodeQa = false;
    bool preservedSourceNodeQa = false;

    int numCreatedNodes() const {
        return static_cast<int>(createdNodes.size());
    }

    int numAbsorbedNodes() const {
        return static_cast<int>(absorbedNodes.size());
    }

    NodeId createdNodeAtLevel(AltitudeType level) const {
        for (const auto &created : createdNodesByLevel) {
            if (created.level == level) {
                return created.nodeId;
            }
        }
        return InvalidNode;
    }
};

class Algorithm2FlatRegionMergeUpdater {
public:
    static Algorithm2FlatRaiseApplyResult updateMaxTreeFlatRaiseSingleOwnerMaterialized(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatRaiseInput(maxTree, adjacency, change);

        RegionLevelChange<AltitudeType> materializedChange = change;
        materializedChange.assumeMaterializedAsProperParts = true;

        const auto eligibility =
            requireAlgorithm2FlatRaiseEligibility(maxTree, adjacency, materializedChange);
        const NodeId sourceNodeQa = eligibility.sourceNodeQa;
        materializedChange.a = eligibility.a;
        materializedChange.direction =
            directionFromLevels(materializedChange.a, materializedChange.b);

        auto workspace = buildRegionMergeWorkspace(maxTree, adjacency, materializedChange);
        return applyMaterializedAlgorithm2(maxTree,
                                           materializedChange,
                                           std::move(workspace),
                                           sourceNodeQa,
                                           false);
    }

    static Algorithm2FlatRaiseApplyResult updateMaxTreeFlatRaiseMaterialized(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        return updateMaxTreeFlatRaiseSingleOwnerMaterialized(maxTree, adjacency, change);
    }

    static Algorithm2FlatRaiseApplyResult updateMaxTreeFlatRaiseSingleOwner(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatRaiseInput(maxTree, adjacency, change);

        RegionLevelChange<AltitudeType> materializedChange = change;
        materializedChange.assumeMaterializedAsProperParts = true;

        const auto eligibility =
            requireAlgorithm2FlatRaiseEligibility(maxTree, adjacency, materializedChange);
        const NodeId sourceNodeQa = eligibility.sourceNodeQa;
        materializedChange.a = eligibility.a;
        materializedChange.direction =
            directionFromLevels(materializedChange.a, materializedChange.b);

        auto workspace = buildRegionMergeWorkspace(maxTree, adjacency, materializedChange);
        return applyMaterializedAlgorithm2(maxTree,
                                           materializedChange,
                                           std::move(workspace),
                                           sourceNodeQa,
                                           false);
    }

    static Algorithm2FlatRaiseApplyResult updateMaxTreeFlatRaise(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        return updateMaxTreeFlatRaiseSingleOwner(maxTree, adjacency, change);
    }

    static Algorithm2FlatRaiseApplyResult updateMaxTreeFlatRaiseAfterMaterializingPartialOwners(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatRaiseInput(maxTree, adjacency, change);

        RegionLevelChange<AltitudeType> materializedChange = change;
        bool materialized = false;

        NodeId sourceNodeQa = InvalidNode;
        try {
            const auto eligibility =
                requireAlgorithm2FlatRaiseEligibility(maxTree, adjacency, materializedChange);
            sourceNodeQa = eligibility.sourceNodeQa;
            materializedChange.a = eligibility.a;
        } catch (const std::runtime_error &) {
            auto materialization = analyzeRegionMaterialization(maxTree, adjacency, change);
            if (materialization.requiresLocalMaterialization()) {
                (void) materializePartialOwners(maxTree, materialization);
                materialized = true;
            }
            const auto eligibility =
                requireAlgorithm2FlatRaiseEligibility(maxTree, adjacency, materializedChange);
            sourceNodeQa = eligibility.sourceNodeQa;
            materializedChange.a = eligibility.a;
        }

        materializedChange.assumeMaterializedAsProperParts = true;
        materializedChange.direction =
            directionFromLevels(materializedChange.a, materializedChange.b);

        auto workspace = buildRegionMergeWorkspace(maxTree, adjacency, materializedChange);
        return applyMaterializedAlgorithm2(maxTree,
                                           materializedChange,
                                           std::move(workspace),
                                           sourceNodeQa,
                                           materialized);
    }

private:
    static std::vector<PixelId> collectProperParts(const DynamicComponentTree &tree,
                                                   NodeId nodeId) {
        std::vector<PixelId> pixels;
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            pixels.push_back(pixelId);
        }
        return pixels;
    }

    static std::vector<NodeId> collectChildren(const DynamicComponentTree &tree,
                                               NodeId nodeId) {
        std::vector<NodeId> children;
        for (NodeId childId : tree.getChildren(nodeId)) {
            children.push_back(childId);
        }
        return children;
    }

    static void recordAttributeDirtySeed(Algorithm2FlatRaiseApplyResult &result,
                                         NodeId nodeId) {
        if (nodeId != InvalidNode) {
            result.attributeDirtySeeds.push_back(nodeId);
        }
    }

    static void validateFlatRaiseInput(const DynamicComponentTree &maxTree,
                                       AdjacencyRelation &adjacency,
                                       const RegionLevelChange<AltitudeType> &change) {
        if (!maxTree.isMaxtree()) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater requires a max-tree");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater requires a non-empty flat region P");
        }
        if (change.b <= change.a) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater currently supports only raise updates b > a");
        }
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater requires connected region P");
        }

        const int domainSize = maxTree.getNumTotalProperParts();
        std::vector<unsigned char> seen(static_cast<size_t>(domainSize), 0);
        auto image = maxTree.reconstructionImage();
        for (PixelId pixelId : change.pixelsC) {
            if (pixelId < 0 || pixelId >= domainSize) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater received an out-of-bounds pixel");
            }
            const size_t index = static_cast<size_t>(pixelId);
            if (seen[index] != 0) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater received duplicated pixels in P");
            }
            seen[index] = 1;

            if ((*image)[pixelId] != change.a) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater requires every pixel in P to have level a");
            }
        }
    }

    static NodeId validateSingleSourceNodeContainingP(
        const DynamicComponentTree &maxTree,
        const RegionLevelChange<AltitudeType> &change) {
        NodeId sourceNode = InvalidNode;

        for (PixelId pixelId : change.pixelsC) {
            const NodeId ownerNode = maxTree.getSmallestComponent(pixelId);
            if (ownerNode == InvalidNode || !maxTree.isAlive(ownerNode)) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater cannot locate a live owner for P");
            }
            if (sourceNode == InvalidNode) {
                sourceNode = ownerNode;
            } else if (sourceNode != ownerNode) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater currently requires P to be owned by a single source node Q_a");
            }
        }

        if (sourceNode == InvalidNode) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater cannot locate Q_a");
        }
        if (maxTree.getAltitude(sourceNode) != change.a) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater requires Q_a altitude to be a");
        }
        for (PixelId pixelId : change.pixelsC) {
            bool foundAsDirectProperPart = false;
            for (PixelId sourcePixel : maxTree.getProperParts(sourceNode)) {
                if (sourcePixel == pixelId) {
                    foundAsDirectProperPart = true;
                    break;
                }
            }
            if (!foundAsDirectProperPart) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater requires every pixel in P to be a direct proper part of Q_a");
            }
        }

        return sourceNode;
    }

    static NodeId detachIfAttached(DynamicComponentTree &tree, NodeId nodeId) {
        if (nodeId == InvalidNode || !tree.isAlive(nodeId)) {
            return InvalidNode;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId != InvalidNode && parentId != nodeId) {
            tree.detachNode(nodeId);
            return parentId;
        }
        return InvalidNode;
    }

    static NodeId attachChildIfNeeded(DynamicComponentTree &tree,
                                      NodeId parentId,
                                      NodeId childId) {
        if (parentId == InvalidNode ||
            childId == InvalidNode ||
            parentId == childId ||
            !tree.isAlive(parentId) ||
            !tree.isAlive(childId)) {
            return InvalidNode;
        }
        if (tree.getNodeParent(childId) == parentId) {
            return InvalidNode;
        }
        const NodeId oldParentId = detachIfAttached(tree, childId);
        tree.attachNode(parentId, childId);
        return oldParentId;
    }

    static void releaseEmptyDetachedOrChildNode(DynamicComponentTree &tree,
                                                NodeId nodeId,
                                                Algorithm2FlatRaiseApplyResult &result) {
        if (nodeId == InvalidNode || !tree.isAlive(nodeId)) {
            return;
        }
        if (tree.getNumProperParts(nodeId) != 0 || tree.getNumChildren(nodeId) != 0) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater tried to release a non-empty node");
        }

        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId != InvalidNode && parentId != nodeId) {
            tree.removeChild(parentId, nodeId, true);
            recordAttributeDirtySeed(result, parentId);
        } else {
            tree.releaseNode(nodeId);
        }
        recordAttributeDirtySeed(result, nodeId);
        result.releasedNodes.push_back(nodeId);
    }

    static void moveAllDirectProperParts(DynamicComponentTree &tree,
                                         NodeId targetNode,
                                         NodeId sourceNode) {
        auto properParts = collectProperParts(tree, sourceNode);
        for (PixelId pixelId : properParts) {
            tree.moveProperPart(targetNode, sourceNode, pixelId);
        }
    }

    static void absorbOldMergeNode(DynamicComponentTree &tree,
                                   NodeId oldNode,
                                   NodeId newNode,
                                   const std::unordered_set<NodeId> &mergeNodeSet,
                                   Algorithm2FlatRaiseApplyResult &result) {
        if (oldNode == InvalidNode ||
            oldNode == newNode ||
            !tree.isAlive(oldNode)) {
            return;
        }

        for (NodeId childId : collectChildren(tree, oldNode)) {
            if (!tree.isAlive(childId) ||
                mergeNodeSet.find(childId) != mergeNodeSet.end()) {
                continue;
            }
            const NodeId oldParentId = attachChildIfNeeded(tree, newNode, childId);
            recordAttributeDirtySeed(result, oldParentId);
            recordAttributeDirtySeed(result, newNode);
            recordAttributeDirtySeed(result, childId);
        }

        moveAllDirectProperParts(tree, newNode, oldNode);
        recordAttributeDirtySeed(result, newNode);
        recordAttributeDirtySeed(result, oldNode);

        if (tree.getNumChildren(oldNode) != 0) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater found an unprocessed child in an absorbed node");
        }

        result.absorbedNodes.push_back(oldNode);
        releaseEmptyDetachedOrChildNode(tree, oldNode, result);
    }

    static void contractSameAltitudeChildren(DynamicComponentTree &tree,
                                             Algorithm2FlatRaiseApplyResult &result) {
        bool changed = true;
        while (changed) {
            changed = false;
            std::vector<NodeId> nodes;
            for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
                nodes.push_back(nodeId);
            }

            for (NodeId parentId : nodes) {
                if (!tree.isAlive(parentId)) {
                    continue;
                }
                for (NodeId childId : collectChildren(tree, parentId)) {
                    if (!tree.isAlive(childId)) {
                        continue;
                    }
                    if (tree.getAltitude(parentId) == tree.getAltitude(childId)) {
                        tree.mergeNodeIntoParent(childId);
                        recordAttributeDirtySeed(result, parentId);
                        recordAttributeDirtySeed(result, childId);
                        changed = true;
                        break;
                    }
                }
                if (changed) {
                    break;
                }
            }
        }
    }

    static void validateNoEmptyLiveNodes(const DynamicComponentTree &tree) {
        for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
            if (tree.getNumProperParts(nodeId) == 0) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater left a live node without direct proper parts");
            }
        }
    }

    static void validateStrictMaxTreeOrder(const DynamicComponentTree &tree) {
        for (NodeId parentId : tree.getNodeSubtree(tree.getRoot())) {
            for (NodeId childId : tree.getChildren(parentId)) {
                if (tree.getAltitude(parentId) >= tree.getAltitude(childId)) {
                    throw std::runtime_error("Algorithm2FlatRegionMergeUpdater left a non-strict max-tree altitude edge");
                }
            }
        }
    }

    static Algorithm2FlatRaiseApplyResult applyMaterializedAlgorithm2(
        DynamicComponentTree &maxTree,
        const RegionLevelChange<AltitudeType> &change,
        RegionMergeWorkspace<AltitudeType> workspace,
        NodeId sourceNodeQa,
        bool materializedBeforeUpdate) {
        Algorithm2FlatRaiseApplyResult result;
        result.workspace = std::move(workspace);
        result.sourceNodeQa = sourceNodeQa;
        result.materializedBeforeUpdate = materializedBeforeUpdate;
        result.rhoQaEqualsPBeforeUpdate =
            maxTree.getNumProperParts(sourceNodeQa) == static_cast<int>(change.pixelsC.size());
        result.xiCircleNodes = result.workspace.frontierNodesAboveB;

        std::unordered_set<NodeId> mergeNodeSet;
        std::map<int, std::vector<NodeId>> bucketsByLevel;
        for (const auto &bucket : result.workspace.mergeNodesByLevel) {
            const int level = static_cast<int>(bucket.level);
            if (level <= static_cast<int>(change.a) ||
                level > static_cast<int>(change.b)) {
                continue;
            }
            for (NodeId nodeId : bucket.mergeNodesVMinus) {
                if (nodeId == InvalidNode ||
                    nodeId == sourceNodeQa ||
                    !maxTree.isAlive(nodeId)) {
                    continue;
                }
                bucketsByLevel[level].push_back(nodeId);
                mergeNodeSet.insert(nodeId);
            }
        }

        for (const auto &[level, nodes] : bucketsByLevel) {
            Algorithm2XiBulletLevel xiLevel;
            xiLevel.level = static_cast<AltitudeType>(level);
            xiLevel.nodes = nodes;
            result.xiBulletByLevel.push_back(xiLevel);
            result.xiBulletNodes.insert(result.xiBulletNodes.end(),
                                        nodes.begin(),
                                        nodes.end());
        }
        std::sort(result.xiBulletByLevel.begin(),
                  result.xiBulletByLevel.end(),
                  [](const auto &lhs, const auto &rhs) {
                      return lhs.level > rhs.level;
                  });

        std::vector<int> levels;
        levels.push_back(static_cast<int>(change.b));
        for (const auto &[level, nodes] : bucketsByLevel) {
            static_cast<void>(nodes);
            levels.push_back(level);
        }
        std::sort(levels.begin(), levels.end(), std::greater<int>());
        levels.erase(std::unique(levels.begin(), levels.end()), levels.end());
        for (int level : levels) {
            result.affectedLevels.push_back(static_cast<AltitudeType>(level));
        }

        NodeId previousHigherNode = InvalidNode;
        for (int level : levels) {
            const NodeId newNode = maxTree.createDetachedNode(level);
            if (newNode == InvalidNode) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater could not create a recomposition node");
            }
            result.createdNodes.push_back(newNode);
            recordAttributeDirtySeed(result, newNode);
            Algorithm2CreatedNodeAtLevel createdAtLevel;
            createdAtLevel.level = static_cast<AltitudeType>(level);
            createdAtLevel.nodeId = newNode;
            createdAtLevel.isQb = level == static_cast<int>(change.b);
            result.createdNodesByLevel.push_back(createdAtLevel);
            if (level == static_cast<int>(change.b)) {
                result.nodeQb = newNode;
            }
            for (auto &xiLevel : result.xiBulletByLevel) {
                if (xiLevel.level == static_cast<AltitudeType>(level)) {
                    xiLevel.createdNodeQh = newNode;
                    break;
                }
            }

            auto bucketIt = bucketsByLevel.find(level);
            if (bucketIt != bucketsByLevel.end()) {
                for (NodeId oldNode : bucketIt->second) {
                    absorbOldMergeNode(maxTree,
                                       oldNode,
                                       newNode,
                                       mergeNodeSet,
                                       result);
                }
            }

            if (level == static_cast<int>(change.b)) {
                for (PixelId pixelId : change.pixelsC) {
                    const NodeId owner = maxTree.getSmallestComponent(pixelId);
                    if (owner != sourceNodeQa) {
                        throw std::runtime_error("Algorithm2FlatRegionMergeUpdater can only move direct proper parts of Q_a into Q_b");
                    }
                    maxTree.moveProperPart(newNode, sourceNodeQa, pixelId);
                }
                recordAttributeDirtySeed(result, newNode);
                recordAttributeDirtySeed(result, sourceNodeQa);

                for (NodeId frontierNode : result.workspace.frontierNodesAboveB) {
                    if (frontierNode == InvalidNode ||
                        !maxTree.isAlive(frontierNode) ||
                        frontierNode == newNode) {
                        continue;
                    }
                    const NodeId oldParentId = attachChildIfNeeded(maxTree, newNode, frontierNode);
                    recordAttributeDirtySeed(result, oldParentId);
                    recordAttributeDirtySeed(result, newNode);
                    recordAttributeDirtySeed(result, frontierNode);
                }
            }

            if (previousHigherNode != InvalidNode &&
                maxTree.isAlive(previousHigherNode)) {
                const NodeId oldParentId = attachChildIfNeeded(maxTree, newNode, previousHigherNode);
                recordAttributeDirtySeed(result, oldParentId);
                recordAttributeDirtySeed(result, newNode);
                recordAttributeDirtySeed(result, previousHigherNode);
            }

            previousHigherNode = newNode;
        }

        if (previousHigherNode == InvalidNode) {
            throw std::runtime_error("Algorithm2FlatRegionMergeUpdater did not create the Q_b node");
        }

        result.recompositionRoot = previousHigherNode;

        const bool removeSourceNodeQa = maxTree.getNumProperParts(sourceNodeQa) == 0;
        if (removeSourceNodeQa) {
            for (NodeId childId : collectChildren(maxTree, sourceNodeQa)) {
                if (!maxTree.isAlive(childId) ||
                    childId == result.recompositionRoot) {
                    continue;
                }
                if (maxTree.getAltitude(childId) <= maxTree.getAltitude(result.recompositionRoot)) {
                    throw std::runtime_error("Algorithm2FlatRegionMergeUpdater cannot attach a remaining child below the recomposition root");
                }
                const NodeId oldParentId = attachChildIfNeeded(maxTree, result.recompositionRoot, childId);
                recordAttributeDirtySeed(result, oldParentId);
                recordAttributeDirtySeed(result, result.recompositionRoot);
                recordAttributeDirtySeed(result, childId);
            }

            if (maxTree.getNumChildren(sourceNodeQa) != 0) {
                throw std::runtime_error("Algorithm2FlatRegionMergeUpdater expected Q_a to be childless before release");
            }

            const bool sourceWasRoot = maxTree.isRoot(sourceNodeQa);
            const NodeId sourceParent = maxTree.getNodeParent(sourceNodeQa);
            if (sourceWasRoot) {
                maxTree.setRoot(result.recompositionRoot);
                recordAttributeDirtySeed(result, result.recompositionRoot);
                recordAttributeDirtySeed(result, sourceNodeQa);
                releaseEmptyDetachedOrChildNode(maxTree,
                                                sourceNodeQa,
                                                result);
            } else {
                if (sourceParent == InvalidNode || sourceParent == sourceNodeQa) {
                    throw std::runtime_error("Algorithm2FlatRegionMergeUpdater cannot locate Q_a parent");
                }
                const NodeId oldParentId = attachChildIfNeeded(maxTree, sourceParent, result.recompositionRoot);
                recordAttributeDirtySeed(result, oldParentId);
                recordAttributeDirtySeed(result, sourceParent);
                recordAttributeDirtySeed(result, result.recompositionRoot);
                recordAttributeDirtySeed(result, sourceNodeQa);
                releaseEmptyDetachedOrChildNode(maxTree,
                                                sourceNodeQa,
                                                result);
            }
            result.releasedSourceNodeQa = true;
        } else {
            const NodeId oldParentId = attachChildIfNeeded(maxTree, sourceNodeQa, result.recompositionRoot);
            recordAttributeDirtySeed(result, oldParentId);
            recordAttributeDirtySeed(result, sourceNodeQa);
            recordAttributeDirtySeed(result, result.recompositionRoot);
            result.preservedSourceNodeQa = true;
        }

        contractSameAltitudeChildren(maxTree, result);
        validateNoEmptyLiveNodes(maxTree);
        validateStrictMaxTreeOrder(maxTree);

        result.reconstruction = maxTree.reconstructionImage();
        return result;
    }
};

} // namespace debug::problem2
