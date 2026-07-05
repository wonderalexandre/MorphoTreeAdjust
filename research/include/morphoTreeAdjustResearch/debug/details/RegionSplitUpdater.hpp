#pragma once

#include <algorithm>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../Region.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::region {

enum class RegionSplitLocalOutcome {
    Removed,
    Contracted,
    Split,
};

enum class RegionSplitBatchInteraction {
    IndependentSingleOwners,
    SharedExistingLevelBTarget,
    RepeatedSourceNode,
    RequiresInsertedLevelBAncestor,
    MixedInsertedAndNonInsertedTargets,
    RequiresMultipleInsertedLevelBAncestors,
    InteractingInsertedLevelBAncestors,
    PreservedAncestorSplits,
};

struct RegionSplitAffectedComponent {
    int level = 0;
    NodeId nodeId = InvalidNode;
    std::vector<PixelId> originalSupport;
    std::vector<PixelId> removedPixels;
    std::vector<std::vector<PixelId>> remainingComponents;
    RegionSplitLocalOutcome outcome = RegionSplitLocalOutcome::Removed;
};

struct RegionSplitPlan {
    RegionLevelChange<AltitudeType> change;
    int firstAffectedLevel = 0;
    int lastAffectedLevel = 0;
    std::vector<RegionSplitAffectedComponent> affectedComponents;

    bool empty() const {
        return affectedComponents.empty();
    }
};

struct RegionSplitApplyResult {
    ImageUInt8Ptr reconstruction;
    NodeId sourceNode = InvalidNode;
    NodeId targetNode = InvalidNode;
    std::vector<NodeId> sourceNodes;
    bool createdTargetNode = false;
    bool releasedSourceNode = false;
};

struct RegionSplitBatchApplyResult {
    ImageUInt8Ptr reconstruction;
    std::vector<RegionSplitApplyResult> componentResults;

    int numComponents() const {
        return static_cast<int>(componentResults.size());
    }

    int numCreatedTargetNodes() const {
        int count = 0;
        for (const auto &result : componentResults) {
            if (result.createdTargetNode) {
                ++count;
            }
        }
        return count;
    }

    int numReleasedSourceNodes() const {
        int count = 0;
        for (const auto &result : componentResults) {
            if (result.releasedSourceNode) {
                ++count;
            }
        }
        return count;
    }
};

struct RegionSplitMaterializationResult {
    ImageUInt8Ptr reconstruction;
    NodeId originalOwner = InvalidNode;
    std::vector<NodeId> originalOwners;
    std::vector<NodeId> createdOutsideNodes;
    bool changed = false;

    int numCreatedNodes() const {
        return static_cast<int>(createdOutsideNodes.size());
    }
};

struct RegionSplitBatchComponentPlan {
    std::vector<PixelId> pixelsC;
    NodeId sourceNode = InvalidNode;
    NodeId parentNode = InvalidNode;
    bool sourceIsRoot = false;
    bool reusesExistingLevelBTarget = false;
    bool requiresInsertedLevelBAncestor = false;
};

struct RegionSplitBatchPlan {
    RegionLevelChange<AltitudeType> change;
    std::vector<RegionSplitBatchComponentPlan> components;
    RegionSplitBatchInteraction interaction = RegionSplitBatchInteraction::IndependentSingleOwners;

    int numComponents() const {
        return static_cast<int>(components.size());
    }

    bool supportedByCurrentUpdater() const {
        return interaction == RegionSplitBatchInteraction::IndependentSingleOwners ||
               interaction == RegionSplitBatchInteraction::SharedExistingLevelBTarget;
    }

    bool candidateForFlatLowerLocalUpdater() const {
        return supportedByCurrentUpdater() ||
               interaction == RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor ||
               interaction == RegionSplitBatchInteraction::MixedInsertedAndNonInsertedTargets ||
               interaction == RegionSplitBatchInteraction::RequiresMultipleInsertedLevelBAncestors ||
               interaction == RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors ||
               interaction == RegionSplitBatchInteraction::PreservedAncestorSplits;
    }
};

struct RegionSplitLowerNativeNodePlan {
    int level = 0;
    std::vector<PixelId> supportPixels;
    std::vector<PixelId> directProperParts;
    std::vector<NodeId> preDirectOwners;
    bool touchesRegion = false;
    bool touchesPreserved = false;
    int parentPlanIndex = -1;
    std::vector<int> childPlanIndices;
};

struct RegionSplitLowerNativeRecompositionPlan {
    RegionLevelChange<AltitudeType> change;
    std::vector<int> affectedLevels;
    std::vector<std::vector<PixelId>> levelBEnvelopeComponents;
    std::vector<RegionSplitLowerNativeNodePlan> nodes;

    int numRootNodes() const {
        return static_cast<int>(std::count_if(
            nodes.begin(),
            nodes.end(),
            [](const auto &node) {
                return node.parentPlanIndex < 0;
            }));
    }

    int numConnectorCandidates() const {
        return static_cast<int>(std::count_if(
            nodes.begin(),
            nodes.end(),
            [](const auto &node) {
                return node.directProperParts.empty();
            }));
    }
};

struct RegionSplitLowerNativeApplyResult {
    ImageUInt8Ptr reconstruction;
    RegionSplitLowerNativeRecompositionPlan plan;
    std::vector<NodeId> createdNodesByPlanIndex;
    NodeId replacedEnvelopeNode = InvalidNode;
    NodeId externalParentNode = InvalidNode;
    std::vector<NodeId> replacedEnvelopeNodes;
    std::vector<NodeId> externalParentNodes;
    int releasedEnvelopeSubtreeNodes = 0;

    int numCreatedNodes() const {
        return static_cast<int>(std::count_if(
            createdNodesByPlanIndex.begin(),
            createdNodesByPlanIndex.end(),
            [](NodeId nodeId) {
                return nodeId != InvalidNode;
            }));
    }

    int numReplacedEnvelopeNodes() const {
        return static_cast<int>(replacedEnvelopeNodes.size());
    }
};

class RegionSplitUpdater {
public:
    static RegionSplitPlan buildFlatLowerPlan(const DynamicComponentTree &maxTree,
                                              AdjacencyRelation &adjacency,
                                              const RegionLevelChange<AltitudeType> &change) {
        validateFlatLowerAssumptions(maxTree, adjacency, change);

        RegionSplitPlan plan;
        plan.change = change;
        plan.firstAffectedLevel = static_cast<int>(change.b) + 1;
        plan.lastAffectedLevel = static_cast<int>(change.a);

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);

        for (int level = plan.firstAffectedLevel; level <= plan.lastAffectedLevel; ++level) {
            for (NodeId nodeId : maxTree.getNodeSubtree(maxTree.getRoot())) {
                if (!nodeContainsMaxTreeLevel(maxTree, nodeId, level)) {
                    continue;
                }

                auto support = maxTree.getPixelsOfCC(nodeId);
                if (!intersectsRegion(support, regionSet)) {
                    continue;
                }

                RegionSplitAffectedComponent affected;
                affected.level = level;
                affected.nodeId = nodeId;
                affected.originalSupport = std::move(support);
                std::sort(affected.originalSupport.begin(), affected.originalSupport.end());
                affected.removedPixels.assign(change.pixelsC.begin(), change.pixelsC.end());
                std::sort(affected.removedPixels.begin(), affected.removedPixels.end());
                affected.remainingComponents = connectedComponentsAfterRemovingRegion(
                    adjacency,
                    affected.originalSupport,
                    regionSet);

                if (affected.remainingComponents.empty()) {
                    affected.outcome = RegionSplitLocalOutcome::Removed;
                } else if (affected.remainingComponents.size() == 1) {
                    affected.outcome = RegionSplitLocalOutcome::Contracted;
                } else {
                    affected.outcome = RegionSplitLocalOutcome::Split;
                }

                plan.affectedComponents.push_back(std::move(affected));
            }
        }

        std::sort(plan.affectedComponents.begin(),
                  plan.affectedComponents.end(),
                  [](const RegionSplitAffectedComponent &lhs,
                     const RegionSplitAffectedComponent &rhs) {
                      if (lhs.level != rhs.level) {
                          return lhs.level < rhs.level;
                      }
                      return lhs.nodeId < rhs.nodeId;
                  });
        return plan;
    }

    static RegionSplitPlan buildLowerPlan(const DynamicComponentTree &maxTree,
                                          AdjacencyRelation &adjacency,
                                          const RegionLevelChange<AltitudeType> &change) {
        validateLowerPlanAssumptions(maxTree, adjacency, change);

        RegionSplitPlan plan;
        plan.change = change;
        plan.firstAffectedLevel = static_cast<int>(change.b) + 1;
        plan.lastAffectedLevel = static_cast<int>(change.a);

        const auto reconstructed = maxTree.reconstructionImage();
        for (int level = plan.firstAffectedLevel; level <= plan.lastAffectedLevel; ++level) {
            auto removedAtLevel = removedPixelsAtOrAboveLevel(
                *reconstructed,
                change.pixelsC,
                level);
            if (removedAtLevel.empty()) {
                continue;
            }

            const std::unordered_set<PixelId> removedSet(
                removedAtLevel.begin(),
                removedAtLevel.end());

            for (NodeId nodeId : maxTree.getNodeSubtree(maxTree.getRoot())) {
                if (!nodeContainsMaxTreeLevel(maxTree, nodeId, level)) {
                    continue;
                }

                auto support = maxTree.getPixelsOfCC(nodeId);
                if (!intersectsRegion(support, removedSet)) {
                    continue;
                }

                RegionSplitAffectedComponent affected;
                affected.level = level;
                affected.nodeId = nodeId;
                affected.originalSupport = std::move(support);
                std::sort(affected.originalSupport.begin(), affected.originalSupport.end());
                affected.removedPixels = removedAtLevel;
                affected.remainingComponents = connectedComponentsAfterRemovingRegion(
                    adjacency,
                    affected.originalSupport,
                    removedSet);

                if (affected.remainingComponents.empty()) {
                    affected.outcome = RegionSplitLocalOutcome::Removed;
                } else if (affected.remainingComponents.size() == 1) {
                    affected.outcome = RegionSplitLocalOutcome::Contracted;
                } else {
                    affected.outcome = RegionSplitLocalOutcome::Split;
                }

                plan.affectedComponents.push_back(std::move(affected));
            }
        }

        std::sort(plan.affectedComponents.begin(),
                  plan.affectedComponents.end(),
                  [](const RegionSplitAffectedComponent &lhs,
                     const RegionSplitAffectedComponent &rhs) {
                      if (lhs.level != rhs.level) {
                          return lhs.level < rhs.level;
                      }
                      return lhs.nodeId < rhs.nodeId;
                  });
        return plan;
    }

    static RegionSplitLowerNativeRecompositionPlan buildLowerNativeRecompositionPlan(
        const DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateLowerPossiblyDisconnectedAssumptions(maxTree, change);
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionSplitUpdater native lower plan requires connected region C");
        }

        const int domainSize = maxTree.getNumTotalProperParts();
        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        std::vector<unsigned char> inRegion(static_cast<size_t>(domainSize), 0);
        for (PixelId pixelId : change.pixelsC) {
            inRegion[static_cast<size_t>(pixelId)] = 1;
        }

        auto updatedImage = maxTree.reconstructionImage()->clone();
        for (PixelId pixelId : change.pixelsC) {
            (*updatedImage)[pixelId] = static_cast<uint8_t>(change.b);
        }

        RegionSplitLowerNativeRecompositionPlan plan;
        plan.change = change;
        plan.affectedLevels.reserve(static_cast<size_t>(
            static_cast<int>(change.a) - static_cast<int>(change.b) + 1));
        for (int level = static_cast<int>(change.b);
             level <= static_cast<int>(change.a);
             ++level) {
            plan.affectedLevels.push_back(level);
        }

        std::vector<PixelId> levelBActivePixels;
        levelBActivePixels.reserve(static_cast<size_t>(domainSize));
        for (PixelId pixelId = 0; pixelId < domainSize; ++pixelId) {
            if (static_cast<int>((*updatedImage)[pixelId]) >= static_cast<int>(change.b)) {
                levelBActivePixels.push_back(pixelId);
            }
        }

        std::vector<unsigned char> inLocalEnvelope(static_cast<size_t>(domainSize), 0);
        auto levelBComponents = connectedComponentsInRegion(
            adjacency,
            std::span<const PixelId>(levelBActivePixels),
            domainSize);
        for (auto &component : levelBComponents) {
            if (!intersectsRegion(component, regionSet)) {
                continue;
            }
            for (PixelId pixelId : component) {
                inLocalEnvelope[static_cast<size_t>(pixelId)] = 1;
            }
            plan.levelBEnvelopeComponents.push_back(std::move(component));
        }
        if (plan.levelBEnvelopeComponents.empty()) {
            throw std::runtime_error("RegionSplitUpdater native lower plan could not locate the level-b envelope of C");
        }

        std::vector<int> recompositionLevels = plan.affectedLevels;
        for (PixelId pixelId = 0; pixelId < domainSize; ++pixelId) {
            if (inLocalEnvelope[static_cast<size_t>(pixelId)] == 0) {
                continue;
            }
            const int value = static_cast<int>((*updatedImage)[pixelId]);
            if (value > static_cast<int>(change.a)) {
                recompositionLevels.push_back(value);
            }
        }
        std::sort(recompositionLevels.begin(), recompositionLevels.end());
        recompositionLevels.erase(
            std::unique(recompositionLevels.begin(), recompositionLevels.end()),
            recompositionLevels.end());

        for (int level : recompositionLevels) {
            std::vector<PixelId> activePixels;
            activePixels.reserve(static_cast<size_t>(domainSize));
            for (PixelId pixelId = 0; pixelId < domainSize; ++pixelId) {
                if (inLocalEnvelope[static_cast<size_t>(pixelId)] == 0) {
                    continue;
                }
                if (static_cast<int>((*updatedImage)[pixelId]) >= level) {
                    activePixels.push_back(pixelId);
                }
            }

            auto components = connectedComponentsInRegion(
                adjacency,
                std::span<const PixelId>(activePixels),
                domainSize);
            for (const auto &component : components) {
                RegionSplitLowerNativeNodePlan node;
                node.level = level;
                node.supportPixels = component;

                for (PixelId pixelId : node.supportPixels) {
                    if (inRegion[static_cast<size_t>(pixelId)] != 0) {
                        node.touchesRegion = true;
                    } else {
                        node.touchesPreserved = true;
                    }

                    if (static_cast<int>((*updatedImage)[pixelId]) != level) {
                        continue;
                    }

                    node.directProperParts.push_back(pixelId);
                    const NodeId ownerId = maxTree.getSmallestComponent(pixelId);
                    if (ownerId != InvalidNode && maxTree.isAlive(ownerId)) {
                        node.preDirectOwners.push_back(ownerId);
                    }
                }

                if (node.directProperParts.empty()) {
                    continue;
                }

                std::sort(node.directProperParts.begin(), node.directProperParts.end());
                std::sort(node.preDirectOwners.begin(), node.preDirectOwners.end());
                node.preDirectOwners.erase(
                    std::unique(node.preDirectOwners.begin(), node.preDirectOwners.end()),
                    node.preDirectOwners.end());
                plan.nodes.push_back(std::move(node));
            }
        }

        std::sort(plan.nodes.begin(),
                  plan.nodes.end(),
                  [](const auto &lhs, const auto &rhs) {
                      if (lhs.level != rhs.level) {
                          return lhs.level < rhs.level;
                      }
                      return lhs.supportPixels < rhs.supportPixels;
                  });
        assignLowerNativeParentIndices(plan);
        return plan;
    }

    static RegionSplitLowerNativeApplyResult updateLowerNativeRecomposition(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        DynamicComponentTree preflightTree = maxTree;
        (void) updateLowerNativeRecompositionInPlace(preflightTree, adjacency, change);
        return updateLowerNativeRecompositionInPlace(maxTree, adjacency, change);
    }

    static RegionSplitApplyResult updateFlatLowerSingleOwner(DynamicComponentTree &maxTree,
                                                             AdjacencyRelation &adjacency,
                                                             const RegionLevelChange<AltitudeType> &change) {
        const NodeId sourceNodeId = validateFlatLowerSingleOwnerEditable(maxTree, adjacency, change);

        RegionSplitApplyResult result;
        result.sourceNode = sourceNodeId;
        result.sourceNodes.push_back(sourceNodeId);

        const bool sourceWasRoot = maxTree.isRoot(sourceNodeId);
        const NodeId oldParentId = maxTree.getNodeParent(sourceNodeId);

        if (!sourceWasRoot && oldParentId != InvalidNode && maxTree.getAltitude(oldParentId) == change.b) {
            result.targetNode = oldParentId;
            moveRegionProperParts(maxTree, result.targetNode, sourceNodeId, change.pixelsC);
            maxTree.moveChildren(result.targetNode, sourceNodeId);
            maxTree.removeChild(result.targetNode, sourceNodeId, true);
            result.releasedSourceNode = true;
            result.reconstruction = maxTree.reconstructionImage();
            return result;
        }

        if (!sourceWasRoot &&
            (oldParentId == InvalidNode || !maxTree.isAlive(oldParentId) || maxTree.getAltitude(oldParentId) > change.b)) {
            throw std::runtime_error("RegionSplitUpdater single-owner update cannot attach b below a higher parent altitude");
        }

        const NodeId targetNodeId = maxTree.createDetachedNode(static_cast<int>(change.b));
        if (targetNodeId == InvalidNode) {
            throw std::runtime_error("RegionSplitUpdater single-owner update could not create target node");
        }
        result.targetNode = targetNodeId;
        result.createdTargetNode = true;

        moveRegionProperParts(maxTree, targetNodeId, sourceNodeId, change.pixelsC);
        maxTree.moveChildren(targetNodeId, sourceNodeId);

        if (sourceWasRoot) {
            maxTree.setRoot(targetNodeId);
            maxTree.releaseNode(sourceNodeId);
        } else {
            maxTree.removeChild(oldParentId, sourceNodeId, true);
            maxTree.attachNode(oldParentId, targetNodeId);
        }
        result.releasedSourceNode = true;
        result.reconstruction = maxTree.reconstructionImage();
        return result;
    }

    static RegionSplitBatchApplyResult updateLowerByFlatLevelBatches(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateLowerPlanAssumptions(maxTree, adjacency, change);
        const auto levelChanges = buildFlatLevelChangesForLowerBatches(maxTree, change);
        if (levelChanges.empty()) {
            throw std::runtime_error("RegionSplitUpdater lower batch update requires at least one affected level");
        }

        DynamicComponentTree preflightTree = maxTree;
        for (const auto &levelChange : levelChanges) {
            (void) updateFlatLowerLocalSupportedBatch(preflightTree, adjacency, levelChange);
        }

        RegionSplitBatchApplyResult result;
        for (const auto &levelChange : levelChanges) {
            auto levelResult = updateFlatLowerLocalSupportedBatch(maxTree, adjacency, levelChange);
            result.componentResults.insert(result.componentResults.end(),
                                           levelResult.componentResults.begin(),
                                           levelResult.componentResults.end());
        }

        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static RegionSplitBatchApplyResult updateLowerByFlatLevelBatchesWithMaterialization(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateLowerPossiblyDisconnectedAssumptions(maxTree, change);
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionSplitUpdater lower materialized batch update requires connected region C");
        }

        auto levelChanges = buildFlatLevelChangesForLowerBatches(maxTree, change);
        if (levelChanges.empty()) {
            throw std::runtime_error("RegionSplitUpdater lower materialized batch update requires at least one affected level");
        }
        for (auto &levelChange : levelChanges) {
            levelChange.assumeConnected = true;
            levelChange.assumeMaterializedAsProperParts = false;
        }

        DynamicComponentTree preflightTree = maxTree;
        for (const auto &levelChange : levelChanges) {
            (void) updateFlatLowerWithMaterialization(preflightTree, adjacency, levelChange);
        }

        RegionSplitBatchApplyResult result;
        for (const auto &levelChange : levelChanges) {
            auto levelResult = updateFlatLowerWithMaterialization(maxTree, adjacency, levelChange);
            result.componentResults.insert(result.componentResults.end(),
                                           levelResult.componentResults.begin(),
                                           levelResult.componentResults.end());
        }

        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static RegionSplitMaterializationResult materializeFlatLowerRegion(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatLowerPossiblyDisconnectedAssumptions(maxTree, change);
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionSplitUpdater materialization requires connected region C");
        }

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        RegionSplitMaterializationResult result;
        if (isMaterializedAsProperParts(maxTree, regionSet)) {
            result.reconstruction = maxTree.reconstructionImage();
            return result;
        }

        const auto materialization = analyzeRegionMaterialization(maxTree, adjacency, change);
        if (materialization.numPartialOwners() == 0) {
            throw std::runtime_error("RegionSplitUpdater materialization found a non-materialized region without partial owners");
        }

        DynamicComponentTree preflightTree = maxTree;
        for (const auto &owner : materialization.affectedOwners) {
            if (!owner.cutsDirectProperParts()) {
                continue;
            }
            (void) materializeFlatLowerSinglePartialOwnerAsSource(
                preflightTree,
                adjacency,
                owner,
                regionSet);
        }
        if (!isMaterializedAsProperParts(preflightTree, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater materialization preflight could not make C a union of complete direct proper parts");
        }

        for (const auto &owner : materialization.affectedOwners) {
            if (!owner.cutsDirectProperParts()) {
                continue;
            }
            auto ownerResult = materializeFlatLowerSinglePartialOwnerAsSource(
                maxTree,
                adjacency,
                owner,
                regionSet);
            appendMaterializationResult(result, ownerResult);
        }

        if (!isMaterializedAsProperParts(maxTree, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater materialization could not make C a union of complete direct proper parts");
        }
        result.reconstruction = maxTree.reconstructionImage();
        return result;
    }

    static RegionSplitBatchApplyResult updateFlatLowerWithMaterialization(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatLowerPossiblyDisconnectedAssumptions(maxTree, change);
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionSplitUpdater materialized update requires connected region C");
        }

        RegionLevelChange<AltitudeType> materializedChange = change;
        materializedChange.assumeMaterializedAsProperParts = false;

        DynamicComponentTree preflightTree = maxTree;
        (void) materializeFlatLowerRegion(preflightTree, adjacency, materializedChange);
        (void) updateFlatLowerMaterializedDispatch(preflightTree, adjacency, materializedChange);

        (void) materializeFlatLowerRegion(maxTree, adjacency, materializedChange);
        return updateFlatLowerMaterializedDispatch(maxTree, adjacency, materializedChange);
    }

    static RegionSplitBatchPlan buildFlatLowerBatchPlan(
        const DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatLowerPossiblyDisconnectedAssumptions(maxTree, change);
        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        if (!change.assumeMaterializedAsProperParts &&
            !isMaterializedAsProperParts(maxTree, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater batch plan requires C to be materialized as complete direct proper parts");
        }

        auto components = connectedComponentsInRegion(
            adjacency,
            change.pixelsC,
            maxTree.getNumTotalProperParts());
        if (components.empty()) {
            throw std::runtime_error("RegionSplitUpdater batch plan requires a non-empty region C");
        }

        RegionSplitBatchPlan plan;
        plan.change = change;
        plan.components.reserve(components.size());

        std::unordered_set<NodeId> seenSources;
        bool hasRepeatedSource = false;
        bool hasInsertedLevelBAncestor = false;
        bool hasSharedExistingTarget = false;
        std::unordered_set<NodeId> seenExistingTargets;

        for (const auto &component : components) {
            const NodeId sourceNodeId = findSingleDirectOwner(maxTree, component);
            if (maxTree.getAltitude(sourceNodeId) != change.a) {
                throw std::runtime_error("RegionSplitUpdater batch plan requires each owner altitude to be a");
            }

            RegionSplitBatchComponentPlan componentPlan;
            componentPlan.pixelsC = component;
            componentPlan.sourceNode = sourceNodeId;
            componentPlan.sourceIsRoot = maxTree.isRoot(sourceNodeId);
            componentPlan.parentNode = componentPlan.sourceIsRoot
                ? InvalidNode
                : maxTree.getNodeParent(sourceNodeId);

            if (!seenSources.insert(sourceNodeId).second) {
                hasRepeatedSource = true;
            }

            if (!componentPlan.sourceIsRoot) {
                if (componentPlan.parentNode == InvalidNode ||
                    !maxTree.isAlive(componentPlan.parentNode)) {
                    throw std::runtime_error("RegionSplitUpdater batch plan requires every non-root source to have a live parent");
                }

                const int parentAltitude = maxTree.getAltitude(componentPlan.parentNode);
                if (parentAltitude == change.b) {
                    componentPlan.reusesExistingLevelBTarget = true;
                    if (!seenExistingTargets.insert(componentPlan.parentNode).second) {
                        hasSharedExistingTarget = true;
                    }
                } else if (parentAltitude > change.b) {
                    componentPlan.requiresInsertedLevelBAncestor = true;
                    hasInsertedLevelBAncestor = true;
                }
            }

            plan.components.push_back(std::move(componentPlan));
        }

        if (hasRepeatedSource) {
            plan.interaction = RegionSplitBatchInteraction::RepeatedSourceNode;
            return plan;
        }

        for (const auto &component : plan.components) {
            const auto componentSet = makeRegionSet(maxTree, component.pixelsC);
            if (!nodeProperPartsEqualRegion(maxTree, component.sourceNode, componentSet)) {
                throw std::runtime_error("RegionSplitUpdater batch plan requires each C component to be the complete direct proper part of one node");
            }
        }

        if (hasInsertedLevelBAncestor) {
            plan.interaction = classifyInsertedLevelBAncestorInteraction(
                maxTree,
                adjacency,
                plan,
                regionSet);
        } else if (hasSharedExistingTarget) {
            plan.interaction = RegionSplitBatchInteraction::SharedExistingLevelBTarget;
        } else {
            plan.interaction = RegionSplitBatchInteraction::IndependentSingleOwners;
        }
        return plan;
    }

    static RegionSplitBatchApplyResult updateFlatLowerIndependentSingleOwners(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerIndependentSingleOwners(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerSharedInsertedLevelBAncestor(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerSharedInsertedLevelBAncestor(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerPreservedAncestorSplit(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerPreservedAncestorSplit(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerMultipleInsertedLevelBAncestors(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerMultipleInsertedLevelBAncestors(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerMixedInsertedAndNonInsertedTargets(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerMixedInsertedAndNonInsertedTargets(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerLocalSupportedBatch(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerBatchFromPlan(maxTree, adjacency, plan);
    }

    static RegionSplitBatchPlan buildFlatLowerSourceGroupedBatchPlan(
        const DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        validateFlatLowerPossiblyDisconnectedAssumptions(maxTree, change);
        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        if (!isMaterializedAsProperParts(maxTree, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater source-grouped plan requires C to be materialized as complete direct proper parts");
        }

        std::vector<NodeId> sourceOrder;
        std::vector<std::vector<PixelId>> pixelsBySource;
        std::unordered_map<NodeId, size_t> sourceIndex;
        sourceIndex.reserve(change.pixelsC.size());

        for (PixelId pixelId : change.pixelsC) {
            const NodeId sourceNodeId = maxTree.getSmallestComponent(pixelId);
            if (sourceNodeId == InvalidNode || !maxTree.isAlive(sourceNodeId)) {
                throw std::runtime_error("RegionSplitUpdater source-grouped plan cannot locate a live source for every pixel in C");
            }

            auto it = sourceIndex.find(sourceNodeId);
            if (it == sourceIndex.end()) {
                const size_t newIndex = sourceOrder.size();
                sourceIndex.emplace(sourceNodeId, newIndex);
                sourceOrder.push_back(sourceNodeId);
                pixelsBySource.emplace_back();
                it = sourceIndex.find(sourceNodeId);
            }
            pixelsBySource[it->second].push_back(pixelId);
        }

        if (sourceOrder.empty()) {
            throw std::runtime_error("RegionSplitUpdater source-grouped plan requires a non-empty region C");
        }

        RegionSplitBatchPlan plan;
        plan.change = change;
        plan.components.reserve(sourceOrder.size());

        bool hasInsertedLevelBAncestor = false;
        bool hasSharedExistingTarget = false;
        std::unordered_set<NodeId> seenExistingTargets;

        for (size_t i = 0; i < sourceOrder.size(); ++i) {
            const NodeId sourceNodeId = sourceOrder[i];
            auto &pixels = pixelsBySource[i];
            std::sort(pixels.begin(), pixels.end());

            if (maxTree.getAltitude(sourceNodeId) != change.a) {
                throw std::runtime_error("RegionSplitUpdater source-grouped plan requires each source altitude to be a");
            }

            const auto sourceRegionSet = makeRegionSet(maxTree, std::span<const PixelId>(pixels));
            if (!nodeProperPartsEqualRegion(maxTree, sourceNodeId, sourceRegionSet)) {
                throw std::runtime_error("RegionSplitUpdater source-grouped plan requires each source group to be the complete direct proper part of its source node");
            }

            RegionSplitBatchComponentPlan componentPlan;
            componentPlan.pixelsC = pixels;
            componentPlan.sourceNode = sourceNodeId;
            componentPlan.sourceIsRoot = maxTree.isRoot(sourceNodeId);
            componentPlan.parentNode = componentPlan.sourceIsRoot
                ? InvalidNode
                : maxTree.getNodeParent(sourceNodeId);

            if (!componentPlan.sourceIsRoot) {
                if (componentPlan.parentNode == InvalidNode ||
                    !maxTree.isAlive(componentPlan.parentNode)) {
                    throw std::runtime_error("RegionSplitUpdater source-grouped plan requires every non-root source to have a live parent");
                }

                const int parentAltitude = maxTree.getAltitude(componentPlan.parentNode);
                if (parentAltitude == change.b) {
                    componentPlan.reusesExistingLevelBTarget = true;
                    if (!seenExistingTargets.insert(componentPlan.parentNode).second) {
                        hasSharedExistingTarget = true;
                    }
                } else if (parentAltitude > change.b) {
                    componentPlan.requiresInsertedLevelBAncestor = true;
                    hasInsertedLevelBAncestor = true;
                }
            }

            plan.components.push_back(std::move(componentPlan));
        }

        if (hasInsertedLevelBAncestor) {
            plan.interaction = classifyInsertedLevelBAncestorInteraction(
                maxTree,
                adjacency,
                plan,
                regionSet);
        } else if (hasSharedExistingTarget) {
            plan.interaction = RegionSplitBatchInteraction::SharedExistingLevelBTarget;
        } else {
            plan.interaction = RegionSplitBatchInteraction::IndependentSingleOwners;
        }
        return plan;
    }

    static RegionSplitBatchApplyResult updateFlatLowerSourceGroupedBatch(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        const auto plan = buildFlatLowerSourceGroupedBatchPlan(maxTree, adjacency, change);
        return updateFlatLowerBatchFromPlan(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerBatchFromPlan(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        switch (plan.interaction) {
        case RegionSplitBatchInteraction::IndependentSingleOwners:
        case RegionSplitBatchInteraction::SharedExistingLevelBTarget:
            return updateFlatLowerIndependentSingleOwners(maxTree, adjacency, plan);
        case RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor:
            return updateFlatLowerSharedInsertedLevelBAncestor(maxTree, adjacency, plan);
        case RegionSplitBatchInteraction::PreservedAncestorSplits:
            return updateFlatLowerPreservedAncestorSplit(maxTree, adjacency, plan);
        case RegionSplitBatchInteraction::RequiresMultipleInsertedLevelBAncestors:
            return updateFlatLowerMultipleInsertedLevelBAncestors(maxTree, adjacency, plan);
        case RegionSplitBatchInteraction::MixedInsertedAndNonInsertedTargets:
            return updateFlatLowerMixedInsertedAndNonInsertedTargets(maxTree, adjacency, plan);
        case RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors:
            return updateFlatLowerInteractingInsertedLevelBAncestors(maxTree, adjacency, plan);
        case RegionSplitBatchInteraction::RepeatedSourceNode:
            throw std::runtime_error(unsupportedBatchInteractionMessage(plan.interaction));
        }

        throw std::runtime_error("RegionSplitUpdater local supported batch encountered an unknown interaction");
    }

    static RegionSplitBatchApplyResult updateFlatLowerIndependentSingleOwners(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        if (!plan.supportedByCurrentUpdater()) {
            throw std::runtime_error(unsupportedBatchInteractionMessage(plan.interaction));
        }

        std::vector<RegionLevelChange<AltitudeType>> componentChanges;
        componentChanges.reserve(plan.components.size());
        for (const auto &component : plan.components) {
            componentChanges.push_back(buildExplicitRegionChange<AltitudeType>(
                std::span<const PixelId>(component.pixelsC),
                plan.change.a,
                plan.change.b,
                true,
                plan.change.assumeMaterializedAsProperParts));
            (void) validateFlatLowerSingleOwnerEditable(maxTree, adjacency, componentChanges.back());
        }

        RegionSplitBatchApplyResult result;
        result.componentResults.reserve(componentChanges.size());
        for (const auto &componentChange : componentChanges) {
            result.componentResults.push_back(updateFlatLowerSingleOwner(
                maxTree,
                adjacency,
                componentChange));
        }
        result.reconstruction = maxTree.reconstructionImage();
        return result;
    }

    static void validateFlatLowerAssumptions(const DynamicComponentTree &maxTree,
                                             AdjacencyRelation &adjacency,
                                             const RegionLevelChange<AltitudeType> &change) {
        validateFlatLowerPossiblyDisconnectedAssumptions(maxTree, change);
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionSplitUpdater requires connected region C");
        }

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        if (!change.assumeMaterializedAsProperParts &&
            !isMaterializedAsProperParts(maxTree, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater flat plan requires C to be materialized as complete direct proper parts");
        }
    }

    static bool isMaterializedAsProperParts(const DynamicComponentTree &tree,
                                            const std::unordered_set<PixelId> &regionSet) {
        for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
            bool touchesRegion = false;
            bool touchesOutside = false;
            for (PixelId pixelId : tree.getProperParts(nodeId)) {
                if (regionSet.find(pixelId) != regionSet.end()) {
                    touchesRegion = true;
                } else {
                    touchesOutside = true;
                }
                if (touchesRegion && touchesOutside) {
                    return false;
                }
            }
        }
        return true;
    }

private:
    static int countSubtreeNodes(const DynamicComponentTree &tree, NodeId subtreeRoot) {
        int count = 0;
        for (NodeId ignored : tree.getNodeSubtree(subtreeRoot)) {
            (void) ignored;
            ++count;
        }
        return count;
    }

    static NodeId findUniqueNodeWithSupport(const DynamicComponentTree &tree,
                                            std::span<const PixelId> supportPixels) {
        NodeId match = InvalidNode;
        for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
            auto nodeSupport = tree.getPixelsOfCC(nodeId);
            std::sort(nodeSupport.begin(), nodeSupport.end());
            if (!std::equal(nodeSupport.begin(),
                            nodeSupport.end(),
                            supportPixels.begin(),
                            supportPixels.end())) {
                continue;
            }
            if (match != InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater native lower update found multiple pre-update nodes with the envelope support");
            }
            match = nodeId;
        }
        return match;
    }

    static std::vector<int> collectLowerNativeRootPlanIndices(
        const RegionSplitLowerNativeRecompositionPlan &plan) {
        std::vector<int> rootIndices;
        for (size_t i = 0; i < plan.nodes.size(); ++i) {
            if (plan.nodes[i].parentPlanIndex >= 0) {
                continue;
            }
            rootIndices.push_back(static_cast<int>(i));
        }
        if (rootIndices.empty()) {
            throw std::runtime_error("RegionSplitUpdater native lower update requires a local root node");
        }
        return rootIndices;
    }

    static int findEnvelopeIndexForLowerNativeRoot(
        const RegionSplitLowerNativeRecompositionPlan &plan,
        int rootPlanIndex,
        std::vector<unsigned char> &usedEnvelopeIndices) {
        if (rootPlanIndex < 0 ||
            static_cast<size_t>(rootPlanIndex) >= plan.nodes.size()) {
            throw std::runtime_error("RegionSplitUpdater native lower update received an invalid root plan index");
        }

        const auto &rootPlanNode = plan.nodes[static_cast<size_t>(rootPlanIndex)];
        for (size_t envelopeIndex = 0;
             envelopeIndex < plan.levelBEnvelopeComponents.size();
             ++envelopeIndex) {
            if (usedEnvelopeIndices[envelopeIndex] != 0) {
                continue;
            }
            if (rootPlanNode.supportPixels != plan.levelBEnvelopeComponents[envelopeIndex]) {
                continue;
            }
            usedEnvelopeIndices[envelopeIndex] = 1;
            return static_cast<int>(envelopeIndex);
        }
        throw std::runtime_error("RegionSplitUpdater native lower update could not match a local root to a level-b envelope");
    }

    static RegionSplitLowerNativeApplyResult updateLowerNativeRecompositionInPlace(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        auto plan = buildLowerNativeRecompositionPlan(maxTree, adjacency, change);
        if (plan.numConnectorCandidates() != 0) {
            throw std::runtime_error("RegionSplitUpdater native lower update refuses connector candidates without direct proper parts");
        }

        const auto rootPlanIndices = collectLowerNativeRootPlanIndices(plan);
        if (rootPlanIndices.size() != plan.levelBEnvelopeComponents.size()) {
            throw std::runtime_error("RegionSplitUpdater native lower update requires one local root per level-b envelope");
        }

        std::vector<int> envelopeIndexByRootPlanIndex(plan.nodes.size(), -1);
        std::vector<unsigned char> usedEnvelopeIndices(plan.levelBEnvelopeComponents.size(), 0);
        for (int rootPlanIndex : rootPlanIndices) {
            const auto &rootPlanNode = plan.nodes[static_cast<size_t>(rootPlanIndex)];
            if (rootPlanNode.level != static_cast<int>(change.b)) {
                throw std::runtime_error("RegionSplitUpdater native lower update requires every local root to be at level b");
            }
            const int envelopeIndex = findEnvelopeIndexForLowerNativeRoot(
                plan,
                rootPlanIndex,
                usedEnvelopeIndices);
            envelopeIndexByRootPlanIndex[static_cast<size_t>(rootPlanIndex)] = envelopeIndex;
        }

        const int domainSize = maxTree.getNumTotalProperParts();
        std::vector<int> directPlanIndexByPixel(static_cast<size_t>(domainSize), -1);
        for (size_t planIndex = 0; planIndex < plan.nodes.size(); ++planIndex) {
            for (PixelId pixelId : plan.nodes[planIndex].directProperParts) {
                auto &owner = directPlanIndexByPixel[static_cast<size_t>(pixelId)];
                if (owner != -1) {
                    throw std::runtime_error("RegionSplitUpdater native lower update found duplicated direct proper parts in the native plan");
                }
                owner = static_cast<int>(planIndex);
            }
        }
        for (const auto &envelopeComponent : plan.levelBEnvelopeComponents) {
            for (PixelId pixelId : envelopeComponent) {
                if (directPlanIndexByPixel[static_cast<size_t>(pixelId)] < 0) {
                    throw std::runtime_error("RegionSplitUpdater native lower update currently requires every envelope pixel to be a direct proper part of the native plan");
                }
            }
        }

        struct LowerNativeEnvelopeEdit {
            int rootPlanIndex = -1;
            int envelopeIndex = -1;
            NodeId envelopeNode = InvalidNode;
            NodeId externalParent = InvalidNode;
            bool envelopeWasRoot = false;
            std::vector<NodeId> oldChildren;
            int releasedSubtreeNodes = 0;
        };

        std::vector<LowerNativeEnvelopeEdit> envelopeEdits;
        envelopeEdits.reserve(rootPlanIndices.size());
        int rootEnvelopeCount = 0;
        for (int rootPlanIndex : rootPlanIndices) {
            const int envelopeIndex = envelopeIndexByRootPlanIndex[static_cast<size_t>(rootPlanIndex)];
            if (envelopeIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater native lower update found a local root without an envelope");
            }

            LowerNativeEnvelopeEdit edit;
            edit.rootPlanIndex = rootPlanIndex;
            edit.envelopeIndex = envelopeIndex;
            edit.envelopeNode = findUniqueNodeWithSupport(
                maxTree,
                std::span<const PixelId>(plan.levelBEnvelopeComponents[static_cast<size_t>(envelopeIndex)]));
            if (edit.envelopeNode == InvalidNode || !maxTree.isAlive(edit.envelopeNode)) {
                throw std::runtime_error("RegionSplitUpdater native lower update could not locate a pre-update envelope node");
            }

            edit.envelopeWasRoot = maxTree.isRoot(edit.envelopeNode);
            if (edit.envelopeWasRoot) {
                ++rootEnvelopeCount;
                edit.externalParent = InvalidNode;
            } else {
                edit.externalParent = maxTree.getNodeParent(edit.envelopeNode);
                if (edit.externalParent == InvalidNode || !maxTree.isAlive(edit.externalParent)) {
                    throw std::runtime_error("RegionSplitUpdater native lower update requires a live external parent");
                }
                if (maxTree.getAltitude(edit.externalParent) >= static_cast<int>(change.b)) {
                    throw std::runtime_error("RegionSplitUpdater native lower update requires every external parent altitude to be below b");
                }
            }

            for (NodeId childId : maxTree.getChildren(edit.envelopeNode)) {
                edit.oldChildren.push_back(childId);
            }
            edit.releasedSubtreeNodes =
                countSubtreeNodes(maxTree, edit.envelopeNode) - (edit.envelopeWasRoot ? 1 : 0);
            envelopeEdits.push_back(std::move(edit));
        }
        if (rootEnvelopeCount > 1 ||
            (rootEnvelopeCount == 1 && envelopeEdits.size() > 1)) {
            throw std::runtime_error("RegionSplitUpdater native lower update cannot combine a root envelope with other envelopes");
        }

        RegionSplitLowerNativeApplyResult result;
        result.plan = std::move(plan);
        result.createdNodesByPlanIndex.assign(result.plan.nodes.size(), InvalidNode);
        for (const auto &edit : envelopeEdits) {
            result.replacedEnvelopeNodes.push_back(edit.envelopeNode);
            result.externalParentNodes.push_back(edit.externalParent);
            result.releasedEnvelopeSubtreeNodes += edit.releasedSubtreeNodes;
        }
        if (envelopeEdits.size() == 1) {
            result.replacedEnvelopeNode = envelopeEdits.front().envelopeNode;
            result.externalParentNode = envelopeEdits.front().externalParent;
        }

        for (size_t planIndex = 0; planIndex < result.plan.nodes.size(); ++planIndex) {
            bool reusesRootEnvelope = false;
            for (const auto &edit : envelopeEdits) {
                if (edit.envelopeWasRoot &&
                    static_cast<int>(planIndex) == edit.rootPlanIndex) {
                    maxTree.setAltitude(edit.envelopeNode, result.plan.nodes[planIndex].level);
                    result.createdNodesByPlanIndex[planIndex] = edit.envelopeNode;
                    reusesRootEnvelope = true;
                    break;
                }
            }
            if (reusesRootEnvelope) {
                continue;
            }

            const NodeId newNodeId = maxTree.createDetachedNode(result.plan.nodes[planIndex].level);
            if (newNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater native lower update could not create a recomposition node");
            }
            result.createdNodesByPlanIndex[planIndex] = newNodeId;
        }

        for (size_t planIndex = 0; planIndex < result.plan.nodes.size(); ++planIndex) {
            const NodeId targetNodeId = result.createdNodesByPlanIndex[planIndex];
            for (PixelId pixelId : result.plan.nodes[planIndex].directProperParts) {
                const NodeId sourceOwnerId = maxTree.getSmallestComponent(pixelId);
                if (sourceOwnerId == InvalidNode || !maxTree.isAlive(sourceOwnerId)) {
                    throw std::runtime_error("RegionSplitUpdater native lower update cannot locate a live source owner");
                }
                if (sourceOwnerId == targetNodeId) {
                    continue;
                }
                maxTree.moveProperPart(targetNodeId, sourceOwnerId, pixelId);
            }
        }

        for (const auto &edit : envelopeEdits) {
            if (edit.envelopeWasRoot) {
                for (NodeId childId : edit.oldChildren) {
                    if (childId != InvalidNode &&
                        maxTree.isAlive(childId) &&
                        maxTree.getNodeParent(childId) == edit.envelopeNode) {
                        maxTree.pruneNode(childId);
                    }
                }
                maxTree.setRoot(edit.envelopeNode);
            } else {
                maxTree.pruneNode(edit.envelopeNode);
                maxTree.attachNode(
                    edit.externalParent,
                    result.createdNodesByPlanIndex[static_cast<size_t>(edit.rootPlanIndex)]);
            }
        }

        for (size_t planIndex = 0; planIndex < result.plan.nodes.size(); ++planIndex) {
            const int parentPlanIndex = result.plan.nodes[planIndex].parentPlanIndex;
            if (parentPlanIndex < 0) {
                continue;
            }
            maxTree.attachNode(
                result.createdNodesByPlanIndex[static_cast<size_t>(parentPlanIndex)],
                result.createdNodesByPlanIndex[planIndex]);
        }

        result.reconstruction = maxTree.reconstructionImage();
        return result;
    }

    static bool containsAllSortedPixels(std::span<const PixelId> container,
                                        std::span<const PixelId> contained) {
        return std::includes(
            container.begin(),
            container.end(),
            contained.begin(),
            contained.end());
    }

    static void assignLowerNativeParentIndices(RegionSplitLowerNativeRecompositionPlan &plan) {
        for (auto &node : plan.nodes) {
            node.parentPlanIndex = -1;
            node.childPlanIndices.clear();
        }

        for (size_t childIndex = 0; childIndex < plan.nodes.size(); ++childIndex) {
            const auto &child = plan.nodes[childIndex];
            int bestParentIndex = -1;
            size_t bestParentSupportSize = 0;
            int bestLevelDistance = 0;

            for (size_t parentIndex = 0; parentIndex < plan.nodes.size(); ++parentIndex) {
                if (parentIndex == childIndex) {
                    continue;
                }

                const auto &candidate = plan.nodes[parentIndex];
                if (candidate.level >= child.level) {
                    continue;
                }
                if (candidate.supportPixels.size() <= child.supportPixels.size()) {
                    continue;
                }
                if (!containsAllSortedPixels(
                        std::span<const PixelId>(candidate.supportPixels),
                        std::span<const PixelId>(child.supportPixels))) {
                    continue;
                }

                const int levelDistance = child.level - candidate.level;
                if (bestParentIndex == -1 ||
                    candidate.supportPixels.size() < bestParentSupportSize ||
                    (candidate.supportPixels.size() == bestParentSupportSize &&
                     levelDistance < bestLevelDistance)) {
                    bestParentIndex = static_cast<int>(parentIndex);
                    bestParentSupportSize = candidate.supportPixels.size();
                    bestLevelDistance = levelDistance;
                }
            }

            plan.nodes[childIndex].parentPlanIndex = bestParentIndex;
        }

        for (size_t childIndex = 0; childIndex < plan.nodes.size(); ++childIndex) {
            const int parentIndex = plan.nodes[childIndex].parentPlanIndex;
            if (parentIndex >= 0) {
                plan.nodes[static_cast<size_t>(parentIndex)].childPlanIndices.push_back(
                    static_cast<int>(childIndex));
            }
        }
    }

    static void validateLowerPlanAssumptions(const DynamicComponentTree &maxTree,
                                             AdjacencyRelation &adjacency,
                                             const RegionLevelChange<AltitudeType> &change) {
        validateLowerPossiblyDisconnectedAssumptions(maxTree, change);
        if (!change.assumeConnected &&
            !isConnectedRegion(adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionSplitUpdater lower plan requires connected region C");
        }

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        if (!change.assumeMaterializedAsProperParts &&
            !isMaterializedAsProperParts(maxTree, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater lower plan requires C to be materialized as complete direct proper parts");
        }
    }

    static void validateLowerPossiblyDisconnectedAssumptions(
        const DynamicComponentTree &maxTree,
        const RegionLevelChange<AltitudeType> &change) {
        if (!maxTree.isMaxtree()) {
            throw std::runtime_error("RegionSplitUpdater lower plan requires a max-tree");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("RegionSplitUpdater lower plan requires a non-empty region C");
        }
        if (change.direction != RegionUpdateDirection::Lower || !(change.b < change.a)) {
            throw std::runtime_error("RegionSplitUpdater lower plan currently supports only lower updates b < a");
        }

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        if (regionSet.size() != change.pixelsC.size()) {
            throw std::runtime_error("RegionSplitUpdater lower plan requires C without duplicate pixels");
        }

        const auto reconstructed = maxTree.reconstructionImage();
        AltitudeType maxRegionValue = static_cast<AltitudeType>((*reconstructed)[change.pixelsC.front()]);
        for (PixelId pixelId : change.pixelsC) {
            const auto value = static_cast<AltitudeType>((*reconstructed)[pixelId]);
            if (!(change.b < value)) {
                throw std::runtime_error("RegionSplitUpdater lower plan requires every pixel in C to be above b");
            }
            if (value > maxRegionValue) {
                maxRegionValue = value;
            }
        }
        if (maxRegionValue != change.a) {
            throw std::runtime_error("RegionSplitUpdater lower plan requires a to be max I(C)");
        }
    }

    static std::vector<PixelId> removedPixelsAtOrAboveLevel(const ImageUInt8 &image,
                                                            std::span<const PixelId> pixelsC,
                                                            int level) {
        std::vector<PixelId> removed;
        removed.reserve(pixelsC.size());
        for (PixelId pixelId : pixelsC) {
            if (static_cast<int>(image[pixelId]) >= level) {
                removed.push_back(pixelId);
            }
        }
        std::sort(removed.begin(), removed.end());
        return removed;
    }

    static RegionLevelChange<AltitudeType> buildChangeForPixelsLike(
        std::vector<PixelId> pixels,
        const RegionLevelChange<AltitudeType> &referenceChange) {
        return buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixels),
            referenceChange.a,
            referenceChange.b,
            false,
            referenceChange.assumeMaterializedAsProperParts);
    }

    static std::vector<RegionLevelChange<AltitudeType>> buildFlatLevelChangesForLowerBatches(
        const DynamicComponentTree &maxTree,
        const RegionLevelChange<AltitudeType> &change) {
        const auto reconstructed = maxTree.reconstructionImage();
        std::vector<std::vector<PixelId>> pixelsByLevel(256);
        for (PixelId pixelId : change.pixelsC) {
            const int value = static_cast<int>((*reconstructed)[pixelId]);
            if (value > static_cast<int>(change.b)) {
                pixelsByLevel[static_cast<size_t>(value)].push_back(pixelId);
            }
        }

        std::vector<RegionLevelChange<AltitudeType>> levelChanges;
        for (int level = static_cast<int>(change.b) + 1;
             level <= static_cast<int>(change.a);
             ++level) {
            auto &pixels = pixelsByLevel[static_cast<size_t>(level)];
            if (pixels.empty()) {
                continue;
            }
            std::sort(pixels.begin(), pixels.end());
            levelChanges.push_back(buildExplicitRegionChange<AltitudeType>(
                std::span<const PixelId>(pixels),
                static_cast<AltitudeType>(level),
                change.b,
                false,
                change.assumeMaterializedAsProperParts));
        }
        return levelChanges;
    }

    static RegionSplitMaterializationResult materializeFlatLowerSinglePartialOwnerAsSource(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionOwnerMaterialization &owner,
        const std::unordered_set<PixelId> &regionSet) {
        if (!owner.cutsDirectProperParts()) {
            throw std::runtime_error("RegionSplitUpdater materialization requires a partial owner");
        }
        if (owner.nodeId == InvalidNode || !maxTree.isAlive(owner.nodeId)) {
            throw std::runtime_error("RegionSplitUpdater materialization requires a live partial owner");
        }
        if (maxTree.getAltitude(owner.nodeId) != owner.altitude) {
            throw std::runtime_error("RegionSplitUpdater materialization received a stale partial-owner plan");
        }

        std::vector<NodeId> originalChildren;
        for (NodeId childId : maxTree.getChildren(owner.nodeId)) {
            originalChildren.push_back(childId);
        }

        auto remainingComponents = connectedComponentsAfterRemovingRegion(
            adjacency,
            maxTree.getPixelsOfCC(owner.nodeId),
            regionSet);
        if (remainingComponents.empty()) {
            throw std::runtime_error("RegionSplitUpdater materialization requires a non-empty complement of C");
        }

        std::vector<int> componentIndexByPixel(
            static_cast<size_t>(maxTree.getNumTotalProperParts()),
            -1);
        for (size_t componentIndex = 0; componentIndex < remainingComponents.size(); ++componentIndex) {
            for (PixelId pixelId : remainingComponents[componentIndex]) {
                componentIndexByPixel[static_cast<size_t>(pixelId)] = static_cast<int>(componentIndex);
            }
        }

        std::vector<int> numOutsideProperPartsByComponent(remainingComponents.size(), 0);
        for (PixelId pixelId : owner.outsideDirectProperParts) {
            const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
            if (componentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater materialization found an outside proper part outside X minus C");
            }
            ++numOutsideProperPartsByComponent[static_cast<size_t>(componentIndex)];
        }

        auto locateRemainingComponent = [&](NodeId nodeId) {
            int selectedComponentIndex = -1;
            for (PixelId pixelId : maxTree.getPixelsOfCC(nodeId)) {
                if (regionSet.find(pixelId) != regionSet.end()) {
                    throw std::runtime_error("RegionSplitUpdater materialization cannot reattach a child that intersects C");
                }
                const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
                if (componentIndex < 0) {
                    throw std::runtime_error("RegionSplitUpdater materialization found a child outside X minus C");
                }
                if (selectedComponentIndex == -1) {
                    selectedComponentIndex = componentIndex;
                } else if (selectedComponentIndex != componentIndex) {
                    throw std::runtime_error("RegionSplitUpdater materialization cannot reattach a child spanning multiple remaining components");
                }
            }
            if (selectedComponentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater materialization cannot reattach an empty child");
            }
            return selectedComponentIndex;
        };

        std::vector<std::pair<NodeId, int>> childMoves;
        childMoves.reserve(originalChildren.size());
        std::vector<int> numChildrenByComponent(remainingComponents.size(), 0);
        for (NodeId childId : originalChildren) {
            const int componentIndex = locateRemainingComponent(childId);
            childMoves.push_back({childId, componentIndex});
            ++numChildrenByComponent[static_cast<size_t>(componentIndex)];
        }

        for (size_t componentIndex = 0; componentIndex < remainingComponents.size(); ++componentIndex) {
            if (numOutsideProperPartsByComponent[componentIndex] == 0 &&
                numChildrenByComponent[componentIndex] != 1) {
                throw std::runtime_error("RegionSplitUpdater materialization requires components without direct proper parts to contain exactly one preserved child");
            }
        }

        std::vector<NodeId> componentTargetNodes(remainingComponents.size(), InvalidNode);
        RegionSplitMaterializationResult result;
        result.originalOwner = owner.nodeId;
        result.originalOwners.push_back(owner.nodeId);
        result.changed = true;

        for (size_t componentIndex = 0; componentIndex < remainingComponents.size(); ++componentIndex) {
            if (numOutsideProperPartsByComponent[componentIndex] == 0) {
                continue;
            }
            const NodeId splitNodeId = maxTree.createDetachedNode(owner.altitude);
            if (splitNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater materialization could not create an outside split node");
            }
            componentTargetNodes[componentIndex] = splitNodeId;
            result.createdOutsideNodes.push_back(splitNodeId);
            maxTree.attachNode(owner.nodeId, splitNodeId);
        }

        for (PixelId pixelId : owner.outsideDirectProperParts) {
            const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
            const NodeId splitNodeId = componentTargetNodes[static_cast<size_t>(componentIndex)];
            if (splitNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater materialization cannot move an outside proper part without a split node");
            }
            maxTree.moveProperPart(splitNodeId, owner.nodeId, pixelId);
        }

        for (const auto &childMove : childMoves) {
            const NodeId splitNodeId = componentTargetNodes[static_cast<size_t>(childMove.second)];
            if (splitNodeId != InvalidNode) {
                maxTree.moveNode(childMove.first, splitNodeId);
            }
        }

        if (maxTree.getNumProperParts(owner.nodeId) != static_cast<int>(owner.insideDirectProperParts.size())) {
            throw std::runtime_error("RegionSplitUpdater materialization must leave the original owner with exactly C as direct proper parts");
        }
        for (PixelId pixelId : maxTree.getProperParts(owner.nodeId)) {
            if (regionSet.find(pixelId) == regionSet.end()) {
                throw std::runtime_error("RegionSplitUpdater materialization left an outside proper part in the original owner");
            }
        }
        return result;
    }

    static void appendMaterializationResult(RegionSplitMaterializationResult &target,
                                            const RegionSplitMaterializationResult &source) {
        if (!source.changed) {
            return;
        }
        if (target.originalOwner == InvalidNode) {
            target.originalOwner = source.originalOwner;
        }
        target.originalOwners.insert(target.originalOwners.end(),
                                     source.originalOwners.begin(),
                                     source.originalOwners.end());
        target.createdOutsideNodes.insert(target.createdOutsideNodes.end(),
                                          source.createdOutsideNodes.begin(),
                                          source.createdOutsideNodes.end());
        target.changed = true;
    }

    static RegionSplitBatchApplyResult wrapSingleOwnerResult(RegionSplitApplyResult singleResult) {
        RegionSplitBatchApplyResult batchResult;
        batchResult.reconstruction = singleResult.reconstruction;
        batchResult.componentResults.push_back(std::move(singleResult));
        return batchResult;
    }

    static bool isCompleteDirectRegionOfSingleOwner(const DynamicComponentTree &tree,
                                                    std::span<const PixelId> pixelsC) {
        if (pixelsC.empty()) {
            return false;
        }

        NodeId ownerNodeId = InvalidNode;
        for (PixelId pixelId : pixelsC) {
            const NodeId currentOwner = tree.getSmallestComponent(pixelId);
            if (currentOwner == InvalidNode || !tree.isAlive(currentOwner)) {
                return false;
            }
            if (ownerNodeId == InvalidNode) {
                ownerNodeId = currentOwner;
            } else if (ownerNodeId != currentOwner) {
                return false;
            }
        }

        const auto regionSet = makeRegionSet(tree, pixelsC);
        return nodeProperPartsEqualRegion(tree, ownerNodeId, regionSet);
    }

    static bool canUseFlatLowerSingleOwnerShortcut(const DynamicComponentTree &tree,
                                                   std::span<const PixelId> pixelsC,
                                                   AltitudeType b) {
        if (!isCompleteDirectRegionOfSingleOwner(tree, pixelsC)) {
            return false;
        }

        const NodeId sourceNodeId = findSingleDirectOwner(tree, pixelsC);
        if (tree.isRoot(sourceNodeId)) {
            return true;
        }

        const NodeId parentNodeId = tree.getNodeParent(sourceNodeId);
        return parentNodeId != InvalidNode &&
               tree.isAlive(parentNodeId) &&
               tree.getAltitude(parentNodeId) <= b;
    }

    static RegionSplitBatchApplyResult updateFlatLowerMaterializedDispatch(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionLevelChange<AltitudeType> &change) {
        if (canUseFlatLowerSingleOwnerShortcut(
                maxTree,
                std::span<const PixelId>(change.pixelsC),
                change.b)) {
            return wrapSingleOwnerResult(updateFlatLowerSingleOwner(maxTree, adjacency, change));
        }

        const auto plan = buildFlatLowerBatchPlan(maxTree, adjacency, change);
        if (plan.interaction == RegionSplitBatchInteraction::RepeatedSourceNode) {
            return updateFlatLowerSourceGroupedBatch(maxTree, adjacency, change);
        }
        return updateFlatLowerBatchFromPlan(maxTree, adjacency, plan);
    }

    static RegionSplitBatchApplyResult updateFlatLowerMixedInsertedAndNonInsertedTargets(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        if (plan.interaction != RegionSplitBatchInteraction::MixedInsertedAndNonInsertedTargets) {
            throw std::runtime_error("RegionSplitUpdater mixed target update requires a mixed-target batch plan");
        }

        std::vector<PixelId> nonInsertedPixels;
        std::vector<RegionLevelChange<AltitudeType>> insertedGroupChanges;
        std::vector<NodeId> insertedGroupParents;

        for (const auto &component : plan.components) {
            if (!component.requiresInsertedLevelBAncestor) {
                nonInsertedPixels.insert(nonInsertedPixels.end(),
                                         component.pixelsC.begin(),
                                         component.pixelsC.end());
                continue;
            }

            auto it = std::find(insertedGroupParents.begin(), insertedGroupParents.end(), component.parentNode);
            if (it == insertedGroupParents.end()) {
                insertedGroupParents.push_back(component.parentNode);
                insertedGroupChanges.push_back(buildChangeForPixelsLike(component.pixelsC, plan.change));
            } else {
                const size_t groupIndex = static_cast<size_t>(it - insertedGroupParents.begin());
                insertedGroupChanges[groupIndex].pixelsC.insert(insertedGroupChanges[groupIndex].pixelsC.end(),
                                                                component.pixelsC.begin(),
                                                                component.pixelsC.end());
            }
        }

        if (nonInsertedPixels.empty() || insertedGroupChanges.empty()) {
            throw std::runtime_error("RegionSplitUpdater mixed target update requires both inserted and non-inserted groups");
        }

        const auto nonInsertedChange = buildChangeForPixelsLike(nonInsertedPixels, plan.change);
        const auto nonInsertedPlan = buildFlatLowerBatchPlan(maxTree, adjacency, nonInsertedChange);
        if (!nonInsertedPlan.supportedByCurrentUpdater()) {
            throw std::runtime_error("RegionSplitUpdater mixed target update requires the non-inserted group to be supported by the single-owner batch updater");
        }

        for (const auto &insertedGroupChange : insertedGroupChanges) {
            const auto insertedGroupPlan = buildFlatLowerBatchPlan(maxTree, adjacency, insertedGroupChange);
            if (insertedGroupPlan.interaction != RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor) {
                throw std::runtime_error("RegionSplitUpdater mixed target update currently requires each inserted group to stay connected after removing C");
            }
        }

        RegionSplitBatchApplyResult result;
        auto nonInsertedResult = updateFlatLowerIndependentSingleOwners(
            maxTree,
            adjacency,
            nonInsertedChange);
        result.componentResults.insert(result.componentResults.end(),
                                       nonInsertedResult.componentResults.begin(),
                                       nonInsertedResult.componentResults.end());

        for (const auto &insertedGroupChange : insertedGroupChanges) {
            auto insertedGroupResult = updateFlatLowerSharedInsertedLevelBAncestor(
                maxTree,
                adjacency,
                insertedGroupChange);
            result.componentResults.insert(result.componentResults.end(),
                                           insertedGroupResult.componentResults.begin(),
                                           insertedGroupResult.componentResults.end());
        }

        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static RegionSplitBatchApplyResult updateFlatLowerMultipleInsertedLevelBAncestors(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        if (plan.interaction != RegionSplitBatchInteraction::RequiresMultipleInsertedLevelBAncestors) {
            throw std::runtime_error("RegionSplitUpdater multiple inserted ancestor update requires a multiple-ancestor batch plan");
        }

        std::vector<RegionLevelChange<AltitudeType>> groupChanges;
        std::vector<NodeId> groupParents;
        for (const auto &component : plan.components) {
            if (!component.requiresInsertedLevelBAncestor) {
                throw std::runtime_error("RegionSplitUpdater multiple inserted ancestor update requires every component to need insertion at b");
            }

            auto it = std::find(groupParents.begin(), groupParents.end(), component.parentNode);
            if (it == groupParents.end()) {
                groupParents.push_back(component.parentNode);
                groupChanges.push_back(buildChangeForPixelsLike(component.pixelsC, plan.change));
            } else {
                const size_t groupIndex = static_cast<size_t>(it - groupParents.begin());
                groupChanges[groupIndex].pixelsC.insert(groupChanges[groupIndex].pixelsC.end(),
                                                        component.pixelsC.begin(),
                                                        component.pixelsC.end());
            }
        }

        if (groupChanges.size() < 2) {
            throw std::runtime_error("RegionSplitUpdater multiple inserted ancestor update requires at least two parent groups");
        }

        for (const auto &groupChange : groupChanges) {
            const auto groupPlan = buildFlatLowerBatchPlan(maxTree, adjacency, groupChange);
            if (groupPlan.interaction != RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor) {
                throw std::runtime_error("RegionSplitUpdater multiple inserted ancestor update currently requires each parent group to stay connected after removing C");
            }
        }

        RegionSplitBatchApplyResult result;
        for (const auto &groupChange : groupChanges) {
            auto groupResult = updateFlatLowerSharedInsertedLevelBAncestor(
                maxTree,
                adjacency,
                groupChange);
            result.componentResults.insert(result.componentResults.end(),
                                           groupResult.componentResults.begin(),
                                           groupResult.componentResults.end());
        }
        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static RegionSplitBatchApplyResult updateFlatLowerPreservedAncestorSplit(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        if (plan.interaction != RegionSplitBatchInteraction::PreservedAncestorSplits) {
            throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires a split-ancestor batch plan");
        }
        if (plan.components.empty()) {
            throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires at least one component");
        }

        const NodeId preservedAncestorId = plan.components.front().parentNode;
        if (preservedAncestorId == InvalidNode || !maxTree.isAlive(preservedAncestorId)) {
            throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires a live preserved ancestor");
        }

        std::unordered_set<NodeId> sourceNodes;
        sourceNodes.reserve(plan.components.size());
        for (const auto &component : plan.components) {
            if (component.parentNode != preservedAncestorId) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update currently requires all sources to share the same preserved ancestor");
            }
            if (!component.requiresInsertedLevelBAncestor) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires every source parent to be above b");
            }
            if (maxTree.getNodeParent(component.sourceNode) != preservedAncestorId) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires sources to remain direct children of the preserved ancestor");
            }
            sourceNodes.insert(component.sourceNode);
        }

        const auto regionSet = makeRegionSet(maxTree, plan.change.pixelsC);
        auto remainingAncestorComponents = connectedComponentsAfterRemovingRegion(
            adjacency,
            maxTree.getPixelsOfCC(preservedAncestorId),
            regionSet);
        if (remainingAncestorComponents.size() < 2) {
            throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires at least two remaining ancestor components");
        }

        std::vector<int> componentIndexByPixel(
            static_cast<size_t>(maxTree.getNumTotalProperParts()),
            -1);
        for (size_t componentIndex = 0; componentIndex < remainingAncestorComponents.size(); ++componentIndex) {
            for (PixelId pixelId : remainingAncestorComponents[componentIndex]) {
                componentIndexByPixel[static_cast<size_t>(pixelId)] = static_cast<int>(componentIndex);
            }
        }

        std::vector<PixelId> preservedAncestorProperParts;
        preservedAncestorProperParts.reserve(static_cast<size_t>(maxTree.getNumProperParts(preservedAncestorId)));
        std::vector<int> numDirectProperPartsByComponent(remainingAncestorComponents.size(), 0);
        for (PixelId pixelId : maxTree.getProperParts(preservedAncestorId)) {
            const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
            if (componentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update found a preserved proper part outside X minus C");
            }
            preservedAncestorProperParts.push_back(pixelId);
            ++numDirectProperPartsByComponent[static_cast<size_t>(componentIndex)];
        }
        auto locateRemainingComponent = [&](NodeId nodeId) {
            int selectedComponentIndex = -1;
            for (PixelId pixelId : maxTree.getPixelsOfCC(nodeId)) {
                if (regionSet.find(pixelId) != regionSet.end()) {
                    throw std::runtime_error("RegionSplitUpdater preserved ancestor split update cannot reattach a child that intersects C");
                }
                const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
                if (componentIndex < 0) {
                    throw std::runtime_error("RegionSplitUpdater preserved ancestor split update found a child outside X minus C");
                }
                if (selectedComponentIndex == -1) {
                    selectedComponentIndex = componentIndex;
                } else if (selectedComponentIndex != componentIndex) {
                    throw std::runtime_error("RegionSplitUpdater preserved ancestor split update cannot reattach a child spanning multiple remaining components");
                }
            }
            if (selectedComponentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update cannot reattach an empty child");
            }
            return selectedComponentIndex;
        };

        std::vector<std::pair<NodeId, int>> childMoves;
        for (NodeId childId : maxTree.getChildren(preservedAncestorId)) {
            if (sourceNodes.find(childId) == sourceNodes.end()) {
                childMoves.push_back({childId, locateRemainingComponent(childId)});
            }
        }
        for (const auto &component : plan.components) {
            for (NodeId childId : maxTree.getChildren(component.sourceNode)) {
                childMoves.push_back({childId, locateRemainingComponent(childId)});
            }
        }

        std::vector<int> numMovedChildrenByComponent(remainingAncestorComponents.size(), 0);
        for (const auto &childMove : childMoves) {
            ++numMovedChildrenByComponent[static_cast<size_t>(childMove.second)];
        }
        for (size_t componentIndex = 0; componentIndex < remainingAncestorComponents.size(); ++componentIndex) {
            if (numDirectProperPartsByComponent[componentIndex] == 0 &&
                numMovedChildrenByComponent[componentIndex] != 1) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires components without direct proper parts to contain exactly one preserved child");
            }
        }

        const bool preservedAncestorWasRoot = maxTree.isRoot(preservedAncestorId);
        const NodeId externalParentId = preservedAncestorWasRoot
            ? InvalidNode
            : maxTree.getNodeParent(preservedAncestorId);
        bool reusesExistingLevelBTarget = false;
        if (!preservedAncestorWasRoot) {
            if (externalParentId == InvalidNode || !maxTree.isAlive(externalParentId)) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update requires a live external parent");
            }
            if (maxTree.getAltitude(externalParentId) == plan.change.b) {
                reusesExistingLevelBTarget = true;
            } else if (maxTree.getAltitude(externalParentId) > plan.change.b) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update currently requires the external parent altitude to be below b");
            }
        }

        NodeId targetNodeId = externalParentId;
        bool createdTargetNode = false;
        if (!reusesExistingLevelBTarget) {
            targetNodeId = maxTree.createDetachedNode(static_cast<int>(plan.change.b));
            if (targetNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update could not create target node");
            }
            createdTargetNode = true;
        }

        std::vector<NodeId> componentTargetNodes(remainingAncestorComponents.size(), InvalidNode);
        for (size_t i = 0; i < remainingAncestorComponents.size(); ++i) {
            if (numDirectProperPartsByComponent[i] == 0) {
                continue;
            }
            const NodeId splitNodeId = maxTree.createDetachedNode(maxTree.getAltitude(preservedAncestorId));
            if (splitNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater preserved ancestor split update could not create a split node");
            }
            componentTargetNodes[i] = splitNodeId;
        }

        if (!reusesExistingLevelBTarget) {
            if (preservedAncestorWasRoot) {
                maxTree.setRoot(targetNodeId);
            } else {
                maxTree.removeChild(externalParentId, preservedAncestorId, false);
                maxTree.attachNode(externalParentId, targetNodeId);
            }
        }

        for (NodeId splitNodeId : componentTargetNodes) {
            if (splitNodeId != InvalidNode) {
                maxTree.attachNode(targetNodeId, splitNodeId);
            }
        }

        for (size_t componentIndex = 0; componentIndex < remainingAncestorComponents.size(); ++componentIndex) {
            const NodeId splitNodeId = componentTargetNodes[componentIndex];
            if (splitNodeId == InvalidNode) {
                continue;
            }
            for (PixelId pixelId : preservedAncestorProperParts) {
                if (componentIndexByPixel[static_cast<size_t>(pixelId)] != static_cast<int>(componentIndex)) {
                    continue;
                }
                maxTree.moveProperPart(splitNodeId, preservedAncestorId, pixelId);
            }
        }

        for (const auto &childMove : childMoves) {
            NodeId targetForChild = componentTargetNodes[static_cast<size_t>(childMove.second)];
            if (targetForChild == InvalidNode) {
                targetForChild = targetNodeId;
            }
            maxTree.moveNode(childMove.first, targetForChild);
        }

        RegionSplitBatchApplyResult result;
        result.componentResults.reserve(plan.components.size());
        bool targetCreationRecorded = false;
        for (const auto &component : plan.components) {
            RegionSplitApplyResult componentResult;
            componentResult.sourceNode = component.sourceNode;
            componentResult.sourceNodes.push_back(component.sourceNode);
            componentResult.targetNode = targetNodeId;
            componentResult.createdTargetNode = createdTargetNode && !targetCreationRecorded;

            moveRegionProperParts(
                maxTree,
                targetNodeId,
                component.sourceNode,
                std::span<const PixelId>(component.pixelsC));
            maxTree.removeChild(preservedAncestorId, component.sourceNode, true);

            componentResult.releasedSourceNode = true;
            result.componentResults.push_back(std::move(componentResult));
            targetCreationRecorded = true;
        }

        if (reusesExistingLevelBTarget) {
            maxTree.removeChild(targetNodeId, preservedAncestorId, true);
        } else {
            maxTree.releaseNode(preservedAncestorId);
        }

        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static RegionSplitBatchApplyResult updateFlatLowerInteractingInsertedLevelBAncestors(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        if (plan.interaction != RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors) {
            throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires an interacting-ancestor batch plan");
        }
        if (plan.components.empty()) {
            throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires at least one component");
        }

        std::vector<NodeId> insertedParentIds;
        insertedParentIds.reserve(plan.components.size());
        std::unordered_set<NodeId> sourceNodes;
        sourceNodes.reserve(plan.components.size());
        for (const auto &component : plan.components) {
            if (!component.requiresInsertedLevelBAncestor) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires every source to need insertion at b");
            }
            if (component.parentNode == InvalidNode ||
                !maxTree.isAlive(component.parentNode)) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires live preserved parents");
            }
            if (maxTree.getNodeParent(component.sourceNode) != component.parentNode) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires sources to remain direct children of their preserved parents");
            }
            if (std::find(insertedParentIds.begin(), insertedParentIds.end(), component.parentNode) == insertedParentIds.end()) {
                insertedParentIds.push_back(component.parentNode);
            }
            sourceNodes.insert(component.sourceNode);
        }

        NodeId outerAncestorId = InvalidNode;
        for (NodeId candidateId : insertedParentIds) {
            bool containsAllParents = true;
            for (NodeId parentId : insertedParentIds) {
                if (candidateId == parentId) {
                    continue;
                }
                if (!isStrictAncestor(maxTree, candidateId, parentId)) {
                    containsAllParents = false;
                    break;
                }
            }
            if (containsAllParents) {
                if (outerAncestorId != InvalidNode) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires a unique outer preserved ancestor");
                }
                outerAncestorId = candidateId;
            }
        }
        if (outerAncestorId == InvalidNode) {
            throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires one preserved parent to contain all inserted parent groups");
        }

        const bool outerWasRoot = maxTree.isRoot(outerAncestorId);
        const NodeId externalParentId = outerWasRoot
            ? InvalidNode
            : maxTree.getNodeParent(outerAncestorId);
        bool reusesExistingLevelBTarget = false;
        if (!outerWasRoot) {
            if (externalParentId == InvalidNode || !maxTree.isAlive(externalParentId)) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires a live external parent");
            }
            if (maxTree.getAltitude(externalParentId) == plan.change.b) {
                reusesExistingLevelBTarget = true;
            } else if (maxTree.getAltitude(externalParentId) > plan.change.b) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update currently requires the external parent altitude to be below b");
            }
        }

        const auto regionSet = makeRegionSet(maxTree, plan.change.pixelsC);
        auto remainingOuterComponents = connectedComponentsAfterRemovingRegion(
            adjacency,
            maxTree.getPixelsOfCC(outerAncestorId),
            regionSet);
        if (remainingOuterComponents.empty()) {
            throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires a non-empty outer complement after removing C");
        }
        if (remainingOuterComponents.size() == 1) {
            NodeId targetNodeId = externalParentId;
            bool createdTargetNode = false;
            if (!reusesExistingLevelBTarget) {
                targetNodeId = maxTree.createDetachedNode(static_cast<int>(plan.change.b));
                if (targetNodeId == InvalidNode) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update could not create target node");
                }
                createdTargetNode = true;

                if (outerWasRoot) {
                    maxTree.setRoot(targetNodeId);
                    maxTree.attachNode(targetNodeId, outerAncestorId);
                } else {
                    maxTree.removeChild(externalParentId, outerAncestorId, false);
                    maxTree.attachNode(externalParentId, targetNodeId);
                    maxTree.attachNode(targetNodeId, outerAncestorId);
                }
            }

            RegionSplitBatchApplyResult result;
            result.componentResults.reserve(plan.components.size());
            bool targetCreationRecorded = false;
            for (const auto &component : plan.components) {
                RegionSplitApplyResult componentResult;
                componentResult.sourceNode = component.sourceNode;
                componentResult.sourceNodes.push_back(component.sourceNode);
                componentResult.targetNode = targetNodeId;
                componentResult.createdTargetNode = createdTargetNode && !targetCreationRecorded;

                moveRegionProperParts(
                    maxTree,
                    targetNodeId,
                    component.sourceNode,
                    std::span<const PixelId>(component.pixelsC));
                maxTree.moveChildren(component.parentNode, component.sourceNode);
                maxTree.removeChild(component.parentNode, component.sourceNode, true);

                componentResult.releasedSourceNode = true;
                result.componentResults.push_back(std::move(componentResult));
                targetCreationRecorded = true;
            }

            result.reconstruction = maxTree.reconstructionImage();
            for (auto &componentResult : result.componentResults) {
                componentResult.reconstruction = result.reconstruction;
            }
            return result;
        }

        std::vector<int> componentIndexByPixel(
            static_cast<size_t>(maxTree.getNumTotalProperParts()),
            -1);
        for (size_t componentIndex = 0; componentIndex < remainingOuterComponents.size(); ++componentIndex) {
            for (PixelId pixelId : remainingOuterComponents[componentIndex]) {
                componentIndexByPixel[static_cast<size_t>(pixelId)] = static_cast<int>(componentIndex);
            }
        }

        std::vector<PixelId> outerProperParts;
        outerProperParts.reserve(static_cast<size_t>(maxTree.getNumProperParts(outerAncestorId)));
        std::vector<int> numOuterProperPartsByComponent(remainingOuterComponents.size(), 0);
        for (PixelId pixelId : maxTree.getProperParts(outerAncestorId)) {
            if (regionSet.find(pixelId) != regionSet.end()) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update does not support C as direct proper parts of the outer ancestor");
            }
            const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
            if (componentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update found an outer proper part outside X minus C");
            }
            outerProperParts.push_back(pixelId);
            ++numOuterProperPartsByComponent[static_cast<size_t>(componentIndex)];
        }

        auto locateRemainingComponentOrSpanning = [&](NodeId nodeId) {
            int selectedComponentIndex = -1;
            for (PixelId pixelId : maxTree.getPixelsOfCC(nodeId)) {
                if (regionSet.find(pixelId) != regionSet.end()) {
                    continue;
                }
                const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
                if (componentIndex < 0) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update found a child outside X minus C");
                }
                if (selectedComponentIndex == -1) {
                    selectedComponentIndex = componentIndex;
                } else if (selectedComponentIndex != componentIndex) {
                    return -2;
                }
            }
            if (selectedComponentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update cannot reattach a child with empty support after removing C");
            }
            return selectedComponentIndex;
        };
        auto locateRemainingComponent = [&](NodeId nodeId) {
            const int componentIndex = locateRemainingComponentOrSpanning(nodeId);
            if (componentIndex < 0) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update cannot reattach a child spanning multiple remaining components");
            }
            return componentIndex;
        };

        struct InteractingChildMove {
            NodeId childNode = InvalidNode;
            int componentIndex = -1;
            int parentSplitIndex = -1;
        };

        struct InteractingSplitChildPlan {
            NodeId childNode = InvalidNode;
            int parentSplitIndex = -1;
            std::vector<std::vector<PixelId>> properPartsByComponent;
            std::vector<NodeId> splitNodesByComponent;
            std::vector<unsigned char> subtreeTouchesComponent;
        };

        std::vector<InteractingChildMove> childMoves;
        std::vector<InteractingSplitChildPlan> splitChildPlans;
        std::vector<int> numSplitChildrenByComponent(remainingOuterComponents.size(), 0);
        std::unordered_set<NodeId> insertedParentSet(insertedParentIds.begin(), insertedParentIds.end());
        std::unordered_map<NodeId, int> splitPlanIndexByNode;
        splitPlanIndexByNode.reserve(insertedParentIds.size());

        auto collectSplitChildPlan = [&](auto &&self,
                                         NodeId childId,
                                         int parentSplitIndex) -> void {
            const int componentIndex = locateRemainingComponentOrSpanning(childId);
            if (componentIndex >= 0) {
                childMoves.push_back({childId, componentIndex, parentSplitIndex});
                return;
            }

            if (insertedParentSet.find(childId) == insertedParentSet.end()) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update can split only inserted-parent children that span multiple remaining components");
            }

            InteractingSplitChildPlan splitChildPlan;
            splitChildPlan.childNode = childId;
            splitChildPlan.parentSplitIndex = parentSplitIndex;
            splitChildPlan.properPartsByComponent.resize(remainingOuterComponents.size());
            splitChildPlan.splitNodesByComponent.assign(remainingOuterComponents.size(), InvalidNode);
            splitChildPlan.subtreeTouchesComponent.assign(remainingOuterComponents.size(), 0);

            for (PixelId pixelId : maxTree.getPixelsOfCC(childId)) {
                if (regionSet.find(pixelId) != regionSet.end()) {
                    continue;
                }
                const int touchedComponentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
                if (touchedComponentIndex < 0) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update found an internal subtree pixel outside X minus C");
                }
                splitChildPlan.subtreeTouchesComponent[static_cast<size_t>(touchedComponentIndex)] = 1;
            }

            for (PixelId pixelId : maxTree.getProperParts(childId)) {
                if (regionSet.find(pixelId) != regionSet.end()) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update does not support C as direct proper parts of an internal preserved parent");
                }
                const int properPartComponentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
                if (properPartComponentIndex < 0) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update found an internal proper part outside X minus C");
                }
                splitChildPlan
                    .properPartsByComponent[static_cast<size_t>(properPartComponentIndex)]
                    .push_back(pixelId);
            }

            const int splitPlanIndex = static_cast<int>(splitChildPlans.size());
            splitPlanIndexByNode.emplace(childId, splitPlanIndex);
            splitChildPlans.push_back(std::move(splitChildPlan));

            for (NodeId grandChildId : maxTree.getChildren(childId)) {
                if (sourceNodes.find(grandChildId) != sourceNodes.end()) {
                    continue;
                }
                self(self, grandChildId, splitPlanIndex);
            }

            if (parentSplitIndex == -1) {
                for (size_t splitComponentIndex = 0;
                     splitComponentIndex < splitChildPlans[static_cast<size_t>(splitPlanIndex)].subtreeTouchesComponent.size();
                     ++splitComponentIndex) {
                    if (splitChildPlans[static_cast<size_t>(splitPlanIndex)].subtreeTouchesComponent[splitComponentIndex] != 0) {
                        ++numSplitChildrenByComponent[splitComponentIndex];
                    }
                }
            }
        };

        for (NodeId childId : maxTree.getChildren(outerAncestorId)) {
            if (sourceNodes.find(childId) != sourceNodes.end()) {
                continue;
            }
            collectSplitChildPlan(collectSplitChildPlan, childId, -1);
        }
        for (const auto &component : plan.components) {
            int parentSplitIndex = -1;
            auto splitPlanIt = splitPlanIndexByNode.find(component.parentNode);
            if (splitPlanIt != splitPlanIndexByNode.end()) {
                parentSplitIndex = splitPlanIt->second;
            }
            for (NodeId childId : maxTree.getChildren(component.sourceNode)) {
                childMoves.push_back({childId, locateRemainingComponent(childId), parentSplitIndex});
            }
        }

        std::vector<int> numMovedChildrenByComponent(remainingOuterComponents.size(), 0);
        for (const auto &childMove : childMoves) {
            if (childMove.parentSplitIndex == -1) {
                ++numMovedChildrenByComponent[static_cast<size_t>(childMove.componentIndex)];
            }
        }
        for (size_t componentIndex = 0; componentIndex < remainingOuterComponents.size(); ++componentIndex) {
            if (numOuterProperPartsByComponent[componentIndex] == 0 &&
                numMovedChildrenByComponent[componentIndex] + numSplitChildrenByComponent[componentIndex] != 1) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update requires components without direct outer proper parts to contain exactly one preserved child");
            }
        }

        NodeId targetNodeId = externalParentId;
        bool createdTargetNode = false;
        if (!reusesExistingLevelBTarget) {
            targetNodeId = maxTree.createDetachedNode(static_cast<int>(plan.change.b));
            if (targetNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update could not create target node");
            }
            createdTargetNode = true;
        }

        std::vector<NodeId> componentTargetNodes(remainingOuterComponents.size(), InvalidNode);
        for (size_t componentIndex = 0; componentIndex < remainingOuterComponents.size(); ++componentIndex) {
            if (numOuterProperPartsByComponent[componentIndex] == 0) {
                continue;
            }
            const NodeId splitNodeId = maxTree.createDetachedNode(maxTree.getAltitude(outerAncestorId));
            if (splitNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update could not create an outer split node");
            }
            componentTargetNodes[componentIndex] = splitNodeId;
        }

        if (!reusesExistingLevelBTarget) {
            if (outerWasRoot) {
                maxTree.setRoot(targetNodeId);
            } else {
                maxTree.removeChild(externalParentId, outerAncestorId, false);
                maxTree.attachNode(externalParentId, targetNodeId);
            }
        }

        for (NodeId splitNodeId : componentTargetNodes) {
            if (splitNodeId != InvalidNode) {
                maxTree.attachNode(targetNodeId, splitNodeId);
            }
        }

        for (PixelId pixelId : outerProperParts) {
            const int componentIndex = componentIndexByPixel[static_cast<size_t>(pixelId)];
            const NodeId splitNodeId = componentTargetNodes[static_cast<size_t>(componentIndex)];
            if (splitNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update cannot move an outer proper part without an outer split node");
            }
            maxTree.moveProperPart(splitNodeId, outerAncestorId, pixelId);
        }

        auto targetForComponentAndSplitParent = [&](size_t componentIndex,
                                                    int parentSplitIndex) {
            int currentSplitIndex = parentSplitIndex;
            while (currentSplitIndex >= 0) {
                const auto &parentSplitPlan = splitChildPlans[static_cast<size_t>(currentSplitIndex)];
                const NodeId splitNodeId = parentSplitPlan.splitNodesByComponent[componentIndex];
                if (splitNodeId != InvalidNode) {
                    return splitNodeId;
                }
                currentSplitIndex = parentSplitPlan.parentSplitIndex;
            }

            NodeId targetForComponent = componentTargetNodes[componentIndex];
            if (targetForComponent == InvalidNode) {
                targetForComponent = targetNodeId;
            }
            return targetForComponent;
        };

        for (size_t splitPlanIndex = 0; splitPlanIndex < splitChildPlans.size(); ++splitPlanIndex) {
            const auto &splitChildPlan = splitChildPlans[splitPlanIndex];
            for (size_t componentIndex = 0;
                 componentIndex < splitChildPlan.properPartsByComponent.size();
                 ++componentIndex) {
                const auto &properParts = splitChildPlan.properPartsByComponent[componentIndex];
                if (properParts.empty()) {
                    continue;
                }

                const NodeId childSplitNodeId = maxTree.createDetachedNode(maxTree.getAltitude(splitChildPlan.childNode));
                if (childSplitNodeId == InvalidNode) {
                    throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update could not create an internal split node");
                }

                splitChildPlans[splitPlanIndex].splitNodesByComponent[componentIndex] = childSplitNodeId;
                const NodeId parentForChildSplit = targetForComponentAndSplitParent(
                    componentIndex,
                    splitChildPlan.parentSplitIndex);
                maxTree.attachNode(parentForChildSplit, childSplitNodeId);

                for (PixelId pixelId : properParts) {
                    maxTree.moveProperPart(childSplitNodeId, splitChildPlan.childNode, pixelId);
                }
            }
        }

        for (const auto &childMove : childMoves) {
            const NodeId targetForChild = targetForComponentAndSplitParent(
                static_cast<size_t>(childMove.componentIndex),
                childMove.parentSplitIndex);
            maxTree.moveNode(childMove.childNode, targetForChild);
        }

        RegionSplitBatchApplyResult result;
        result.componentResults.reserve(plan.components.size());
        bool targetCreationRecorded = false;
        for (const auto &component : plan.components) {
            RegionSplitApplyResult componentResult;
            componentResult.sourceNode = component.sourceNode;
            componentResult.sourceNodes.push_back(component.sourceNode);
            componentResult.targetNode = targetNodeId;
            componentResult.createdTargetNode = createdTargetNode && !targetCreationRecorded;

            moveRegionProperParts(
                maxTree,
                targetNodeId,
                component.sourceNode,
                std::span<const PixelId>(component.pixelsC));
            maxTree.removeChild(component.parentNode, component.sourceNode, true);

            componentResult.releasedSourceNode = true;
            result.componentResults.push_back(std::move(componentResult));
            targetCreationRecorded = true;
        }

        for (size_t reverseIndex = splitChildPlans.size(); reverseIndex > 0; --reverseIndex) {
            const auto &splitChildPlan = splitChildPlans[reverseIndex - 1];
            if (maxTree.getNumProperParts(splitChildPlan.childNode) != 0 ||
                maxTree.getNumChildren(splitChildPlan.childNode) != 0) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update could not empty an internal split child");
            }
            const NodeId splitChildParentId = maxTree.getNodeParent(splitChildPlan.childNode);
            if (splitChildParentId == InvalidNode ||
                !maxTree.isAlive(splitChildParentId) ||
                splitChildParentId == splitChildPlan.childNode) {
                throw std::runtime_error("RegionSplitUpdater interacting inserted ancestor update cannot locate an internal split child's parent");
            }
            maxTree.removeChild(splitChildParentId, splitChildPlan.childNode, true);
        }

        if (reusesExistingLevelBTarget) {
            maxTree.removeChild(targetNodeId, outerAncestorId, true);
        } else {
            maxTree.releaseNode(outerAncestorId);
        }

        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static RegionSplitBatchApplyResult updateFlatLowerSharedInsertedLevelBAncestor(
        DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan) {
        if (plan.interaction != RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor) {
            throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires an inserted-ancestor batch plan");
        }
        if (plan.components.empty()) {
            throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires at least one component");
        }

        const NodeId preservedAncestorId = plan.components.front().parentNode;
        if (preservedAncestorId == InvalidNode || !maxTree.isAlive(preservedAncestorId)) {
            throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires a live preserved ancestor");
        }

        for (const auto &component : plan.components) {
            if (component.parentNode != preservedAncestorId) {
                throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update currently requires all sources to share the same preserved ancestor");
            }
            if (!component.requiresInsertedLevelBAncestor) {
                throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires every source parent to be above b");
            }
            if (maxTree.getNodeParent(component.sourceNode) != preservedAncestorId) {
                throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires sources to remain direct children of the preserved ancestor");
            }
        }

        const auto regionSet = makeRegionSet(maxTree, plan.change.pixelsC);
        auto remainingAncestorComponents = connectedComponentsAfterRemovingRegion(
            adjacency,
            maxTree.getPixelsOfCC(preservedAncestorId),
            regionSet);
        if (remainingAncestorComponents.size() != 1) {
            throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires the preserved ancestor to stay connected after removing C");
        }

        const bool preservedAncestorWasRoot = maxTree.isRoot(preservedAncestorId);
        const NodeId externalParentId = preservedAncestorWasRoot
            ? InvalidNode
            : maxTree.getNodeParent(preservedAncestorId);
        bool reusesExistingLevelBTarget = false;
        if (!preservedAncestorWasRoot) {
            if (externalParentId == InvalidNode || !maxTree.isAlive(externalParentId)) {
                throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update requires a live external parent");
            }
            if (maxTree.getAltitude(externalParentId) == plan.change.b) {
                reusesExistingLevelBTarget = true;
            } else if (maxTree.getAltitude(externalParentId) > plan.change.b) {
                throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update currently requires the external parent altitude to be below b");
            }
        }

        NodeId targetNodeId = externalParentId;
        bool createdTargetNode = false;
        if (!reusesExistingLevelBTarget) {
            targetNodeId = maxTree.createDetachedNode(static_cast<int>(plan.change.b));
            if (targetNodeId == InvalidNode) {
                throw std::runtime_error("RegionSplitUpdater shared inserted ancestor update could not create target node");
            }
            createdTargetNode = true;

            if (preservedAncestorWasRoot) {
                maxTree.setRoot(targetNodeId);
                maxTree.attachNode(targetNodeId, preservedAncestorId);
            } else {
                maxTree.removeChild(externalParentId, preservedAncestorId, false);
                maxTree.attachNode(externalParentId, targetNodeId);
                maxTree.attachNode(targetNodeId, preservedAncestorId);
            }
        }

        RegionSplitBatchApplyResult result;
        result.componentResults.reserve(plan.components.size());
        bool targetCreationRecorded = false;
        for (const auto &component : plan.components) {
            RegionSplitApplyResult componentResult;
            componentResult.sourceNode = component.sourceNode;
            componentResult.sourceNodes.push_back(component.sourceNode);
            componentResult.targetNode = targetNodeId;
            componentResult.createdTargetNode = createdTargetNode && !targetCreationRecorded;

            moveRegionProperParts(
                maxTree,
                targetNodeId,
                component.sourceNode,
                std::span<const PixelId>(component.pixelsC));
            maxTree.moveChildren(preservedAncestorId, component.sourceNode);
            maxTree.removeChild(preservedAncestorId, component.sourceNode, true);

            componentResult.releasedSourceNode = true;
            result.componentResults.push_back(std::move(componentResult));
            targetCreationRecorded = true;
        }

        result.reconstruction = maxTree.reconstructionImage();
        for (auto &componentResult : result.componentResults) {
            componentResult.reconstruction = result.reconstruction;
        }
        return result;
    }

    static const char *unsupportedBatchInteractionMessage(
        RegionSplitBatchInteraction interaction) {
        switch (interaction) {
        case RegionSplitBatchInteraction::RepeatedSourceNode:
            return "RegionSplitUpdater batch update does not yet support multiple C components owned by the same source node";
        case RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor:
            return "RegionSplitUpdater batch update does not yet support inserting a shared level-b ancestor below a higher preserved component";
        case RegionSplitBatchInteraction::MixedInsertedAndNonInsertedTargets:
            return "RegionSplitUpdater batch update does not yet support mixing inserted level-b ancestors with independent or existing level-b targets";
        case RegionSplitBatchInteraction::RequiresMultipleInsertedLevelBAncestors:
            return "RegionSplitUpdater batch update does not yet support inserting level-b ancestors under multiple preserved components";
        case RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors:
            return "RegionSplitUpdater batch update does not yet support interacting inserted level-b ancestor groups";
        case RegionSplitBatchInteraction::PreservedAncestorSplits:
            return "RegionSplitUpdater batch update does not yet support the general preserved-ancestor split case";
        case RegionSplitBatchInteraction::IndependentSingleOwners:
        case RegionSplitBatchInteraction::SharedExistingLevelBTarget:
            break;
        }
        return "RegionSplitUpdater batch update encountered an unsupported interaction";
    }

    static RegionSplitBatchInteraction classifyInsertedLevelBAncestorInteraction(
        const DynamicComponentTree &maxTree,
        AdjacencyRelation &adjacency,
        const RegionSplitBatchPlan &plan,
        const std::unordered_set<PixelId> &regionSet) {
        NodeId insertedParentId = InvalidNode;
        bool hasInserted = false;
        bool hasNonInserted = false;
        std::vector<NodeId> insertedParentIds;

        for (const auto &component : plan.components) {
            if (!component.requiresInsertedLevelBAncestor) {
                hasNonInserted = true;
                continue;
            }
            hasInserted = true;
            if (std::find(insertedParentIds.begin(), insertedParentIds.end(), component.parentNode) == insertedParentIds.end()) {
                insertedParentIds.push_back(component.parentNode);
            }
            if (insertedParentId == InvalidNode) {
                insertedParentId = component.parentNode;
            }
        }

        if (!hasInserted) {
            return RegionSplitBatchInteraction::IndependentSingleOwners;
        }
        if (hasNonInserted) {
            return RegionSplitBatchInteraction::MixedInsertedAndNonInsertedTargets;
        }
        if (insertedParentIds.size() > 1) {
            for (size_t i = 0; i < insertedParentIds.size(); ++i) {
                for (size_t j = i + 1; j < insertedParentIds.size(); ++j) {
                    if (isStrictAncestor(maxTree, insertedParentIds[i], insertedParentIds[j]) ||
                        isStrictAncestor(maxTree, insertedParentIds[j], insertedParentIds[i])) {
                        return RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors;
                    }
                }
            }
            return RegionSplitBatchInteraction::RequiresMultipleInsertedLevelBAncestors;
        }
        if (insertedParentId == InvalidNode || !maxTree.isAlive(insertedParentId)) {
            throw std::runtime_error("RegionSplitUpdater batch plan requires a live inserted-ancestor parent");
        }

        auto remainingAncestorComponents = connectedComponentsAfterRemovingRegion(
            adjacency,
            maxTree.getPixelsOfCC(insertedParentId),
            regionSet);
        if (remainingAncestorComponents.size() != 1) {
            return RegionSplitBatchInteraction::PreservedAncestorSplits;
        }
        return RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor;
    }

    static bool isStrictAncestor(const DynamicComponentTree &tree,
                                 NodeId ancestorId,
                                 NodeId descendantId) {
        if (ancestorId == InvalidNode || descendantId == InvalidNode ||
            !tree.isAlive(ancestorId) || !tree.isAlive(descendantId) ||
            ancestorId == descendantId) {
            return false;
        }

        NodeId current = descendantId;
        while (!tree.isRoot(current)) {
            current = tree.getNodeParent(current);
            if (current == ancestorId) {
                return true;
            }
            if (current == InvalidNode || !tree.isAlive(current)) {
                return false;
            }
        }
        return false;
    }

    static void validateFlatLowerPossiblyDisconnectedAssumptions(
        const DynamicComponentTree &maxTree,
        const RegionLevelChange<AltitudeType> &change) {
        if (!maxTree.isMaxtree()) {
            throw std::runtime_error("RegionSplitUpdater requires a max-tree");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("RegionSplitUpdater requires a non-empty region C");
        }
        if (change.direction != RegionUpdateDirection::Lower || !(change.b < change.a)) {
            throw std::runtime_error("RegionSplitUpdater currently supports only lower updates b < a");
        }

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        if (regionSet.size() != change.pixelsC.size()) {
            throw std::runtime_error("RegionSplitUpdater requires C without duplicate pixels");
        }

        const auto reconstructed = maxTree.reconstructionImage();
        for (PixelId pixelId : change.pixelsC) {
            if ((*reconstructed)[pixelId] != change.a) {
                throw std::runtime_error("RegionSplitUpdater flat plan requires every pixel in C to have value a");
            }
        }
    }

    static NodeId validateFlatLowerSingleOwnerEditable(const DynamicComponentTree &maxTree,
                                                       AdjacencyRelation &adjacency,
                                                       const RegionLevelChange<AltitudeType> &change) {
        (void) buildFlatLowerPlan(maxTree, adjacency, change);

        const auto regionSet = makeRegionSet(maxTree, change.pixelsC);
        const NodeId sourceNodeId = findSingleDirectOwner(maxTree, change.pixelsC);
        if (maxTree.getAltitude(sourceNodeId) != change.a) {
            throw std::runtime_error("RegionSplitUpdater single-owner update requires the owner altitude to be a");
        }
        if (!nodeProperPartsEqualRegion(maxTree, sourceNodeId, regionSet)) {
            throw std::runtime_error("RegionSplitUpdater single-owner update requires C to be the complete direct proper part of one node");
        }

        const bool sourceWasRoot = maxTree.isRoot(sourceNodeId);
        const NodeId oldParentId = maxTree.getNodeParent(sourceNodeId);
        if (!sourceWasRoot &&
            (oldParentId == InvalidNode || !maxTree.isAlive(oldParentId) || maxTree.getAltitude(oldParentId) > change.b)) {
            throw std::runtime_error("RegionSplitUpdater single-owner update cannot attach b below a higher parent altitude");
        }
        return sourceNodeId;
    }

    static NodeId findSingleDirectOwner(const DynamicComponentTree &tree,
                                        std::span<const PixelId> pixelsC) {
        NodeId ownerNodeId = InvalidNode;
        for (PixelId pixelId : pixelsC) {
            const NodeId currentOwner = tree.getSmallestComponent(pixelId);
            if (currentOwner == InvalidNode || !tree.isAlive(currentOwner)) {
                throw std::runtime_error("RegionSplitUpdater cannot locate a live owner for every pixel in C");
            }
            if (ownerNodeId == InvalidNode) {
                ownerNodeId = currentOwner;
            } else if (ownerNodeId != currentOwner) {
                throw std::runtime_error("RegionSplitUpdater single-owner update requires all pixels in C to have the same direct owner");
            }
        }
        if (ownerNodeId == InvalidNode) {
            throw std::runtime_error("RegionSplitUpdater single-owner update requires a non-empty region C");
        }
        return ownerNodeId;
    }

    static bool nodeProperPartsEqualRegion(const DynamicComponentTree &tree,
                                           NodeId nodeId,
                                           const std::unordered_set<PixelId> &regionSet) {
        if (tree.getNumProperParts(nodeId) != static_cast<int>(regionSet.size())) {
            return false;
        }
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            if (regionSet.find(pixelId) == regionSet.end()) {
                return false;
            }
        }
        return true;
    }

    static void moveRegionProperParts(DynamicComponentTree &tree,
                                      NodeId targetNodeId,
                                      NodeId sourceNodeId,
                                      std::span<const PixelId> pixelsC) {
        for (PixelId pixelId : pixelsC) {
            if (tree.getSmallestComponent(pixelId) != sourceNodeId) {
                throw std::runtime_error("RegionSplitUpdater can only move direct proper parts from the source node");
            }
            tree.moveProperPart(targetNodeId, sourceNodeId, pixelId);
        }
    }

    static std::vector<std::vector<PixelId>> connectedComponentsInRegion(
        AdjacencyRelation &adjacency,
        std::span<const PixelId> pixels,
        int domainSize) {
        std::vector<unsigned char> allowed(static_cast<size_t>(domainSize), 0);
        std::vector<unsigned char> visited(static_cast<size_t>(domainSize), 0);
        for (PixelId pixelId : pixels) {
            if (pixelId < 0 || pixelId >= domainSize) {
                throw std::runtime_error("RegionSplitUpdater received an out-of-bounds region pixel");
            }
            allowed[static_cast<size_t>(pixelId)] = 1;
        }

        std::vector<std::vector<PixelId>> components;
        std::vector<PixelId> stack;
        stack.reserve(pixels.size());

        for (PixelId seed : pixels) {
            if (visited[static_cast<size_t>(seed)] != 0) {
                continue;
            }

            std::vector<PixelId> component;
            stack.clear();
            stack.push_back(seed);
            visited[static_cast<size_t>(seed)] = 1;

            while (!stack.empty()) {
                const PixelId current = stack.back();
                stack.pop_back();
                component.push_back(current);

                for (PixelId neighbor : adjacency.getNeighborPixels(current)) {
                    if (neighbor < 0 || neighbor >= domainSize) {
                        continue;
                    }
                    const size_t neighborIndex = static_cast<size_t>(neighbor);
                    if (allowed[neighborIndex] == 0 || visited[neighborIndex] != 0) {
                        continue;
                    }
                    visited[neighborIndex] = 1;
                    stack.push_back(neighbor);
                }
            }

            std::sort(component.begin(), component.end());
            components.push_back(std::move(component));
        }

        std::sort(components.begin(),
                  components.end(),
                  [](const std::vector<PixelId> &lhs,
                     const std::vector<PixelId> &rhs) {
                      if (lhs.empty() || rhs.empty()) {
                          return lhs.size() < rhs.size();
                      }
                      return lhs.front() < rhs.front();
                  });
        return components;
    }

    static std::unordered_set<PixelId> makeRegionSet(const DynamicComponentTree &tree,
                                                     std::span<const PixelId> pixelsC) {
        std::unordered_set<PixelId> regionSet;
        regionSet.reserve(pixelsC.size());
        for (PixelId pixelId : pixelsC) {
            if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
                throw std::runtime_error("RegionSplitUpdater received an out-of-bounds region pixel");
            }
            regionSet.insert(pixelId);
        }
        return regionSet;
    }

    static bool nodeContainsMaxTreeLevel(const DynamicComponentTree &tree,
                                         NodeId nodeId,
                                         int level) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            return false;
        }

        const int altitude = tree.getAltitude(nodeId);
        if (level > altitude) {
            return false;
        }

        if (tree.isRoot(nodeId)) {
            return true;
        }

        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || !tree.isAlive(parentId)) {
            return false;
        }

        return level > tree.getAltitude(parentId);
    }

    static bool intersectsRegion(const std::vector<PixelId> &support,
                                 const std::unordered_set<PixelId> &regionSet) {
        for (PixelId pixelId : support) {
            if (regionSet.find(pixelId) != regionSet.end()) {
                return true;
            }
        }
        return false;
    }

    static std::vector<std::vector<PixelId>> connectedComponentsAfterRemovingRegion(
        AdjacencyRelation &adjacency,
        const std::vector<PixelId> &support,
        const std::unordered_set<PixelId> &regionSet) {
        std::unordered_set<PixelId> remaining;
        remaining.reserve(support.size());
        for (PixelId pixelId : support) {
            if (regionSet.find(pixelId) == regionSet.end()) {
                remaining.insert(pixelId);
            }
        }

        std::vector<std::vector<PixelId>> components;
        std::unordered_set<PixelId> visited;
        visited.reserve(remaining.size());
        std::vector<PixelId> stack;

        for (PixelId seed : remaining) {
            if (visited.find(seed) != visited.end()) {
                continue;
            }

            std::vector<PixelId> component;
            stack.clear();
            stack.push_back(seed);
            visited.insert(seed);

            while (!stack.empty()) {
                const PixelId current = stack.back();
                stack.pop_back();
                component.push_back(current);

                for (PixelId neighbor : adjacency.getNeighborPixels(current)) {
                    if (remaining.find(neighbor) == remaining.end()) {
                        continue;
                    }
                    if (visited.insert(neighbor).second) {
                        stack.push_back(neighbor);
                    }
                }
            }

            std::sort(component.begin(), component.end());
            components.push_back(std::move(component));
        }

        std::sort(components.begin(),
                  components.end(),
                  [](const std::vector<PixelId> &lhs,
                     const std::vector<PixelId> &rhs) {
                      if (lhs.empty() || rhs.empty()) {
                          return lhs.size() < rhs.size();
                      }
                      return lhs.front() < rhs.front();
                  });
        return components;
    }
};

} // namespace debug::region
