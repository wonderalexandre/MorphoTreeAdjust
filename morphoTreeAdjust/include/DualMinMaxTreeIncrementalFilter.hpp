#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <span>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "DynamicComponentTree.hpp"

/**
 * @brief Incremental updater for the dynamic min-tree / max-tree pair.
 *
 * This class implements the "update rather than rebuild" strategy for the
 * main per-pixel line. Given a set of nodes to prune in a primal tree, it
 * updates the dual tree in place so that both remain consistent with the
 * filtered image, avoiding a full rebuild of the dual tree after each change.
 *
 * Algorithmic intuition:
 *
 * - a subtree is removed from the primal tree;
 * - the affected proper parts define a set `C` in the dual tree;
 * - the nodes adjacent to `C` in the dual tree are collected;
 * - those nodes are organized by altitude in `mergeNodesByLevel`;
 * - a level sweep merges and reconnects the region until a consistent
 *   hierarchy is restored.
 *
 * Structures maintained by the class:
 *
 * - `mintree_`, `maxtree_`: the two mutable dynamic trees;
 * - `graph_`: shared adjacency of the image domain;
 * - `attrComputerMin_`, `attrComputerMax_`, and external buffers:
 *   incremental attribute maintenance after edits;
 * - `mergeNodesByLevel_`: merge buckets and frontier nodes of the current step;
 * - `properPartSetC_`, `removedMarks_`, and auxiliary buffers: collections and
 *   temporary marks for the current update.
 *
 * Important internal conventions:
 *
 * - `nodeCa`: extremal node that represents set `C` at level `a`;
 * - `b`: altitude of the parent of the removed subtree in the primal tree;
 * - `frontierNodesAboveB`: roots that remain above `b` and must be
 *   reattached when the sweep reaches their altitude;
 * - `properPartSetC`: pixels directly removed/relocated in the dual tree.
 *
 * Operational invariants:
 *
 * - at the beginning of the step, both trees are consistent with the same image;
 * - at the end of the step, the updated dual tree remains consistent with the
 *   new implicit image after pruning the primal tree;
 * - the adjuster never depends on rebuilding the entire dual tree;
 * - incremental attributes, when configured, remain synchronized with the
 *   edited structure.
 *
 * Granularity:
 *
 * This variant works in `subtree` mode: the caller provides subtree roots to
 * prune, and the adjuster processes the entire affected region in one shot.
 * The `leaf` variant lives in `DualMinMaxTreeIncrementalFilterLeaf`.
 *
 * Complexity:
 *
 * Use:
 *
 * - `P_C`: number of pixels in `properPartSetC`;
 * - `A_C`: number of adjacent nodes collected around `C`;
 * - `M_C`: total number of nodes placed in `mergeNodesByLevel`;
 * - `K_C`: number of nodes actually removed/merged/reconnected in the step.
 *
 * The cost of an adjustment step is local and depends on the affected region.
 * In asymptotic terms, it is dominated by:
 *
 * - collecting `properPartSetC` and removed nodes;
 * - computing adjacent nodes;
 * - building the paths up to `nodeCa` and the per-level buckets;
 * - the merge-and-reconnect sweep.
 *
 * In asymptotic notation, a typical step looks like:
 *
 * - `O(P_C + A_C + M_C + K_C)`,
 *
 * without depending directly on the total number of pixels in the image,
 * except when the affected region grows to the whole domain. This local nature
 * is precisely what makes incremental adjustment advantageous relative to a
 * global rebuild.
 *
 * Practical note:
 *
 * - this class is one of the main bottlenecks of the per-pixel `subtree`
 *   pipeline;
 * - its real cost depends heavily on how the threshold fragments or merges the
 *   affected region;
 * - the quality of the incremental update comes from trading a global cost for
 *   a cost proportional to the modified region.
 *
 * This implementation does not own the attribute buffers. The caller is
 * responsible for keeping them alive and indexed in the same global id space
 * as the associated `DynamicComponentTree`.
 */
template<typename PixelType = AltitudeType>
class DualMinMaxTreeIncrementalFilter {
public:
    /**
     * @brief Collection of merged and nested nodes grouped by altitude.
     * @details Over a single adjustment step, the structure stores:
     * - the `mergeNodesByLevel` buckets, indexed by altitude;
     * - the `frontierNodesAboveB` list, used when an adjacent node is above
     *   the sweep limit `b`;
     * - generation-based marks to avoid duplicate insertions during collection.
     */
    class MergedNodesCollection {
    private:
        static_assert(std::numeric_limits<PixelType>::is_integer, "MergedNodesCollection requires integral PixelType");
        static_assert(std::numeric_limits<PixelType>::lowest() == 0, "MergedNodesCollection expects non-negative contiguous levels");
        static constexpr std::size_t NumLevels = static_cast<std::size_t>(std::numeric_limits<PixelType>::max()) + 1;

        std::array<std::vector<NodeId>, NumLevels> mergeNodesByLevelStorage_;
        std::vector<PixelType> mergeLevels_;
        std::vector<NodeId> frontierNodesAboveB_;
        GenerationStampSet collectedNodeMarks_;
        GenerationStampSet mergeBucketNodeMarks_;
        GenerationStampSet adjacentSeedMarks_;
        std::size_t maxBucketSize_ = 0;
        int currentMergeLevelIndex_ = 0;
        bool isMaxtree_ = false;

    public:
        /**
         * @brief Creates the merged-node collection for a known node-id space.
         * @param maxNodes Expected maximum size of the node id space.
         */
        explicit MergedNodesCollection(int maxNodes = 0)
            : collectedNodeMarks_(std::max(0, maxNodes)),
              mergeBucketNodeMarks_(std::max(0, maxNodes)),
              adjacentSeedMarks_(std::max(0, maxNodes)) {}

        /**
         * @brief Reinitializes all temporary state for a new adjustment step.
         * @param isMaxtree `true` for a max-tree update and `false` for a min-tree update.
         */
        void resetCollection(bool isMaxtree) {
            isMaxtree_ = isMaxtree;
            for (auto level : mergeLevels_) {
                // The current backend uses one preallocated bucket for every
                // possible altitude, so it is enough to clear the levels that became active.
                mergeNodesByLevelStorage_[static_cast<std::size_t>(level)].clear();
            }
            mergeLevels_.clear();
            frontierNodesAboveB_.clear();
            collectedNodeMarks_.resetAll();
            mergeBucketNodeMarks_.resetAll();
            adjacentSeedMarks_.resetAll();
            maxBucketSize_ = 0;
            currentMergeLevelIndex_ = 0;
        }

        /**
         * @brief Returns the bucket associated with an altitude value.
         * @param level Queried level.
         * @return Mutable reference to the bucket corresponding to `level`.
         */
        std::vector<NodeId> &getMergedNodes(const PixelType &level) {
            return mergeNodesByLevelStorage_[static_cast<std::size_t>(level)];
        }

        /**
         * @brief Returns the frontier roots collected above the limit `b`.
         * @details Each stored node is the root of the first branch that left
         * the merge interval while climbing from a valid adjacent seed.
         */
        std::vector<NodeId> &getFrontierNodesAboveB() {
            return frontierNodesAboveB_;
        }

        /**
         * @brief Marks an adjacent seed only once in the current step.
         * @return `true` if the seed was accepted now; `false` if it had already been seen.
         */
        bool markAdjacentSeed(NodeId nodeId) {
            if (adjacentSeedMarks_.isMarked(nodeId)) {
                return false;
            }

            adjacentSeedMarks_.mark(nodeId);
            return true;
        }

        /**
         * @brief Returns the largest bucket observed during the current collection.
         * @details This value is used to reserve the sweep worklist without
         * reintroducing reallocations on the hot path.
         */
        std::size_t getMaxBucketSize() const {
            return maxBucketSize_;
        }

        /**
         * @brief Records a frontier root above `b` only once.
         */
        void addFrontierNodeAboveB(NodeId nodeId) {
            if (!collectedNodeMarks_.isMarked(nodeId)) {
                frontierNodesAboveB_.push_back(nodeId);
                collectedNodeMarks_.mark(nodeId);
            }
        }

        /**
         * @brief Inserts a single node into the bucket of its own level.
         * @details This operation corresponds more directly to partition
         * `{\Gamma_v^-}_{v \in I_C^*}` from Algorithm 2 of the paper: each
         * visited node enters only the bucket of its own level, without
         * immediately materializing the full path up to `nodeCa`.
         */
        void addMergeNode(const DynamicComponentTree &tree, NodeId nodeId) {
            ensureNodeCapacity(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
            if (collectedNodeMarks_.isMarked(nodeId)) {
                return;
            }

            auto &bucket = getMergedNodes(static_cast<PixelType>(tree.getAltitude(nodeId)));
            bucket.push_back(nodeId);
            maxBucketSize_ = std::max(maxBucketSize_, bucket.size());
            collectedNodeMarks_.mark(nodeId);
            mergeBucketNodeMarks_.mark(nodeId);
        }

        /**
         * @brief Inserts a detached synthetic node into the bucket of `level`.
         */
        void addSyntheticMergeNode(PixelType level, NodeId nodeId) {
            if (nodeId == InvalidNode) {
                return;
            }
            ensureNodeCapacity(static_cast<std::size_t>(nodeId + 1));
            auto &bucket = getMergedNodes(level);
            bucket.push_back(nodeId);
            maxBucketSize_ = std::max(maxBucketSize_, bucket.size());
            collectedNodeMarks_.mark(static_cast<std::size_t>(nodeId));
            mergeBucketNodeMarks_.mark(static_cast<std::size_t>(nodeId));
        }

        /**
         * @brief Tests whether a node was placed in any merge bucket of the step.
         */
        bool isMergeNode(NodeId nodeId) const {
            return mergeBucketNodeMarks_.isMarked(nodeId);
        }

        /**
         * @brief Grows generation-backed marks to cover the current node-id space.
         * @param maxNodes Minimum number of node slots that may be addressed.
         */
        void ensureNodeCapacity(std::size_t maxNodes) {
            collectedNodeMarks_.ensureSize(maxNodes);
            mergeBucketNodeMarks_.ensureSize(maxNodes);
            adjacentSeedMarks_.ensureSize(maxNodes);
        }


        /**
         * @brief Builds the ordered list of active levels and returns the first one.
         * @details Because the current implementation uses dense altitude-indexed
         * buckets, the collection scans the full discrete domain and keeps only
         * the levels that remained non-empty in the current step.
         */
        PixelType firstMergeLevel() {
            mergeLevels_.clear();
            // The altitude domain is discrete and dense in this backend, so the
            // scan visits all levels and keeps only the non-empty ones.
            for (std::size_t i = 0; i < mergeNodesByLevelStorage_.size(); ++i) {
                if (!mergeNodesByLevelStorage_[i].empty()) {
                    mergeLevels_.push_back(static_cast<PixelType>(i));
                }
            }

            if (mergeLevels_.empty()) {
                return PixelType{};
            }

            currentMergeLevelIndex_ = isMaxtree_ ? (int) mergeLevels_.size() - 1 : 0;
            return mergeLevels_[currentMergeLevelIndex_];
        }

        /**
         * @brief Returns `true` while there is still an active bucket to visit in the sweep.
         * @details The internal iterator does not traverse empty levels; it
         * navigates only the compact `mergeLevels_` list built in
         * `firstMergeLevel()`. This helper is the main guard of the sweep loop.
         */
        bool hasMergeLevel() const {
            return !mergeLevels_.empty() && currentMergeLevelIndex_ >= 0 && currentMergeLevelIndex_ < (int) mergeLevels_.size();
        }

        /**
         * @brief Advances the ordered level iterator and returns the next level.
         * @details In a max-tree the order is descending; in a min-tree it is
         * ascending. When the iterator leaves the valid range, the method
         * returns `PixelType{}` as a sentinel, and `hasMergeLevel()` starts
         * returning `false`.
         */
        PixelType nextMergeLevel() {
            currentMergeLevelIndex_ = isMaxtree_ ? currentMergeLevelIndex_ - 1 : currentMergeLevelIndex_ + 1;
            if (!hasMergeLevel()) {
                return PixelType{};
            }
            return mergeLevels_[currentMergeLevelIndex_];
        }
    };

private:
    // Dynamic primal/dual trees and shared adjacency of the domain.
    DynamicComponentTree *mintree_ = nullptr;
    DynamicComponentTree *maxtree_ = nullptr;
    AdjacencyRelation *graph_ = nullptr;

    // Incremental computers and external buffers used after local edits.
    DynamicAttributeComputer *attrComputerMin_ = nullptr;
    DynamicAttributeComputer *attrComputerMax_ = nullptr;
    std::span<float> bufferMin_;
    std::span<float> bufferMax_;

    // Temporary state for the current step.
    MergedNodesCollection mergeNodesByLevel_;
    GenerationStampSet removedMarks_;
    std::vector<NodeId> removedNodesPendingAbsorption_;
    GenerationStampSet pixelsInCMarks_;
    GenerationStampSet climbedNodeMarks_;
    GenerationStampSet attributeUpdateMarks_;
    mutable GenerationStampSet directComponentInRegionMarks_;
    mutable GenerationStampSet directComponentVisitedMarks_;
    std::vector<PixelId> properPartSetC_;
    std::vector<NodeId> nodesPendingRemoval_;
    mutable std::vector<PixelId> directComponentSeeds_;
    mutable std::vector<PixelId> directComponentStack_;
    mutable std::vector<PixelId> directComponent_;
    PixelType altitudeCa = PixelType{};

    /**
     * @brief Recomputes the incremental attribute of a single node after a local edit.
     * @details The recomputation uses the registered incremental computer and
     * applies the pre-processing, child-merge, and post-processing phases
     * directly on the buffer corresponding to the tree.
     * @param tree Tree that contains the edited node.
     * @param nodeId Node whose attribute must be updated.
     * @param enforceDualUpdateInterval When true, recomputes only nodes whose
     * altitude is inside the active dual-update interval anchored at
     * `altitudeCa`; when false, recomputes the marked node immediately. The
     * false path is used for direct primal updates that must not be filtered by
     * the dual sweep interval.
     */
    void computeAttributeOnTreeNode(DynamicComponentTree *tree, NodeId nodeId, bool enforceDualUpdateInterval = true) {
        DynamicAttributeComputer *computer = nullptr;
        std::span<float> buffer;
        if (tree == mintree_) {
            computer = attrComputerMin_;
            buffer = bufferMin_;
        } else if (tree == maxtree_) {
            computer = attrComputerMax_;
            buffer = bufferMax_;
        }
        if (tree == nullptr || computer == nullptr || buffer.empty() || nodeId == InvalidNode || !tree->isNode(nodeId) || !tree->isAlive(nodeId)) {
            return;
        }

        if (!attributeUpdateMarks_.isMarked(static_cast<std::size_t>(nodeId))) {
            return;
        }

        if (enforceDualUpdateInterval) {
            const PixelType level = static_cast<PixelType>(tree->getAltitude(nodeId));
            if ((tree == maxtree_ && level < altitudeCa) || (tree == mintree_ && level > altitudeCa)) {
                attributeUpdateMarks_.unmark(static_cast<std::size_t>(nodeId));
                return;
            }
        }

        computer->preProcessing(nodeId, buffer);
        for (NodeId childId : tree->getChildren(nodeId)) {
            computer->mergeProcessing(nodeId, childId, buffer);
        }
        computer->postProcessing(nodeId, buffer);

        attributeUpdateMarks_.unmark(static_cast<std::size_t>(nodeId));
    }

    /**
     * @brief Resolves the attribute computer associated with the given tree.
     * @details Because the adjuster maintains a min-tree/max-tree pair, this
     * helper centralizes the choice of the correct incremental runtime without
     * exposing that logic at local mutation points.
     */
    DynamicAttributeComputer *getAttributeComputer(DynamicComponentTree *tree) const {
        if (tree == mintree_) {
            return attrComputerMin_;
        }
        if (tree == maxtree_) {
            return attrComputerMax_;
        }
        return nullptr;
    }

    /**
     * @brief Forwards the definitive removal of a node to the attribute computer.
     * @details The hook is used by attributes with a persistent cache per
     * `NodeId`, allowing internal summaries to be cleared when the node slot is
     * released.
     */
    void notifyNodeRemoved(DynamicComponentTree *tree, NodeId nodeId) {
        auto *computer = getAttributeComputer(tree);
        if (computer != nullptr && nodeId != InvalidNode) {
            computer->onNodeRemoved(nodeId);
        }
    }

    /**
     * @brief Forwards a full merge of direct proper parts to the attribute.
     * @details The structural adjustment first performs the notification and
     * only then applies the real mutation in the tree backend, preserving
     * access to the previous state for the attribute when necessary.
     */
    void notifyMoveProperParts(DynamicComponentTree *tree, NodeId targetNodeId, NodeId sourceNodeId) {
        auto *computer = getAttributeComputer(tree);
        if (computer != nullptr) {
            computer->onMoveProperParts(targetNodeId, sourceNodeId);
        }
    }

    /**
     * @brief Forwards the transfer of a single proper part to the attribute.
     * @details This is the fine-grained hook used by geometric attributes for
     * incremental local-state updates.
     */
    void notifyMoveProperPart(DynamicComponentTree *tree, NodeId targetNodeId, NodeId sourceNodeId, PixelId pixelId) {
        auto *computer = getAttributeComputer(tree);
        if (computer != nullptr) {
            computer->onMoveProperPart(targetNodeId, sourceNodeId, pixelId);
        }
    }

    /**
     * @brief Marks a node for attribute update in the current step.
     * @details Marking is idempotent and uses the global id space of the tree.
     * Marking does not recompute immediately; it only schedules the node for
     * when the flow reaches the appropriate point of structural stabilization.
     */
    void markAttributeUpdate(NodeId nodeId) {
        if (nodeId != InvalidNode) {
            attributeUpdateMarks_.mark(static_cast<std::size_t>(nodeId));
        }
    }

    /**
     * @brief Merges a primal node into its parent and refreshes the impacted primal attribute.
     * @details This helper is used after the paired dual updates have succeeded,
     * so the primal mutation and its local attribute cache update remain a single
     * final step.
     */
    void mergePrimalNodeIntoParentAndUpdateAttribute(DynamicComponentTree *primalTree, NodeId nodeId, NodeId parentId) {
        notifyMoveProperParts(primalTree, parentId, nodeId);
        primalTree->mergeNodeIntoParent(nodeId);
        notifyNodeRemoved(primalTree, nodeId);

        markAttributeUpdate(parentId);
        computeAttributeOnTreeNode(primalTree, parentId, false);
    }

    /**
     * @brief Grows temporary node-indexed marks to the current tree id space.
     */
    void ensureUpdateStepNodeCapacity() {
        const std::size_t minSlots = static_cast<std::size_t>(std::max(mintree_ ? mintree_->getNumInternalNodeSlots() : 0, maxtree_ ? maxtree_->getNumInternalNodeSlots() : 0));
        mergeNodesByLevel_.ensureNodeCapacity(minSlots);
        removedMarks_.ensureSize(minSlots);
        climbedNodeMarks_.ensureSize(minSlots);
        attributeUpdateMarks_.ensureSize(minSlots);
    }

    /**
     * @brief Initializes the temporary state shared by one `updateTree` step.
     * @details This routine concentrates all resets that do not depend on the
     * size of the merge buckets not yet built. The goal is to make it explicit
     * that each step begins with no marks, no residual `C` collection, and no
     * pending nodes inherited from the previous step.
     */
    void initializeUpdateStepState(bool isMaxtree) {
        ensureUpdateStepNodeCapacity();
        attributeUpdateMarks_.resetAll();
        altitudeCa = PixelType{};
        pixelsInCMarks_.resetAll();
        properPartSetC_.clear();
        properPartSetC_.reserve(64);
        mergeNodesByLevel_.resetCollection(isMaxtree);
        climbedNodeMarks_.resetAll();
        removedMarks_.resetAll();
        removedNodesPendingAbsorption_.clear();
        nodesPendingRemoval_.clear();
    }

    /**
     * @brief Tests whether a marked node can still be absorbed during closure.
     * @details An absorbable node must exist, remain alive, stay marked in
     * `removedMarks_`, and have no direct proper parts. These guards filter
     * obsolete seeds and nodes already stabilized by later mutations.
     */
    bool canAbsorbRemovedNode(DynamicComponentTree *dualTree, NodeId nodeId) const {
        return nodeId != InvalidNode && dualTree->isNode(nodeId) && dualTree->isAlive(nodeId) && removedMarks_.isMarked(nodeId) && dualTree->getNumProperParts(nodeId) == 0;
    }

    /**
     * @brief Absorbs a removed root by choosing a new root representative.
     * @details The method selects, among the remaining immediate children, the
     * representative compatible with the altitude order of the current tree and
     * reattaches the others under it before releasing the removed root.
     */
    void absorbRemovedRootNode(DynamicComponentTree *dualTree, NodeId removedNodeId, bool preserveCanonicalForm) {
        const NodeId firstChild = dualTree->getFirstChild(removedNodeId);
        if (firstChild == InvalidNode) {
            return;
        }

        NodeId newRoot = firstChild;
        bool newRootChanged = false;
        for (NodeId childId = firstChild; childId != InvalidNode; childId = dualTree->getNextSibling(childId)) {
            if ((dualTree->isMaxtree() && dualTree->getAltitude(childId) < dualTree->getAltitude(newRoot)) || (!dualTree->isMaxtree() && dualTree->getAltitude(childId) > dualTree->getAltitude(newRoot))) {
                newRoot = childId;
            }
        }

        for (NodeId childId = firstChild; childId != InvalidNode;) {
            const NodeId next = dualTree->getNextSibling(childId);
            if (childId != newRoot && !dualTree->hasChild(newRoot, childId)) {
                if (preserveCanonicalForm) {
                    attachChildOrMergeSameLevelStructure(dualTree, newRoot, childId);
                } else {
                    dualTree->detachNode(childId);
                    dualTree->attachNode(newRoot, childId);
                }
                newRootChanged = true;
            }
            childId = next;
        }

        dualTree->setRoot(newRoot);
        dualTree->releaseNode(removedNodeId);
        notifyNodeRemoved(dualTree, removedNodeId);
        if (newRootChanged) {
            markAttributeUpdate(newRoot);
        }
        computeAttributeOnTreeNode(dualTree, newRoot);
    }

    /**
     * @brief Absorbs a removed non-root node by merging its contents into the parent.
     * @details After moving children and proper parts to `parentId`, the method
     * recomputes the parent's incremental attribute and queues it again for
     * absorption if it has itself become empty of proper parts.
     * @return The parent queued again for absorption, or `InvalidNode` when
     * there is no new propagation to perform.
     */
    NodeId absorbRemovedNonRootNode(DynamicComponentTree *dualTree, NodeId removedNodeId, bool preserveCanonicalForm) {
        const NodeId parentId = dualTree->getNodeParent(removedNodeId);
        if (parentId == InvalidNode || parentId == removedNodeId || !dualTree->isAlive(parentId)) {
            return InvalidNode;
        }

        const bool parentChanged = dualTree->getFirstChild(removedNodeId) != InvalidNode;
        if (preserveCanonicalForm) {
            mergeNodeIntoParentKeepingCanonicalForm(dualTree, parentId, removedNodeId);
        } else {
            mergedParentAndChildren(dualTree, parentId, removedNodeId);
        }
        disconnect(dualTree, removedNodeId, true);
        if (parentChanged) {
            markAttributeUpdate(parentId);
        }
        computeAttributeOnTreeNode(dualTree, parentId);

        if (dualTree->isAlive(parentId) && dualTree->getNumProperParts(parentId) == 0) {
            removedMarks_.mark(parentId);
            return parentId;
        }

        return InvalidNode;
    }

    /**
     * @brief Closes the final reconnection of the step and eliminates the remaining empty nodes.
     * @details This phase runs after the level sweep, when the local topology
     * is already stable. The explicit stack traverses in post-order only the
     * nodes still marked and empty of proper parts, ignoring seeds that became
     * obsolete during earlier mutations of the step.
     */
    void absorbRemovedNodes(DynamicComponentTree *dualTree, const std::vector<NodeId> &removedNodeIds, bool preserveCanonicalForm = false) {
        struct Frame {
            NodeId nodeId = InvalidNode;
            NodeId nextChildId = InvalidNode;
        };

        if (dualTree == nullptr || removedNodeIds.empty()) {
            return;
        }

        std::vector<Frame> stack;
        stack.reserve(std::max<std::size_t>(64, removedNodeIds.size()));
        const auto makeFrame = [dualTree](NodeId nodeId) {
            return Frame{nodeId, dualTree->getFirstChild(nodeId)};
        };
        for (auto it = removedNodeIds.rbegin(); it != removedNodeIds.rend(); ++it) {
            if (*it != InvalidNode) {
                stack.push_back(makeFrame(*it));
            }
        }

        while (!stack.empty()) {
            Frame &frame = stack.back();
            NodeId currentNodeId = InvalidNode;

            if (!canAbsorbRemovedNode(dualTree, frame.nodeId)) {
                stack.pop_back();
            } else if (frame.nextChildId != InvalidNode) {
                const NodeId childId = frame.nextChildId;
                frame.nextChildId = dualTree->getNextSibling(childId);
                stack.push_back(makeFrame(childId));
            } else {
                currentNodeId = frame.nodeId;
                stack.pop_back();
            }

            if (canAbsorbRemovedNode(dualTree, currentNodeId)) {
                if (dualTree->isRoot(currentNodeId)) {
                    absorbRemovedRootNode(dualTree, currentNodeId, preserveCanonicalForm);
                } else {
                    const NodeId parentId = absorbRemovedNonRootNode(dualTree, currentNodeId, preserveCanonicalForm);
                    if (parentId != InvalidNode) {
                        stack.push_back(makeFrame(parentId));
                    }
                }
            }
        }
    }

    /**
     * @brief Normalizes a child branch of `rootId` until reaching its surviving representative.
     * @details While the top of that branch remains marked in `removedMarks_`
     * and has no proper parts, the best child is promoted to sit directly under
     * `rootId`. The goal is to expose, already at this level, the
     * representative that will actually survive the closure of the step.
     *
     * This is used in the special case where `nodeCa` is the root: before
     * choosing the new tree root, each child branch is reduced to its immediate
     * representative, avoiding the selection of a node that would still be
     * contracted immediately afterward.
     *
     * @param dualTree Dual tree being updated.
     * @param rootId Current root under which the branch will be normalized.
     * @param childId Initial direct child of the branch.
     * @return The immediate surviving representative of that branch under `rootId`.
     */
    NodeId collapseRemovedRootBranch(DynamicComponentTree *dualTree, NodeId rootId, NodeId childId, bool preserveCanonicalForm) {
        assert(dualTree != nullptr);
        assert(rootId != InvalidNode);

        NodeId current = childId;
        while (current != InvalidNode && dualTree->isNode(current) && dualTree->isAlive(current) && dualTree->getNodeParent(current) == rootId && removedMarks_.isMarked(current) && dualTree->getNumProperParts(current) == 0) {
            const NodeId firstGrandchild = dualTree->getFirstChild(current);
            if (firstGrandchild == InvalidNode) {
                break; // A branch without grandchildren has already reached the final representative of this path.
            }

            NodeId promoted = firstGrandchild;
            for (NodeId grandchildId = dualTree->getFirstChild(current); grandchildId != InvalidNode; grandchildId = dualTree->getNextSibling(grandchildId)) {
                if ((dualTree->isMaxtree() && dualTree->getAltitude(grandchildId) < dualTree->getAltitude(promoted)) || (!dualTree->isMaxtree() && dualTree->getAltitude(grandchildId) > dualTree->getAltitude(promoted))) {
                    promoted = grandchildId;
                }
            }

            if (!dualTree->isRoot(promoted)) {
                disconnect(dualTree, promoted, false);
            }
            dualTree->attachNode(rootId, promoted);

            for (NodeId grandchildId = dualTree->getFirstChild(current); grandchildId != InvalidNode;) {
                const NodeId next = dualTree->getNextSibling(grandchildId);
                if (grandchildId != promoted && !dualTree->hasChild(promoted, grandchildId)) {
                    if (preserveCanonicalForm) {
                        attachChildOrMergeSameLevelStructure(dualTree, promoted, grandchildId);
                    } else {
                        if (!dualTree->isRoot(grandchildId)) {
                            disconnect(dualTree, grandchildId, false);
                        }
                        dualTree->attachNode(promoted, grandchildId);
                    }
                }
                grandchildId = next;
            }

            disconnect(dualTree, current, true);
            markAttributeUpdate(promoted);
            computeAttributeOnTreeNode(dualTree, promoted);
            current = promoted;
        }

        return current;
    }

    /**
     * @brief Closes the final reconnection among `nodeCa`, `finalUnionNode`, and the removed nodes.
     * @details It first resolves the final topological position of the union
     * node produced by the level sweep. It then consumes
     * `removedNodesPendingAbsorption_` to contract, in post-order, the nodes
     * that remained alive but empty until the step was closed.
     */
    void finalizeUpdateTreeAndContractRemovedNodes(DynamicComponentTree *dualTree,
                                                   NodeId nodeCa,
                                                   NodeId finalUnionNode,
                                                   bool preserveCanonicalForm = false) {
        if (dualTree == nullptr) {
            return;
        }

        if (finalUnionNode != InvalidNode && dualTree->isAlive(finalUnionNode)) {
            // First resolve the final topological position of the union node
            // produced by the level sweep.
            if (removedMarks_.isMarked(nodeCa)) {
                // If `nodeCa` was emptied, the final union node must occupy its
                // topological position in the updated hierarchy.
                if (!dualTree->isRoot(nodeCa)) {
                    const NodeId nodeCaParentId = dualTree->getNodeParent(nodeCa);
                    bool finalUnionNodeChanged = false;

                    if (!dualTree->isRoot(finalUnionNode)) {
                        disconnect(dualTree, finalUnionNode, false);
                    }
                    NodeId attachmentTarget = finalUnionNode;
                    if (preserveCanonicalForm) {
                        attachChildOrMergeSameLevelStructure(dualTree, nodeCaParentId, finalUnionNode);
                        if (!dualTree->isAlive(finalUnionNode)) {
                            attachmentTarget = nodeCaParentId;
                        }
                    } else {
                        dualTree->attachNode(nodeCaParentId, finalUnionNode);
                    }

                    for (NodeId n = dualTree->getFirstChild(nodeCa); n != InvalidNode;) {
                        const NodeId next = dualTree->getNextSibling(n);
                        if (n != attachmentTarget && dualTree->isAlive(attachmentTarget) && !dualTree->hasChild(attachmentTarget, n)) {
                            if (preserveCanonicalForm) {
                                attachChildOrMergeSameLevelStructure(dualTree, attachmentTarget, n);
                            } else {
                                if (!dualTree->isRoot(n)) {
                                    disconnect(dualTree, n, false);
                                }
                                dualTree->attachNode(attachmentTarget, n);
                            }
                            finalUnionNodeChanged = true;
                        }
                        n = next;
                    }

                    disconnect(dualTree, nodeCa, true);
                    if (finalUnionNodeChanged) {
                        markAttributeUpdate(attachmentTarget);
                    }
                    computeAttributeOnTreeNode(dualTree, attachmentTarget);
                    markAttributeUpdate(nodeCaParentId);
                    computeAttributeOnTreeNode(dualTree, nodeCaParentId);
                } else {
                    // When `nodeCa` is the root, we first normalize each direct
                    // child branch to its surviving representative. The new root
                    // must be chosen among those representatives, not among
                    // possibly transient immediate children.
                    NodeId survivingFinalUnionNode = finalUnionNode;
                    for (NodeId childId = dualTree->getFirstChild(nodeCa); childId != InvalidNode;) {
                        const NodeId next = dualTree->getNextSibling(childId);
                        const NodeId normalizedChild = collapseRemovedRootBranch(dualTree, nodeCa, childId, preserveCanonicalForm);
                        if (childId == finalUnionNode && normalizedChild != InvalidNode) {
                            survivingFinalUnionNode = normalizedChild;
                        }
                        childId = next;
                    }

                    NodeId candidateRootId = survivingFinalUnionNode;
                    bool candidateRootChanged = false;
                    // When `nodeCa` was the root, a new root is chosen
                    // consistently with the altitude order of the current tree.
                    for (NodeId n = dualTree->getFirstChild(nodeCa); n != InvalidNode; n = dualTree->getNextSibling(n)) {
                        if ((dualTree->isMaxtree() && dualTree->getAltitude(n) < dualTree->getAltitude(candidateRootId)) ||
                            (!dualTree->isMaxtree() && dualTree->getAltitude(n) > dualTree->getAltitude(candidateRootId))) {
                            candidateRootId = n;
                        }
                    }

                    if (candidateRootId != survivingFinalUnionNode) {
                        if (preserveCanonicalForm) {
                            attachChildOrMergeSameLevelStructure(dualTree, candidateRootId, survivingFinalUnionNode);
                        } else {
                            if (!dualTree->isRoot(survivingFinalUnionNode)) {
                                disconnect(dualTree, survivingFinalUnionNode, false);
                            }
                            dualTree->attachNode(candidateRootId, survivingFinalUnionNode);
                        }
                        candidateRootChanged = true;
                    }

                    for (NodeId n = dualTree->getFirstChild(nodeCa); n != InvalidNode;) {
                        const NodeId next = dualTree->getNextSibling(n);
                        if (n != candidateRootId && !dualTree->hasChild(candidateRootId, n)) {
                            if (preserveCanonicalForm) {
                                attachChildOrMergeSameLevelStructure(dualTree, candidateRootId, n);
                            } else {
                                if (!dualTree->isRoot(n)) {
                                    disconnect(dualTree, n, false);
                                }
                                dualTree->attachNode(candidateRootId, n);
                            }
                            candidateRootChanged = true;
                        }
                        n = next;
                    }

                    dualTree->setRoot(candidateRootId);
                    dualTree->releaseNode(nodeCa);
                    notifyNodeRemoved(dualTree, nodeCa);
                    if (candidateRootChanged) {
                        markAttributeUpdate(candidateRootId);
                    }
                    computeAttributeOnTreeNode(dualTree, candidateRootId);
                }
            } else if (finalUnionNode != nodeCa) {
                // If `nodeCa` survived, the final union node becomes its child.
                if (!dualTree->isRoot(finalUnionNode)) {
                    disconnect(dualTree, finalUnionNode, false);
                }
                if (preserveCanonicalForm) {
                    attachChildOrMergeSameLevelStructure(dualTree, nodeCa, finalUnionNode);
                } else {
                    dualTree->attachNode(nodeCa, finalUnionNode);
                }
                markAttributeUpdate(nodeCa);
                computeAttributeOnTreeNode(dualTree, nodeCa);
            }
        }

        absorbRemovedNodes(dualTree, removedNodesPendingAbsorption_, preserveCanonicalForm);
    }

    /**
     * @brief Checks the structural post-condition of `updateTree`.
     * @details In debug builds, guarantees that no live node remains without
     * proper parts at the end of the step. In release builds this routine has
     * no runtime cost.
     */
    void assertAllAliveNodesHaveProperParts(DynamicComponentTree *tree) const {
#ifndef NDEBUG
        assert(tree != nullptr);
        for (NodeId nodeId = 0; nodeId < tree->getNumInternalNodeSlots(); ++nodeId) {
            if (!tree->isNode(nodeId) || !tree->isAlive(nodeId)) {
                continue; // Nonexistent or dead nodes do not participate in the post-condition.
            }
            assert(tree->getNumProperParts(nodeId) > 0);
        }
#else
        (void) tree;
#endif
    }

    /**
     * @brief Disconnects a node from its parent, optionally releasing its slot.
     * @details The operation delegates to `DynamicComponentTree::removeChild`,
     * preserving the root case and ignoring nodes with an invalid parent.
     *
     * When `releaseNode` is `true`, the removal is definitive: the node slot
     * returns to the tree's free pool and the attribute computer is notified to
     * discard any persistent state associated with the `NodeId`.
     *
     * @param tree Tree that contains `nodeId`.
     * @param nodeId Node to disconnect.
     * @param releaseNode `true` to release the node slot after detaching.
     */
    void disconnect(DynamicComponentTree *tree, NodeId nodeId, bool releaseNode) {
        assert(tree != nullptr);
        if (tree->isRoot(nodeId)) {
            return;
        }

        const NodeId parentId = tree->getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId) {
            return;
        }

        tree->removeChild(parentId, nodeId, releaseNode);
        if (releaseNode) {
            notifyNodeRemoved(tree, nodeId);
        }
    }

    /**
     * @brief Reattaches to `parentId` all direct children and proper parts of `childId`.
     * @details This is the basic structural merge operation used during the
     * level sweep. Node `childId` transfers its children and proper parts to
     * `parentId` before being disconnected from the hierarchy.
     *
     * On the attribute side, the full proper-part merge is notified before the
     * actual `moveProperParts` so the incremental runtime can exploit the full
     * transfer of the local summary.
     *
     * @param tree Updated tree.
     * @param parentId Destination node of the merge.
     * @param childId Node absorbed by `parentId`.
     */
    void mergedParentAndChildren(DynamicComponentTree *tree, NodeId parentId, NodeId childId) {
        assert(tree != nullptr);
        tree->moveChildren(parentId, childId);
        notifyMoveProperParts(tree, parentId, childId);
        tree->moveProperParts(parentId, childId);
    }

    /**
     * @brief Releases a node that has already transferred its children and proper parts.
     * @details The node may still be attached to its parent or may already be
     * detached as an isolated slot; both cases notify the attribute runtime.
     */
    void releaseMergedNode(DynamicComponentTree *tree, NodeId nodeId) {
        assert(tree != nullptr);
        if (nodeId == InvalidNode || !tree->isNode(nodeId) || !tree->isAlive(nodeId)) {
            return;
        }

        const NodeId parentId = tree->getNodeParent(nodeId);
        if (parentId != InvalidNode && parentId != nodeId) {
            disconnect(tree, nodeId, true);
        } else {
            tree->releaseNode(nodeId);
            notifyNodeRemoved(tree, nodeId);
        }
    }

    /**
     * @brief Merges `childId` into `parentId` without creating same-level child edges.
     * @details Direct children of `childId` are moved one by one through
     * `attachChildOrMergeSameLevelStructure`, so a child with the same altitude
     * as `parentId` is absorbed immediately instead of being left for a later
     * canonicalization pass.
     */
    void mergeNodeIntoParentKeepingCanonicalForm(DynamicComponentTree *tree, NodeId parentId, NodeId childId) {
        assert(tree != nullptr);
        assert(parentId != InvalidNode);
        assert(childId != InvalidNode);

        for (NodeId grandchildId = tree->getFirstChild(childId); grandchildId != InvalidNode;) {
            const NodeId next = tree->getNextSibling(grandchildId);
            if (grandchildId != parentId && tree->isAlive(grandchildId)) {
                if (!tree->isRoot(grandchildId)) {
                    disconnect(tree, grandchildId, false);
                }
                attachChildOrMergeSameLevelStructure(tree, parentId, grandchildId);
            }
            grandchildId = next;
        }

        notifyMoveProperParts(tree, parentId, childId);
        tree->moveProperParts(parentId, childId);
        markAttributeUpdate(parentId);
    }

    /**
     * @brief Attaches a child to a parent while preserving canonical level edges.
     * @details If parent and child have the same altitude, the child is merged
     * into the parent immediately. Otherwise the child is detached from any old
     * parent and attached normally. Attribute recomputation is only marked here;
     * callers decide when the marked node is recomputed.
     * @return `true` when the topology or proper-part ownership changed.
     */
    bool attachChildOrMergeSameLevelStructure(DynamicComponentTree *tree, NodeId parentId, NodeId childId) {
        assert(tree != nullptr);
        if (parentId == InvalidNode || childId == InvalidNode || parentId == childId) {
            return false;
        }
        if (!tree->isNode(parentId) || !tree->isNode(childId) ||
            !tree->isAlive(parentId) || !tree->isAlive(childId)) {
            return false;
        }

        if (tree->getAltitude(parentId) == tree->getAltitude(childId)) {
            if (tree->isRoot(childId)) {
                throw std::runtime_error("attachChildOrMergeSameLevel cannot merge the current root into a same-level parent");
            }

            mergeNodeIntoParentKeepingCanonicalForm(tree, parentId, childId);
            releaseMergedNode(tree, childId);
            return true;
        }

        if (tree->hasChild(parentId, childId)) {
            return false;
        }
        if (tree->isRoot(childId)) {
            throw std::runtime_error("attachChildOrMergeSameLevel cannot attach the current root under a non-root parent");
        }

        const NodeId oldParentId = tree->getNodeParent(childId);
        if (oldParentId != InvalidNode && oldParentId != childId) {
            disconnect(tree, childId, false);
        }
        tree->attachNode(parentId, childId);
        markAttributeUpdate(parentId);
        return true;
    }

    /**
     * @brief Canonical attach/merge helper that also recomputes the destination attribute.
     * @return `true` when the structure changed.
     */
    bool attachChildOrMergeSameLevel(DynamicComponentTree *tree, NodeId parentId, NodeId childId) {
        const bool changed = attachChildOrMergeSameLevelStructure(tree, parentId, childId);
        if (changed) {
            computeAttributeOnTreeNode(tree, parentId);
        }
        return changed;
    }

    /**
     * @brief Reattaches to `targetNodeId` only the direct children of `sourceNodeId` that lie outside `Γ[a,b]`.
     * @details Children still marked as belonging to the merge interval are
     * detached and temporarily isolated so they can be handled at the
     * corresponding sweep level. The remaining children are reattached to
     * `targetNodeId`.
     *
     * This routine separates two responsibilities:
     * - preserving, outside the interval, the connectivity that must rise
     *   together with the current union node;
     * - keeping, inside the interval, only the branches that will still be
     *   handled by the sweep buckets.
     */
    void reattachOutsideIntervalChildren(DynamicComponentTree *tree, NodeId targetNodeId, NodeId sourceNodeId, bool preserveCanonicalForm = false) {
        assert(tree != nullptr);
        if (sourceNodeId == targetNodeId) {
            for (NodeId childId = tree->getFirstChild(sourceNodeId); childId != InvalidNode;) {
                const NodeId next = tree->getNextSibling(childId);
                if (mergeNodesByLevel_.isMergeNode(childId)) {
                    tree->detachNode(childId);
                }
                childId = next;
            }
            return;
        }

        for (NodeId childId = tree->getFirstChild(sourceNodeId); childId != InvalidNode;) {
            const NodeId next = tree->getNextSibling(childId);
            if (mergeNodesByLevel_.isMergeNode(childId)) {
                tree->detachNode(childId);
            }
            childId = next;
        }

        if (preserveCanonicalForm) {
            for (NodeId childId = tree->getFirstChild(sourceNodeId); childId != InvalidNode;) {
                const NodeId next = tree->getNextSibling(childId);
                if (childId != targetNodeId && tree->isAlive(childId)) {
                    disconnect(tree, childId, false);
                    attachChildOrMergeSameLevelStructure(tree, targetNodeId, childId);
                }
                childId = next;
            }
        } else {
            tree->moveChildren(targetNodeId, sourceNodeId);
        }
    }

    /**
     * @brief Moves the removed proper-part set to the current union node.
     * @details Any node left without proper parts after the transfer is marked
     * in `removedMarks_` so it can be eliminated in later sweep stages without
     * being processed twice.
     *
     * This is the point in the algorithm where set `C` effectively enters the
     * updated dual tree. Each pixel is moved from the smallest live component
     * that still contains it to the level-`b` union node.
     *
     * @param tree Updated tree.
     * @param unionNode Union node that receives the proper parts.
     * @param properPartSetC Set `C` collected from the subtree removed in the primal tree.
     */
    void moveSelectedProperPartsToNode(DynamicComponentTree *tree, NodeId unionNode, const std::vector<PixelId> &properPartSetC) {
        for (auto pixelId : properPartSetC) {
            const NodeId ownerId = tree->getSmallestComponent(pixelId);
            if (ownerId == InvalidNode || ownerId == unionNode) {
                continue; // Ignore pixels without a current owner or already belonging to the union node.
            }

            tree->moveProperPart(unionNode, ownerId, pixelId);
            notifyMoveProperPart(tree, unionNode, ownerId, pixelId);

            if (tree->isAlive(ownerId) && tree->getNumProperParts(ownerId) == 0 && ownerId != unionNode && !removedMarks_.isMarked(ownerId)) {
                removedMarks_.mark(ownerId);
                removedNodesPendingAbsorption_.push_back(ownerId);
            }
        }
    }

    /**
     * @brief Creates the synthetic level-`b` union node required by direct merge updates.
     * @details The node is created only when the sweep interval reaches `b` and
     * no bucket already exists for that level. Subtree prune updates pass
     * `allowTargetNodeCreation=false`, so this path is used by the clean
     * merge-node API.
     * @return The created detached node, or `InvalidNode` when no synthetic node is needed.
     */
    NodeId createTargetLevelUnionNodeIfMissing(DynamicComponentTree *dualTree, PixelType b, bool isMaxtree) {
        if (dualTree == nullptr) {
            return InvalidNode;
        }

        const bool reachesTargetLevel = (isMaxtree && b > altitudeCa) ||
                                        (!isMaxtree && b < altitudeCa);
        if (!reachesTargetLevel) {
            return InvalidNode;
        }

        auto &targetBucket = mergeNodesByLevel_.getMergedNodes(b);
        if (!targetBucket.empty()) {
            return InvalidNode;
        }

        const NodeId targetNode = dualTree->createDetachedNode(static_cast<int>(b));
        if (targetNode == InvalidNode) {
            throw std::runtime_error("explicit region update could not create a level-b union node");
        }

        ensureUpdateStepNodeCapacity();
        mergeNodesByLevel_.addSyntheticMergeNode(b, targetNode);
        markAttributeUpdate(targetNode);
        return targetNode;
    }

    /**
     * @brief Runs the explicit region update on the dual tree.
     * @details The method builds the affected buckets, sweeps levels between
     * `b` and `altitudeCa`, inserts set `C` at level `b`, reconnects frontier
     * branches, and finalizes the updated hierarchy in canonical form.
     */
    void runExplicitRegionUpdate(DynamicComponentTree *dualTree, PixelType b, bool isMaxtree, NodeId nodeCa, bool allowTargetNodeCreation) {
        buildMergedAndNestedCollections(dualTree, properPartSetC_, b, isMaxtree);
        const NodeId syntheticTargetNode = allowTargetNodeCreation ? createTargetLevelUnionNodeIfMissing(dualTree, b, isMaxtree) : InvalidNode;

        PixelType currentMergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId currentUnionNode = InvalidNode;
        NodeId previousLevelUnionNode = InvalidNode;
        nodesPendingRemoval_.reserve(mergeNodesByLevel_.getMaxBucketSize());

        while (mergeNodesByLevel_.hasMergeLevel() &&
               ((isMaxtree && currentMergeLevel > altitudeCa) ||
                (!isMaxtree && currentMergeLevel < altitudeCa))) {
            auto &nodesAtCurrentLevel = mergeNodesByLevel_.getMergedNodes(currentMergeLevel);
            currentUnionNode = InvalidNode;
            nodesPendingRemoval_.clear();

            for (auto nodeId : nodesAtCurrentLevel) {
                if (!dualTree->isAlive(nodeId)) {
                    continue;
                } else if (currentUnionNode == InvalidNode && !removedMarks_.isMarked(nodeId)) {
                    currentUnionNode = nodeId;
                    if (!dualTree->isRoot(currentUnionNode)) {
                        disconnect(dualTree, currentUnionNode, false);
                    }

                    for (auto pendingNodeId : nodesPendingRemoval_) {
                        reattachOutsideIntervalChildren(dualTree, currentUnionNode, pendingNodeId, true);
                        notifyMoveProperParts(dualTree, currentUnionNode, pendingNodeId);
                        dualTree->moveProperParts(currentUnionNode, pendingNodeId);
                        disconnect(dualTree, pendingNodeId, true);
                    }
                    nodesPendingRemoval_.clear();
                } else if (nodeId != currentUnionNode) {
                    if (currentUnionNode == InvalidNode) {
                        if (removedMarks_.isMarked(nodeId)) {
                            nodesPendingRemoval_.push_back(nodeId);
                        }
                    } else {
                        reattachOutsideIntervalChildren(dualTree, currentUnionNode, nodeId, true);
                        notifyMoveProperParts(dualTree, currentUnionNode, nodeId);
                        dualTree->moveProperParts(currentUnionNode, nodeId);
                        disconnect(dualTree, nodeId, true);
                    }
                }
            }

            if (currentUnionNode == InvalidNode) {
                for (auto nodeId : nodesPendingRemoval_) {
                    if (!dualTree->isAlive(nodeId) || dualTree->isRoot(nodeId)) {
                        continue;
                    }
                    mergeNodeIntoParentKeepingCanonicalForm(dualTree, dualTree->getNodeParent(nodeId), nodeId);
                    disconnect(dualTree, nodeId, true);
                }
                currentMergeLevel = mergeNodesByLevel_.nextMergeLevel();
                currentUnionNode = previousLevelUnionNode;
            } else {
                if (currentMergeLevel == b) {
                    moveSelectedProperPartsToNode(dualTree, currentUnionNode, properPartSetC_);

                    for (auto nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                        attachChildOrMergeSameLevel(dualTree, currentUnionNode, nodeId);
                    }
                }

                if (previousLevelUnionNode != InvalidNode &&
                    dualTree->isAlive(previousLevelUnionNode) &&
                    !dualTree->hasChild(currentUnionNode, previousLevelUnionNode)) {
                    attachChildOrMergeSameLevel(dualTree, currentUnionNode, previousLevelUnionNode);
                }

                markAttributeUpdate(currentUnionNode);
                computeAttributeOnTreeNode(dualTree, currentUnionNode);

                previousLevelUnionNode = currentUnionNode;
                currentMergeLevel = mergeNodesByLevel_.nextMergeLevel();
            }
        }

        finalizeUpdateTreeAndContractRemovedNodes(dualTree, nodeCa, previousLevelUnionNode, true);
        (void) syntheticTargetNode;
        assertAllAliveNodesHaveProperParts(dualTree);
    }

    /**
     * @brief Initializes and applies one explicit connected-region update.
     * @details The input pixels form set `C` directly. This path is used by
     * `mergeNodeIntoParentAndUpdateDual`, where the primal node's direct proper
     * parts are processed component by component.
     */
    void applyExplicitRegionUpdate(DynamicComponentTree *dualTree, std::span<const PixelId> pixelsC, PixelType b, bool allowTargetNodeCreation, const char *apiName) {
        const std::string api = apiName == nullptr ? "explicit region update" : std::string(apiName);
        if (dualTree == nullptr) {
            throw std::invalid_argument(api + " requires a valid dual tree");
        }
        if (pixelsC.empty()) {
            throw std::invalid_argument(api + " requires a non-empty connected region");
        }
        if (dualTree != mintree_ && dualTree != maxtree_) {
            throw std::invalid_argument(api + " requires one of the adjuster's paired trees");
        }
        if (graph_ == nullptr) {
            throw std::runtime_error(api + " requires a valid adjacency relation");
        }

        const bool isMaxtree = dualTree == maxtree_;
        initializeUpdateStepState(isMaxtree);

        NodeId nodeCa = InvalidNode;
        for (PixelId pixelId : pixelsC) {
            if (pixelId < 0 || pixelId >= dualTree->getNumTotalProperParts()) {
                throw std::invalid_argument(api + " received an out-of-bounds pixel");
            }
            if (pixelsInCMarks_.isMarked(static_cast<std::size_t>(pixelId))) {
                continue;
            }

            properPartSetC_.push_back(pixelId);
            pixelsInCMarks_.mark(static_cast<std::size_t>(pixelId));

            const NodeId ownerNodeId = dualTree->getSmallestComponent(pixelId);
            if (ownerNodeId == InvalidNode || !dualTree->isAlive(ownerNodeId)) {
                throw std::runtime_error(api + " cannot locate a live owner for the region");
            }
            const PixelType altitudeP = static_cast<PixelType>(dualTree->getAltitude(ownerNodeId));
            if (nodeCa == InvalidNode ||
                ((isMaxtree && altitudeP < altitudeCa) ||
                 (!isMaxtree && altitudeP > altitudeCa))) {
                altitudeCa = altitudeP;
                nodeCa = ownerNodeId;
            }
        }

        if (properPartSetC_.empty() || nodeCa == InvalidNode) {
            throw std::invalid_argument(api + " requires a non-empty connected region");
        }

        runExplicitRegionUpdate(dualTree, b, isMaxtree, nodeCa, allowTargetNodeCreation);
    }

    /**
     * @brief Builds set `C` from a primal subtree and locates the corresponding dual anchor.
     * @details The subtree must be connected to its parent by a strict level
     * edge. The returned `nodeCa` is the extremal live dual owner among pixels in
     * `C`, and `b` is set to the primal parent altitude.
     */
    NodeId collectSubtreeRegionCForDualUpdate(DynamicComponentTree *primalTree, DynamicComponentTree *dualTree, NodeId subtreeRoot, PixelType &b, bool dualIsMaxtree) {
        const NodeId subtreeParentId = primalTree->getNodeParent(subtreeRoot);
        if (subtreeParentId == InvalidNode ||
            subtreeParentId == subtreeRoot ||
            !primalTree->isAlive(subtreeParentId)) {
            throw std::runtime_error("pruneSubtreeAndUpdateDual cannot locate the subtree parent");
        }

        b = static_cast<PixelType>(primalTree->getAltitude(subtreeParentId));
        const PixelType subtreeRootLevel = static_cast<PixelType>(primalTree->getAltitude(subtreeRoot));
        const bool strictPrimalEdge = primalTree->isMaxtree() ? (b < subtreeRootLevel) : (b > subtreeRootLevel);
        if (!strictPrimalEdge) {
            throw std::runtime_error(primalTree->isMaxtree()
                ? "pruneSubtreeAndUpdateDual requires a strict primal max-tree parent/subtree edge"
                : "pruneSubtreeAndUpdateDual requires a strict primal min-tree parent/subtree edge");
        }

        NodeId nodeCa = InvalidNode;
        for (auto subtreeNodeId : primalTree->getNodeSubtree(subtreeRoot)) {
            for (auto pixelId : primalTree->getProperParts(subtreeNodeId)) {
                properPartSetC_.push_back(pixelId);
                pixelsInCMarks_.mark(static_cast<std::size_t>(pixelId));

                const NodeId ownerNodeId = dualTree->getSmallestComponent(pixelId);
                if (ownerNodeId == InvalidNode) {
                    throw std::runtime_error("pruneSubtreeAndUpdateDual found a pixel without a dual owner");
                }
                assert(dualTree->isAlive(ownerNodeId));

                const PixelType altitudeP = static_cast<PixelType>(dualTree->getAltitude(ownerNodeId));
                if (nodeCa == InvalidNode ||
                    ((dualIsMaxtree && altitudeP < altitudeCa) ||
                     (!dualIsMaxtree && altitudeP > altitudeCa))) {
                    altitudeCa = altitudeP;
                    nodeCa = ownerNodeId;
                }
            }
        }

        if (properPartSetC_.empty() || nodeCa == InvalidNode) {
            throw std::runtime_error("pruneSubtreeAndUpdateDual received an empty region C");
        }

        return nodeCa;
    }

    /**
     * @brief Visits connected components induced by the direct proper parts of one node.
     * @details Components are discovered with the shared image adjacency and
     * emitted immediately through `onComponent`; no list of all components is
     * materialized.
     * @return Number of components emitted to the callback.
     */
    template<class ComponentCallback>
    std::size_t forEachConnectedDirectProperPartComponent(const DynamicComponentTree *tree, NodeId nodeId, ComponentCallback &&onComponent) const {
        if (tree == nullptr || graph_ == nullptr) {
            throw std::invalid_argument("mergeNodeIntoParentAndUpdateDual requires valid tree and adjacency relation");
        }

        const int numPixels = tree->getNumTotalProperParts();
        directComponentInRegionMarks_.ensureSize(static_cast<std::size_t>(numPixels));
        directComponentVisitedMarks_.ensureSize(static_cast<std::size_t>(numPixels));
        directComponentInRegionMarks_.resetAll();
        directComponentVisitedMarks_.resetAll();
        directComponentSeeds_.clear();
        directComponentSeeds_.reserve(static_cast<std::size_t>(tree->getNumProperParts(nodeId)));
        directComponent_.reserve(static_cast<std::size_t>(tree->getNumProperParts(nodeId)));

        for (PixelId pixelId : tree->getProperParts(nodeId)) {
            if (pixelId < 0 || pixelId >= numPixels) {
                throw std::runtime_error("mergeNodeIntoParentAndUpdateDual found an out-of-bounds proper part");
            }
            const auto pixelIndex = static_cast<std::size_t>(pixelId);
            if (directComponentInRegionMarks_.isMarked(pixelIndex)) {
                continue;
            }
            directComponentInRegionMarks_.mark(pixelIndex);
            directComponentSeeds_.push_back(pixelId);
        }

        std::size_t componentCount = 0;
        directComponentStack_.clear();
        for (PixelId seed : directComponentSeeds_) {
            const auto seedIndex = static_cast<std::size_t>(seed);
            if (directComponentVisitedMarks_.isMarked(seedIndex)) {
                continue;
            }

            directComponent_.clear();
            directComponentStack_.clear();
            directComponentStack_.push_back(seed);
            directComponentVisitedMarks_.mark(seedIndex);

            while (!directComponentStack_.empty()) {
                const PixelId pixelId = directComponentStack_.back();
                directComponentStack_.pop_back();
                directComponent_.push_back(pixelId);

                for (PixelId neighbor : graph_->getNeighborPixels(pixelId)) {
                    if (neighbor < 0 || neighbor >= numPixels) {
                        continue;
                    }
                    const auto neighborIndex = static_cast<std::size_t>(neighbor);
                    if (!directComponentInRegionMarks_.isMarked(neighborIndex) ||
                        directComponentVisitedMarks_.isMarked(neighborIndex)) {
                        continue;
                    }
                    directComponentVisitedMarks_.mark(neighborIndex);
                    directComponentStack_.push_back(neighbor);
                }
            }

            onComponent(std::span<const PixelId>(directComponent_.data(), directComponent_.size()));
            ++componentCount;
        }

        return componentCount;
    }

    /**
     * @brief Validates that a requested update uses this adjuster's primal/dual pair.
     */
    void requirePairedTrees(DynamicComponentTree *primalTree, DynamicComponentTree *dualTree, const char *apiName) const {
        if (primalTree == nullptr || dualTree == nullptr) {
            const char *api = apiName == nullptr ? "operation" : apiName;
            throw std::invalid_argument(std::string(api) + " requires valid primal and dual trees");
        }
        if (!((primalTree == maxtree_ && dualTree == mintree_) ||
              (primalTree == mintree_ && dualTree == maxtree_))) {
            const char *api = apiName == nullptr ? "operation" : apiName;
            throw std::invalid_argument(std::string(api) + " requires the adjuster's paired primal/dual trees");
        }
    }

    /**
     * @brief Validates that a public update target is a live non-root node.
     */
    void requireValidNonRootNode(const DynamicComponentTree *tree, NodeId nodeId, const char *apiName) const {
        if (tree == nullptr ||
            nodeId == InvalidNode ||
            nodeId < 0 ||
            nodeId >= tree->getNumInternalNodeSlots() ||
            !tree->isNode(nodeId) ||
            !tree->isAlive(nodeId) ||
            tree->isRoot(nodeId)) {
            const char *api = apiName == nullptr ? "operation" : apiName;
            throw std::invalid_argument(std::string(api) + " requires a live non-root node");
        }
    }

public:
    /**
     * @brief Creates the adjuster from a fixed dynamic min-tree / max-tree pair.
     * @details The pair passed here defines the lifetime association of the
     * adjuster. Subsequent prune/update calls operate on this same pair
     * and do not accept per-call tree overrides.
     * @param mintree Pointer to the dynamic min-tree owned externally.
     * @param maxtree Pointer to the dynamic max-tree owned externally.
     * @param graph Adjacency relation shared by both trees.
     */
    DualMinMaxTreeIncrementalFilter(DynamicComponentTree *mintree, DynamicComponentTree *maxtree, AdjacencyRelation &graph)
        : mintree_(mintree), maxtree_(maxtree), graph_(&graph),
          mergeNodesByLevel_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)),
          removedMarks_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)),
          pixelsInCMarks_(std::max(mintree ? mintree->getNumTotalProperParts() : 0, maxtree ? maxtree->getNumTotalProperParts() : 0)),
          climbedNodeMarks_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)),
          attributeUpdateMarks_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        assert(graph_ != nullptr);
    }

    /**
     * @brief Returns the frontier roots stored in the latest `frontierNodesAboveB` collection.
     * @return Mutable reference to the frontier-root list.
     */
    std::vector<NodeId> &getFrontierNodesAboveB() {
        return mergeNodesByLevel_.getFrontierNodesAboveB();
    }

    /**
     * @brief Returns the merge bucket associated with an altitude value.
     * @param level Queried level.
     * @return Mutable reference to the bucket for that level.
     */
    std::vector<NodeId> &getMergedNodes(const PixelType &level) {
        return mergeNodesByLevel_.getMergedNodes(level);
    }

    /**
     * @brief Registers the incremental attribute computers and their external buffers.
     * @details Buffer sizing and lifetime are the caller's responsibility. Each
     * buffer must cover the full global id space of the corresponding tree.
     * The registered computers and buffers must correspond to the same fixed
     * min-tree / max-tree pair provided at construction time.
     * @param computerMin Incremental attribute computer for the min-tree.
     * @param computerMax Incremental attribute computer for the max-tree.
     * @param bufferMin External buffer associated with the min-tree.
     * @param bufferMax External buffer associated with the max-tree.
     */
    void setAttributeComputer(DynamicAttributeComputer &computerMin, DynamicAttributeComputer &computerMax, std::span<float> bufferMin, std::span<float> bufferMax) {
        attrComputerMin_ = &computerMin;
        attrComputerMax_ = &computerMax;
        bufferMin_ = bufferMin;
        bufferMax_ = bufferMax;
    }

    /**
     * @brief Prunes one subtree from a primal tree and updates the paired dual tree.
     * @details This is the clean public path for the explicit prune/update
     * engine. It does not expose diagnostics and reports invalid input or
     * rejected events through exceptions.
     */
    void pruneSubtreeAndUpdateDual(DynamicComponentTree *primalTree, DynamicComponentTree *dualTree, NodeId subtreeRoot) {
        requirePairedTrees(primalTree, dualTree, "pruneSubtreeAndUpdateDual");
        requireValidNonRootNode(primalTree, subtreeRoot, "pruneSubtreeAndUpdateDual");

        const bool isMaxtree = dualTree == maxtree_;
        initializeUpdateStepState(isMaxtree);
        PixelType b{};
        const NodeId nodeCa = collectSubtreeRegionCForDualUpdate(
            primalTree,
            dualTree,
            subtreeRoot,
            b,
            isMaxtree);
        runExplicitRegionUpdate(dualTree, b, isMaxtree, nodeCa, false);
        primalTree->pruneNode(subtreeRoot);
    }

    /**
     * @brief Prunes a batch of subtree roots from a primal tree and updates the paired dual tree.
     * @details Each root is processed as an independent explicit update step.
     * Invalid roots throw instead of being silently skipped.
     */
    void pruneSubtreesAndUpdateDual(DynamicComponentTree *primalTree, DynamicComponentTree *dualTree, std::span<const NodeId> nodesToPrune) {
        requirePairedTrees(primalTree, dualTree, "pruneSubtreesAndUpdateDual");
        const bool isMaxtree = dualTree == maxtree_;
        for (NodeId rootSubtree : nodesToPrune) {
            requireValidNonRootNode(primalTree, rootSubtree, "pruneSubtreesAndUpdateDual");
            initializeUpdateStepState(isMaxtree);
            PixelType b{};
            const NodeId nodeCa = collectSubtreeRegionCForDualUpdate(
                primalTree,
                dualTree,
                rootSubtree,
                b,
                isMaxtree);
            runExplicitRegionUpdate(dualTree, b, isMaxtree, nodeCa, false);
            primalTree->pruneNode(rootSubtree);
        }
    }

    /**
     * @brief Merges a primal node into its parent and updates the paired dual tree.
     * @details The direct proper parts of `nodeId` are split into connected
     * components and each component is updated through the explicit region
     * engine. The primal merge is applied only after the dual updates succeed.
     */
    void mergeNodeIntoParentAndUpdateDual(DynamicComponentTree *primalTree, DynamicComponentTree *dualTree, NodeId nodeId) {
        requirePairedTrees(primalTree, dualTree, "mergeNodeIntoParentAndUpdateDual");
        requireValidNonRootNode(primalTree, nodeId, "mergeNodeIntoParentAndUpdateDual");

        const NodeId parentId = primalTree->getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !primalTree->isAlive(parentId)) {
            throw std::runtime_error("mergeNodeIntoParentAndUpdateDual cannot locate the primal parent");
        }

        const PixelType b = static_cast<PixelType>(primalTree->getAltitude(parentId));
        const std::size_t componentCount = forEachConnectedDirectProperPartComponent(
            primalTree,
            nodeId,
            [&](std::span<const PixelId> component) {
                applyExplicitRegionUpdate(
                    dualTree,
                    component,
                    b,
                    true,
                    "mergeNodeIntoParentAndUpdateDual");
            });
        if (componentCount == 0) {
            throw std::runtime_error("mergeNodeIntoParentAndUpdateDual requires a non-empty direct proper-part region");
        }

        mergePrimalNodeIntoParentAndUpdateAttribute(primalTree, nodeId, parentId);
    }

    /**
     * @brief Builds the merged and nested collections of one subtree-adjustment step.
     * @details The routine computes the sets equivalent to `mergeNodesByLevel`
     * and `frontierNodesAboveB` from the nodes adjacent to set `C` in the tree
     * being updated.
     *
     * The flow follows Algorithm 2 from the paper: adjacent nodes serve as
     * seeds and each relevant ancestor is visited at most once through a
     * marking mechanism. Each visited node enters only the bucket of its own
     * level or, if it is above `b`, into `frontierNodesAboveB`.
     *
     * @param dualTree Target tree updated in place.
     * @param properPartSetC Set `C` formed by the proper parts of the subtree removed in the primal tree.
     * @param b Altitude of the parent of the removed subtree in the primal tree.
     * @param isMaxtree `true` for the max-tree direction and `false` for the min-tree direction.
     */
    void buildMergedAndNestedCollections(DynamicComponentTree *dualTree, std::vector<PixelId> &properPartSetC, PixelType b, bool isMaxtree) {
        assert(dualTree != nullptr);
        assert(!properPartSetC.empty());

        // Phase 1: step collection already arrives clean thanks to `initializeUpdateStepState`.
        assert(graph_ != nullptr);

        // Phase 2: for each pixel of C, locate valid adjacent seeds and climb
        // each relevant ancestral path only once.
        for (PixelId p : properPartSetC) {
            for (PixelId q : graph_->getNeighborPixels(p)) {
                if (pixelsInCMarks_.isMarked(static_cast<std::size_t>(q))) {
                    continue; // Neighbors internal to C do not generate adjacent seeds.
                }

                const NodeId nodeQ = dualTree->getSmallestComponent(q);
                if (nodeQ == InvalidNode) {
                    continue; // Pixels without a corresponding live component do not enter the collection.
                }

                const PixelType altitudeQ = static_cast<PixelType>(dualTree->getAltitude(nodeQ));
                const bool validSeed = (isMaxtree && altitudeQ >= altitudeCa) || (!isMaxtree && altitudeQ <= altitudeCa);
                if (!validSeed) {
                    continue; // Seeds outside the valid interval do not participate in the merge.
                }

                if (mergeNodesByLevel_.markAdjacentSeed(nodeQ)) {
                    NodeId nodeSubtree = nodeQ;
                    NodeId n = nodeQ;
                    while (n != InvalidNode && dualTree->isAlive(n) && !climbedNodeMarks_.isMarked(static_cast<std::size_t>(n))) {
                        const PixelType levelCurrent = static_cast<PixelType>(dualTree->getAltitude(n));
                        if (!((isMaxtree && levelCurrent >= altitudeCa) || (!isMaxtree && levelCurrent <= altitudeCa))) {
                            break; // The climb left the level interval relevant to C.
                        }

                        climbedNodeMarks_.mark(static_cast<std::size_t>(n));
                        nodeSubtree = n;

                        if ((isMaxtree && levelCurrent <= b) || (!isMaxtree && levelCurrent >= b)) {
                            mergeNodesByLevel_.addMergeNode(*dualTree, nodeSubtree);
                        } else {
                            NodeId parentId = dualTree->getNodeParent(nodeSubtree);
                            if (parentId == nodeSubtree) {
                                parentId = InvalidNode;
                            }
                            if (!(parentId != InvalidNode && ((isMaxtree && dualTree->getAltitude(parentId) > b) || (!isMaxtree && dualTree->getAltitude(parentId) < b)))) {
                                mergeNodesByLevel_.addFrontierNodeAboveB(nodeSubtree);
                            }
                        }

                        const NodeId parentId = dualTree->getNodeParent(n);
                        if (parentId == n) {
                            break; // Reached the structural top of this path.
                        }
                        n = parentId;
                    }
                }
            }
        }
    }



    /**
     * @brief Updates the dual tree after removing a subtree from the primal tree.
     * @details The method executes the full incremental step:
     * 1. collects in `properPartSetC_` the proper parts of the subtree still alive in the primal tree;
     * 2. locates `nodeCa` and the altitude interval to traverse;
     * 3. builds `mergeNodesByLevel_` and `frontierNodesAboveB`;
     * 4. merges nodes level by level and moves set `C` when the sweep reaches `b`;
     * 5. reconnects the final union node in the updated hierarchy.
     *
     * The method assumes that `subtreeRoot` is still alive in the primal tree
     * at call time. Actual pruning in the primal occurs only after the dual has
     * been adjusted, in the public `prune*AndUpdate*` wrappers.
     *
     * @param dualTree Dual hierarchy updated in place.
     * @param subtreeRoot Root of the subtree to be removed in the primal tree.
     */
    void updateTree(DynamicComponentTree *dualTree, NodeId subtreeRoot) {
        assert(dualTree != nullptr);
        assert(subtreeRoot != InvalidNode);
        const bool isMaxtree = dualTree == maxtree_;
        initializeUpdateStepState(isMaxtree);
        DynamicComponentTree *primalTree = isMaxtree ? mintree_ : maxtree_;

        assert(primalTree != nullptr);
        assert(subtreeRoot >= 0 && subtreeRoot < primalTree->getNumInternalNodeSlots());
        assert(primalTree->isNode(subtreeRoot) && primalTree->isAlive(subtreeRoot));

        const NodeId subtreeParentId = primalTree->getNodeParent(subtreeRoot);
        assert(subtreeParentId != InvalidNode && subtreeParentId != subtreeRoot);
        const PixelType b = static_cast<PixelType>(primalTree->getAltitude(subtreeParentId));

        NodeId nodeCa = InvalidNode;
        // Phase 1: collect set C in the primal tree and locate `nodeCa`
        // as the extremal representative of C in the dual tree.
        for (auto subtreeNodeId : primalTree->getNodeSubtree(subtreeRoot)) {
            // Set `C` is formed by the proper parts of all nodes in the subtree removed from the primal tree.
            for (auto p : primalTree->getProperParts(subtreeNodeId)) {
                properPartSetC_.push_back(p);
                pixelsInCMarks_.mark(static_cast<std::size_t>(p));

                const NodeId ownerNodeId = dualTree->getSmallestComponent(p);
                if (ownerNodeId == InvalidNode) {
                    continue; // Pixels without a live dual component do not contribute to nodeCa.
                }
                const PixelType altitudeP = static_cast<PixelType>(dualTree->getAltitude(ownerNodeId));
                if (nodeCa == InvalidNode || ((isMaxtree && altitudeP < altitudeCa) || (!isMaxtree && altitudeP > altitudeCa))) {
                    altitudeCa = altitudeP;
                    nodeCa = ownerNodeId;
                }
            }
        }

        if (properPartSetC_.empty()) {
            return;
        }
        assert(nodeCa != InvalidNode);
        assert(!properPartSetC_.empty());

        // Phase 2: build the per-level merge buckets and the frontier roots above b.
        buildMergedAndNestedCollections(dualTree, properPartSetC_, b, isMaxtree);

        // Phase 3: level sweep between b and a, merging active buckets and
        // propagating the union node built at each level.
        PixelType currentMergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId currentUnionNode = InvalidNode;
        NodeId previousLevelUnionNode = InvalidNode;
        nodesPendingRemoval_.reserve(mergeNodesByLevel_.getMaxBucketSize());
        while (mergeNodesByLevel_.hasMergeLevel() && ((isMaxtree && currentMergeLevel > altitudeCa) || (!isMaxtree && currentMergeLevel < altitudeCa))) {
            auto &nodesAtCurrentLevel = mergeNodesByLevel_.getMergedNodes((int) currentMergeLevel);
            currentUnionNode = InvalidNode;
            nodesPendingRemoval_.clear();

            for (auto nodeId : nodesAtCurrentLevel) {
                if (!dualTree->isAlive(nodeId)) {
                    continue; // Nodes already removed before this level do not need to be revisited.
                } else if (currentUnionNode == InvalidNode && !removedMarks_.isMarked(nodeId)) {
                    currentUnionNode = nodeId;

                    if (!dualTree->isRoot(currentUnionNode)) {
                        disconnect(dualTree, currentUnionNode, false);
                    }

                    for (auto pendingNodeId : nodesPendingRemoval_) {
                        reattachOutsideIntervalChildren(dualTree, currentUnionNode, pendingNodeId);
                        notifyMoveProperParts(dualTree, currentUnionNode, pendingNodeId);
                        dualTree->moveProperParts(currentUnionNode, pendingNodeId);
                        disconnect(dualTree, pendingNodeId, true);
                    }
                    nodesPendingRemoval_.clear();
                } else if (nodeId != currentUnionNode) {
                    if (currentUnionNode == InvalidNode) {
                        // If the first candidate at this level was already
                        // emptied while moving proper parts, it stays pending
                        // until the effective union node for this level appears.
                        if (removedMarks_.isMarked(nodeId)) {
                            nodesPendingRemoval_.push_back(nodeId);
                        }
                    } else {
                        reattachOutsideIntervalChildren(dualTree, currentUnionNode, nodeId);
                        notifyMoveProperParts(dualTree, currentUnionNode, nodeId);
                        dualTree->moveProperParts(currentUnionNode, nodeId);
                        disconnect(dualTree, nodeId, true);
                    }
                }
            }

            if (currentUnionNode == InvalidNode) {
                // Rare case: all nodes at this level were already marked for
                // removal. In that case, each one is absorbed by its current parent.
                for (auto nodeId : nodesPendingRemoval_) {
                    if (!dualTree->isAlive(nodeId) || dualTree->isRoot(nodeId)) {
                        continue; // Dead seeds or the root can no longer be absorbed through this shortcut.
                    }
                    mergedParentAndChildren(dualTree, dualTree->getNodeParent(nodeId), nodeId);
                    disconnect(dualTree, nodeId, true);
                }
                currentMergeLevel = mergeNodesByLevel_.nextMergeLevel();
                currentUnionNode = previousLevelUnionNode;
            } else {
                if (currentMergeLevel == b) {
                    // Level `b` is the point where set `C` effectively enters
                    // the union node and where frontiers above `b` are
                    // reattached under that same node.
                    moveSelectedProperPartsToNode(dualTree, currentUnionNode, properPartSetC_);

                    for (auto nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                        disconnect(dualTree, nodeId, false);
                        dualTree->attachNode(currentUnionNode, nodeId);
                    }

                }

                if (previousLevelUnionNode != InvalidNode && dualTree->isAlive(previousLevelUnionNode) && !dualTree->hasChild(currentUnionNode, previousLevelUnionNode)) {
                    if (!dualTree->isRoot(previousLevelUnionNode)) {
                        disconnect(dualTree, previousLevelUnionNode, false);
                    }
                    dualTree->attachNode(currentUnionNode, previousLevelUnionNode);
                }

                markAttributeUpdate(currentUnionNode);
                computeAttributeOnTreeNode(dualTree, currentUnionNode);

                previousLevelUnionNode = currentUnionNode;
                currentMergeLevel = mergeNodesByLevel_.nextMergeLevel();
            }
        }

        // Phase 4: position the final union node and contract the removed
        // nodes that became temporarily empty during the sweep.
        finalizeUpdateTreeAndContractRemovedNodes(dualTree, nodeCa, previousLevelUnionNode);
        assertAllAliveNodesHaveProperParts(dualTree);
    }

    /**
     * @brief Propagates max-tree prunes to the dual min-tree of this adjuster.
     * @details For each valid root in `nodesToPrune`, the wrapper first calls
     * `updateTree(mintree_, rootSubtree)` and only then executes
     * `maxtree_->pruneNode(rootSubtree)`. This order is mandatory because the
     * incremental update still needs to inspect the live subtree in the primal
     * tree to build `C`, the frontiers, and the corresponding merges.
     *
     * The cost is the sum of `updateTree` and `pruneNode` for each valid root.
     *
     * @param nodesToPrune Roots of the subtrees removed from the max-tree.
     */
    void pruneMaxTreeAndUpdateMinTree(std::vector<NodeId> &nodesToPrune) {
        assert(removedMarks_.stamp.size() >= static_cast<std::size_t>(std::max(mintree_ ? mintree_->getNumInternalNodeSlots() : 0, maxtree_ ? maxtree_->getNumInternalNodeSlots() : 0)));
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == maxtree_->getRoot() || !maxtree_->isNode(rootSubtree) || !maxtree_->isAlive(rootSubtree)) {
                continue; // Ignore invalid roots, the global root, and nodes already removed.
            }
            updateTree(mintree_, rootSubtree);
            maxtree_->pruneNode(rootSubtree);
        }
    }

    /**
     * @brief Propagates min-tree prunes to the dual max-tree of this adjuster.
     * @details This is the min->max counterpart of the previous wrapper. For
     * each valid root in `nodesToPrune`, the method first executes
     * `updateTree(maxtree_, rootSubtree)` and then
     * `mintree_->pruneNode(rootSubtree)`. The order remains the same: the dual
     * tree is updated while the primal-tree subtree is still alive.
     *
     * The cost is the sum of `updateTree` and `pruneNode` for each valid root.
     *
     * @param nodesToPrune Roots of the subtrees removed from the min-tree.
     */
    void pruneMinTreeAndUpdateMaxTree(std::vector<NodeId> &nodesToPrune) {
        assert(removedMarks_.stamp.size() >= static_cast<std::size_t>(std::max(mintree_ ? mintree_->getNumInternalNodeSlots() : 0, maxtree_ ? maxtree_->getNumInternalNodeSlots() : 0)));
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == mintree_->getRoot() || !mintree_->isNode(rootSubtree) || !mintree_->isAlive(rootSubtree)) {
                continue; // Ignore invalid roots, the global root, and nodes already removed.
            }
            updateTree(maxtree_, rootSubtree);
            mintree_->pruneNode(rootSubtree);
        }
    }

    /**
     * @brief Returns the min-tree currently associated with the adjuster.
     */
    DynamicComponentTree *getMinTree() const { return mintree_; }

    /**
     * @brief Returns the max-tree currently associated with the adjuster.
     */
    DynamicComponentTree *getMaxTree() const { return maxtree_; }
};
