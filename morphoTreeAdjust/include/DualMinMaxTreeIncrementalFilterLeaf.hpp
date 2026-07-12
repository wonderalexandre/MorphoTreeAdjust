#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "DynamicComponentTree.hpp"

/**
 * @brief Leaf-by-leaf incremental updater for the dynamic min-tree / max-tree pair.
 *
 * This class is the `leaf` variant of `DualMinMaxTreeIncrementalFilter`. The
 * main difference is not in the final CASF semantics, but in the granularity
 * of the update:
 *
 * - here, the basic step is the removal of one live leaf;
 * - the dual tree is updated after each removed leaf;
 * - the full sequence of leaves removed at one threshold may be long.
 *
 * Consequences of this granularity:
 *
 * - more events are produced over the course of CASF;
 * - each event tends to have smaller support;
 * - direct CASF filtering is usually more expensive than in the `subtree` variant;
 * - the resulting intermediate hierarchy usually has more nodes.
 *
 * Structures maintained by the class:
 *
 * - `mintree_`, `maxtree_`: mutable dynamic trees;
 * - incremental attribute computers and external buffers;
 * - `mergeNodesByLevel_`: temporary buckets and collections for the current step.
 *
 * Algorithmic intuition:
 *
 * - while a live primal leaf is still present, its proper parts define set `C`;
 * - the relevant neighbors of `C` are collected;
 * - the dual tree is updated by a level sweep, analogously to the `subtree`
 *   adjuster, but usually over a smaller region.
 * - only after the dual update succeeds is the leaf pruned from the primal tree.
 *
 * Operational invariants:
 *
 * - before and after each step, both trees represent the same image;
 * - the removed leaf is non-root and live at the time of the operation;
 * - the dual tree is adjusted locally without a global rebuild;
 * - incremental attributes remain consistent when configured.
 *
 * Complexity:
 *
 * Use:
 *
 * - `P_L`: number of proper parts in the removed leaf;
 * - `E_L`: total adjacency-stencil positions examined for those pixels
 *   (`E_L <= P_L * D` for maximum stencil size `D`);
 * - `A_L`: number of adjacent nodes collected in the step;
 * - `M_L`: number of nodes placed in the merge buckets;
 * - `K_L`: number of nodes effectively edited in the dual tree.
 *
 * Let `L = max(PixelType) + 1`, the number of representable altitude levels.
 * `PixelType` must be integral and have a lowest value of zero because this
 * implementation allocates one dense merge bucket per level. An individual
 * step costs approximately:
 *
 * - `O(L + P_L + E_L + A_L + M_L + K_L)` time, because the active-level scan
 *   traverses all `L` buckets;
 * - `O(L + N + P)` persistent workspace, where `N` is the paired trees'
 *   node-id capacity and `P` is the full image domain used by pixel-indexed
 *   generation stamps, plus the collections for the current affected region.
 *
 * The practical difference relative to `DualMinMaxTreeIncrementalFilter` lies in
 * the number of steps: the total cost of a `leaf` threshold can be the sum of
 * many small steps, whereas `subtree` aggregates many of those changes into
 * fewer, larger steps.
 *
 * In terms of project usage:
 *
 * - this class is the basis of per-pixel `leaf` CASF;
 * - its advantage is fine-grained granularity;
 * - its cost is the large number of dual updates.
 *
 * The adjuster owns none of its trees, adjacency, registered attribute
 * computers, or attribute-buffer storage. All borrowed objects and allocations
 * must remain alive and at stable addresses. A reallocation of a span's backing
 * container, or tree growth beyond a buffer's current coverage, requires
 * resizing and a new `setAttributeComputer()` call before the next update.
 *
 * @tparam PixelType Integral altitude type whose minimum is zero. Wide types
 * are generally impractical because storage is proportional to
 * `max(PixelType) + 1`.
 */
template<typename PixelType = AltitudeType>
class DualMinMaxTreeIncrementalFilterLeaf {
public:
    class MergedNodesCollection {
    private:
        static_assert(std::numeric_limits<PixelType>::is_integer, "MergedNodesCollection requires integral PixelType");
        static_assert(std::numeric_limits<PixelType>::lowest() == 0, "MergedNodesCollection expects non-negative contiguous levels");
        static constexpr std::size_t NumLevels = static_cast<std::size_t>(std::numeric_limits<PixelType>::max()) + 1;

        std::array<std::vector<NodeId>, NumLevels> mergeNodesByLevelStorage_;
        std::vector<PixelType> mergeLevels_;
        std::vector<NodeId> adjacentNodes_;
        std::vector<NodeId> frontierNodesAboveB_;
        GenerationStampSet visited_;
        GenerationStampSet visitedAdj_;
        int currentMergeLevelIndex_ = 0;
        bool isMaxtree_ = false;

    public:
        explicit MergedNodesCollection(int maxNodes = 0)
            : visited_(std::max(0, maxNodes)),
              visitedAdj_(std::max(0, maxNodes)) {}

        void resetCollection(bool isMaxtree) {
            isMaxtree_ = isMaxtree;
            for (auto level : mergeLevels_) {
                mergeNodesByLevelStorage_[static_cast<std::size_t>(level)].clear();
            }
            mergeLevels_.clear();
            adjacentNodes_.clear();
            frontierNodesAboveB_.clear();
            visited_.resetAll();
            visitedAdj_.resetAll();
            currentMergeLevelIndex_ = 0;
        }

        std::vector<NodeId> &getMergedNodes(const PixelType &level) {
            return mergeNodesByLevelStorage_[static_cast<std::size_t>(level)];
        }

        std::vector<NodeId> &getAdjacentNodes() {
            return adjacentNodes_;
        }

        std::vector<NodeId> &getFrontierNodesAboveB() {
            return frontierNodesAboveB_;
        }

        void addFrontierNodeAboveB(NodeId nodeId) {
            if (!visited_.isMarked(nodeId)) {
                frontierNodesAboveB_.push_back(nodeId);
                visited_.mark(nodeId);
            }
        }

        void addNodesOfPath(const DynamicComponentTree &tree, NodeId adjacentNode, NodeId nodeCa) {
            if (visited_.isMarked(adjacentNode)) {
                return;
            }

            NodeId nodeId = adjacentNode;
            while (nodeId != InvalidNode) {
                if (!visited_.isMarked(nodeId)) {
                    auto &bucket = getMergedNodes(static_cast<PixelType>(tree.getAltitude(nodeId)));
                    bucket.push_back(nodeId);
                    visited_.mark(nodeId);
                } else {
                    break;
                }

                if (nodeId == nodeCa) {
                    break;
                }

                const NodeId parentId = tree.getNodeParent(nodeId);
                nodeId = (parentId == nodeId) ? InvalidNode : parentId;
            }
        }

        void addMergeNode(const DynamicComponentTree &tree, NodeId nodeId) {
            if (visited_.isMarked(nodeId)) {
                return;
            }

            auto &bucket = getMergedNodes(static_cast<PixelType>(tree.getAltitude(nodeId)));
            bucket.push_back(nodeId);
            visited_.mark(nodeId);
        }

        bool markAdjacentSeed(NodeId nodeId) {
            if (visitedAdj_.isMarked(nodeId)) {
                return false;
            }

            adjacentNodes_.push_back(nodeId);
            visitedAdj_.mark(nodeId);
            return true;
        }

        PixelType firstMergeLevel() {
            mergeLevels_.clear();
            for (std::size_t i = 0; i < mergeNodesByLevelStorage_.size(); ++i) {
                if (!mergeNodesByLevelStorage_[i].empty()) {
                    mergeLevels_.push_back(static_cast<PixelType>(i));
                }
            }

            if (mergeLevels_.empty()) {
                return PixelType{};
            }

            currentMergeLevelIndex_ = isMaxtree_ ? static_cast<int>(mergeLevels_.size()) - 1 : 0;
            return mergeLevels_[currentMergeLevelIndex_];
        }

        bool hasMergeLevel() const {
            return !mergeLevels_.empty() && currentMergeLevelIndex_ >= 0 && currentMergeLevelIndex_ < static_cast<int>(mergeLevels_.size());
        }

        PixelType nextMergeLevel() {
            currentMergeLevelIndex_ = isMaxtree_ ? currentMergeLevelIndex_ - 1 : currentMergeLevelIndex_ + 1;
            if (!hasMergeLevel()) {
                return PixelType{};
            }
            return mergeLevels_[currentMergeLevelIndex_];
        }
    };

private:
    DynamicComponentTree *mintree_ = nullptr;
    DynamicComponentTree *maxtree_ = nullptr;
    AdjacencyRelation *graph_ = nullptr;
    DynamicAttributeComputer *attrComputerMin_ = nullptr;
    DynamicAttributeComputer *attrComputerMax_ = nullptr;
    std::span<float> bufferMin_;
    std::span<float> bufferMax_;
    MergedNodesCollection mergeNodesByLevel_;
    GenerationStampSet removedMarks_;
    std::vector<NodeId> removedNodesPendingAbsorption_;
    GenerationStampSet pixelsInLeafMarks_;
    GenerationStampSet climbedNodeMarks_;
    GenerationStampSet attributeUpdateMarks_;
    PixelType altitudeCa_ = PixelType{};

    /** @brief Returns the primal tree corresponding to the current adjustment direction. */
    DynamicComponentTree *getPrimalTree(bool isMaxtree) {
        return isMaxtree ? mintree_ : maxtree_;
    }

    /** @brief Selects the incremental attribute computer associated with the given tree. */
    DynamicAttributeComputer *getAttributeComputer(DynamicComponentTree *tree) {
        if (tree == mintree_) {
            return attrComputerMin_;
        }
        if (tree == maxtree_) {
            return attrComputerMax_;
        }
        return nullptr;
    }

    /** @brief Selects the external attribute buffer associated with the given tree. */
    std::span<float> getAttributeBuffer(DynamicComponentTree *tree) {
        if (tree == mintree_) {
            return bufferMin_;
        }
        if (tree == maxtree_) {
            return bufferMax_;
        }
        return {};
    }

    /**
     * @brief Incrementally recomputes the attribute of a node after a local edit.
     *
     * The method reapplies the `pre/merge/postProcessing` protocol only to the
     * touched node, assuming that its children are already consistent.
     */
    void computeAttributeOnTreeNode(DynamicComponentTree *tree, NodeId nodeId) {
        auto *computer = getAttributeComputer(tree);
        auto buffer = getAttributeBuffer(tree);
        if (tree == nullptr || computer == nullptr || buffer.empty() ||
            nodeId == InvalidNode || !tree->isNode(nodeId) || !tree->isAlive(nodeId)) {
            return;
        }

        if (!attributeUpdateMarks_.isMarked(static_cast<size_t>(nodeId))) {
            return;
        }

        const PixelType level = static_cast<PixelType>(tree->getAltitude(nodeId));
        if ((tree == maxtree_ && level < altitudeCa_) || (tree == mintree_ && level > altitudeCa_)) {
            attributeUpdateMarks_.unmark(static_cast<size_t>(nodeId));
            return;
        }

        computer->preProcessing(nodeId, buffer);
        for (NodeId childId : tree->getChildren(nodeId)) {
            computer->mergeProcessing(nodeId, childId, buffer);
        }
        computer->postProcessing(nodeId, buffer);
        attributeUpdateMarks_.unmark(static_cast<size_t>(nodeId));
    }

    DynamicAttributeComputer *getAttributeComputer(DynamicComponentTree *tree) const {
        if (tree == mintree_) {
            return attrComputerMin_;
        }
        if (tree == maxtree_) {
            return attrComputerMax_;
        }
        return nullptr;
    }

    void notifyNodeRemoved(DynamicComponentTree *tree, NodeId nodeId) {
        auto *computer = getAttributeComputer(tree);
        if (computer != nullptr && nodeId != InvalidNode) {
            computer->onNodeRemoved(nodeId);
        }
    }

    void notifyMoveProperParts(DynamicComponentTree *tree, NodeId targetNodeId, NodeId sourceNodeId) {
        auto *computer = getAttributeComputer(tree);
        if (computer != nullptr) {
            computer->onMoveProperParts(targetNodeId, sourceNodeId);
        }
    }

    void notifyMoveProperPart(DynamicComponentTree *tree, NodeId targetNodeId, NodeId sourceNodeId, PixelId pixelId) {
        auto *computer = getAttributeComputer(tree);
        if (computer != nullptr) {
            computer->onMoveProperPart(targetNodeId, sourceNodeId, pixelId);
        }
    }

    void markAttributeUpdate(NodeId nodeId) {
        if (nodeId != InvalidNode) {
            attributeUpdateMarks_.mark(static_cast<size_t>(nodeId));
        }
    }

    void initializeUpdateStepState(bool isMaxtree) {
        mergeNodesByLevel_.resetCollection(isMaxtree);
        removedMarks_.resetAll();
        removedNodesPendingAbsorption_.clear();
        pixelsInLeafMarks_.resetAll();
        climbedNodeMarks_.resetAll();
        attributeUpdateMarks_.resetAll();
        altitudeCa_ = PixelType{};
    }

    /** @brief Detaches a node from its current parent, optionally releasing its slot. */
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

    /** @brief Merges the structural contents of `childId` into `parentId`. */
    void mergedParentAndChildren(DynamicComponentTree *tree, NodeId parentId, NodeId childId) {
        assert(tree != nullptr);
        tree->moveChildren(parentId, childId);
        notifyMoveProperParts(tree, parentId, childId);
        tree->moveProperParts(parentId, childId);
    }

    bool canAbsorbRemovedNode(DynamicComponentTree *tree, NodeId nodeId) const {
        return tree != nullptr &&
               nodeId != InvalidNode &&
               tree->isNode(nodeId) &&
               tree->isAlive(nodeId) &&
               removedMarks_.isMarked(nodeId) &&
               tree->getNumProperParts(nodeId) == 0;
    }

    void absorbRemovedRootNode(DynamicComponentTree *tree, NodeId removedNodeId) {
        const NodeId firstChild = tree->getFirstChild(removedNodeId);
        if (firstChild == InvalidNode) {
            return;
        }

        NodeId newRoot = firstChild;
        bool newRootChanged = false;
        for (NodeId childId = firstChild; childId != InvalidNode; childId = tree->getNextSibling(childId)) {
            if ((tree->isMaxtree() && tree->getAltitude(childId) < tree->getAltitude(newRoot)) ||
                (!tree->isMaxtree() && tree->getAltitude(childId) > tree->getAltitude(newRoot))) {
                newRoot = childId;
            }
        }

        for (NodeId childId = firstChild; childId != InvalidNode;) {
            const NodeId nextChildId = tree->getNextSibling(childId);
            if (childId != newRoot && !tree->hasChild(newRoot, childId)) {
                tree->detachNode(childId);
                tree->attachNode(newRoot, childId);
                newRootChanged = true;
            }
            childId = nextChildId;
        }

        tree->setRoot(newRoot);
        tree->releaseNode(removedNodeId);
        notifyNodeRemoved(tree, removedNodeId);
        if (newRootChanged) {
            markAttributeUpdate(newRoot);
        }
        computeAttributeOnTreeNode(tree, newRoot);
    }

    NodeId absorbRemovedNonRootNode(DynamicComponentTree *tree, NodeId removedNodeId) {
        const NodeId parentId = tree->getNodeParent(removedNodeId);
        if (parentId == InvalidNode || parentId == removedNodeId || !tree->isAlive(parentId)) {
            return InvalidNode;
        }

        const bool parentChanged = tree->getFirstChild(removedNodeId) != InvalidNode;
        mergedParentAndChildren(tree, parentId, removedNodeId);
        disconnect(tree, removedNodeId, true);
        if (parentChanged) {
            markAttributeUpdate(parentId);
        }
        computeAttributeOnTreeNode(tree, parentId);

        if (tree->isAlive(parentId) && tree->getNumProperParts(parentId) == 0) {
            removedMarks_.mark(parentId);
            return parentId;
        }

        return InvalidNode;
    }

    void absorbRemovedNodes(DynamicComponentTree *tree, const std::vector<NodeId> &removedNodeIds) {
        struct Frame {
            NodeId nodeId = InvalidNode;
            NodeId nextChildId = InvalidNode;
        };

        if (tree == nullptr || removedNodeIds.empty()) {
            return;
        }

        std::vector<Frame> stack;
        stack.reserve(std::max<std::size_t>(64, removedNodeIds.size()));
        const auto makeFrame = [tree](NodeId nodeId) {
            return Frame{nodeId, tree->getFirstChild(nodeId)};
        };

        for (auto it = removedNodeIds.rbegin(); it != removedNodeIds.rend(); ++it) {
            if (*it != InvalidNode) {
                stack.push_back(makeFrame(*it));
            }
        }

        while (!stack.empty()) {
            Frame &frame = stack.back();
            NodeId currentNodeId = InvalidNode;

            if (!canAbsorbRemovedNode(tree, frame.nodeId)) {
                stack.pop_back();
            } else if (frame.nextChildId != InvalidNode) {
                const NodeId childId = frame.nextChildId;
                frame.nextChildId = tree->getNextSibling(childId);
                stack.push_back(makeFrame(childId));
            } else {
                currentNodeId = frame.nodeId;
                stack.pop_back();
            }

            if (canAbsorbRemovedNode(tree, currentNodeId)) {
                if (tree->isRoot(currentNodeId)) {
                    absorbRemovedRootNode(tree, currentNodeId);
                } else {
                    const NodeId parentId = absorbRemovedNonRootNode(tree, currentNodeId);
                    if (parentId != InvalidNode) {
                        stack.push_back(makeFrame(parentId));
                    }
                }
            }
        }
    }

    /** @brief Collects the subtree in post-order to ensure leaves are removed first. */
    static void collectPostOrderNodes(const DynamicComponentTree &tree, NodeId nodeId, std::vector<NodeId> &out) {
        for (NodeId childId : tree.getChildren(nodeId)) {
            collectPostOrderNodes(tree, childId, out);
        }
        out.push_back(nodeId);
    }

public:
    /**
     * @brief Builds the `leaf` adjuster over a fixed dynamic min-tree / max-tree pair.
     * @details The arguments are borrowed and remain the operational pair of
     * the adjuster. They must outlive it and remain at stable addresses. The
     * trees must be a min-tree and max-tree for the same image dimensions,
     * current image state, and adjacency. The public `prune*AndUpdate*` methods
     * operate on this same pair and do not rebind trees per call.
     *
     * Null tree pointers and matching image dimensions are assertion-checked.
     * Tree polarity, shared image state, and adjacency consistency remain caller
     * responsibilities.
     *
     * @param mintree Non-null pointer to the externally owned dynamic min-tree.
     * @param maxtree Non-null pointer to the externally owned dynamic max-tree.
     * @param graph Externally owned adjacency relation shared by both trees.
     */
    DualMinMaxTreeIncrementalFilterLeaf(DynamicComponentTree *mintree, DynamicComponentTree *maxtree, AdjacencyRelation &graph)
        : mintree_(mintree), maxtree_(maxtree), graph_(&graph),
          mergeNodesByLevel_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)),
          removedMarks_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)),
          pixelsInLeafMarks_([mintree, maxtree]() {
              const DynamicComponentTree *tree = mintree != nullptr ? mintree : maxtree;
              if (tree == nullptr) {
                  return 0;
              }
              return tree->getNumRowsOfImage() * tree->getNumColsOfImage();
          }()),
          climbedNodeMarks_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)),
          attributeUpdateMarks_(std::max(mintree ? mintree->getNumInternalNodeSlots() : 0, maxtree ? maxtree->getNumInternalNodeSlots() : 0)) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        assert(graph_ != nullptr);
        assert(mintree_->getNumRowsOfImage() == maxtree_->getNumRowsOfImage());
        assert(mintree_->getNumColsOfImage() == maxtree_->getNumColsOfImage());
    }

    /** @brief Exposes the adjacent nodes collected in the current step. */
    std::vector<NodeId> &getAdjacentNodes() {
        return mergeNodesByLevel_.getAdjacentNodes();
    }

    /** @brief Exposes the frontiers above `b` collected in the current step. */
    std::vector<NodeId> &getFrontierNodesAboveB() {
        return mergeNodesByLevel_.getFrontierNodesAboveB();
    }

    /** @brief Exposes the merge bucket associated with a level. */
    std::vector<NodeId> &getMergedNodes(const PixelType &level) {
        return mergeNodesByLevel_.getMergedNodes(level);
    }

    /**
     * @brief Registers the incremental attribute computers and external buffers.
     * @details All four arguments are borrowed. The computers must be bound to
     * the same fixed tree pair supplied to the constructor and must outlive the
     * adjuster (or remain alive until replaced by another registration). Each
     * writable span must cover every global `NodeId` slot its corresponding
     * tree may address during an update.
     *
     * Reallocating a backing container leaves the stored span dangling. Tree
     * growth can also make an existing span too small. In either case, resize
     * the buffers and call this method again before another update.
     *
     * @param computerMin Externally owned computer bound to the min-tree.
     * @param computerMax Externally owned computer bound to the max-tree.
     * @param bufferMin Borrowed writable buffer indexed by min-tree `NodeId`.
     * @param bufferMax Borrowed writable buffer indexed by max-tree `NodeId`.
     */
    void setAttributeComputer(DynamicAttributeComputer &computerMin, DynamicAttributeComputer &computerMax, std::span<float> bufferMin, std::span<float> bufferMax) {
        attrComputerMin_ = &computerMin;
        attrComputerMax_ = &computerMax;
        bufferMin_ = bufferMin;
        bufferMax_ = bufferMax;
    }

    /**
     * @brief Builds the merge buckets and frontiers for the removed leaf.
     *
     * Starting from set `C` formed by the proper parts of the leaf in
     * `primalTree`, the method
     * collects adjacent nodes in the dual tree and distributes the relevant
     * paths up to `nodeCa` into altitude buckets.
     *
     * @param dualTree Exact dual member of the fixed constructor pair.
     * @param primalTree Opposite member of that pair containing `leafId`.
     * @param leafId Live non-root leaf whose non-empty direct proper parts form `C`.
     * @param nodeCa Live extremal owner of `C` in `dualTree`.
     * @param b Altitude of the live primal parent of `leafId`.
     * @param isMaxtree Whether `dualTree` is the max-tree; must match its kind.
     * @pre Both trees are synchronized, all pixels of the leaf have valid dual
     * ownership, and `b` plus every visited altitude fit in `PixelType`.
     * @warning This low-level phase validates its contract only with assertions.
     * Invalid inputs can access invalid storage or corrupt collection state.
     */
    void buildMergedAndNestedCollections(DynamicComponentTree *dualTree, DynamicComponentTree *primalTree, NodeId leafId, NodeId nodeCa, PixelType b, bool isMaxtree) {
        assert(dualTree != nullptr);
        assert(primalTree != nullptr);
        mergeNodesByLevel_.resetCollection(isMaxtree);
        assert(nodeCa != InvalidNode);
        const PixelType altitudeCa = static_cast<PixelType>(dualTree->getAltitude(nodeCa));
        pixelsInLeafMarks_.resetAll();
        climbedNodeMarks_.resetAll();
        for (PixelId pixelId : primalTree->getProperParts(leafId)) {
            pixelsInLeafMarks_.mark(static_cast<size_t>(pixelId));
        }

        for (PixelId p : primalTree->getProperParts(leafId)) {
            for (PixelId q : graph_->getNeighborPixels(p)) {
                if (pixelsInLeafMarks_.isMarked(static_cast<size_t>(q))) {
                    continue;
                }

                const NodeId nodeQ = dualTree->getSmallestComponent(q);
                if (nodeQ == InvalidNode) {
                    continue;
                }

                const PixelType altitudeQ = static_cast<PixelType>(dualTree->getAltitude(nodeQ));
                const bool validSeed = (isMaxtree && altitudeQ >= altitudeCa) ||
                                       (!isMaxtree && altitudeQ <= altitudeCa);
                if (!validSeed || !mergeNodesByLevel_.markAdjacentSeed(nodeQ)) {
                    continue;
                }

                NodeId nodeSubtree = nodeQ;
                NodeId n = nodeQ;
                while (n != InvalidNode && dualTree->isAlive(n) && !climbedNodeMarks_.isMarked(static_cast<size_t>(n))) {
                    const PixelType levelCurrent = static_cast<PixelType>(dualTree->getAltitude(n));
                    if (!((isMaxtree && levelCurrent >= altitudeCa) ||
                          (!isMaxtree && levelCurrent <= altitudeCa))) {
                        break;
                    }

                    climbedNodeMarks_.mark(static_cast<size_t>(n));
                    nodeSubtree = n;

                    if ((isMaxtree && levelCurrent <= b) || (!isMaxtree && levelCurrent >= b)) {
                        mergeNodesByLevel_.addMergeNode(*dualTree, nodeSubtree);
                    } else {
                        NodeId parentId = dualTree->getNodeParent(nodeSubtree);
                        if (parentId == nodeSubtree) {
                            parentId = InvalidNode;
                        }
                        if (!(parentId != InvalidNode &&
                              ((isMaxtree && dualTree->getAltitude(parentId) > b) ||
                               (!isMaxtree && dualTree->getAltitude(parentId) < b)))) {
                            mergeNodesByLevel_.addFrontierNodeAboveB(nodeSubtree);
                        }
                    }

                    const NodeId parentId = dualTree->getNodeParent(n);
                    if (parentId == n) {
                        break;
                    }
                    n = parentId;
                }
            }
        }
    }

    /**
     * @brief Updates the dual tree in preparation for removing one primal leaf.
     *
     * This is the algorithmic core of `leaf` mode. The method:
     *
     * - identifies `b` and `nodeCa`;
     * - builds the merge buckets;
     * - executes the altitude sweep;
     * - moves proper parts and reattaches frontiers;
     * - handles cases where `nodeCa` disappears or the root must change.
     *
     * `dualTree` must be exactly one member of the fixed constructor pair; the
     * opposite member is inferred as the primal tree. `leafId` must still be a
     * live non-root leaf in that primal tree, with at least one direct proper
     * part and a live dual owner for every relevant pixel. The leaf/parent edge
     * must be strict and canonical: the parent is lower than the leaf in a
     * max-tree and higher than the leaf in a min-tree. Both trees must be
     * synchronized before entry, all altitudes must fit in `PixelType`, and any
     * registered attribute runtime must satisfy its lifetime and sizing contract.
     *
     * This method mutates only the dual tree. It deliberately reads the leaf
     * before pruning, so callers must prune the primal leaf only after this
     * method succeeds; the `prune*AndUpdate*` wrappers enforce that order.
     *
     * @param dualTree Dual member of the fixed pair, updated in place.
     * @param leafId Live non-root leaf in the inferred primal tree.
     * @warning Membership of `dualTree` in the fixed pair is not checked: any
     * pointer other than the stored max-tree is treated as the stored min-tree.
     * Other preconditions are primarily assertion-only. Invalid input can cause
     * undefined behavior or structural corruption. The operation is not
     * transactional and does not roll back partial dual-tree mutations if an
     * internal operation fails.
     */
    void updateTree(DynamicComponentTree *dualTree, NodeId leafId) {
        assert(dualTree != nullptr);
        assert(leafId != InvalidNode);

        const bool isMaxtree = dualTree == maxtree_;
        initializeUpdateStepState(isMaxtree);
        DynamicComponentTree *primalTree = getPrimalTree(isMaxtree);

        assert(primalTree != nullptr);
        assert(primalTree->isNode(leafId) && primalTree->isAlive(leafId));
        assert(primalTree->isLeaf(leafId));

        const NodeId leafParent = primalTree->getNodeParent(leafId);
        assert(leafParent != InvalidNode && leafParent != leafId);
        const PixelType b = static_cast<PixelType>(primalTree->getAltitude(leafParent));

        assert(primalTree->getNumProperParts(leafId) > 0);
        NodeId nodeCa = InvalidNode;
        PixelType altitudeCa = PixelType{};
        for (PixelId pixelId : primalTree->getProperParts(leafId)) {
            const NodeId ownerId = dualTree->getSmallestComponent(pixelId);
            if (ownerId == InvalidNode) {
                continue;
            }

            const PixelType ownerAltitude = static_cast<PixelType>(dualTree->getAltitude(ownerId));
            if (nodeCa == InvalidNode ||
                (isMaxtree && ownerAltitude < altitudeCa) ||
                (!isMaxtree && ownerAltitude > altitudeCa)) {
                nodeCa = ownerId;
                altitudeCa = ownerAltitude;
            }
        }
        assert(nodeCa != InvalidNode);
        altitudeCa_ = altitudeCa;

        buildMergedAndNestedCollections(dualTree, primalTree, leafId, nodeCa, b, isMaxtree);

        PixelType mergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId unionNode = InvalidNode;
        NodeId previousUnionNode = InvalidNode;
        bool nodeCaRemoved = false;

        while (mergeNodesByLevel_.hasMergeLevel() && ((isMaxtree && mergeLevel > altitudeCa) || (!isMaxtree && mergeLevel < altitudeCa))) {
            auto &mergeNodesAtLevel = mergeNodesByLevel_.getMergedNodes(static_cast<PixelType>(mergeLevel));
            unionNode = InvalidNode;

            for (auto nodeId : mergeNodesAtLevel) {
                if (!dualTree->isAlive(nodeId)) {
                    continue;
                }
                if (unionNode == InvalidNode || nodeId < unionNode) {
                    unionNode = nodeId;
                }
            }

            for (auto nodeId : mergeNodesAtLevel) {
                if (!dualTree->isAlive(nodeId)) {
                    continue;
                }

                if (nodeId == unionNode) {
                    disconnect(dualTree, unionNode, false);
                    continue;
                }

                if (unionNode == InvalidNode) {
                    continue;
                }

                mergedParentAndChildren(dualTree, unionNode, nodeId);
                disconnect(dualTree, nodeId, true);
            }

            if (unionNode == InvalidNode) {
                mergeLevel = mergeNodesByLevel_.nextMergeLevel();
                unionNode = previousUnionNode;
                continue;
            }

            if (mergeLevel == b) {
                for (auto pixelId : primalTree->getProperParts(leafId)) {
                    const NodeId ownerId = dualTree->getSmallestComponent(pixelId);
                    if (ownerId == InvalidNode || ownerId == unionNode) {
                        continue;
                    }

                    dualTree->moveProperPart(unionNode, ownerId, pixelId);
                    notifyMoveProperPart(dualTree, unionNode, ownerId, pixelId);
                    if (dualTree->isAlive(ownerId) && dualTree->getNumProperParts(ownerId) == 0 && !removedMarks_.isMarked(ownerId)) {
                        removedMarks_.mark(ownerId);
                        removedNodesPendingAbsorption_.push_back(ownerId);
                        if (ownerId == nodeCa) {
                            nodeCaRemoved = true;
                        }
                    }
                }

                for (auto nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                    disconnect(dualTree, nodeId, false);
                    dualTree->attachNode(unionNode, nodeId);
                }
            }

            if (previousUnionNode != InvalidNode && dualTree->isAlive(previousUnionNode) && !dualTree->hasChild(unionNode, previousUnionNode)) {
                if (!dualTree->isRoot(previousUnionNode)) {
                    disconnect(dualTree, previousUnionNode, false);
                }
                dualTree->attachNode(unionNode, previousUnionNode);
            }

            markAttributeUpdate(unionNode);
            computeAttributeOnTreeNode(dualTree, unionNode);
            previousUnionNode = unionNode;
            mergeLevel = mergeNodesByLevel_.nextMergeLevel();
        }

        const NodeId finalUnionNode = previousUnionNode;
        if (finalUnionNode == InvalidNode || !dualTree->isAlive(finalUnionNode)) {
            return;
        }

        if (nodeCaRemoved) {
            if (!dualTree->isRoot(nodeCa)) {
                const NodeId parentIdNodeCa = dualTree->getNodeParent(nodeCa);
                if (!dualTree->isRoot(finalUnionNode)) {
                    disconnect(dualTree, finalUnionNode, false);
                }
                dualTree->attachNode(parentIdNodeCa, finalUnionNode);

                for (NodeId n = dualTree->getFirstChild(nodeCa); n != InvalidNode;) {
                    const NodeId next = dualTree->getNextSibling(n);
                    if (n != finalUnionNode && !dualTree->hasChild(finalUnionNode, n)) {
                        if (!dualTree->isRoot(n)) {
                            disconnect(dualTree, n, false);
                        }
                        dualTree->attachNode(finalUnionNode, n);
                    }
                    n = next;
                }

                disconnect(dualTree, nodeCa, true);
                markAttributeUpdate(finalUnionNode);
                computeAttributeOnTreeNode(dualTree, finalUnionNode);
                markAttributeUpdate(parentIdNodeCa);
                computeAttributeOnTreeNode(dualTree, parentIdNodeCa);
            } else {
                NodeId newRoot = finalUnionNode;
                bool newRootChanged = false;
                for (NodeId n = dualTree->getFirstChild(nodeCa); n != InvalidNode; n = dualTree->getNextSibling(n)) {
                    if ((isMaxtree && dualTree->getAltitude(n) < dualTree->getAltitude(newRoot)) || (!isMaxtree && dualTree->getAltitude(n) > dualTree->getAltitude(newRoot))) {
                        newRoot = n;
                    }
                }

                if (newRoot != finalUnionNode) {
                    if (!dualTree->isRoot(finalUnionNode)) {
                        disconnect(dualTree, finalUnionNode, false);
                    }
                    dualTree->attachNode(newRoot, finalUnionNode);
                    newRootChanged = true;
                }

                for (NodeId n = dualTree->getFirstChild(nodeCa); n != InvalidNode;) {
                    const NodeId next = dualTree->getNextSibling(n);
                    if (n != newRoot && !dualTree->hasChild(newRoot, n)) {
                        if (!dualTree->isRoot(n)) {
                            disconnect(dualTree, n, false);
                        }
                        dualTree->attachNode(newRoot, n);
                        newRootChanged = true;
                    }
                    n = next;
                }

                dualTree->setRoot(newRoot);
                dualTree->releaseNode(nodeCa);
                notifyNodeRemoved(dualTree, nodeCa);
                if (newRootChanged) {
                    markAttributeUpdate(newRoot);
                }
                computeAttributeOnTreeNode(dualTree, newRoot);
            }
        } else if (finalUnionNode != nodeCa) {
            if (!dualTree->isRoot(finalUnionNode)) {
                disconnect(dualTree, finalUnionNode, false);
            }
            dualTree->attachNode(nodeCa, finalUnionNode);
            markAttributeUpdate(nodeCa);
            computeAttributeOnTreeNode(dualTree, nodeCa);
        }

        absorbRemovedNodes(dualTree, removedNodesPendingAbsorption_);
    }

    /**
     * @brief Prunes max-tree subtrees and updates this adjuster's min-tree leaf by leaf.
     *
     * Each subtree is converted to post-order and processed only while the
     * nodes remain alive, ensuring that the elementary step remains a leaf
     * removal. For every accepted leaf, `updateTree(mintree_, leafId)` runs
     * while the max-tree leaf is still alive; `maxtree_->pruneNode(leafId)` runs
     * only afterward.
     *
     * @param nodesToPrune Max-tree subtree roots. Invalid ids, the global root,
     * dead nodes, and nodes removed by an earlier entry are silently skipped.
     * @pre Each accepted subtree is structurally valid and the fixed pair and
     * optional attribute runtime satisfy the `updateTree()` contract.
     * @warning Leafness and other structural requirements are assertion-only.
     * The batch is not transactional; completed updates are not rolled back.
     */
    void pruneMaxTreeAndUpdateMinTree(std::vector<NodeId> &nodesToPrune) {
        std::vector<NodeId> postOrderNodes;
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == maxtree_->getRoot() || !maxtree_->isNode(rootSubtree) || !maxtree_->isAlive(rootSubtree)) {
                continue;
            }

            postOrderNodes.clear();
            collectPostOrderNodes(*maxtree_, rootSubtree, postOrderNodes);
            for (NodeId leafId : postOrderNodes) {
                if (leafId == InvalidNode || leafId == maxtree_->getRoot() || !maxtree_->isNode(leafId) || !maxtree_->isAlive(leafId)) {
                    continue;
                }
                assert(maxtree_->isLeaf(leafId));
                updateTree(mintree_, leafId);
                maxtree_->pruneNode(leafId);
            }
        }
    }

    /**
     * @brief Updates the min-tree, then prunes one max-tree leaf.
     * @param leafId Live non-root max-tree leaf. Invalid ids, the root, and dead
     * nodes are silently ignored; leafness is checked only by assertion.
     * @warning The update is not transactional and does not roll back a partial
     * dual mutation on failure.
     */
    void pruneMaxTreeAndUpdateMinTree(NodeId leafId) {
        if (leafId == InvalidNode || leafId == maxtree_->getRoot() || !maxtree_->isNode(leafId) || !maxtree_->isAlive(leafId)) {
            return;
        }
        assert(maxtree_->isLeaf(leafId));
        updateTree(mintree_, leafId);
        maxtree_->pruneNode(leafId);
    }

    /**
     * @brief Prunes min-tree subtrees and updates this adjuster's max-tree leaf by leaf.
     * @details Each accepted leaf is first passed to
     * `updateTree(maxtree_, leafId)` while still alive, and is pruned from the
     * min-tree only after the dual update completes.
     * @param nodesToPrune Min-tree subtree roots. Invalid ids, the global root,
     * dead nodes, and nodes removed by an earlier entry are silently skipped.
     * @pre Each accepted subtree is structurally valid and the fixed pair and
     * optional attribute runtime satisfy the `updateTree()` contract.
     * @warning Leafness and other structural requirements are assertion-only.
     * The batch is not transactional; completed updates are not rolled back.
     */
    void pruneMinTreeAndUpdateMaxTree(std::vector<NodeId> &nodesToPrune) {
        std::vector<NodeId> postOrderNodes;
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == mintree_->getRoot() || !mintree_->isNode(rootSubtree) || !mintree_->isAlive(rootSubtree)) {
                continue;
            }

            postOrderNodes.clear();
            collectPostOrderNodes(*mintree_, rootSubtree, postOrderNodes);
            for (NodeId leafId : postOrderNodes) {
                if (leafId == InvalidNode || leafId == mintree_->getRoot() || !mintree_->isNode(leafId) || !mintree_->isAlive(leafId)) {
                    continue;
                }
                assert(mintree_->isLeaf(leafId));
                updateTree(maxtree_, leafId);
                mintree_->pruneNode(leafId);
            }
        }
    }

    /**
     * @brief Updates the max-tree, then prunes one min-tree leaf.
     * @param leafId Live non-root min-tree leaf. Invalid ids, the root, and dead
     * nodes are silently ignored; leafness is checked only by assertion.
     * @warning The update is not transactional and does not roll back a partial
     * dual mutation on failure.
     */
    void pruneMinTreeAndUpdateMaxTree(NodeId leafId) {
        if (leafId == InvalidNode || leafId == mintree_->getRoot() || !mintree_->isNode(leafId) || !mintree_->isAlive(leafId)) {
            return;
        }
        assert(mintree_->isLeaf(leafId));
        updateTree(maxtree_, leafId);
        mintree_->pruneNode(leafId);
    }

    /**
     * @brief Returns mutable access to the borrowed min-tree.
     * @return Non-owning pointer valid only while the external tree remains alive.
     * @warning External edits can desynchronize the fixed pair and invalidate
     * subsequent leaf-update assumptions.
     */
    DynamicComponentTree *getMinTree() const { return mintree_; }

    /**
     * @brief Returns mutable access to the borrowed max-tree.
     * @return Non-owning pointer valid only while the external tree remains alive.
     * @warning External edits can desynchronize the fixed pair and invalidate
     * subsequent leaf-update assumptions.
     */
    DynamicComponentTree *getMaxTree() const { return maxtree_; }
};
