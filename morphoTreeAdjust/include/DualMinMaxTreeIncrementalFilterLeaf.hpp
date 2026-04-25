#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <limits>
#include <span>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "DynamicComponentTree.hpp"

#ifndef PRINT_LOG
#define PRINT_LOG 0
#endif

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
 * - a leaf is removed from the primal tree;
 * - the proper parts of that leaf define set `C` in the dual tree;
 * - the relevant neighbors of `C` are collected;
 * - the dual tree is updated by a level sweep, analogously to the `subtree`
 *   adjuster, but usually over a smaller region.
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
 * - `A_L`: number of adjacent nodes collected in the step;
 * - `M_L`: number of nodes placed in the merge buckets;
 * - `K_L`: number of nodes effectively edited in the dual tree.
 *
 * An individual step has local cost, approximately:
 *
 * - `O(P_L + A_L + M_L + K_L)`.
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
    std::ostringstream outputLog_;

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

    /** @brief Writes a compact textual description of a node to the log. */
    void appendNodeState(DynamicComponentTree *tree, NodeId nodeId) {
        outputLog_ << "id:" << nodeId
                   << ", level:" << tree->getAltitude(nodeId)
                   << ", |cnps|:" << tree->getNumProperParts(nodeId)
                   << ", |children|:" << tree->getNumChildren(nodeId);
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
     * @details The pair passed here remains the operational pair of the adjuster.
     * The public `prune*AndUpdate*` methods act on this same pair and do not
     * rebind trees per call.
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

    /** @brief Clears the textual log used when `PRINT_LOG` is active. */
    void clearOutputLog() {
        outputLog_.str("");
        outputLog_.clear();
    }

    /** @brief Returns the current contents of the debug log. */
    std::string getOutputLog() const {
        return outputLog_.str();
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
     * @details The registered computers and buffers must match the same fixed
     * tree pair supplied to the constructor.
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
     * @brief Updates the dual tree after the removal of a single leaf.
     *
     * This is the algorithmic core of `leaf` mode. The method:
     *
     * - identifies `b` and `nodeCa`;
     * - builds the merge buckets;
     * - executes the altitude sweep;
     * - moves proper parts and reattaches frontiers;
     * - handles cases where `nodeCa` disappears or the root must change.
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

        if (PRINT_LOG) {
            clearOutputLog();
            outputLog_ << "Atributo(leaf)= " << primalTree->getNumProperParts(leafId) << ", nivel(leaf)= " << static_cast<double>(primalTree->getAltitude(leafId)) << ", nivel(parent(leaf))= " << static_cast<double>(primalTree->getAltitude(leafParent)) << "\n";
            outputLog_ << "b: " << static_cast<double>(b) << "\n";
            outputLog_ << "Proper parts (tau_L): [";
            bool first = true;
            for (auto rep : primalTree->getProperParts(leafId)) {
                const NodeId nodeTau = dualTree->getSmallestComponent(rep);
                if (nodeTau == InvalidNode) {
                    continue;
                }
                if (!first) {
                    outputLog_ << "\t";
                }
                first = false;
                outputLog_ << "\t(id:" << nodeTau << ", level:" << static_cast<double>(dualTree->getAltitude(nodeTau)) << ", |cnps|:" << dualTree->getNumProperParts(nodeTau) << ", rep:" << rep << "), \n";
            }
            outputLog_ << "]\n";
            if (isMaxtree) {
                outputLog_ << "Intervalo: [" << static_cast<double>(altitudeCa) << ", " << static_cast<double>(b) << "]\n";
            } else {
                outputLog_ << "Intervalo: [" << static_cast<double>(b) << ", " << static_cast<double>(altitudeCa) << "]\n";
            }
            outputLog_ << "nodeCa: Id:" << nodeCa << "; level:" << static_cast<double>(altitudeCa) << "; |cnps|:" << dualTree->getNumProperParts(nodeCa) << "\n";
        }

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

                if (PRINT_LOG) {
                    outputLog_ << "frontier attach @b: unionNode=(";
                    appendNodeState(dualTree, unionNode);
                    outputLog_ << "), frontier_count=" << mergeNodesByLevel_.getFrontierNodesAboveB().size() << "\n";
                    for (auto nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                        outputLog_ << "\tfrontier->union: frontier=(";
                        appendNodeState(dualTree, nodeId);
                        outputLog_ << ") -> unionNode=" << unionNode << "\n";
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
     * removal.
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

    /** @brief Prunes one max-tree leaf and updates this adjuster's min-tree. */
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

    /** @brief Prunes one min-tree leaf and updates this adjuster's max-tree. */
    void pruneMinTreeAndUpdateMaxTree(NodeId leafId) {
        if (leafId == InvalidNode || leafId == mintree_->getRoot() || !mintree_->isNode(leafId) || !mintree_->isAlive(leafId)) {
            return;
        }
        assert(mintree_->isLeaf(leafId));
        updateTree(maxtree_, leafId);
        mintree_->pruneNode(leafId);
    }

    DynamicComponentTree *getMinTree() const { return mintree_; }
    DynamicComponentTree *getMaxTree() const { return maxtree_; }
};
