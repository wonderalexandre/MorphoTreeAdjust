#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

#include "AdjacencyRelation.hpp"
#include "Common.hpp"

/**
 * @brief Dynamic component tree built directly from the image pixels.
 *
 * This is the central structure of the project's per-pixel dynamic line. It
 * represents a mutable max-tree or min-tree of the current image and was
 * designed to support frequent local updates without a full rebuild after each
 * pruning step.
 *
 * Data model:
 *
 * - image pixels are indexed by `PixelId`;
 * - hierarchy nodes are indexed by `NodeId`, in a dense id space from `0` to
 *   `numInternalNodeSlots_ - 1`;
 * - each live node represents a connected component at one level;
 * - each pixel belongs directly to exactly one live node;
 * - the full support of a node is given by its subtree.
 *
 * Hierarchy encoding:
 *
 * - the parent relation is encoded by `nodeParent_`, which stores the direct
 *   parent of each node and therefore has size `numInternalNodeSlots_`;
 * - the child relation is encoded by the doubly linked child backend
 *   `firstChild_`, `lastChild_`, `nextSibling_`, and `prevSibling_`; each of
 *   these arrays has size `numInternalNodeSlots_`;
 * - the node-to-pixel relation for direct proper parts is encoded by
 *   `properHead_`, `properTail_`, `nextProperPart_`, and `prevProperPart_`;
 *   `properHead_` and `properTail_` are node-indexed and have size
 *   `numInternalNodeSlots_`, while `nextProperPart_` and `prevProperPart_` are
 *   pixel-indexed and have size `numTotalProperParts_`;
 * - the pixel-to-node relation is encoded by `properPartOwner_`, which stores
 *   the current direct owner of each pixel and has size
 *   `numTotalProperParts_`;
 * - `altitude_` stores the canonical level associated with each node and has
 *   size `numInternalNodeSlots_`;
 * - `freeNodeIds_` stores released node slots that may later be reused; its
 *   current size is dynamic and never exceeds `numInternalNodeSlots_`.

 *
 * Fundamental invariants:
 *
 * - the `rootNodeId_` root is a live node and points to itself as parent;
 * - every live node other than the root has a live parent;
 * - every pixel in the domain has a single direct owner;
 * - the direct pixels of a node are uniform at level `altitude_[nodeId]`;
 * - concatenating the proper parts of a node subtree reconstructs its full
 *   support;
 * - the id of a live node never changes during its lifetime;
 * - released nodes leave the topology and their slots can only reappear
 *   through reuse via `allocateNode`.
 *
 * Role in the pipeline:
 *
 * - serves as the primal and dual structure for the dynamic adjusters;
 * - supports CASF execution in `leaf` and `subtree` modes;
 * - provides the proper parts and topology used by the dynamic adjusters and
 *   CASF routines built on top of this representation.
 *
 * Operation types for which the class was optimized:
 *
 * - moving proper parts between nodes;
 * - merging neighboring nodes;
 * - detaching nodes and pruning subtrees;
 * - reattaching children and updating the root;
 * - iterating quickly over children, proper parts, and subtrees.
 *
 * Structural complexity:
 *
 * Use:
 *
 * - `P`: number of pixels;
 * - `N`: number of nodes materialized in the global id space.
 *
 * Main costs:
 *
 * - space: `Theta(P + N)`;
 * - local topology and metadata queries (`getParent`, `getAltitude`,
 *   `isLeaf`, `isAlive`, `getProperPartOwner`): `O(1)`;
 * - direct-child iteration: linear in the number of children;
 * - direct proper-part iteration: linear in the number of owned pixels;
 * - subtree iteration: linear in the number of visited nodes;
 * - moving or attaching proper parts: linear in the number of pixels actually
 *   moved;
 * - merges and prunes: cost proportional to the size of the affected
 *   region/locality, not to the total tree size.
 *
 * In practice, the class does not by itself determine the cost of the CASF
 * pipeline: that cost still depends heavily on the dynamic adjusters. Even so,
 * several representation choices here directly affect how efficiently those
 * adjusters can perform local edits:
 *
 * - node ids give a stable array-based address space while a node is live,
 *   enabling `O(1)` metadata access, reusable external mark/buffer structures,
 *   and slot reuse without pointer allocation churn;
 * - storing both parent links and explicit child links avoids global scans
 *   when the adjuster needs to detach, reattach, promote, or traverse local
 *   subtrees;
 * - doubly linked child lists make local unlinks and reattachments cheap once
 *   the participating nodes are known, which is a common pattern in merge and
 *   contraction steps;
 * - keeping both pixel-to-node ownership and node-to-pixel proper-part lists
 *   supports the two directions required by the update algorithms: locate the
 *   current owner of a pixel, and move one or many direct proper parts between
 *   nodes without rebuilding component supports from scratch;
 * - using doubly linked proper-part lists keeps single-pixel removals,
 *   insertions, and concatenations local to the edited nodes.
 *
 * The tradeoff is extra `Theta(P + N)` storage to maintain these parallel
 * backends, but the benefit is that the adjusters can operate on the affected
 * region instead of paying for a global rebuild after each pruning step.
 */
class DynamicComponentTree {
private:
    // Permanent metadata of the base image and the adjacency used in construction.
    AdjacencyRelationPtr adj_;
    int numRows_ = 0;
    int numCols_ = 0;
    bool isMaxtree_ = true;

    // Size of the pixel domain and the node domain.
    int numTotalProperParts_ = 0;
    int numInternalNodeSlots_ = 0;
    NodeId rootNodeId_ = InvalidNode;

    // Main hierarchy structure indexed only by nodes.
    std::vector<int> altitude_; //This is unique attribute of the node, so it is stored in a separate vector for better cache performance.
    std::vector<NodeId> freeNodeIds_;
    std::vector<NodeId> nodeParent_;
    
    // Backend for the doubly linked child lists.
    std::vector<NodeId> firstChild_;
    std::vector<NodeId> lastChild_;
    std::vector<NodeId> nextSibling_;
    std::vector<NodeId> prevSibling_;
    std::vector<int> numChildrenByNode_;

    // Backend for the doubly linked lists of direct proper parts.
    std::vector<PixelId> properHead_;
    std::vector<PixelId> properTail_;
    std::vector<int> numProperPartsByNode_;
    std::vector<NodeId> properPartOwner_;
    std::vector<PixelId> nextProperPart_;
    std::vector<PixelId> prevProperPart_;

    // Mutation counters used in internal checks and tests.
    std::size_t nodeStructureVersion_ = 0;
    std::size_t topologyVersion_ = 0;
    std::size_t properPartVersion_ = 0;

    /**
     * @brief Initializes the tree backend vectors.
     * @param numProperParts Number of image pixels.
     * @param numInternalNodeSlots Number of nodes materialized in the initial build.
     *
     * The node vectors and pixel vectors are kept separate.
     */
    void initializeStorage(int numProperParts, int numInternalNodeSlots) {
        numTotalProperParts_ = numProperParts;
        numInternalNodeSlots_ = numInternalNodeSlots;
        rootNodeId_ = InvalidNode;

        altitude_.assign((size_t) std::max(0, numInternalNodeSlots), 0);
        freeNodeIds_.clear();
        nodeParent_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        firstChild_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        lastChild_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        nextSibling_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        prevSibling_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        numChildrenByNode_.assign((size_t) std::max(0, numInternalNodeSlots), 0);
        properHead_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        properTail_.assign((size_t) std::max(0, numInternalNodeSlots), InvalidNode);
        numProperPartsByNode_.assign((size_t) std::max(0, numInternalNodeSlots), 0);

        properPartOwner_.assign((size_t) std::max(0, numTotalProperParts_), InvalidNode);
        nextProperPart_.assign((size_t) std::max(0, numTotalProperParts_), InvalidNode);
        prevProperPart_.assign((size_t) std::max(0, numTotalProperParts_), InvalidNode);

        nodeStructureVersion_ = 0;
        topologyVersion_ = 0;
        properPartVersion_ = 0;
    }

    /**
     * @brief Appends a pixel to the list of direct proper parts of a node.
     * @param nodeId Node that will directly own the pixel.
     * @param pixelId Pixel to append.
     */
    void appendProperPartToNode(NodeId nodeId, PixelId pixelId) {
        properPartOwner_[pixelId] = nodeId;
        if (properHead_[nodeId] == InvalidNode) {
            properHead_[nodeId] = pixelId;
            properTail_[nodeId] = pixelId;
        } else {
            const PixelId tail = properTail_[nodeId];
            nextProperPart_[tail] = pixelId;
            prevProperPart_[pixelId] = tail;
            properTail_[nodeId] = pixelId;
        }
        numProperPartsByNode_[nodeId]++;
    }

    /**
     * @brief Inserts a child at the end of a node's child list.
     * @param parentId Parent node.
     * @param childId Direct child node to append.
     */
    void linkChildBack(NodeId parentId, NodeId childId) {
        const NodeId tail = lastChild_[parentId];
        if (tail == InvalidNode) {
            firstChild_[parentId] = childId;
            lastChild_[parentId] = childId;
        } else {
            nextSibling_[tail] = childId;
            prevSibling_[childId] = tail;
            lastChild_[parentId] = childId;
        }
        numChildrenByNode_[parentId]++;
    }

    /**
     * @brief Removes a node from the child list of its current parent.
     * @param childId Child to detach.
     */
    void unlinkChild(NodeId childId) {
        const NodeId parentId = nodeParent_[childId];
        const NodeId prev = prevSibling_[childId];
        const NodeId next = nextSibling_[childId];

        if (prev != InvalidNode) {
            nextSibling_[prev] = next;
        } else {
            firstChild_[parentId] = next;
        }

        if (next != InvalidNode) {
            prevSibling_[next] = prev;
        } else {
            lastChild_[parentId] = prev;
        }

        prevSibling_[childId] = InvalidNode;
        nextSibling_[childId] = InvalidNode;
        numChildrenByNode_[parentId]--;
    }

    /**
     * @brief Moves all direct proper parts from one node to another.
     *
     * The backend only relinks the lists and updates the owner node of each pixel.
     */
    void moveProperPartsInBackend(NodeId targetNodeId, NodeId sourceNodeId) {
        const PixelId movedHead = properHead_[sourceNodeId];
        if (movedHead == InvalidNode) {
            return;
        }

        for (PixelId pixelId = movedHead; pixelId != InvalidNode; pixelId = nextProperPart_[pixelId]) {
            properPartOwner_[pixelId] = targetNodeId;
        }

        if (numProperPartsByNode_[targetNodeId] == 0) {
            properHead_[targetNodeId] = movedHead;
            properTail_[targetNodeId] = properTail_[sourceNodeId];
            numProperPartsByNode_[targetNodeId] = numProperPartsByNode_[sourceNodeId];
        } else {
            prevProperPart_[movedHead] = properTail_[targetNodeId];
            nextProperPart_[properTail_[targetNodeId]] = movedHead;
            properTail_[targetNodeId] = properTail_[sourceNodeId];
            numProperPartsByNode_[targetNodeId] += numProperPartsByNode_[sourceNodeId];
        }

        properHead_[sourceNodeId] = InvalidNode;
        properTail_[sourceNodeId] = InvalidNode;
        numProperPartsByNode_[sourceNodeId] = 0;
    }

    /**
     * @brief Detaches a node from its parent and turns it into the root of an isolated subtree.
     */
    void detachNodeInBackend(NodeId nodeId) {
        unlinkChild(nodeId);
        nodeParent_[nodeId] = nodeId;
    }

    /**
     * @brief Releases a node completely.
     *
     * The node index becomes available for reuse by `allocateNode`.
     */
    void releaseNodeSlot(NodeId nodeId) {
        firstChild_[nodeId] = InvalidNode;
        lastChild_[nodeId] = InvalidNode;
        nextSibling_[nodeId] = InvalidNode;
        prevSibling_[nodeId] = InvalidNode;
        numChildrenByNode_[nodeId] = 0;
        properHead_[nodeId] = InvalidNode;
        properTail_[nodeId] = InvalidNode;
        numProperPartsByNode_[nodeId] = 0;
        nodeParent_[nodeId] = InvalidNode;
        altitude_[nodeId] = 0;
        freeNodeIds_.push_back(nodeId);
    }

public:
    /**
     * @brief Lightweight range of the direct children of a node.
     *
     * The range encapsulates the singly traversable list via `nextSibling_`,
     * exposing only a minimal forward interface.
     */
    class ChildrenRange;
    /**
     * @brief Lightweight range of a node's direct proper parts.
     *
     * Each proper part is a leaf/pixel of the hierarchy. The range traverses
     * the doubly linked list whose head lives in `properHead_[nodeId]`.
     */
    class ProperPartsRange;
    /**
     * @brief Range for traversing a subtree in depth-first order.
     *
     * The traversal uses iterative DFS and returns only live internal nodes.
     */
    class SubtreeNodeRange;

    /**
     * @brief Empty constructor; the tree can be built later with `build`.
     */
    DynamicComponentTree() = default;

    /**
     * @brief Range for breadth-first traversal of a subtree.
     */
    class BreadthFirstNodeRange;

    /**
     * @brief Builds the tree directly from an image.
     * @param image Grayscale image.
     * @param isMaxtree `true` for max-tree, `false` for min-tree.
     * @param adj Image adjacency relation.
     */
    DynamicComponentTree(ImageUInt8Ptr image, bool isMaxtree, AdjacencyRelationPtr adj) {
        build(image, isMaxtree, adj);
    }

    /**
     * @brief Builds or rebuilds the entire hierarchy from the input image.
     * @param image Base grayscale image.
     * @param isMaxtree `true` for max-tree, `false` for min-tree.
     * @param adj Adjacency relation used during construction.
     */
    void build(ImageUInt8Ptr image, bool isMaxtree, AdjacencyRelationPtr adj) {
        assert(image != nullptr);
        assert(adj != nullptr);
        adj_ = adj;
        numRows_ = image->getNumRows();
        numCols_ = image->getNumCols();
        isMaxtree_ = isMaxtree;
        auto orderedPixels = countingSort(image);
        createTreeByUnionFind(orderedPixels, image);
    }

    /**
     * @brief Sorts pixels for union-find construction.
     * @return Vector of pixels in min-tree/max-tree-compatible order.
     *
     * For a max-tree, the order is increasing in gray level; for a min-tree,
     * the order is inverted by the `maxvalue - gray` transform.
     */
    std::vector<PixelId> countingSort(ImageUInt8Ptr image) const {
        assert(image != nullptr);
        const uint8_t *img = image->rawData();
        const int n = numRows_ * numCols_;
        if (n == 0) {
            return {};
        }

        int maxvalue = img[0];
        for (int i = 1; i < n; ++i) {
            if (maxvalue < img[i]) {
                maxvalue = img[i];
            }
        }

        std::vector<uint32_t> counter(maxvalue + 1, 0);
        std::vector<PixelId> orderedPixels(n);

        if (isMaxtree_) {
            for (int i = 0; i < n; ++i) {
                counter[img[i]]++;
            }
            for (int i = 1; i < maxvalue; ++i) {
                counter[i] += counter[i - 1];
            }
            if (maxvalue > 0) {
                counter[maxvalue] += counter[maxvalue - 1];
            }
            for (int i = n - 1; i >= 0; --i) {
                orderedPixels[--counter[img[i]]] = i;
            }
        } else {
            for (int i = 0; i < n; ++i) {
                counter[maxvalue - img[i]]++;
            }
            for (int i = 1; i < maxvalue; ++i) {
                counter[i] += counter[i - 1];
            }
            if (maxvalue > 0) {
                counter[maxvalue] += counter[maxvalue - 1];
            }
            for (int i = n - 1; i >= 0; --i) {
                orderedPixels[--counter[maxvalue - img[i]]] = i;
            }
        }

        return orderedPixels;
    }

    /**
     * @brief Materializes the hierarchy by union-find from already sorted pixels.
     * @param orderedPixels Pixels in min/max-tree-compatible order.
     *
     * The result is produced directly in the mutable format used by the rest of
     * the class: explicit nodes, explicit proper parts, and linked child lists.
     */
    void createTreeByUnionFind(std::vector<PixelId> &orderedPixels, ImageUInt8Ptr image) {
        assert(image != nullptr);
        assert(adj_ != nullptr);
        const uint8_t *img = image->rawData();

        const int numPixels = numRows_ * numCols_;
        std::vector<int> zPar(numPixels, -1);
        std::vector<int> parent(numPixels, -1);
        std::vector<NodeId> pixelToNodeId(numPixels, InvalidNode);
        auto findRoot = [&](NodeId p) {
            while (zPar[p] != p) {
                zPar[p] = zPar[zPar[p]];
                p = zPar[p];
            }
            return p;
        };
        for (int i = numPixels - 1; i >= 0; --i) {
            const int p = orderedPixels[i];
            parent[p] = p;
            zPar[p] = p;
            for (int q : adj_->getNeighborPixels(p)) {
                if (zPar[q] != -1) {
                    const int r = findRoot(q);
                    if (p != r) {
                        parent[r] = p;
                        zPar[r] = p;
                    }
                }
            }
        }

        int numNodes = 0;
        for (int i = 0; i < numPixels; ++i) {
            const int p = orderedPixels[i];
            const int q = parent[p];
            if (img[parent[q]] == img[q]) {
                parent[p] = parent[q];
            }
            if (parent[p] == p || img[parent[p]] != img[p]) {
                ++numNodes;
            }
        }

        initializeStorage(numPixels, numNodes);

        int nextNodeId = 0;
        for (int i = 0; i < numPixels; ++i) {
            const int p = orderedPixels[i];
            if (p == parent[p]) {
                const NodeId dynamicNodeId = nextNodeId++;
                rootNodeId_ = dynamicNodeId;
                nodeParent_[dynamicNodeId] = dynamicNodeId;
                pixelToNodeId[p] = dynamicNodeId;
                altitude_[dynamicNodeId] = img[p];
            } else if (img[p] != img[parent[p]]) {
                const NodeId dynamicNodeId = nextNodeId++;
                const NodeId dynamicParentId = pixelToNodeId[parent[p]];
                nodeParent_[dynamicNodeId] = dynamicParentId;
                linkChildBack(dynamicParentId, dynamicNodeId);
                pixelToNodeId[p] = dynamicNodeId;
                altitude_[dynamicNodeId] = img[p];
            } else {
                pixelToNodeId[p] = pixelToNodeId[parent[p]];
            }
            
        }

        for (int i = 0; i < numPixels; ++i) {
            const int p = orderedPixels[i];
            appendProperPartToNode(pixelToNodeId[p], p);
        }
    }

    /**
     * @brief Total number of image pixels.
     */
    int getNumTotalProperParts() const { return numTotalProperParts_; }
    /**
     * @brief Total number of node indices maintained by the instance.
     *
     * Some of these nodes may be live, while others may already have been
     * released for reuse.
     */
    int getNumInternalNodeSlots() const { return numInternalNodeSlots_; }
    /**
     * @brief Size of the node domain.
     */
    int getGlobalIdSpaceSize() const { return numInternalNodeSlots_; }
    /**
     * @brief Id of the current hierarchy root.
     */
    NodeId getRoot() const { return rootNodeId_; }
    /**
     * @brief Number of rows in the base image.
     */
    int getNumRowsOfImage() const { return numRows_; }
    /**
     * @brief Explicit alias for the number of columns in the base image.
     */
    int getNumColsOfImage() const { return numCols_; }
    /**
     * @brief Indicates whether the instance represents a max-tree.
     */
    bool isMaxtree() const { return isMaxtree_; }
    /**
     * @brief Adjacency relation used in construction and auxiliary queries.
     */
    AdjacencyRelationPtr getAdjacencyRelation() const { return adj_; }
    /**
     * @brief Altitude/level of an internal node.
     */
    int getAltitude(NodeId nodeId) const { return altitude_[nodeId]; }

    /**
     * @brief Tests whether an index belongs to the node domain.
     *
     * This does not imply that the node is live; for that, use `isAlive`.
     */
    bool isNode(NodeId nodeId) const {
        return nodeId >= 0 && nodeId < getGlobalIdSpaceSize();
    }

    /**
     * @brief Tests whether `nodeId` is the current root of the tree.
     */
    bool isRoot(NodeId nodeId) const {
        return rootNodeId_ == nodeId;
    }

    /**
     * @brief Counts how many nodes are currently live.
     *
     * Unlike `getNumInternalNodeSlots()`, this function ignores nodes that have
     * already been released.
     */
    int getNumNodes() const {
        int count = 0;
        for (NodeId nodeId = 0; nodeId < getGlobalIdSpaceSize(); ++nodeId) {
            if (isAlive(nodeId)) {
                ++count;
            }
        }
        return count;
    }

    /**
     * @brief Tests whether a node is currently live.
     */
    bool isAlive(NodeId nodeId) const {
        return isNode(nodeId) && nodeParent_[nodeId] != InvalidNode;
    }

    /**
     * @brief Direct parent of an internal node.
     */
    NodeId getNodeParent(NodeId nodeId) const {
        return nodeParent_[nodeId];
    }

    /**
     * @brief Number of direct children of a node.
     */
    int getNumChildren(NodeId nodeId) const {
        return numChildrenByNode_[nodeId];
    }

    /**
     * @brief Tests whether a node has no children.
     */
    bool isLeaf(NodeId nodeId) const {
        return firstChild_[nodeId] == InvalidNode;
    }

    /**
     * @brief First direct child of `nodeId`, or `InvalidNode` if none exists.
     */
    NodeId getFirstChild(NodeId nodeId) const {
        return firstChild_[nodeId];
    }

    /**
     * @brief Next sibling of `nodeId` in the parent's child list.
     */
    NodeId getNextSibling(NodeId nodeId) const {
        return nextSibling_[nodeId];
    }

    /**
     * @brief Number of direct proper parts of a node.
     */
    int getNumProperParts(NodeId nodeId) const {
        return numProperPartsByNode_[nodeId];
    }

    /**
     * @brief Node that currently owns the given pixel.
     */
    NodeId getSmallestComponent(PixelId pixelId) const {
        return properPartOwner_[pixelId];
    }

    /**
     * @brief Tests whether `childId` is a direct child of `parentId`.
     */
    bool hasChild(NodeId parentId, NodeId childId) const {
        return getNodeParent(childId) == parentId;
    }

    /**
     * @brief Returns a range for traversing direct children.
     */
    ChildrenRange getChildren(NodeId nodeId) const {
        return ChildrenRange(this, firstChild_[nodeId]);
    }

    /**
     * @brief Returns a range for traversing the node's direct proper parts.
     */
    ProperPartsRange getProperParts(NodeId nodeId) const {
        return ProperPartsRange(this, properHead_[nodeId]);
    }

    /**
     * @brief Returns a range for traversing the subtree rooted at `nodeId`.
     */
    SubtreeNodeRange getNodeSubtree(NodeId nodeId) const {
        return SubtreeNodeRange(this, nodeId);
    }

    /**
     * @brief Lists all pixels that belong to the connected component of `nodeId`.
     */
    std::vector<PixelId> getPixelsOfCC(NodeId nodeId) const {
        std::vector<PixelId> pixels;
        for (NodeId subtreeNodeId : getNodeSubtree(nodeId)) {
            for (PixelId pixelId : getProperParts(subtreeNodeId)) {
                pixels.push_back(pixelId);
            }
        }
        return pixels;
    }

    /**
     * @brief Counts the internal descendant nodes of `nodeId`.
     */
    int countDescendants(NodeId nodeId) const {
        int count = 0;
        for (NodeId subtreeNodeId : getNodeSubtree(nodeId)) {
            if (subtreeNodeId != nodeId) {
                ++count;
            }
        }
        return count;
    }

    /**
     * @brief Counts the proper parts present only in the descendants of `nodeId`.
     */
    int countProperPartsInDescendants(NodeId nodeId) const {
        int count = 0;
        for (NodeId childId : getChildren(nodeId)) {
            count += getNumProperParts(childId);
            count += countProperPartsInDescendants(childId);
        }
        return count;
    }

    /**
     * @brief Moves a single proper part from one node to another.
     * @param targetNodeId Node that will receive the proper part.
     * @param sourceNodeId Node that currently owns the proper part.
     * @param pixelId Proper part/pixel to move.
     */
    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, PixelId pixelId) {
        const PixelId prevPixel = prevProperPart_[pixelId];
        const PixelId nextPixel = nextProperPart_[pixelId];

        if (prevPixel == InvalidNode) {
            properHead_[sourceNodeId] = nextPixel;
        } else {
            nextProperPart_[prevPixel] = nextPixel;
        }

        if (nextPixel == InvalidNode) {
            properTail_[sourceNodeId] = prevPixel;
        } else {
            prevProperPart_[nextPixel] = prevPixel;
        }

        numProperPartsByNode_[sourceNodeId]--;

        prevProperPart_[pixelId] = properTail_[targetNodeId];
        nextProperPart_[pixelId] = InvalidNode;

        if (numProperPartsByNode_[targetNodeId] == 0) {
            properHead_[targetNodeId] = pixelId;
            properTail_[targetNodeId] = pixelId;
            numProperPartsByNode_[targetNodeId] = 1;
        } else {
            nextProperPart_[properTail_[targetNodeId]] = pixelId;
            properTail_[targetNodeId] = pixelId;
            numProperPartsByNode_[targetNodeId] += 1;
        }
        properPartOwner_[pixelId] = targetNodeId;
        properPartVersion_++;
    }

    /**
     * @brief Moves all direct proper parts from `sourceNodeId` to `targetNodeId`.
     */
    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        moveProperPartsInBackend(targetNodeId, sourceNodeId);
        properPartVersion_++;
    }

    /**
     * @brief Sets a new root for the hierarchy.
     *
     * The given node becomes its own parent reference and the previous root,
     * if any, is detached from its old position.
     */
    void setRoot(NodeId nodeId) {
        const NodeId oldRoot = getRoot();
        const NodeId oldParentId = getNodeParent(nodeId);
        if (oldParentId != InvalidNode && oldParentId != nodeId) {
            unlinkChild(nodeId);
        }
        if (oldRoot != InvalidNode && oldRoot != nodeId && isNode(oldRoot) && isAlive(oldRoot)) {
            nodeParent_[oldRoot] = getNodeParent(oldRoot);
        }
        rootNodeId_ = nodeId;
        nodeParent_[nodeId] = nodeId;
        topologyVersion_++;
    }

    /**
     * @brief Reuses a previously released internal id.
     * @return Id of the newly allocated node, or `InvalidNode` if no free ids remain.
     */
    NodeId allocateNode() {
        if (freeNodeIds_.empty()) {
            return InvalidNode;
        }
        const NodeId nodeId = freeNodeIds_.back();
        freeNodeIds_.pop_back();
        nodeParent_[nodeId] = nodeId;
        firstChild_[nodeId] = InvalidNode;
        lastChild_[nodeId] = InvalidNode;
        nextSibling_[nodeId] = InvalidNode;
        prevSibling_[nodeId] = InvalidNode;
        numChildrenByNode_[nodeId] = 0;
        properHead_[nodeId] = InvalidNode;
        properTail_[nodeId] = InvalidNode;
        numProperPartsByNode_[nodeId] = 0;
        altitude_[nodeId] = 0;
        nodeStructureVersion_++;
        return nodeId;
    }

    /**
     * @brief Releases an isolated node with no children and no proper parts.
     */
    void releaseNode(NodeId nodeId) {
        assert(isAlive(nodeId));
        assert(getNodeParent(nodeId) == nodeId);
        assert(firstChild_[nodeId] == InvalidNode);
        assert(numProperPartsByNode_[nodeId] == 0);
        releaseNodeSlot(nodeId);
        nodeStructureVersion_++;
    }

    /**
     * @brief Removes a child from its current parent, optionally releasing the node.
     * @param parentId Expected parent of the node.
     * @param childId Child to remove.
     * @param releaseNodeFlag If `true`, the node slot is returned to the free pool.
     */
    void removeChild(NodeId parentId, NodeId childId, bool releaseNodeFlag) {
        const bool wasDirectChild = isNode(childId) && isAlive(childId) && getNodeParent(childId) == parentId;
        if (!wasDirectChild) {
            return;
        }
        unlinkChild(childId);
        if (releaseNodeFlag) {
            nodeParent_[childId] = InvalidNode;
            altitude_[childId] = 0;
            freeNodeIds_.push_back(childId);
            nodeStructureVersion_++;
        } else {
            nodeParent_[childId] = childId;
        }
        topologyVersion_++;
    }

    /**
     * @brief Attaches `nodeId` as the last child of `parentId`.
     */
    void attachNode(NodeId parentId, NodeId nodeId) {
        linkChildBack(parentId, nodeId);
        nodeParent_[nodeId] = parentId;
        topologyVersion_++;
    }

    /**
     * @brief Detaches a node from its current parent without releasing it.
     */
    void detachNode(NodeId nodeId) {
        unlinkChild(nodeId);
        nodeParent_[nodeId] = nodeId;
        topologyVersion_++;
    }

    /**
     * @brief Moves an entire node to a new parent.
     */
    void moveNode(NodeId nodeId, NodeId newParentId) {
        const NodeId oldParentId = getNodeParent(nodeId);
        if (oldParentId != newParentId) {
            unlinkChild(nodeId);
            linkChildBack(newParentId, nodeId);
        }
        nodeParent_[nodeId] = newParentId;
        topologyVersion_++;
    }

    /**
     * @brief Moves all children of `sourceId` to `parentId`.
     *
     * The relative order of the moved children is preserved in the final linkage.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        const NodeId firstChildId = firstChild_[sourceId];
        if (firstChildId == InvalidNode) {
            return;
        }

        const NodeId lastChildId = lastChild_[sourceId];
        const int movedCount = numChildrenByNode_[sourceId];
        const NodeId tail = lastChild_[parentId];

        firstChild_[sourceId] = InvalidNode;
        lastChild_[sourceId] = InvalidNode;
        numChildrenByNode_[sourceId] = 0;
        prevSibling_[firstChildId] = InvalidNode;
        nextSibling_[lastChildId] = InvalidNode;

        if (tail == InvalidNode) {
            firstChild_[parentId] = firstChildId;
            lastChild_[parentId] = lastChildId;
        } else {
            nextSibling_[tail] = firstChildId;
            prevSibling_[firstChildId] = tail;
            lastChild_[parentId] = lastChildId;
        }
        numChildrenByNode_[parentId] += movedCount;

        for (NodeId childId = firstChildId; childId != InvalidNode; childId = nextSibling_[childId]) {
            nodeParent_[childId] = parentId;
        }
        topologyVersion_++;
    }

    /**
     * @brief Removes an entire subtree, promoting its proper parts to the parent.
     *
     * Each released node has its proper parts transferred directly to the
     * parent of the pruned root before the slot is released.
     */
    void pruneNode(NodeId nodeId) {
        const NodeId parentId = getNodeParent(nodeId);
        auto descendToDeepestLastChild = [this](NodeId startId) {
            NodeId currentId = startId;
            while (!isLeaf(currentId)) {
                currentId = lastChild_[currentId];
            }
            return currentId;
        };

        NodeId currentId = descendToDeepestLastChild(nodeId);
        while (true) {
            if (currentId == nodeId) {
                detachNodeInBackend(nodeId);
                moveProperPartsInBackend(parentId, nodeId);
                releaseNodeSlot(nodeId);
                topologyVersion_++;
                nodeStructureVersion_++;
                properPartVersion_++;
                break;
            }

            const NodeId currentParentId = getNodeParent(currentId);
            detachNodeInBackend(currentId);
            moveProperPartsInBackend(parentId, currentId);
            releaseNodeSlot(currentId);
            topologyVersion_++;
            nodeStructureVersion_++;
            properPartVersion_++;

            currentId = isLeaf(currentParentId) ? currentParentId : descendToDeepestLastChild(currentParentId);
        }
    }

    /**
     * @brief Merges a node into its parent, preserving children and proper parts.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        const NodeId parentId = getNodeParent(nodeId);
        detachNodeInBackend(nodeId);
        moveChildren(parentId, nodeId);
        moveProperPartsInBackend(parentId, nodeId);
        releaseNodeSlot(nodeId);
        topologyVersion_++;
        nodeStructureVersion_++;
        properPartVersion_++;
    }

    /**
     * @brief Reconstructs the current image from the dynamic hierarchy.
     */
    ImageUInt8Ptr reconstructionImage() const {
        auto image = ImageUInt8::create(numRows_, numCols_, 0);
        auto *data = image->rawData();
        for (PixelId pixelId = 0; pixelId < numTotalProperParts_; ++pixelId) {
            const NodeId ownerId = properPartOwner_[pixelId];
            if (ownerId == InvalidNode) {
                continue;
            }
            data[pixelId] = static_cast<uint8_t>(altitude_[ownerId]);
        }
        return image;
    }

    static std::vector<NodeId> getNodesThreshold(DynamicComponentTree *tree,
                                                 int areaThreshold,
                                                 bool enableLog = false) {
        std::vector<NodeId> nodes;
        if (tree == nullptr || tree->getRoot() == InvalidNode) {
            return nodes;
        }

        auto computeArea = [&](auto &&self, NodeId nodeId) -> int {
            int area = tree->getNumProperParts(nodeId);
            for (NodeId childId : tree->getChildren(nodeId)) {
                area += self(self, childId);
            }
            return area;
        };

        FastQueue<NodeId> queue;
        queue.push(tree->getRoot());
        while (!queue.empty()) {
            const NodeId nodeId = queue.pop();
            if (computeArea(computeArea, nodeId) > areaThreshold) {
                for (NodeId childId : tree->getChildren(nodeId)) {
                    queue.push(childId);
                }
            } else {
                nodes.push_back(nodeId);
            }
        }

        if (enableLog) {
            std::cout << "\tArea threshold: " << areaThreshold
                      << ", #InputTreeNodes: " << tree->getNumNodes()
                      << ", #NodesToPrune: " << nodes.size()
                      << std::endl;
        }
        return nodes;
    }

    /**
     * @brief Iterable range of a node's direct children.
     *
     * The class does not materialize copies. It only stores the first child and
     * uses the tree's linked backend to advance.
     */
    class ChildrenRange {
    private:
        const DynamicComponentTree *tree_ = nullptr;
        NodeId firstChildId_ = InvalidNode;

    public:
        /**
         * @brief Forward iterator over direct children.
         *
         * The iterator is stable as long as the topology of the iterated node is not mutated.
         */
        class iterator {
        private:
            const DynamicComponentTree *tree_ = nullptr;
            NodeId currentChildId_ = InvalidNode;

        public:
            iterator(const DynamicComponentTree *tree, NodeId currentChildId) : tree_(tree), currentChildId_(currentChildId) {}

            NodeId operator*() const { return currentChildId_; }

            iterator &operator++() {
                currentChildId_ = (currentChildId_ == InvalidNode) ? InvalidNode : tree_->nextSibling_[currentChildId_];
                return *this;
            }

            bool operator!=(const iterator &other) const {
                return currentChildId_ != other.currentChildId_;
            }
        };

        ChildrenRange(const DynamicComponentTree *tree, NodeId firstChildId) : tree_(tree), firstChildId_(firstChildId) {}

        iterator begin() const { return iterator(tree_, firstChildId_); }
        iterator end() const { return iterator(tree_, InvalidNode); }
    };

    /**
     * @brief Iterable range of a node's direct proper parts.
     *
     * The range exposes the pixels directly owned by the node, without
     * including pixels inherited from descendants.
     */
    class ProperPartsRange {
    private:
        const DynamicComponentTree *tree_ = nullptr;
        PixelId firstPixelId_ = InvalidNode;

    public:
        /**
         * @brief Forward iterator over the node's pixels/proper parts.
         *
         * The order matches the current order of the internal linked list.
         */
        class iterator {
        private:
            const DynamicComponentTree *tree_ = nullptr;
            PixelId currentPixelId_ = InvalidNode;

        public:
            iterator(const DynamicComponentTree *tree, PixelId currentPixelId) : tree_(tree), currentPixelId_(currentPixelId) {}

            PixelId operator*() const { return currentPixelId_; }

            iterator &operator++() {
                currentPixelId_ = (currentPixelId_ == InvalidNode) ? InvalidNode : tree_->nextProperPart_[currentPixelId_];
                return *this;
            }

            bool operator!=(const iterator &other) const {
                return currentPixelId_ != other.currentPixelId_;
            }
        };

        ProperPartsRange(const DynamicComponentTree *tree, PixelId firstPixelId) : tree_(tree), firstPixelId_(firstPixelId) {}

        iterator begin() const { return iterator(tree_, firstPixelId_); }
        iterator end() const { return iterator(tree_, InvalidNode); }
    };

    /**
     * @brief Range for depth-first traversal of a subtree.
     *
     * The subtree root is included in the traversal. Children are pushed in
     * reverse order to preserve the natural linkage order when visiting.
     */
    class SubtreeNodeRange {
    private:
        const DynamicComponentTree *tree_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;

    public:
        /**
         * @brief Non-recursive DFS iterator over internal nodes.
         *
         * This iterator avoids explicit recursion and traverses only live
         * internal nodes reachable from the root provided to the range.
         */
        class iterator {
        private:
            const DynamicComponentTree *tree_ = nullptr;
            FastStack<int> stack_;
            NodeId currentNodeId_ = InvalidNode;

            void advance() {
                if (stack_.empty()) {
                    currentNodeId_ = InvalidNode;
                    return;
                }
                currentNodeId_ = stack_.pop();
                for (NodeId childId = tree_->lastChild_[currentNodeId_];
                     childId != InvalidNode;
                     childId = tree_->prevSibling_[childId]) {
                    stack_.push(childId);
                }
            }

        public:
            iterator(const DynamicComponentTree *tree, NodeId rootNodeId, bool isEnd) : tree_(tree) {
                if (!isEnd && tree_ != nullptr && rootNodeId != InvalidNode) {
                    stack_.push(rootNodeId);
                    advance();
                }
            }

            NodeId operator*() const { return currentNodeId_; }

            iterator &operator++() {
                advance();
                return *this;
            }

            bool operator!=(const iterator &other) const {
                return currentNodeId_ != other.currentNodeId_;
            }
        };

        SubtreeNodeRange(const DynamicComponentTree *tree, NodeId rootNodeId) : tree_(tree), rootNodeId_(rootNodeId) {}

        iterator begin() const { return iterator(tree_, rootNodeId_, false); }
        iterator end() const { return iterator(tree_, InvalidNode, true); }
    };

    class BreadthFirstNodeRange {
    private:
        const DynamicComponentTree *tree_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;

    public:
        /**
         * @brief Breadth-first forward iterator over internal nodes.
         */
        class iterator {
        private:
            const DynamicComponentTree *tree_ = nullptr;
            FastQueue<NodeId> queue_;

        public:
            iterator(const DynamicComponentTree *tree, NodeId rootNodeId, bool isEnd) : tree_(tree) {
                if (!isEnd && tree_ != nullptr && rootNodeId != InvalidNode) {
                    queue_.push(rootNodeId);
                }
            }

            NodeId operator*() const {
                return queue_.front();
            }

            iterator &operator++() {
                if (!queue_.empty()) {
                    const NodeId nodeId = queue_.pop();
                    for (NodeId childId : tree_->getChildren(nodeId)) {
                        queue_.push(childId);
                    }
                }
                return *this;
            }

            bool operator!=(const iterator &other) const {
                return queue_.empty() != other.queue_.empty();
            }
        };

        /**
         * @brief Builds the BFS range rooted at `rootNodeId`.
         */
        BreadthFirstNodeRange(const DynamicComponentTree *tree, NodeId rootNodeId) : tree_(tree), rootNodeId_(rootNodeId) {}

        iterator begin() const { return iterator(tree_, rootNodeId_, false); }
        iterator end() const { return iterator(tree_, InvalidNode, true); }
    };

    /**
     * @brief Returns a range for breadth-first traversal from the current root.
     */
    BreadthFirstNodeRange getIteratorBreadthFirstTraversal() const {
        return BreadthFirstNodeRange(this, getRoot());
    }

    /**
     * @brief Returns a range for breadth-first traversal from `rootNodeId`.
     */
    BreadthFirstNodeRange getIteratorBreadthFirstTraversal(NodeId rootNodeId) const {
        return BreadthFirstNodeRange(this, rootNodeId);
    }
};
