#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "Common.hpp"
#include "DynamicComponentTree.hpp"

enum class Attribute {
    AREA,
    BOX_WIDTH,
    BOX_HEIGHT,
    DIAGONAL_LENGTH,
};
using enum Attribute;

struct DynamicBoundingBoxState {
    int xmin = 0;
    int xmax = -1;
    int ymin = 0;
    int ymax = -1;
    bool empty = true;
};

/**
 * @brief Common interface for incremental attribute computers on dynamic trees.
 *
 * A local node update always follows the same protocol:
 * - initialize the local state in `preProcessing`;
 * - merge child contributions through `mergeProcessing`;
 * - materialize the final scalar value in `postProcessing`.
 *
 * `DynamicComponentTreeAdjustment` depends only on this interface, which makes
 * it possible to swap the incremental attribute without changing the
 * structural-adjustment logic.
 */
class DynamicAttributeComputer {
public:
    virtual ~DynamicAttributeComputer() = default;

    /**
     * @brief Initializes the direct contribution of `nodeId` in `buffer`.
     * @details This phase prepares the node's local state before any merge with
     * children. For simple attributes, this usually means writing only the
     * node's direct contribution into `buffer`. For attributes with persistent
     * summaries, this step may also synchronize internal caches before
     * structural aggregation.
     */
    virtual void preProcessing(NodeId nodeId, std::span<float> buffer) const = 0;

    /**
     * @brief Aggregates into `parentId` the already consolidated contribution of `childId`.
     * @details The contract assumes that `preProcessing(childId, ...)` and the
     * full recomputation of the child subtree have already happened. This
     * method therefore only needs to merge the child state into the parent
     * accumulator.
     */
    virtual void mergeProcessing(NodeId parentId, NodeId childId, std::span<float> buffer) const = 0;

    /**
     * @brief Materializes the final attribute value of `nodeId` in `buffer`.
     * @details This phase closes the node recomputation after all child
     * contributions have already been aggregated. For attributes whose
     * authoritative value lives directly in `buffer`, this method can be a
     * no-op.
     */
    virtual void postProcessing(NodeId nodeId, std::span<float> buffer) const = 0;

    /**
     * @brief Notifies a full merge of the direct proper parts of `sourceId` into `targetId`.
     * @details The default implementation does nothing. Attributes that keep
     * additional incremental state can override this hook to avoid reprocessing
     * pixels after `moveProperParts(targetId, sourceId)`.
     */
    virtual void onMoveProperParts(NodeId, NodeId) const {
        /* no-op */
    }

    /**
     * @brief Notifies the transfer of a single proper part from `sourceId` to `targetId`.
     * @details The default implementation does nothing. Geometric attributes
     * can use this hook to incrementally update the local state of the receiver
     * node and mark the donor as locally inconsistent.
     */
    virtual void onMoveProperPart(NodeId, NodeId, PixelId) const {
        /* no-op */
    }

    /**
     * @brief Notifies the definitive removal of a node from the hierarchy.
     * @details The default implementation does nothing. Attributes with
     * persistent summaries can use this hook to discard auxiliary state.
     */
    virtual void onNodeRemoved(NodeId) const {
        /* no-op */
    }

};

/**
 * @brief Generic infrastructure for attributes with persistent per-node summaries.
 *
 * @details This base explicitly separates two summary levels:
 * - `LocalSummary`: state associated only with the node's direct proper parts;
 * - `SubtreeSummary`: aggregated state of the node subtree after merging children.
 *
 * The `Policy` template parameter defines how those summaries are built,
 * combined, and converted into the final scalar attribute value. The adjuster
 * sees only the `DynamicAttributeComputer` interface; all attribute-specific
 * logic remains encapsulated in the policy.
 */
template<class Policy>
class DynamicSummaryComputer : public DynamicAttributeComputer {
public:
    using LocalSummary = typename Policy::LocalSummary;
    using SubtreeSummary = typename Policy::SubtreeSummary;

private:
    DynamicComponentTree *tree_ = nullptr;
    Policy policy_;
    mutable std::vector<LocalSummary> local_;
    mutable std::vector<SubtreeSummary> subtree_;

    /**
     * @brief Recursively recomputes the subtree of `nodeId` and writes the result into `buffer`.
     * @details This helper is used for full computations outside the adjuster's
     * incremental flow, for example when creating reference buffers or
     * explicitly bootstrapping the attribute.
     */
    void computeNode(NodeId nodeId, std::span<float> buffer) const {
        preProcessing(nodeId, buffer);
        for (NodeId childId : tree_->getChildren(nodeId)) {
            computeNode(childId, buffer);
            mergeProcessing(nodeId, childId, buffer);
        }
        postProcessing(nodeId, buffer);
    }

    /**
     * @brief Initializes the summaries of `nodeId` and of its entire subtree.
     * @details The routine first ensures the local summary of the current node
     * and then propagates that state into the subtree summary. It then visits
     * the children and incorporates their already initialized summaries into
     * the parent node.
     */
    void initializeSummariesFromNode(NodeId nodeId) const {
        policy_.ensureLocal(tree_, nodeId, local_[(size_t) nodeId]);
        policy_.copyLocalToSubtree(local_[(size_t) nodeId], subtree_[(size_t) nodeId]);
        for (NodeId childId : tree_->getChildren(nodeId)) {
            initializeSummariesFromNode(childId);
            policy_.mergeSubtree(subtree_[(size_t) nodeId], subtree_[(size_t) childId]);
        }
    }

    /**
     * @brief Bootstraps the persistent summaries from the current tree root.
     * @details The method is called only once in the constructor. The premise
     * is that the tree already represents a consistent state and that the
     * computer only needs to build its internal runtime to track future
     * updates.
     */
    void initializeSummaries() const {
        if (tree_ == nullptr || tree_->getRoot() == InvalidNode) {
            return;
        }
        initializeSummariesFromNode(tree_->getRoot());
    }

protected:
    /**
     * @brief Builds the generic infrastructure and initializes the persistent summaries.
     * @param tree Dynamic tree observed by the computer.
     * @param policy Policy that implements the attribute semantics.
     */
    DynamicSummaryComputer(DynamicComponentTree *tree, Policy policy)
        : tree_(tree),
          policy_(std::move(policy)),
          local_((size_t) (tree ? tree->getGlobalIdSpaceSize() : 0)),
          subtree_((size_t) (tree ? tree->getGlobalIdSpaceSize() : 0)) {
        initializeSummaries();
    }

    /**
     * @brief Exposes the tree associated with the derived computer.
     * @details Used only by concrete subclasses that need to clear or inspect
     * the persistent state of a specific node.
     */
    DynamicComponentTree *tree() const {
        return tree_;
    }

    /**
     * @brief Accesses the persistent local summary of `nodeId`.
     * @details This helper avoids exposing the internal vectors directly and
     * makes it clear in derived code that the access is to the summary of the
     * node's direct proper parts.
     */
    LocalSummary &localSummary(NodeId nodeId) const {
        return local_[(size_t) nodeId];
    }

    /**
     * @brief Accesses the persistent subtree summary of `nodeId`.
     * @details The returned value represents the aggregated state after merging
     * children, not just the node's direct contribution.
     */
    SubtreeSummary &subtreeSummary(NodeId nodeId) const {
        return subtree_[(size_t) nodeId];
    }

public:
    /**
     * @brief Generic base for attributes whose authoritative state lives in internal summaries.
     * @details This infrastructure is appropriate for attribute families where
     * the cost of a persistent summary compensates for the extra indirection on
     * the hot path, such as `bbox_*`. It should not be used indiscriminately
     * for simple additive attributes like `AREA`, whose direct path based only
     * on the external buffer is cheaper.
     */
    void preProcessing(NodeId nodeId, std::span<float>) const override {
        policy_.ensureLocal(tree_, nodeId, local_[(size_t) nodeId]);
        policy_.copyLocalToSubtree(local_[(size_t) nodeId], subtree_[(size_t) nodeId]);
    }

    /**
     * @brief Merges into the parent's subtree summary the already consolidated child contribution.
     * @details The merge is fully delegated to `Policy`, which allows the
     * generic infrastructure to remain independent of the attribute semantics.
     */
    void mergeProcessing(NodeId parentId, NodeId childId, std::span<float>) const override {
        policy_.mergeSubtree(subtree_[(size_t) parentId], subtree_[(size_t) childId]);
    }

    /**
     * @brief Writes into `buffer` the final scalar value derived from the subtree summary.
     * @details `buffer` is treated as the external view of the attribute, while
     * the internal state remains in the `local_` and `subtree_` vectors.
     */
    void postProcessing(NodeId nodeId, std::span<float> buffer) const override {
        buffer[(size_t) nodeId] = policy_.value(subtree_[(size_t) nodeId]);
    }

    /**
     * @brief Updates the local summaries after a full proper-part merge.
     * @details This operation is used when the entire local set of `sourceId`
     * is transferred to `targetId`, allowing a cheaper merge than reprocessing
     * pixel by pixel.
     */
    void onMoveProperParts(NodeId targetId, NodeId sourceId) const override {
        policy_.moveLocalSummary(tree_,
                                 targetId,
                                 sourceId,
                                 local_[(size_t) targetId],
                                 local_[(size_t) sourceId]);
    }

    /**
     * @brief Updates the local summaries after transferring a single pixel.
     * @details The policy decides how to incorporate the pixel into the
     * receiver and how to handle the possible invalidation of the donor's local
     * summary.
     */
    void onMoveProperPart(NodeId targetId, NodeId sourceId, PixelId pixelId) const override {
        policy_.movePixel(tree_,
                          targetId,
                          sourceId,
                          pixelId,
                          local_[(size_t) targetId],
                          local_[(size_t) sourceId]);
    }

    /**
     * @brief Computes the full attribute and returns a new buffer with the values.
     * @details Convenience method used in benchmarks, tools, and tests when the
     * entire attribute needs to be recomputed without relying on the
     * incremental flow.
     */
    std::vector<float> compute() const {
        std::vector<float> buffer((size_t) (tree_ ? tree_->getGlobalIdSpaceSize() : 0), 0.0f);
        compute(std::span<float>(buffer));
        return buffer;
    }

    /**
     * @brief Recomputes the full tree attribute over a caller-provided `buffer`.
     * @details Unlike the adjuster's incremental flow, this routine visits the
     * entire tree starting from the current root.
     */
    void compute(std::span<float> buffer) const {
        if (tree_ == nullptr || tree_->getRoot() == InvalidNode) {
            return;
        }
        computeNode(tree_->getRoot(), buffer);
    }
};

/**
 * @brief Incremental computer for the area attribute on the pixel-based dynamic tree.
 *
 * The area of a node is defined as the sum of:
 * - the node's direct proper parts;
 * - the areas of all children in the subtree.
 *
 * The class also exposes the `preProcessing` / `mergeProcessing` /
 * `postProcessing` interface used by `DynamicComponentTreeAdjustment` to
 * recompute only the nodes affected by a local mutation.
 *
 * The implementation remains intentionally explicit: for `AREA`, the
 * incremental value is already handled directly in the external buffer and
 * does not benefit from the persistent-summary infrastructure used by
 * `bbox_*`.
 */
class DynamicAreaComputer : public DynamicAttributeComputer {
private:
    DynamicComponentTree *tree_ = nullptr;

    /**
     * @brief Recursively recomputes the accumulated area of the subtree of `nodeId`.
     * @details This helper provides the full version of the area computation,
     * independent of the adjuster's incremental flow.
     */
    void computeNodeArea(NodeId nodeId, std::span<float> buffer) {
        preProcessing(nodeId, buffer);
        for (NodeId childId : tree_->getChildren(nodeId)) {
            computeNodeArea(childId, buffer);
            mergeProcessing(nodeId, childId, buffer);
        }
        postProcessing(nodeId, buffer);
    }

public:
    /**
     * @brief Builds the incremental area computer for a dynamic tree.
     */
    explicit DynamicAreaComputer(DynamicComponentTree *tree) : tree_(tree) {}

    /**
     * @brief Initializes the local contribution of `nodeId` with the number of direct proper parts.
     * @details For `AREA`, the local node value is exactly the cardinality of
     * its set of direct proper parts.
     */
    void preProcessing(NodeId nodeId, std::span<float> buffer) const override {
        buffer[(size_t) nodeId] = static_cast<float>(tree_->getNumProperParts(nodeId));
    }

    /**
     * @brief Adds the child's already consolidated area to the parent.
     * @details `AREA` aggregation is purely additive, with no need for
     * auxiliary summaries beyond `buffer` itself.
     */
    void mergeProcessing(NodeId parentId, NodeId childId, std::span<float> buffer) const override {
        buffer[(size_t) parentId] += buffer[(size_t) childId];
    }

    /**
     * @brief Empty finalization step for the `AREA` incremental protocol.
     * @details All attribute logic has already been resolved in
     * `preProcessing` and `mergeProcessing`, so there is no extra
     * post-processing.
     */
    void postProcessing(NodeId, std::span<float>) const override {
        /* no-op */
    }

    /**
     * @brief Computes the full tree area and returns a new buffer.
     */
    std::vector<float> compute() {
        std::vector<float> buffer((size_t) tree_->getGlobalIdSpaceSize(), 0.0f);
        compute(std::span<float>(buffer));
        return buffer;
    }

    /**
     * @brief Computes the full tree area on a caller-provided buffer.
     */
    void compute(std::span<float> buffer) {
        if (tree_ == nullptr || tree_->getRoot() == InvalidNode) {
            return;
        }
        computeNodeArea(tree_->getRoot(), buffer);
    }
};

/**
 * @brief Local summary used by the `bbox_*` family.
 *
 * @details The local summary describes only the node's direct proper parts. In
 * addition to the four rectangle extremes, it stores counts on the current
 * borders to detect when a removal makes the summary inaccurate and requires a
 * rebuild.
 */
struct DynamicBoundingBoxLocalSummary {
    int xmin = 0;
    int xmax = -1;
    int ymin = 0;
    int ymax = -1;
    int xminCount = 0;
    int xmaxCount = 0;
    int yminCount = 0;
    int ymaxCount = 0;
    bool empty = true;
    bool dirty = true;
};

/**
 * @brief Policy that implements the incremental semantics of the `bbox_*` family.
 *
 * @details This policy defines:
 * - how to maintain the local bounding box of direct proper parts;
 * - how to derive the subtree summary from the local state and the children;
 * - how to convert the final summary into `width`, `height`, or `diagonal`.
 *
 * The current strategy uses a lightweight local summary with border counts and
 * a `dirty` fallback when removing a pixel exhausts some extremum.
 */
class DynamicBoundingBoxPolicy {
private:
    Attribute attribute_ = AREA;
    int numCols_ = 0;
    int numRows_ = 0;

    /**
     * @brief Places the local summary in the canonical empty state.
     * @details The empty state uses image-compatible sentinels so that the
     * first pixel expansion rewrites every extremum correctly.
     */
    void resetLocalBox(DynamicBoundingBoxLocalSummary &local) const {
        local.xmin = numCols_;
        local.xmax = -1;
        local.ymin = numRows_;
        local.ymax = -1;
        local.xminCount = 0;
        local.xmaxCount = 0;
        local.yminCount = 0;
        local.ymaxCount = 0;
        local.empty = true;
    }

    /**
     * @brief Expands the local bounding box with a single pixel.
     * @details In addition to updating the four extrema, the method keeps the
     * counts of the pixels that support each current border.
     */
    void expandLocalBoxWithPixel(DynamicBoundingBoxLocalSummary &local, PixelId pixelId) const {
        auto [row, col] = ImageUtils::to2D(pixelId, numCols_);
        if (local.empty) {
            local.xmin = col;
            local.xmax = col;
            local.ymin = row;
            local.ymax = row;
            local.xminCount = 1;
            local.xmaxCount = 1;
            local.yminCount = 1;
            local.ymaxCount = 1;
            local.empty = false;
            return;
        }

        if (col < local.xmin) {
            local.xmin = col;
            local.xminCount = 1;
        } else if (col == local.xmin) {
            ++local.xminCount;
        }

        if (col > local.xmax) {
            local.xmax = col;
            local.xmaxCount = 1;
        } else if (col == local.xmax) {
            ++local.xmaxCount;
        }

        if (row < local.ymin) {
            local.ymin = row;
            local.yminCount = 1;
        } else if (row == local.ymin) {
            ++local.yminCount;
        }

        if (row > local.ymax) {
            local.ymax = row;
            local.ymaxCount = 1;
        } else if (row == local.ymax) {
            ++local.ymaxCount;
        }
    }

    /**
     * @brief Unions two complete local summaries.
     * @details This operation is used when all local proper parts of one node
     * are moved to another. Border counts are preserved whenever the extrema
     * coincide.
     */
    void mergeLocalBoxes(DynamicBoundingBoxLocalSummary &target,
                         const DynamicBoundingBoxLocalSummary &source) const {
        if (source.empty) {
            return;
        }

        if (target.empty) {
            target = source;
            return;
        }

        if (source.xmin < target.xmin) {
            target.xmin = source.xmin;
            target.xminCount = source.xminCount;
        } else if (source.xmin == target.xmin) {
            target.xminCount += source.xminCount;
        }

        if (source.xmax > target.xmax) {
            target.xmax = source.xmax;
            target.xmaxCount = source.xmaxCount;
        } else if (source.xmax == target.xmax) {
            target.xmaxCount += source.xmaxCount;
        }

        if (source.ymin < target.ymin) {
            target.ymin = source.ymin;
            target.yminCount = source.yminCount;
        } else if (source.ymin == target.ymin) {
            target.yminCount += source.yminCount;
        }

        if (source.ymax > target.ymax) {
            target.ymax = source.ymax;
            target.ymaxCount = source.ymaxCount;
        } else if (source.ymax == target.ymax) {
            target.ymaxCount += source.ymaxCount;
        }
    }

    /**
     * @brief Rebuilds the local summary of `nodeId` from scratch.
     * @details This is the fallback used when the node enters the `dirty`
     * state, that is, when a removal makes the current incremental summary
     * insufficient.
     */
    void rebuildLocalBox(DynamicComponentTree *tree,
                         NodeId nodeId,
                         DynamicBoundingBoxLocalSummary &local) const {
        resetLocalBox(local);
        for (PixelId pixelId : tree->getProperParts(nodeId)) {
            expandLocalBoxWithPixel(local, pixelId);
        }
        local.dirty = false;
    }

    /**
     * @brief Unions two subtree summaries.
     * @details For `bbox_*`, subtree merging is simply the union of the minimum
     * rectangles of the nodes involved.
     */
    void mergeSubtreeStates(DynamicBoundingBoxState &state, const DynamicBoundingBoxState &other) const {
        if (other.empty) {
            return;
        }

        if (state.empty) {
            state = other;
            return;
        }

        state.xmin = std::min(state.xmin, other.xmin);
        state.xmax = std::max(state.xmax, other.xmax);
        state.ymin = std::min(state.ymin, other.ymin);
        state.ymax = std::max(state.ymax, other.ymax);
        state.empty = false;
    }

public:
    using LocalSummary = DynamicBoundingBoxLocalSummary;
    using SubtreeSummary = DynamicBoundingBoxState;

    /**
     * @brief Default constructor required by the generic infrastructure.
     */
    DynamicBoundingBoxPolicy() = default;

    /**
     * @brief Builds the policy for a tree and a specific geometric measure.
     * @param tree Tree whose geometry defines the image domain.
     * @param attribute Derived measure to produce: width, height, or diagonal.
     */
    DynamicBoundingBoxPolicy(DynamicComponentTree *tree, Attribute attribute)
        : attribute_(attribute),
          numCols_(tree ? tree->getNumColsOfImage() : 0),
          numRows_(tree ? tree->getNumRowsOfImage() : 0) {}

    /**
     * @brief Ensures that the local summary of `nodeId` is consistent.
     * @details If the node is marked as `dirty`, it performs a full rebuild
     * from the current direct proper parts; otherwise, it reuses the cache.
     */
    void ensureLocal(DynamicComponentTree *tree,
                     NodeId nodeId,
                     LocalSummary &local) const {
        if (nodeId == InvalidNode || tree == nullptr || !tree->isNode(nodeId)) {
            return;
        }
        if (!local.dirty) {
            return;
        }
        rebuildLocalBox(tree, nodeId, local);
    }

    /**
     * @brief Copies the local summary into the initial subtree summary.
     * @details The subtree always starts with the node's direct contribution
     * before incorporating the contribution of the children.
     */
    void copyLocalToSubtree(const LocalSummary &local, SubtreeSummary &subtree) const {
        subtree.xmin = local.xmin;
        subtree.xmax = local.xmax;
        subtree.ymin = local.ymin;
        subtree.ymax = local.ymax;
        subtree.empty = local.empty;
    }

    /**
     * @brief Aggregates into the parent's subtree summary the child's already consolidated contribution.
     */
    void mergeSubtree(SubtreeSummary &parent, const SubtreeSummary &child) const {
        mergeSubtreeStates(parent, child);
    }

    /**
     * @brief Converts the final subtree summary into the scalar attribute value.
     * @details The same policy serves `bbox_width`, `bbox_height`, and
     * `bbox_diagonal`; only the final projection on the rectangle changes.
     */
    float value(const SubtreeSummary &subtree) const {
        if (subtree.empty) {
            return 0.0f;
        }

        const float width = static_cast<float>(subtree.xmax - subtree.xmin + 1);
        const float height = static_cast<float>(subtree.ymax - subtree.ymin + 1);
        if (attribute_ == BOX_WIDTH) {
            return width;
        }
        if (attribute_ == BOX_HEIGHT) {
            return height;
        }
        return std::sqrt(width * width + height * height);
    }

    /**
     * @brief Updates the local summaries after a full proper-part merge.
     * @details The destination receives the union of the local rectangles and
     * the source is zeroed because it no longer owns direct proper parts.
     */
    void moveLocalSummary(DynamicComponentTree *tree,
                          NodeId targetId,
                          NodeId sourceId,
                          LocalSummary &target,
                          LocalSummary &source) const {
        if (tree == nullptr || targetId == InvalidNode || sourceId == InvalidNode ||
            !tree->isNode(targetId) || !tree->isNode(sourceId)) {
            return;
        }

        ensureLocal(tree, targetId, target);
        ensureLocal(tree, sourceId, source);
        mergeLocalBoxes(target, source);
        resetLocalBox(source);
        source.dirty = false;
    }

    /**
     * @brief Incrementally updates the local summaries after moving a single pixel.
     * @details The receiver can always be updated exactly. The donor, however,
     * may enter `dirty` if the removed pixel was the last one supporting some
     * local extremum and the next extremum cannot be inferred without a rebuild.
     */
    void movePixel(DynamicComponentTree *tree,
                   NodeId targetId,
                   NodeId sourceId,
                   PixelId pixelId,
                   LocalSummary &target,
                   LocalSummary &source) const {
        if (tree == nullptr || targetId == InvalidNode || sourceId == InvalidNode ||
            !tree->isNode(targetId) || !tree->isNode(sourceId)) {
            return;
        }

        ensureLocal(tree, targetId, target);
        expandLocalBoxWithPixel(target, pixelId);
        target.dirty = false;

        if (!tree->isAlive(sourceId) || tree->getNumProperParts(sourceId) == 0) {
            resetLocalBox(source);
            source.dirty = false;
        } else if (!source.dirty) {
            const auto [row, col] = ImageUtils::to2D(pixelId, numCols_);
            bool exhaustsXmin = false;
            bool exhaustsXmax = false;
            bool exhaustsYmin = false;
            bool exhaustsYmax = false;

            if (col == source.xmin) {
                if (source.xminCount <= 1) {
                    exhaustsXmin = true;
                } else {
                    --source.xminCount;
                }
            }

            if (col == source.xmax) {
                if (source.xmaxCount <= 1) {
                    exhaustsXmax = true;
                } else {
                    --source.xmaxCount;
                }
            }

            if (row == source.ymin) {
                if (source.yminCount <= 1) {
                    exhaustsYmin = true;
                } else {
                    --source.yminCount;
                }
            }

            if (row == source.ymax) {
                if (source.ymaxCount <= 1) {
                    exhaustsYmax = true;
                } else {
                    --source.ymaxCount;
                }
            }

            const bool needsDirty = exhaustsXmin || exhaustsXmax || exhaustsYmin || exhaustsYmax;
            if (needsDirty) {
                source.dirty = true;
            }
        } else {
            source.dirty = true;
        }
    }
};

/**
 * @brief Computer for geometric attributes derived from the bounding box of each node.
 *
 * The class traverses the tree and maintains, for each node, the bounding box
 * that covers all proper parts in its subtree. From that box, it produces one
 * of the attributes listed in `Attribute`:
 * - width;
 * - height;
 * - diagonal.
 */
class DynamicBoundingBoxComputer : public DynamicSummaryComputer<DynamicBoundingBoxPolicy> {
private:
    using Base = DynamicSummaryComputer<DynamicBoundingBoxPolicy>;

    /**
     * @brief Clears the persistent local summary of a removed node.
     * @details After the node slot is returned to the tree, the local attribute
     * cache must also return to the empty state to avoid contamination if the
     * `NodeId` is reused in the future.
     */
    void resetLocalSummary(NodeId nodeId) const {
        auto &local = localSummary(nodeId);
        local.xmin = 0;
        local.xmax = -1;
        local.ymin = 0;
        local.ymax = -1;
        local.xminCount = 0;
        local.xmaxCount = 0;
        local.yminCount = 0;
        local.ymaxCount = 0;
        local.empty = true;
        local.dirty = false;
    }

    /**
     * @brief Clears the persistent subtree summary of a removed node.
     * @details This reset guarantees that a recycled `NodeId` does not inherit
     * the rectangle of an old subtree.
     */
    void resetSubtreeSummary(NodeId nodeId) const {
        auto &subtree = subtreeSummary(nodeId);
        subtree.xmin = 0;
        subtree.xmax = -1;
        subtree.ymin = 0;
        subtree.ymax = -1;
        subtree.empty = true;
    }

public:
    /**
     * @brief Builds the incremental computer for the `bbox_*` family.
     * @param tree Observed dynamic tree.
     * @param attribute Desired derived geometric measure.
     */
    DynamicBoundingBoxComputer(DynamicComponentTree *tree, Attribute attribute)
        : Base(tree, DynamicBoundingBoxPolicy{tree, attribute}) {}

    /**
     * @brief Discards the persistent state associated with a removed node.
     * @details The method is called by the adjuster when a node ceases to
     * exist in the dynamic hierarchy. The cleanup keeps the summary vectors
     * coherent if the same `NodeId` is allocated again in the future.
     */
    void onNodeRemoved(NodeId nodeId) const override {
        if (tree() == nullptr || nodeId == InvalidNode) {
            return;
        }
        resetLocalSummary(nodeId);
        resetSubtreeSummary(nodeId);
    }
};
