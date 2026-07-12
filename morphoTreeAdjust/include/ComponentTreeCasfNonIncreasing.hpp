#pragma once

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "DynamicComponentTree.hpp"
#include "DualMinMaxTreeIncrementalFilter.hpp"

/**
 * @brief Local-merge CASF for supported non-order-preserving flat attributes.
 *
 * This class accepts `AREA_PER_BOX_PERIMETER`, `RECTANGULARITY`,
 * `BOX_ASPECT_RATIO`, `BOX_ELONGATION`, and
 * `AREA_PER_BOX_DIAGONAL_SQUARED`. Unlike conventional subtree-pruning CASF,
 * this path contracts each selected node locally into its parent and uses
 * `DualMinMaxTreeIncrementalFilter` to keep the paired tree synchronized.
 *
 * Candidate selection takes one deterministic snapshot at the beginning of
 * each directional phase. A non-root node is selected when its value in that
 * snapshot is less than or equal to the threshold; `NaN` is not selected.
 * Candidates invalidated by an earlier contraction in the same phase are
 * skipped, but surviving candidates are not threshold-tested again if their
 * attribute changes. Every threshold processes the max-tree phase before
 * taking the min-tree snapshot and processing the min-tree phase.
 *
 * Thresholds and repeated `filter()` calls are cumulative. This class neither
 * exposes a reset operation nor rebuilds the pair between thresholds.
 *
 * @tparam PixelType Unsigned integral altitude type used by the incremental
 * adjuster. Its level domain must start at zero.
 *
 * @warning Mutating either tree through a non-const getter bypasses paired-tree
 * and attribute maintenance and can desynchronize later operations.
 */
template<typename PixelType = AltitudeType>
class ComponentTreeCasfNonIncreasing {
private:
    /**
     * @brief Stable ordering signature for one node selected by a threshold.
     */
    struct CandidateOrderEntry {
        NodeId nodeId = InvalidNode; ///< Candidate node in the current primal tree.
        int altitude = 0; ///< Node altitude, used as the first deterministic ordering key.
        PixelId supportMinPixel = std::numeric_limits<PixelId>::max(); ///< Smallest support pixel.
        int supportSize = 0; ///< Candidate support size, used after altitude and support pixel.
    };

    /**
     * @brief Reusable buffers used while collecting merge candidates.
     *
     * The buffers are kept as object state to avoid reallocating vectors on
     * every threshold iteration.
     */
    struct CandidateScratch {
        FastQueue<NodeId> queue; ///< Breadth-first traversal queue over the current tree.
        std::vector<CandidateOrderEntry> ordered; ///< Candidate entries before projection to node ids.
        std::vector<NodeId> output; ///< Sorted candidate node ids.
    };

    AdjacencyRelationPtr adjacency_; ///< Adjacency shared by the trees and adjuster.
    Attribute attribute_ = AREA_PER_BOX_PERIMETER; ///< Attribute used to select merge candidates.
    std::unique_ptr<DynamicComponentTree> maxtree_; ///< Dynamic max-tree for upper level sets.
    std::unique_ptr<DynamicComponentTree> mintree_; ///< Dynamic min-tree for lower level sets.
    std::unique_ptr<DynamicBoundingBoxComputer> maxAttributeComputer_; ///< Max-tree attribute computer.
    std::unique_ptr<DynamicBoundingBoxComputer> minAttributeComputer_; ///< Min-tree attribute computer.
    std::vector<float> maxAttribute_; ///< External attribute buffer for `maxtree_`.
    std::vector<float> minAttribute_; ///< External attribute buffer for `mintree_`.
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<PixelType>> adjust_; ///< Primal/dual update engine.
    CandidateScratch candidateScratch_; ///< Scratch buffers reused while scanning thresholds.

    /**
     * @brief Returns whether `attribute` is accepted by this public CASF path.
     */
    static bool isSupportedAttribute(Attribute attribute) {
        switch (attribute) {
            case AREA_PER_BOX_PERIMETER:
            case RECTANGULARITY:
            case BOX_ASPECT_RATIO:
            case BOX_ELONGATION:
            case AREA_PER_BOX_DIAGONAL_SQUARED:
                return true;
            default:
                return false;
        }
    }

    /**
     * @brief Throws when `attribute` is not a supported non-order-preserving attribute.
     */
    static void requireSupportedAttribute(Attribute attribute) {
        if (!isSupportedAttribute(attribute)) {
            throw std::runtime_error(
                "ComponentTreeCasfNonIncreasing supports only AREA_PER_BOX_PERIMETER, "
                "RECTANGULARITY, BOX_ASPECT_RATIO, BOX_ELONGATION, and "
                "AREA_PER_BOX_DIAGONAL_SQUARED.");
        }
    }

    /**
     * @brief Calculates an attribute-buffer capacity safe for one local merge.
     *
     * Dynamic contractions can temporarily allocate ids beyond the current
     * internal node slots. Reserving the pixel count margin keeps the span used
     * by the adjuster valid during the merge.
     */
    static std::size_t bufferCapacityForSingleMerge(const DynamicComponentTree &tree) {
        const auto slots = static_cast<std::size_t>(
            std::max(0, tree.getNumInternalNodeSlots()));
        const auto pixels = static_cast<std::size_t>(
            std::max(1, tree.getNumTotalProperParts()));
        return slots + pixels;
    }

    /**
     * @brief Clears only the active node slots in an attribute buffer.
     */
    static void clearActiveAttributeSlots(const DynamicComponentTree &tree,
                                          std::vector<float> &attribute) {
        const auto slots = static_cast<std::size_t>(
            std::max(0, tree.getNumInternalNodeSlots()));
        if (attribute.size() < slots) {
            attribute.resize(slots, 0.0f);
        }
        std::fill_n(attribute.begin(), slots, 0.0f);
    }

    /**
     * @brief Snapshots nodes whose current attribute value is at or below `threshold`.
     *
     * The returned vector is owned by `work` and remains valid until the next
     * call. Candidates are sorted by altitude and by support signatures so the
     * merge order is stable even when child traversal order changes. A `NaN`
     * attribute value compares false and is not included.
     */
    static std::vector<NodeId> &collectMergeCandidates(const DynamicComponentTree &tree,
                                                       const std::vector<float> &attribute,
                                                       const DynamicBoundingBoxComputer &signatureComputer,
                                                       float threshold,
                                                       CandidateScratch &work) {
        work.queue.clear();
        work.ordered.clear();
        work.queue.push(tree.getRoot());

        while (!work.queue.empty()) {
            const NodeId nodeId = work.queue.pop();
            if (nodeId == InvalidNode ||
                !tree.isNode(nodeId) ||
                !tree.isAlive(nodeId)) {
                continue;
            }

            if (nodeId != tree.getRoot() &&
                static_cast<std::size_t>(nodeId) < attribute.size() &&
                attribute[static_cast<std::size_t>(nodeId)] <= threshold) {
                CandidateOrderEntry entry;
                entry.nodeId = nodeId;
                entry.altitude = tree.getAltitude(nodeId);
                entry.supportMinPixel = signatureComputer.supportMinPixel(nodeId);
                entry.supportSize = signatureComputer.supportSize(nodeId);
                work.ordered.push_back(entry);
            }

            for (NodeId childId : tree.getChildren(nodeId)) {
                work.queue.push(childId);
            }
        }

        std::sort(work.ordered.begin(),
                  work.ordered.end(),
                  [](const CandidateOrderEntry &lhs,
                     const CandidateOrderEntry &rhs) {
                      return std::tie(lhs.altitude, lhs.supportMinPixel, lhs.supportSize, lhs.nodeId) <
                             std::tie(rhs.altitude, rhs.supportMinPixel, rhs.supportSize, rhs.nodeId);
                  });

        work.output.clear();
        work.output.reserve(work.ordered.size());
        for (const auto &entry : work.ordered) {
            work.output.push_back(entry.nodeId);
        }

        return work.output;
    }

    /**
     * @brief Checks whether `nodeId` is still mergeable after previous contractions.
     */
    static bool isCurrentMergeCandidate(const DynamicComponentTree &tree, NodeId nodeId) {
        return nodeId != InvalidNode &&
               tree.isNode(nodeId) &&
               tree.isAlive(nodeId) &&
               nodeId != tree.getRoot() &&
               tree.getNumProperParts(nodeId) > 0;
    }

    /**
     * @brief Registers the current attribute computers and buffers in the adjuster.
     */
    void registerAttributeComputers() {
        adjust_->setAttributeComputer(*minAttributeComputer_,
                                      *maxAttributeComputer_,
                                      std::span<float>(minAttribute_),
                                      std::span<float>(maxAttribute_));
    }

    /**
     * @brief Ensures attribute buffers are large enough for the next local merge.
     */
    void ensureAttributeRuntimeCapacityForSingleMerge() {
        bool resized = false;
        const auto maxCapacity = bufferCapacityForSingleMerge(*maxtree_);
        const auto minCapacity = bufferCapacityForSingleMerge(*mintree_);
        if (maxAttribute_.size() < maxCapacity) {
            maxAttribute_.resize(maxCapacity, 0.0f);
            resized = true;
        }
        if (minAttribute_.size() < minCapacity) {
            minAttribute_.resize(minCapacity, 0.0f);
            resized = true;
        }
        if (resized) {
            registerAttributeComputers();
        }
    }

    /**
     * @brief Reserves candidate scratch space according to the current tree sizes.
     */
    void reserveCandidateScratch() {
        const auto slots = static_cast<std::size_t>(
            std::max(maxtree_ ? maxtree_->getNumInternalNodeSlots() : 0,
                     mintree_ ? mintree_->getNumInternalNodeSlots() : 0));
        candidateScratch_.queue.reserve(slots);
        candidateScratch_.ordered.reserve(slots);
        candidateScratch_.output.reserve(slots);
    }

    /**
     * @brief Recomputes both attribute buffers from the current min/max trees.
     *
     * This is the bootstrap path used after tree construction. Subsequent local
     * contractions are handled incrementally by the registered computers.
     */
    void refreshAttributeBuffers() {
        ensureAttributeRuntimeCapacityForSingleMerge();
        clearActiveAttributeSlots(*maxtree_, maxAttribute_);
        clearActiveAttributeSlots(*mintree_, minAttribute_);
        maxAttributeComputer_->compute(std::span<float>(maxAttribute_));
        minAttributeComputer_->compute(std::span<float>(minAttribute_));
        registerAttributeComputers();
    }

    /**
     * @brief Builds the paired dynamic trees and initializes their attribute state.
     */
    void rebuildFromImage(const ImageUInt8Ptr &image) {
        if (image == nullptr) {
            throw std::runtime_error("ComponentTreeCasfNonIncreasing requires a valid image.");
        }

        maxtree_ = std::make_unique<DynamicComponentTree>(image->clone(), true, adjacency_);
        mintree_ = std::make_unique<DynamicComponentTree>(image->clone(), false, adjacency_);
        maxAttributeComputer_ = std::make_unique<DynamicBoundingBoxComputer>(
            maxtree_.get(),
            attribute_);
        minAttributeComputer_ = std::make_unique<DynamicBoundingBoxComputer>(
            mintree_.get(),
            attribute_);
        adjust_ = std::make_unique<DualMinMaxTreeIncrementalFilter<PixelType>>(
            mintree_.get(),
            maxtree_.get(),
            *adjacency_);
        reserveCandidateScratch();
        refreshAttributeBuffers();
    }

    /**
     * @brief Applies one snapshotted threshold phase to a primal/dual direction.
     *
     * `primalTree` supplies the candidates, while `dualTree` is updated by the
     * adjuster so the max-tree/min-tree pair remains synchronized. Candidate
     * attributes are not re-evaluated after every contraction in this phase;
     * only candidates that cease to be live or mergeable are skipped.
     */
    void applyMergePhase(DynamicComponentTree *primalTree,
                         DynamicComponentTree *dualTree,
                         std::vector<float> &primalAttribute,
                         const DynamicBoundingBoxComputer &primalSignatureComputer,
                         float threshold) {
        const auto &candidates = collectMergeCandidates(
            *primalTree,
            primalAttribute,
            primalSignatureComputer,
            threshold,
            candidateScratch_);

        for (NodeId nodeId : candidates) {
            if (!isCurrentMergeCandidate(*primalTree, nodeId)) {
                continue;
            }

            ensureAttributeRuntimeCapacityForSingleMerge();
            adjust_->mergeNodeIntoParentAndUpdateDual(primalTree, dualTree, nodeId);
        }
    }

    /**
     * @brief Applies one CASF threshold to the max-tree direction and then the min-tree direction.
     */
    void applyThreshold(float threshold) {
        applyMergePhase(
            maxtree_.get(),
            mintree_.get(),
            maxAttribute_,
            *maxAttributeComputer_,
            threshold);
        applyMergePhase(
            mintree_.get(),
            maxtree_.get(),
            minAttribute_,
            *minAttributeComputer_,
            threshold);
    }

public:
    /**
     * @brief Builds the filter from an image and an adjacency radius.
     * @param image Input image used to initialize the paired component trees.
     * The image is cloned and is not retained.
     * @param radiusAdj Non-negative radius used to build the internally owned
     * adjacency relation.
     * @param attribute Non-order-preserving flat attribute used for local
     * node selection and contraction. Accepted values are
     * `AREA_PER_BOX_PERIMETER`, `RECTANGULARITY`, `BOX_ASPECT_RATIO`,
     * `BOX_ELONGATION`, and `AREA_PER_BOX_DIAGONAL_SQUARED`.
     * @throws std::runtime_error If `image` is null or `attribute` is unsupported.
     * @throws std::invalid_argument If the image dimensions or radius cannot
     * define a valid adjacency relation.
     */
    ComponentTreeCasfNonIncreasing(ImageUInt8Ptr image,
                                   double radiusAdj,
                                   Attribute attribute = AREA_PER_BOX_PERIMETER)
        : adjacency_(image
              ? std::make_shared<AdjacencyRelation>(
                    image->getNumRows(),
                    image->getNumCols(),
                    radiusAdj)
              : nullptr),
          attribute_(attribute) {
        requireSupportedAttribute(attribute_);
        rebuildFromImage(image);
    }

    /**
     * @brief Builds the filter from an image and an existing adjacency relation.
     * @param image Input image used to initialize the paired component trees.
     * The image is cloned and is not retained.
     * @param adjacency Non-null adjacency retained through shared ownership by
     * this class and both dynamic trees. It must describe the image domain.
     * @param attribute Non-order-preserving flat attribute used for local
     * node selection and contraction. Accepted values are
     * `AREA_PER_BOX_PERIMETER`, `RECTANGULARITY`, `BOX_ASPECT_RATIO`,
     * `BOX_ELONGATION`, and `AREA_PER_BOX_DIAGONAL_SQUARED`.
     * @throws std::runtime_error If `image` or `adjacency` is null, or if
     * `attribute` is unsupported.
     */
    ComponentTreeCasfNonIncreasing(ImageUInt8Ptr image,
                                   AdjacencyRelationPtr adjacency,
                                   Attribute attribute = AREA_PER_BOX_PERIMETER)
        : adjacency_(std::move(adjacency)),
          attribute_(attribute) {
        requireSupportedAttribute(attribute_);
        if (adjacency_ == nullptr) {
            throw std::runtime_error(
                "ComponentTreeCasfNonIncreasing requires a valid adjacency relation.");
        }
        rebuildFromImage(image);
    }

    /**
     * @brief Applies all thresholds and returns the reconstructed filtered image.
     *
     * Thresholds are applied sequentially to the same mutable tree pair, so each
     * value continues from the state produced by the previous one. For each
     * threshold, max-tree candidates are snapshotted and merged first; min-tree
     * candidates are snapshotted from the resulting state and merged second.
     * A live non-root candidate is selected when its value at snapshot time is
     * `<=` the threshold. `NaN` values are not selected at that point. Surviving
     * candidates are not tested against the threshold again after earlier merges
     * in the same phase, even if their recomputed value would differ. An empty
     * sequence leaves the current state unchanged.
     *
     * @param thresholds Ordered floating-point thresholds.
     * @return A newly reconstructed image representing the current min-tree.
     * @throws std::invalid_argument If the paired update engine rejects an
     * invalid tree, node, pixel region, or other violated precondition.
     * @throws std::runtime_error If a required parent/owner cannot be found or
     * a local paired-tree update cannot be completed.
     * @warning Filtering is not transactional. If an exception is raised after
     * an earlier contraction, that successful mutation is not rolled back.
     */
    ImageUInt8Ptr filter(const std::vector<float> &thresholds) {
        for (float threshold : thresholds) {
            applyThreshold(threshold);
        }
        return reconstructImage();
    }

    /**
     * @brief Applies all thresholds from a span and returns the reconstructed image.
     * @param thresholds Ordered floating-point thresholds applied cumulatively.
     * Candidate snapshots, comparison rules, exception behavior, and mutation
     * order are identical to the `std::vector<float>` overload.
     * @return A newly reconstructed image representing the current min-tree.
     * @throws std::invalid_argument If the paired update engine rejects an
     * invalid tree, node, pixel region, or other violated precondition.
     * @throws std::runtime_error If a local paired-tree update cannot be completed.
     * @warning Filtering is not transactional; completed contractions are not
     * rolled back if a later operation throws.
     */
    ImageUInt8Ptr filter(std::span<const float> thresholds) {
        for (float threshold : thresholds) {
            applyThreshold(threshold);
        }
        return reconstructImage();
    }

    /**
     * @brief Reconstructs the current filtered image from the min-tree.
     * @return A newly allocated image for the current cumulative state.
     */
    ImageUInt8Ptr reconstructImage() const {
        return mintree_->reconstructionImage();
    }

    /**
     * @brief Returns the mutable max-tree after the current filtering state.
     * @return Reference valid until this CASF object is destroyed.
     * @warning External structural edits bypass the dual update engine and
     * attribute computers and can desynchronize subsequent filtering.
     */
    DynamicComponentTree &getMaxTree() {
        return *maxtree_;
    }

    /**
     * @brief Returns the read-only max-tree after the current filtering state.
     * @return Reference valid until this CASF object is destroyed.
     */
    const DynamicComponentTree &getMaxTree() const {
        return *maxtree_;
    }

    /**
     * @brief Returns the mutable min-tree after the current filtering state.
     * @return Reference valid until this CASF object is destroyed.
     * @warning External structural edits bypass the dual update engine and
     * attribute computers and can desynchronize subsequent filtering.
     */
    DynamicComponentTree &getMinTree() {
        return *mintree_;
    }

    /**
     * @brief Returns the read-only min-tree after the current filtering state.
     * @return Reference valid until this CASF object is destroyed.
     */
    const DynamicComponentTree &getMinTree() const {
        return *mintree_;
    }
};
