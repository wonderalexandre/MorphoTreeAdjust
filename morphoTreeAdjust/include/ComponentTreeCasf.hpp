#pragma once

#include <cctype>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "DynamicComponentTree.hpp"
#include "DualMinMaxTreeIncrementalFilter.hpp"

/**
 * @brief Conventional CASF over a synchronized dynamic max-tree/min-tree pair.
 *
 * This path is intended for the order-preserving attributes `AREA`,
 * `BOX_WIDTH`, `BOX_HEIGHT`, and `DIAGONAL_LENGTH`. For the supported flat
 * non-order-preserving bounding-box ratios, use
 * `ComponentTreeCasfNonIncreasing` instead.
 *
 * A threshold performs a breadth-first frontier selection: nodes whose current
 * attribute is greater than the threshold are expanded, and the first node at
 * or below the threshold on each path becomes a candidate while its descendants
 * are skipped. The root is never removed. If the root itself qualifies, it is
 * the only candidate and that direction performs no pruning. For every
 * threshold, the max-tree direction is processed before the min-tree direction.
 * Thresholds and repeated calls to `filter()` are cumulative: the object keeps
 * the mutable pair produced by all preceding operations.
 *
 * `Mode::Updating` edits the existing tree objects in place. `Mode::Naive`
 * rebuilds and replaces both tree objects after every threshold, while
 * `Mode::Hybrid` does so for the first threshold of each non-empty `filter()`
 * call and then switches to the incremental path. Consequently, references returned by `getMaxTree()` or
 * `getMinTree()` are invalidated whenever a naive rebuild occurs (and when the
 * CASF object is destroyed).
 *
 * @tparam PixelType Unsigned integral altitude type used by the incremental
 * adjuster. Its level domain must start at zero.
 *
 * @warning Mutating a tree obtained from a non-const getter bypasses CASF's
 * paired-tree and attribute maintenance and can desynchronize subsequent
 * filtering operations.
 */
template<typename PixelType = AltitudeType>
class ComponentTreeCasf {
public:
    /** @brief Strategy used to execute a threshold sequence. */
    enum class Mode {
        Updating, ///< Incrementally update the existing paired trees.
        Naive, ///< Reconstruct and replace both trees after every threshold.
        Hybrid, ///< Rebuild for the first threshold of each non-empty call, then update incrementally.
    };

private:
    AdjacencyRelationPtr adjacency_;
    Attribute attribute_ = AREA;
    std::unique_ptr<DynamicComponentTree> maxtree_;
    std::unique_ptr<DynamicComponentTree> mintree_;
    std::unique_ptr<DynamicAttributeComputer> maxAttributeComputer_;
    std::unique_ptr<DynamicAttributeComputer> minAttributeComputer_;
    std::vector<float> maxAttribute_;
    std::vector<float> minAttribute_;
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<PixelType>> adjust_;

    static std::string normalizeToken(std::string_view token) {
        std::string normalized(token.begin(), token.end());
        for (char &c : normalized) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        return normalized;
    }

    static std::vector<NodeId> getNodesToPrune(const DynamicComponentTree &tree, const std::vector<float> &attribute, int threshold) {
        std::vector<NodeId> out;
        FastQueue<NodeId> queue;
        queue.push(tree.getRoot());
        while (!queue.empty()) {
            const NodeId nodeId = queue.pop();
            if (attribute[static_cast<std::size_t>(nodeId)] > threshold) {
                for (NodeId childId : tree.getChildren(nodeId)) {
                    queue.push(childId);
                }
            } else {
                out.push_back(nodeId);
            }
        }
        return out;
    }

    static std::vector<float> computeAttribute(DynamicComponentTree &tree, Attribute attribute) {
        if (attribute == AREA) {
            DynamicAreaComputer computer(&tree);
            return computer.compute();
        }

        DynamicBoundingBoxComputer computer(&tree, attribute);
        return computer.compute();
    }

    static std::unique_ptr<DynamicAttributeComputer> makeAttributeComputer(DynamicComponentTree *tree, Attribute attribute) {
        if (attribute == AREA) {
            return std::make_unique<DynamicAreaComputer>(tree);
        }

        return std::make_unique<DynamicBoundingBoxComputer>(tree, attribute);
    }

    void rebuildFromImage(const ImageUInt8Ptr &image) {
        if (image == nullptr) {
            throw std::runtime_error("ComponentTreeCasf requires a valid image.");
        }

        maxtree_ = std::make_unique<DynamicComponentTree>(image->clone(), true, adjacency_);
        mintree_ = std::make_unique<DynamicComponentTree>(image->clone(), false, adjacency_);
        maxAttributeComputer_ = makeAttributeComputer(maxtree_.get(), attribute_);
        minAttributeComputer_ = makeAttributeComputer(mintree_.get(), attribute_);
        adjust_ = std::make_unique<DualMinMaxTreeIncrementalFilter<PixelType>>(mintree_.get(), maxtree_.get(), *adjacency_);
        refreshAttributeBuffers();
    }

    void refreshAttributeBuffers() {
        maxAttribute_ = computeAttribute(*maxtree_, attribute_);
        minAttribute_ = computeAttribute(*mintree_, attribute_);
        adjust_->setAttributeComputer(*minAttributeComputer_, *maxAttributeComputer_, std::span<float>(minAttribute_), std::span<float>(maxAttribute_));
    }

    void applyUpdatingThreshold(int threshold) {
        auto maxNodes = getNodesToPrune(*maxtree_, maxAttribute_, threshold);
        adjust_->pruneMaxTreeAndUpdateMinTree(maxNodes);

        auto minNodes = getNodesToPrune(*mintree_, minAttribute_, threshold);
        adjust_->pruneMinTreeAndUpdateMaxTree(minNodes);
    }

    void applyNaiveThreshold(int threshold) {
        auto current = mintree_->reconstructionImage();

        DynamicComponentTree maxTree(current, true, adjacency_);
        const auto maxAttribute = computeAttribute(maxTree, attribute_);
        for (NodeId nodeId : getNodesToPrune(maxTree, maxAttribute, threshold)) {
            if (nodeId != maxTree.getRoot()) {
                maxTree.pruneNode(nodeId);
            }
        }

        current = maxTree.reconstructionImage();
        DynamicComponentTree minTree(current, false, adjacency_);
        const auto minAttribute = computeAttribute(minTree, attribute_);
        for (NodeId nodeId : getNodesToPrune(minTree, minAttribute, threshold)) {
            if (nodeId != minTree.getRoot()) {
                minTree.pruneNode(nodeId);
            }
        }

        rebuildFromImage(minTree.reconstructionImage());
    }

public:
    /**
     * @brief Builds a conventional CASF and creates its adjacency from a radius.
     * @param image Input image. Its pixels are cloned into the two dynamic trees;
     * the input image itself is not retained.
     * @param radiusAdj Non-negative Euclidean radius used by the internally
     * owned adjacency relation.
     * @param attribute Order-preserving attribute. Supported values are `AREA`,
     * `BOX_WIDTH`, `BOX_HEIGHT`, and `DIAGONAL_LENGTH`.
     * @pre `attribute` is one of the supported values; this constructor does
     * not reject non-order-preserving values at runtime.
     * @throws std::runtime_error If `image` is null.
     * @throws std::invalid_argument If the image dimensions or radius cannot
     * define a valid adjacency relation.
     */
    ComponentTreeCasf(ImageUInt8Ptr image, double radiusAdj, Attribute attribute = AREA)
        : adjacency_(image ? std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radiusAdj) : nullptr), attribute_(attribute) {
        rebuildFromImage(image);
    }

    /**
     * @brief Builds a conventional CASF with an existing shared adjacency.
     * @param image Input image. Its pixels are cloned into the two dynamic trees;
     * the input image itself is not retained.
     * @param adjacency Non-null adjacency retained through shared ownership. It
     * must describe the same image domain as `image`.
     * @param attribute Order-preserving attribute. Supported values are `AREA`,
     * `BOX_WIDTH`, `BOX_HEIGHT`, and `DIAGONAL_LENGTH`.
     * @pre `attribute` is one of the supported values and `adjacency` describes
     * the same domain as `image`; these conditions are not fully runtime-checked.
     * @throws std::runtime_error If `image` or `adjacency` is null.
     */
    ComponentTreeCasf(ImageUInt8Ptr image, AdjacencyRelationPtr adjacency, Attribute attribute = AREA)
        : adjacency_(std::move(adjacency)), attribute_(attribute) {
        if (adjacency_ == nullptr) {
            throw std::runtime_error("ComponentTreeCasf requires a valid adjacency relation.");
        }
        rebuildFromImage(image);
    }

    /**
     * @brief Parses a case-insensitive execution-mode token.
     * @param mode `"updating"`, `"naive"`, or `"hybrid"`. The legacy
     * misspelling `"updaing"` is also accepted as `Mode::Updating`.
     * @return The corresponding execution mode.
     * @throws std::runtime_error If the normalized token is not recognized.
     */
    static Mode parseMode(std::string_view mode) {
        const std::string normalized = normalizeToken(mode);
        if (normalized == "updating" || normalized == "updaing") {
            return Mode::Updating;
        }
        if (normalized == "naive") {
            return Mode::Naive;
        }
        if (normalized == "hybrid") {
            return Mode::Hybrid;
        }

        throw std::runtime_error("Unknown ComponentTreeCasf mode. Expected one of: updating, naive, hybrid.");
    }

    /**
     * @brief Applies a sequence of conventional CASF thresholds.
     *
     * Each threshold continues from the result of the preceding threshold and
     * from any earlier `filter()` call. A node is eligible when its current
     * attribute is less than or equal to the threshold, and descendants of an
     * eligible node are not examined independently. The global root is never
     * removed; if it qualifies, that direction stops without pruning a
     * descendant. Within each threshold the max-tree direction runs before the
     * min-tree direction. An empty sequence leaves the current state unchanged.
     *
     * In `Mode::Naive`, every threshold replaces the tree objects. In
     * `Mode::Hybrid`, the first threshold of every non-empty call replaces them.
     * References previously returned by the tree getters are invalid after
     * either event.
     *
     * @param thresholds Ordered integer thresholds.
     * @param mode Execution strategy; defaults to incremental updating.
     * @return A newly reconstructed image representing the current min-tree.
     * @warning Filtering is not transactional. An exception does not roll back
     * thresholds, subtree operations, or tree replacements already completed;
     * references can therefore be invalidated even when the call fails.
     */
    ImageUInt8Ptr filter(const std::vector<int> &thresholds, Mode mode = Mode::Updating) {
        switch (mode) {
            case Mode::Updating:
                for (int threshold : thresholds) {
                    applyUpdatingThreshold(threshold);
                }
                break;
            case Mode::Naive:
                for (int threshold : thresholds) {
                    applyNaiveThreshold(threshold);
                }
                break;
            case Mode::Hybrid:
                if (!thresholds.empty()) {
                    applyNaiveThreshold(thresholds.front());
                    for (std::size_t i = 1; i < thresholds.size(); ++i) {
                        applyUpdatingThreshold(thresholds[i]);
                    }
                }
                break;
        }

        return mintree_->reconstructionImage();
    }

    /**
     * @brief Applies thresholds using a case-insensitive mode token.
     * @param thresholds Ordered integer thresholds applied cumulatively.
     * @param mode Token accepted by `parseMode()`.
     * @return A newly reconstructed image representing the current min-tree.
     * @throws std::runtime_error If `mode` is unknown.
     */
    ImageUInt8Ptr filter(const std::vector<int> &thresholds, std::string_view mode) {
        return filter(thresholds, parseMode(mode));
    }

    /**
     * @brief Returns mutable access to the current max-tree.
     * @return Reference valid until the next naive rebuild or object destruction.
     * @warning External structural edits bypass paired-tree and attribute
     * synchronization and can make later filtering invalid.
     */
    DynamicComponentTree &getMaxTree() {
        return *maxtree_;
    }

    /**
     * @brief Returns read-only access to the current max-tree.
     * @return Reference valid until the next naive rebuild or object destruction.
     */
    const DynamicComponentTree &getMaxTree() const {
        return *maxtree_;
    }

    /**
     * @brief Returns mutable access to the current min-tree.
     * @return Reference valid until the next naive rebuild or object destruction.
     * @warning External structural edits bypass paired-tree and attribute
     * synchronization and can make later filtering invalid.
     */
    DynamicComponentTree &getMinTree() {
        return *mintree_;
    }

    /**
     * @brief Returns read-only access to the current min-tree.
     * @return Reference valid until the next naive rebuild or object destruction.
     */
    const DynamicComponentTree &getMinTree() const {
        return *mintree_;
    }
};
