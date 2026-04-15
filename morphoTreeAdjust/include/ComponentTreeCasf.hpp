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
#include "DynamicComponentTreeAdjustment.hpp"

template<typename PixelType = AltitudeType>
class ComponentTreeCasf {
public:
    enum class Mode {
        Updating,
        Naive,
        Hybrid,
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
    std::unique_ptr<DynamicComponentTreeAdjustment<PixelType>> adjust_;

    static std::string normalizeToken(std::string_view token) {
        std::string normalized(token.begin(), token.end());
        for (char &c : normalized) {
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        }
        return normalized;
    }

    static std::vector<NodeId> getNodesToPrune(const DynamicComponentTree &tree,
                                               const std::vector<float> &attribute,
                                               int threshold) {
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
        adjust_ = std::make_unique<DynamicComponentTreeAdjustment<PixelType>>(mintree_.get(), maxtree_.get(), *adjacency_);
        refreshAttributeBuffers();
    }

    void refreshAttributeBuffers() {
        maxAttribute_ = computeAttribute(*maxtree_, attribute_);
        minAttribute_ = computeAttribute(*mintree_, attribute_);
        adjust_->setAttributeComputer(*minAttributeComputer_,
                                      *maxAttributeComputer_,
                                      std::span<float>(minAttribute_),
                                      std::span<float>(maxAttribute_));
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
    ComponentTreeCasf(ImageUInt8Ptr image,
                      double radiusAdj,
                      Attribute attribute = AREA)
        : adjacency_(image ? std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radiusAdj) : nullptr),
          attribute_(attribute) {
        rebuildFromImage(image);
    }

    ComponentTreeCasf(ImageUInt8Ptr image,
                      AdjacencyRelationPtr adjacency,
                      Attribute attribute = AREA)
        : adjacency_(std::move(adjacency)),
          attribute_(attribute) {
        if (adjacency_ == nullptr) {
            throw std::runtime_error("ComponentTreeCasf requires a valid adjacency relation.");
        }
        rebuildFromImage(image);
    }

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

    ImageUInt8Ptr filter(const std::vector<int> &thresholds, std::string_view mode) {
        return filter(thresholds, parseMode(mode));
    }

    DynamicComponentTree &getMaxTree() {
        return *maxtree_;
    }

    const DynamicComponentTree &getMaxTree() const {
        return *maxtree_;
    }

    DynamicComponentTree &getMinTree() {
        return *mintree_;
    }

    const DynamicComponentTree &getMinTree() const {
        return *mintree_;
    }
};
