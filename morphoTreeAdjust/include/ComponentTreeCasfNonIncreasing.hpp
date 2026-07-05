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

template<typename PixelType = AltitudeType>
class ComponentTreeCasfNonIncreasing {
private:
    struct CandidateOrderEntry {
        NodeId nodeId = InvalidNode;
        int altitude = 0;
        PixelId supportMinPixel = std::numeric_limits<PixelId>::max();
        int supportSize = 0;
    };

    struct CandidateScratch {
        FastQueue<NodeId> queue;
        std::vector<NodeId> traversalOrder;
        std::vector<NodeId> scanOrder;
        std::vector<CandidateOrderEntry> ordered;
        std::vector<NodeId> output;
        std::vector<int> supportSizeByNode;
        std::vector<PixelId> supportMinPixelByNode;
    };

    AdjacencyRelationPtr adjacency_;
    Attribute attribute_ = AREA_PER_BOX_PERIMETER;
    std::unique_ptr<DynamicComponentTree> maxtree_;
    std::unique_ptr<DynamicComponentTree> mintree_;
    std::unique_ptr<DynamicBoundingBoxComputer> maxAttributeComputer_;
    std::unique_ptr<DynamicBoundingBoxComputer> minAttributeComputer_;
    std::vector<float> maxAttribute_;
    std::vector<float> minAttribute_;
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<PixelType>> adjust_;
    CandidateScratch candidateScratch_;

    static void requireSupportedAttribute(Attribute attribute) {
        if (attribute != AREA_PER_BOX_PERIMETER) {
            throw std::runtime_error(
                "ComponentTreeCasfNonIncreasing currently supports only AREA_PER_BOX_PERIMETER.");
        }
    }

    static std::size_t bufferCapacityForSingleMerge(const DynamicComponentTree &tree) {
        const auto slots = static_cast<std::size_t>(
            std::max(0, tree.getNumInternalNodeSlots()));
        const auto pixels = static_cast<std::size_t>(
            std::max(1, tree.getNumTotalProperParts()));
        return slots + pixels;
    }

    static std::vector<NodeId> &collectMergeCandidates(const DynamicComponentTree &tree,
                                                       const std::vector<float> &attribute,
                                                       float threshold,
                                                       CandidateScratch &work) {
        work.queue.clear();
        work.traversalOrder.clear();
        work.scanOrder.clear();
        work.queue.push(tree.getRoot());

        while (!work.queue.empty()) {
            const NodeId nodeId = work.queue.pop();
            if (nodeId == InvalidNode ||
                !tree.isNode(nodeId) ||
                !tree.isAlive(nodeId)) {
                continue;
            }

            work.traversalOrder.push_back(nodeId);
            if (nodeId != tree.getRoot() &&
                static_cast<std::size_t>(nodeId) < attribute.size() &&
                attribute[static_cast<std::size_t>(nodeId)] <= threshold) {
                work.scanOrder.push_back(nodeId);
            }

            for (NodeId childId : tree.getChildren(nodeId)) {
                work.queue.push(childId);
            }
        }

        const auto nodeSlots = static_cast<std::size_t>(
            tree.getNumInternalNodeSlots());
        if (work.supportSizeByNode.size() < nodeSlots) {
            work.supportSizeByNode.resize(nodeSlots, 0);
        }
        if (work.supportMinPixelByNode.size() < nodeSlots) {
            work.supportMinPixelByNode.resize(
                nodeSlots,
                std::numeric_limits<PixelId>::max());
        }

        for (auto it = work.traversalOrder.rbegin(); it != work.traversalOrder.rend(); ++it) {
            const NodeId nodeId = *it;
            const auto index = static_cast<std::size_t>(nodeId);
            int supportSize = tree.getNumProperParts(nodeId);
            PixelId supportMinPixel = std::numeric_limits<PixelId>::max();

            for (PixelId pixelId : tree.getProperParts(nodeId)) {
                supportMinPixel = std::min(supportMinPixel, pixelId);
            }

            for (NodeId childId : tree.getChildren(nodeId)) {
                if (childId == InvalidNode ||
                    childId < 0 ||
                    static_cast<std::size_t>(childId) >= work.supportSizeByNode.size()) {
                    continue;
                }

                const auto childIndex = static_cast<std::size_t>(childId);
                supportSize += work.supportSizeByNode[childIndex];
                supportMinPixel =
                    std::min(supportMinPixel, work.supportMinPixelByNode[childIndex]);
            }

            work.supportSizeByNode[index] = supportSize;
            work.supportMinPixelByNode[index] = supportMinPixel;
        }

        work.ordered.resize(work.scanOrder.size());
        for (std::size_t index = 0; index < work.scanOrder.size(); ++index) {
            const NodeId nodeId = work.scanOrder[index];
            auto &entry = work.ordered[index];
            entry.nodeId = nodeId;
            entry.altitude = tree.getAltitude(nodeId);
            const auto nodeIndex = static_cast<std::size_t>(nodeId);
            entry.supportMinPixel = work.supportMinPixelByNode[nodeIndex];
            entry.supportSize = work.supportSizeByNode[nodeIndex];
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

    static bool isCurrentMergeCandidate(const DynamicComponentTree &tree, NodeId nodeId) {
        return nodeId != InvalidNode &&
               tree.isNode(nodeId) &&
               tree.isAlive(nodeId) &&
               nodeId != tree.getRoot() &&
               tree.getNumProperParts(nodeId) > 0;
    }

    void registerAttributeComputers() {
        adjust_->setAttributeComputer(*minAttributeComputer_,
                                      *maxAttributeComputer_,
                                      std::span<float>(minAttribute_),
                                      std::span<float>(maxAttribute_));
    }

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

    void refreshAttributeBuffers() {
        ensureAttributeRuntimeCapacityForSingleMerge();
        std::fill(maxAttribute_.begin(), maxAttribute_.end(), 0.0f);
        std::fill(minAttribute_.begin(), minAttribute_.end(), 0.0f);
        maxAttributeComputer_->compute(std::span<float>(maxAttribute_));
        minAttributeComputer_->compute(std::span<float>(minAttribute_));
        registerAttributeComputers();
    }

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
        refreshAttributeBuffers();
    }

    void applyMergePhase(DynamicComponentTree *primalTree,
                         DynamicComponentTree *dualTree,
                         std::vector<float> &primalAttribute,
                         float threshold) {
        const auto candidates = collectMergeCandidates(
            *primalTree,
            primalAttribute,
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

    void applyThreshold(float threshold) {
        applyMergePhase(maxtree_.get(), mintree_.get(), maxAttribute_, threshold);
        applyMergePhase(mintree_.get(), maxtree_.get(), minAttribute_, threshold);
    }

public:
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

    ImageUInt8Ptr filter(const std::vector<float> &thresholds) {
        for (float threshold : thresholds) {
            applyThreshold(threshold);
        }
        return reconstructImage();
    }

    ImageUInt8Ptr filter(std::span<const float> thresholds) {
        for (float threshold : thresholds) {
            applyThreshold(threshold);
        }
        return reconstructImage();
    }

    ImageUInt8Ptr reconstructImage() const {
        return mintree_->reconstructionImage();
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
