#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Common.hpp"
#include "DynamicComponentTree.hpp"
#include "details/NonIncreasingCasfKernel.hpp"

class ComponentTreeCasfNonIncreasing {
public:
    enum class Attribute {
        BoxAreaPerimeterRatio,
    };

    explicit ComponentTreeCasfNonIncreasing(
        ImageUInt8Ptr image,
        double radiusAdj,
        Attribute attribute = Attribute::BoxAreaPerimeterRatio)
        : casf_(std::make_unique<morphotree::details::non_increasing::NonIncreasingCasfKernel>(
              std::move(image),
              radiusAdj,
              morphotree::details::non_increasing::NonIncreasingDualMaintenance::MergeNodeCleanUnchecked,
              toKernelAttribute(attribute))) {}

    ~ComponentTreeCasfNonIncreasing() = default;

    ComponentTreeCasfNonIncreasing(const ComponentTreeCasfNonIncreasing &) = delete;
    ComponentTreeCasfNonIncreasing &operator=(const ComponentTreeCasfNonIncreasing &) = delete;
    ComponentTreeCasfNonIncreasing(ComponentTreeCasfNonIncreasing &&) noexcept = default;
    ComponentTreeCasfNonIncreasing &operator=(ComponentTreeCasfNonIncreasing &&) noexcept = default;

    ImageUInt8Ptr filter(const std::vector<int> &thresholds) {
        return filter(toFloatThresholds(thresholds));
    }

    ImageUInt8Ptr filter(const std::vector<float> &thresholds) {
        for (float threshold : thresholds) {
            casf_->applyThreshold(threshold);
            requireAcceptedCleanThreshold(threshold);
        }
        return casf_->reconstructImage();
    }

    ImageUInt8Ptr reconstructImage() const {
        return casf_->reconstructImage();
    }

    DynamicComponentTree &getMaxTree() {
        return casf_->maxTree();
    }

    const DynamicComponentTree &getMaxTree() const {
        return casf_->maxTree();
    }

    DynamicComponentTree &getMinTree() {
        return casf_->minTree();
    }

    const DynamicComponentTree &getMinTree() const {
        return casf_->minTree();
    }

private:
    std::unique_ptr<morphotree::details::non_increasing::NonIncreasingCasfKernel> casf_;

    static morphotree::details::non_increasing::NonIncreasingAttributeKind toKernelAttribute(
        Attribute attribute) {
        switch (attribute) {
            case Attribute::BoxAreaPerimeterRatio:
                return morphotree::details::non_increasing::NonIncreasingAttributeKind::BoxAreaPerimeterRatio;
        }
        return morphotree::details::non_increasing::NonIncreasingAttributeKind::BoxAreaPerimeterRatio;
    }

    static std::vector<float> toFloatThresholds(const std::vector<int> &thresholds) {
        std::vector<float> out;
        out.reserve(thresholds.size());
        for (int threshold : thresholds) {
            out.push_back(static_cast<float>(threshold));
        }
        return out;
    }

    void requireAcceptedCleanThreshold(float threshold) const {
        if (casf_->lastLocalUpdateRejectionStopReached()) {
            throw std::runtime_error(
                "ComponentTreeCasfNonIncreasing clean local update rejected threshold " +
                std::to_string(threshold));
        }
        if (casf_->lastPhaseFallbackStopReached() || casf_->lastPhaseDualRebuilds() != 0) {
            throw std::runtime_error(
                "ComponentTreeCasfNonIncreasing clean local update required a phase rebuild at threshold " +
                std::to_string(threshold));
        }
    }
};
