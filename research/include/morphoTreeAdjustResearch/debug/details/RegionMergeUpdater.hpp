#pragma once

#include <chrono>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "Algorithm2FlatRaiseAreaAttributeUpdater.hpp"
#include "Algorithm2FlatRegionMergeUpdater.hpp"
#include "RegionLevelChange.hpp"
#include "RegionLevelChangeRebuildUpdater.hpp"
#include "RegionMaterialization.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../AttributeComputer.hpp"
#include "../../Common.hpp"
#include "../../DualMinMaxTreeIncrementalFilter.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::region {

using ::debug::problem2::Algorithm2FlatRaiseAreaAttributeUpdater;
using ::debug::problem2::Algorithm2FlatRegionMergeUpdater;

struct RegionMergeAlgorithm2AreaOptInDiagnostics {
    bool attemptedAlgorithm2 = false;
    bool acceptedAlgorithm2 = false;
    bool rejectedByAttribute = false;
    bool rejectedByContract = false;
    bool fallbackUsed = false;
    bool areaAttributeUpdated = false;
    std::string rejectionReason;
    int createdNodes = 0;
    int absorbedNodes = 0;
    int releasedNodes = 0;
    int affectedLevels = 0;
    int xiBulletNodes = 0;
    int xiCircleNodes = 0;
    int areaDirtySeeds = 0;
    int areaDirtyClosure = 0;
    int areaRecomputedNodes = 0;
    long long algorithm2StructuralMicros = 0;
    long long areaAttributeMicros = 0;
    long long primalRebuildMicros = 0;
    long long fallbackRebuildMicros = 0;
};

class RegionMergeUpdater {
private:
    static void requireAreaAttributeForAlgorithm2Raise(Attribute attribute) {
        if (attribute != AREA) {
            throw std::runtime_error("RegionMergeUpdater Algorithm 2 raise backend is currently validated only for AREA attributes");
        }
    }

    static long long elapsedMicros(std::chrono::steady_clock::time_point start,
                                   std::chrono::steady_clock::time_point end) {
        return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr rebuildPairFromImageUpdate(
        DynamicComponentTree &minTree,
        DynamicComponentTree &maxTree,
        const AdjacencyRelationPtr &adjacency,
        const ImageUInt8Ptr &inputImage,
        const RegionLevelChange<PixelType> &change) {
        auto updatedImage = RegionLevelChangeRebuildUpdater::applyToImage(
            inputImage,
            change);
        maxTree.build(updatedImage->clone(), true, adjacency);
        minTree.build(updatedImage->clone(), false, adjacency);
        return updatedImage;
    }

public:
    template<typename PixelType = AltitudeType>
    static DynamicComponentTree &selectDualTree(DynamicComponentTree &minTree,
                                                DynamicComponentTree &maxTree,
                                                const RegionLevelChange<PixelType> &change) {
        if (change.direction == RegionUpdateDirection::Raise) {
            static_cast<void>(minTree);
            return maxTree;
        }
        if (change.direction == RegionUpdateDirection::Lower) {
            return minTree;
        }
        throw std::runtime_error("RegionMergeUpdater requires a non-trivial raise or lower update");
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr updateDualTreeOnly(DynamicComponentTree &minTree,
                                            DynamicComponentTree &maxTree,
                                            const AdjacencyRelationPtr &adjacency,
                                            const RegionLevelChange<PixelType> &change,
                                            Attribute attribute = AREA) {
        if (!adjacency) {
            throw std::runtime_error("RegionMergeUpdater requires a valid adjacency relation");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("RegionMergeUpdater requires a non-empty region C");
        }
        if (!change.assumeConnected &&
            !isConnectedRegion(*adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionMergeUpdater requires connected region C");
        }
        if (!change.assumeMaterializedAsProperParts) {
            throw std::runtime_error("RegionMergeUpdater currently requires C to be materialized as proper parts");
        }

        if (change.direction == RegionUpdateDirection::Raise) {
            requireAreaAttributeForAlgorithm2Raise(attribute);
            auto result = Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseMaterialized(
                maxTree,
                *adjacency,
                change);
            return result.reconstruction;
        }

        DynamicComponentTree &dualTree = selectDualTree(minTree, maxTree, change);
        const auto materialization = analyzeRegionMaterialization(dualTree, *adjacency, change);
        if (materialization.requiresLocalMaterialization()) {
            throw std::runtime_error("RegionMergeUpdater cannot yet process C that cuts a node's direct proper parts; run local materialization first");
        }

        DualMinMaxTreeIncrementalFilter<PixelType> adjust(&minTree, &maxTree, *adjacency);
        adjust.updateTreeWithExplicitRegion(&dualTree, std::span<const PixelId>(change.pixelsC), change.b);
        return dualTree.reconstructionImage();
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr updateDualTreeOnlyAfterMaterializingSinglePartialOwner(
        DynamicComponentTree &minTree,
        DynamicComponentTree &maxTree,
        const AdjacencyRelationPtr &adjacency,
        const RegionLevelChange<PixelType> &change,
        Attribute attribute = AREA) {
        if (!adjacency) {
            throw std::runtime_error("RegionMergeUpdater requires a valid adjacency relation");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("RegionMergeUpdater requires a non-empty region C");
        }
        if (!change.assumeConnected &&
            !isConnectedRegion(*adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionMergeUpdater requires connected region C");
        }

        if (change.direction == RegionUpdateDirection::Raise) {
            requireAreaAttributeForAlgorithm2Raise(attribute);
            const auto materialization = analyzeRegionMaterialization(maxTree, *adjacency, change);
            if (materialization.requiresLocalMaterialization() &&
                materialization.numPartialOwners() != 1) {
                throw std::runtime_error("RegionMergeUpdater can only materialize one partial owner in this debug path");
            }
            auto result = Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseAfterMaterializingPartialOwners(
                maxTree,
                *adjacency,
                change);
            return result.reconstruction;
        }

        DynamicComponentTree &dualTree = selectDualTree(minTree, maxTree, change);
        auto materialization = analyzeRegionMaterialization(dualTree, *adjacency, change);
        if (materialization.requiresLocalMaterialization()) {
            if (materialization.numPartialOwners() != 1) {
                throw std::runtime_error("RegionMergeUpdater can only materialize one partial owner in this debug path");
            }
            (void) materializeSinglePartialOwner(dualTree, materialization);
            materialization = analyzeRegionMaterialization(dualTree, *adjacency, change);
        }

        if (materialization.requiresLocalMaterialization()) {
            throw std::runtime_error("RegionMergeUpdater could not materialize C as complete direct proper parts");
        }

        RegionLevelChange<PixelType> materializedChange = change;
        materializedChange.assumeMaterializedAsProperParts = true;

        DualMinMaxTreeIncrementalFilter<PixelType> adjust(&minTree, &maxTree, *adjacency);
        adjust.updateTreeWithExplicitRegion(&dualTree, std::span<const PixelId>(materializedChange.pixelsC), materializedChange.b);
        return dualTree.reconstructionImage();
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr updateDualTreeOnlyAfterMaterializingPartialOwners(
        DynamicComponentTree &minTree,
        DynamicComponentTree &maxTree,
        const AdjacencyRelationPtr &adjacency,
        const RegionLevelChange<PixelType> &change,
        Attribute attribute = AREA) {
        if (!adjacency) {
            throw std::runtime_error("RegionMergeUpdater requires a valid adjacency relation");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("RegionMergeUpdater requires a non-empty region C");
        }
        if (!change.assumeConnected &&
            !isConnectedRegion(*adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionMergeUpdater requires connected region C");
        }

        if (change.direction == RegionUpdateDirection::Raise) {
            requireAreaAttributeForAlgorithm2Raise(attribute);
            auto result = Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseAfterMaterializingPartialOwners(
                maxTree,
                *adjacency,
                change);
            return result.reconstruction;
        }

        DynamicComponentTree &dualTree = selectDualTree(minTree, maxTree, change);
        auto materialization = analyzeRegionMaterialization(dualTree, *adjacency, change);
        if (materialization.requiresLocalMaterialization()) {
            (void) materializePartialOwners(dualTree, materialization);
            materialization = analyzeRegionMaterialization(dualTree, *adjacency, change);
        }

        if (materialization.requiresLocalMaterialization()) {
            throw std::runtime_error("RegionMergeUpdater could not materialize C as complete direct proper parts");
        }

        RegionLevelChange<PixelType> materializedChange = change;
        materializedChange.assumeMaterializedAsProperParts = true;

        DualMinMaxTreeIncrementalFilter<PixelType> adjust(&minTree, &maxTree, *adjacency);
        adjust.updateTreeWithExplicitRegion(&dualTree, std::span<const PixelId>(materializedChange.pixelsC), materializedChange.b);
        return dualTree.reconstructionImage();
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr updatePairWithLocalDualAndRebuiltPrimal(
        DynamicComponentTree &minTree,
        DynamicComponentTree &maxTree,
        const AdjacencyRelationPtr &adjacency,
        const RegionLevelChange<PixelType> &change,
        Attribute attribute = AREA) {
        DynamicComponentTree &dualTree = selectDualTree(minTree, maxTree, change);
        const bool dualIsMaxTree = &dualTree == &maxTree;

        auto updatedImage = updateDualTreeOnlyAfterMaterializingPartialOwners(
            minTree,
            maxTree,
            adjacency,
            change,
            attribute);

        if (dualIsMaxTree) {
            minTree.build(updatedImage->clone(), false, adjacency);
        } else {
            maxTree.build(updatedImage->clone(), true, adjacency);
        }
        return updatedImage;
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr updatePairWithAlgorithm2AreaOptIn(
        DynamicComponentTree &minTree,
        DynamicComponentTree &maxTree,
        const AdjacencyRelationPtr &adjacency,
        const RegionLevelChange<PixelType> &change,
        Attribute attribute,
        std::vector<float> *maxAreaBuffer = nullptr,
        RegionMergeAlgorithm2AreaOptInDiagnostics *diagnostics = nullptr) {
        if (!adjacency) {
            throw std::runtime_error("RegionMergeUpdater requires a valid adjacency relation");
        }
        if (diagnostics != nullptr) {
            *diagnostics = RegionMergeAlgorithm2AreaOptInDiagnostics{};
        }
        const auto inputImageBeforeAttempt = maxTree.reconstructionImage()->clone();

        auto useFallback = [&](const std::string &reason,
                               bool rejectedByAttribute,
                               bool rejectedByContract) {
            if (diagnostics != nullptr) {
                diagnostics->fallbackUsed = true;
                diagnostics->rejectedByAttribute = rejectedByAttribute;
                diagnostics->rejectedByContract = rejectedByContract;
                diagnostics->rejectionReason = reason;
            }
            const auto start = std::chrono::steady_clock::now();
            auto updatedImage = rebuildPairFromImageUpdate(minTree,
                                                           maxTree,
                                                           adjacency,
                                                           inputImageBeforeAttempt,
                                                           change);
            if (diagnostics != nullptr) {
                diagnostics->fallbackRebuildMicros =
                    elapsedMicros(start, std::chrono::steady_clock::now());
            }
            return updatedImage;
        };

        if (attribute != AREA) {
            return useFallback(
                "Algorithm 2 area opt-in backend only supports AREA",
                true,
                false);
        }

        if (change.direction != RegionUpdateDirection::Raise) {
            return useFallback(
                "Algorithm 2 area opt-in backend only supports raise updates",
                false,
                true);
        }

        if (diagnostics != nullptr) {
            diagnostics->attemptedAlgorithm2 = true;
        }

        try {
            const auto structuralStart = std::chrono::steady_clock::now();
            auto structuralResult =
                Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
                    maxTree,
                    *adjacency,
                    change);
            if (diagnostics != nullptr) {
                diagnostics->algorithm2StructuralMicros =
                    elapsedMicros(structuralStart, std::chrono::steady_clock::now());
            }

            if (maxAreaBuffer != nullptr) {
                const auto areaStart = std::chrono::steady_clock::now();
                const auto areaResult = Algorithm2FlatRaiseAreaAttributeUpdater::recomputeArea(
                    maxTree,
                    structuralResult,
                    *maxAreaBuffer);
                if (diagnostics != nullptr) {
                    diagnostics->areaAttributeMicros =
                        elapsedMicros(areaStart, std::chrono::steady_clock::now());
                    diagnostics->areaAttributeUpdated = true;
                    diagnostics->areaDirtySeeds = static_cast<int>(areaResult.dirtySeeds.size());
                    diagnostics->areaDirtyClosure = static_cast<int>(areaResult.dirtyClosure.size());
                    diagnostics->areaRecomputedNodes = areaResult.numRecomputedNodes();
                }
            }

            const auto primalRebuildStart = std::chrono::steady_clock::now();
            minTree.build(structuralResult.reconstruction->clone(), false, adjacency);
            if (diagnostics != nullptr) {
                diagnostics->primalRebuildMicros =
                    elapsedMicros(primalRebuildStart, std::chrono::steady_clock::now());
            }

            if (diagnostics != nullptr) {
                diagnostics->acceptedAlgorithm2 = true;
                diagnostics->createdNodes = structuralResult.numCreatedNodes();
                diagnostics->absorbedNodes = structuralResult.numAbsorbedNodes();
                diagnostics->releasedNodes = static_cast<int>(structuralResult.releasedNodes.size());
                diagnostics->affectedLevels = static_cast<int>(structuralResult.affectedLevels.size());
                diagnostics->xiBulletNodes = static_cast<int>(structuralResult.xiBulletNodes.size());
                diagnostics->xiCircleNodes = static_cast<int>(structuralResult.xiCircleNodes.size());
            }

            return structuralResult.reconstruction;
        } catch (const std::runtime_error &e) {
            return useFallback(e.what(), false, true);
        }
    }
};

} // namespace debug::region
