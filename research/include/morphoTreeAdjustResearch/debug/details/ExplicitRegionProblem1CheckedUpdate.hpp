#pragma once

#include <span>
#include <stdexcept>
#include <string>

#include "../Region.hpp"
#include "ExplicitRegionProblem1Contract.hpp"
#include "ExplicitRegionProblem1LocalCorrectnessValidator.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DualMinMaxTreeIncrementalFilter.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem1 {

/**
 * @brief Checked diagnostic wrapper for the direct explicit-region Problem 1 update.
 * @details This function is intentionally kept outside
 * `DualMinMaxTreeIncrementalFilter`: it owns the expensive proof/debug layer
 * around the direct updater. It runs H1-H8/local-correctness diagnostics,
 * materializes partial owners when accepted, validates the prepared contract,
 * and then calls the adjuster's direct local applier with explicit cleanups
 * enabled.
 */
template<typename PixelType = AltitudeType>
void updateTreeWithExplicitRegionChecked(
    DualMinMaxTreeIncrementalFilter<PixelType> &adjust,
    DynamicComponentTree *dualTree,
    AdjacencyRelation &graph,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool allowTargetNodeCreation) {
    if (dualTree == nullptr) {
        throw std::runtime_error("updateTreeWithExplicitRegionChecked requires a valid dual tree");
    }
    if (pixelsC.empty()) {
        adjust.lastSameLevelEdgesCollapsed_ = 0;
        adjust.lastPostEmptyDetachedChildlessNodesRequired_ = 0;
        adjust.lastPostSameLevelParentChildEdgesRequired_ = 0;
        return;
    }
    if (dualTree != adjust.mintree_ && dualTree != adjust.maxtree_) {
        throw std::runtime_error("updateTreeWithExplicitRegionChecked requires one of the adjuster's paired trees");
    }
    if (adjust.graph_ == nullptr || adjust.graph_ != &graph) {
        throw std::runtime_error("updateTreeWithExplicitRegionChecked requires the adjuster's adjacency relation");
    }

    const bool isMaxtree = dualTree == adjust.maxtree_;
    const ImageUInt8Ptr preUpdateImage = dualTree->reconstructionImage();
    AdjacencyRelationPtr adjacencyAlias(&graph, [](AdjacencyRelation *) {});

    ExplicitRegionProblem1LocalCorrectnessValidationOptions validationOptions;
    validationOptions.assumeConnected = false;
    validationOptions.allowLocalMaterialization = true;
    validationOptions.allowTargetNodeCreation = allowTargetNodeCreation;
    const auto diagnostic =
        diagnoseExplicitRegionProblem1LocalCorrectness<PixelType>(
            preUpdateImage,
            adjacencyAlias,
            isMaxtree,
            pixelsC,
            b,
            validationOptions);
    if (!diagnostic.validation.satisfiesIntegratedReleaseCriteria()) {
        const std::string message = diagnostic.validation.message.empty()
            ? diagnostic.analysis.message
            : diagnostic.validation.message;
        throw std::runtime_error(
            "updateTreeWithExplicitRegionChecked rejected the Problem 1 event before mutation: " +
            message);
    }

    const auto materializationPlan =
        analyzeRegionMaterialization(
            *dualTree,
            graph,
            pixelsC);
    if (materializationPlan.requiresLocalMaterialization()) {
        materializePartialOwners(*dualTree, materializationPlan);
    }

    const auto preparedContract =
        checkExplicitRegionProblem1Contract<PixelType>(
            *dualTree,
            graph,
            pixelsC,
            b,
            true);
    switch (preparedContract.preflight.status) {
        case Problem1PreflightStatus::Accepted:
            break;
        case Problem1PreflightStatus::RequiresTargetNodeCreation:
            if (!allowTargetNodeCreation) {
                throw std::runtime_error(
                    "updateTreeWithExplicitRegionChecked requires level-b node creation for this event");
            }
            break;
        case Problem1PreflightStatus::RequiresLocalMaterialization:
            throw std::runtime_error(
                "updateTreeWithExplicitRegionChecked could not materialize all partial owners");
        case Problem1PreflightStatus::UnsupportedLocalTopology:
            throw std::runtime_error(
                "updateTreeWithExplicitRegionChecked rejected an unsupported post-materialization Problem 1 event: " +
                preparedContract.currentUpdaterLimitation);
        case Problem1PreflightStatus::EmptyRegion:
        case Problem1PreflightStatus::OutOfBoundsPixel:
        case Problem1PreflightStatus::DuplicatePixel:
        case Problem1PreflightStatus::MissingLiveOwner:
        case Problem1PreflightStatus::DisconnectedRegion:
        case Problem1PreflightStatus::InvalidTargetLevel:
            throw std::runtime_error(
                "updateTreeWithExplicitRegionChecked rejected the post-materialization Problem 1 event: " +
                preparedContract.rejectionReason);
    }

    adjust.applyExplicitRegionUpdate(
        dualTree,
        pixelsC,
        b,
        allowTargetNodeCreation,
        true,
        true,
        "updateTreeWithExplicitRegionChecked");
}

/**
 * @brief Direct Problem 1 updater with expensive post-canonicalization diagnostics.
 * @details This is a debug-only opt-in wrapper. The public direct updater keeps
 * the local mutation path clean; callers that need the post-update diagnostic
 * counters should route through this function instead.
 */
template<typename PixelType = AltitudeType>
void updateTreeWithExplicitRegionWithPostDiagnostics(
    DualMinMaxTreeIncrementalFilter<PixelType> &adjust,
    DynamicComponentTree *dualTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool allowTargetNodeCreation) {
    adjust.applyExplicitRegionUpdate(
        dualTree,
        pixelsC,
        b,
        allowTargetNodeCreation,
        false,
        true,
        "debug::problem1::updateTreeWithExplicitRegionWithPostDiagnostics");
}

} // namespace debug::problem1
