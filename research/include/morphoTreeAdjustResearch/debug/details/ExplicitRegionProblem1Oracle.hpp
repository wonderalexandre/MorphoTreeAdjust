#pragma once

#include <exception>
#include <span>
#include <string>

#include "../Region.hpp"
#include "../StructuralComparison.hpp"
#include "ExplicitRegionProblem1CheckedUpdate.hpp"
#include "ExplicitRegionProblem1Contract.hpp"
#include "ExplicitRegionProblem1LocalCorrectnessValidator.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DualMinMaxTreeIncrementalFilter.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem1 {

using ::debug::structural::compareDynamicComponentTreeStructure;

enum class ExplicitRegionProblem1OracleStatus {
    RejectedByContract,
    KnownUpdaterLimitation,
    LocalMaterializationThrew,
    LocalUpdateMatchedRebuild,
    LocalUpdateMismatchedRebuild,
    LocalUpdateThrew,
};

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1OracleResult {
    ExplicitRegionProblem1OracleStatus status =
        ExplicitRegionProblem1OracleStatus::RejectedByContract;
    ExplicitRegionProblem1Contract<PixelType> contract;
    ExplicitRegionProblem1Contract<PixelType> postMaterializationContract;
    ImageUInt8Ptr expectedImage;
    ImageUInt8Ptr localReconstruction;
    bool hasPostMaterializationContract = false;
    bool materializationAttempted = false;
    bool materializationApplied = false;
    bool targetNodeCreationRequired = false;
    bool targetNodeCreationAllowed = false;
    bool localCorrectnessValidationPerformed = false;
    bool localCorrectnessValidationPassed = false;
    bool rejectedByLocalCorrectnessValidation = false;
    ExplicitRegionProblem1LocalCorrectnessValidation<PixelType>
        localCorrectnessValidation;
    bool structuralValidationPerformed = false;
    bool structuralMatchedRebuild = false;
    int sameLevelEdgesCollapsed = 0;
    std::string structuralValidationMessage;
    int materializedOwners = 0;
    int materializedNodes = 0;
    std::string message;

    bool ranLocalUpdate() const {
        return status == ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild ||
               status == ExplicitRegionProblem1OracleStatus::LocalUpdateMismatchedRebuild ||
               status == ExplicitRegionProblem1OracleStatus::LocalUpdateThrew;
    }

    bool matchedRebuild() const {
        return status == ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild;
    }
};

template<typename PixelType = AltitudeType>
bool explicitRegionProblem1OracleNeedsOnlyTargetNodeCreation(
    const ExplicitRegionProblem1Contract<PixelType> &contract) {
    return contract.requiresTargetNodeCreation();
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1OracleResult<PixelType> runExplicitRegionProblem1Oracle(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false,
    bool allowLocalMaterialization = false,
    bool allowTargetNodeCreation = false,
    bool validateStructureAgainstRebuild = false,
    bool enforceLocalCorrectnessValidation = false) {
    ExplicitRegionProblem1OracleResult<PixelType> out;

    if (!input) {
        out.message = "oracle requires a valid input image";
        return out;
    }
    if (!adjacency) {
        out.message = "oracle requires a valid adjacency relation";
        return out;
    }

    DynamicComponentTree minTree(input->clone(), false, adjacency);
    DynamicComponentTree maxTree(input->clone(), true, adjacency);
    DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;

    out.contract = checkExplicitRegionProblem1Contract<PixelType>(
        dualTree,
        *adjacency,
        pixelsC,
        b,
        assumeConnected);

    if (!out.contract.acceptsPublishedProblem()) {
        out.status = ExplicitRegionProblem1OracleStatus::RejectedByContract;
        out.message = out.contract.rejectionReason;
        return out;
    }

    RegionLevelChange<PixelType> change;
    change.pixelsC.assign(pixelsC.begin(), pixelsC.end());
    change.a = out.contract.a;
    change.b = b;
    change.direction = out.contract.direction;
    change.assumeConnected = true;
    change.assumeMaterializedAsProperParts =
        !out.contract.requiresLocalMaterialization();

    out.expectedImage = RegionLevelChangeRebuildUpdater::applyToImage(input, change);

    if (enforceLocalCorrectnessValidation) {
        ExplicitRegionProblem1LocalCorrectnessValidationOptions validationOptions;
        validationOptions.assumeConnected = assumeConnected;
        validationOptions.allowLocalMaterialization = allowLocalMaterialization;
        validationOptions.allowTargetNodeCreation = allowTargetNodeCreation;
        const auto diagnostic =
            diagnoseExplicitRegionProblem1LocalCorrectness<PixelType>(
                input,
                adjacency,
                updateMaxTree,
                pixelsC,
                b,
                validationOptions);
        out.localCorrectnessValidationPerformed = true;
        out.localCorrectnessValidation = diagnostic.validation;
        out.localCorrectnessValidationPassed =
            diagnostic.validation.satisfiesIntegratedReleaseCriteria();
        if (!out.localCorrectnessValidationPassed) {
            out.status = ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation;
            out.rejectedByLocalCorrectnessValidation = true;
            out.message = diagnostic.validation.message;
            return out;
        }
    }

    const ExplicitRegionProblem1Contract<PixelType> *activeContract = &out.contract;

    if (activeContract->requiresLocalMaterialization() &&
        allowLocalMaterialization) {
        out.materializationAttempted = true;
        try {
            const auto materializationPlan = analyzeRegionMaterialization(
                dualTree,
                *adjacency,
                pixelsC);
            const auto materializationResult = materializePartialOwners(
                dualTree,
                materializationPlan);
            out.materializationApplied = true;
            out.materializedOwners = materializationResult.numMaterializedOwners();
            out.materializedNodes = materializationResult.numCreatedNodes();
            out.postMaterializationContract = checkExplicitRegionProblem1Contract<PixelType>(
                dualTree,
                *adjacency,
                pixelsC,
                b,
                true);
            out.hasPostMaterializationContract = true;
            activeContract = &out.postMaterializationContract;
        } catch (const std::exception &error) {
            out.status = ExplicitRegionProblem1OracleStatus::LocalMaterializationThrew;
            out.message = error.what();
            return out;
        }
    }

    out.targetNodeCreationRequired =
        explicitRegionProblem1OracleNeedsOnlyTargetNodeCreation(*activeContract);
    out.targetNodeCreationAllowed =
        out.targetNodeCreationRequired && allowTargetNodeCreation;

    const bool canRunLocalUpdate =
        activeContract->readyForLocalUpdater() ||
        out.targetNodeCreationAllowed;
    if (!canRunLocalUpdate) {
        out.status = ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation;
        out.message = activeContract->currentUpdaterLimitation;
        return out;
    }

    try {
        DualMinMaxTreeIncrementalFilter<PixelType> adjust(&minTree, &maxTree, *adjacency);
        adjust.setRuntimePostConditionValidationEnabled(true);
        updateTreeWithExplicitRegionChecked(
            adjust,
            &dualTree,
            *adjacency,
            pixelsC,
            b,
            allowTargetNodeCreation);
        out.sameLevelEdgesCollapsed = adjust.getLastSameLevelEdgesCollapsed();
        out.localReconstruction = dualTree.reconstructionImage();
    } catch (const std::exception &error) {
        out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateThrew;
        out.message = error.what();
        return out;
    }

    if (out.localReconstruction->isEqual(out.expectedImage)) {
        if (validateStructureAgainstRebuild) {
            DynamicComponentTree expectedTree(
                out.expectedImage->clone(),
                updateMaxTree,
                adjacency);
            const auto structuralComparison =
                compareDynamicComponentTreeStructure(dualTree, expectedTree);
            out.structuralValidationPerformed = true;
            out.structuralMatchedRebuild = structuralComparison.equivalent;
            out.structuralValidationMessage = structuralComparison.message;
            if (!structuralComparison.equivalent) {
                out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateMismatchedRebuild;
                out.message = "local update reconstruction matches rebuild, but structure differs: " +
                              structuralComparison.message;
                return out;
            }
        }

        out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild;
        out.message = "local update reconstruction matches rebuild oracle";
    } else {
        out.status = ExplicitRegionProblem1OracleStatus::LocalUpdateMismatchedRebuild;
        out.message = "local update reconstruction differs from rebuild oracle";
    }

    return out;
}

} // namespace debug::problem1
