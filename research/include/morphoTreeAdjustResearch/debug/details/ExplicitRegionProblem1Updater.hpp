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

enum class ExplicitRegionProblem1UpdateStatus {
    RejectedByContract,
    KnownUpdaterLimitation,
    LocalMaterializationThrew,
    LocalUpdateApplied,
    LocalUpdateMismatchedRebuild,
    LocalUpdateThrew,
    RebuildFallbackApplied,
    RebuildFallbackThrew,
};

enum class ExplicitRegionProblem1FallbackPolicy {
    NoFallback,
    RebuildPair,
};

struct ExplicitRegionProblem1UpdateOptions {
    bool assumeConnected = false;
    bool allowLocalMaterialization = false;
    bool allowTargetNodeCreation = false;
    ExplicitRegionProblem1FallbackPolicy fallbackPolicy =
        ExplicitRegionProblem1FallbackPolicy::NoFallback;
    bool enforceLocalCorrectnessValidation = false;
    bool validateAgainstRebuild = false;
    bool validateStructureAgainstRebuild = false;
    bool runtimePostConditionValidation = true;
    bool requireCanonicalByConstruction = false;
};

inline ExplicitRegionProblem1UpdateOptions makeValidatedSafeExplicitRegionProblem1Options(
    ExplicitRegionProblem1FallbackPolicy fallbackPolicy =
        ExplicitRegionProblem1FallbackPolicy::NoFallback) {
    ExplicitRegionProblem1UpdateOptions options;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.fallbackPolicy = fallbackPolicy;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;
    options.requireCanonicalByConstruction = true;
    return options;
}

/**
 * @brief Rich debug payload for an explicit-region Problem 1 update attempt.
 * @details Clean callers should use `updateTreeWithExplicitRegion(...)` or the
 * public checked wrapper. This diagnostic payload is for tests, surveys,
 * oracles, and benchmark tooling that need contract, materialization, rebuild,
 * structural, or canonical-by-construction details.
 */
template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1UpdateDiagnostics {
    ExplicitRegionProblem1UpdateStatus status =
        ExplicitRegionProblem1UpdateStatus::RejectedByContract;
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
    int materializedOwners = 0;
    int materializedNodes = 0;
    bool rebuildValidationPerformed = false;
    bool rebuildMatched = false;
    bool structuralValidationPerformed = false;
    bool structuralMatchedRebuild = false;
    int sameLevelEdgesCollapsed = 0;
    int postEmptyDetachedChildlessNodesRequired = 0;
    int postSameLevelParentChildEdgesRequired = 0;
    int postSameLevelSiblingMergesRequired = 0;
    int postComparableSiblingMovesRequired = 0;
    int postCanonicalizationChangesRequired = 0;
    bool canonicalByConstructionPreflightPerformed = false;
    bool canonicalByConstructionAccepted = false;
    std::string canonicalByConstructionRejectReason;
    bool fallbackAttempted = false;
    bool fallbackApplied = false;
    ImageUInt8Ptr fallbackReconstruction;
    std::string structuralValidationMessage;
    std::string message;

    bool appliedLocalUpdate() const {
        return status == ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied ||
               status == ExplicitRegionProblem1UpdateStatus::LocalUpdateMismatchedRebuild;
    }

    bool appliedRebuildFallback() const {
        return status == ExplicitRegionProblem1UpdateStatus::RebuildFallbackApplied;
    }

    bool appliedAnyUpdate() const {
        return appliedLocalUpdate() || appliedRebuildFallback();
    }
};

enum class ExplicitRegionProblem1PreflightStatus {
    RejectedByContract,
    RejectedByLocalCorrectnessValidation,
    KnownUpdaterLimitation,
    LocalCorrectnessAnalysisThrew,
    ReadyForLocalUpdate,
};

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1PreflightResult {
    ExplicitRegionProblem1PreflightStatus status =
        ExplicitRegionProblem1PreflightStatus::RejectedByContract;
    ExplicitRegionProblem1Contract<PixelType> contract;
    ExplicitRegionProblem1Contract<PixelType> postMaterializationContract;
    bool hasPostMaterializationContract = false;
    bool materializationRequired = false;
    bool materializationAllowed = false;
    bool targetNodeCreationRequired = false;
    bool targetNodeCreationAllowed = false;
    bool localCorrectnessValidationPerformed = false;
    bool localCorrectnessValidationPassed = false;
    bool rejectedByLocalCorrectnessValidation = false;
    ExplicitRegionProblem1LocalCorrectnessValidation<PixelType>
        localCorrectnessValidation;
    std::string message;

    bool readyForLocalUpdate() const {
        return status == ExplicitRegionProblem1PreflightStatus::ReadyForLocalUpdate;
    }
};

template<typename PixelType = AltitudeType>
bool explicitRegionProblem1NeedsOnlyTargetNodeCreation(
    const ExplicitRegionProblem1Contract<PixelType> &contract) {
    return contract.requiresTargetNodeCreation();
}

template<typename PixelType = AltitudeType>
bool explicitRegionProblem1CanRunLocalUpdater(
    const ExplicitRegionProblem1Contract<PixelType> &contract,
    const ExplicitRegionProblem1UpdateOptions &options) {
    return contract.readyForLocalUpdater() ||
           (options.allowTargetNodeCreation &&
            explicitRegionProblem1NeedsOnlyTargetNodeCreation(contract));
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1PreflightResult<PixelType> preflightExplicitRegionProblem1Update(
    const DynamicComponentTree &minTree,
    const DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    const ExplicitRegionProblem1UpdateOptions &options = {}) {
    ExplicitRegionProblem1PreflightResult<PixelType> out;

    if (!adjacency) {
        out.message = "Problem 1 preflight requires a valid adjacency relation";
        return out;
    }

    const DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;
    const ImageUInt8Ptr preUpdateImage = dualTree.reconstructionImage();

    out.contract = checkExplicitRegionProblem1Contract<PixelType>(
        dualTree,
        *adjacency,
        pixelsC,
        b,
        options.assumeConnected);

    if (!out.contract.acceptsPublishedProblem()) {
        out.status = ExplicitRegionProblem1PreflightStatus::RejectedByContract;
        out.message = out.contract.rejectionReason;
        return out;
    }

    out.materializationRequired =
        out.contract.requiresLocalMaterialization();
    out.materializationAllowed = options.allowLocalMaterialization;

    ExplicitRegionProblem1LocalCorrectnessDiagnostic<PixelType> diagnostic;
    try {
        ExplicitRegionProblem1LocalCorrectnessValidationOptions validationOptions;
        validationOptions.assumeConnected = options.assumeConnected;
        validationOptions.allowLocalMaterialization = options.allowLocalMaterialization;
        validationOptions.allowTargetNodeCreation = options.allowTargetNodeCreation;
        diagnostic = diagnoseExplicitRegionProblem1LocalCorrectness<PixelType>(
            preUpdateImage,
            adjacency,
            updateMaxTree,
            pixelsC,
            b,
            validationOptions);
    } catch (const std::exception &error) {
        out.status =
            ExplicitRegionProblem1PreflightStatus::LocalCorrectnessAnalysisThrew;
        out.message = error.what();
        return out;
    }

    out.localCorrectnessValidationPerformed = true;
    out.localCorrectnessValidation = diagnostic.validation;
    out.localCorrectnessValidationPassed =
        diagnostic.validation.satisfiesIntegratedReleaseCriteria();
    out.rejectedByLocalCorrectnessValidation =
        !out.localCorrectnessValidationPassed;

    if (diagnostic.analysis.hasPostMaterializationContract) {
        out.postMaterializationContract =
            diagnostic.analysis.postMaterializationContract;
        out.hasPostMaterializationContract = true;
    }

    if (!out.localCorrectnessValidationPassed) {
        out.status =
            ExplicitRegionProblem1PreflightStatus::RejectedByLocalCorrectnessValidation;
        out.message = diagnostic.validation.message;
        return out;
    }

    const ExplicitRegionProblem1Contract<PixelType> &activeContract =
        diagnostic.analysis.activeContract();
    out.targetNodeCreationRequired =
        explicitRegionProblem1NeedsOnlyTargetNodeCreation(activeContract);
    out.targetNodeCreationAllowed =
        out.targetNodeCreationRequired && options.allowTargetNodeCreation;

    if (!explicitRegionProblem1CanRunLocalUpdater(activeContract, options)) {
        out.status = ExplicitRegionProblem1PreflightStatus::KnownUpdaterLimitation;
        out.message = activeContract.currentUpdaterLimitation;
        return out;
    }

    out.status = ExplicitRegionProblem1PreflightStatus::ReadyForLocalUpdate;
    out.message = "Problem 1 event passed the non-mutating local update preflight";
    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1UpdateDiagnostics<PixelType> updateExplicitRegionProblem1WithDiagnostics(
    DynamicComponentTree &minTree,
    DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    const ExplicitRegionProblem1UpdateOptions &options = {}) {
    ExplicitRegionProblem1UpdateDiagnostics<PixelType> out;

    if (!adjacency) {
        out.message = "Problem 1 updater requires a valid adjacency relation";
        return out;
    }

    DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;
    const ImageUInt8Ptr preUpdateImage = dualTree.reconstructionImage();

    out.contract = checkExplicitRegionProblem1Contract<PixelType>(
        dualTree,
        *adjacency,
        pixelsC,
        b,
        options.assumeConnected);

    if (!out.contract.acceptsPublishedProblem()) {
        out.status = ExplicitRegionProblem1UpdateStatus::RejectedByContract;
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

    if (options.validateAgainstRebuild ||
        options.validateStructureAgainstRebuild ||
        options.requireCanonicalByConstruction) {
        out.expectedImage = RegionLevelChangeRebuildUpdater::applyToImage(preUpdateImage, change);
        out.rebuildValidationPerformed = true;
    }

    if (options.enforceLocalCorrectnessValidation ||
        options.requireCanonicalByConstruction) {
        ExplicitRegionProblem1LocalCorrectnessValidationOptions validationOptions;
        validationOptions.assumeConnected = options.assumeConnected;
        validationOptions.allowLocalMaterialization = options.allowLocalMaterialization;
        validationOptions.allowTargetNodeCreation = options.allowTargetNodeCreation;
        const auto diagnostic =
            diagnoseExplicitRegionProblem1LocalCorrectness<PixelType>(
                preUpdateImage,
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
            out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
            out.rejectedByLocalCorrectnessValidation = true;
            out.message = diagnostic.validation.message;
            return out;
        }
    }

    const ExplicitRegionProblem1Contract<PixelType> *activeContract = &out.contract;

    if (options.requireCanonicalByConstruction) {
        out.canonicalByConstructionPreflightPerformed = true;
        if (activeContract->requiresLocalMaterialization() &&
            !options.allowLocalMaterialization) {
            out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
            out.canonicalByConstructionRejectReason =
                "canonical-by-construction update requires C-only local materialization for partial owners";
            out.message =
                "Problem 1 canonical-by-construction preflight rejected before mutation: " +
                out.canonicalByConstructionRejectReason;
            return out;
        }
    }

    if (activeContract->requiresLocalMaterialization() &&
        options.allowLocalMaterialization) {
        out.materializationAttempted = true;
        try {
            const auto materializationPlan = analyzeRegionMaterialization(
                dualTree,
                *adjacency,
                pixelsC);
            const auto materializationResult = options.requireCanonicalByConstruction
                ? materializePartialOwnersInsideOnly(dualTree, materializationPlan)
                : materializePartialOwners(dualTree, materializationPlan);
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
            out.status = ExplicitRegionProblem1UpdateStatus::LocalMaterializationThrew;
            out.message = error.what();
            return out;
        }

        if (options.requireCanonicalByConstruction &&
            activeContract->requiresLocalMaterialization()) {
            out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
            out.canonicalByConstructionRejectReason =
                "C-only materialization left partial owners in the prepared contract";
            out.message =
                "Problem 1 canonical-by-construction preparation failed after C-only materialization: " +
                out.canonicalByConstructionRejectReason;
            return out;
        }
    }

    out.targetNodeCreationRequired =
        explicitRegionProblem1NeedsOnlyTargetNodeCreation(*activeContract);
    out.targetNodeCreationAllowed =
        out.targetNodeCreationRequired && options.allowTargetNodeCreation;

    if (!explicitRegionProblem1CanRunLocalUpdater(*activeContract, options)) {
        out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
        if (options.requireCanonicalByConstruction) {
            out.canonicalByConstructionRejectReason =
                activeContract->currentUpdaterLimitation;
        }
        out.message = activeContract->currentUpdaterLimitation;
        return out;
    }

    if (options.requireCanonicalByConstruction) {
        out.canonicalByConstructionAccepted = true;
    }

    try {
        DualMinMaxTreeIncrementalFilter<PixelType> adjust(&minTree, &maxTree, *adjacency);
        adjust.setRuntimePostConditionValidationEnabled(options.runtimePostConditionValidation);
        if (options.requireCanonicalByConstruction) {
            updateTreeWithExplicitRegionWithPostDiagnostics(
                adjust,
                &dualTree,
                pixelsC,
                b,
                options.allowTargetNodeCreation);
        } else {
            updateTreeWithExplicitRegionChecked(
                adjust,
                &dualTree,
                *adjacency,
                pixelsC,
                b,
                options.allowTargetNodeCreation);
        }
        out.sameLevelEdgesCollapsed = adjust.getLastSameLevelEdgesCollapsed();
        out.postEmptyDetachedChildlessNodesRequired =
            adjust.getLastPostEmptyDetachedChildlessNodesRequired();
        out.postSameLevelParentChildEdgesRequired =
            adjust.getLastPostSameLevelParentChildEdgesRequired();
        out.postCanonicalizationChangesRequired =
            adjust.getLastPostCanonicalizationChangesRequired();
        out.localReconstruction = dualTree.reconstructionImage();
    } catch (const std::exception &error) {
        out.status = ExplicitRegionProblem1UpdateStatus::LocalUpdateThrew;
        out.message = error.what();
        return out;
    }

    if (options.requireCanonicalByConstruction &&
        (out.sameLevelEdgesCollapsed != 0 ||
         out.postCanonicalizationChangesRequired != 0)) {
        out.status = ExplicitRegionProblem1UpdateStatus::LocalUpdateMismatchedRebuild;
        out.canonicalByConstructionRejectReason =
            "accepted canonical-by-construction plan still required post-update canonicalization";
        out.message =
            "Problem 1 canonical-by-construction contract failed after mutation: " +
            out.canonicalByConstructionRejectReason;
        return out;
    }

    if (out.rebuildValidationPerformed &&
        !out.localReconstruction->isEqual(out.expectedImage)) {
        out.status = ExplicitRegionProblem1UpdateStatus::LocalUpdateMismatchedRebuild;
        out.message = "Problem 1 local update reconstruction differs from rebuild oracle";
        return out;
    }

    if (options.validateStructureAgainstRebuild) {
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
            out.status = ExplicitRegionProblem1UpdateStatus::LocalUpdateMismatchedRebuild;
            out.message = "Problem 1 local update reconstruction matches rebuild, but structure differs: " +
                          structuralComparison.message;
            return out;
        }
    }

    out.status = ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied;
    out.rebuildMatched = out.rebuildValidationPerformed;
    out.message = out.rebuildValidationPerformed
        ? "Problem 1 local update reconstruction matches rebuild oracle"
        : "Problem 1 local update applied";
    return out;
}

template<typename PixelType = AltitudeType>
void copyExplicitRegionProblem1PreflightToUpdateDiagnostics(
    const ExplicitRegionProblem1PreflightResult<PixelType> &preflight,
    ExplicitRegionProblem1UpdateDiagnostics<PixelType> &out) {
    out.contract = preflight.contract;
    if (preflight.hasPostMaterializationContract) {
        out.postMaterializationContract = preflight.postMaterializationContract;
        out.hasPostMaterializationContract = true;
    }
    out.targetNodeCreationRequired = preflight.targetNodeCreationRequired;
    out.targetNodeCreationAllowed = preflight.targetNodeCreationAllowed;
    out.localCorrectnessValidationPerformed =
        preflight.localCorrectnessValidationPerformed;
    out.localCorrectnessValidationPassed =
        preflight.localCorrectnessValidationPassed;
    out.rejectedByLocalCorrectnessValidation =
        preflight.rejectedByLocalCorrectnessValidation;
    out.localCorrectnessValidation = preflight.localCorrectnessValidation;
}

template<typename PixelType = AltitudeType>
RegionLevelChange<PixelType> makeExplicitRegionProblem1ChangeForFallback(
    const ExplicitRegionProblem1Contract<PixelType> &contract,
    std::span<const PixelId> pixelsC,
    PixelType b) {
    RegionLevelChange<PixelType> change;
    change.pixelsC.assign(pixelsC.begin(), pixelsC.end());
    change.a = contract.a;
    change.b = b;
    change.direction = contract.direction;
    change.assumeConnected = true;
    change.assumeMaterializedAsProperParts =
        !contract.requiresLocalMaterialization();
    return change;
}

template<typename PixelType = AltitudeType>
bool explicitRegionProblem1DescriptorMatchesContract(
    const RegionLevelChange<PixelType> &change,
    const ExplicitRegionProblem1Contract<PixelType> &contract,
    std::string *message = nullptr) {
    if (change.a != contract.a) {
        if (message != nullptr) {
            *message = "RegionLevelChange level a does not match the current Problem 1 contract";
        }
        return false;
    }
    if (change.b != contract.b) {
        if (message != nullptr) {
            *message = "RegionLevelChange level b does not match the current Problem 1 contract";
        }
        return false;
    }
    if (change.direction != contract.direction) {
        if (message != nullptr) {
            *message = "RegionLevelChange direction does not match the current Problem 1 contract";
        }
        return false;
    }
    if (message != nullptr) {
        *message = "RegionLevelChange descriptor matches the current Problem 1 contract";
    }
    return true;
}

template<typename PixelType = AltitudeType>
bool explicitRegionProblem1AnchoredDescriptorMatchesContract(
    const AnchoredRegionLevelChange<PixelType> &anchoredChange,
    const ExplicitRegionProblem1Contract<PixelType> &contract,
    std::string *message = nullptr) {
    if (!explicitRegionProblem1DescriptorMatchesContract(
            anchoredChange.change,
            contract,
            message)) {
        return false;
    }
    if (anchoredChange.nodeCa != contract.nodeCa) {
        if (message != nullptr) {
            *message = "AnchoredRegionLevelChange nodeCa does not match the current Problem 1 contract";
        }
        return false;
    }
    if (message != nullptr) {
        *message = "AnchoredRegionLevelChange descriptor matches the current Problem 1 contract";
    }
    return true;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1PreflightResult<PixelType> preflightExplicitRegionProblem1Update(
    const DynamicComponentTree &minTree,
    const DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    const RegionLevelChange<PixelType> &change,
    const ExplicitRegionProblem1UpdateOptions &options = {}) {
    auto out = preflightExplicitRegionProblem1Update<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        std::span<const PixelId>(change.pixelsC),
        change.b,
        options);

    if (!out.contract.acceptsPublishedProblem()) {
        return out;
    }

    std::string descriptorMessage;
    if (!explicitRegionProblem1DescriptorMatchesContract(
            change,
            out.contract,
            &descriptorMessage)) {
        out.status = ExplicitRegionProblem1PreflightStatus::RejectedByContract;
        out.message = descriptorMessage;
        out.rejectedByLocalCorrectnessValidation = false;
    }

    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1PreflightResult<PixelType> preflightExplicitRegionProblem1Update(
    const DynamicComponentTree &minTree,
    const DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    const AnchoredRegionLevelChange<PixelType> &anchoredChange,
    const ExplicitRegionProblem1UpdateOptions &options = {}) {
    auto out = preflightExplicitRegionProblem1Update<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        std::span<const PixelId>(anchoredChange.change.pixelsC),
        anchoredChange.change.b,
        options);

    if (!out.contract.acceptsPublishedProblem()) {
        return out;
    }

    std::string descriptorMessage;
    if (!explicitRegionProblem1AnchoredDescriptorMatchesContract(
            anchoredChange,
            out.contract,
            &descriptorMessage)) {
        out.status = ExplicitRegionProblem1PreflightStatus::RejectedByContract;
        out.message = descriptorMessage;
        out.rejectedByLocalCorrectnessValidation = false;
    }

    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1UpdateDiagnostics<PixelType> safeUpdateExplicitRegionProblem1WithDiagnostics(
    DynamicComponentTree &minTree,
    DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    const ExplicitRegionProblem1UpdateOptions &options = {}) {
    const auto preflight = preflightExplicitRegionProblem1Update<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        pixelsC,
        b,
        options);

    if (!preflight.readyForLocalUpdate()) {
        ExplicitRegionProblem1UpdateDiagnostics<PixelType> out;
        copyExplicitRegionProblem1PreflightToUpdateDiagnostics(preflight, out);
        out.message = preflight.message;

        if (options.fallbackPolicy == ExplicitRegionProblem1FallbackPolicy::RebuildPair &&
            preflight.contract.acceptsPublishedProblem()) {
            out.fallbackAttempted = true;
            try {
                const auto change = makeExplicitRegionProblem1ChangeForFallback<PixelType>(
                    preflight.contract,
                    pixelsC,
                    b);
                out.fallbackReconstruction =
                    RegionLevelChangeRebuildUpdater::updatePair(
                        minTree,
                        maxTree,
                        adjacency,
                        change);
                out.fallbackApplied = true;
                out.status = ExplicitRegionProblem1UpdateStatus::RebuildFallbackApplied;
                out.message =
                    "Problem 1 safe update used rebuild fallback after preflight rejection: " +
                    preflight.message;
                return out;
            } catch (const std::exception &error) {
                out.status = ExplicitRegionProblem1UpdateStatus::RebuildFallbackThrew;
                out.message = error.what();
                return out;
            }
        }

        switch (preflight.status) {
            case ExplicitRegionProblem1PreflightStatus::RejectedByContract:
                out.status = ExplicitRegionProblem1UpdateStatus::RejectedByContract;
                break;
            case ExplicitRegionProblem1PreflightStatus::LocalCorrectnessAnalysisThrew:
                out.status = ExplicitRegionProblem1UpdateStatus::LocalMaterializationThrew;
                break;
            case ExplicitRegionProblem1PreflightStatus::RejectedByLocalCorrectnessValidation:
            case ExplicitRegionProblem1PreflightStatus::KnownUpdaterLimitation:
            case ExplicitRegionProblem1PreflightStatus::ReadyForLocalUpdate:
                out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
                break;
        }
        return out;
    }

    ExplicitRegionProblem1UpdateOptions updateOptions = options;
    updateOptions.enforceLocalCorrectnessValidation = false;
    auto out = updateExplicitRegionProblem1WithDiagnostics<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        pixelsC,
        b,
        updateOptions);
    copyExplicitRegionProblem1PreflightToUpdateDiagnostics(preflight, out);
    return out;
}

/**
 * @brief Validated diagnostic facade for broader Problem 1 events.
 * @details This is the richest current debug entry point for explicit Problem
 * 1 experiments. It keeps the production `updateTree` path untouched, but wraps
 * `updateTreeWithExplicitRegion` with non-mutating H1-H8 preflight, opt-in
 * local materialization, opt-in level-`b` node creation, and reconstruction
 * plus structural validation against rebuild.
 *
 * For connected Problem 1 regions, the H6-H8 portion of the preflight is the
 * current debug release contract beyond the direct single-`nodeCa`
 * sweep:
 *
 * - H6 requires one general workspace root (`closureRoot`);
 * - H7 requires `nodeCa` to reach `closureRoot` through a same-level chain;
 * - H8 requires every owner-interval bucket gap to be explained by that
 *   same-level closure path or by direct `C`-only support.
 *
 * These checks allow multi-owner and multi-extremal connected events to be
 * represented by one `closureRoot`; they do not promote a multi-workspace
 * backend and they do not change the published CASF path.
 *
 * Events outside the published Problem 1 contract are rejected without
 * fallback. Events inside the published contract but not released by the local
 * preflight can optionally use `RebuildPair` fallback through `fallbackPolicy`.
 */
template<typename PixelType = AltitudeType>
ExplicitRegionProblem1UpdateDiagnostics<PixelType> validatedSafeUpdateExplicitRegionProblem1WithDiagnostics(
    DynamicComponentTree &minTree,
    DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    ExplicitRegionProblem1FallbackPolicy fallbackPolicy =
        ExplicitRegionProblem1FallbackPolicy::NoFallback) {
    const auto options =
        makeValidatedSafeExplicitRegionProblem1Options(fallbackPolicy);
    return safeUpdateExplicitRegionProblem1WithDiagnostics<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        pixelsC,
        b,
        options);
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1UpdateDiagnostics<PixelType> validatedSafeUpdateExplicitRegionProblem1WithDiagnostics(
    DynamicComponentTree &minTree,
    DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    const RegionLevelChange<PixelType> &change,
    ExplicitRegionProblem1FallbackPolicy fallbackPolicy =
        ExplicitRegionProblem1FallbackPolicy::NoFallback) {
    const auto options =
        makeValidatedSafeExplicitRegionProblem1Options(fallbackPolicy);
    const auto preflight = preflightExplicitRegionProblem1Update<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        change,
        options);

    if (!preflight.readyForLocalUpdate()) {
        ExplicitRegionProblem1UpdateDiagnostics<PixelType> out;
        copyExplicitRegionProblem1PreflightToUpdateDiagnostics(preflight, out);
        out.message = preflight.message;

        if (fallbackPolicy == ExplicitRegionProblem1FallbackPolicy::RebuildPair &&
            preflight.contract.acceptsPublishedProblem() &&
            preflight.status != ExplicitRegionProblem1PreflightStatus::RejectedByContract) {
            out.fallbackAttempted = true;
            try {
                out.fallbackReconstruction =
                    RegionLevelChangeRebuildUpdater::updatePair(
                        minTree,
                        maxTree,
                        adjacency,
                        change);
                out.fallbackApplied = true;
                out.status = ExplicitRegionProblem1UpdateStatus::RebuildFallbackApplied;
                out.message =
                    "Problem 1 descriptor safe update used rebuild fallback after preflight rejection: " +
                    preflight.message;
                return out;
            } catch (const std::exception &error) {
                out.status = ExplicitRegionProblem1UpdateStatus::RebuildFallbackThrew;
                out.message = error.what();
                return out;
            }
        }

        switch (preflight.status) {
            case ExplicitRegionProblem1PreflightStatus::RejectedByContract:
                out.status = ExplicitRegionProblem1UpdateStatus::RejectedByContract;
                break;
            case ExplicitRegionProblem1PreflightStatus::LocalCorrectnessAnalysisThrew:
                out.status = ExplicitRegionProblem1UpdateStatus::LocalMaterializationThrew;
                break;
            case ExplicitRegionProblem1PreflightStatus::RejectedByLocalCorrectnessValidation:
            case ExplicitRegionProblem1PreflightStatus::KnownUpdaterLimitation:
            case ExplicitRegionProblem1PreflightStatus::ReadyForLocalUpdate:
                out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
                break;
        }
        return out;
    }

    auto updateOptions = options;
    updateOptions.enforceLocalCorrectnessValidation = false;
    auto out = updateExplicitRegionProblem1WithDiagnostics<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        std::span<const PixelId>(change.pixelsC),
        change.b,
        updateOptions);
    copyExplicitRegionProblem1PreflightToUpdateDiagnostics(preflight, out);
    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1UpdateDiagnostics<PixelType> validatedSafeUpdateExplicitRegionProblem1WithDiagnostics(
    DynamicComponentTree &minTree,
    DynamicComponentTree &maxTree,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    const AnchoredRegionLevelChange<PixelType> &anchoredChange,
    ExplicitRegionProblem1FallbackPolicy fallbackPolicy =
        ExplicitRegionProblem1FallbackPolicy::NoFallback) {
    const auto options =
        makeValidatedSafeExplicitRegionProblem1Options(fallbackPolicy);
    const auto preflight = preflightExplicitRegionProblem1Update<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        anchoredChange,
        options);

    if (!preflight.readyForLocalUpdate()) {
        ExplicitRegionProblem1UpdateDiagnostics<PixelType> out;
        copyExplicitRegionProblem1PreflightToUpdateDiagnostics(preflight, out);
        out.message = preflight.message;

        if (fallbackPolicy == ExplicitRegionProblem1FallbackPolicy::RebuildPair &&
            preflight.contract.acceptsPublishedProblem() &&
            preflight.status != ExplicitRegionProblem1PreflightStatus::RejectedByContract) {
            out.fallbackAttempted = true;
            try {
                out.fallbackReconstruction =
                    RegionLevelChangeRebuildUpdater::updatePair(
                        minTree,
                        maxTree,
                        adjacency,
                        anchoredChange.change);
                out.fallbackApplied = true;
                out.status = ExplicitRegionProblem1UpdateStatus::RebuildFallbackApplied;
                out.message =
                    "Problem 1 anchored descriptor safe update used rebuild fallback after preflight rejection: " +
                    preflight.message;
                return out;
            } catch (const std::exception &error) {
                out.status = ExplicitRegionProblem1UpdateStatus::RebuildFallbackThrew;
                out.message = error.what();
                return out;
            }
        }

        switch (preflight.status) {
            case ExplicitRegionProblem1PreflightStatus::RejectedByContract:
                out.status = ExplicitRegionProblem1UpdateStatus::RejectedByContract;
                break;
            case ExplicitRegionProblem1PreflightStatus::LocalCorrectnessAnalysisThrew:
                out.status = ExplicitRegionProblem1UpdateStatus::LocalMaterializationThrew;
                break;
            case ExplicitRegionProblem1PreflightStatus::RejectedByLocalCorrectnessValidation:
            case ExplicitRegionProblem1PreflightStatus::KnownUpdaterLimitation:
            case ExplicitRegionProblem1PreflightStatus::ReadyForLocalUpdate:
                out.status = ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation;
                break;
        }
        return out;
    }

    auto updateOptions = options;
    updateOptions.enforceLocalCorrectnessValidation = false;
    auto out = updateExplicitRegionProblem1WithDiagnostics<PixelType>(
        minTree,
        maxTree,
        adjacency,
        updateMaxTree,
        std::span<const PixelId>(anchoredChange.change.pixelsC),
        anchoredChange.change.b,
        updateOptions);
    copyExplicitRegionProblem1PreflightToUpdateDiagnostics(preflight, out);
    return out;
}

} // namespace debug::problem1
