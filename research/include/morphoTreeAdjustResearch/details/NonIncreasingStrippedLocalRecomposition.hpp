#pragma once

#include <stdexcept>
#include <utility>

#include "../AdjacencyRelation.hpp"
#include "../DynamicComponentTree.hpp"
#include "NonIncreasingLocalOracle.hpp"
#include "NonIncreasingLocalOracleDependency.hpp"
#include "NonIncreasingLocalPartition.hpp"
#include "NonIncreasingLocalPartitionPlanner.hpp"
#include "NonIncreasingLocalRecomposer.hpp"

namespace morphotree::details::non_increasing {

struct StrippedNonOracleLocalRecompositionBuild {
    LocalPartitionPlan partitionPlan;
    LocalRecompositionPlan recompositionPlan;
    LocalOracleDependencyReport dependencyReport;

    bool isOracleFree() const {
        return dependencyReport.isOracleFree();
    }
};

struct StrippedNonOracleLocalRecompositionMaterialization {
    StrippedNonOracleLocalRecompositionBuild build;
    LocalMergeOracle materializationContext;
    LocalMaterializationOracleDependencyReport materializationDependencyReport;
    LocalRecompositionMaterializationResult materialization;

    bool isOracleFree() const {
        return build.isOracleFree() &&
               materializationDependencyReport.isOracleFreeMaterializationInput();
    }
};

inline StrippedNonOracleLocalRecompositionBuild buildStrippedNonOracleLocalRecomposition(
    const LocalPartitionInput &input,
    AdjacencyRelation &adjacency) {
    StrippedNonOracleLocalRecompositionBuild result;
    result.partitionPlan = buildLocalPartitionPlanWithoutOracle(input, adjacency);
    result.recompositionPlan = buildLocalRecompositionPlanWithoutOracle(
        input,
        result.partitionPlan);
    result.dependencyReport = auditLocalPlanOracleDependencies(
        &result.partitionPlan,
        &result.recompositionPlan);

    if (!result.dependencyReport.isOracleFree()) {
        throw std::runtime_error("stripped non-oracle recomposition build exposed oracle-guided fields");
    }
    return result;
}

inline LocalMergeOracle makeStrippedNonOracleMaterializationContext(
    const DynamicComponentTree &dualTree,
    const LocalPartitionInput &input,
    bool includePreImage = true) {
    if (!localPartitionHasUpdatedTarget(input)) {
        throw std::runtime_error("stripped non-oracle materialization requires target values");
    }
    if (dualTree.isMaxtree() != input.contraction.dualIsMaxtree) {
        throw std::runtime_error("stripped non-oracle materialization received the wrong dual-tree polarity");
    }

    LocalMergeOracle context;
    context.contraction = input.contraction;
    if (includePreImage) {
        context.preImage = dualTree.reconstructionImage();
    }
    context.updatedImage = input.updatedImage;
    context.domainSize = localPartitionDomainSize(input);
    context.localPixels = input.localPixels;
    context.updatedLocalPixelValues = input.updatedLocalPixelValues;
    context.preDualNodes = input.preDualNodes;
    return context;
}

inline StrippedNonOracleLocalRecompositionMaterialization materializeStrippedNonOracleLocalRecomposition(
    DynamicComponentTree &dualTree,
    const LocalPartitionInput &input,
    AdjacencyRelation &partitionAdjacency,
    AdjacencyRelation *materializationAdjacency = nullptr) {
    StrippedNonOracleLocalRecompositionMaterialization result;
    result.build = buildStrippedNonOracleLocalRecomposition(
        input,
        partitionAdjacency);
    result.materializationContext = makeStrippedNonOracleMaterializationContext(
        dualTree,
        input);
    result.materializationDependencyReport = auditLocalMaterializationOracleDependencies(
        result.materializationContext,
        result.build.recompositionPlan);

    if (!result.materializationDependencyReport.isOracleFreeMaterializationInput()) {
        throw std::runtime_error("stripped non-oracle materialization exposed oracle-guided fields");
    }

    result.materialization = materializeLocalRecompositionNodes(
        dualTree,
        result.materializationContext,
        result.build.recompositionPlan,
        materializationAdjacency);
    return result;
}

inline StrippedNonOracleLocalRecompositionMaterialization materializeStrippedNonOracleLocalRecompositionPlan(
    DynamicComponentTree &dualTree,
    const LocalPartitionInput &input,
    LocalRecompositionPlan recompositionPlan,
    AdjacencyRelation *materializationAdjacency = nullptr) {
    StrippedNonOracleLocalRecompositionMaterialization result;
    result.build.recompositionPlan = std::move(recompositionPlan);
    result.build.dependencyReport = auditLocalPlanOracleDependencies(
        nullptr,
        &result.build.recompositionPlan);
    if (!result.build.dependencyReport.isOracleFree()) {
        throw std::runtime_error("stripped non-oracle recomposition plan exposed oracle-guided fields");
    }

    result.materializationContext = makeStrippedNonOracleMaterializationContext(
        dualTree,
        input);
    result.materializationDependencyReport = auditLocalMaterializationOracleDependencies(
        result.materializationContext,
        result.build.recompositionPlan);
    if (!result.materializationDependencyReport.isOracleFreeMaterializationInput()) {
        throw std::runtime_error("stripped non-oracle materialization plan exposed oracle-guided fields");
    }

    result.materialization = materializeLocalRecompositionNodes(
        dualTree,
        result.materializationContext,
        result.build.recompositionPlan,
        materializationAdjacency);
    return result;
}

} // namespace morphotree::details::non_increasing
