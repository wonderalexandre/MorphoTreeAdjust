#pragma once

#include <algorithm>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../Region.hpp"
#include "ExplicitRegionProblem1Contract.hpp"
#include "ExplicitRegionProblem1Survey.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem1 {

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1LevelMaterializationRecord {
    PixelType level{};
    int componentSize = 0;
    bool isIntermediateLevel = false;
    bool hasWorkspaceCandidate = false;
    bool hasCanonicalRebuildNode = false;

    bool missingWorkspaceCandidate() const {
        return !hasWorkspaceCandidate;
    }

    bool missingCandidateButMaterializedInRebuild() const {
        return missingWorkspaceCandidate() && hasCanonicalRebuildNode;
    }

    bool missingCandidateImplicitInRebuild() const {
        return missingWorkspaceCandidate() && !hasCanonicalRebuildNode;
    }
};

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1LevelMaterializationAnalysis {
    bool analyzed = false;
    ExplicitRegionProblem1Contract<PixelType> contract;
    ExplicitRegionProblem1Contract<PixelType> postMaterializationContract;
    std::vector<ExplicitRegionProblem1LevelMaterializationRecord<PixelType>> levels;
    std::string message;

    bool hasPostMaterializationContract = false;
    bool materializationAttempted = false;
    bool materializationApplied = false;
    int materializedOwners = 0;
    int materializedNodes = 0;
    int intervalLevels = 0;
    int workspaceCandidateLevels = 0;
    int missingWorkspaceCandidateLevels = 0;
    int canonicalRebuildNodeLevels = 0;
    int missingCandidateMaterializedInRebuildLevels = 0;
    int missingIntermediateCandidateMaterializedInRebuildLevels = 0;
    int missingCandidateImplicitInRebuildLevels = 0;
    int missingIntermediateCandidateImplicitInRebuildLevels = 0;
};

struct ExplicitRegionProblem1LevelMaterializationSurveySummary {
    int connectedRegions = 0;
    int totalEvents = 0;
    int rejectedByContract = 0;
    int acceptedByPublishedProblem = 0;
    long long materializationApplied = 0;
    long long materializedOwners = 0;
    long long materializedNodes = 0;
    long long intervalLevels = 0;
    long long workspaceCandidateLevels = 0;
    long long missingWorkspaceCandidateLevels = 0;
    long long canonicalRebuildNodeLevels = 0;
    long long missingCandidateMaterializedInRebuildLevels = 0;
    long long missingIntermediateCandidateMaterializedInRebuildLevels = 0;
    long long missingCandidateImplicitInRebuildLevels = 0;
    long long missingIntermediateCandidateImplicitInRebuildLevels = 0;
};

namespace detail {

inline bool pixelPassesTreeThreshold(int value, int level, bool isMaxtree) {
    return isMaxtree ? value >= level : value <= level;
}

inline std::vector<PixelId> thresholdComponentContainingRegion(
    const ImageUInt8Ptr &image,
    AdjacencyRelation &adjacency,
    std::span<const PixelId> pixelsC,
    int level,
    bool isMaxtree) {
    if (!image || pixelsC.empty()) {
        return {};
    }

    const int domainSize = image->getSize();
    const PixelId seed = pixelsC.front();
    if (seed < 0 || seed >= domainSize ||
        !pixelPassesTreeThreshold(static_cast<int>((*image)[seed]), level, isMaxtree)) {
        return {};
    }

    std::vector<unsigned char> visited(static_cast<size_t>(domainSize), 0);
    std::vector<PixelId> stack;
    std::vector<PixelId> component;
    stack.push_back(seed);
    visited[static_cast<size_t>(seed)] = 1;

    while (!stack.empty()) {
        const PixelId current = stack.back();
        stack.pop_back();
        component.push_back(current);

        for (PixelId neighbor : adjacency.getNeighborPixels(current)) {
            if (neighbor < 0 || neighbor >= domainSize) {
                continue;
            }
            const size_t neighborIndex = static_cast<size_t>(neighbor);
            if (visited[neighborIndex] != 0) {
                continue;
            }
            if (!pixelPassesTreeThreshold(
                    static_cast<int>((*image)[neighbor]), level, isMaxtree)) {
                continue;
            }
            visited[neighborIndex] = 1;
            stack.push_back(neighbor);
        }
    }

    std::sort(component.begin(), component.end());
    return component;
}

inline bool treeHasNodeWithLevelAndSupport(const DynamicComponentTree &tree,
                                           int level,
                                           const std::vector<PixelId> &support) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isAlive(nodeId) || tree.getAltitude(nodeId) != level) {
            continue;
        }

        std::vector<PixelId> nodeSupport = tree.getPixelsOfCC(nodeId);
        std::sort(nodeSupport.begin(), nodeSupport.end());
        if (nodeSupport == support) {
            return true;
        }
    }
    return false;
}

template<typename PixelType = AltitudeType>
void accumulateLevelMaterializationAnalysis(
    ExplicitRegionProblem1LevelMaterializationSurveySummary &summary,
    const ExplicitRegionProblem1LevelMaterializationAnalysis<PixelType> &analysis) {
    if (!analysis.contract.acceptsPublishedProblem()) {
        ++summary.rejectedByContract;
        return;
    }

    ++summary.acceptedByPublishedProblem;
    if (analysis.materializationApplied) {
        ++summary.materializationApplied;
        summary.materializedOwners += analysis.materializedOwners;
        summary.materializedNodes += analysis.materializedNodes;
    }
    summary.intervalLevels += analysis.intervalLevels;
    summary.workspaceCandidateLevels += analysis.workspaceCandidateLevels;
    summary.missingWorkspaceCandidateLevels += analysis.missingWorkspaceCandidateLevels;
    summary.canonicalRebuildNodeLevels += analysis.canonicalRebuildNodeLevels;
    summary.missingCandidateMaterializedInRebuildLevels +=
        analysis.missingCandidateMaterializedInRebuildLevels;
    summary.missingIntermediateCandidateMaterializedInRebuildLevels +=
        analysis.missingIntermediateCandidateMaterializedInRebuildLevels;
    summary.missingCandidateImplicitInRebuildLevels +=
        analysis.missingCandidateImplicitInRebuildLevels;
    summary.missingIntermediateCandidateImplicitInRebuildLevels +=
        analysis.missingIntermediateCandidateImplicitInRebuildLevels;
}

} // namespace detail

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1LevelMaterializationAnalysis<PixelType>
analyzeExplicitRegionProblem1LevelMaterialization(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    std::span<const PixelId> pixelsC,
    PixelType b,
    bool assumeConnected = false) {
    ExplicitRegionProblem1LevelMaterializationAnalysis<PixelType> out;

    if (!input) {
        out.message = "level materialization analysis requires a valid input image";
        return out;
    }
    if (!adjacency) {
        out.message = "level materialization analysis requires a valid adjacency relation";
        return out;
    }

    DynamicComponentTree preDualTree(input->clone(), updateMaxTree, adjacency);
    out.contract = checkExplicitRegionProblem1Contract<PixelType>(
        preDualTree,
        *adjacency,
        pixelsC,
        b,
        assumeConnected);

    if (!out.contract.acceptsPublishedProblem()) {
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

    const ImageUInt8Ptr updatedImage =
        RegionLevelChangeRebuildUpdater::applyToImage(input, change);

    const ExplicitRegionProblem1Contract<PixelType> *activeContract = &out.contract;
    if (out.contract.requiresLocalMaterialization()) {
        out.materializationAttempted = true;
        try {
            const auto materializationPlan = analyzeRegionMaterialization(
                preDualTree,
                *adjacency,
                pixelsC);
            const auto materializationResult = materializePartialOwners(
                preDualTree,
                materializationPlan);
            out.materializationApplied = true;
            out.materializedOwners = materializationResult.numMaterializedOwners();
            out.materializedNodes = materializationResult.numCreatedNodes();
            out.postMaterializationContract =
                checkExplicitRegionProblem1Contract<PixelType>(
                    preDualTree,
                    *adjacency,
                    pixelsC,
                    b,
                    true);
            out.hasPostMaterializationContract = true;
            activeContract = &out.postMaterializationContract;
        } catch (const std::exception &error) {
            out.message = error.what();
            return out;
        }
    }

    change.a = activeContract->a;
    change.direction = activeContract->direction;
    change.assumeMaterializedAsProperParts =
        !activeContract->requiresLocalMaterialization();

    const RegionMergeWorkspace<PixelType> workspace =
        buildRegionMergeWorkspace(preDualTree, *adjacency, change);
    DynamicComponentTree rebuiltDualTree(updatedImage->clone(), updateMaxTree, adjacency);

    const bool isMaxtree = preDualTree.isMaxtree();
    const int firstAffectedLevel = isMaxtree
        ? static_cast<int>(activeContract->a) + 1
        : static_cast<int>(activeContract->a) - 1;
    const int lastAffectedLevel = static_cast<int>(b);
    const int step = isMaxtree ? 1 : -1;

    for (int level = firstAffectedLevel;
         isMaxtree ? level <= lastAffectedLevel : level >= lastAffectedLevel;
         level += step) {
        ExplicitRegionProblem1LevelMaterializationRecord<PixelType> record;
        record.level = static_cast<PixelType>(level);
        record.isIntermediateLevel = level != lastAffectedLevel;
        const auto *bucket = workspace.mergeNodesVMinus(record.level);
        record.hasWorkspaceCandidate = bucket != nullptr && !bucket->empty();

        const std::vector<PixelId> component =
            detail::thresholdComponentContainingRegion(
                updatedImage,
                *adjacency,
                pixelsC,
                level,
                isMaxtree);
        record.componentSize = static_cast<int>(component.size());
        record.hasCanonicalRebuildNode =
            detail::treeHasNodeWithLevelAndSupport(
                rebuiltDualTree,
                level,
                component);

        ++out.intervalLevels;
        if (record.hasWorkspaceCandidate) {
            ++out.workspaceCandidateLevels;
        } else {
            ++out.missingWorkspaceCandidateLevels;
        }
        if (record.hasCanonicalRebuildNode) {
            ++out.canonicalRebuildNodeLevels;
        }
        if (record.missingCandidateButMaterializedInRebuild()) {
            ++out.missingCandidateMaterializedInRebuildLevels;
            if (record.isIntermediateLevel) {
                ++out.missingIntermediateCandidateMaterializedInRebuildLevels;
            }
        }
        if (record.missingCandidateImplicitInRebuild()) {
            ++out.missingCandidateImplicitInRebuildLevels;
            if (record.isIntermediateLevel) {
                ++out.missingIntermediateCandidateImplicitInRebuildLevels;
            }
        }

        out.levels.push_back(record);
    }

    out.analyzed = true;
    out.message = "level materialization analysis completed";
    return out;
}

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1LevelMaterializationSurveySummary
surveySeededExplicitRegionProblem1LevelMaterializationEvents(
    const ImageUInt8Ptr &input,
    const AdjacencyRelationPtr &adjacency,
    bool updateMaxTree,
    int maxRegionSize,
    int maxRegions,
    PixelType minB,
    PixelType maxB) {
    ExplicitRegionProblem1LevelMaterializationSurveySummary summary;
    if (!input || !adjacency || minB > maxB) {
        return summary;
    }

    auto connectedRegions = detail::enumerateSeededConnectedRegions(
        *adjacency,
        input->getSize(),
        maxRegionSize,
        maxRegions);
    summary.connectedRegions = static_cast<int>(connectedRegions.size());

    for (const auto &region : connectedRegions) {
        for (int bValue = static_cast<int>(minB);
             bValue <= static_cast<int>(maxB);
             ++bValue) {
            ++summary.totalEvents;
            const auto analysis =
                analyzeExplicitRegionProblem1LevelMaterialization<PixelType>(
                    input,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    static_cast<PixelType>(bValue),
                    true);
            detail::accumulateLevelMaterializationAnalysis(summary, analysis);
        }
    }

    return summary;
}

} // namespace debug::problem1
