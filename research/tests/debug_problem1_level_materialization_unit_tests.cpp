#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1LevelMaterialization.hpp"

namespace problem1_debug = debug::problem1;

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr make_demo_image() {
    auto image = ImageUInt8::create(4, 4);
    const std::vector<uint8_t> raw = {
        4, 4, 2, 1,
        4, 3, 2, 1,
        5, 5, 2, 0,
        5, 6, 6, 0,
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

NodeId firstNodeWithAtLeastProperParts(const DynamicComponentTree &tree,
                                       int minProperParts,
                                       int minAltitudeExclusive = -1) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (tree.getNumProperParts(nodeId) >= minProperParts &&
            tree.getAltitude(nodeId) > minAltitudeExclusive) {
            return nodeId;
        }
    }
    return InvalidNode;
}

std::vector<PixelId> collectProperParts(const DynamicComponentTree &tree, NodeId nodeId) {
    std::vector<PixelId> pixels;
    for (PixelId pixelId : tree.getProperParts(nodeId)) {
        pixels.push_back(pixelId);
    }
    return pixels;
}

void requireConsistentAnalysis(
    const problem1_debug::ExplicitRegionProblem1LevelMaterializationAnalysis<AltitudeType> &analysis) {
    require(analysis.analyzed, "analysis should run for accepted Problem 1 events");
    require(analysis.contract.acceptedByPublishedProblem,
            "analysis fixture should satisfy the published Problem 1 contract");
    const auto &activeContract = analysis.hasPostMaterializationContract
        ? analysis.postMaterializationContract
        : analysis.contract;
    require(analysis.intervalLevels == activeContract.affectedIntervalLevels,
            "analysis interval must match the active contract interval");
    require(static_cast<int>(analysis.levels.size()) == analysis.intervalLevels,
            "analysis must expose one record per interval level");
    require(analysis.workspaceCandidateLevels +
                analysis.missingWorkspaceCandidateLevels ==
            analysis.intervalLevels,
            "workspace candidates and missing candidates must partition levels");
    require(analysis.missingCandidateMaterializedInRebuildLevels +
                analysis.missingCandidateImplicitInRebuildLevels ==
            analysis.missingWorkspaceCandidateLevels,
            "missing candidates must be materialized or implicit in the rebuild");
    require(analysis.missingIntermediateCandidateMaterializedInRebuildLevels <=
                analysis.missingCandidateMaterializedInRebuildLevels,
            "intermediate materialized gaps must be a subset of materialized gaps");
    require(analysis.missingIntermediateCandidateImplicitInRebuildLevels <=
                analysis.missingCandidateImplicitInRebuildLevels,
            "intermediate implicit gaps must be a subset of implicit gaps");
}

void testAnalyzesAcceptedMaxTreeEvent() {
    auto image = make_demo_image();
    auto adjacency =
        std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const std::vector<PixelId> pixelsC = {0, 1};

    const auto analysis =
        problem1_debug::analyzeExplicitRegionProblem1LevelMaterialization<AltitudeType>(
            image,
            adjacency,
            true,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(7),
            true);

    requireConsistentAnalysis(analysis);
    require(analysis.intervalLevels == 3,
            "max-tree fixture should analyze levels 5, 6, and 7");
    require(analysis.canonicalRebuildNodeLevels > 0,
            "analysis should identify at least one canonical rebuild node");
}

void testAnalyzesAcceptedMinTreeDualEvent() {
    auto image = make_demo_image();
    auto adjacency =
        std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const std::vector<PixelId> pixelsC = {13, 14};

    const auto analysis =
        problem1_debug::analyzeExplicitRegionProblem1LevelMaterialization<AltitudeType>(
            image,
            adjacency,
            false,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4),
            true);

    requireConsistentAnalysis(analysis);
    require(analysis.intervalLevels == 2,
            "min-tree fixture should analyze levels 5 and 4");
    require(analysis.canonicalRebuildNodeLevels > 0,
            "dual analysis should identify at least one canonical rebuild node");
}

void testAnalyzesAfterMaterializingPartialOwners() {
    auto image = make_demo_image();
    auto adjacency =
        std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2, 0);
    require(owner != InvalidNode,
            "test image must have a nonzero min-tree node with at least two proper parts");
    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};
    const auto b = static_cast<AltitudeType>(minTree.getAltitude(owner) - 1);

    const auto analysis =
        problem1_debug::analyzeExplicitRegionProblem1LevelMaterialization<AltitudeType>(
            image,
            adjacency,
            false,
            std::span<const PixelId>(pixelsC),
            b,
            true);

    requireConsistentAnalysis(analysis);
    require(analysis.materializationAttempted,
            "partial-owner event should attempt materialization before analysis");
    require(analysis.materializationApplied,
            "partial-owner event should be materialized before workspace analysis");
    require(analysis.hasPostMaterializationContract,
            "analysis should expose the post-materialization contract");
    require(analysis.contract.partialOwners > 0,
            "fixture should start with a partial owner");
    require(analysis.postMaterializationContract.partialOwners == 0,
            "post-materialization contract should remove partial owners");
    require(analysis.materializedOwners == 1,
            "fixture should materialize one owner");
    require(analysis.materializedNodes >= 1,
            "materialization should create at least one node");
}

void testSeededSurveyFindsImplicitIntermediateGaps() {
    auto image = make_demo_image();
    auto adjacency =
        std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto summary =
        problem1_debug::surveySeededExplicitRegionProblem1LevelMaterializationEvents<AltitudeType>(
            image,
            adjacency,
            true,
            3,
            30,
            static_cast<AltitudeType>(0),
            static_cast<AltitudeType>(7));

    require(summary.connectedRegions > 0,
            "level-materialization survey should enumerate connected regions");
    require(summary.totalEvents ==
                summary.rejectedByContract + summary.acceptedByPublishedProblem,
            "level-materialization survey must classify every event");
    require(summary.acceptedByPublishedProblem > 0,
            "level-materialization survey should find accepted events");
    require(summary.workspaceCandidateLevels +
                summary.missingWorkspaceCandidateLevels ==
            summary.intervalLevels,
            "survey candidates and missing candidates must partition interval levels");
    require(summary.missingCandidateMaterializedInRebuildLevels +
                summary.missingCandidateImplicitInRebuildLevels ==
            summary.missingWorkspaceCandidateLevels,
            "survey missing candidates must be materialized or implicit in rebuild");
    require(summary.missingIntermediateCandidateImplicitInRebuildLevels > 0,
            "fixture should expose at least one implicit intermediate gap");
}

} // namespace

int main() {
    testAnalyzesAcceptedMaxTreeEvent();
    testAnalyzesAcceptedMinTreeDualEvent();
    testAnalyzesAfterMaterializingPartialOwners();
    testSeededSurveyFindsImplicitIntermediateGaps();
    return 0;
}
