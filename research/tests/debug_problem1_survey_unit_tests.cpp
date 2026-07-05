#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1LocalCorrectnessValidator.hpp"
#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1Survey.hpp"

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

ImageUInt8Ptr make_pattern_image(int rows, int cols, int variant) {
    auto image = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            const int pixelId = row * cols + col;
            int value = 0;
            if (variant == 0) {
                value = (3 * row + 5 * col + ((row * col) % 4)) % 10;
            } else if (variant == 1) {
                value = ((row < rows / 2 ? 2 : 7) +
                         (col % 3) +
                         ((row + col) % 2)) % 10;
            } else {
                value = (pixelId * 7 +
                         row * 3 +
                         col * 11 +
                         ((pixelId * pixelId) % 5)) % 10;
            }
            (*image)[pixelId] = static_cast<uint8_t>(value);
        }
    }
    return image;
}

ImageUInt8Ptr make_prng_image(int rows, int cols, uint32_t seed) {
    auto image = ImageUInt8::create(rows, cols);
    uint32_t state = seed;
    for (int pixelId = 0; pixelId < rows * cols; ++pixelId) {
        state = state * 1664525u + 1013904223u;
        const int row = pixelId / cols;
        const int col = pixelId % cols;
        const uint32_t mixed =
            state ^
            static_cast<uint32_t>(row * 97 + col * 193 + pixelId * 389);
        (*image)[pixelId] = static_cast<uint8_t>(mixed % 16);
    }
    return image;
}

void requireConsistentSummary(const problem1_debug::ExplicitRegionProblem1SurveySummary &summary) {
    require(summary.connectedRegions > 0, "survey must enumerate connected regions");
    require(summary.totalEvents > 0, "survey must classify events");
    require(summary.totalEvents == summary.rejectedByContract + summary.acceptedByPublishedProblem,
            "every event must be rejected by contract or accepted by the published problem");
    require(summary.acceptedByPublishedProblem ==
                summary.knownUpdaterLimitations +
                summary.localUpdateMatchedRebuild +
                summary.localUpdateMismatchedRebuild +
                summary.localUpdateThrew +
                summary.localMaterializationThrew +
                summary.rebuildFallbackMatchedRebuild +
                summary.rebuildFallbackMismatchedRebuild +
                summary.rebuildFallbackThrew,
            "accepted events must be classified by oracle status");
    require(summary.knownUpdaterLimitations ==
                summary.requiresMaterialization +
                summary.requiresTargetNodeCreation +
                summary.otherKnownLimitations,
            "policy/updater limitations must be partitioned by cause");
    if (summary.acceptedByPublishedProblem > 0) {
        require(summary.affectedIntervalLevels >= summary.acceptedByPublishedProblem,
                "each accepted event must contribute at least one interval level");
        require(summary.missingLevelUnionCandidates <= summary.affectedIntervalLevels,
                "missing level candidates cannot exceed interval levels");
        require(summary.missingIntermediateLevelUnionCandidates <=
                    summary.missingLevelUnionCandidates,
                "missing intermediate candidates must be a subset of all missing candidates");
        require(summary.missingIntermediateLevelUnionCandidates <=
                    summary.affectedIntervalLevels - summary.acceptedByPublishedProblem,
                "intermediate gaps must exclude the target level of each accepted event");
        require(summary.safeUpdateApplied <= summary.safeUpdateAttempted,
                "safe applied events must be attempted events");
        require(summary.safeUpdateMatchedRebuild + summary.safeUpdateMismatchedRebuild <=
                    summary.safeUpdateApplied,
                "safe update rebuild checks must be applied events");
        require(summary.directUpdateApplied <= summary.directUpdateAttempted,
                "direct applied events must be attempted events");
        require(summary.directUpdateMatchedRebuild +
                    summary.directUpdateMismatchedRebuild <=
                    summary.directUpdateApplied,
                "direct rebuild checks must be applied events");
        require(summary.directUpdateThrew <= summary.directUpdateAttempted,
                "direct throws must be attempted direct updates");
        require(summary.rebuildFallbackAttempted <= summary.safeUpdateAttempted,
                "fallback attempts must be safe update attempts");
        require(summary.rebuildFallbackApplied <= summary.rebuildFallbackAttempted,
                "fallback applications must be attempted fallbacks");
        require(summary.rebuildFallbackMatchedRebuild +
                    summary.rebuildFallbackMismatchedRebuild <=
                    summary.rebuildFallbackApplied,
                "fallback rebuild checks must be applied fallbacks");
        require(summary.rebuildFallbackThrew <= summary.rebuildFallbackAttempted,
                "fallback throws must be attempted fallbacks");
        require(summary.localCorrectnessValidationPassed <=
                    summary.localCorrectnessValidationPerformed,
                "passed local-correctness checks must be checked events");
        require(summary.rejectedByLocalCorrectnessValidation <=
                    summary.localCorrectnessValidationPerformed,
                "local-correctness rejections must be checked events");
        require(summary.localCorrectnessRejectedH2 <=
                    summary.rejectedByLocalCorrectnessValidation,
                "H2 local-correctness rejections must be gate rejections");
        require(summary.localCorrectnessRejectedH5 <=
                    summary.rejectedByLocalCorrectnessValidation,
                "H5 local-correctness rejections must be gate rejections");
        require(summary.localCorrectnessRejectedH6 <=
                    summary.rejectedByLocalCorrectnessValidation,
                "H6 local-correctness rejections must be gate rejections");
        require(summary.localCorrectnessRejectedH7 <=
                    summary.rejectedByLocalCorrectnessValidation,
                "H7 local-correctness rejections must be gate rejections");
        require(summary.localCorrectnessRejectedH8 <=
                    summary.rejectedByLocalCorrectnessValidation,
                "H8 local-correctness rejections must be gate rejections");
    }
}

void requireConsistentMultiAnchorSummary(
    const problem1_debug::ExplicitRegionProblem1MultiAnchorSurveySummary &summary) {
    require(summary.connectedRegions > 0,
            "multi-anchor survey must enumerate connected regions");
    require(summary.totalEvents > 0,
            "multi-anchor survey must classify events");
    require(summary.totalEvents == summary.rejectedByContract + summary.acceptedByPublishedProblem,
            "multi-anchor survey must partition events by contract status");
    require(summary.analyzedEvents <= summary.acceptedByPublishedProblem,
            "only accepted events can be analyzed for multi-anchor signals");
    if (summary.analyzedEvents > 0) {
        require(summary.owners >= summary.analyzedEvents,
                "each analyzed event must have at least one owner");
        require(summary.extremalOwners >= summary.analyzedEvents,
                "each analyzed event must have at least one extremal owner");
        require(summary.ownerPartitionComponents >= summary.owners,
                "owner partitions must contain at least one component per owner");
        require(summary.sparseWorkspaceAnchors >= summary.analyzedEvents,
                "each analyzed event must have at least one sparse workspace anchor");
        require(summary.ownerIntervalWorkspaceAnchors >= summary.analyzedEvents,
                "each analyzed event must have at least one owner interval workspace anchor");
        require(summary.generalWorkspaceAnchors >= summary.analyzedEvents,
                "each analyzed event must have at least one general workspace anchor");
        require(summary.multipleOwnerEvents <= summary.analyzedEvents,
                "multiple-owner events must be analyzed events");
        require(summary.multiOwnerSingleWorkspaceEvents <= summary.multipleOwnerEvents,
                "multi-owner single-workspace events must have multiple owners");
        require(summary.multiAnchorCandidateEvents <= summary.analyzedEvents,
                "multi-anchor candidates must be analyzed events");
        require(summary.multipleExtremalOwnerEvents <= summary.multipleOwnerEvents,
                "multiple extremal owners imply multiple owners");
        require(summary.multipleExtremalOwnerEvents <= summary.multiAnchorCandidateEvents,
                "multiple extremal owners imply a multi-anchor candidate");
        require(summary.disconnectedOwnerPartitionEvents <= summary.multiAnchorCandidateEvents,
                "disconnected owner partitions imply a multi-anchor candidate");
        require(summary.multipleSparseWorkspaceAnchorEvents <= summary.multiAnchorCandidateEvents,
                "multiple sparse workspace anchors imply a multi-anchor candidate");
        require(summary.multipleGeneralWorkspaceAnchorEvents <= summary.multiAnchorCandidateEvents,
                "multiple general workspace anchors imply a multi-anchor candidate");
        require(summary.ownerIntervalWorkspaceAnchors == summary.analyzedEvents,
                "the owner-interval workspace should have exactly one anchor per analyzed event");
        require(summary.generalWorkspaceAnchors == summary.analyzedEvents,
                "the general workspace should have exactly one anchor per analyzed event");
        require(summary.maxOwnerIntervalWorkspaceAnchors == 1,
                "the owner-interval workspace should not expose multiple local roots");
        require(summary.maxGeneralWorkspaceAnchors == 1,
                "the general workspace should not expose multiple local roots");
        require(summary.closureRootEvents == summary.analyzedEvents,
                "each analyzed event should expose one closure root");
        require(summary.closureRootEvents ==
                    summary.closureRootMatchesNodeCaEvents +
                    summary.closureRootDiffersFromNodeCaEvents,
                "closure root events must be partitioned by nodeCa agreement");
        require(summary.nodeCaReachesClosureRootEvents == summary.closureRootEvents,
                "nodeCa should reach closureRoot along the parent chain");
        require(summary.nodeCaToClosureRootSameLevelPathEvents == summary.closureRootEvents,
                "nodeCa-to-closureRoot paths should stay at the same level");
        require(summary.closureRootDiffersButSameLevelPathEvents ==
                    summary.closureRootDiffersFromNodeCaEvents,
                "closure-root mismatches should still be same-level paths");
        require(summary.multiOwnerClosureRootDiffersFromNodeCaEvents <=
                    summary.closureRootDiffersFromNodeCaEvents,
                "multi-owner closure-root mismatches are a subset of all mismatches");
        require(summary.ownerIntervalBucketMissingNodes ==
                    summary.ownerIntervalBucketClosurePathMissingNodes +
                    summary.ownerIntervalBucketCOnlyProperPartMissingNodes +
                    summary.ownerIntervalBucketUnexplainedMissingNodes,
                "owner-interval bucket gaps must be partitioned by explanation");
        require(summary.ownerIntervalBucketUnexplainedGapEvents <=
                    summary.ownerIntervalBucketGapEvents,
                "unexplained bucket-gap events are a subset of all bucket-gap events");
        require(summary.ownerIntervalBucketUnexplainedGapEvents == 0,
                "the current Problem 1 survey corpus should not expose unexplained bucket gaps");
        require(summary.ownerIntervalBucketUnexplainedMissingNodes == 0,
                "the current Problem 1 survey corpus should explain every missing bucket node");
        require(summary.maxOwnerIntervalBucketUnexplainedMissingNodes <=
                    summary.maxOwnerIntervalBucketMissingNodes,
                "unexplained missing bucket nodes are a subset of all missing bucket nodes");
        require(summary.h6UniqueGeneralWorkspaceRootEvents == summary.analyzedEvents,
                "all analyzed events should satisfy H6 in the current survey corpus");
        require(summary.h7SameLevelClosureEvents == summary.analyzedEvents,
                "all analyzed events should satisfy H7 in the current survey corpus");
        require(summary.h8BucketCompletenessEvents == summary.analyzedEvents,
                "all analyzed events should satisfy H8 in the current survey corpus");
        require(summary.targetLevelPresentOrCreatedEvents == summary.analyzedEvents,
                "all analyzed events should either have or allow the target level");
        require(summary.localCorrectnessCandidateEvents == summary.analyzedEvents,
                "all analyzed events should satisfy the integrated local correctness criteria");
    }
}

void testSurveyMapsMaxTreeProblem1Limitations() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
        image,
        adjacency,
        true,
        2,
        static_cast<AltitudeType>(0),
        static_cast<AltitudeType>(7));

    requireConsistentSummary(summary);
    require(summary.acceptedByPublishedProblem > 0,
            "max-tree survey must find events accepted by the published Problem 1");
    require(summary.requiresMaterialization > 0,
            "max-tree survey must expose C regions requiring materialization");
    require(summary.requiresTargetNodeCreation > 0,
            "max-tree survey must expose events requiring a new level-b node");
    require(!summary.hasLocalUpdateFailure(),
            "current compatible max-tree events should not mismatch or throw in the survey");
}

void testSurveyMapsMinTreeDualProblem1LimitationsAndCompatibleEvents() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
        image,
        adjacency,
        false,
        2,
        static_cast<AltitudeType>(0),
        static_cast<AltitudeType>(7));

    requireConsistentSummary(summary);
    require(summary.acceptedByPublishedProblem > 0,
            "min-tree dual survey must find events accepted by the dual Problem 1");
    require(summary.localUpdateMatchedRebuild > 0,
            "min-tree dual survey must find compatible events already handled by updateTreeWithExplicitRegion");
    require(summary.requiresMaterialization > 0,
            "min-tree dual survey must expose C regions requiring materialization");
    require(summary.requiresTargetNodeCreation > 0,
            "min-tree dual survey must expose events requiring a new level-b node");
    require(!summary.hasLocalUpdateFailure(),
            "current compatible min-tree events should not mismatch or throw in the survey");
}

void testSurveyWithMaterializationResolvesSomePartialOwnerEvents() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
        image,
        adjacency,
        false,
        1,
        static_cast<AltitudeType>(0),
        static_cast<AltitudeType>(7),
        true);

    requireConsistentSummary(summary);
    require(summary.materializationAttempted > 0,
            "survey should attempt materialization when the option is enabled");
    require(summary.materializationApplied > 0,
            "survey should apply materialization to at least one accepted event");
    require(summary.materializedLocalUpdateMatchedRebuild > 0,
            "survey should find events resolved by materialize-then-update");
    require(!summary.hasLocalUpdateFailure(),
            "materialized compatible events should not mismatch, throw, or fail materialization");
}

void testSurveyWithTargetNodeCreationResolvesSomeLevelBEvents() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    for (bool updateMaxTree : {true, false}) {
        const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
            image,
            adjacency,
            updateMaxTree,
            2,
            static_cast<AltitudeType>(0),
            static_cast<AltitudeType>(7),
            false,
            true);

        requireConsistentSummary(summary);
        require(summary.targetNodeCreationLocalUpdateMatchedRebuild > 0,
                "survey should find target-node-creation events resolved by the opt-in updater");
        require(!summary.hasLocalUpdateFailure(),
                "target-node-creation events should not mismatch or throw in the survey");
    }
}

void testSurveyCanGateLocalUpdatesByIntegratedCorrectness() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
        image,
        adjacency,
        true,
        2,
        static_cast<AltitudeType>(0),
        static_cast<AltitudeType>(7),
        true,
        false,
        false,
        true);

    requireConsistentSummary(summary);
    require(summary.acceptedByPublishedProblem > 0,
            "gated survey should find accepted Problem 1 events");
    require(summary.localCorrectnessValidationPerformed == summary.acceptedByPublishedProblem,
            "gated survey should validate every accepted event");
    require(summary.rejectedByLocalCorrectnessValidation > 0,
            "gated survey should reject level-b creation events when target creation is disabled");
    require(summary.localCorrectnessRejectedH5 ==
                summary.rejectedByLocalCorrectnessValidation,
            "disabled target-node creation should explain these policy rejections by H5");
    require(summary.localCorrectnessValidationPassed == summary.localUpdateMatchedRebuild,
            "events passing the gate should be exactly the events run to rebuild match");
    require(!summary.hasLocalUpdateFailure(),
            "gated survey should not introduce local update failures");
}

void testSurveyWithMaterializationAndTargetNodeCreationResolvesCombinedEvents() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
        image,
        adjacency,
        false,
        2,
        static_cast<AltitudeType>(0),
        static_cast<AltitudeType>(7),
        true,
        true);

    requireConsistentSummary(summary);
    require(summary.materializedTargetNodeCreationLocalUpdateMatchedRebuild > 0,
            "survey should find events resolved by materialization plus target-node creation");
    require(!summary.hasLocalUpdateFailure(),
            "combined materialization plus target-node-creation events should not mismatch or throw");
}

void testSurveyWithBothOptInsCoversAllSmallPublishedEvents() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    for (bool updateMaxTree : {true, false}) {
        const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
            image,
            adjacency,
            updateMaxTree,
            3,
            static_cast<AltitudeType>(0),
            static_cast<AltitudeType>(7),
            true,
            true);

        requireConsistentSummary(summary);
        require(summary.acceptedByPublishedProblem > 0,
                "survey should find published Problem 1 events");
        require(summary.knownUpdaterLimitations == 0,
                "materialization and target creation enabled should handle small published events");
        require(summary.acceptedByPublishedProblem == summary.localUpdateMatchedRebuild,
                "all small published events should match rebuild with both opt-ins enabled");
        require(!summary.hasLocalUpdateFailure(),
                "both opt-ins should not introduce mismatch, throw, or materialization failure");
    }
}

void testSafeSurveyFallbackResolvesGateRejectedTargetNodeEvents() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    int totalFallbackMatched = 0;
    int totalLocalMatched = 0;
    for (bool updateMaxTree : {true, false}) {
        const auto summary = problem1_debug::surveyExplicitRegionProblem1Events<AltitudeType>(
            image,
            adjacency,
            updateMaxTree,
            2,
            static_cast<AltitudeType>(0),
            static_cast<AltitudeType>(7),
            true,
            false,
            true,
            true,
            true,
            true);

        requireConsistentSummary(summary);
        require(summary.safeUpdateAttempted == summary.acceptedByPublishedProblem,
                "safe survey should attempt every accepted Problem 1 event");
        require(summary.rebuildFallbackAttempted > 0,
                "safe survey should attempt fallback for events blocked by disabled target-node creation");
        require(summary.rebuildFallbackApplied == summary.rebuildFallbackMatchedRebuild,
                "every applied fallback should match rebuild");
        require(summary.rebuildFallbackMismatchedRebuild == 0,
                "fallback path should not mismatch rebuild");
        require(summary.rebuildFallbackThrew == 0,
                "fallback path should not throw");
        require(summary.knownUpdaterLimitations == 0,
                "enabled rebuild fallback should cover policy gate rejections");
        require(summary.acceptedByPublishedProblem ==
                    summary.localUpdateMatchedRebuild +
                    summary.rebuildFallbackMatchedRebuild,
                "safe local plus fallback should cover all accepted events");
        require(summary.acceptedByPublishedProblem == summary.safeUpdateMatchedRebuild,
                "all safe updates should match rebuild");
        require(summary.structuralValidationPerformed == summary.acceptedByPublishedProblem,
                "safe survey should structurally validate local and fallback events");
        require(summary.structuralMatchedRebuild == summary.acceptedByPublishedProblem,
                "safe survey should structurally match rebuild for all accepted events");
        require(!summary.hasLocalUpdateFailure(),
                "safe fallback survey should not report update failures");

        totalFallbackMatched += summary.rebuildFallbackMatchedRebuild;
        totalLocalMatched += summary.localUpdateMatchedRebuild;
    }

    require(totalFallbackMatched > 0,
            "combined safe survey should validate at least one rebuild fallback");
    require(totalLocalMatched > 0,
            "combined safe survey should still validate at least one local update");
}

void testSeededSurveyCoversLargerSyntheticImagesWithBothOptIns() {
    for (int variant = 0; variant < 3; ++variant) {
        for (const auto [rows, cols] : {std::pair<int, int>{5, 5},
                                        std::pair<int, int>{6, 6}}) {
            auto image = make_pattern_image(rows, cols, variant);
            auto adjacency = std::make_shared<AdjacencyRelation>(
                image->getNumRows(),
                image->getNumCols(),
                1.5);

            for (bool updateMaxTree : {true, false}) {
                const auto summary =
                    problem1_debug::surveySeededExplicitRegionProblem1Events<AltitudeType>(
                        image,
                        adjacency,
                        updateMaxTree,
                        6,
                        240,
                        static_cast<AltitudeType>(0),
                        static_cast<AltitudeType>(10),
                        true,
                        true);

                requireConsistentSummary(summary);
                require(summary.connectedRegions > 0,
                        "seeded survey should generate connected regions");
                require(summary.acceptedByPublishedProblem > 0,
                        "seeded survey should find accepted Problem 1 events");
                require(summary.knownUpdaterLimitations == 0,
                        "materialization and target creation enabled should handle the seeded survey");
                require(summary.acceptedByPublishedProblem == summary.localUpdateMatchedRebuild,
                        "all seeded accepted events should match rebuild with both opt-ins enabled");
                require(!summary.hasLocalUpdateFailure(),
                        "seeded survey should not mismatch, throw, or fail materialization");
            }
        }
    }
}

void testSeededSurveyCoversPseudoRandomImagesWithBothOptIns() {
    int totalAccepted = 0;
    int totalMatched = 0;

    for (uint32_t seed : {17u, 101u, 4099u}) {
        for (const auto [rows, cols] : {std::pair<int, int>{5, 7},
                                        std::pair<int, int>{7, 5},
                                        std::pair<int, int>{7, 7}}) {
            auto image = make_prng_image(rows, cols, seed);
            auto adjacency = std::make_shared<AdjacencyRelation>(
                image->getNumRows(),
                image->getNumCols(),
                1.5);

            for (bool updateMaxTree : {true, false}) {
                const auto summary =
                    problem1_debug::surveySeededExplicitRegionProblem1Events<AltitudeType>(
                        image,
                        adjacency,
                        updateMaxTree,
                        7,
                        260,
                        static_cast<AltitudeType>(0),
                        static_cast<AltitudeType>(15),
                        true,
                        true);

                requireConsistentSummary(summary);
                require(summary.connectedRegions > 0,
                        "pseudo-random seeded survey should generate connected regions");
                require(summary.acceptedByPublishedProblem > 0,
                        "pseudo-random seeded survey should find accepted Problem 1 events");
                require(summary.knownUpdaterLimitations == 0,
                        "materialization and target creation enabled should handle pseudo-random surveys");
                require(summary.acceptedByPublishedProblem == summary.localUpdateMatchedRebuild,
                        "all pseudo-random accepted events should match rebuild");
                require(!summary.hasLocalUpdateFailure(),
                        "pseudo-random survey should not mismatch, throw, or fail materialization");

                totalAccepted += summary.acceptedByPublishedProblem;
                totalMatched += summary.localUpdateMatchedRebuild;
            }
        }
    }

    require(totalAccepted >= 15000,
            "pseudo-random seeded survey should cover a large event set");
    require(totalAccepted == totalMatched,
            "all pseudo-random accepted events should match rebuild in aggregate");
}

void testSeededSurveyCoversEightByEightPseudoRandomImagesWithBothOptIns() {
    int totalAccepted = 0;
    int totalMatched = 0;

    for (uint32_t seed : {23u, 911u}) {
        auto image = make_prng_image(8, 8, seed);
        auto adjacency = std::make_shared<AdjacencyRelation>(
            image->getNumRows(),
            image->getNumCols(),
            1.5);

        for (bool updateMaxTree : {true, false}) {
            const auto summary =
                problem1_debug::surveySeededExplicitRegionProblem1Events<AltitudeType>(
                    image,
                    adjacency,
                    updateMaxTree,
                    8,
                    340,
                    static_cast<AltitudeType>(0),
                    static_cast<AltitudeType>(19),
                    true,
                    true);

            requireConsistentSummary(summary);
            require(summary.connectedRegions > 0,
                    "8x8 pseudo-random seeded survey should generate connected regions");
            require(summary.acceptedByPublishedProblem > 0,
                    "8x8 pseudo-random seeded survey should find accepted events");
            require(summary.knownUpdaterLimitations == 0,
                    "materialization and target creation enabled should handle 8x8 surveys");
            require(summary.acceptedByPublishedProblem == summary.localUpdateMatchedRebuild,
                    "all 8x8 pseudo-random accepted events should match rebuild");
            require(!summary.hasLocalUpdateFailure(),
                    "8x8 pseudo-random survey should not mismatch, throw, or fail materialization");

            totalAccepted += summary.acceptedByPublishedProblem;
            totalMatched += summary.localUpdateMatchedRebuild;
        }
    }

    require(totalAccepted >= 6500,
            "8x8 pseudo-random seeded survey should cover a large event set");
    require(totalAccepted == totalMatched,
            "all 8x8 pseudo-random accepted events should match rebuild in aggregate");
}

void testMultiAnchorSurveyRecordsAcceptedEventGeometry() {
    auto image = make_prng_image(6, 6, 17u);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.5);

    int totalAnalyzed = 0;
    int totalMultiAnchorCandidates = 0;
    int totalMultipleOwnerEvents = 0;
    int totalMultiOwnerSingleWorkspaceEvents = 0;
    int totalMultipleGeneralWorkspaceAnchorEvents = 0;
    int totalClosureRootDiffersFromNodeCaEvents = 0;
    for (bool updateMaxTree : {true, false}) {
        const auto summary =
            problem1_debug::surveySeededExplicitRegionProblem1MultiAnchorEvents<AltitudeType>(
                image,
                adjacency,
                updateMaxTree,
                6,
                180,
                static_cast<AltitudeType>(0),
                static_cast<AltitudeType>(15),
                true);

        requireConsistentMultiAnchorSummary(summary);
        require(summary.acceptedByPublishedProblem > 0,
                "multi-anchor survey must find accepted Problem 1 events");
        require(summary.analyzedEvents == summary.acceptedByPublishedProblem,
                "all accepted events should be analyzed by the diagnostic");
        totalAnalyzed += summary.analyzedEvents;
        totalMultiAnchorCandidates += summary.multiAnchorCandidateEvents;
        totalMultipleOwnerEvents += summary.multipleOwnerEvents;
        totalMultiOwnerSingleWorkspaceEvents += summary.multiOwnerSingleWorkspaceEvents;
        totalMultipleGeneralWorkspaceAnchorEvents += summary.multipleGeneralWorkspaceAnchorEvents;
        totalClosureRootDiffersFromNodeCaEvents += summary.closureRootDiffersFromNodeCaEvents;
    }

    require(totalAnalyzed > 0,
            "multi-anchor diagnostic should analyze at least one event");
    require(totalMultiAnchorCandidates > 0,
            "seeded synthetic survey should expose at least one multi-anchor signal");
    require(totalMultipleOwnerEvents > 0,
            "seeded synthetic survey should expose events with multiple owners");
    require(totalMultiOwnerSingleWorkspaceEvents == totalMultipleOwnerEvents,
            "current multi-owner evidence should collapse to one general workspace");
    require(totalMultipleGeneralWorkspaceAnchorEvents == 0,
            "current synthetic corpus should not expose multiple general workspace roots");
    require(totalClosureRootDiffersFromNodeCaEvents > 0,
            "survey should keep covering closureRoot != nodeCa accepted events");
}

void testSafeSurveyWithGeneralProblem1OptInsMatchesStructure() {
    auto image = make_prng_image(6, 6, 23u);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.5);

    int totalAccepted = 0;
    for (bool updateMaxTree : {true, false}) {
        const auto summary =
            problem1_debug::surveySeededExplicitRegionProblem1Events<AltitudeType>(
                image,
                adjacency,
                updateMaxTree,
                5,
                160,
                static_cast<AltitudeType>(0),
                static_cast<AltitudeType>(15),
                true,
                true,
                true,
                true,
                true,
                false);

        requireConsistentSummary(summary);
        require(summary.acceptedByPublishedProblem > 0,
                "safe general survey should find accepted Problem 1 events");
        require(summary.knownUpdaterLimitations == 0,
                "materialization plus target-node creation should leave no updater limitations in this corpus");
        require(summary.localCorrectnessValidationPerformed ==
                    summary.acceptedByPublishedProblem,
                "safe general survey should validate H1-H8 for every accepted event");
        require(summary.localCorrectnessValidationPassed ==
                    summary.acceptedByPublishedProblem,
                "all accepted events should pass H1-H8 in this regression corpus");
        require(summary.rejectedByLocalCorrectnessValidation == 0,
                "H1-H8 should not reject events in this regression corpus");
        require(summary.safeUpdateAttempted == summary.acceptedByPublishedProblem,
                "safe updater should attempt every accepted event");
        require(summary.safeUpdateApplied == summary.acceptedByPublishedProblem,
                "safe updater should apply every accepted event");
        require(summary.safeUpdateMatchedRebuild == summary.acceptedByPublishedProblem,
                "safe updater should match rebuild for every accepted event");
        require(summary.structuralValidationPerformed == summary.acceptedByPublishedProblem,
                "safe updater should structurally validate every accepted event");
        require(summary.structuralMatchedRebuild == summary.acceptedByPublishedProblem,
                "safe updater should structurally match rebuild for every accepted event");
        require(!summary.hasLocalUpdateFailure(),
                "safe general survey should not mismatch, throw, or fail materialization");

        totalAccepted += summary.acceptedByPublishedProblem;
    }

    require(totalAccepted >= 1000,
            "safe general regression should cover a non-trivial accepted event set");
}

void testDirectMethodSurveyWithGeneralProblem1OptInsMatchesStructure() {
    auto image = make_prng_image(6, 6, 23u);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.5);

    int totalAccepted = 0;
    for (bool updateMaxTree : {true, false}) {
        const auto summary =
            problem1_debug::surveySeededExplicitRegionProblem1Events<AltitudeType>(
                image,
                adjacency,
                updateMaxTree,
                5,
                160,
                static_cast<AltitudeType>(0),
                static_cast<AltitudeType>(15),
                true,
                true,
                true,
                false,
                false,
                false,
                true);

        requireConsistentSummary(summary);
        require(summary.acceptedByPublishedProblem > 0,
                "direct-method survey should find accepted Problem 1 events");
        require(summary.safeUpdateAttempted == 0,
                "direct-method survey must not route through the safe facade");
        require(summary.directUpdateAttempted == summary.acceptedByPublishedProblem,
                "direct-method survey should attempt every accepted event released by H1-H8");
        require(summary.directUpdateApplied == summary.acceptedByPublishedProblem,
                "direct-method survey should apply every accepted event");
        require(summary.directUpdateMatchedRebuild == summary.acceptedByPublishedProblem,
                "direct updateTreeWithExplicitRegion should match rebuild for every accepted event");
        require(summary.directUpdateMismatchedRebuild == 0,
                "direct updateTreeWithExplicitRegion should not mismatch rebuild");
        require(summary.directUpdateThrew == 0,
                "direct updateTreeWithExplicitRegion should not throw in the released corpus");
        require(summary.localCorrectnessValidationPerformed ==
                    summary.acceptedByPublishedProblem,
                "direct-method survey should record H1-H8 validation for every accepted event");
        require(summary.localCorrectnessValidationPassed ==
                    summary.acceptedByPublishedProblem,
                "every accepted direct-method event should pass H1-H8 in this corpus");
        require(summary.structuralValidationPerformed == summary.acceptedByPublishedProblem,
                "direct-method survey should structurally validate every accepted event");
        require(summary.structuralMatchedRebuild == summary.acceptedByPublishedProblem,
                "direct-method survey should structurally match rebuild for every accepted event");
        require(summary.materializedTargetNodeCreationLocalUpdateMatchedRebuild > 0,
                "direct-method survey should cover materialization plus target-node creation");
        require(!summary.hasLocalUpdateFailure(),
                "direct-method survey should not mismatch, throw, or fail materialization");

        totalAccepted += summary.acceptedByPublishedProblem;
    }

    require(totalAccepted >= 1000,
            "direct-method regression should cover a non-trivial accepted event set");
}

void testLocalCorrectnessValidatorReportsTheoremHypotheses() {
    auto image = make_prng_image(6, 6, 17u);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.5);

    const std::vector<PixelId> region = {0, 1, 7};
    int validatedEvents = 0;
    int rejectedByTargetPolicy = 0;

    for (bool updateMaxTree : {true, false}) {
        for (int b = 0; b <= 15; ++b) {
            problem1_debug::ExplicitRegionProblem1LocalCorrectnessValidationOptions options;
            options.assumeConnected = true;
            options.allowLocalMaterialization = true;
            options.allowTargetNodeCreation = true;

            const auto diagnostic =
                problem1_debug::diagnoseExplicitRegionProblem1LocalCorrectness<AltitudeType>(
                    image,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    static_cast<AltitudeType>(b),
                    options);

            if (!diagnostic.validation.acceptedByPublishedProblem) {
                continue;
            }

            ++validatedEvents;
            require(diagnostic.validation.satisfiesH6ToH8(),
                    "accepted event should satisfy reusable H6-H8 validation");
            require(diagnostic.validation.satisfiesIntegratedReleaseCriteria(),
                    "accepted event should satisfy integrated release validation with opt-ins");

            options.allowTargetNodeCreation = false;
            const auto noTargetCreationDiagnostic =
                problem1_debug::diagnoseExplicitRegionProblem1LocalCorrectness<AltitudeType>(
                    image,
                    adjacency,
                    updateMaxTree,
                    std::span<const PixelId>(region),
                    static_cast<AltitudeType>(b),
                    options);
            if (noTargetCreationDiagnostic.validation.acceptedByPublishedProblem &&
                !noTargetCreationDiagnostic.validation.targetLevelPresent) {
                ++rejectedByTargetPolicy;
                require(!noTargetCreationDiagnostic.validation.satisfiesIntegratedReleaseCriteria(),
                        "missing target level with target creation disabled must fail integrated release validation");
            }
        }
    }

    require(validatedEvents > 0,
            "local correctness validator should find accepted events in the fixture");
    require(rejectedByTargetPolicy > 0,
            "local correctness validator should expose the target-level policy gate");
}

} // namespace

int main() {
    testSurveyMapsMaxTreeProblem1Limitations();
    testSurveyMapsMinTreeDualProblem1LimitationsAndCompatibleEvents();
    testSurveyWithMaterializationResolvesSomePartialOwnerEvents();
    testSurveyWithTargetNodeCreationResolvesSomeLevelBEvents();
    testSurveyCanGateLocalUpdatesByIntegratedCorrectness();
    testSurveyWithMaterializationAndTargetNodeCreationResolvesCombinedEvents();
    testSurveyWithBothOptInsCoversAllSmallPublishedEvents();
    testSafeSurveyFallbackResolvesGateRejectedTargetNodeEvents();
    testSeededSurveyCoversLargerSyntheticImagesWithBothOptIns();
    testSeededSurveyCoversPseudoRandomImagesWithBothOptIns();
    testSeededSurveyCoversEightByEightPseudoRandomImagesWithBothOptIns();
    testMultiAnchorSurveyRecordsAcceptedEventGeometry();
    testSafeSurveyWithGeneralProblem1OptInsMatchesStructure();
    testDirectMethodSurveyWithGeneralProblem1OptInsMatchesStructure();
    testLocalCorrectnessValidatorReportsTheoremHypotheses();
    return 0;
}
