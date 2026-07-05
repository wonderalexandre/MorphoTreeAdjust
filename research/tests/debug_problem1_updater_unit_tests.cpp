#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1Updater.hpp"

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

ImageUInt8Ptr make_empty_detached_regression_image() {
    auto image = ImageUInt8::create(6, 6);
    const std::vector<uint8_t> raw = {
        0, 5, 0, 5, 0, 5,
        3, 9, 5, 1, 3, 9,
        6, 3, 6, 3, 6, 3,
        9, 7, 1, 5, 9, 7,
        2, 7, 2, 7, 2, 7,
        5, 1, 7, 3, 5, 1,
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_complex_image() {
    auto image = ImageUInt8::create(12, 12);
    const std::vector<uint8_t> raw = {
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
        2, 5, 5, 5, 3, 3, 3, 2, 2, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 5, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 6, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 5, 2, 3, 1,
        2, 5, 5, 5, 3, 3, 3, 2, 2, 2, 3, 1,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1,
        2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 1,
        6, 6, 6, 6, 6, 2, 4, 4, 4, 4, 4, 1,
        7, 3, 6, 4, 6, 2, 4, 8, 4, 8, 4, 1,
        7, 8, 6, 6, 6, 2, 4, 4, 4, 4, 4, 1,
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
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
            static_cast<uint32_t>(row * 97 + col * 193 + pixelId * 389 + seed * 17);
        (*image)[pixelId] = static_cast<uint8_t>(mixed % 20);
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

bool regionEqualsSubtreeFootprint(const DynamicComponentTree &tree,
                                  NodeId subtreeRoot,
                                  std::span<const PixelId> pixelsC) {
    std::vector<unsigned char> inRegion(
        static_cast<size_t>(tree.getNumTotalProperParts()),
        0);
    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
            return false;
        }
        inRegion[static_cast<size_t>(pixelId)] = 1;
    }

    int footprintSize = 0;
    for (NodeId subtreeNodeId : tree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            if (inRegion[static_cast<size_t>(pixelId)] == 0) {
                return false;
            }
            ++footprintSize;
        }
    }

    return footprintSize == static_cast<int>(pixelsC.size());
}

bool regionMatchesAnyNonRootSubtreeFootprint(const DynamicComponentTree &tree,
                                             std::span<const PixelId> pixelsC) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (tree.isRoot(nodeId)) {
            continue;
        }
        if (regionEqualsSubtreeFootprint(tree, nodeId, pixelsC)) {
            return true;
        }
    }
    return false;
}

void requireNoAliveEmptyNodes(const DynamicComponentTree &tree) {
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isNode(nodeId) &&
            tree.isAlive(nodeId) &&
            tree.getNumProperParts(nodeId) == 0) {
            throw std::runtime_error("tree should not contain alive nodes without proper parts");
        }
    }
}

struct TargetNodeCreationEvent {
    std::vector<PixelId> pixelsC;
    AltitudeType b = 0;
    bool updateMaxTree = false;
    bool found = false;
};

struct MaterializationEvent {
    std::vector<PixelId> pixelsC;
    AltitudeType b = 0;
    bool found = false;
};

void enumeratePixelCombinations(int domainSize,
                                int targetSize,
                                int start,
                                std::vector<PixelId> &current,
                                std::vector<std::vector<PixelId>> &out) {
    if (static_cast<int>(current.size()) == targetSize) {
        out.push_back(current);
        return;
    }

    const int remaining = targetSize - static_cast<int>(current.size());
    for (int pixelId = start; pixelId <= domainSize - remaining; ++pixelId) {
        current.push_back(static_cast<PixelId>(pixelId));
        enumeratePixelCombinations(domainSize, targetSize, pixelId + 1, current, out);
        current.pop_back();
    }
}

std::vector<std::vector<PixelId>> enumerateConnectedRegions(const ImageUInt8Ptr &image,
                                                           const AdjacencyRelationPtr &adjacency,
                                                           int maxRegionSize) {
    std::vector<std::vector<PixelId>> regions;
    const int domainSize = image->getSize();
    for (int regionSize = 1; regionSize <= maxRegionSize; ++regionSize) {
        std::vector<std::vector<PixelId>> candidates;
        std::vector<PixelId> current;
        current.reserve(static_cast<size_t>(regionSize));
        enumeratePixelCombinations(domainSize, regionSize, 0, current, candidates);

        for (auto &candidate : candidates) {
            if (problem1_debug::isConnectedRegion(*adjacency, std::span<const PixelId>(candidate))) {
                regions.push_back(std::move(candidate));
            }
        }
    }
    return regions;
}

TargetNodeCreationEvent findTargetNodeCreationEvent(const ImageUInt8Ptr &image,
                                                    const AdjacencyRelationPtr &adjacency,
                                                    bool updateMaxTree) {
    const auto regions = enumerateConnectedRegions(image, adjacency, 2);

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;

    for (const auto &region : regions) {
        for (int bValue = 0; bValue <= 7; ++bValue) {
            const auto contract =
                problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
                    dualTree,
                    *adjacency,
                    std::span<const PixelId>(region),
                    static_cast<AltitudeType>(bValue),
                    true);
            if (contract.acceptedByPublishedProblem &&
                contract.partialOwners == 0 &&
                !contract.hasTargetLevelUnionCandidate) {
                return TargetNodeCreationEvent{
                    region,
                    static_cast<AltitudeType>(bValue),
                    updateMaxTree,
                    true};
            }
        }
    }

    return TargetNodeCreationEvent{};
}

TargetNodeCreationEvent findMaterializationThenTargetNodeCreationEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    const auto regions = enumerateConnectedRegions(image, adjacency, 2);

    for (bool updateMaxTree : {true, false}) {
        for (const auto &region : regions) {
            for (int bValue = 0; bValue <= 7; ++bValue) {
                DynamicComponentTree minTree(image->clone(), false, adjacency);
                DynamicComponentTree maxTree(image->clone(), true, adjacency);
                DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;
                const auto contract =
                    problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
                        dualTree,
                        *adjacency,
                        std::span<const PixelId>(region),
                        static_cast<AltitudeType>(bValue),
                        true);
                if (!contract.acceptedByPublishedProblem || contract.partialOwners == 0) {
                    continue;
                }

                problem1_debug::ExplicitRegionProblem1UpdateOptions options;
                options.assumeConnected = true;
                options.allowLocalMaterialization = true;
                options.allowTargetNodeCreation = false;
                options.validateAgainstRebuild = true;

                const auto withoutTargetCreation =
                    problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
                        minTree,
                        maxTree,
                        adjacency,
                        updateMaxTree,
                        std::span<const PixelId>(region),
                        static_cast<AltitudeType>(bValue),
                        options);

                if (withoutTargetCreation.status ==
                        problem1_debug::ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation &&
                    withoutTargetCreation.materializationApplied &&
                    withoutTargetCreation.hasPostMaterializationContract &&
                    withoutTargetCreation.postMaterializationContract.partialOwners == 0 &&
                    !withoutTargetCreation.postMaterializationContract.hasTargetLevelUnionCandidate) {
                    return TargetNodeCreationEvent{
                        region,
                        static_cast<AltitudeType>(bValue),
                        updateMaxTree,
                        true};
                }
            }
        }
    }

    return TargetNodeCreationEvent{};
}

MaterializationEvent findMultiplePartialOwnerMaterializationEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    for (NodeId firstOwner : minTree.getNodeSubtree(minTree.getRoot())) {
        if (minTree.getAltitude(firstOwner) <= 0 ||
            minTree.getNumProperParts(firstOwner) < 2) {
            continue;
        }

        const auto firstOwnerPixels = collectProperParts(minTree, firstOwner);
        for (PixelId firstPixel : firstOwnerPixels) {
            for (PixelId neighbor : adjacency->getNeighborPixels(firstPixel)) {
                const NodeId secondOwner = minTree.getSmallestComponent(neighbor);
                if (secondOwner == InvalidNode ||
                    secondOwner == firstOwner ||
                    !minTree.isAlive(secondOwner) ||
                    minTree.getAltitude(secondOwner) <= 0 ||
                    minTree.getNumProperParts(secondOwner) < 2) {
                    continue;
                }

                std::vector<PixelId> pixelsC = {firstPixel, neighbor};
                std::sort(pixelsC.begin(), pixelsC.end());
                pixelsC.erase(std::unique(pixelsC.begin(), pixelsC.end()), pixelsC.end());
                if (pixelsC.size() != 2 ||
                    !problem1_debug::isConnectedRegion(*adjacency, std::span<const PixelId>(pixelsC))) {
                    continue;
                }

                const auto plan = problem1_debug::analyzeRegionMaterialization(
                    minTree,
                    *adjacency,
                    std::span<const PixelId>(pixelsC));
                if (plan.numPartialOwners() < 2) {
                    continue;
                }

                int minAltitude = 255;
                for (const auto &owner : plan.affectedOwners) {
                    minAltitude = std::min(minAltitude, owner.altitude);
                }
                if (minAltitude <= 0) {
                    continue;
                }

                return MaterializationEvent{
                    std::move(pixelsC),
                    static_cast<AltitudeType>(minAltitude - 1),
                    true};
            }
        }
    }

    return MaterializationEvent{};
}

MaterializationEvent findMultipleInsideComponentsMaterializationEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    for (NodeId partialOwner : minTree.getNodeSubtree(minTree.getRoot())) {
        if (minTree.getAltitude(partialOwner) <= 0 ||
            minTree.getNumProperParts(partialOwner) < 3) {
            continue;
        }

        const auto ownerPixels = collectProperParts(minTree, partialOwner);
        for (size_t i = 0; i < ownerPixels.size(); ++i) {
            for (size_t j = i + 1; j < ownerPixels.size(); ++j) {
                const std::vector<PixelId> selectedPartialPixels = {
                    ownerPixels[i],
                    ownerPixels[j],
                };
                const auto selectedComponents =
                    problem1_debug::collectLocalConnectedComponents(
                        *adjacency,
                        std::span<const PixelId>(selectedPartialPixels),
                        minTree.getNumTotalProperParts());
                if (selectedComponents.size() < 2) {
                    continue;
                }

                for (PixelId selectedPixel : selectedPartialPixels) {
                    for (PixelId neighbor : adjacency->getNeighborPixels(selectedPixel)) {
                        const NodeId fullOwner = minTree.getSmallestComponent(neighbor);
                        if (fullOwner == InvalidNode ||
                            fullOwner == partialOwner ||
                            !minTree.isAlive(fullOwner) ||
                            minTree.getAltitude(fullOwner) <= 0 ||
                            minTree.getNumProperParts(fullOwner) == 0) {
                            continue;
                        }

                        std::vector<PixelId> pixelsC = selectedPartialPixels;
                        const auto fullOwnerPixels = collectProperParts(minTree, fullOwner);
                        pixelsC.insert(
                            pixelsC.end(),
                            fullOwnerPixels.begin(),
                            fullOwnerPixels.end());
                        std::sort(pixelsC.begin(), pixelsC.end());
                        pixelsC.erase(std::unique(pixelsC.begin(), pixelsC.end()), pixelsC.end());

                        if (!problem1_debug::isConnectedRegion(
                                *adjacency,
                                std::span<const PixelId>(pixelsC))) {
                            continue;
                        }

                        const auto plan = problem1_debug::analyzeRegionMaterialization(
                            minTree,
                            *adjacency,
                            std::span<const PixelId>(pixelsC));
                        bool hasTargetOwner = false;
                        int minAltitude = 255;
                        for (const auto &owner : plan.affectedOwners) {
                            minAltitude = std::min(minAltitude, owner.altitude);
                            if (owner.nodeId == partialOwner &&
                                owner.cutsDirectProperParts() &&
                                owner.insideComponents.size() > 1) {
                                hasTargetOwner = true;
                            }
                        }
                        if (!hasTargetOwner ||
                            plan.numFullyCoveredOwners() < 1 ||
                            minAltitude <= 0) {
                            continue;
                        }

                        return MaterializationEvent{
                            std::move(pixelsC),
                            static_cast<AltitudeType>(minAltitude - 1),
                            true};
                    }
                }
            }
        }
    }

    return MaterializationEvent{};
}

TargetNodeCreationEvent findCompatibleMaterializedNonSubtreeEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    const auto regions = enumerateConnectedRegions(image, adjacency, 3);

    for (bool updateMaxTree : {true, false}) {
        DynamicComponentTree minTree(image->clone(), false, adjacency);
        DynamicComponentTree maxTree(image->clone(), true, adjacency);
        DynamicComponentTree &primalTree = updateMaxTree ? minTree : maxTree;
        DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;

        for (const auto &region : regions) {
            if (regionMatchesAnyNonRootSubtreeFootprint(
                    primalTree,
                    std::span<const PixelId>(region))) {
                continue;
            }

            for (int bValue = 0; bValue <= 9; ++bValue) {
                const auto b = static_cast<AltitudeType>(bValue);
                const auto contract =
                    problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
                        dualTree,
                        *adjacency,
                        std::span<const PixelId>(region),
                        b,
                        true);
                if (!contract.acceptedByPublishedProblem ||
                    !contract.compatibleWithCurrentUpdater ||
                    contract.partialOwners != 0 ||
                    !contract.hasTargetLevelUnionCandidate) {
                    continue;
                }

                problem1_debug::ExplicitRegionProblem1UpdateOptions options;
                options.assumeConnected = true;
                const auto preflight =
                    problem1_debug::preflightExplicitRegionProblem1Update<AltitudeType>(
                        minTree,
                        maxTree,
                        adjacency,
                        updateMaxTree,
                        std::span<const PixelId>(region),
                        b,
                        options);
                if (preflight.readyForLocalUpdate()) {
                    return TargetNodeCreationEvent{region, b, updateMaxTree, true};
                }
            }
        }
    }

    return TargetNodeCreationEvent{};
}

TargetNodeCreationEvent findMultiExtremalSingleWorkspaceEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    const auto regions = enumerateConnectedRegions(image, adjacency, 4);

    for (bool updateMaxTree : {true, false}) {
        for (const auto &region : regions) {
            for (int bValue = 0; bValue <= 15; ++bValue) {
                const auto b = static_cast<AltitudeType>(bValue);
                const auto analysis =
                    problem1_debug::analyzeExplicitRegionProblem1MultiAnchor<AltitudeType>(
                        image,
                        adjacency,
                        updateMaxTree,
                        std::span<const PixelId>(region),
                        b,
                        true,
                        true);
                if (!analysis.acceptedByPublishedProblem ||
                    !analysis.analyzed ||
                    !analysis.hasMultipleExtremalOwners() ||
                    analysis.generalWorkspaceAnchorCount != 1) {
                    continue;
                }

                const auto validation =
                    problem1_debug::validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
                        analysis,
                        true);
                if (validation.satisfiesIntegratedReleaseCriteria()) {
                    return TargetNodeCreationEvent{
                        region,
                        b,
                        updateMaxTree,
                        true};
                }
            }
        }
    }

    return TargetNodeCreationEvent{};
}

TargetNodeCreationEvent findClosureRootMismatchEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    const auto regions = enumerateConnectedRegions(image, adjacency, 4);

    for (bool updateMaxTree : {true, false}) {
        for (const auto &region : regions) {
            for (int bValue = 0; bValue <= 15; ++bValue) {
                const auto b = static_cast<AltitudeType>(bValue);
                const auto analysis =
                    problem1_debug::analyzeExplicitRegionProblem1MultiAnchor<AltitudeType>(
                        image,
                        adjacency,
                        updateMaxTree,
                        std::span<const PixelId>(region),
                        b,
                        true,
                        true);
                if (!analysis.acceptedByPublishedProblem ||
                    !analysis.analyzed ||
                    !analysis.closureRootDiffersFromActiveNodeCa() ||
                    !analysis.nodeCaToClosureRootSameLevelPath) {
                    continue;
                }

                const auto validation =
                    problem1_debug::validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
                        analysis,
                        true);
                if (validation.satisfiesIntegratedReleaseCriteria()) {
                    return TargetNodeCreationEvent{
                        region,
                        b,
                        updateMaxTree,
                        true};
                }
            }
        }
    }

    return TargetNodeCreationEvent{};
}

bool containsFailedHypothesis(
    const problem1_debug::ExplicitRegionProblem1LocalCorrectnessValidation<AltitudeType> &validation,
    const std::string &hypothesis) {
    return std::find(
        validation.failedHypotheses.begin(),
        validation.failedHypotheses.end(),
        hypothesis) != validation.failedHypotheses.end();
}

void testUpdaterRejectsInvalidProblem1WithoutMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 1};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.validateAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RejectedByContract,
            "non-extensive max-tree event must be rejected by the Problem 1 updater");
    require(!result.appliedLocalUpdate(), "rejected event must not be applied");
    require(maxTree.reconstructionImage()->isEqual(image),
            "rejected update must not mutate the max-tree reconstruction");
}

void testSafeUpdaterRejectsDisconnectedRegionBeforeMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 15};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.fallbackPolicy = problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(7),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RejectedByContract,
            "safe updater should reject disconnected C at the Problem 1 contract boundary");
    require(result.message.find("connected region C") != std::string::npos,
            "safe updater should report the connectivity contract violation");
    require(!result.fallbackAttempted,
            "safe updater must not fallback for events outside the published contract");
    require(!result.appliedAnyUpdate(),
            "disconnected C rejection must not apply any update");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "disconnected C rejection must leave both trees unchanged");
}

void testUpdaterAppliesCompatibleLowerEventAndMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        false,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "compatible lower event must be applied");
    require(result.appliedLocalUpdate(), "compatible event must run the local update");
    require(result.rebuildValidationPerformed, "test should enable rebuild validation");
    require(result.rebuildMatched, "local update must match rebuild oracle");
    require(result.structuralValidationPerformed,
            "test should enable structural validation");
    require(result.structuralMatchedRebuild,
            "local update tree must be structurally equivalent to rebuild tree");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated min-tree reconstruction must match expected image");
}

void testSafeUpdaterAppliesMaterializedNonSubtreeEventAndMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findCompatibleMaterializedNonSubtreeEvent(image, adjacency);
    require(event.found,
            "fixture must expose a materialized Problem 1 event that is not a primal subtree footprint");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree &primalTree = event.updateMaxTree ? minTree : maxTree;
    DynamicComponentTree &dualTree = event.updateMaxTree ? maxTree : minTree;

    require(!regionMatchesAnyNonRootSubtreeFootprint(
                primalTree,
                std::span<const PixelId>(event.pixelsC)),
            "selected explicit region C must not be representable as a Problem 2 subtree footprint");

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "safe updater should apply the materialized non-subtree Problem 1 event");
    require(result.contract.compatibleWithCurrentUpdater,
            "selected event should be compatible with the materialized/no-target-node subcase");
    require(result.contract.partialOwners == 0,
            "selected event should not require local materialization");
    require(!result.targetNodeCreationRequired,
            "selected event should not require target-node creation");
    require(result.localCorrectnessValidationPerformed,
            "safe updater should run non-mutating preflight before the local update");
    require(result.localCorrectnessValidationPassed,
            "selected event should pass preflight");
    require(result.rebuildMatched,
            "materialized non-subtree Problem 1 event must match rebuild");
    require(result.structuralMatchedRebuild,
            "materialized non-subtree Problem 1 event must match rebuild structurally");
    require(dualTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated dual-tree reconstruction must match expected image");
}

void testSafeUpdaterHandlesMultiExtremalSingleWorkspaceEvent() {
    auto image = make_prng_image(6, 6, 17u);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findMultiExtremalSingleWorkspaceEvent(image, adjacency);
    require(event.found,
            "fixture must expose a multi-extremal Problem 1 event with one general workspace root");

    const auto analysis =
        problem1_debug::analyzeExplicitRegionProblem1MultiAnchor<AltitudeType>(
            image,
            adjacency,
            event.updateMaxTree,
            std::span<const PixelId>(event.pixelsC),
            event.b,
            true,
            true);
    require(analysis.hasMultipleExtremalOwners(),
            "selected event should expose multiple extremal owners before update");
    require(analysis.generalWorkspaceAnchorCount == 1,
            "selected event should remain inside one general workspace root");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree &dualTree = event.updateMaxTree ? maxTree : minTree;

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "safe updater should apply multi-extremal events contained in one workspace");
    require(result.localCorrectnessValidationPerformed,
            "safe updater should validate the multi-extremal event before mutation");
    require(result.localCorrectnessValidationPassed,
            "multi-extremal single-workspace event should pass local correctness validation");
    require(result.localCorrectnessValidation.generalWorkspaceAnchorCount == 1,
            "safe updater should expose the single-workspace condition");
    require(result.rebuildMatched,
            "multi-extremal single-workspace event must match rebuild");
    require(result.structuralMatchedRebuild,
            "multi-extremal single-workspace event must match rebuild structurally");
    require(dualTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated dual-tree reconstruction must match expected image");
}

void testSafeUpdaterHandlesClosureRootDifferentFromNodeCa() {
    auto image = make_prng_image(6, 6, 17u);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findClosureRootMismatchEvent(image, adjacency);
    require(event.found,
            "fixture must expose a Problem 1 event where closureRoot differs from nodeCa");

    const auto analysis =
        problem1_debug::analyzeExplicitRegionProblem1MultiAnchor<AltitudeType>(
            image,
            adjacency,
            event.updateMaxTree,
            std::span<const PixelId>(event.pixelsC),
            event.b,
            true,
            true);
    require(analysis.closureRootDiffersFromActiveNodeCa(),
            "selected event should expose closureRoot != nodeCa before update");
    require(analysis.nodeCaToClosureRootSameLevelPath,
            "selected event should keep nodeCa and closureRoot on a same-level path");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree &dualTree = event.updateMaxTree ? maxTree : minTree;

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "safe updater should apply closure-root mismatch events that pass H6-H8");
    require(result.localCorrectnessValidationPerformed,
            "safe updater should validate closure-root mismatch events before mutation");
    require(result.localCorrectnessValidationPassed,
            "closure-root mismatch event should pass local correctness validation");
    require(result.localCorrectnessValidation.nodeCa != result.localCorrectnessValidation.closureRoot,
            "safe updater should expose closureRoot != nodeCa in the validation result");
    require(result.rebuildMatched,
            "closure-root mismatch event must match rebuild");
    require(result.structuralMatchedRebuild,
            "closure-root mismatch event must match rebuild structurally");
    require(dualTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated dual-tree reconstruction must match expected image");
}

void testLocalCorrectnessRejectsMultipleGeneralWorkspaceRoots() {
    problem1_debug::ExplicitRegionProblem1MultiAnchorAnalysis<AltitudeType> analysis;
    analysis.acceptedByPublishedProblem = true;
    analysis.analyzed = true;
    analysis.contract.acceptedByPublishedProblem = true;
    analysis.contract.compatibleWithCurrentUpdater = true;
    analysis.contract.partialOwners = 0;
    analysis.contract.hasTargetLevelUnionCandidate = true;
    analysis.contract.nodeCa = 3;
    analysis.closureRoot = 5;
    analysis.nodeCaToClosureRootPathLength = 1;
    analysis.nodeCaToClosureRootSameLevelPath = true;
    analysis.generalWorkspaceAnchorCount = 2;
    analysis.ownerIntervalBucketUnexplainedMissingNodeCount = 0;

    const auto validation =
        problem1_debug::validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
            analysis,
            true);

    require(!validation.satisfiesIntegratedReleaseCriteria(),
            "events with multiple general workspace roots must not be released");
    require(!validation.h6UniqueGeneralWorkspaceRoot,
            "multiple general workspace roots must fail H6");
    require(containsFailedHypothesis(validation, "H6"),
            "multiple general workspace roots must be reported as an H6 rejection");
    require(validation.h7SameLevelClosure,
            "synthetic H6 regression fixture should keep H7 satisfied");
    require(validation.h8BucketCompleteness,
            "synthetic H6 regression fixture should keep H8 satisfied");
}

void testLocalCorrectnessRejectsNonSameLevelClosurePath() {
    problem1_debug::ExplicitRegionProblem1MultiAnchorAnalysis<AltitudeType> analysis;
    analysis.acceptedByPublishedProblem = true;
    analysis.analyzed = true;
    analysis.contract.acceptedByPublishedProblem = true;
    analysis.contract.compatibleWithCurrentUpdater = true;
    analysis.contract.partialOwners = 0;
    analysis.contract.hasTargetLevelUnionCandidate = true;
    analysis.contract.nodeCa = 3;
    analysis.closureRoot = 5;
    analysis.nodeCaToClosureRootPathLength = 1;
    analysis.nodeCaToClosureRootSameLevelPath = false;
    analysis.generalWorkspaceAnchorCount = 1;
    analysis.ownerIntervalBucketUnexplainedMissingNodeCount = 0;

    const auto validation =
        problem1_debug::validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
            analysis,
            true);

    require(!validation.satisfiesIntegratedReleaseCriteria(),
            "events with a non-same-level closure path must not be released");
    require(validation.h6UniqueGeneralWorkspaceRoot,
            "synthetic H7 regression fixture should keep H6 satisfied");
    require(!validation.h7SameLevelClosure,
            "non-same-level closure paths must fail H7");
    require(containsFailedHypothesis(validation, "H7"),
            "non-same-level closure paths must be reported as an H7 rejection");
    require(validation.h8BucketCompleteness,
            "synthetic H7 regression fixture should keep H8 satisfied");
}

void testLocalCorrectnessRejectsUnexplainedBucketGap() {
    problem1_debug::ExplicitRegionProblem1MultiAnchorAnalysis<AltitudeType> analysis;
    analysis.acceptedByPublishedProblem = true;
    analysis.analyzed = true;
    analysis.contract.acceptedByPublishedProblem = true;
    analysis.contract.compatibleWithCurrentUpdater = true;
    analysis.contract.partialOwners = 0;
    analysis.contract.hasTargetLevelUnionCandidate = true;
    analysis.contract.nodeCa = 3;
    analysis.closureRoot = 5;
    analysis.nodeCaToClosureRootPathLength = 1;
    analysis.nodeCaToClosureRootSameLevelPath = true;
    analysis.generalWorkspaceAnchorCount = 1;
    analysis.ownerIntervalBucketMissingNodeCount = 1;
    analysis.ownerIntervalBucketUnexplainedMissingNodeCount = 1;

    const auto validation =
        problem1_debug::validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(
            analysis,
            true);

    require(!validation.satisfiesIntegratedReleaseCriteria(),
            "events with unexplained owner-interval bucket gaps must not be released");
    require(validation.h6UniqueGeneralWorkspaceRoot,
            "synthetic H8 regression fixture should keep H6 satisfied");
    require(validation.h7SameLevelClosure,
            "synthetic H8 regression fixture should keep H7 satisfied");
    require(!validation.h8BucketCompleteness,
            "unexplained owner-interval bucket gaps must fail H8");
    require(containsFailedHypothesis(validation, "H8"),
            "unexplained owner-interval bucket gaps must be reported as an H8 rejection");
}

void testUpdaterCreatesMissingTargetLevelNodeThenMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    for (bool updateMaxTree : {true, false}) {
        const auto event = findTargetNodeCreationEvent(image, adjacency, updateMaxTree);
        require(event.found,
                "fixture must expose a Problem 1 event that only lacks a level-b union node");

        DynamicComponentTree minTree(image->clone(), false, adjacency);
        DynamicComponentTree maxTree(image->clone(), true, adjacency);

        problem1_debug::ExplicitRegionProblem1UpdateOptions options;
        options.assumeConnected = true;
        options.allowTargetNodeCreation = true;
        options.validateAgainstRebuild = true;

        const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            updateMaxTree,
            std::span<const PixelId>(event.pixelsC),
            event.b,
            options);

        require(result.targetNodeCreationRequired,
                "selected event should require a new level-b union node");
        require(result.targetNodeCreationAllowed,
                "target node creation should be enabled for this test");
        require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
                "target-node-creation event must be applied");
        require(result.rebuildMatched,
                "target-node-creation event must match rebuild oracle");
        DynamicComponentTree &updatedTree = updateMaxTree ? maxTree : minTree;
        require(updatedTree.reconstructionImage()->isEqual(result.expectedImage),
                "updated dual-tree reconstruction must match expected image");
    }
}

void testUpdaterLocalCorrectnessGateRejectsMissingTargetPolicy() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findTargetNodeCreationEvent(image, adjacency, true);
    require(event.found,
            "fixture must expose a Problem 1 event that requires target-node creation");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = false;
    options.enforceLocalCorrectnessValidation = true;
    options.validateAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation,
            "updater should reject through the local-correctness gate when target creation is disabled");
    require(result.localCorrectnessValidationPerformed,
            "updater should record local-correctness validation");
    require(!result.localCorrectnessValidationPassed,
            "disabled target creation should fail local-correctness validation for this policy test");
    require(result.rejectedByLocalCorrectnessValidation,
            "updater should report gate rejection");
    require(!result.localCorrectnessValidation.targetLevelPresentOrCreated,
            "validation should expose the failed target-level gate");
    require(!result.appliedLocalUpdate(),
            "updater must not mutate through local update after gate rejection");
    DynamicComponentTree &dualTree = event.updateMaxTree ? maxTree : minTree;
    require(dualTree.reconstructionImage()->isEqual(image),
            "gate rejection must leave the target tree reconstruction unchanged");
}

void testPreflightDoesNotMutateMaterializationCandidate() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2, 0);
    require(owner != InvalidNode,
            "test image must have a min-tree partial-owner candidate");

    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};
    const auto b = static_cast<AltitudeType>(minTree.getAltitude(owner) - 1);
    const int numNodesBefore = minTree.getNumNodes();
    const int numSlotsBefore = minTree.getNumInternalNodeSlots();

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;

    const auto preflight =
        problem1_debug::preflightExplicitRegionProblem1Update<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            std::span<const PixelId>(pixelsC),
            b,
            options);

    require(preflight.readyForLocalUpdate(),
            "preflight should release the materialization candidate");
    require(preflight.localCorrectnessValidationPerformed,
            "preflight should run the integrated local-correctness validation");
    require(preflight.localCorrectnessValidationPassed,
            "preflight should pass the integrated local-correctness validation");
    require(preflight.hasPostMaterializationContract,
            "preflight should expose the post-materialization contract computed on a clone");
    require(minTree.getNumNodes() == numNodesBefore,
            "preflight must not create or release live nodes in the target tree");
    require(minTree.getNumInternalNodeSlots() == numSlotsBefore,
            "preflight must not allocate node slots in the target tree");
    require(minTree.reconstructionImage()->isEqual(image),
            "preflight must not alter the target tree reconstruction");
}

void testSafeUpdaterRejectsMissingTargetBeforeMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findTargetNodeCreationEvent(image, adjacency, true);
    require(event.found,
            "fixture must expose a Problem 1 event that requires target-node creation");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = false;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::KnownUpdaterLimitation,
            "safe updater should reject before mutation when target creation is disabled");
    require(result.localCorrectnessValidationPerformed,
            "safe updater should always run the local-correctness preflight");
    require(result.rejectedByLocalCorrectnessValidation,
            "safe updater should expose the local-correctness rejection");
    require(!result.localCorrectnessValidation.targetLevelPresentOrCreated,
            "safe updater should report the failed H5 target-level gate");
    require(!result.appliedLocalUpdate(),
            "safe updater must not run the local update after preflight rejection");
    DynamicComponentTree &dualTree = event.updateMaxTree ? maxTree : minTree;
    require(dualTree.reconstructionImage()->isEqual(image),
            "safe updater rejection must leave the target tree reconstruction unchanged");
}

void testSafeUpdaterFallsBackToRebuildForPreflightRejection() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findTargetNodeCreationEvent(image, adjacency, true);
    require(event.found,
            "fixture must expose a Problem 1 event that requires target-node creation");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = false;
    options.fallbackPolicy = problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RebuildFallbackApplied,
            "safe updater should apply rebuild fallback after preflight rejection");
    require(result.fallbackAttempted,
            "safe updater should report fallback attempt");
    require(result.fallbackApplied,
            "safe updater should report fallback application");
    require(result.rejectedByLocalCorrectnessValidation,
            "fallback should preserve the original preflight rejection reason");
    require(!result.appliedLocalUpdate(),
            "fallback path must not report a local update");
    require(result.appliedRebuildFallback(),
            "fallback path should expose rebuild application");
    require(result.appliedAnyUpdate(),
            "fallback path should count as an applied update");
    require(result.fallbackReconstruction &&
                minTree.reconstructionImage()->isEqual(result.fallbackReconstruction) &&
                maxTree.reconstructionImage()->isEqual(result.fallbackReconstruction),
            "rebuild fallback should rebuild both trees to the same updated image");
    require(!result.fallbackReconstruction->isEqual(image),
            "rebuild fallback should apply the requested level change");
}

void testSafeUpdaterDoesNotFallbackOutsidePublishedContract() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 1};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.fallbackPolicy = problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RejectedByContract,
            "fallback should not apply to events outside the published Problem 1 contract");
    require(!result.fallbackAttempted,
            "fallback should not be attempted for contract rejection");
    require(!result.appliedAnyUpdate(),
            "contract rejection must not apply any update");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "contract rejection with fallback enabled must leave both trees unchanged");
}

void testSafeUpdaterAppliesReleasedEventAfterPreflight() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::safeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        false,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "safe updater should apply an event released by preflight");
    require(result.localCorrectnessValidationPerformed,
            "safe updater should carry the preflight validation into the update result");
    require(result.localCorrectnessValidationPassed,
            "safe updater should expose the passing local-correctness validation");
    require(!result.rejectedByLocalCorrectnessValidation,
            "released safe update must not be marked as rejected by the gate");
    require(result.rebuildMatched,
            "safe update should match rebuild when rebuild validation is enabled");
    require(result.structuralMatchedRebuild,
            "safe update should match the rebuild structure when structural validation is enabled");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "safe update reconstruction should match expected image");
}

void testValidatedSafeFacadeOptionsEnableFullSafeBoundary() {
    const auto options =
        problem1_debug::makeValidatedSafeExplicitRegionProblem1Options();
    require(!options.assumeConnected,
            "validated safe defaults should check connectedness instead of assuming it");
    require(options.allowLocalMaterialization,
            "validated safe defaults should allow local materialization");
    require(options.allowTargetNodeCreation,
            "validated safe defaults should allow target-level node creation");
    require(options.fallbackPolicy ==
                problem1_debug::ExplicitRegionProblem1FallbackPolicy::NoFallback,
            "validated safe defaults should not fallback unless explicitly requested");
    require(options.validateAgainstRebuild,
            "validated safe defaults should compare reconstruction against rebuild");
    require(options.validateStructureAgainstRebuild,
            "validated safe defaults should compare structure against rebuild");
    require(options.runtimePostConditionValidation,
            "validated safe defaults should keep runtime post-condition validation enabled");

    const auto fallbackOptions =
        problem1_debug::makeValidatedSafeExplicitRegionProblem1Options(
            problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair);
    require(fallbackOptions.fallbackPolicy ==
                problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair,
            "validated safe options should preserve the requested fallback policy");
}

void testValidatedSafeFacadeAppliesMaterializationAndTargetCreation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findMaterializationThenTargetNodeCreationEvent(image, adjacency);
    require(event.found,
            "fixture must expose a Problem 1 event requiring materialization and target-node creation");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const auto result =
        problem1_debug::validatedSafeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            event.updateMaxTree,
            std::span<const PixelId>(event.pixelsC),
            event.b);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "validated safe facade should apply events released by H1-H8");
    require(result.localCorrectnessValidationPerformed,
            "validated safe facade should run preflight validation before mutation");
    require(result.localCorrectnessValidationPassed,
            "validated safe facade should only mutate after H1-H8 release");
    require(result.materializationAttempted,
            "validated safe facade should attempt required local materialization");
    require(result.materializationApplied,
            "validated safe facade should materialize C before the local update");
    require(result.hasPostMaterializationContract,
            "validated safe facade should expose the post-materialization contract");
    require(result.postMaterializationContract.partialOwners == 0,
            "validated safe facade should remove partial owners before the local update");
    require(result.targetNodeCreationRequired,
            "validated safe facade should detect required level-b node creation");
    require(result.targetNodeCreationAllowed,
            "validated safe facade should allow level-b node creation");
    require(result.rebuildValidationPerformed,
            "validated safe facade should validate reconstruction against rebuild");
    require(result.rebuildMatched,
            "validated safe facade local reconstruction must match rebuild");
    require(result.structuralValidationPerformed,
            "validated safe facade should validate structure against rebuild");
    require(result.structuralMatchedRebuild,
            "validated safe facade local structure must match rebuild");

    DynamicComponentTree &updatedTree = event.updateMaxTree ? maxTree : minTree;
    require(updatedTree.reconstructionImage()->isEqual(result.expectedImage),
            "validated safe facade should leave the updated dual tree at the rebuild image");
    requireNoAliveEmptyNodes(updatedTree);
}

void testValidatedSafeFacadeRejectsOutsideContractWithoutMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> disconnectedC = {0, 15};
    const auto result =
        problem1_debug::validatedSafeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            true,
            std::span<const PixelId>(disconnectedC),
            static_cast<AltitudeType>(7),
            problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RejectedByContract,
            "validated safe facade should reject events outside the published Problem 1 contract");
    require(!result.fallbackAttempted,
            "validated safe facade must not fallback for events outside the published contract");
    require(!result.appliedAnyUpdate(),
            "validated safe facade rejection must not apply local or rebuild updates");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "validated safe facade contract rejection must leave both trees unchanged");
}

void testDescriptorPreflightAcceptsConsistentRegionLevelChangeWithoutMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    const auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));
    const int minNodesBefore = minTree.getNumNodes();
    const int maxNodesBefore = maxTree.getNumNodes();
    const auto options =
        problem1_debug::makeValidatedSafeExplicitRegionProblem1Options();

    const auto preflight =
        problem1_debug::preflightExplicitRegionProblem1Update<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored.change,
            options);

    require(preflight.readyForLocalUpdate(),
            "descriptor preflight should release a consistent RegionLevelChange");
    require(preflight.contract.a == anchored.change.a,
            "descriptor preflight should validate source level a");
    require(preflight.contract.b == anchored.change.b,
            "descriptor preflight should validate target level b");
    require(preflight.contract.direction == anchored.change.direction,
            "descriptor preflight should validate direction");
    require(preflight.localCorrectnessValidationPerformed,
            "descriptor preflight should run H1-H8 validation");
    require(preflight.localCorrectnessValidationPassed,
            "descriptor preflight should pass H1-H8 validation");
    require(minTree.getNumNodes() == minNodesBefore &&
                maxTree.getNumNodes() == maxNodesBefore,
            "descriptor preflight must not add or remove live nodes");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "descriptor preflight must not mutate either tree reconstruction");
}

void testDescriptorPreflightRejectsInconsistentRegionLevelChangeWithoutMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));
    anchored.change.direction = problem1_debug::RegionUpdateDirection::Raise;
    const int minNodesBefore = minTree.getNumNodes();
    const int maxNodesBefore = maxTree.getNumNodes();
    const auto options =
        problem1_debug::makeValidatedSafeExplicitRegionProblem1Options(
            problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair);

    const auto preflight =
        problem1_debug::preflightExplicitRegionProblem1Update<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored.change,
            options);

    require(preflight.status == problem1_debug::ExplicitRegionProblem1PreflightStatus::RejectedByContract,
            "descriptor preflight should reject inconsistent RegionLevelChange metadata");
    require(!preflight.readyForLocalUpdate(),
            "descriptor preflight rejection must not release the local update");
    require(preflight.message.find("direction") != std::string::npos,
            "descriptor preflight should report the inconsistent descriptor field");
    require(minTree.getNumNodes() == minNodesBefore &&
                maxTree.getNumNodes() == maxNodesBefore,
            "descriptor preflight rejection must not add or remove live nodes");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "descriptor preflight rejection must not mutate either tree reconstruction");
}

void testAnchoredDescriptorPreflightAcceptsConsistentAnchorWithoutMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    const auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));
    const int minNodesBefore = minTree.getNumNodes();
    const int maxNodesBefore = maxTree.getNumNodes();
    const auto options =
        problem1_debug::makeValidatedSafeExplicitRegionProblem1Options();

    const auto preflight =
        problem1_debug::preflightExplicitRegionProblem1Update<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored,
            options);

    require(preflight.readyForLocalUpdate(),
            "anchored descriptor preflight should release a consistent descriptor");
    require(preflight.contract.nodeCa == anchored.nodeCa,
            "anchored descriptor preflight should validate nodeCa");
    require(preflight.contract.a == anchored.change.a,
            "anchored descriptor preflight should validate source level a");
    require(preflight.contract.b == anchored.change.b,
            "anchored descriptor preflight should validate target level b");
    require(preflight.contract.direction == anchored.change.direction,
            "anchored descriptor preflight should validate direction");
    require(minTree.getNumNodes() == minNodesBefore &&
                maxTree.getNumNodes() == maxNodesBefore,
            "anchored descriptor preflight must not add or remove live nodes");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "anchored descriptor preflight must not mutate either tree reconstruction");
}

void testAnchoredDescriptorPreflightRejectsInconsistentAnchorWithoutMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));
    anchored.nodeCa = InvalidNode;
    const int minNodesBefore = minTree.getNumNodes();
    const int maxNodesBefore = maxTree.getNumNodes();
    const auto options =
        problem1_debug::makeValidatedSafeExplicitRegionProblem1Options(
            problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair);

    const auto preflight =
        problem1_debug::preflightExplicitRegionProblem1Update<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored,
            options);

    require(preflight.status == problem1_debug::ExplicitRegionProblem1PreflightStatus::RejectedByContract,
            "anchored descriptor preflight should reject inconsistent nodeCa");
    require(!preflight.readyForLocalUpdate(),
            "anchored descriptor preflight rejection must not release the local update");
    require(preflight.message.find("nodeCa") != std::string::npos,
            "anchored descriptor preflight should report the inconsistent anchor field");
    require(minTree.getNumNodes() == minNodesBefore &&
                maxTree.getNumNodes() == maxNodesBefore,
            "anchored descriptor preflight rejection must not add or remove live nodes");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "anchored descriptor preflight rejection must not mutate either tree reconstruction");
}

void testValidatedSafeFacadeAcceptsRegionLevelChangeDescriptor() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));

    const auto result =
        problem1_debug::validatedSafeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored.change);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "validated safe facade should accept a consistent RegionLevelChange descriptor");
    require(result.contract.a == anchored.change.a,
            "accepted descriptor should expose the validated source level a");
    require(result.contract.b == anchored.change.b,
            "accepted descriptor should expose the validated target level b");
    require(result.contract.direction == anchored.change.direction,
            "accepted descriptor should expose the validated update direction");
    require(result.localCorrectnessValidationPassed,
            "descriptor facade should keep the H1-H8 gate enabled");
    require(result.rebuildMatched,
            "descriptor facade should validate reconstruction against rebuild");
    require(result.structuralMatchedRebuild,
            "descriptor facade should validate structure against rebuild");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "descriptor facade should leave the updated tree at the rebuild image");
}

void testValidatedSafeFacadeRejectsInconsistentRegionLevelChangeDescriptor() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));
    anchored.change.a = static_cast<AltitudeType>(anchored.change.a + 1);

    const auto result =
        problem1_debug::validatedSafeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored.change,
            problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RejectedByContract,
            "descriptor facade should reject inconsistent RegionLevelChange metadata");
    require(result.message.find("level a") != std::string::npos,
            "descriptor facade should report the inconsistent descriptor field");
    require(!result.fallbackAttempted,
            "descriptor mismatch must not fall back to rebuild");
    require(!result.appliedAnyUpdate(),
            "descriptor mismatch must not apply local or rebuild updates");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "descriptor mismatch must leave both trees unchanged");
}

void testValidatedSafeFacadeAcceptsAnchoredRegionLevelChangeDescriptor() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    const auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));

    const auto result =
        problem1_debug::validatedSafeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "validated safe facade should accept a consistent AnchoredRegionLevelChange descriptor");
    require(result.contract.nodeCa == anchored.nodeCa,
            "accepted anchored descriptor should expose the validated nodeCa");
    require(result.contract.a == anchored.change.a,
            "accepted anchored descriptor should expose the validated source level a");
    require(result.contract.b == anchored.change.b,
            "accepted anchored descriptor should expose the validated target level b");
    require(result.contract.direction == anchored.change.direction,
            "accepted anchored descriptor should expose the validated update direction");
    require(result.localCorrectnessValidationPassed,
            "anchored descriptor facade should keep the H1-H8 gate enabled");
    require(result.rebuildMatched,
            "anchored descriptor facade should validate reconstruction against rebuild");
    require(result.structuralMatchedRebuild,
            "anchored descriptor facade should validate structure against rebuild");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "anchored descriptor facade should leave the updated tree at the rebuild image");
}

void testValidatedSafeFacadeRejectsInconsistentAnchoredDescriptor() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto anchored =
        problem1_debug::buildExplicitAnchoredRegionChange<AltitudeType>(
            minTree,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(4));
    anchored.nodeCa = InvalidNode;

    const auto result =
        problem1_debug::validatedSafeUpdateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
            minTree,
            maxTree,
            adjacency,
            false,
            anchored,
            problem1_debug::ExplicitRegionProblem1FallbackPolicy::RebuildPair);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::RejectedByContract,
            "anchored descriptor facade should reject inconsistent nodeCa");
    require(result.message.find("nodeCa") != std::string::npos,
            "anchored descriptor facade should report the inconsistent anchor field");
    require(!result.fallbackAttempted,
            "anchored descriptor mismatch must not fall back to rebuild");
    require(!result.appliedAnyUpdate(),
            "anchored descriptor mismatch must not apply local or rebuild updates");
    require(minTree.reconstructionImage()->isEqual(image) &&
                maxTree.reconstructionImage()->isEqual(image),
            "anchored descriptor mismatch must leave both trees unchanged");
}

void testUpdaterMaterializesPartialOwnerThenMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2, 0);
    require(owner != InvalidNode,
            "test image must have a nonzero min-tree node with at least two direct proper parts");

    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};
    const auto b = static_cast<AltitudeType>(minTree.getAltitude(owner) - 1);

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.validateAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        false,
        std::span<const PixelId>(pixelsC),
        b,
        options);

    require(result.materializationAttempted,
            "partial-owner event should attempt materialization");
    require(result.materializationApplied,
            "partial-owner event should be materialized");
    require(result.materializedOwners == 1,
            "fixture should materialize one owner");
    require(result.materializedNodes >= 1,
            "materialization should create at least one local node");
    require(result.hasPostMaterializationContract,
            "updater should expose post-materialization contract");
    require(result.postMaterializationContract.partialOwners == 0,
            "post-materialization contract should no longer have partial owners");
    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "materialize-then-update must be applied");
    require(result.rebuildMatched,
            "materialize-then-update must match rebuild oracle");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated min-tree reconstruction must match expected image");
}

void testUpdaterMaterializesMultiplePartialOwnersThenMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findMultiplePartialOwnerMaterializationEvent(image, adjacency);
    require(event.found,
            "fixture must expose a Problem 1 event with multiple partial owners");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const auto prePlan = problem1_debug::analyzeRegionMaterialization(
        minTree,
        *adjacency,
        std::span<const PixelId>(event.pixelsC));
    require(prePlan.numPartialOwners() >= 2,
            "selected event should start with multiple partial owners");

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        false,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.materializationAttempted,
            "multiple-partial-owner event should attempt materialization");
    require(result.materializationApplied,
            "multiple-partial-owner event should be materialized");
    require(result.materializedOwners >= 2,
            "multiple-partial-owner event should materialize at least two owners");
    require(result.hasPostMaterializationContract,
            "updater should expose post-materialization contract");
    require(result.postMaterializationContract.partialOwners == 0,
            "post-materialization contract should remove all partial owners");
    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "multiple-partial-owner materialize-then-update must be applied");
    require(result.rebuildMatched,
            "multiple-partial-owner materialize-then-update must match rebuild oracle");
    require(result.structuralMatchedRebuild,
            "multiple-partial-owner materialize-then-update must match rebuild structurally");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated min-tree reconstruction must match expected image");
    requireNoAliveEmptyNodes(minTree);
}

void testUpdaterMaterializesMultipleInsideComponentsThenMatchesRebuild() {
    auto image = make_complex_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findMultipleInsideComponentsMaterializationEvent(image, adjacency);
    require(event.found,
            "fixture must expose a partial owner with multiple inside components");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const auto prePlan = problem1_debug::analyzeRegionMaterialization(
        minTree,
        *adjacency,
        std::span<const PixelId>(event.pixelsC));
    bool hasMultipleInsideComponents = false;
    for (const auto &owner : prePlan.affectedOwners) {
        if (owner.cutsDirectProperParts() &&
            owner.insideComponents.size() > 1) {
            hasMultipleInsideComponents = true;
        }
    }
    require(hasMultipleInsideComponents,
            "selected event should force one partial owner to produce multiple inside components");

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        false,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.materializationAttempted,
            "multiple-inside-components event should attempt materialization");
    require(result.materializationApplied,
            "multiple-inside-components event should be materialized");
    require(result.materializedNodes >= 2,
            "multiple-inside-components event should create at least two local nodes");
    require(result.hasPostMaterializationContract,
            "updater should expose post-materialization contract");
    require(result.postMaterializationContract.partialOwners == 0,
            "post-materialization contract should remove all partial owners");
    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "multiple-inside-components materialize-then-update must be applied");
    require(result.rebuildMatched,
            "multiple-inside-components materialize-then-update must match rebuild oracle");
    require(result.structuralMatchedRebuild,
            "multiple-inside-components materialize-then-update must match rebuild structurally");
    require(minTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated min-tree reconstruction must match expected image");
    requireNoAliveEmptyNodes(minTree);
}

void testUpdaterMaterializesThenCreatesMissingTargetNodeAndMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const auto event = findMaterializationThenTargetNodeCreationEvent(image, adjacency);
    require(event.found,
            "fixture must expose a Problem 1 event requiring materialization and target-node creation");

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        options);

    require(result.materializationAttempted,
            "combined event should attempt materialization");
    require(result.materializationApplied,
            "combined event should be materialized before update");
    require(result.hasPostMaterializationContract,
            "combined event should expose post-materialization contract");
    require(result.postMaterializationContract.partialOwners == 0,
            "combined event should remove partial owners after materialization");
    require(result.targetNodeCreationRequired,
            "combined event should still require target-node creation after materialization");
    require(result.targetNodeCreationAllowed,
            "target-node creation should be enabled for the combined event");
    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "materialize-then-create-target-node event must be applied");
    require(result.rebuildMatched,
            "materialize-then-create-target-node event must match rebuild oracle");
    DynamicComponentTree &updatedTree = event.updateMaxTree ? maxTree : minTree;
    require(updatedTree.reconstructionImage()->isEqual(result.expectedImage),
            "updated dual-tree reconstruction must match expected image");
}

void testUpdaterReleasesEmptyDetachedNodesAfterExplicitMaterialization() {
    auto image = make_empty_detached_regression_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {11, 17, 23, 28, 29, 35};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.runtimePostConditionValidation = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(10),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "regression event should be applied");
    require(result.materializationApplied,
            "regression event should require materialization");
    require(result.targetNodeCreationAllowed,
            "regression event should require target-node creation");
    require(result.rebuildMatched,
            "regression event should match rebuild oracle");
    require(maxTree.reconstructionImage()->isEqual(result.expectedImage),
            "regression event reconstruction should match expected image");
    requireNoAliveEmptyNodes(maxTree);
}

void testUpdaterBuildsCanonicalTreeForSameLevelMaterializationRegression() {
    auto image = make_prng_image(8, 8, 23u);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0};
    problem1_debug::ExplicitRegionProblem1UpdateOptions options;
    options.assumeConnected = true;
    options.allowLocalMaterialization = true;
    options.allowTargetNodeCreation = true;
    options.validateAgainstRebuild = true;
    options.validateStructureAgainstRebuild = true;
    options.runtimePostConditionValidation = true;
    options.requireCanonicalByConstruction = true;

    const auto result = problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
        minTree,
        maxTree,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(10),
        options);

    require(result.status == problem1_debug::ExplicitRegionProblem1UpdateStatus::LocalUpdateApplied,
            "same-level edge regression event should be applied canonically");
    require(result.canonicalByConstructionPreflightPerformed,
            "same-level edge regression should run canonical-by-construction preflight");
    require(result.canonicalByConstructionAccepted,
            "same-level edge regression should be released by the canonical-by-construction preparation");
    require(result.materializationApplied,
            "same-level edge regression should use C-only local materialization");
    require(result.structuralValidationPerformed,
            "same-level edge regression should run structural validation");
    require(result.structuralMatchedRebuild,
            "same-level edge regression should be structurally equivalent to rebuild");
    require(result.sameLevelEdgesCollapsed == 0,
            "same-level edge regression must not collapse same-level edges after the update");
    require(result.postCanonicalizationChangesRequired == 0,
            "same-level edge regression must not require corrective post-canonicalization");
    require(result.postEmptyDetachedChildlessNodesRequired == 0,
            "same-level edge regression must not require releasing empty detached nodes");
    require(result.postSameLevelParentChildEdgesRequired == 0,
            "same-level edge regression must not leave same-level parent-child edges");
    require(maxTree.reconstructionImage()->isEqual(result.expectedImage),
            "same-level edge regression reconstruction should match expected image");
    requireNoAliveEmptyNodes(maxTree);
}

} // namespace

int main() {
    testUpdaterRejectsInvalidProblem1WithoutMutation();
    testSafeUpdaterRejectsDisconnectedRegionBeforeMutation();
    testUpdaterAppliesCompatibleLowerEventAndMatchesRebuild();
    testSafeUpdaterAppliesMaterializedNonSubtreeEventAndMatchesRebuild();
    testSafeUpdaterHandlesMultiExtremalSingleWorkspaceEvent();
    testSafeUpdaterHandlesClosureRootDifferentFromNodeCa();
    testLocalCorrectnessRejectsMultipleGeneralWorkspaceRoots();
    testLocalCorrectnessRejectsNonSameLevelClosurePath();
    testLocalCorrectnessRejectsUnexplainedBucketGap();
    testUpdaterCreatesMissingTargetLevelNodeThenMatchesRebuild();
    testUpdaterLocalCorrectnessGateRejectsMissingTargetPolicy();
    testPreflightDoesNotMutateMaterializationCandidate();
    testSafeUpdaterRejectsMissingTargetBeforeMutation();
    testSafeUpdaterFallsBackToRebuildForPreflightRejection();
    testSafeUpdaterDoesNotFallbackOutsidePublishedContract();
    testSafeUpdaterAppliesReleasedEventAfterPreflight();
    testValidatedSafeFacadeOptionsEnableFullSafeBoundary();
    testValidatedSafeFacadeAppliesMaterializationAndTargetCreation();
    testValidatedSafeFacadeRejectsOutsideContractWithoutMutation();
    testDescriptorPreflightAcceptsConsistentRegionLevelChangeWithoutMutation();
    testDescriptorPreflightRejectsInconsistentRegionLevelChangeWithoutMutation();
    testAnchoredDescriptorPreflightAcceptsConsistentAnchorWithoutMutation();
    testAnchoredDescriptorPreflightRejectsInconsistentAnchorWithoutMutation();
    testValidatedSafeFacadeAcceptsRegionLevelChangeDescriptor();
    testValidatedSafeFacadeRejectsInconsistentRegionLevelChangeDescriptor();
    testValidatedSafeFacadeAcceptsAnchoredRegionLevelChangeDescriptor();
    testValidatedSafeFacadeRejectsInconsistentAnchoredDescriptor();
    testUpdaterMaterializesPartialOwnerThenMatchesRebuild();
    testUpdaterMaterializesMultiplePartialOwnersThenMatchesRebuild();
    testUpdaterMaterializesMultipleInsideComponentsThenMatchesRebuild();
    testUpdaterMaterializesThenCreatesMissingTargetNodeAndMatchesRebuild();
    testUpdaterReleasesEmptyDetachedNodesAfterExplicitMaterialization();
    testUpdaterBuildsCanonicalTreeForSameLevelMaterializationRegression();
    return 0;
}
