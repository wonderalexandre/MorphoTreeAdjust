#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DualMinMaxTreeIncrementalFilter.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Problem1.hpp"

namespace problem1_debug = debug::problem1;

namespace {

using namespace debug::problem1;

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

ImageUInt8Ptr make_large_demo_image() {
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

struct DebugCheckedUpdaterEvent {
    std::vector<PixelId> pixelsC;
    AltitudeType b = 0;
    bool updateMaxTree = false;
    Problem1PreflightStatus status = Problem1PreflightStatus::EmptyRegion;
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
                                                            AdjacencyRelation &adjacency,
                                                            int maxRegionSize) {
    std::vector<std::vector<PixelId>> regions;
    const int domainSize = image->getSize();
    for (int regionSize = 1; regionSize <= maxRegionSize; ++regionSize) {
        std::vector<std::vector<PixelId>> candidates;
        std::vector<PixelId> current;
        current.reserve(static_cast<size_t>(regionSize));
        enumeratePixelCombinations(domainSize, regionSize, 0, current, candidates);

        for (auto &candidate : candidates) {
            if (isConnectedExplicitRegionProblem1Region(
                    adjacency,
                    std::span<const PixelId>(candidate))) {
                regions.push_back(std::move(candidate));
            }
        }
    }
    return regions;
}

DebugCheckedUpdaterEvent findDebugCheckedUpdaterReadyEvent(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency) {
    const auto regions = enumerateConnectedRegions(image, *adjacency, 3);
    DebugCheckedUpdaterEvent acceptedFallback;

    for (bool updateMaxTree : {true, false}) {
        DynamicComponentTree minTree(image->clone(), false, adjacency);
        DynamicComponentTree maxTree(image->clone(), true, adjacency);
        DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;

        for (const auto &region : regions) {
            for (int bValue = 0; bValue <= 9; ++bValue) {
                const auto b = static_cast<AltitudeType>(bValue);
                const auto preflight = checkExplicitRegionProblem1Preflight<AltitudeType>(
                    dualTree,
                    *adjacency,
                    std::span<const PixelId>(region),
                    b,
                    true);

                if (preflight.status == Problem1PreflightStatus::RequiresTargetNodeCreation) {
                    return DebugCheckedUpdaterEvent{
                        region,
                        b,
                        updateMaxTree,
                        preflight.status,
                        true};
                }

                if (!acceptedFallback.found &&
                    preflight.status == Problem1PreflightStatus::Accepted) {
                    acceptedFallback = DebugCheckedUpdaterEvent{
                        region,
                        b,
                        updateMaxTree,
                        preflight.status,
                        true};
                }
            }
        }
    }

    return acceptedFallback;
}

void requireLightPreflightMatchesContract(const DynamicComponentTree &tree,
                                          AdjacencyRelation &adjacency,
                                          const std::vector<PixelId> &pixelsC,
                                          AltitudeType b,
                                          const std::string &message) {
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        tree,
        adjacency,
        std::span<const PixelId>(pixelsC),
        b);
    const auto preflight = checkExplicitRegionProblem1Preflight<AltitudeType>(
        tree,
        adjacency,
        std::span<const PixelId>(pixelsC),
        b);

    require(preflight.status == contract.preflight.status, message);
    if (contract.acceptsPublishedProblem()) {
        require(preflight.a == contract.a,
                "light preflight must match contract level a");
        require(preflight.b == contract.b,
                "light preflight must match contract level b");
        require(preflight.direction == contract.direction,
                "light preflight must match contract direction");
        require(preflight.nodeCa == contract.nodeCa,
                "light preflight must match contract nodeCa");
    }
}

void testAcceptsPublishedProblem1MaxTreeRaise() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 1};
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        maxTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(7));

    require(contract.acceptedByPublishedProblem,
            "max-tree raise event must satisfy published Problem 1");
    require(contract.preflight.acceptedByPublishedProblem(),
            "typed preflight must accept the published Problem 1 event");
    require(contract.preflight.a == contract.a,
            "legacy a mirror must match the typed preflight result");
    require(contract.preflight.direction == contract.direction,
            "legacy direction mirror must match the typed preflight result");
    const auto preflight = checkExplicitRegionProblem1Preflight<AltitudeType>(
        maxTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(7));
    require(preflight.status == contract.preflight.status,
            "typed preflight helper must match the diagnostic contract status");
    require(preflight.nodeCa == contract.nodeCa,
            "typed preflight helper must return the same nodeCa");
    requireLightPreflightMatchesContract(
        maxTree,
        *adjacency,
        pixelsC,
        static_cast<AltitudeType>(7),
        "light preflight must match the full contract for max-tree raise");
    require(contract.isMaxtree, "contract must identify max-tree orientation");
    require(contract.direction == problem1_debug::RegionUpdateDirection::Raise,
            "max-tree Problem 1 event must be a raise");
    require(contract.a == 4, "a must be min I-(C)");
    require(contract.maxC == 4, "max I-(C) must be tracked");
    require(contract.regionSize == 2, "region size must be tracked");
    require(contract.affectedIntervalLevels == 3,
            "max-tree interval should track levels a+1 through b");
    require(contract.missingLevelUnionCandidates >= 0 &&
                contract.missingLevelUnionCandidates <= contract.affectedIntervalLevels,
            "missing candidate count must stay inside the affected interval");
    require(contract.missingIntermediateLevelUnionCandidates >= 0 &&
                contract.missingIntermediateLevelUnionCandidates <= 2,
            "missing intermediate candidate count must exclude the target level");
}

void testRejectsNonExtensiveMaxTreeUpdate() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 1};
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        maxTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4));

    require(!contract.acceptedByPublishedProblem,
            "max-tree event with b <= max I-(C) must be rejected");
    require(contract.preflight.status == problem1_debug::Problem1PreflightStatus::InvalidTargetLevel,
            "typed preflight must report an invalid target level");
    requireLightPreflightMatchesContract(
        maxTree,
        *adjacency,
        pixelsC,
        static_cast<AltitudeType>(4),
        "light preflight must match the full contract for invalid target level");
}

void testRejectsDisconnectedRegion() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 15};
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        maxTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(7));

    require(!contract.acceptedByPublishedProblem,
            "disconnected region C must be rejected");
    require(contract.preflight.status == problem1_debug::Problem1PreflightStatus::DisconnectedRegion,
            "typed preflight must report a disconnected region");
    requireLightPreflightMatchesContract(
        maxTree,
        *adjacency,
        pixelsC,
        static_cast<AltitudeType>(7),
        "light preflight must match the full contract for disconnected C");
}

void testRejectsDuplicatedPixels() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {0, 0};
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        maxTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(7));

    require(!contract.acceptedByPublishedProblem,
            "duplicated pixels in C must be rejected");
    require(contract.preflight.status == problem1_debug::Problem1PreflightStatus::DuplicatePixel,
            "typed preflight must report duplicated pixels");
    requireLightPreflightMatchesContract(
        maxTree,
        *adjacency,
        pixelsC,
        static_cast<AltitudeType>(7),
        "light preflight must match the full contract for duplicated pixels");
}

void testAcceptsDualMinTreeLower() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const std::vector<PixelId> pixelsC = {13, 14};
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        minTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4));

    require(contract.acceptedByPublishedProblem,
            "min-tree dual lower event must satisfy the dual contract");
    require(contract.preflight.acceptedByPublishedProblem(),
            "typed preflight must accept the published min-tree dual event");
    require(contract.preflight.a == contract.a,
            "legacy dual a mirror must match the typed preflight result");
    requireLightPreflightMatchesContract(
        minTree,
        *adjacency,
        pixelsC,
        static_cast<AltitudeType>(4),
        "light preflight must match the full contract for min-tree dual lower");
    require(!contract.isMaxtree, "contract must identify min-tree orientation");
    require(contract.direction == problem1_debug::RegionUpdateDirection::Lower,
            "min-tree dual event must be a lower update");
    require(contract.a == 6, "dual a must be max I-(C)");
    require(contract.minC == 6, "min I-(C) must be tracked");
    require(contract.affectedIntervalLevels == 2,
            "min-tree interval should track levels a-1 down to b");
    require(contract.missingLevelUnionCandidates >= 0 &&
                contract.missingLevelUnionCandidates <= contract.affectedIntervalLevels,
            "dual missing candidate count must stay inside the affected interval");
    require(contract.missingIntermediateLevelUnionCandidates >= 0 &&
                contract.missingIntermediateLevelUnionCandidates <= 1,
            "dual missing intermediate candidate count must exclude the target level");
}

void testAcceptsPublishedProblem1WithMultipleOwners() {
    auto image = ImageUInt8::create(2, 2);
    const std::vector<uint8_t> raw = {
        12, 10,
        11, 10,
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }

    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);

    const std::vector<PixelId> pixelsC = {1, 2, 3};
    const auto contract = problem1_debug::checkExplicitRegionProblem1Contract<AltitudeType>(
        maxTree,
        *adjacency,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(12));

    require(contract.acceptedByPublishedProblem,
            "max-tree explicit region with multiple owners must satisfy published Problem 1 contract");
    require(contract.affectedOwners > 1,
            "published contract should track more than one affected owner");
    requireLightPreflightMatchesContract(
        maxTree,
        *adjacency,
        pixelsC,
        static_cast<AltitudeType>(12),
        "light preflight must match the full contract for multiple owners");
}

void testDebugCheckedUpdaterRejectsBeforeMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adjacency);

    const auto before = maxTree.reconstructionImage();
    const std::vector<PixelId> pixelsC = {0, 15};
    bool rejected = false;
    try {
        updateTreeWithExplicitRegionChecked<AltitudeType>(
            adjust,
            &maxTree,
            *adjacency,
            std::span<const PixelId>(pixelsC),
            static_cast<AltitudeType>(7),
            true);
    } catch (const std::runtime_error &) {
        rejected = true;
    }

    require(rejected,
            "debug checked updater must reject invalid preflight cases before mutation");
    require(maxTree.reconstructionImage()->isEqual(before),
            "debug checked updater rejection must leave the tree reconstruction unchanged");
}

void testDebugCheckedUpdaterAppliesWhenTargetCreationIsAllowed() {
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    const auto event = findDebugCheckedUpdaterReadyEvent(image, adjacency);
    require(event.found,
            "fixture search must find a debug checked updater-ready Problem 1 event");

    DynamicComponentTree &dualTree = event.updateMaxTree ? maxTree : minTree;
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adjacency);

    const auto preflight = checkExplicitRegionProblem1Preflight<AltitudeType>(
        dualTree,
        *adjacency,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        true);
    require(preflight.status == event.status,
            "fixture search must preserve the debug checked updater preflight status");

    if (preflight.requiresTargetNodeCreation()) {
        bool rejected = false;
        try {
            updateTreeWithExplicitRegionChecked<AltitudeType>(
                adjust,
                &dualTree,
                *adjacency,
                std::span<const PixelId>(event.pixelsC),
                event.b,
                false);
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        require(rejected,
                "debug checked updater must reject missing target creation unless it is allowed");
    }

    updateTreeWithExplicitRegionChecked<AltitudeType>(
        adjust,
        &dualTree,
        *adjacency,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        true);
}

} // namespace

int main() {
    testAcceptsPublishedProblem1MaxTreeRaise();
    testRejectsNonExtensiveMaxTreeUpdate();
    testRejectsDisconnectedRegion();
    testRejectsDuplicatedPixels();
    testAcceptsDualMinTreeLower();
    testAcceptsPublishedProblem1WithMultipleOwners();
    testDebugCheckedUpdaterRejectsBeforeMutation();
    testDebugCheckedUpdaterAppliesWhenTargetCreationIsAllowed();
    return 0;
}
