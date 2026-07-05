#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1Oracle.hpp"

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

struct TargetNodeCreationEvent {
    std::vector<PixelId> pixelsC;
    AltitudeType b = 0;
    bool updateMaxTree = false;
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

TargetNodeCreationEvent findTargetNodeCreationEvent(const ImageUInt8Ptr &image,
                                                    const AdjacencyRelationPtr &adjacency,
                                                    bool updateMaxTree) {
    std::vector<std::vector<PixelId>> candidates;
    std::vector<PixelId> current;
    enumeratePixelCombinations(image->getSize(), 2, 0, current, candidates);

    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree &dualTree = updateMaxTree ? maxTree : minTree;

    for (const auto &region : candidates) {
        if (!problem1_debug::isConnectedRegion(*adjacency, std::span<const PixelId>(region))) {
            continue;
        }
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

void testOracleMatchesRebuildForCurrentCompatibleLowerEvent() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const std::vector<PixelId> pixelsC = {13, 14};

    const auto result = problem1_debug::runExplicitRegionProblem1Oracle<AltitudeType>(
        image,
        adjacency,
        false,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        true);

    require(result.status == problem1_debug::ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild,
            "compatible lower event must match rebuild oracle");
    require(result.ranLocalUpdate(), "oracle should run local update for compatible event");
    require(result.matchedRebuild(), "oracle should report rebuild match");
    require(result.expectedImage != nullptr, "oracle should keep expected rebuild image");
    require(result.localReconstruction != nullptr, "oracle should keep local reconstruction");
}

void testOracleClassifiesInvalidContractBeforeMutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const std::vector<PixelId> pixelsC = {0, 1};

    const auto result = problem1_debug::runExplicitRegionProblem1Oracle<AltitudeType>(
        image,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        true);

    require(result.status == problem1_debug::ExplicitRegionProblem1OracleStatus::RejectedByContract,
            "non-extensive max-tree event must be rejected by Problem 1 contract");
    require(!result.ranLocalUpdate(), "oracle must not run local update for rejected contract");
    require(result.expectedImage == nullptr, "rejected contract should not need rebuild image");
}

void testOracleClassifiesKnownUpdaterLimitation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    const std::vector<PixelId> pixelsC = {0, 1};

    const auto result = problem1_debug::runExplicitRegionProblem1Oracle<AltitudeType>(
        image,
        adjacency,
        true,
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(7),
        true);

    require(result.contract.acceptedByPublishedProblem,
            "raise event should satisfy the published Problem 1 hypothesis");
    require(result.status == problem1_debug::ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation,
            "oracle should classify unsupported accepted event as known updater limitation");
    require(!result.ranLocalUpdate(), "oracle must not run local update for known limitations");
    require(result.expectedImage != nullptr, "known limitation should still build rebuild oracle image");
    require(!result.message.empty(), "known limitation should explain why the updater is not run");
}

void testOracleMaterializesPartialOwnerThenMatchesRebuild() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2, 0);
    require(owner != InvalidNode,
            "test image must have a nonzero min-tree node with at least two direct proper parts");

    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};
    const auto b = static_cast<AltitudeType>(minTree.getAltitude(owner) - 1);

    const auto withoutMaterialization = problem1_debug::runExplicitRegionProblem1Oracle<AltitudeType>(
        image,
        adjacency,
        false,
        std::span<const PixelId>(pixelsC),
        b,
        true,
        false);

    require(withoutMaterialization.status == problem1_debug::ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation,
            "partial-owner event should be a known limitation before materialization");
    require(withoutMaterialization.contract.partialOwners == 1,
            "fixture should expose exactly one partial owner before materialization");

    const auto withMaterialization = problem1_debug::runExplicitRegionProblem1Oracle<AltitudeType>(
        image,
        adjacency,
        false,
        std::span<const PixelId>(pixelsC),
        b,
        true,
        true);

    require(withMaterialization.materializationAttempted,
            "oracle should attempt local materialization when enabled");
    require(withMaterialization.materializationApplied,
            "oracle should apply local materialization");
    require(withMaterialization.materializedOwners == 1,
            "oracle should materialize exactly one owner in this fixture");
    require(withMaterialization.materializedNodes >= 1,
            "oracle should create at least one local node during materialization");
    require(withMaterialization.hasPostMaterializationContract,
            "oracle should expose the post-materialization contract");
    require(withMaterialization.postMaterializationContract.partialOwners == 0,
            "post-materialization contract should no longer require materialization");
    require(withMaterialization.status == problem1_debug::ExplicitRegionProblem1OracleStatus::LocalUpdateMatchedRebuild,
            "materialize-then-update should match rebuild oracle");
}

void testOracleLocalCorrectnessGateRejectsMissingTargetPolicy() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    const auto event = findTargetNodeCreationEvent(image, adjacency, true);
    require(event.found,
            "fixture must expose a target-node-creation event");

    const auto result = problem1_debug::runExplicitRegionProblem1Oracle<AltitudeType>(
        image,
        adjacency,
        event.updateMaxTree,
        std::span<const PixelId>(event.pixelsC),
        event.b,
        true,
        true,
        false,
        false,
        true);

    require(result.status == problem1_debug::ExplicitRegionProblem1OracleStatus::KnownUpdaterLimitation,
            "enforced local-correctness gate should reject missing target creation");
    require(result.localCorrectnessValidationPerformed,
            "oracle should record local-correctness validation");
    require(!result.localCorrectnessValidationPassed,
            "missing target creation should fail local-correctness validation");
    require(result.rejectedByLocalCorrectnessValidation,
            "oracle should report rejection by the local-correctness gate");
    require(!result.localCorrectnessValidation.targetLevelPresentOrCreated,
            "validation should expose the failed target-level gate");
    require(!result.ranLocalUpdate(),
            "oracle must not run the local update after gate rejection");
}

} // namespace

int main() {
    testOracleMatchesRebuildForCurrentCompatibleLowerEvent();
    testOracleClassifiesInvalidContractBeforeMutation();
    testOracleClassifiesKnownUpdaterLimitation();
    testOracleMaterializesPartialOwnerThenMatchesRebuild();
    testOracleLocalCorrectnessGateRejectsMissingTargetPolicy();
    return 0;
}
