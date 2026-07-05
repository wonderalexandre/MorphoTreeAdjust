#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionMaterialization.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionMergeWorkspace.hpp"

namespace region_debug = debug::region;

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

std::vector<PixelId> collectProperParts(const DynamicComponentTree &tree, NodeId nodeId) {
    std::vector<PixelId> pixels;
    for (PixelId pixelId : tree.getProperParts(nodeId)) {
        pixels.push_back(pixelId);
    }
    return pixels;
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

void requirePaperLabels(const std::string &debug) {
    require(debug.find("MergeNodes(C, -)") != std::string::npos,
            "debug string must expose MergeNodes(C, -)");
    require(debug.find("MergeNodes([a,b], -)") != std::string::npos,
            "debug string must expose MergeNodes([a,b], -)");
    require(debug.find("MergeNodes(v, -)") != std::string::npos,
            "debug string must expose MergeNodes(v, -)");
    require(debug.find("MergedNode(v, +)") != std::string::npos,
            "debug string must expose MergedNode(v, +)");
    require(debug.find("affectedThresholds") != std::string::npos,
            "debug string must expose affected thresholds");
}

void requireNonemptyLevelWorkspace(
    const region_debug::RegionMergeWorkspace<AltitudeType> &workspace) {
    require(!workspace.mergeNodesCMinus.empty(),
            "workspace should record adjacent seeds in MergeNodes(C, -)");
    require(!workspace.mergeNodesIntervalMinus.empty(),
            "workspace should record MergeNodes([a,b], -)");
    require(!workspace.mergeNodesByLevel.empty(),
            "workspace should record per-threshold buckets");
    require(workspace.affectedThresholds.size() == workspace.mergeNodesByLevel.size(),
            "affected thresholds should match per-level buckets");

    for (const auto &bucket : workspace.mergeNodesByLevel) {
        require(!bucket.mergeNodesVMinus.empty(), "per-level merge bucket should not be empty");
        require(bucket.mergedNodeVPlusCandidate == bucket.mergeNodesVMinus.front(),
                "read-only MergedNode(v, +) candidate should name the first bucket node");
        const auto *bucketByLevel = workspace.mergeNodesVMinus(bucket.level);
        require(bucketByLevel != nullptr, "level lookup should locate MergeNodes(v,-)");
        require(*bucketByLevel == bucket.mergeNodesVMinus,
                "level lookup should match the stored bucket");
        require(workspace.mergedNodeVPlusCandidate(bucket.level) == bucket.mergedNodeVPlusCandidate,
                "candidate lookup should match the stored MergedNode(v, +) candidate");
    }
}

void testMaterializedLowerRegionWorkspaceUsesPaperNames() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(6),
        static_cast<AltitudeType>(4),
        true,
        true);

    auto workspace = region_debug::buildRegionMergeWorkspace(minTree, *adj, change);

    require(workspace.nodeCa != InvalidNode, "workspace should locate nodeCa");
    require(workspace.requestedA == 6, "workspace should keep the requested a");
    require(workspace.a == 6, "workspace should expose the actual nodeCa altitude as a");
    require(workspace.b == 4, "workspace should expose b");
    require(workspace.direction == region_debug::RegionUpdateDirection::Lower,
            "workspace direction should be computed from actual a and b");
    require(!workspace.isMaxtree, "demo workspace should target the min-tree");
    requireNonemptyLevelWorkspace(workspace);
    requirePaperLabels(workspace.debugString());
}

void testWorkspaceAfterMaterializingPartialOwnerUsesNewNodeCa() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2, 0);
    require(owner != InvalidNode, "test image must have a nonzero node with at least two direct proper parts");

    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};
    const auto a = static_cast<AltitudeType>(minTree.getAltitude(owner));
    const auto b = static_cast<AltitudeType>(a - 1);
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        a,
        b,
        true,
        false);

    const auto prePlan = region_debug::analyzeRegionMaterialization(minTree, *adj, change);
    require(prePlan.requiresLocalMaterialization(),
            "single-pixel region should cut its original direct owner");
    const auto materialized = region_debug::materializePartialOwners(minTree, prePlan);
    require(materialized.numCreatedNodes() >= 1,
            "partial-owner materialization should create a node for C");

    const auto postPlan = region_debug::analyzeRegionMaterialization(minTree, *adj, change);
    require(!postPlan.requiresLocalMaterialization(),
            "materialized C should be represented by complete direct owners");
    require(postPlan.affectedOwners.size() == 1,
            "single-pixel materialized C should have one direct owner");

    auto workspace = region_debug::buildRegionMergeWorkspace(minTree, *adj, change);

    require(workspace.nodeCa == postPlan.affectedOwners.front().nodeId,
            "workspace nodeCa should be the materialized owner of C");
    require(workspace.a == a, "materialization should preserve a before the lower update");
    require(workspace.b == b, "workspace should preserve the requested target level b");
    requireNonemptyLevelWorkspace(workspace);
}

} // namespace

int main() {
    testMaterializedLowerRegionWorkspaceUsesPaperNames();
    testWorkspaceAfterMaterializingPartialOwnerUsesNewNodeCa();
    return 0;
}
