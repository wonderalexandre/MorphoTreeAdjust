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
#include "../include/morphoTreeAdjustResearch/debug/RegionMergeUpdater.hpp"

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

int countAliveNodes(const DynamicComponentTree &tree) {
    int count = 0;
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId)) {
            ++count;
        }
    }
    return count;
}

void requireTreeConsistency(const DynamicComponentTree &tree) {
    require(tree.getRoot() != InvalidNode, "root must be valid");
    require(tree.isAlive(tree.getRoot()), "root must be alive");
    require(tree.getNodeParent(tree.getRoot()) == tree.getRoot(), "root parent must be itself");

    std::vector<unsigned char> reachable(static_cast<size_t>(tree.getNumInternalNodeSlots()), 0);
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        reachable[static_cast<size_t>(nodeId)] = 1;
    }

    int totalProperParts = 0;
    std::vector<unsigned char> seen(static_cast<size_t>(tree.getNumTotalProperParts()), 0);
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }
        require(reachable[static_cast<size_t>(nodeId)] != 0, "alive node must be reachable from root");

        int countedChildren = 0;
        for (NodeId childId : tree.getChildren(nodeId)) {
            require(tree.isAlive(childId), "child must be alive");
            require(tree.getNodeParent(childId) == nodeId, "child parent must match");
            ++countedChildren;
        }
        require(countedChildren == tree.getNumChildren(nodeId), "child count mismatch");

        int countedProperParts = 0;
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            require(pixelId >= 0 && pixelId < tree.getNumTotalProperParts(), "proper part must be in bounds");
            require(tree.getSmallestComponent(pixelId) == nodeId, "proper part owner mismatch");
            require(seen[static_cast<size_t>(pixelId)] == 0, "proper part visited twice");
            seen[static_cast<size_t>(pixelId)] = 1;
            ++countedProperParts;
        }
        require(countedProperParts == tree.getNumProperParts(nodeId), "proper part count mismatch");
        totalProperParts += countedProperParts;
    }

    require(totalProperParts == tree.getNumTotalProperParts(), "global proper part count mismatch");
    require(std::all_of(seen.begin(), seen.end(), [](unsigned char value) { return value != 0; }),
            "proper part partition is incomplete");
}

void requireThrowsLocalMaterialization(DynamicComponentTree &minTree,
                                       DynamicComponentTree &maxTree,
                                       const AdjacencyRelationPtr &adj,
                                       const region_debug::RegionLevelChange<AltitudeType> &change) {
    try {
        (void) region_debug::RegionMergeUpdater::updateDualTreeOnly(minTree, maxTree, adj, change);
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error("partial owner region should require local materialization");
}

void testFullDirectOwnerDoesNotRequireMaterialization() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2);
    require(owner != InvalidNode, "test image must have a node with at least two direct proper parts");
    const auto pixelsC = collectProperParts(minTree, owner);

    auto plan = region_debug::analyzeRegionMaterialization(
        minTree,
        *adj,
        std::span<const PixelId>(pixelsC));

    require(plan.affectedOwners.size() == 1, "full direct owner should affect exactly one owner");
    require(plan.affectedOwners.front().nodeId == owner, "materialization plan should report the selected owner");
    require(plan.numFullyCoveredOwners() == 1, "selected owner should be fully covered");
    require(plan.numPartialOwners() == 0, "full direct owner should not be partial");
    require(!plan.requiresLocalMaterialization(), "full direct owner should not require local materialization");
    require(plan.affectedOwners.front().insideDirectProperParts.size() == pixelsC.size(),
            "inside direct proper parts should match C");
    require(plan.affectedOwners.front().outsideDirectProperParts.empty(),
            "full direct owner should not have outside direct proper parts");
}

void testPartialDirectOwnerRequiresMaterialization() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2);
    require(owner != InvalidNode, "test image must have a node with at least two direct proper parts");
    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};

    auto plan = region_debug::analyzeRegionMaterialization(
        minTree,
        *adj,
        std::span<const PixelId>(pixelsC));

    require(plan.affectedOwners.size() == 1, "partial direct owner should affect exactly one owner");
    require(plan.affectedOwners.front().nodeId == owner, "materialization plan should report the partial owner");
    require(plan.numFullyCoveredOwners() == 0, "partial owner should not be fully covered");
    require(plan.numPartialOwners() == 1, "partial owner should be counted");
    require(plan.requiresLocalMaterialization(), "partial owner should require local materialization");
    require(plan.affectedOwners.front().insideDirectProperParts.size() == 1,
            "partial owner should contain the selected direct proper part");
    require(plan.affectedOwners.front().outsideDirectProperParts.size() == ownerPixels.size() - 1,
            "partial owner should preserve the unselected direct proper parts outside C");
    require(plan.affectedOwners.front().insideComponents.size() == 1,
            "a single selected pixel should form one inside component");
}

void testDuplicatePixelsAreRejected() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsC = {0, 0};
    try {
        (void) region_debug::analyzeRegionMaterialization(minTree, *adj, std::span<const PixelId>(pixelsC));
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error("duplicated region pixels should be rejected");
}

void testMergeUpdaterRejectsMarkedMaterializedPartialOwner() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
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
        true);

    requireThrowsLocalMaterialization(minTree, maxTree, adj, change);
}

void testMaterializeSinglePartialOwnerPreservesReconstructionAndMakesRegionWhole() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree minTree(input->clone(), false, adj);
    DynamicComponentTree rebuildTree(input->clone(), false, adj);

    const NodeId owner = firstNodeWithAtLeastProperParts(minTree, 2);
    require(owner != InvalidNode, "test image must have a node with at least two direct proper parts");
    const auto ownerPixels = collectProperParts(minTree, owner);
    const std::vector<PixelId> pixelsC = {ownerPixels.front()};

    auto plan = region_debug::analyzeRegionMaterialization(
        minTree,
        *adj,
        std::span<const PixelId>(pixelsC));
    const int oldSlots = minTree.getNumInternalNodeSlots();
    const int oldAlive = countAliveNodes(minTree);

    const auto result = region_debug::materializeSinglePartialOwner(minTree, plan);

    require(result.originalOwner == owner, "materialization must report the original owner");
    require(result.createdInsideNodes.size() == plan.affectedOwners.front().insideComponents.size(),
            "materialization should create one inside node per inside component");
    require(result.numCreatedNodes() >= 1, "partial materialization should create at least one node");
    require(countAliveNodes(minTree) == oldAlive + result.numCreatedNodes(),
            "alive node count must increase by the number of created nodes");
    require(minTree.getNumInternalNodeSlots() >= oldSlots + result.numCreatedNodes(),
            "node slot domain should grow when no reusable slots are available");

    const NodeId materializedCNode = result.createdInsideNodes.front();
    require(minTree.isAlive(materializedCNode), "materialized C node must be alive");
    require(minTree.getNodeParent(materializedCNode) == owner, "materialized C node must be attached under the original owner");
    require(minTree.getAltitude(materializedCNode) == minTree.getAltitude(owner),
            "materialized node must preserve the owner's altitude before the level update");
    require(minTree.getSmallestComponent(pixelsC.front()) == materializedCNode,
            "selected pixel must move to the materialized C node");
    require(minTree.reconstructionImage()->isEqual(rebuildTree.reconstructionImage()),
            "local materialization must preserve reconstruction before the level update");
    requireTreeConsistency(minTree);

    auto postPlan = region_debug::analyzeRegionMaterialization(
        minTree,
        *adj,
        std::span<const PixelId>(pixelsC));
    require(!postPlan.requiresLocalMaterialization(),
            "after local materialization, C should be a union of complete direct owners");
    require(postPlan.numFullyCoveredOwners() == 1,
            "after local materialization, the selected C node should be fully covered");
    require(postPlan.affectedOwners.front().nodeId == materializedCNode,
            "post-materialization plan should point C to the new materialized node");
}

} // namespace

int main() {
    testFullDirectOwnerDoesNotRequireMaterialization();
    testPartialDirectOwnerRequiresMaterialization();
    testDuplicatePixelsAreRejected();
    testMergeUpdaterRejectsMarkedMaterializedPartialOwner();
    testMaterializeSinglePartialOwnerPreservesReconstructionAndMakesRegionWhole();
    return 0;
}
