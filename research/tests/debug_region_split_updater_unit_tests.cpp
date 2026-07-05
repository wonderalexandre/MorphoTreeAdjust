#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChangeRebuildUpdater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionSplitUpdater.hpp"

namespace region_debug = debug::region;

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr makeImage(int rows, int cols, const std::vector<uint8_t> &raw) {
    auto image = ImageUInt8::create(rows, cols);
    require(static_cast<int>(raw.size()) == image->getSize(), "raw image size must match dimensions");
    for (int i = 0; i < image->getSize(); ++i) {
        (*image)[i] = raw[static_cast<size_t>(i)];
    }
    return image;
}

region_debug::RegionLevelChange<AltitudeType> makeLowerChange(const std::vector<PixelId> &pixelsC,
                                                              AltitudeType a,
                                                              AltitudeType b,
                                                              bool assumeConnected = false,
                                                              bool assumeMaterialized = false) {
    return region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        a,
        b,
        assumeConnected,
        assumeMaterialized);
}

void requireAllEntriesHaveOutcome(const region_debug::RegionSplitPlan &plan,
                                  region_debug::RegionSplitLocalOutcome outcome,
                                  const std::string &message) {
    require(!plan.affectedComponents.empty(), "plan must contain affected components");
    for (const auto &entry : plan.affectedComponents) {
        require(entry.outcome == outcome, message);
    }
}

bool componentsEqual(const std::vector<std::vector<PixelId>> &lhs,
                     const std::vector<std::vector<PixelId>> &rhs) {
    return lhs == rhs;
}

std::vector<PixelId> removedPixelsAtLevel(const region_debug::RegionSplitPlan &plan,
                                          int level) {
    for (const auto &entry : plan.affectedComponents) {
        if (entry.level == level) {
            return entry.removedPixels;
        }
    }
    throw std::runtime_error("expected an affected entry at the requested level");
}

std::unordered_set<PixelId> makePixelSet(const std::vector<PixelId> &pixels) {
    return std::unordered_set<PixelId>(pixels.begin(), pixels.end());
}

void requireNoEmptyLiveNodes(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        require(tree.getNumProperParts(nodeId) > 0,
                "local split update must not leave live nodes without proper parts");
    }
}

void requireStrictMaxTreeAltitudeOrder(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        for (NodeId childId : tree.getChildren(nodeId)) {
            require(tree.getAltitude(nodeId) < tree.getAltitude(childId),
                    "max-tree parent altitude must be strictly lower than child altitude");
        }
    }
}

void requireMatchesRebuildMaxTree(const ImageUInt8Ptr &input,
                                  const AdjacencyRelationPtr &adj,
                                  const region_debug::RegionLevelChange<AltitudeType> &change,
                                  const DynamicComponentTree &updatedMaxTree) {
    auto expectedImage = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    DynamicComponentTree rebuiltMaxTree(expectedImage->clone(), true, adj);

    require(updatedMaxTree.reconstructionImage()->isEqual(expectedImage),
            "local split update reconstruction must match the globally updated image");
    require(updatedMaxTree.reconstructionImage()->isEqual(rebuiltMaxTree.reconstructionImage()),
            "local split update reconstruction must match rebuilt max-tree reconstruction");
    require(updatedMaxTree.getNumNodes() == rebuiltMaxTree.getNumNodes(),
            "restricted local split update should match rebuilt max-tree node count on synthetic fixtures");
    requireNoEmptyLiveNodes(updatedMaxTree);
    requireStrictMaxTreeAltitudeOrder(updatedMaxTree);
}

template<typename Fn>
void requireThrows(Fn &&fn, const std::string &message) {
    try {
        fn();
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error(message);
}

int findNativeNodePlanIndex(const region_debug::RegionSplitLowerNativeRecompositionPlan &plan,
                            int level,
                            const std::vector<PixelId> &supportPixels) {
    for (size_t i = 0; i < plan.nodes.size(); ++i) {
        if (plan.nodes[i].level == level &&
            plan.nodes[i].supportPixels == supportPixels) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void testFlatLowerRemovalPlan() {
    auto input = makeImage(1, 3, {1, 5, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 1);
    auto plan = region_debug::RegionSplitUpdater::buildFlatLowerPlan(maxTree, *adj, change);

    require(plan.firstAffectedLevel == 2, "removal plan first affected level must be b+1");
    require(plan.lastAffectedLevel == 5, "removal plan last affected level must be a");
    require(static_cast<int>(plan.affectedComponents.size()) == 4,
            "single-pixel removal should affect one component at each level 2..5");
    requireAllEntriesHaveOutcome(plan,
                                 region_debug::RegionSplitLocalOutcome::Removed,
                                 "single-pixel high component should be removed at every affected level");
    for (const auto &entry : plan.affectedComponents) {
        require(entry.remainingComponents.empty(), "removed entries must have no remaining components");
    }
}

void testFlatLowerRemovalUpdate() {
    auto input = makeImage(1, 3, {1, 5, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 1);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerSingleOwner(maxTree, *adj, change);

    require(!result.createdTargetNode, "removal update should reuse the existing parent at level b");
    require(result.releasedSourceNode, "removal update must release the emptied source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerContractionPlan() {
    auto input = makeImage(1, 2, {5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0}, 5, 3);
    auto plan = region_debug::RegionSplitUpdater::buildFlatLowerPlan(maxTree, *adj, change);

    require(plan.firstAffectedLevel == 4, "contraction plan first affected level must be b+1");
    require(plan.lastAffectedLevel == 5, "contraction plan last affected level must be a");
    require(static_cast<int>(plan.affectedComponents.size()) == 2,
            "contraction fixture should affect one component at each level 4..5");
    requireAllEntriesHaveOutcome(plan,
                                 region_debug::RegionSplitLocalOutcome::Contracted,
                                 "removing the shell should leave one preserved child component");
    for (const auto &entry : plan.affectedComponents) {
        require(componentsEqual(entry.remainingComponents, {{1}}),
                "contracted entries must keep the high child as the only component");
    }
}

void testFlatLowerContractionUpdate() {
    auto input = makeImage(1, 2, {5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0}, 5, 3);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerSingleOwner(maxTree, *adj, change);

    require(result.createdTargetNode, "contraction update should create a new node at b");
    require(result.releasedSourceNode, "contraction update must release the emptied source node");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "contraction update should create the new root at b");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerSplitPlan() {
    auto input = makeImage(1, 3, {6, 5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3);
    auto plan = region_debug::RegionSplitUpdater::buildFlatLowerPlan(maxTree, *adj, change);

    require(plan.firstAffectedLevel == 4, "split plan first affected level must be b+1");
    require(plan.lastAffectedLevel == 5, "split plan last affected level must be a");
    require(static_cast<int>(plan.affectedComponents.size()) == 2,
            "split fixture should affect one component at each level 4..5");
    requireAllEntriesHaveOutcome(plan,
                                 region_debug::RegionSplitLocalOutcome::Split,
                                 "removing the bridge shell should split the component");
    for (const auto &entry : plan.affectedComponents) {
        require(componentsEqual(entry.remainingComponents, {{0}, {2}}),
                "split entries must keep the two high children as separate components");
    }
}

void testFlatLowerSplitUpdate() {
    auto input = makeImage(1, 3, {6, 5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerSingleOwner(maxTree, *adj, change);

    require(result.createdTargetNode, "split update should create a new node at b");
    require(result.releasedSourceNode, "split update must release the emptied source node");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "split update should create the new root at b");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerPlanNonFlatThresholds() {
    auto input = makeImage(1, 3, {6, 5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 1}, 6, 3);
    const auto plan = region_debug::RegionSplitUpdater::buildLowerPlan(maxTree, *adj, change);

    require(plan.firstAffectedLevel == 4,
            "non-flat lower plan first affected level must be b+1");
    require(plan.lastAffectedLevel == 6,
            "non-flat lower plan last affected level must be max I(C)");
    require(static_cast<int>(plan.affectedComponents.size()) == 3,
            "non-flat fixture should affect root at levels 4..5 and the high child at level 6");
    require(removedPixelsAtLevel(plan, 4) == std::vector<PixelId>({0, 1}),
            "D_4 must contain every C pixel with value at least 4");
    require(removedPixelsAtLevel(plan, 5) == std::vector<PixelId>({0, 1}),
            "D_5 must contain every C pixel with value at least 5");
    require(removedPixelsAtLevel(plan, 6) == std::vector<PixelId>({0}),
            "D_6 must contain only the C pixels with value at least 6");
    require(plan.affectedComponents[0].outcome == region_debug::RegionSplitLocalOutcome::Contracted,
            "removing D_4 should contract the root component to the untouched high child");
    require(plan.affectedComponents[1].outcome == region_debug::RegionSplitLocalOutcome::Contracted,
            "removing D_5 should contract the root component to the untouched high child");
    require(plan.affectedComponents[2].outcome == region_debug::RegionSplitLocalOutcome::Removed,
            "removing D_6 should remove the high child inside C");
}

void testLowerPlanRejectsWrongMaxA() {
    auto input = makeImage(1, 3, {6, 5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 1}, 5, 3);
    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::buildLowerPlan(maxTree, *adj, change);
        },
        "non-flat lower plan must require a to be max I(C)");
}

void testLowerUpdateByFlatLevelBatches() {
    auto input = makeImage(1, 3, {6, 5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 1}, 6, 3);
    const auto result = region_debug::RegionSplitUpdater::updateLowerByFlatLevelBatches(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "non-flat lower update should process one flat batch per original level");
    require(result.numCreatedTargetNodes() == 1,
            "non-flat lower update should create the b-node in the first flat batch");
    require(result.numReleasedSourceNodes() == 2,
            "non-flat lower update should release the level-5 and level-6 source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "non-flat lower update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerUpdateByFlatLevelBatchesRejectsBeforePartialMutation() {
    auto input = makeImage(1, 8, {0, 5, 4, 5, 3, 5, 4, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2, 3, 5}, 5, 2, true);
    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::updateLowerByFlatLevelBatches(
                maxTree,
                *adj,
                change);
        },
        "non-flat lower update must reject unsupported later batches during preflight");
    require(maxTree.reconstructionImage()->isEqual(input),
            "non-flat lower update preflight failure must leave the tree unchanged");
}

void testFlatLowerMaterializesPartialRootOwnerAsSource() {
    auto input = makeImage(1, 3, {5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3);
    const auto result = region_debug::RegionSplitUpdater::materializeFlatLowerRegion(
        maxTree,
        *adj,
        change);

    require(result.changed,
            "partial-owner materialization should report that it changed the tree");
    require(result.numCreatedNodes() == 2,
            "partial root materialization should create one outside node per remaining component");
    require(result.reconstruction->isEqual(input),
            "partial-owner materialization must preserve the current reconstruction");
    require(region_debug::RegionSplitUpdater::isMaterializedAsProperParts(maxTree, makePixelSet({1})),
            "partial-owner materialization must make C a complete direct proper part");
    require(maxTree.getSmallestComponent(1) == result.originalOwner,
            "the original owner should become the source node that directly owns C");
    require(maxTree.getNumProperParts(result.originalOwner) == 1,
            "the original owner should keep only C as direct proper part");
    requireNoEmptyLiveNodes(maxTree);
}

void testFlatLowerUpdateAfterMaterializingPartialRootOwner() {
    auto input = makeImage(1, 3, {5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "materialized partial-owner update should process the single connected C component");
    require(result.numCreatedTargetNodes() == 1,
            "materialized partial-owner update should create the b-node");
    require(result.numReleasedSourceNodes() == 1,
            "materialized partial-owner update should release the temporary source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "materialized partial-owner result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerUpdateAfterMaterializingPartialOwnerWithChild() {
    auto input = makeImage(1, 4, {5, 6, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2}, 5, 3);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "materialized partial-owner update should process the bridge component with a preserved child");
    require(result.numCreatedTargetNodes() == 1,
            "materialized partial-owner update with child should create the b-node");
    require(result.numReleasedSourceNodes() == 1,
            "materialized partial-owner update with child should release the temporary source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "materialized partial-owner child result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerMaterializesDisconnectedInsideComponentsAsOneSource() {
    auto input = makeImage(1, 3, {5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 2}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::materializeFlatLowerRegion(
        maxTree,
        *adj,
        change);

    require(result.changed,
            "disconnected inside materialization should report that it changed the tree");
    require(result.numCreatedNodes() == 1,
            "disconnected inside materialization should create one outside node for the middle component");
    require(result.reconstruction->isEqual(input),
            "disconnected inside materialization must preserve the current reconstruction");
    require(region_debug::RegionSplitUpdater::isMaterializedAsProperParts(maxTree, makePixelSet({0, 2})),
            "disconnected inside materialization must make C a complete direct proper part");
    require(maxTree.getSmallestComponent(0) == result.originalOwner,
            "the first inside component should remain owned by the original source");
    require(maxTree.getSmallestComponent(2) == result.originalOwner,
            "the second inside component should remain owned by the original source");
    require(maxTree.getNumProperParts(result.originalOwner) == 2,
            "the original owner should keep both disconnected C pixels as direct proper parts");
    requireNoEmptyLiveNodes(maxTree);
}

void testFlatLowerUpdateAfterMaterializingDisconnectedInsideComponents() {
    auto input = makeImage(1, 3, {5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 2}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "disconnected inside update should dispatch as one materialized source node");
    require(result.numCreatedTargetNodes() == 1,
            "disconnected inside update should create the b-node");
    require(result.numReleasedSourceNodes() == 1,
            "disconnected inside update should release the materialized source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "disconnected inside result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerMaterializesMultiplePartialOwnersAsSources() {
    auto input = makeImage(1, 5, {5, 5, 0, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::materializeFlatLowerRegion(
        maxTree,
        *adj,
        change);

    require(result.changed,
            "multiple partial-owner materialization should report that it changed the tree");
    require(static_cast<int>(result.originalOwners.size()) == 2,
            "multiple partial-owner materialization should record both original owners");
    require(result.numCreatedNodes() == 2,
            "multiple partial-owner materialization should create one outside node per owner");
    require(result.reconstruction->isEqual(input),
            "multiple partial-owner materialization must preserve the current reconstruction");
    require(region_debug::RegionSplitUpdater::isMaterializedAsProperParts(maxTree, makePixelSet({1, 3})),
            "multiple partial-owner materialization must make C a union of complete direct proper parts");
    require(maxTree.getNumProperParts(maxTree.getSmallestComponent(1)) == 1,
            "first partial source should keep only its C pixel");
    require(maxTree.getNumProperParts(maxTree.getSmallestComponent(3)) == 1,
            "second partial source should keep only its C pixel");
    requireNoEmptyLiveNodes(maxTree);
}

void testFlatLowerUpdateAfterMaterializingMultiplePartialOwners() {
    auto input = makeImage(1, 5, {5, 5, 0, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "multiple partial-owner update should process both disconnected C components");
    require(result.numCreatedTargetNodes() == 2,
            "multiple partial-owner update should create one b-node per component");
    require(result.numReleasedSourceNodes() == 2,
            "multiple partial-owner update should release both materialized source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "multiple partial-owner result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerUpdateAfterMaterializingMixedPartialAndCoveredOwners() {
    auto input = makeImage(1, 5, {5, 5, 0, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "mixed partial and covered update should process both disconnected C components");
    require(result.numCreatedTargetNodes() == 2,
            "mixed partial and covered update should create one b-node per component");
    require(result.numReleasedSourceNodes() == 2,
            "mixed partial and covered update should release both source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "mixed partial and covered result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerUpdateAfterMaterializingMixedSourcesWithRepeatedOwner() {
    auto input = makeImage(1, 6, {5, 5, 0, 5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 5}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "source-grouped update should process the two materialized source nodes");
    require(result.numCreatedTargetNodes() == 2,
            "source-grouped update should create one b-node per independent source group");
    require(result.numReleasedSourceNodes() == 2,
            "source-grouped update should release both materialized source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "source-grouped result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerSourceGroupedPreservedAncestorSplitAfterMaterialization() {
    auto input = makeImage(1, 7, {0, 5, 5, 5, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 5}, 5, 2, true);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "source-grouped preserved-ancestor split should process the two materialized source nodes");
    require(result.numCreatedTargetNodes() == 1,
            "source-grouped preserved-ancestor split should create one shared b-node");
    require(result.numReleasedSourceNodes() == 2,
            "source-grouped preserved-ancestor split should release both materialized source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "source-grouped preserved-ancestor split result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerMaterializedSourceGroupedInteractingInsertedAncestorsSplitsInternalChildUpdate() {
    auto input = makeImage(1, 8, {0, 5, 5, 5, 4, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 6}, 5, 2, true);
    DynamicComponentTree materializedTree = maxTree;
    (void) region_debug::RegionSplitUpdater::materializeFlatLowerRegion(
        materializedTree,
        *adj,
        change);
    const auto groupedPlan = region_debug::RegionSplitUpdater::buildFlatLowerSourceGroupedBatchPlan(
        materializedTree,
        *adj,
        change);

    require(groupedPlan.numComponents() == 2,
            "source-grouped interacting fixture should group the two direct sources");
    require(groupedPlan.interaction == region_debug::RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors,
            "source-grouped plan should classify nested inserted ancestors as interacting");
    require(groupedPlan.candidateForFlatLowerLocalUpdater(),
            "interacting source groups should now be attempted by the restricted joint updater");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "materialized interacting split-child update should process the two materialized source nodes");
    require(result.numCreatedTargetNodes() == 1,
            "materialized interacting split-child update should create one shared b-node");
    require(result.numReleasedSourceNodes() == 2,
            "materialized interacting split-child update should release both source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "materialized interacting split-child result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerMaterializedSourceGroupedInteractingInsertedAncestorsUpdate() {
    auto input = makeImage(1, 8, {0, 5, 5, 5, 4, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2, 6}, 5, 2, true);
    DynamicComponentTree materializedTree = maxTree;
    (void) region_debug::RegionSplitUpdater::materializeFlatLowerRegion(
        materializedTree,
        *adj,
        change);
    const auto groupedPlan = region_debug::RegionSplitUpdater::buildFlatLowerSourceGroupedBatchPlan(
        materializedTree,
        *adj,
        change);

    require(groupedPlan.numComponents() == 2,
            "source-grouped materialized interacting update fixture should group the two direct sources");
    require(groupedPlan.interaction == region_debug::RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors,
            "source-grouped materialized update fixture should classify nested inserted ancestors as interacting");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "materialized interacting inserted ancestor update should process the two materialized source nodes");
    require(result.numCreatedTargetNodes() == 1,
            "materialized interacting inserted ancestor update should create one shared b-node");
    require(result.numReleasedSourceNodes() == 2,
            "materialized interacting inserted ancestor update should release both source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "materialized interacting inserted ancestor result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerInteractingInsertedAncestorsRecursiveInternalSplitUpdate() {
    auto input = makeImage(1, 10, {0, 7, 6, 7, 6, 5, 7, 4, 7, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 6, 8}, 7, 3, true);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 4,
            "recursive interacting fixture should identify the four high source nodes");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors,
            "recursive interacting fixture should classify nested inserted ancestors as interacting");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 4,
            "recursive interacting update should process the four source nodes");
    require(result.numCreatedTargetNodes() == 1,
            "recursive interacting update should create one shared b-node");
    require(result.numReleasedSourceNodes() == 4,
            "recursive interacting update should release all high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "recursive interacting result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerUpdateByFlatLevelBatchesWithMaterialization() {
    auto input = makeImage(1, 4, {6, 5, 5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 1}, 6, 3);
    const auto result = region_debug::RegionSplitUpdater::updateLowerByFlatLevelBatchesWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "non-flat materialized update should process the level-5 and level-6 batches");
    require(result.numCreatedTargetNodes() == 1,
            "non-flat materialized update should create the b-node in the first batch");
    require(result.numReleasedSourceNodes() == 2,
            "non-flat materialized update should release both level source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "non-flat materialized result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerUpdateByFlatLevelBatchesWithSourceGroupedMaterialization() {
    auto input = makeImage(1, 6, {5, 5, 0, 5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 5}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerByFlatLevelBatchesWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "source-grouped level-batch update should process the two materialized source nodes");
    require(result.numCreatedTargetNodes() == 2,
            "source-grouped level-batch update should create one b-node per independent source group");
    require(result.numReleasedSourceNodes() == 2,
            "source-grouped level-batch update should release both materialized source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "source-grouped level-batch result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionPlanForNestedNonFlatSource() {
    auto input = makeImage(1, 5, {0, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2}, 6, 3, true);
    const auto plan = region_debug::RegionSplitUpdater::buildLowerNativeRecompositionPlan(
        maxTree,
        *adj,
        change);

    require(plan.affectedLevels == std::vector<int>({3, 4, 5, 6}),
            "native non-flat plan should include b and every threshold up to a");
    require(plan.levelBEnvelopeComponents.size() == 1,
            "native non-flat plan should find one level-b envelope component");
    require(plan.levelBEnvelopeComponents.front() == std::vector<PixelId>({1, 2, 3}),
            "native non-flat level-b envelope should be the updated component containing C");
    require(plan.nodes.size() == 2,
            "native non-flat plan should keep only nodes with direct proper parts");
    require(plan.numConnectorCandidates() == 0,
            "native non-flat plan must not create connector candidates");
    require(plan.numRootNodes() == 1,
            "native non-flat plan should have one local root at b");

    const int targetIndex = findNativeNodePlanIndex(plan, 3, {1, 2, 3});
    const int childIndex = findNativeNodePlanIndex(plan, 5, {3});
    require(targetIndex >= 0, "native non-flat plan should contain the b target node");
    require(childIndex >= 0, "native non-flat plan should contain the preserved level-5 child");

    const auto &targetNode = plan.nodes[static_cast<size_t>(targetIndex)];
    require(targetNode.directProperParts == std::vector<PixelId>({1, 2}),
            "native b target should directly own the lowered pixels");
    require(targetNode.touchesRegion && targetNode.touchesPreserved,
            "native b target should touch both C and the preserved child support");
    require(targetNode.parentPlanIndex == -1,
            "native b target should be a local root");

    const auto &childNode = plan.nodes[static_cast<size_t>(childIndex)];
    require(childNode.directProperParts == std::vector<PixelId>({3}),
            "native preserved child should keep its direct proper part");
    require(childNode.parentPlanIndex == targetIndex,
            "native preserved child should attach below the b target");
}

void testLowerNativeRecompositionPlanForRecursiveSplitFixture() {
    auto input = makeImage(1, 10, {0, 7, 6, 7, 6, 5, 7, 4, 7, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 6, 8}, 7, 3, true);
    const auto plan = region_debug::RegionSplitUpdater::buildLowerNativeRecompositionPlan(
        maxTree,
        *adj,
        change);

    require(plan.affectedLevels == std::vector<int>({3, 4, 5, 6, 7}),
            "recursive native plan should include every threshold from b to a");
    require(plan.levelBEnvelopeComponents.size() == 1,
            "recursive native plan should find one updated level-b envelope");
    require(plan.levelBEnvelopeComponents.front() == std::vector<PixelId>({1, 2, 3, 4, 5, 6, 7, 8}),
            "recursive native envelope should cover the full local updated component");
    require(plan.nodes.size() == 5,
            "recursive native plan should contain the b target and four preserved nodes");
    require(plan.numConnectorCandidates() == 0,
            "recursive native plan must not keep empty connector nodes");
    require(plan.numRootNodes() == 1,
            "recursive native plan should have one local root at b");

    const int targetIndex = findNativeNodePlanIndex(plan, 3, {1, 2, 3, 4, 5, 6, 7, 8});
    const int isolatedHighIndex = findNativeNodePlanIndex(plan, 6, {2});
    const int branchRootIndex = findNativeNodePlanIndex(plan, 5, {4, 5});
    const int branchHighIndex = findNativeNodePlanIndex(plan, 6, {4});
    const int sideIndex = findNativeNodePlanIndex(plan, 4, {7});

    require(targetIndex >= 0, "recursive native plan should contain the b target");
    require(isolatedHighIndex >= 0, "recursive native plan should contain the isolated level-6 child");
    require(branchRootIndex >= 0, "recursive native plan should contain the level-5 branch");
    require(branchHighIndex >= 0, "recursive native plan should contain the level-6 child inside the branch");
    require(sideIndex >= 0, "recursive native plan should contain the level-4 side component");

    require(plan.nodes[static_cast<size_t>(targetIndex)].directProperParts ==
                std::vector<PixelId>({1, 3, 6, 8}),
            "recursive native b target should directly own all lowered pixels");
    require(plan.nodes[static_cast<size_t>(isolatedHighIndex)].parentPlanIndex == targetIndex,
            "isolated high preserved component should attach directly to b");
    require(plan.nodes[static_cast<size_t>(branchRootIndex)].parentPlanIndex == targetIndex,
            "level-5 preserved branch should attach directly to b");
    require(plan.nodes[static_cast<size_t>(branchHighIndex)].parentPlanIndex == branchRootIndex,
            "level-6 branch child should attach below the level-5 branch");
    require(plan.nodes[static_cast<size_t>(sideIndex)].parentPlanIndex == targetIndex,
            "level-4 side component should attach directly to b");
}

void testLowerNativeRecompositionUpdateForNestedNonFlatSource() {
    auto input = makeImage(1, 5, {0, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2}, 6, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numCreatedNodes() == 2,
            "native non-flat update should create the b target and preserved child");
    require(result.releasedEnvelopeSubtreeNodes == 2,
            "native non-flat update should replace the old level-5/level-6 envelope subtree");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native non-flat update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionUpdateForRecursiveSplitFixture() {
    auto input = makeImage(1, 10, {0, 7, 6, 7, 6, 5, 7, 4, 7, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 6, 8}, 7, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numCreatedNodes() == 5,
            "recursive native update should create the b target and four preserved nodes");
    require(result.releasedEnvelopeSubtreeNodes == 7,
            "recursive native update should replace the full old envelope subtree");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "recursive native update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionPreservesPixelsAboveA() {
    auto input = makeImage(1, 5, {0, 5, 8, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3, true);
    const auto plan = region_debug::RegionSplitUpdater::buildLowerNativeRecompositionPlan(
        maxTree,
        *adj,
        change);

    require(plan.affectedLevels == std::vector<int>({3, 4, 5}),
            "native high-preserved plan should keep affected levels as b..a");
    const int targetIndex = findNativeNodePlanIndex(plan, 3, {1, 2, 3});
    const int levelFiveIndex = findNativeNodePlanIndex(plan, 5, {2, 3});
    const int levelEightIndex = findNativeNodePlanIndex(plan, 8, {2});
    require(targetIndex >= 0,
            "native high-preserved plan should contain the b target");
    require(levelFiveIndex >= 0,
            "native high-preserved plan should contain the preserved level-5 component");
    require(levelEightIndex >= 0,
            "native high-preserved plan should contain the preserved level-8 child");
    require(plan.nodes[static_cast<size_t>(levelFiveIndex)].parentPlanIndex == targetIndex,
            "preserved level-5 component should attach below b");
    require(plan.nodes[static_cast<size_t>(levelEightIndex)].parentPlanIndex == levelFiveIndex,
            "preserved level-8 child should attach below the level-5 component");

    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numCreatedNodes() == 3,
            "native high-preserved update should create b, level-5, and level-8 nodes");
    require(result.releasedEnvelopeSubtreeNodes == 2,
            "native high-preserved update should replace the old level-5/level-8 envelope subtree");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native high-preserved result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionReusesRootEnvelope() {
    auto input = makeImage(1, 3, {5, 6, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);
    const NodeId oldRoot = maxTree.getRoot();

    const auto change = makeLowerChange({0, 1}, 6, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.replacedEnvelopeNode == oldRoot,
            "native root-envelope update should report the old root as the replaced envelope");
    require(result.externalParentNode == InvalidNode,
            "native root-envelope update should not report an external parent");
    require(maxTree.getRoot() == oldRoot,
            "native root-envelope update should reuse the old root as the new b-node");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "native root-envelope update should retag the reused root at b");
    require(result.releasedEnvelopeSubtreeNodes == 1,
            "native root-envelope update should release the old child subtree");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native root-envelope result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionReusesRootEnvelopeWithPreservedPixelAboveA() {
    auto input = makeImage(1, 3, {5, 8, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);
    const NodeId oldRoot = maxTree.getRoot();

    const auto change = makeLowerChange({0}, 5, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.replacedEnvelopeNode == oldRoot,
            "native root-envelope high-preserved update should reuse the old root");
    require(maxTree.getRoot() == oldRoot,
            "native root-envelope high-preserved update should keep the old root id");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "native root-envelope high-preserved update should retag the root at b");
    require(result.numCreatedNodes() == 3,
            "native root-envelope high-preserved update should materialize b, level-5, and level-8 nodes");
    require(result.releasedEnvelopeSubtreeNodes == 1,
            "native root-envelope high-preserved update should release the old level-8 child");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native root-envelope high-preserved result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionHandlesMultipleEnvelopes() {
    auto input = makeImage(1, 9, {0, 5, 6, 5, 0, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2, 5, 6}, 6, 3, true);
    const auto plan = region_debug::RegionSplitUpdater::buildLowerNativeRecompositionPlan(
        maxTree,
        *adj,
        change);

    require(plan.levelBEnvelopeComponents.size() == 2,
            "native multi-envelope plan should find two level-b envelopes");
    require(plan.numRootNodes() == 2,
            "native multi-envelope plan should have one local root per envelope");
    require(plan.numConnectorCandidates() == 0,
            "native multi-envelope plan should not create connector candidates");

    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numReplacedEnvelopeNodes() == 2,
            "native multi-envelope update should replace both old envelope nodes");
    require(result.numCreatedNodes() == 4,
            "native multi-envelope update should create two b nodes and two preserved children");
    require(result.releasedEnvelopeSubtreeNodes == 4,
            "native multi-envelope update should release both old level-5/level-6 subtrees");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native multi-envelope result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionUpdatesPartialRootOwnerWithoutMaterialization() {
    auto input = makeImage(1, 3, {5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3, true);
    require(!region_debug::RegionSplitUpdater::isMaterializedAsProperParts(maxTree, makePixelSet({1})),
            "native partial-root fixture should start non-materialized");

    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numReplacedEnvelopeNodes() == 1,
            "native partial-root update should replace the single root envelope");
    require(result.numCreatedNodes() == 3,
            "native partial-root update should create the b root and two preserved level-5 children");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native partial-root result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionUpdatesPartialOwnerWithPreservedChildWithoutMaterialization() {
    auto input = makeImage(1, 4, {5, 6, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2}, 5, 3, true);
    require(!region_debug::RegionSplitUpdater::isMaterializedAsProperParts(maxTree, makePixelSet({2})),
            "native partial-owner child fixture should start non-materialized");

    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numReplacedEnvelopeNodes() == 1,
            "native partial-owner child update should replace the single root envelope");
    require(result.numCreatedNodes() == 4,
            "native partial-owner child update should create b, preserved level-5 parts, and level-6 child");
    require(result.releasedEnvelopeSubtreeNodes == 1,
            "native partial-owner child update should release the old level-6 child");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native partial-owner child result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerNativeRecompositionUpdatesNestedPartialOwnerWithoutMaterialization() {
    auto input = makeImage(1, 5, {0, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2}, 6, 3, true);
    require(!region_debug::RegionSplitUpdater::isMaterializedAsProperParts(maxTree, makePixelSet({1, 2})),
            "native nested partial-owner fixture should start non-materialized");

    const auto result = region_debug::RegionSplitUpdater::updateLowerNativeRecomposition(
        maxTree,
        *adj,
        change);

    require(result.numReplacedEnvelopeNodes() == 1,
            "native nested partial-owner update should replace the old level-5 envelope");
    require(result.numCreatedNodes() == 2,
            "native nested partial-owner update should create b and the preserved level-5 child");
    require(result.releasedEnvelopeSubtreeNodes == 2,
            "native nested partial-owner update should release the old level-5/level-6 subtree");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "native nested partial-owner result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerMaterializedBatchesReuseExistingLevelBAncestorForNestedSource() {
    auto input = makeImage(1, 5, {0, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 2}, 6, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerByFlatLevelBatchesWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "nested non-flat update should process the level-5 and level-6 sources");
    require(result.numCreatedTargetNodes() == 1,
            "nested non-flat update should create b once and reuse it for the higher source");
    require(result.numReleasedSourceNodes() == 2,
            "nested non-flat update should release both updated source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "nested non-flat result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testLowerMaterializedBatchesReuseExistingLevelBAncestorForMultipleParents() {
    auto input = makeImage(1, 7, {0, 5, 6, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2, 3, 4}, 6, 3, true);
    const auto result = region_debug::RegionSplitUpdater::updateLowerByFlatLevelBatchesWithMaterialization(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 3,
            "nested multiple-parent update should process one low source and two high sources");
    require(result.numCreatedTargetNodes() == 1,
            "nested multiple-parent update should reuse the b-node created by the low source");
    require(result.numReleasedSourceNodes() == 3,
            "nested multiple-parent update should release all updated source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "nested multiple-parent result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerPreservedAncestorSplitReusesExistingLevelBAncestor() {
    auto input = makeImage(1, 7, {0, 5, 6, 5, 6, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto firstChange = makeLowerChange({1}, 5, 3, true);
    auto firstResult = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        firstChange);
    require(firstResult.numCreatedTargetNodes() == 1,
            "staged preserved-split fixture should create the level-b ancestor first");

    const auto secondChange = makeLowerChange({2, 4}, 6, 3, true);
    const auto secondPlan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        secondChange);
    require(secondPlan.interaction == region_debug::RegionSplitBatchInteraction::PreservedAncestorSplits,
            "staged high update should split the preserved ancestor below an existing b-node");

    const auto secondResult = region_debug::RegionSplitUpdater::updateFlatLowerWithMaterialization(
        maxTree,
        *adj,
        secondChange);
    require(secondResult.numComponents() == 2,
            "staged preserved-split update should process both high sources");
    require(secondResult.numCreatedTargetNodes() == 0,
            "staged preserved-split update should reuse the existing level-b ancestor");
    require(secondResult.numReleasedSourceNodes() == 2,
            "staged preserved-split update should release both high source nodes");
    require(secondResult.reconstruction->isEqual(maxTree.reconstructionImage()),
            "staged preserved-split result image must be the current tree reconstruction");

    const auto combinedChange = makeLowerChange({1, 2, 4}, 6, 3, true);
    requireMatchesRebuildMaxTree(input, adj, combinedChange, maxTree);
}

void testBatchPlanClassifiesSharedExistingLevelBTarget() {
    auto input = makeImage(1, 5, {1, 5, 1, 5, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 1);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 2,
            "batch plan should identify the two disconnected C components");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::SharedExistingLevelBTarget,
            "batch plan should classify reuse of the same existing level-b parent");
    require(plan.supportedByCurrentUpdater(),
            "shared existing level-b target is supported by the current sequential updater");
    for (const auto &component : plan.components) {
        require(component.reusesExistingLevelBTarget,
                "each source should reuse its existing level-b parent");
    }
}

void testBatchPlanClassifiesInsertedLevelBAncestorNeed() {
    auto input = makeImage(1, 5, {0, 5, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 2,
            "batch plan should identify the two disconnected C components");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::RequiresInsertedLevelBAncestor,
            "batch plan should classify the need to insert a shared level-b ancestor");
    require(!plan.supportedByCurrentUpdater(),
            "inserting a shared level-b ancestor is not supported by the current updater");
    for (const auto &component : plan.components) {
        require(component.requiresInsertedLevelBAncestor,
                "each component should be below a preserved parent with altitude greater than b");
    }

    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::updateFlatLowerIndependentSingleOwners(
                maxTree,
                *adj,
                change);
        },
        "batch update must reject cases that need an inserted level-b ancestor");
    require(maxTree.reconstructionImage()->isEqual(input),
            "unsupported inserted-ancestor case must leave the tree unchanged");
}

void testBatchPlanClassifiesMultipleInsertedLevelBAncestors() {
    auto input = makeImage(1, 7, {0, 5, 3, 0, 5, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 4}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 2,
            "batch plan should identify two disconnected C components under different ancestors");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::RequiresMultipleInsertedLevelBAncestors,
            "batch plan should classify insertion under multiple preserved ancestors");
    require(plan.candidateForFlatLowerLocalUpdater(),
            "multiple inserted level-b ancestors should be a candidate for the strict local updater");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "multiple inserted ancestor update should process both disconnected C components");
    require(result.numCreatedTargetNodes() == 2,
            "multiple inserted ancestor update should create one b-node per preserved ancestor");
    require(result.numReleasedSourceNodes() == 2,
            "multiple inserted ancestor update should release both high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "multiple inserted ancestor result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testBatchPlanClassifiesInteractingInsertedLevelBAncestors() {
    auto input = makeImage(1, 7, {0, 5, 4, 5, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3, 5}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 3,
            "batch plan should identify three disconnected C components");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::InteractingInsertedLevelBAncestors,
            "batch plan should classify nested inserted level-b ancestor groups as interacting");
    require(plan.candidateForFlatLowerLocalUpdater(),
            "restricted interacting inserted ancestor groups should be candidates for joint recomposition");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 3,
            "interacting inserted ancestor update should process the three source nodes");
    require(result.numCreatedTargetNodes() == 1,
            "interacting inserted ancestor update should create one shared b-node");
    require(result.numReleasedSourceNodes() == 3,
            "interacting inserted ancestor update should release all high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "interacting inserted ancestor result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testBatchPlanClassifiesMixedInsertedAndIndependentTargets() {
    auto input = makeImage(1, 6, {0, 5, 0, 5, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 2,
            "batch plan should identify two disconnected C components with different target types");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::MixedInsertedAndNonInsertedTargets,
            "batch plan should classify mixed inserted and non-inserted target types");
    require(plan.candidateForFlatLowerLocalUpdater(),
            "mixed inserted and non-inserted targets should be a candidate for the strict local updater");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "mixed target update should process both disconnected C components");
    require(result.numCreatedTargetNodes() == 2,
            "mixed target update should create one b-node for each local target");
    require(result.numReleasedSourceNodes() == 2,
            "mixed target update should release both high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "mixed target result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testBatchPlanClassifiesPreservedAncestorSplit() {
    auto input = makeImage(1, 5, {0, 3, 5, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 1,
            "batch plan should identify the bridge component C");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::PreservedAncestorSplits,
            "batch plan should classify a preserved ancestor split after removing C");
    require(plan.candidateForFlatLowerLocalUpdater(),
            "preserved ancestor split should be a candidate for the strict local split updater");
}

void testFlatLowerPreservedAncestorSplitUpdate() {
    auto input = makeImage(1, 5, {0, 3, 5, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "preserved ancestor split update should process the single bridge component C");
    require(result.numCreatedTargetNodes() == 1,
            "preserved ancestor split update should create one shared b-node");
    require(result.numReleasedSourceNodes() == 1,
            "preserved ancestor split update should release the high bridge source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "preserved ancestor split result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerPreservedAncestorSplitRootUpdate() {
    auto input = makeImage(1, 3, {3, 5, 3});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "root preserved ancestor split update should process the single bridge component C");
    require(result.numCreatedTargetNodes() == 1,
            "root preserved ancestor split update should create one new root at b");
    require(result.numReleasedSourceNodes() == 1,
            "root preserved ancestor split update should release the high bridge source node");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "root preserved ancestor split update should make the b-node the root");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "root preserved ancestor split result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerPreservedAncestorSplitReattachesSourceChildrenUpdate() {
    auto input = makeImage(1, 6, {0, 3, 5, 6, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.interaction == region_debug::RegionSplitBatchInteraction::PreservedAncestorSplits,
            "batch plan should still classify source-child fixtures as preserved ancestor splits");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "preserved ancestor split update should process the bridge source with a child");
    require(result.numCreatedTargetNodes() == 1,
            "preserved ancestor split update should create one b-node when reattaching source children");
    require(result.numReleasedSourceNodes() == 1,
            "preserved ancestor split update should release the source after moving its child");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "preserved ancestor split source-child result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerPreservedAncestorSplitPromotesEffectiveChildUpdate() {
    auto input = makeImage(1, 5, {0, 6, 5, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({2}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.interaction == region_debug::RegionSplitBatchInteraction::PreservedAncestorSplits,
            "batch plan should classify fixtures with effective child promotion as preserved ancestor splits");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "preserved ancestor split update should process the bridge source");
    require(result.numCreatedTargetNodes() == 1,
            "preserved ancestor split update should create one b-node");
    require(result.numReleasedSourceNodes() == 1,
            "preserved ancestor split update should release the source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "preserved ancestor split effective-child result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerPreservedAncestorSplitReattachesUnaffectedChildrenUpdate() {
    auto input = makeImage(1, 7, {0, 3, 6, 3, 5, 3, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({4}, 5, 2);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.interaction == region_debug::RegionSplitBatchInteraction::PreservedAncestorSplits,
            "batch plan should classify fixtures with unaffected children as preserved ancestor splits");

    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 1,
            "preserved ancestor split update should process the bridge source with unaffected siblings");
    require(result.numCreatedTargetNodes() == 1,
            "preserved ancestor split update should create one b-node when reattaching unaffected children");
    require(result.numReleasedSourceNodes() == 1,
            "preserved ancestor split update should release the source node");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "preserved ancestor split unaffected-child result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerSharedInsertedLevelBAncestorUpdate() {
    auto input = makeImage(1, 5, {0, 5, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerSharedInsertedLevelBAncestor(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "shared inserted ancestor update should process both disconnected C components");
    require(result.numCreatedTargetNodes() == 1,
            "shared inserted ancestor update should create a single shared b-node");
    require(result.numReleasedSourceNodes() == 2,
            "shared inserted ancestor update should release both high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "shared inserted ancestor result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerSharedInsertedLevelBRootAncestorUpdate() {
    auto input = makeImage(1, 3, {5, 3, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 2}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerSharedInsertedLevelBAncestor(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "root inserted ancestor update should process both disconnected C components");
    require(result.numCreatedTargetNodes() == 1,
            "root inserted ancestor update should create one new root at b");
    require(result.numReleasedSourceNodes() == 2,
            "root inserted ancestor update should release both high source nodes");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "root inserted ancestor update should make the new b-node the root");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "root inserted ancestor result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerLocalSupportedBatchDispatchesIndependentUpdate() {
    auto input = makeImage(1, 5, {0, 5, 0, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "local supported batch should dispatch independent components");
    require(result.numCreatedTargetNodes() == 2,
            "local supported batch should preserve independent target creation");
    require(result.numReleasedSourceNodes() == 2,
            "local supported batch should release both independent source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "local supported batch independent result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerLocalSupportedBatchDispatchesSharedInsertedAncestorUpdate() {
    auto input = makeImage(1, 5, {0, 5, 3, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "local supported batch should dispatch shared inserted ancestor components");
    require(result.numCreatedTargetNodes() == 1,
            "local supported batch should create a single shared b-node");
    require(result.numReleasedSourceNodes() == 2,
            "local supported batch should release both high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "local supported batch inserted-ancestor result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testBatchPlanClassifiesRepeatedSourceNode() {
    auto input = makeImage(1, 3, {5, 6, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 2}, 5, 3);
    const auto plan = region_debug::RegionSplitUpdater::buildFlatLowerBatchPlan(
        maxTree,
        *adj,
        change);

    require(plan.numComponents() == 2,
            "batch plan should identify two disconnected parts of C");
    require(plan.interaction == region_debug::RegionSplitBatchInteraction::RepeatedSourceNode,
            "batch plan should classify multiple C components with the same direct owner");
    require(!plan.supportedByCurrentUpdater(),
            "repeated source node is not supported by the current updater");

    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::updateFlatLowerIndependentSingleOwners(
                maxTree,
                *adj,
                change);
        },
        "batch update must reject repeated source nodes before mutating");
    require(maxTree.reconstructionImage()->isEqual(input),
            "unsupported repeated-source case must leave the tree unchanged");
}

void testFlatLowerLocalSupportedBatchRejectsRepeatedSourceNode() {
    auto input = makeImage(1, 3, {5, 6, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 2}, 5, 3);
    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::updateFlatLowerLocalSupportedBatch(
                maxTree,
                *adj,
                change);
        },
        "local supported batch must reject repeated source nodes");
    require(maxTree.reconstructionImage()->isEqual(input),
            "local supported batch repeated-source rejection must leave the tree unchanged");
}

void testFlatLowerIndependentSingleOwnersReuseParentUpdate() {
    auto input = makeImage(1, 5, {1, 5, 1, 5, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 1);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerIndependentSingleOwners(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "batch update should process two disconnected single-owner components");
    require(result.numCreatedTargetNodes() == 0,
            "batch update should reuse the common parent when b is already materialized");
    require(result.numReleasedSourceNodes() == 2,
            "batch update should release both high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "batch update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testFlatLowerIndependentSingleOwnersCreateTargetsUpdate() {
    auto input = makeImage(1, 5, {0, 5, 0, 5, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1, 3}, 5, 2);
    const auto result = region_debug::RegionSplitUpdater::updateFlatLowerIndependentSingleOwners(
        maxTree,
        *adj,
        change);

    require(result.numComponents() == 2,
            "batch update should process two disconnected single-owner components");
    require(result.numCreatedTargetNodes() == 2,
            "batch update should create one b-node per disconnected component");
    require(result.numReleasedSourceNodes() == 2,
            "batch update should release both high source nodes");
    require(result.reconstruction->isEqual(maxTree.reconstructionImage()),
            "batch update result image must be the current tree reconstruction");
    requireMatchesRebuildMaxTree(input, adj, change, maxTree);
}

void testBatchRejectsBeforePartialMutation() {
    auto input = makeImage(1, 4, {5, 1, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 2}, 5, 1);
    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::updateFlatLowerIndependentSingleOwners(
                maxTree,
                *adj,
                change);
        },
        "batch update must reject if any component is not a complete direct owner");
    require(maxTree.reconstructionImage()->isEqual(input),
            "batch preflight failure must leave the tree unchanged");
}

void testRejectsNonFlatRegion() {
    auto input = makeImage(1, 2, {5, 6});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({0, 1}, 5, 3, true, true);
    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::buildFlatLowerPlan(maxTree, *adj, change);
        },
        "flat lower split plan must reject non-flat C");
}

void testRejectsNonMaterializedRegion() {
    auto input = makeImage(1, 3, {5, 5, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.0);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const auto change = makeLowerChange({1}, 5, 3);
    requireThrows(
        [&]() {
            (void) region_debug::RegionSplitUpdater::buildFlatLowerPlan(maxTree, *adj, change);
        },
        "flat lower split plan must reject C that cuts direct proper parts");
}

} // namespace

int main() {
    testFlatLowerRemovalPlan();
    testFlatLowerRemovalUpdate();
    testFlatLowerContractionPlan();
    testFlatLowerContractionUpdate();
    testFlatLowerSplitPlan();
    testFlatLowerSplitUpdate();
    testLowerPlanNonFlatThresholds();
    testLowerPlanRejectsWrongMaxA();
    testLowerUpdateByFlatLevelBatches();
    testLowerUpdateByFlatLevelBatchesRejectsBeforePartialMutation();
    testFlatLowerMaterializesPartialRootOwnerAsSource();
    testFlatLowerUpdateAfterMaterializingPartialRootOwner();
    testFlatLowerUpdateAfterMaterializingPartialOwnerWithChild();
    testFlatLowerMaterializesDisconnectedInsideComponentsAsOneSource();
    testFlatLowerUpdateAfterMaterializingDisconnectedInsideComponents();
    testFlatLowerMaterializesMultiplePartialOwnersAsSources();
    testFlatLowerUpdateAfterMaterializingMultiplePartialOwners();
    testFlatLowerUpdateAfterMaterializingMixedPartialAndCoveredOwners();
    testFlatLowerUpdateAfterMaterializingMixedSourcesWithRepeatedOwner();
    testFlatLowerSourceGroupedPreservedAncestorSplitAfterMaterialization();
    testFlatLowerMaterializedSourceGroupedInteractingInsertedAncestorsSplitsInternalChildUpdate();
    testFlatLowerMaterializedSourceGroupedInteractingInsertedAncestorsUpdate();
    testFlatLowerInteractingInsertedAncestorsRecursiveInternalSplitUpdate();
    testLowerUpdateByFlatLevelBatchesWithMaterialization();
    testLowerUpdateByFlatLevelBatchesWithSourceGroupedMaterialization();
    testLowerNativeRecompositionPlanForNestedNonFlatSource();
    testLowerNativeRecompositionPlanForRecursiveSplitFixture();
    testLowerNativeRecompositionUpdateForNestedNonFlatSource();
    testLowerNativeRecompositionUpdateForRecursiveSplitFixture();
    testLowerNativeRecompositionPreservesPixelsAboveA();
    testLowerNativeRecompositionReusesRootEnvelope();
    testLowerNativeRecompositionReusesRootEnvelopeWithPreservedPixelAboveA();
    testLowerNativeRecompositionHandlesMultipleEnvelopes();
    testLowerNativeRecompositionUpdatesPartialRootOwnerWithoutMaterialization();
    testLowerNativeRecompositionUpdatesPartialOwnerWithPreservedChildWithoutMaterialization();
    testLowerNativeRecompositionUpdatesNestedPartialOwnerWithoutMaterialization();
    testLowerMaterializedBatchesReuseExistingLevelBAncestorForNestedSource();
    testLowerMaterializedBatchesReuseExistingLevelBAncestorForMultipleParents();
    testFlatLowerPreservedAncestorSplitReusesExistingLevelBAncestor();
    testBatchPlanClassifiesSharedExistingLevelBTarget();
    testBatchPlanClassifiesInsertedLevelBAncestorNeed();
    testBatchPlanClassifiesMultipleInsertedLevelBAncestors();
    testBatchPlanClassifiesInteractingInsertedLevelBAncestors();
    testBatchPlanClassifiesMixedInsertedAndIndependentTargets();
    testBatchPlanClassifiesPreservedAncestorSplit();
    testFlatLowerPreservedAncestorSplitUpdate();
    testFlatLowerPreservedAncestorSplitRootUpdate();
    testFlatLowerPreservedAncestorSplitReattachesSourceChildrenUpdate();
    testFlatLowerPreservedAncestorSplitPromotesEffectiveChildUpdate();
    testFlatLowerPreservedAncestorSplitReattachesUnaffectedChildrenUpdate();
    testFlatLowerSharedInsertedLevelBAncestorUpdate();
    testFlatLowerSharedInsertedLevelBRootAncestorUpdate();
    testFlatLowerLocalSupportedBatchDispatchesIndependentUpdate();
    testFlatLowerLocalSupportedBatchDispatchesSharedInsertedAncestorUpdate();
    testBatchPlanClassifiesRepeatedSourceNode();
    testFlatLowerLocalSupportedBatchRejectsRepeatedSourceNode();
    testFlatLowerIndependentSingleOwnersReuseParentUpdate();
    testFlatLowerIndependentSingleOwnersCreateTargetsUpdate();
    testBatchRejectsBeforePartialMutation();
    testRejectsNonFlatRegion();
    testRejectsNonMaterializedRegion();
    return 0;
}
