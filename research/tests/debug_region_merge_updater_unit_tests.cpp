#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/AttributeComputer.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DualMinMaxTreeIncrementalFilter.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseAreaAttributeUpdater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseEligibility.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Algorithm2FlatRegionMergeUpdater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/AttributeDirtyRecomputer.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChangeRebuildUpdater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionMergeUpdater.hpp"

namespace problem2_debug = debug::problem2;
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
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

NodeId firstNonRootNode(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isRoot(nodeId)) {
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

struct RegionCandidate {
    std::vector<PixelId> pixelsC;
    AltitudeType a = 0;
    AltitudeType b = 0;
};

RegionCandidate findOnePartialPlusFullOwnerRegion(const DynamicComponentTree &tree,
                                                  AdjacencyRelation &adj) {
    for (NodeId partialOwner : tree.getNodeSubtree(tree.getRoot())) {
        if (tree.getAltitude(partialOwner) <= 0 || tree.getNumProperParts(partialOwner) < 2) {
            continue;
        }

        const auto partialOwnerPixels = collectProperParts(tree, partialOwner);
        for (PixelId partialPixel : partialOwnerPixels) {
            for (PixelId neighbor : adj.getNeighborPixels(partialPixel)) {
                const NodeId fullOwner = tree.getSmallestComponent(neighbor);
                if (fullOwner == InvalidNode ||
                    fullOwner == partialOwner ||
                    !tree.isAlive(fullOwner) ||
                    tree.getAltitude(fullOwner) <= 0 ||
                    tree.getNumProperParts(fullOwner) == 0) {
                    continue;
                }

                std::vector<PixelId> pixelsC = {partialPixel};
                const auto fullOwnerPixels = collectProperParts(tree, fullOwner);
                pixelsC.insert(pixelsC.end(), fullOwnerPixels.begin(), fullOwnerPixels.end());
                std::sort(pixelsC.begin(), pixelsC.end());
                pixelsC.erase(std::unique(pixelsC.begin(), pixelsC.end()), pixelsC.end());

                if (!region_debug::isConnectedRegion(adj, std::span<const PixelId>(pixelsC))) {
                    continue;
                }

                auto plan = region_debug::analyzeRegionMaterialization(
                    tree,
                    adj,
                    std::span<const PixelId>(pixelsC));
                if (plan.numPartialOwners() != 1 || plan.numFullyCoveredOwners() < 1) {
                    continue;
                }

                int minAltitude = 255;
                int maxAltitude = 0;
                for (const auto &owner : plan.affectedOwners) {
                    minAltitude = std::min(minAltitude, owner.altitude);
                    maxAltitude = std::max(maxAltitude, owner.altitude);
                }
                if (minAltitude <= 0 || maxAltitude <= 0) {
                    continue;
                }

                RegionCandidate candidate;
                candidate.pixelsC = std::move(pixelsC);
                candidate.a = static_cast<AltitudeType>(maxAltitude);
                candidate.b = static_cast<AltitudeType>(minAltitude - 1);
                return candidate;
            }
        }
    }

    throw std::runtime_error("could not find a connected region with one partial owner and one fully covered owner");
}

RegionCandidate findMultiplePartialOwnerRegion(const DynamicComponentTree &tree,
                                               AdjacencyRelation &adj) {
    for (NodeId firstOwner : tree.getNodeSubtree(tree.getRoot())) {
        if (tree.getAltitude(firstOwner) <= 0 || tree.getNumProperParts(firstOwner) < 2) {
            continue;
        }

        const auto firstOwnerPixels = collectProperParts(tree, firstOwner);
        for (PixelId firstPixel : firstOwnerPixels) {
            for (PixelId neighbor : adj.getNeighborPixels(firstPixel)) {
                const NodeId secondOwner = tree.getSmallestComponent(neighbor);
                if (secondOwner == InvalidNode ||
                    secondOwner == firstOwner ||
                    !tree.isAlive(secondOwner) ||
                    tree.getAltitude(secondOwner) <= 0 ||
                    tree.getNumProperParts(secondOwner) < 2) {
                    continue;
                }

                std::vector<PixelId> pixelsC = {firstPixel, neighbor};
                std::sort(pixelsC.begin(), pixelsC.end());
                pixelsC.erase(std::unique(pixelsC.begin(), pixelsC.end()), pixelsC.end());
                if (pixelsC.size() != 2 ||
                    !region_debug::isConnectedRegion(adj, std::span<const PixelId>(pixelsC))) {
                    continue;
                }

                auto plan = region_debug::analyzeRegionMaterialization(
                    tree,
                    adj,
                    std::span<const PixelId>(pixelsC));
                if (plan.numPartialOwners() < 2) {
                    continue;
                }

                int minAltitude = 255;
                int maxAltitude = 0;
                for (const auto &owner : plan.affectedOwners) {
                    minAltitude = std::min(minAltitude, owner.altitude);
                    maxAltitude = std::max(maxAltitude, owner.altitude);
                }
                if (minAltitude <= 0 || maxAltitude <= 0) {
                    continue;
                }

                RegionCandidate candidate;
                candidate.pixelsC = std::move(pixelsC);
                candidate.a = static_cast<AltitudeType>(maxAltitude);
                candidate.b = static_cast<AltitudeType>(minAltitude - 1);
                return candidate;
            }
        }
    }

    throw std::runtime_error("could not find a connected region with multiple partial owners");
}

RegionCandidate findPartialOwnerWithMultipleInsideComponentsRegion(const DynamicComponentTree &tree,
                                                                   AdjacencyRelation &adj) {
    for (NodeId partialOwner : tree.getNodeSubtree(tree.getRoot())) {
        if (tree.getAltitude(partialOwner) <= 0 || tree.getNumProperParts(partialOwner) < 3) {
            continue;
        }

        const auto ownerPixels = collectProperParts(tree, partialOwner);
        for (size_t i = 0; i < ownerPixels.size(); ++i) {
            for (size_t j = i + 1; j < ownerPixels.size(); ++j) {
                const std::vector<PixelId> selectedPartialPixels = {ownerPixels[i], ownerPixels[j]};
                const auto selectedComponents = region_debug::collectLocalConnectedComponents(
                    adj,
                    std::span<const PixelId>(selectedPartialPixels),
                    tree.getNumTotalProperParts());
                if (selectedComponents.size() < 2) {
                    continue;
                }

                for (PixelId selectedPixel : selectedPartialPixels) {
                    for (PixelId neighbor : adj.getNeighborPixels(selectedPixel)) {
                        const NodeId fullOwner = tree.getSmallestComponent(neighbor);
                        if (fullOwner == InvalidNode ||
                            fullOwner == partialOwner ||
                            !tree.isAlive(fullOwner) ||
                            tree.getAltitude(fullOwner) <= 0 ||
                            tree.getNumProperParts(fullOwner) == 0) {
                            continue;
                        }

                        std::vector<PixelId> pixelsC = selectedPartialPixels;
                        const auto fullOwnerPixels = collectProperParts(tree, fullOwner);
                        pixelsC.insert(pixelsC.end(), fullOwnerPixels.begin(), fullOwnerPixels.end());
                        std::sort(pixelsC.begin(), pixelsC.end());
                        pixelsC.erase(std::unique(pixelsC.begin(), pixelsC.end()), pixelsC.end());

                        if (!region_debug::isConnectedRegion(adj, std::span<const PixelId>(pixelsC))) {
                            continue;
                        }

                        auto plan = region_debug::analyzeRegionMaterialization(
                            tree,
                            adj,
                            std::span<const PixelId>(pixelsC));
                        bool hasTargetOwner = false;
                        int minAltitude = 255;
                        int maxAltitude = 0;
                        for (const auto &owner : plan.affectedOwners) {
                            minAltitude = std::min(minAltitude, owner.altitude);
                            maxAltitude = std::max(maxAltitude, owner.altitude);
                            if (owner.nodeId == partialOwner && owner.insideComponents.size() > 1) {
                                hasTargetOwner = true;
                            }
                        }
                        if (!hasTargetOwner ||
                            plan.numFullyCoveredOwners() < 1 ||
                            minAltitude <= 0 ||
                            maxAltitude <= 0) {
                            continue;
                        }

                        RegionCandidate candidate;
                        candidate.pixelsC = std::move(pixelsC);
                        candidate.a = static_cast<AltitudeType>(maxAltitude);
                        candidate.b = static_cast<AltitudeType>(minAltitude - 1);
                        return candidate;
                    }
                }
            }
        }
    }

    throw std::runtime_error("could not find a region with multiple inside components in one partial owner");
}

void requireThrowsUnmaterializedRegion(DynamicComponentTree &minTree,
                                       DynamicComponentTree &maxTree,
                                       const AdjacencyRelationPtr &adj,
                                       const region_debug::RegionLevelChange<AltitudeType> &change) {
    try {
        (void) region_debug::RegionMergeUpdater::updateDualTreeOnly(minTree, maxTree, adj, change);
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error("unmaterialized region update should throw");
}

ImageUInt8Ptr makeImage(int rows, int cols, const std::vector<uint8_t> &raw) {
    auto image = ImageUInt8::create(rows, cols);
    require(static_cast<int>(raw.size()) == image->getSize(), "raw image size must match dimensions");
    for (int i = 0; i < image->getSize(); ++i) {
        (*image)[i] = raw[static_cast<size_t>(i)];
    }
    return image;
}

void requireNoEmptyLiveNodes(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        require(tree.getNumProperParts(nodeId) > 0,
                "local merge update must not leave live nodes without direct proper parts");
    }
}

void requireStrictMaxTreeAltitudeOrder(const DynamicComponentTree &tree) {
    for (NodeId parentId : tree.getNodeSubtree(tree.getRoot())) {
        for (NodeId childId : tree.getChildren(parentId)) {
            require(tree.getAltitude(parentId) < tree.getAltitude(childId),
                    "max-tree parent altitude must be strictly lower than child altitude");
        }
    }
}

bool containsNode(const std::vector<NodeId> &nodes, NodeId nodeId) {
    return std::find(nodes.begin(), nodes.end(), nodeId) != nodes.end();
}

bool containsPixel(const std::vector<PixelId> &pixels, PixelId pixelId) {
    return std::find(pixels.begin(), pixels.end(), pixelId) != pixels.end();
}

void requireAreaBufferMatchesFullRecompute(DynamicComponentTree &tree,
                                           const std::vector<float> &areaBuffer,
                                           const std::string &context) {
    DynamicAreaComputer areaComputer(&tree);
    const std::vector<float> expected = areaComputer.compute();
    require(areaBuffer.size() >= expected.size(),
            context + ": area buffer must cover every node slot");

    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        const float actualArea = areaBuffer[static_cast<size_t>(nodeId)];
        const float expectedArea = expected[static_cast<size_t>(nodeId)];
        require(std::fabs(actualArea - expectedArea) < 0.5f,
                context + ": area mismatch at node " + std::to_string(nodeId));
    }
}

std::vector<PixelId> collectAdjacentInteriorPixels(int sideLength,
                                                   const std::vector<PixelId> &pixelsP) {
    std::vector<PixelId> adjacent;
    for (PixelId pixelId : pixelsP) {
        const int row = pixelId / sideLength;
        const int col = pixelId % sideLength;
        for (int dr = -1; dr <= 1; ++dr) {
            for (int dc = -1; dc <= 1; ++dc) {
                if (dr == 0 && dc == 0) {
                    continue;
                }
                const int nr = row + dr;
                const int nc = col + dc;
                if (nr <= 0 || nr >= sideLength - 1 ||
                    nc <= 0 || nc >= sideLength - 1) {
                    continue;
                }
                const PixelId neighbor = nr * sideLength + nc;
                if (!containsPixel(pixelsP, neighbor) &&
                    !containsPixel(adjacent, neighbor)) {
                    adjacent.push_back(neighbor);
                }
            }
        }
    }
    std::sort(adjacent.begin(), adjacent.end());
    return adjacent;
}

const problem2_debug::Algorithm2XiBulletLevel *findXiBulletLevel(
    const problem2_debug::Algorithm2FlatRaiseApplyResult &result,
    AltitudeType level) {
    for (const auto &xiLevel : result.xiBulletByLevel) {
        if (xiLevel.level == level) {
            return &xiLevel;
        }
    }
    return nullptr;
}

void requireCreatedChainEdge(const DynamicComponentTree &tree,
                             const problem2_debug::Algorithm2FlatRaiseApplyResult &result,
                             AltitudeType parentLevel,
                             AltitudeType childLevel,
                             const std::string &message) {
    const NodeId parent = result.createdNodeAtLevel(parentLevel);
    const NodeId child = result.createdNodeAtLevel(childLevel);
    require(parent != InvalidNode && child != InvalidNode, message + ": missing created node");
    require(tree.getNodeParent(child) == parent, message);
}

void requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
    DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    const region_debug::RegionLevelChange<AltitudeType> &change,
    const std::string &message) {
    const auto before = tree.reconstructionImage();
    try {
        (void) problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
            tree,
            adjacency,
            change);
    } catch (const std::runtime_error &) {
        require(tree.reconstructionImage()->isEqual(before),
                message + ": rejected update must preserve reconstruction");
        return;
    }
    throw std::runtime_error(message + ": update should have been rejected");
}

void requireEligibilityRejected(
    const DynamicComponentTree &tree,
    AdjacencyRelation &adjacency,
    const region_debug::RegionLevelChange<AltitudeType> &change,
    const std::string &message) {
    const auto eligibility = problem2_debug::checkAlgorithm2FlatRaiseEligibility(
        tree,
        adjacency,
        change);
    require(!eligibility.accepted, message + ": eligibility should reject");
    require(!eligibility.reason.empty(), message + ": rejection should report a reason");
}

struct DeterministicRng {
    uint32_t state = 0xC0FFEEu;

    int nextInt(int minValue, int maxValue) {
        state = state * 1664525u + 1013904223u;
        const uint32_t range = static_cast<uint32_t>(maxValue - minValue + 1);
        return minValue + static_cast<int>((state >> 8) % range);
    }
};

void requireAlgorithm2FlatRaiseResultConsistent(
    const DynamicComponentTree &tree,
    const problem2_debug::Algorithm2FlatRaiseApplyResult &result,
    const region_debug::RegionLevelChange<AltitudeType> &change,
    const std::string &context) {
    require(result.nodeQb != InvalidNode &&
            tree.isAlive(result.nodeQb) &&
            tree.getAltitude(result.nodeQb) == change.b,
            context + ": Q_b must be alive at level b");
    require(result.createdNodeAtLevel(change.b) == result.nodeQb,
            context + ": created level-b node must be Q_b");
    require(!result.affectedLevels.empty() &&
            result.affectedLevels.front() == change.b,
            context + ": affected levels must start at b");

    for (size_t i = 0; i < result.affectedLevels.size(); ++i) {
        const AltitudeType level = result.affectedLevels[i];
        require(level > change.a && level <= change.b,
                context + ": affected level must lie in (a,b]");
        if (i > 0) {
            require(result.affectedLevels[i - 1] > level,
                    context + ": affected levels must be strictly decreasing");
        }
        const NodeId created = result.createdNodeAtLevel(level);
        require(created != InvalidNode &&
                tree.isAlive(created) &&
                tree.getAltitude(created) == level,
                context + ": every affected level must have an alive created node");
    }

    require(result.createdNodesByLevel.size() == result.affectedLevels.size(),
            context + ": created-node levels must match affected levels");
    for (const auto &xiLevel : result.xiBulletByLevel) {
        require(xiLevel.level > change.a && xiLevel.level <= change.b,
                context + ": Xi^bullet level must lie in (a,b]");
        require(xiLevel.createdNodeQh == result.createdNodeAtLevel(xiLevel.level),
                context + ": Xi^bullet level must point to the created node at the same level");
        for (NodeId nodeId : xiLevel.nodes) {
            require(containsNode(result.xiBulletNodes, nodeId),
                    context + ": Xi^bullet per-level node must also be in flat Xi^bullet list");
            require(containsNode(result.absorbedNodes, nodeId),
                    context + ": Xi^bullet node must be recorded as absorbed");
        }
    }

    for (NodeId nodeId : result.xiCircleNodes) {
        require(nodeId != InvalidNode &&
                tree.isAlive(nodeId) &&
                tree.getAltitude(nodeId) > change.b,
                context + ": Xi^circle node must remain alive above b");
        require(tree.getNodeParent(nodeId) == result.nodeQb,
                context + ": Xi^circle node must be attached under Q_b");
    }

    require(result.releasedSourceNodeQa != result.preservedSourceNodeQa,
            context + ": Q_a must be either released or preserved");
    if (result.preservedSourceNodeQa) {
        require(tree.isAlive(result.sourceNodeQa),
                context + ": preserved Q_a must remain alive");
        require(tree.getNodeParent(result.recompositionRoot) == result.sourceNodeQa,
                context + ": recomposition root must be attached under preserved Q_a");
    }

    for (PixelId pixelId : change.pixelsC) {
        require(tree.getSmallestComponent(pixelId) == result.nodeQb,
                context + ": every pixel of P must move to Q_b");
    }
}

void testSubtreeDescriptorMatchesProductionDualUpdate() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree productionMaxTree(input->clone(), true, adj);
    DynamicComponentTree productionMinTree(input->clone(), false, adj);
    DynamicComponentTree debugMaxTree(input->clone(), true, adj);
    DynamicComponentTree debugMinTree(input->clone(), false, adj);

    const NodeId subtreeRoot = firstNonRootNode(productionMaxTree);
    require(subtreeRoot != InvalidNode, "test image must have a non-root subtree");

    auto anchored = region_debug::buildSubtreeRegionChange<AltitudeType>(
        debugMaxTree,
        debugMinTree,
        subtreeRoot);

    DualMinMaxTreeIncrementalFilter<AltitudeType> productionAdjust(
        &productionMinTree,
        &productionMaxTree,
        *adj);
    productionAdjust.updateTree(&productionMinTree, subtreeRoot);

    auto debugImage = region_debug::RegionMergeUpdater::updateDualTreeOnly(
        debugMinTree,
        debugMaxTree,
        adj,
        anchored.change);

    require(debugImage->isEqual(productionMinTree.reconstructionImage()),
            "explicit-region dual update must match production subtree dual update");
}

void testSubtreeDescriptorMatchesRebuildOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const NodeId subtreeRoot = firstNonRootNode(maxTree);
    require(subtreeRoot != InvalidNode, "test image must have a non-root subtree");

    auto anchored = region_debug::buildSubtreeRegionChange<AltitudeType>(maxTree, minTree, subtreeRoot);
    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, anchored.change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnly(minTree, maxTree, adj, anchored.change);

    require(updated->isEqual(expected), "explicit-region dual update must match rebuild oracle for subtree C");
}

void testRejectsRegionThatIsNotMarkedMaterialized() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(6),
        static_cast<AltitudeType>(4),
        true,
        false);

    requireThrowsUnmaterializedRegion(minTree, maxTree, adj, change);
}

void testExplicitMaterializedLowerRegionMatchesRebuildOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsC = {13, 14};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(6),
        static_cast<AltitudeType>(4),
        true,
        true);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnly(minTree, maxTree, adj, change);

    require(updated->isEqual(expected),
            "explicit materialized lower region must match rebuild oracle");
}

void testExplicitMaterializedRaiseRegionUsesAlgorithm2AndMatchesRebuildOracle() {
    auto input = makeImage(1, 3, {5, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsC = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        true);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnly(minTree, maxTree, adj, change);

    require(updated->isEqual(expected),
            "Algorithm 2 flat raise path must match rebuild oracle");
    require(maxTree.reconstructionImage()->isEqual(expected),
            "updated max-tree reconstruction must match rebuild oracle");
    require(maxTree.getAltitude(maxTree.getRoot()) == change.b,
            "when rho(Q_a)=P, Q_b should become the max-tree root in this fixture");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaiseRejectsMinTreeBeforeMutation() {
    auto input = makeImage(1, 3, {1, 1, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree minTree(input->clone(), false, adj);
    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
        minTree,
        *adj,
        change,
        "Algorithm 2 flat raise must reject min-tree input");
}

void testAlgorithm2FlatRaiseEligibilityAcceptsPreserveQa() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 3, 1, 0,
        0, 2, 1, 6, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const std::vector<PixelId> pixelsP = {12};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        false);

    const auto eligibility = problem2_debug::checkAlgorithm2FlatRaiseEligibility(
        maxTree,
        *adj,
        change);
    require(eligibility.accepted,
            "Algorithm 2 eligibility should accept preserve-Q_a fixture");
    require(eligibility.sourceNodeQa != InvalidNode,
            "Algorithm 2 eligibility should report Q_a");
    require(eligibility.a == 1 && eligibility.b == 4,
            "Algorithm 2 eligibility should report levels a and b");
    require(eligibility.pixelsP == 1,
            "Algorithm 2 eligibility should report |P|");
    require(!eligibility.releasesQa,
            "Algorithm 2 eligibility should classify preserved Q_a");
    require(eligibility.higherBoundaryNeighbors > 0,
            "Algorithm 2 eligibility should count higher boundary neighbors");

    const auto inferred = problem2_debug::checkAlgorithm2FlatRaiseEligibility(
        maxTree,
        *adj,
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(4),
        true);
    require(inferred.accepted &&
            inferred.sourceNodeQa == eligibility.sourceNodeQa &&
            inferred.a == eligibility.a,
            "Algorithm 2 inferred-level eligibility should match explicit change");
}

void testAlgorithm2FlatRaiseEligibilityAcceptsReleaseQa() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 1, 0, 0,
        0, 1, 1, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const std::vector<PixelId> pixelsP = {6, 7, 11, 12};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    const auto eligibility = problem2_debug::checkAlgorithm2FlatRaiseEligibility(
        maxTree,
        *adj,
        change);
    require(eligibility.accepted,
            "Algorithm 2 eligibility should accept release-Q_a fixture");
    require(eligibility.pixelsP == 4 &&
            eligibility.directProperPartsOfQa == 4,
            "Algorithm 2 eligibility should report full Q_a direct support");
    require(eligibility.releasesQa,
            "Algorithm 2 eligibility should classify released Q_a");
}

void testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts() {
    {
        auto input = makeImage(1, 3, {1, 1, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree minTree(input->clone(), false, adj);
        const std::vector<PixelId> pixelsP = {1};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            true,
            false);
        requireEligibilityRejected(minTree, *adj, change, "eligibility must reject min-tree");
    }
    {
        auto input = makeImage(1, 3, {1, 1, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP;
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            true,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject empty P");
    }
    {
        auto input = makeImage(1, 3, {1, 1, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP = {3};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            true,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject out-of-bounds P");
    }
    {
        auto input = makeImage(1, 3, {1, 1, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP = {1, 1};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            true,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject duplicated pixels in P");
    }
    {
        auto input = makeImage(1, 3, {1, 1, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP = {1};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(1),
            true,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject b <= a");
    }
    {
        auto input = makeImage(1, 5, {1, 0, 1, 0, 0});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP = {0, 2};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            false,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject disconnected P");
    }
    {
        auto input = makeImage(1, 3, {1, 2, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP = {1};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            true,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject level mismatch");
    }
    {
        auto input = makeImage(1, 3, {1, 0, 1});
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);
        const std::vector<PixelId> pixelsP = {0, 2};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            static_cast<AltitudeType>(1),
            static_cast<AltitudeType>(3),
            true,
            false);
        requireEligibilityRejected(maxTree, *adj, change, "eligibility must reject multiple Q_a owners");
    }
}

void testAlgorithm2FlatRaiseRejectsEmptyRegionBeforeMutation() {
    auto input = makeImage(1, 3, {1, 1, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP;
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
        maxTree,
        *adj,
        change,
        "Algorithm 2 flat raise must reject empty P");
}

void testAlgorithm2FlatRaiseRejectsNonRaiseBeforeMutation() {
    auto input = makeImage(1, 3, {1, 1, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(1),
        true,
        false);

    requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
        maxTree,
        *adj,
        change,
        "Algorithm 2 flat raise must reject b <= a");
}

void testAlgorithm2FlatRaiseRejectsDisconnectedRegionBeforeMutation() {
    auto input = makeImage(1, 5, {1, 0, 1, 0, 0});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {0, 2};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        false,
        false);

    requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
        maxTree,
        *adj,
        change,
        "Algorithm 2 flat raise must reject disconnected P");
}

void testAlgorithm2FlatRaiseRejectsLevelMismatchBeforeMutation() {
    auto input = makeImage(1, 3, {1, 2, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
        maxTree,
        *adj,
        change,
        "Algorithm 2 flat raise must reject pixels whose current level differs from a");
}

void testAlgorithm2FlatRaiseRejectsMultipleSourceNodesBeforeMutation() {
    auto input = makeImage(1, 3, {1, 0, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {0, 2};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    requireAlgorithm2FlatRaiseThrowsAndPreservesTree(
        maxTree,
        *adj,
        change,
        "Algorithm 2 flat raise must reject P distributed over multiple Q_a nodes");
}

void testAlgorithm2FlatRaiseCreatesIntermediateLevelNode() {
    auto input = makeImage(1, 3, {3, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        true);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwnerMaterialized(
        maxTree,
        *adj,
        change);

    require(result.reconstruction->isEqual(expected),
            "Algorithm 2 flat raise with an affected lower neighbor must match rebuild oracle");
    require(result.numCreatedNodes() == 2,
            "fixture should create Q_b and one lower recomposition node");
    require(result.numAbsorbedNodes() == 1,
            "fixture should absorb the pre-existing level-3 adjacent node");
    require(result.rhoQaEqualsPBeforeUpdate,
            "fixture should exercise the rho(Q_a)=P branch");
    require(result.affectedLevels == std::vector<AltitudeType>({4, 3}),
            "affected levels should be recorded in decreasing max-tree order");
    require(result.createdNodeAtLevel(4) == result.nodeQb,
            "Q_b should be recorded as the created node at level b");
    require(result.createdNodeAtLevel(3) == result.recompositionRoot,
            "lower created node should be the recomposition root");
    const auto *xi3 = findXiBulletLevel(result, static_cast<AltitudeType>(3));
    require(xi3 != nullptr,
            "Xi^bullet should expose the affected level-3 bucket");
    require(xi3->createdNodeQh == result.recompositionRoot,
            "Xi^bullet level-3 bucket should point to the new level-3 node");
    require(xi3->nodes.size() == 1 &&
            containsNode(result.absorbedNodes, xi3->nodes.front()),
            "level-3 Xi^bullet node should be absorbed by the recomposition node");
    require(maxTree.getAltitude(maxTree.getRoot()) == 3,
            "absorbed level-3 node should become the recomposition root");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaisePreservesQaWhenRhoQaDiffersFromP() {
    auto input = makeImage(1, 3, {1, 1, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    require(result.reconstruction->isEqual(expected),
            "Algorithm 2 flat raise with rho(Q_a) != P must match rebuild oracle");
    require(maxTree.reconstructionImage()->isEqual(expected),
            "partial-Q_a raise max-tree reconstruction must match rebuild oracle");
    require(!result.materializedBeforeUpdate,
            "rho(Q_a) != P should be handled directly without local materialization");
    require(!result.rhoQaEqualsPBeforeUpdate,
            "fixture should exercise the rho(Q_a) != P branch");
    require(!result.releasedSourceNodeQa,
            "Q_a must be preserved when it keeps proper parts outside P");
    require(result.preservedSourceNodeQa,
            "result should report that Q_a was preserved");
    require(result.xiBulletNodes.empty(),
            "constant-level fixture should not collect Xi^bullet nodes");
    require(result.xiCircleNodes.empty(),
            "constant-level fixture should not collect Xi^circle frontier nodes");
    require(result.affectedLevels == std::vector<AltitudeType>({3}),
            "only Q_b level should be affected in the constant-level fixture");
    require(result.createdNodeAtLevel(3) == result.nodeQb,
            "created level-b node should be Q_b");
    require(maxTree.getAltitude(maxTree.getRoot()) == 1,
            "outside level-1 pixels should keep Q_a as the root");
    require(maxTree.getNumProperParts(result.sourceNodeQa) == 2,
            "Q_a should keep exactly the two outside proper parts");
    require(maxTree.getNodeParent(result.recompositionRoot) == result.sourceNodeQa,
            "the recomposition root must be attached below the preserved Q_a");
    require(maxTree.getSmallestComponent(1) == result.nodeQb,
            "P should move to Q_b");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaisePreservesQaAndMovesAdjacentFrontierUnderQb() {
    auto input = makeImage(1, 3, {1, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    const NodeId highNode = maxTree.getSmallestComponent(2);
    require(result.reconstruction->isEqual(expected),
            "Algorithm 2 flat raise with a frontier above b must match rebuild oracle");
    require(result.preservedSourceNodeQa,
            "Q_a must survive when outside level-a proper parts remain");
    require(!result.rhoQaEqualsPBeforeUpdate,
            "frontier fixture should preserve Q_a because rho(Q_a) != P");
    require(result.xiBulletNodes.empty(),
            "frontier fixture should not collect merge nodes inside [a+1,b]");
    require(containsNode(result.xiCircleNodes, highNode),
            "Xi^circle should expose the preserved frontier above b");
    require(result.affectedLevels == std::vector<AltitudeType>({3}),
            "frontier fixture should create only Q_b");
    require(highNode != InvalidNode && highNode != result.nodeQb,
            "high neighbor should remain as a distinct frontier node");
    require(maxTree.getNodeParent(highNode) == result.nodeQb,
            "adjacent frontier above b must be reattached under Q_b");
    require(maxTree.getNodeParent(result.nodeQb) == result.sourceNodeQa,
            "Q_b must be attached under preserved Q_a");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaise2DMultiPixelPreservesQaAndMatchesRebuildOracle() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 1, 1, 0,
        0, 1, 1, 1, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {12, 13};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    require(result.reconstruction->isEqual(expected),
            "2D multi-pixel flat raise must match rebuild oracle");
    require(maxTree.reconstructionImage()->isEqual(expected),
            "2D multi-pixel flat raise max-tree reconstruction must match rebuild oracle");
    require(!maxTree.isRoot(result.sourceNodeQa),
            "2D multi-pixel fixture should make Q_a a non-root node");
    require(result.preservedSourceNodeQa,
            "Q_a must be preserved when P is a strict subset of rho(Q_a)");
    require(!result.rhoQaEqualsPBeforeUpdate,
            "2D multi-pixel fixture should exercise rho(Q_a) != P");
    require(maxTree.getNumProperParts(result.nodeQb) == static_cast<int>(pixelsP.size()),
            "Q_b should own all multi-pixel P proper parts");
    require(maxTree.getNumProperParts(result.sourceNodeQa) == 7,
            "preserved Q_a should keep the seven outside level-a proper parts");
    requireAlgorithm2FlatRaiseResultConsistent(maxTree,
                                               result,
                                               change,
                                               "2D multi-pixel preserve fixture");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaise2DReleasesNonRootQaWhenPIsWholeFlatZone() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 1, 0, 0,
        0, 1, 1, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const std::vector<PixelId> pixelsP = {6, 7, 11, 12};
    const NodeId sourceNodeBefore = maxTree.getSmallestComponent(pixelsP.front());
    require(sourceNodeBefore != InvalidNode &&
            !maxTree.isRoot(sourceNodeBefore) &&
            maxTree.getNumProperParts(sourceNodeBefore) == static_cast<int>(pixelsP.size()),
            "2D whole-flatzone fixture must start with non-root Q_a and rho(Q_a)=P");

    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    require(result.reconstruction->isEqual(expected),
            "2D whole-flatzone flat raise must match rebuild oracle");
    require(maxTree.reconstructionImage()->isEqual(expected),
            "2D whole-flatzone flat raise max-tree reconstruction must match rebuild oracle");
    require(result.sourceNodeQa == sourceNodeBefore,
            "result should report the pre-update Q_a");
    require(result.rhoQaEqualsPBeforeUpdate,
            "whole-flatzone fixture should exercise rho(Q_a)=P");
    require(result.releasedSourceNodeQa && !result.preservedSourceNodeQa,
            "Q_a must be released when all its proper parts move to Q_b");
    require(!maxTree.isAlive(sourceNodeBefore),
            "released Q_a should no longer be alive");
    require(result.recompositionRoot == result.nodeQb,
            "without intermediate levels, Q_b should be the recomposition root");
    require(maxTree.getNodeParent(result.nodeQb) == maxTree.getRoot(),
            "released non-root Q_a should be replaced by Q_b below the old parent");
    require(maxTree.getNumProperParts(result.nodeQb) == static_cast<int>(pixelsP.size()),
            "Q_b should own all whole-flatzone proper parts");
    requireAlgorithm2FlatRaiseResultConsistent(maxTree,
                                               result,
                                               change,
                                               "2D whole-flatzone release fixture");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaise2DWithMultipleXiCircleFrontiersAndNonRootQa() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 1, 1, 0,
        0, 5, 1, 6, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const PixelId p = 12;
    const PixelId leftFrontierPixel = 11;
    const PixelId rightFrontierPixel = 13;
    const std::vector<PixelId> pixelsP = {p};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    const NodeId leftFrontier = maxTree.getSmallestComponent(leftFrontierPixel);
    const NodeId rightFrontier = maxTree.getSmallestComponent(rightFrontierPixel);
    require(result.reconstruction->isEqual(expected),
            "2D multi-frontier raise must match rebuild oracle");
    require(!maxTree.isRoot(result.sourceNodeQa),
            "2D fixture should make Q_a a non-root node");
    require(result.preservedSourceNodeQa,
            "Q_a should be preserved because rho(Q_a) != P");
    require(result.xiBulletNodes.empty(),
            "multi-frontier fixture should not contain Xi^bullet nodes");
    require(containsNode(result.xiCircleNodes, leftFrontier) &&
            containsNode(result.xiCircleNodes, rightFrontier),
            "Xi^circle should contain both preserved frontiers above b");
    require(result.affectedLevels == std::vector<AltitudeType>({3}),
            "multi-frontier fixture should create only Q_b");
    require(maxTree.getNodeParent(result.nodeQb) == result.sourceNodeQa,
            "Q_b should be attached under preserved non-root Q_a");
    require(maxTree.getNodeParent(leftFrontier) == result.nodeQb &&
            maxTree.getNodeParent(rightFrontier) == result.nodeQb,
            "both frontiers above b must be reattached under Q_b");
    require(maxTree.getSmallestComponent(p) == result.nodeQb,
            "P should move to Q_b in the 2D multi-frontier fixture");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaise2DWithMultipleXiBulletLevelsAndNonRootQa() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 1, 1, 0,
        0, 2, 1, 3, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const PixelId p = 12;
    const PixelId level2Pixel = 11;
    const PixelId level3Pixel = 13;
    const std::vector<PixelId> pixelsP = {p};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    const NodeId level2Node = maxTree.getSmallestComponent(level2Pixel);
    const NodeId level3Node = maxTree.getSmallestComponent(level3Pixel);
    const auto *xi2 = findXiBulletLevel(result, static_cast<AltitudeType>(2));
    const auto *xi3 = findXiBulletLevel(result, static_cast<AltitudeType>(3));

    require(result.reconstruction->isEqual(expected),
            "2D multi-level Xi^bullet raise must match rebuild oracle");
    require(!maxTree.isRoot(result.sourceNodeQa),
            "2D multi-level fixture should make Q_a a non-root node");
    require(result.preservedSourceNodeQa,
            "Q_a should be preserved in the 2D multi-level fixture");
    require(result.xiCircleNodes.empty(),
            "multi-level fixture should not contain Xi^circle frontiers");
    require(result.affectedLevels == std::vector<AltitudeType>({4, 3, 2}),
            "affected levels should include Q_b and both Xi^bullet levels");
    require(result.createdNodeAtLevel(4) == result.nodeQb,
            "level 4 created node should be Q_b");
    require(result.createdNodeAtLevel(2) == result.recompositionRoot,
            "lowest affected created node should be the recomposition root");
    require(xi2 != nullptr && xi3 != nullptr,
            "Xi^bullet should expose both level-2 and level-3 buckets");
    require(xi2->createdNodeQh == result.createdNodeAtLevel(2) &&
            xi3->createdNodeQh == result.createdNodeAtLevel(3),
            "each Xi^bullet bucket should point to its created replacement node");
    require(xi2->nodes.size() == 1 && xi3->nodes.size() == 1,
            "each 2D Xi^bullet level should absorb one pre-update node");
    require(containsNode(result.absorbedNodes, xi2->nodes.front()) &&
            containsNode(result.absorbedNodes, xi3->nodes.front()),
            "both Xi^bullet nodes should be absorbed");
    requireCreatedChainEdge(maxTree,
                            result,
                            static_cast<AltitudeType>(2),
                            static_cast<AltitudeType>(3),
                            "created level-3 node should be attached under created level-2 node");
    requireCreatedChainEdge(maxTree,
                            result,
                            static_cast<AltitudeType>(3),
                            static_cast<AltitudeType>(4),
                            "Q_b should be attached under created level-3 node");
    require(maxTree.getNodeParent(result.recompositionRoot) == result.sourceNodeQa,
            "created level-2 recomposition root should be attached under preserved Q_a");
    require(maxTree.getAltitude(level2Node) == 2 &&
            maxTree.getAltitude(level3Node) == 3,
            "absorbed pixels should keep their original gray levels in replacement nodes");
    require(maxTree.getSmallestComponent(p) == result.nodeQb,
            "P should move to Q_b in the 2D multi-level fixture");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaise2DWithXiBulletAndXiCircleTogether() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 3, 1, 0,
        0, 2, 1, 6, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    const PixelId p = 12;
    const PixelId level2Pixel = 11;
    const PixelId level3Pixel = 7;
    const PixelId highPixel = 13;
    const NodeId oldLevel2Node = maxTree.getSmallestComponent(level2Pixel);
    const NodeId oldLevel3Node = maxTree.getSmallestComponent(level3Pixel);
    const NodeId highNode = maxTree.getSmallestComponent(highPixel);
    require(oldLevel2Node != InvalidNode &&
            oldLevel3Node != InvalidNode &&
            highNode != InvalidNode,
            "joint Xi fixture must expose all pre-update frontier nodes");

    const std::vector<PixelId> pixelsP = {p};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        false);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    const auto *xi2 = findXiBulletLevel(result, static_cast<AltitudeType>(2));
    const auto *xi3 = findXiBulletLevel(result, static_cast<AltitudeType>(3));
    require(result.reconstruction->isEqual(expected),
            "joint Xi^bullet/Xi^circle fixture must match rebuild oracle");
    require(result.preservedSourceNodeQa,
            "joint Xi fixture should preserve Q_a because rho(Q_a) != P");
    require(xi2 != nullptr && xi3 != nullptr,
            "joint Xi fixture should expose both Xi^bullet levels");
    require(containsNode(result.xiBulletNodes, oldLevel2Node) &&
            containsNode(result.xiBulletNodes, oldLevel3Node),
            "joint Xi fixture should record both old Xi^bullet nodes");
    require(containsNode(result.absorbedNodes, oldLevel2Node) &&
            containsNode(result.absorbedNodes, oldLevel3Node),
            "joint Xi fixture should absorb both Xi^bullet nodes");
    require(containsNode(result.xiCircleNodes, highNode),
            "joint Xi fixture should preserve the high frontier as Xi^circle");
    require(result.affectedLevels == std::vector<AltitudeType>({4, 3, 2}),
            "joint Xi fixture should create Q_b and both lower affected levels");
    require(result.createdNodeAtLevel(4) == result.nodeQb,
            "joint Xi fixture should record Q_b at level 4");
    require(result.createdNodeAtLevel(2) == result.recompositionRoot,
            "joint Xi fixture should use the lowest affected level as recomposition root");
    requireCreatedChainEdge(maxTree,
                            result,
                            static_cast<AltitudeType>(2),
                            static_cast<AltitudeType>(3),
                            "joint Xi fixture should attach level 3 under level 2");
    requireCreatedChainEdge(maxTree,
                            result,
                            static_cast<AltitudeType>(3),
                            static_cast<AltitudeType>(4),
                            "joint Xi fixture should attach Q_b under level 3");
    require(maxTree.getNodeParent(highNode) == result.nodeQb,
            "joint Xi fixture should reattach Xi^circle under Q_b");
    requireAlgorithm2FlatRaiseResultConsistent(maxTree,
                                               result,
                                               change,
                                               "joint Xi^bullet/Xi^circle fixture");
    requireNoEmptyLiveNodes(maxTree);
    requireStrictMaxTreeAltitudeOrder(maxTree);
}

void testAlgorithm2FlatRaiseAreaDirtyClosureMatchesFullRecomputePreserveQa() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 3, 1, 0,
        0, 2, 1, 6, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicAreaComputer areaComputer(&maxTree);
    std::vector<float> areaBuffer = areaComputer.compute();

    const std::vector<PixelId> pixelsP = {12};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        false);

    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);
    const auto areaUpdate = problem2_debug::Algorithm2FlatRaiseAreaAttributeUpdater::recomputeArea(
        maxTree,
        result,
        areaBuffer);

    require(result.preservedSourceNodeQa,
            "area dirty closure preserve fixture should preserve Q_a");
    require(!result.attributeDirtySeeds.empty(),
            "Algorithm 2 should expose attribute dirty seeds in preserve fixture");
    require(areaUpdate.numRecomputedNodes() > 0,
            "area dirty closure preserve fixture should recompute at least one node");
    require(containsNode(areaUpdate.recomputedNodes, result.nodeQb),
            "area dirty closure preserve fixture should recompute Q_b");
    require(containsNode(areaUpdate.dirtyClosure, result.sourceNodeQa),
            "area dirty closure preserve fixture should include preserved Q_a");
    requireAreaBufferMatchesFullRecompute(maxTree,
                                          areaBuffer,
                                          "area dirty closure preserve fixture");
}

void testAlgorithm2FlatRaiseAreaDirtyClosureMatchesFullRecomputeReleaseQa() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 1, 0, 0,
        0, 1, 1, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicAreaComputer areaComputer(&maxTree);
    std::vector<float> areaBuffer = areaComputer.compute();

    const std::vector<PixelId> pixelsP = {6, 7, 11, 12};
    const NodeId sourceNodeBefore = maxTree.getSmallestComponent(pixelsP.front());
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(3),
        true,
        false);

    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);
    const auto areaUpdate = problem2_debug::Algorithm2FlatRaiseAreaAttributeUpdater::recomputeArea(
        maxTree,
        result,
        areaBuffer);

    require(result.releasedSourceNodeQa,
            "area dirty closure release fixture should release Q_a");
    require(!maxTree.isAlive(sourceNodeBefore),
            "area dirty closure release fixture should leave source Q_a dead");
    require(!result.attributeDirtySeeds.empty(),
            "Algorithm 2 should expose attribute dirty seeds in release fixture");
    require(areaUpdate.numRecomputedNodes() > 0,
            "area dirty closure release fixture should recompute at least one node");
    require(containsNode(areaUpdate.recomputedNodes, result.nodeQb),
            "area dirty closure release fixture should recompute Q_b");
    requireAreaBufferMatchesFullRecompute(maxTree,
                                          areaBuffer,
                                          "area dirty closure release fixture");
}

void testAttributeDirtyRecomputerFullTreeFallbackMatchesAreaRecompute() {
    auto input = makeImage(5, 5, {
        0, 0, 0, 0, 0,
        0, 1, 3, 1, 0,
        0, 2, 1, 6, 0,
        0, 1, 1, 1, 0,
        0, 0, 0, 0, 0,
    });
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicAreaComputer areaComputerBefore(&maxTree);
    std::vector<float> areaBuffer = areaComputerBefore.compute();

    const std::vector<PixelId> pixelsP = {12};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        false);

    const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
        maxTree,
        *adj,
        change);

    DynamicAreaComputer areaComputerAfter(&maxTree);
    const auto attributeUpdate = region_debug::AttributeDirtyRecomputer::recomputeFullTree(
        maxTree,
        areaComputerAfter,
        std::span<const NodeId>(result.releasedNodes.data(), result.releasedNodes.size()),
        areaBuffer);

    int liveNodeCount = 0;
    for (NodeId nodeId : maxTree.getNodeSubtree(maxTree.getRoot())) {
        static_cast<void>(nodeId);
        ++liveNodeCount;
    }
    require(attributeUpdate.numRecomputedNodes() == liveNodeCount,
            "full-tree attribute fallback should recompute every live node");
    require(containsNode(attributeUpdate.dirtyClosure, maxTree.getRoot()),
            "full-tree attribute fallback should include the root");
    requireAreaBufferMatchesFullRecompute(maxTree,
                                          areaBuffer,
                                          "full-tree attribute fallback");
}

void testAlgorithm2FlatRaiseRandomSmallFixturesMatchRebuildOracle() {
    DeterministicRng rng;
    const std::vector<PixelId> highNeighborPixels = {7, 11, 13, 17};
    constexpr PixelId centerPixel = 12;

    for (int caseIndex = 0; caseIndex < 96; ++caseIndex) {
        const auto a = static_cast<AltitudeType>(1 + rng.nextInt(0, 2));
        const auto b = static_cast<AltitudeType>(static_cast<int>(a) + 1 + rng.nextInt(0, 6 - static_cast<int>(a) - 1));

        std::vector<uint8_t> raw(25, 0);
        for (int row = 1; row <= 3; ++row) {
            for (int col = 1; col <= 3; ++col) {
                raw[static_cast<size_t>(row * 5 + col)] = static_cast<uint8_t>(a);
            }
        }

        bool hasRaisedNeighbor = false;
        for (PixelId pixelId : highNeighborPixels) {
            if (rng.nextInt(0, 2) != 0) {
                raw[static_cast<size_t>(pixelId)] =
                    static_cast<uint8_t>(static_cast<int>(a) + 1 + rng.nextInt(0, 6 - static_cast<int>(a) - 1));
                hasRaisedNeighbor = true;
            }
        }
        if (!hasRaisedNeighbor) {
            raw[static_cast<size_t>(highNeighborPixels[static_cast<size_t>(caseIndex % static_cast<int>(highNeighborPixels.size()))])] =
                static_cast<uint8_t>(static_cast<int>(a) + 1);
        }

        auto input = makeImage(5, 5, raw);
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);

        const std::vector<PixelId> pixelsP = {centerPixel};
        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            a,
            b,
            true,
            false);

        const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
        const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
            maxTree,
            *adj,
            change);

        const std::string context = "random flat raise case " + std::to_string(caseIndex);
        require(result.reconstruction->isEqual(expected),
                context + ": reconstruction must match rebuild oracle");
        require(maxTree.reconstructionImage()->isEqual(expected),
                context + ": max-tree reconstruction must match rebuild oracle");
        require(!maxTree.isRoot(result.sourceNodeQa),
                context + ": Q_a should be non-root because the fixture has background below a");
        require(result.preservedSourceNodeQa,
                context + ": Q_a should be preserved because outside level-a proper parts remain");
        require(!result.rhoQaEqualsPBeforeUpdate,
                context + ": P is a strict subset of rho(Q_a)");
        requireAlgorithm2FlatRaiseResultConsistent(maxTree, result, change, context);
        requireNoEmptyLiveNodes(maxTree);
        requireStrictMaxTreeAltitudeOrder(maxTree);
    }
}

void testAlgorithm2FlatRaiseRandomPreserveAndReleaseBranchesMatchRebuildOracle() {
    DeterministicRng rng;
    const std::vector<std::vector<PixelId>> preserveShapes = {
        {12, 13},
        {12, 17},
        {11, 12},
        {7, 12},
    };
    const std::vector<std::vector<PixelId>> releaseShapes = {
        {12, 13},
        {12, 17},
        {12, 13, 17, 18},
        {7, 12, 17},
    };

    for (int caseIndex = 0; caseIndex < 64; ++caseIndex) {
        const bool expectRelease = caseIndex % 2 == 0;
        const auto a = static_cast<AltitudeType>(1 + rng.nextInt(0, 2));
        const auto b = static_cast<AltitudeType>(
            static_cast<int>(a) + 1 + rng.nextInt(0, 6 - static_cast<int>(a) - 1));

        std::vector<uint8_t> raw(25, 0);
        std::vector<PixelId> pixelsP;
        if (expectRelease) {
            pixelsP = releaseShapes[static_cast<size_t>(
                rng.nextInt(0, static_cast<int>(releaseShapes.size()) - 1))];
            for (PixelId pixelId : pixelsP) {
                raw[static_cast<size_t>(pixelId)] = static_cast<uint8_t>(a);
            }
        } else {
            for (int row = 1; row <= 3; ++row) {
                for (int col = 1; col <= 3; ++col) {
                    raw[static_cast<size_t>(row * 5 + col)] = static_cast<uint8_t>(a);
                }
            }
            pixelsP = preserveShapes[static_cast<size_t>(
                rng.nextInt(0, static_cast<int>(preserveShapes.size()) - 1))];
        }

        auto input = makeImage(5, 5, raw);
        auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
        DynamicComponentTree maxTree(input->clone(), true, adj);

        const NodeId sourceNodeBefore = maxTree.getSmallestComponent(pixelsP.front());
        const std::string context =
            std::string("random flat raise ") +
            (expectRelease ? "release" : "preserve") +
            " branch case " +
            std::to_string(caseIndex);
        require(sourceNodeBefore != InvalidNode && !maxTree.isRoot(sourceNodeBefore),
                context + ": fixture must start from a non-root Q_a");
        if (expectRelease) {
            require(maxTree.getNumProperParts(sourceNodeBefore) == static_cast<int>(pixelsP.size()),
                    context + ": release fixture must satisfy rho(Q_a)=P before update");
        } else {
            require(maxTree.getNumProperParts(sourceNodeBefore) > static_cast<int>(pixelsP.size()),
                    context + ": preserve fixture must satisfy rho(Q_a)!=P before update");
        }

        auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
            std::span<const PixelId>(pixelsP),
            a,
            b,
            true,
            false);

        const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
        const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
            maxTree,
            *adj,
            change);

        require(result.reconstruction->isEqual(expected),
                context + ": reconstruction must match rebuild oracle");
        require(maxTree.reconstructionImage()->isEqual(expected),
                context + ": max-tree reconstruction must match rebuild oracle");
        require(result.sourceNodeQa == sourceNodeBefore,
                context + ": result must report the original Q_a");
        require(result.rhoQaEqualsPBeforeUpdate == expectRelease,
                context + ": rho(Q_a)=P flag must match the fixture branch");
        require(result.releasedSourceNodeQa == expectRelease,
                context + ": released-source flag must match the fixture branch");
        require(result.preservedSourceNodeQa != expectRelease,
                context + ": preserved-source flag must match the fixture branch");
        if (expectRelease) {
            require(!maxTree.isAlive(sourceNodeBefore),
                    context + ": Q_a should be released");
        } else {
            require(maxTree.isAlive(sourceNodeBefore),
                    context + ": Q_a should be preserved");
        }
        require(maxTree.getNumProperParts(result.nodeQb) == static_cast<int>(pixelsP.size()),
                context + ": Q_b should own every pixel of P as direct proper parts");
        requireAlgorithm2FlatRaiseResultConsistent(maxTree, result, change, context);
        requireNoEmptyLiveNodes(maxTree);
        requireStrictMaxTreeAltitudeOrder(maxTree);
    }
}

void testAlgorithm2FlatRaiseRandomLargerFixturesMatchRebuildOracle() {
    DeterministicRng rng;
    const std::vector<int> sideLengths = {7, 9};

    for (int sideLength : sideLengths) {
        const int center = (sideLength / 2) * sideLength + (sideLength / 2);
        const std::vector<std::vector<PixelId>> localShapes = {
            {center},
            {center, center + 1},
            {center, center + sideLength},
            {center, center + 1, center + sideLength, center + sideLength + 1},
        };

        for (int caseIndex = 0; caseIndex < 40; ++caseIndex) {
            const bool expectRelease = caseIndex % 2 == 0;
            const auto a = static_cast<AltitudeType>(1 + rng.nextInt(0, 2));
            const auto b = static_cast<AltitudeType>(
                static_cast<int>(a) + 1 + rng.nextInt(0, 6 - static_cast<int>(a) - 1));
            std::vector<uint8_t> raw(static_cast<size_t>(sideLength * sideLength), 0);
            const std::vector<PixelId> pixelsP =
                localShapes[static_cast<size_t>(rng.nextInt(0, static_cast<int>(localShapes.size()) - 1))];

            if (expectRelease) {
                for (PixelId pixelId : pixelsP) {
                    raw[static_cast<size_t>(pixelId)] = static_cast<uint8_t>(a);
                }
            } else {
                for (int row = 1; row < sideLength - 1; ++row) {
                    for (int col = 1; col < sideLength - 1; ++col) {
                        raw[static_cast<size_t>(row * sideLength + col)] = static_cast<uint8_t>(a);
                    }
                }
            }

            const auto adjacent = collectAdjacentInteriorPixels(sideLength, pixelsP);
            require(adjacent.size() >= 2,
                    "larger random fixture must provide adjacent pixels for Xi buckets");

            raw[static_cast<size_t>(adjacent[0])] =
                static_cast<uint8_t>(static_cast<int>(a) + 1 + rng.nextInt(0, static_cast<int>(b) - static_cast<int>(a) - 1));
            raw[static_cast<size_t>(adjacent[1])] =
                static_cast<uint8_t>(static_cast<int>(b) + 1);

            for (size_t i = 2; i < adjacent.size(); ++i) {
                if (rng.nextInt(0, 4) == 0) {
                    raw[static_cast<size_t>(adjacent[i])] =
                        static_cast<uint8_t>(static_cast<int>(a) + 1 + rng.nextInt(0, 6 - static_cast<int>(a) - 1));
                }
            }

            auto input = makeImage(sideLength, sideLength, raw);
            auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
            DynamicComponentTree maxTree(input->clone(), true, adj);

            const NodeId sourceNodeBefore = maxTree.getSmallestComponent(pixelsP.front());
            const std::string context =
                "random larger flat raise " +
                std::to_string(sideLength) +
                "x" +
                std::to_string(sideLength) +
                " case " +
                std::to_string(caseIndex);
            require(sourceNodeBefore != InvalidNode && !maxTree.isRoot(sourceNodeBefore),
                    context + ": fixture must start with a non-root Q_a");
            if (expectRelease) {
                require(maxTree.getNumProperParts(sourceNodeBefore) == static_cast<int>(pixelsP.size()),
                        context + ": release fixture must satisfy rho(Q_a)=P before update");
            } else {
                require(maxTree.getNumProperParts(sourceNodeBefore) > static_cast<int>(pixelsP.size()),
                        context + ": preserve fixture must satisfy rho(Q_a)!=P before update");
            }

            auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
                std::span<const PixelId>(pixelsP),
                a,
                b,
                true,
                false);

            const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
            const auto result = problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
                maxTree,
                *adj,
                change);

            require(result.reconstruction->isEqual(expected),
                    context + ": reconstruction must match rebuild oracle");
            require(maxTree.reconstructionImage()->isEqual(expected),
                    context + ": max-tree reconstruction must match rebuild oracle");
            require(result.sourceNodeQa == sourceNodeBefore,
                    context + ": result must report the original Q_a");
            require(result.rhoQaEqualsPBeforeUpdate == expectRelease,
                    context + ": rho(Q_a)=P flag must match the generated branch");
            require(result.releasedSourceNodeQa == expectRelease,
                    context + ": released-source flag must match the generated branch");
            require(result.preservedSourceNodeQa != expectRelease,
                    context + ": preserved-source flag must match the generated branch");
            require(!result.xiBulletNodes.empty(),
                    context + ": generated fixture should include Xi^bullet");
            require(!result.xiCircleNodes.empty(),
                    context + ": generated fixture should include Xi^circle");
            require(maxTree.getNumProperParts(result.nodeQb) >= static_cast<int>(pixelsP.size()),
                    context + ": Q_b should own at least all pixels of P");
            requireAlgorithm2FlatRaiseResultConsistent(maxTree, result, change, context);
            requireNoEmptyLiveNodes(maxTree);
            requireStrictMaxTreeAltitudeOrder(maxTree);
        }
    }
}

void testSinglePartialOwnerMaterializeThenLowerMatchesRebuildOracle() {
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
        false);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingSinglePartialOwner(
        minTree,
        maxTree,
        adj,
        change);

    require(updated->isEqual(expected),
            "single-partial-owner materialize-then-lower path must match rebuild oracle");
    require(minTree.reconstructionImage()->isEqual(expected),
            "updated dual tree reconstruction must match rebuild oracle");
}

void testSinglePartialOwnerPlusFullOwnerMaterializeThenLowerMatchesRebuildOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const auto candidate = findOnePartialPlusFullOwnerRegion(minTree, *adj);
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(candidate.pixelsC),
        candidate.a,
        candidate.b,
        true,
        false);

    const auto prePlan = region_debug::analyzeRegionMaterialization(minTree, *adj, change);
    require(prePlan.numPartialOwners() == 1,
            "test region must have exactly one partial owner before materialization");
    require(prePlan.numFullyCoveredOwners() >= 1,
            "test region must include at least one fully covered owner");

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingSinglePartialOwner(
        minTree,
        maxTree,
        adj,
        change);

    require(updated->isEqual(expected),
            "one-partial-plus-full materialize-then-lower path must match rebuild oracle");
    require(minTree.reconstructionImage()->isEqual(expected),
            "updated dual tree reconstruction must match rebuild oracle");
}

void testMultiplePartialOwnersMaterializeThenLowerMatchesRebuildOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const auto candidate = findMultiplePartialOwnerRegion(minTree, *adj);
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(candidate.pixelsC),
        candidate.a,
        candidate.b,
        true,
        false);

    const auto prePlan = region_debug::analyzeRegionMaterialization(minTree, *adj, change);
    require(prePlan.numPartialOwners() >= 2,
            "test region must have multiple partial owners before materialization");

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingPartialOwners(
        minTree,
        maxTree,
        adj,
        change);

    require(updated->isEqual(expected),
            "multiple-partial-owners materialize-then-lower path must match rebuild oracle");
    require(minTree.reconstructionImage()->isEqual(expected),
            "updated dual tree reconstruction must match rebuild oracle");
}

void testPartialOwnerWithMultipleInsideComponentsMaterializeThenLowerMatchesRebuildOracle() {
    auto input = make_complex_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const auto candidate = findPartialOwnerWithMultipleInsideComponentsRegion(minTree, *adj);
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(candidate.pixelsC),
        candidate.a,
        candidate.b,
        true,
        false);

    const auto prePlan = region_debug::analyzeRegionMaterialization(minTree, *adj, change);
    bool hasMultipleInsideComponents = false;
    for (const auto &owner : prePlan.affectedOwners) {
        if (owner.cutsDirectProperParts() && owner.insideComponents.size() > 1) {
            hasMultipleInsideComponents = true;
        }
    }
    require(hasMultipleInsideComponents,
            "test region must force one partial owner to produce multiple inside components");
    require(prePlan.numFullyCoveredOwners() >= 1,
            "test region must include a fully covered owner connecting the partial components");

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    auto updated = region_debug::RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingPartialOwners(
        minTree,
        maxTree,
        adj,
        change);

    require(updated->isEqual(expected),
            "multiple-inside-components materialize-then-lower path must match rebuild oracle");
    require(minTree.reconstructionImage()->isEqual(expected),
            "updated dual tree reconstruction must match rebuild oracle");
}

void testPairUpdateWithLocalDualAndRebuiltPrimalMatchesRebuildOracle() {
    auto input = make_complex_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree localMaxTree(input->clone(), true, adj);
    DynamicComponentTree localMinTree(input->clone(), false, adj);
    DynamicComponentTree oracleMaxTree(input->clone(), true, adj);
    DynamicComponentTree oracleMinTree(input->clone(), false, adj);

    const auto candidate = findPartialOwnerWithMultipleInsideComponentsRegion(localMinTree, *adj);
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(candidate.pixelsC),
        candidate.a,
        candidate.b,
        true,
        false);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::updatePair(
        oracleMinTree,
        oracleMaxTree,
        adj,
        change);
    auto updated = region_debug::RegionMergeUpdater::updatePairWithLocalDualAndRebuiltPrimal(
        localMinTree,
        localMaxTree,
        adj,
        change,
        AREA);

    require(updated->isEqual(expected),
            "debug pair update must return the rebuild-oracle image");
    require(localMinTree.reconstructionImage()->isEqual(expected),
            "local dual min-tree reconstruction must match rebuild oracle");
    require(localMaxTree.reconstructionImage()->isEqual(expected),
            "rebuilt primal max-tree reconstruction must match rebuild oracle");
    require(localMinTree.reconstructionImage()->isEqual(oracleMinTree.reconstructionImage()),
            "local min-tree reconstruction must match oracle min-tree reconstruction");
    require(localMaxTree.reconstructionImage()->isEqual(oracleMaxTree.reconstructionImage()),
            "local max-tree reconstruction must match oracle max-tree reconstruction");
}

void testPairUpdateWithAlgorithm2RaiseAndRebuiltPrimalMatchesRebuildOracle() {
    auto input = makeImage(1, 3, {3, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree localMaxTree(input->clone(), true, adj);
    DynamicComponentTree localMinTree(input->clone(), false, adj);
    DynamicComponentTree oracleMaxTree(input->clone(), true, adj);
    DynamicComponentTree oracleMinTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        true);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::updatePair(
        oracleMinTree,
        oracleMaxTree,
        adj,
        change);
    auto updated = region_debug::RegionMergeUpdater::updatePairWithLocalDualAndRebuiltPrimal(
        localMinTree,
        localMaxTree,
        adj,
        change);

    require(updated->isEqual(expected),
            "pair update with Algorithm 2 raise must return the rebuild-oracle image");
    require(localMaxTree.reconstructionImage()->isEqual(expected),
            "local max-tree after Algorithm 2 raise must match rebuild oracle");
    require(localMinTree.reconstructionImage()->isEqual(expected),
            "rebuilt min-tree after Algorithm 2 raise must match rebuild oracle");
    require(localMaxTree.reconstructionImage()->isEqual(oracleMaxTree.reconstructionImage()),
            "local max-tree reconstruction must match oracle max-tree reconstruction after raise");
    require(localMinTree.reconstructionImage()->isEqual(oracleMinTree.reconstructionImage()),
            "local min-tree reconstruction must match oracle min-tree reconstruction after raise");
}

void testPairUpdateWithAlgorithm2RaiseRejectsNonAreaAttributeBeforeMutation() {
    auto input = makeImage(1, 3, {3, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree localMaxTree(input->clone(), true, adj);
    DynamicComponentTree localMinTree(input->clone(), false, adj);
    const auto maxBefore = localMaxTree.reconstructionImage();
    const auto minBefore = localMinTree.reconstructionImage();

    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        true);

    try {
        (void) region_debug::RegionMergeUpdater::updatePairWithLocalDualAndRebuiltPrimal(
            localMinTree,
            localMaxTree,
            adj,
            change,
            BOX_WIDTH);
    } catch (const std::runtime_error &) {
        require(localMaxTree.reconstructionImage()->isEqual(maxBefore),
                "non-AREA Algorithm 2 raise rejection must preserve max-tree reconstruction");
        require(localMinTree.reconstructionImage()->isEqual(minBefore),
                "non-AREA Algorithm 2 raise rejection must preserve min-tree reconstruction");
        return;
    }

    throw std::runtime_error("Algorithm 2 raise backend should reject non-AREA attributes");
}

void testAlgorithm2AreaOptInAcceptsAreaAndUpdatesAreaBuffer() {
    auto input = makeImage(1, 3, {3, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree localMaxTree(input->clone(), true, adj);
    DynamicComponentTree localMinTree(input->clone(), false, adj);
    DynamicComponentTree oracleMaxTree(input->clone(), true, adj);
    DynamicComponentTree oracleMinTree(input->clone(), false, adj);

    DynamicAreaComputer maxAreaComputer(&localMaxTree);
    std::vector<float> maxArea = maxAreaComputer.compute();

    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        true);

    auto expected = region_debug::RegionLevelChangeRebuildUpdater::updatePair(
        oracleMinTree,
        oracleMaxTree,
        adj,
        change);
    region_debug::RegionMergeAlgorithm2AreaOptInDiagnostics diagnostics;
    auto updated = region_debug::RegionMergeUpdater::updatePairWithAlgorithm2AreaOptIn(
        localMinTree,
        localMaxTree,
        adj,
        change,
        AREA,
        &maxArea,
        &diagnostics);

    require(updated->isEqual(expected),
            "Algorithm 2 AREA opt-in accepted path must match rebuild oracle");
    require(localMaxTree.reconstructionImage()->isEqual(expected),
            "Algorithm 2 AREA opt-in accepted path must update max-tree");
    require(localMinTree.reconstructionImage()->isEqual(expected),
            "Algorithm 2 AREA opt-in accepted path must rebuild min-tree");
    require(diagnostics.attemptedAlgorithm2 && diagnostics.acceptedAlgorithm2,
            "Algorithm 2 AREA opt-in accepted path must report acceptance");
    require(!diagnostics.fallbackUsed,
            "Algorithm 2 AREA opt-in accepted path must not report fallback");
    require(diagnostics.areaAttributeUpdated,
            "Algorithm 2 AREA opt-in accepted path must update AREA buffer");
    require(diagnostics.areaRecomputedNodes > 0,
            "Algorithm 2 AREA opt-in accepted path must report recomputed AREA nodes");
    requireAreaBufferMatchesFullRecompute(localMaxTree,
                                          maxArea,
                                          "Algorithm 2 AREA opt-in accepted path");
}

void testAlgorithm2AreaOptInFallsBackForNonAreaAttribute() {
    auto input = makeImage(1, 3, {3, 1, 5});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree localMaxTree(input->clone(), true, adj);
    DynamicComponentTree localMinTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsP = {1};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(4),
        true,
        true);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    region_debug::RegionMergeAlgorithm2AreaOptInDiagnostics diagnostics;
    auto updated = region_debug::RegionMergeUpdater::updatePairWithAlgorithm2AreaOptIn(
        localMinTree,
        localMaxTree,
        adj,
        change,
        BOX_WIDTH,
        nullptr,
        &diagnostics);

    require(updated->isEqual(expected),
            "Algorithm 2 AREA opt-in non-AREA path must fallback to rebuild output");
    require(localMaxTree.reconstructionImage()->isEqual(expected),
            "Algorithm 2 AREA opt-in non-AREA fallback must rebuild max-tree");
    require(localMinTree.reconstructionImage()->isEqual(expected),
            "Algorithm 2 AREA opt-in non-AREA fallback must rebuild min-tree");
    require(!diagnostics.attemptedAlgorithm2 &&
            !diagnostics.acceptedAlgorithm2 &&
            diagnostics.rejectedByAttribute &&
            diagnostics.fallbackUsed,
            "Algorithm 2 AREA opt-in non-AREA path must report attribute fallback");
}

void testAlgorithm2AreaOptInFallsBackForContractRejection() {
    auto input = makeImage(1, 3, {1, 0, 1});
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    DynamicComponentTree localMaxTree(input->clone(), true, adj);
    DynamicComponentTree localMinTree(input->clone(), false, adj);

    const std::vector<PixelId> pixelsP = {0, 2};
    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsP),
        static_cast<AltitudeType>(1),
        static_cast<AltitudeType>(2),
        true,
        true);

    const auto expected = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(input, change);
    region_debug::RegionMergeAlgorithm2AreaOptInDiagnostics diagnostics;
    auto updated = region_debug::RegionMergeUpdater::updatePairWithAlgorithm2AreaOptIn(
        localMinTree,
        localMaxTree,
        adj,
        change,
        AREA,
        nullptr,
        &diagnostics);

    require(updated->isEqual(expected),
            "Algorithm 2 AREA opt-in contract rejection must fallback to rebuild output");
    require(localMaxTree.reconstructionImage()->isEqual(expected),
            "Algorithm 2 AREA opt-in contract fallback must rebuild max-tree");
    require(localMinTree.reconstructionImage()->isEqual(expected),
            "Algorithm 2 AREA opt-in contract fallback must rebuild min-tree");
    require(diagnostics.attemptedAlgorithm2 &&
            !diagnostics.acceptedAlgorithm2 &&
            diagnostics.rejectedByContract &&
            diagnostics.fallbackUsed,
            "Algorithm 2 AREA opt-in contract rejection must report fallback");
}

} // namespace

int main() {
    testSubtreeDescriptorMatchesProductionDualUpdate();
    testSubtreeDescriptorMatchesRebuildOracle();
    testRejectsRegionThatIsNotMarkedMaterialized();
    testExplicitMaterializedLowerRegionMatchesRebuildOracle();
    testExplicitMaterializedRaiseRegionUsesAlgorithm2AndMatchesRebuildOracle();
    testAlgorithm2FlatRaiseRejectsMinTreeBeforeMutation();
    testAlgorithm2FlatRaiseEligibilityAcceptsPreserveQa();
    testAlgorithm2FlatRaiseEligibilityAcceptsReleaseQa();
    testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts();
    testAlgorithm2FlatRaiseRejectsEmptyRegionBeforeMutation();
    testAlgorithm2FlatRaiseRejectsNonRaiseBeforeMutation();
    testAlgorithm2FlatRaiseRejectsDisconnectedRegionBeforeMutation();
    testAlgorithm2FlatRaiseRejectsLevelMismatchBeforeMutation();
    testAlgorithm2FlatRaiseRejectsMultipleSourceNodesBeforeMutation();
    testAlgorithm2FlatRaiseCreatesIntermediateLevelNode();
    testAlgorithm2FlatRaisePreservesQaWhenRhoQaDiffersFromP();
    testAlgorithm2FlatRaisePreservesQaAndMovesAdjacentFrontierUnderQb();
    testAlgorithm2FlatRaise2DMultiPixelPreservesQaAndMatchesRebuildOracle();
    testAlgorithm2FlatRaise2DReleasesNonRootQaWhenPIsWholeFlatZone();
    testAlgorithm2FlatRaise2DWithMultipleXiCircleFrontiersAndNonRootQa();
    testAlgorithm2FlatRaise2DWithMultipleXiBulletLevelsAndNonRootQa();
    testAlgorithm2FlatRaise2DWithXiBulletAndXiCircleTogether();
    testAlgorithm2FlatRaiseAreaDirtyClosureMatchesFullRecomputePreserveQa();
    testAlgorithm2FlatRaiseAreaDirtyClosureMatchesFullRecomputeReleaseQa();
    testAttributeDirtyRecomputerFullTreeFallbackMatchesAreaRecompute();
    testAlgorithm2FlatRaiseRandomSmallFixturesMatchRebuildOracle();
    testAlgorithm2FlatRaiseRandomPreserveAndReleaseBranchesMatchRebuildOracle();
    testAlgorithm2FlatRaiseRandomLargerFixturesMatchRebuildOracle();
    testSinglePartialOwnerMaterializeThenLowerMatchesRebuildOracle();
    testSinglePartialOwnerPlusFullOwnerMaterializeThenLowerMatchesRebuildOracle();
    testMultiplePartialOwnersMaterializeThenLowerMatchesRebuildOracle();
    testPartialOwnerWithMultipleInsideComponentsMaterializeThenLowerMatchesRebuildOracle();
    testPairUpdateWithLocalDualAndRebuiltPrimalMatchesRebuildOracle();
    testPairUpdateWithAlgorithm2RaiseAndRebuiltPrimalMatchesRebuildOracle();
    testPairUpdateWithAlgorithm2RaiseRejectsNonAreaAttributeBeforeMutation();
    testAlgorithm2AreaOptInAcceptsAreaAndUpdatesAreaBuffer();
    testAlgorithm2AreaOptInFallsBackForNonAreaAttribute();
    testAlgorithm2AreaOptInFallsBackForContractRejection();
    return 0;
}
