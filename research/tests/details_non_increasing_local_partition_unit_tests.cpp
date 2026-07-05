#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingCasfKernel.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingLocalOracle.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingLocalOracleDependency.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingLocalPartition.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingLocalPartitionPlanner.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingLocalRecomposer.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingStrippedLocalRecomposition.hpp"

namespace noninc = morphotree::details::non_increasing;

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr make_reduced_split_fixture_image() {
    auto image = ImageUInt8::create(6, 6);
    const std::vector<uint8_t> raw = {
        20, 44, 68, 44, 68, 44,
        56, 140, 104, 116, 104, 116,
        68, 140, 104, 116, 140, 128,
        56, 104, 56, 20, 32, 104,
        80, 104, 152, 104, 56, 128,
        68, 116, 140, 68, 20, 92,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_monotone_line_image() {
    auto image = ImageUInt8::create(1, 4);
    (*image)[0] = 10;
    (*image)[1] = 20;
    (*image)[2] = 30;
    (*image)[3] = 40;
    return image;
}

ImageUInt8Ptr make_split_child_line_image() {
    auto image = ImageUInt8::create(1, 4);
    (*image)[0] = 10;
    (*image)[1] = 10;
    (*image)[2] = 30;
    (*image)[3] = 40;
    return image;
}

ImageUInt8Ptr clone_with_values(ImageUInt8Ptr input, const std::vector<uint8_t> &values) {
    auto image = input->clone();
    require(static_cast<int>(values.size()) == image->getSize(),
            "replacement value count must match image size");
    for (std::size_t i = 0; i < values.size(); ++i) {
        (*image)[static_cast<int>(i)] = values[i];
    }
    return image;
}

NodeId findReducedFixtureMergeNode(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (nodeId == tree.getRoot() || !tree.isAlive(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
            continue;
        }
        if (tree.getAltitude(nodeId) == 92 &&
            tree.getAltitude(parentId) == 80 &&
            tree.getNumChildren(nodeId) == 1 &&
            tree.getNumProperParts(nodeId) == 1) {
            return nodeId;
        }
    }
    return InvalidNode;
}

noninc::LocalMergeOracle buildReducedFixtureOracle(AdjacencyRelationPtr adjacency) {
    auto image = make_reduced_split_fixture_image();
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "fixture should not fail before the diagnostic threshold");
    }

    const NodeId nodeId = findReducedFixtureMergeNode(runner.maxTree());
    require(nodeId != InvalidNode, "reduced fixture merge node must be present before threshold 10");
    return noninc::buildLocalMergeOracle(runner.maxTree(), runner.minTree(), adjacency, nodeId);
}

bool componentContainsPixel(const noninc::LocalThresholdComponent &component, PixelId pixelId) {
    return std::find(component.pixels.begin(), component.pixels.end(), pixelId) != component.pixels.end();
}

bool vectorContainsPixel(const std::vector<PixelId> &pixels, PixelId pixelId) {
    return std::find(pixels.begin(), pixels.end(), pixelId) != pixels.end();
}

bool vectorContainsIndex(const std::vector<int> &values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

bool componentsMatchIgnoringOracleLabels(const noninc::LocalThresholdComponent &lhs,
                                         const noninc::LocalThresholdComponent &rhs) {
    return lhs.level == rhs.level &&
           lhs.pixels == rhs.pixels &&
           lhs.touchesModified == rhs.touchesModified &&
           lhs.touchesPreserved == rhs.touchesPreserved &&
           lhs.touchesBoundary == rhs.touchesBoundary &&
           lhs.preDualOwners == rhs.preDualOwners &&
           lhs.preThresholdOwners == rhs.preThresholdOwners;
}

bool partitionPlansMatchIgnoringOracleLabels(const noninc::LocalPartitionPlan &lhs,
                                             const noninc::LocalPartitionPlan &rhs) {
    if (lhs.affectedLevels != rhs.affectedLevels ||
        lhs.levels.size() != rhs.levels.size() ||
        lhs.recomposedComponents.size() != rhs.recomposedComponents.size()) {
        return false;
    }

    for (std::size_t levelIndex = 0; levelIndex < lhs.levels.size(); ++levelIndex) {
        const auto &lhsLevel = lhs.levels[levelIndex];
        const auto &rhsLevel = rhs.levels[levelIndex];
        if (lhsLevel.level != rhsLevel.level ||
            lhsLevel.components.size() != rhsLevel.components.size()) {
            return false;
        }
        for (std::size_t componentIndex = 0; componentIndex < lhsLevel.components.size(); ++componentIndex) {
            if (!componentsMatchIgnoringOracleLabels(
                    lhsLevel.components[componentIndex],
                    rhsLevel.components[componentIndex])) {
                return false;
            }
        }
    }

    for (std::size_t componentIndex = 0;
         componentIndex < lhs.recomposedComponents.size();
         ++componentIndex) {
        if (!componentsMatchIgnoringOracleLabels(
                lhs.recomposedComponents[componentIndex],
                rhs.recomposedComponents[componentIndex])) {
            return false;
        }
    }
    return true;
}

bool recompositionNodesMatchIgnoringOracleLabels(const noninc::LocalRecompositionNodePlan &lhs,
                                                 const noninc::LocalRecompositionNodePlan &rhs) {
    return lhs.level == rhs.level &&
           lhs.supportPixels == rhs.supportPixels &&
           lhs.directProperParts == rhs.directProperParts &&
           lhs.touchesModified == rhs.touchesModified &&
           lhs.touchesPreserved == rhs.touchesPreserved &&
           lhs.touchesBoundary == rhs.touchesBoundary &&
           lhs.preDualOwners == rhs.preDualOwners &&
           lhs.parentPlanIndex == rhs.parentPlanIndex &&
           lhs.childPlanIndices == rhs.childPlanIndices;
}

bool recompositionPlansMatchIgnoringOracleLabels(const noninc::LocalRecompositionPlan &lhs,
                                                 const noninc::LocalRecompositionPlan &rhs) {
    if (lhs.affectedLevels != rhs.affectedLevels ||
        lhs.nodes.size() != rhs.nodes.size()) {
        return false;
    }

    for (std::size_t nodeIndex = 0; nodeIndex < lhs.nodes.size(); ++nodeIndex) {
        if (!recompositionNodesMatchIgnoringOracleLabels(lhs.nodes[nodeIndex], rhs.nodes[nodeIndex])) {
            return false;
        }
    }
    return true;
}

void requireTreeConsistency(const DynamicComponentTree &tree) {
    require(tree.getRoot() != InvalidNode, "root must be valid");
    require(tree.isAlive(tree.getRoot()), "root must be alive");
    require(tree.getNodeParent(tree.getRoot()) == tree.getRoot(), "root parent must be itself");

    std::vector<unsigned char> reachable(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0);
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        reachable[static_cast<std::size_t>(nodeId)] = 1;
    }

    int totalProperParts = 0;
    std::vector<unsigned char> seen(static_cast<std::size_t>(tree.getNumTotalProperParts()), 0);
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }
        require(reachable[static_cast<std::size_t>(nodeId)] != 0,
                "alive node must be reachable from root");
        int countedChildren = 0;
        for (NodeId childId : tree.getChildren(nodeId)) {
            require(tree.isAlive(childId), "child must be alive");
            require(tree.getNodeParent(childId) == nodeId, "child parent must match");
            ++countedChildren;
        }
        require(countedChildren == tree.getNumChildren(nodeId), "child count mismatch");

        int countedProperParts = 0;
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            require(pixelId >= 0 && pixelId < tree.getNumTotalProperParts(),
                    "proper part must be in bounds");
            require(tree.getSmallestComponent(pixelId) == nodeId,
                    "proper part owner mismatch");
            require(seen[static_cast<std::size_t>(pixelId)] == 0,
                    "proper part visited twice");
            seen[static_cast<std::size_t>(pixelId)] = 1;
            ++countedProperParts;
        }
        require(countedProperParts == tree.getNumProperParts(nodeId),
                "proper part count mismatch");
        totalProperParts += countedProperParts;
    }

    require(totalProperParts == tree.getNumTotalProperParts(),
            "global proper part count mismatch");
    require(std::all_of(seen.begin(), seen.end(), [](unsigned char value) {
        return value != 0;
    }), "proper part partition is incomplete");
}

void testAffectedLevelsFollowDualInsertionIntervals() {
    noninc::LocalMergeNodeContractionDescriptor lowerIntoMin;
    lowerIntoMin.a = 92;
    lowerIntoMin.b = 80;
    lowerIntoMin.dualIsMaxtree = false;
    lowerIntoMin.direction = ::debug::region::RegionUpdateDirection::Lower;
    const auto lowerLevels = noninc::affectedLevelsForLocalMerge(lowerIntoMin);
    require(lowerLevels.size() == 12, "lower min-tree interval must contain b..a-1");
    require(lowerLevels.front() == 80 && lowerLevels.back() == 91,
            "lower min-tree interval endpoints must be b and a-1");

    noninc::LocalMergeNodeContractionDescriptor raiseIntoMax;
    raiseIntoMax.a = 4;
    raiseIntoMax.b = 7;
    raiseIntoMax.dualIsMaxtree = true;
    raiseIntoMax.direction = ::debug::region::RegionUpdateDirection::Raise;
    const auto raiseLevels = noninc::affectedLevelsForLocalMerge(raiseIntoMax);
    require((raiseLevels == std::vector<AltitudeType>{5, 6, 7}),
            "raise max-tree interval must contain a+1..b");
}

void testDualThresholdPredicateUsesTreePolarity() {
    require(noninc::pixelPassesDualThreshold(92, 80, true),
            "max-tree threshold must accept pixels above the level");
    require(!noninc::pixelPassesDualThreshold(79, 80, true),
            "max-tree threshold must reject pixels below the level");
    require(noninc::pixelPassesDualThreshold(80, 92, false),
            "min-tree threshold must accept pixels below the level");
    require(!noninc::pixelPassesDualThreshold(93, 92, false),
            "min-tree threshold must reject pixels above the level");
}

void testReducedFixtureBuildsLocalThresholdPartitions() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const auto oracle = buildReducedFixtureOracle(adjacency);
    const auto plan = noninc::buildLocalPartitionPlan(oracle, *adjacency);

    require(plan.affectedLevels.size() == 12, "fixture must affect thresholds 80..91");
    require(plan.affectedLevels.front() == 80 && plan.affectedLevels.back() == 91,
            "fixture affected-threshold endpoints must match J_C");
    require(plan.levels.size() == plan.affectedLevels.size(),
            "local partition must keep one entry per affected threshold");
    require(plan.recomposedComponents.size() == plan.affectedLevels.size(),
            "fixture must expose one component containing C for every affected threshold");

    const PixelId modifiedPixel = oracle.contraction.modifiedProperPart.front();
    bool recomposedTouchesStableContext = false;
    for (const auto &partition : plan.levels) {
        require(!partition.components.empty(),
                "every affected threshold must have at least one local component");

        int modifiedComponentCount = 0;
        for (const auto &component : partition.components) {
            require(!component.pixels.empty(), "local threshold component must not be empty");
            require(!component.oracleDualOwners.empty(),
                    "local threshold component must carry oracle owner labels");
            if (component.touchesModified) {
                ++modifiedComponentCount;
                require(componentContainsPixel(component, modifiedPixel),
                        "the recomposed component must contain the modified pixel");
                recomposedTouchesStableContext =
                    recomposedTouchesStableContext ||
                    component.touchesPreserved ||
                    component.touchesBoundary;
            }
        }

        require(modifiedComponentCount == 1,
                "reduced fixture must have exactly one recomposed local component per level");
    }
    require(recomposedTouchesStableContext,
            "at least one recomposed component must meet preserved or boundary context");

    const auto serialized = noninc::debugString(plan);
    require(serialized.find("LocalPartitionPlan") != std::string::npos,
            "local partition must provide a deterministic debug serialization");
    require(serialized.find("level=80") != std::string::npos,
            "local partition serialization must include affected levels");
    require(serialized.find("touches_modified=true") != std::string::npos,
            "local partition serialization must expose the recomposed component");
}

void testNonOracleLocalPartitionMatchesReducedFixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const auto oracle = buildReducedFixtureOracle(adjacency);
    const auto oraclePlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto nonOraclePlan = noninc::buildLocalPartitionPlanWithoutOracle(input, *adjacency);

    require(partitionPlansMatchIgnoringOracleLabels(nonOraclePlan, oraclePlan),
            "non-oracle local partition must match the reduced fixture topology");

    for (const auto &partition : nonOraclePlan.levels) {
        for (const auto &component : partition.components) {
            require(component.oracleDualOwners.empty(),
                    "non-oracle local partition must not populate oracle dual-owner labels");
            require(component.oracleThresholdOwners.empty(),
                    "non-oracle local partition must not populate oracle threshold-owner labels");
        }
    }

    const auto report = noninc::auditLocalPlanOracleDependencies(&nonOraclePlan);
    require(report.partitionLevels == static_cast<int>(oraclePlan.levels.size()),
            "non-oracle dependency audit must count local partition levels");
    require(report.partitionComponents > 0,
            "non-oracle dependency audit must count partition components");
    require(report.isOracleFree(),
            "non-oracle partition-only plan must be classified as oracle-free");
}

void testNonOracleLocalRecompositionMatchesReducedFixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const auto oracle = buildReducedFixtureOracle(adjacency);
    const auto oraclePartition = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto oracleRecomposition = noninc::buildLocalRecompositionPlan(oracle, oraclePartition);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto strippedBuild = noninc::buildStrippedNonOracleLocalRecomposition(
        input,
        *adjacency);

    require(recompositionPlansMatchIgnoringOracleLabels(strippedBuild.recompositionPlan, oracleRecomposition),
            "non-oracle recomposition must match the reduced fixture topology");
    require(strippedBuild.isOracleFree(),
            "non-oracle recomposition for the reduced fixture must be oracle-free");
}

void testReducedFixtureBuildsEffectiveRecompositionPlan() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const auto oracle = buildReducedFixtureOracle(adjacency);
    const auto partitionPlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto recompositionPlan = noninc::buildLocalRecompositionPlan(oracle, partitionPlan);

    require(partitionPlan.recomposedComponents.size() == 12,
            "fixture must still expose one threshold component for each level in J_C");
    require(recompositionPlan.nodes.size() == 1,
            "fixture threshold components should condense to one effective recomposition node");

    const auto &node = recompositionPlan.nodes.front();
    const PixelId modifiedPixel = oracle.contraction.modifiedProperPart.front();
    require(node.level == oracle.contraction.b,
            "the effective recomposition node must be created at level b");
    require(vectorContainsPixel(node.supportPixels, modifiedPixel),
            "the effective recomposition support must contain the modified pixel");
    require(vectorContainsPixel(node.directProperParts, modifiedPixel),
            "the modified pixel must be a direct proper part of the effective node");
    require(node.matchesOracleNode(),
            "the effective recomposition node must match an oracle local node");
    require(recompositionPlan.allNodesMatchOracle(),
            "every effective recomposition node in the fixture must match the oracle");
    require(node.parentPlanIndex == -1,
            "single-node fixture plan must not invent a local parent");

    const auto serialized = noninc::debugString(recompositionPlan);
    require(serialized.find("LocalRecompositionPlan") != std::string::npos,
            "recomposition plan must provide a deterministic debug serialization");
    require(serialized.find("level=80") != std::string::npos,
            "recomposition plan serialization must include the effective level");
    require(serialized.find("all_nodes_match_oracle=true") != std::string::npos,
            "recomposition plan serialization must expose oracle matching");
}

noninc::LocalMergeOracle makeExternalConnectivityPartitionOracle() {
    auto updated = ImageUInt8::create(3, 5);
    for (int pixelId = 0; pixelId < updated->getSize(); ++pixelId) {
        (*updated)[pixelId] = 30;
    }

    const std::vector<PixelId> level13Pixels = {1, 2, 3, 6, 9, 11, 12, 13};
    for (PixelId pixelId : level13Pixels) {
        (*updated)[pixelId] = 13;
    }
    (*updated)[5] = 11;
    (*updated)[7] = 21;
    (*updated)[8] = 21;

    noninc::LocalMergeOracle oracle;
    oracle.preImage = updated->clone();
    oracle.updatedImage = updated;
    oracle.contraction.a = 20;
    oracle.contraction.b = 13;
    oracle.contraction.dualIsMaxtree = false;
    oracle.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    oracle.contraction.modifiedProperPart = {6};
    oracle.contraction.modifiedComponents = {{6}};
    oracle.preservedPixels = {7, 8};
    oracle.boundaryPixels = {1, 2, 3, 5, 9, 11, 12, 13};
    oracle.localPixels = {1, 2, 3, 5, 6, 7, 8, 9, 11, 12, 13};

    constexpr NodeId oracleLevel11 = 101;
    constexpr NodeId oracleLevel12 = 102;
    constexpr NodeId oracleLevel13 = 103;
    constexpr NodeId oracleLevel14 = 104;
    constexpr NodeId oraclePreserved = 201;

    auto addOwner = [&](PixelId pixelId,
                        noninc::LocalMergePixelRole role,
                        NodeId oracleOwner,
                        int ownerLevel) {
        noninc::LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = role;
        owner.preDualOwner = oracleOwner;
        owner.oracleDualOwner = oracleOwner;
        owner.preDualOwnerLevel = ownerLevel;
        owner.oracleDualOwnerLevel = ownerLevel;
        oracle.pixelOwners.push_back(owner);
    };

    for (PixelId pixelId : level13Pixels) {
        addOwner(pixelId,
                 pixelId == 6 ? noninc::LocalMergePixelRole::Modified
                              : noninc::LocalMergePixelRole::Boundary,
                 oracleLevel13,
                 13);
    }
    addOwner(5, noninc::LocalMergePixelRole::Boundary, oracleLevel11, 11);
    addOwner(7, noninc::LocalMergePixelRole::Preserved, oraclePreserved, 21);
    addOwner(8, noninc::LocalMergePixelRole::Preserved, oraclePreserved, 21);
    std::sort(oracle.pixelOwners.begin(), oracle.pixelOwners.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.pixelId < rhs.pixelId;
    });

    const std::vector<PixelId> level13Support = {1, 2, 3, 5, 6, 9, 11, 12, 13};

    noninc::LocalDualNodeSnapshot node11;
    node11.nodeId = oracleLevel11;
    node11.parentId = oracleLevel12;
    node11.parentInSnapshot = true;
    node11.altitude = 11;
    node11.directLocalProperParts = {5};
    node11.supportLocalPixels = {5};
    node11.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot node12;
    node12.nodeId = oracleLevel12;
    node12.parentId = oracleLevel13;
    node12.parentInSnapshot = true;
    node12.altitude = 12;
    node12.supportLocalPixels = {5};
    node12.childrenInSnapshot = {oracleLevel11};
    node12.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot node13;
    node13.nodeId = oracleLevel13;
    node13.parentId = oracleLevel14;
    node13.parentInSnapshot = true;
    node13.altitude = 13;
    node13.directLocalProperParts = level13Pixels;
    node13.supportLocalPixels = level13Support;
    node13.childrenInSnapshot = {oracleLevel12};
    node13.supportTouchesModified = true;
    node13.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot node14;
    node14.nodeId = oracleLevel14;
    node14.parentId = InvalidNode;
    node14.parentInSnapshot = false;
    node14.altitude = 14;
    node14.supportLocalPixels = level13Support;
    node14.childrenInSnapshot = {oracleLevel13};
    node14.supportTouchesModified = true;
    node14.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preserved;
    preserved.nodeId = oraclePreserved;
    preserved.parentId = InvalidNode;
    preserved.parentInSnapshot = false;
    preserved.altitude = 21;
    preserved.directLocalProperParts = {7, 8};
    preserved.supportLocalPixels = {7, 8};
    preserved.supportTouchesPreserved = true;

    oracle.preDualNodes = {node11, node12, node13, node14, preserved};
    oracle.oracleDualNodes = oracle.preDualNodes;
    oracle.preDualMatchesPrimal = true;
    oracle.oracleDualMatchesUpdatedImage = true;
    return oracle;
}

void testExternalConnectivityMergesLocalThresholdComponentsByOracleOwner() {
    const auto oracle = makeExternalConnectivityPartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto partitionPlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    require(!partitionPlan.levels.empty(),
            "external-connectivity fixture must produce affected levels");
    const auto &level13Partition = partitionPlan.levels.front();
    require(level13Partition.level == 13,
            "external-connectivity fixture must start at level 13");
    require(level13Partition.components.size() == 1,
            "oracle threshold owner must merge local pieces connected outside the local window");

    const auto &component = level13Partition.components.front();
    require((component.pixels == std::vector<PixelId>{1, 2, 3, 5, 6, 9, 11, 12, 13}),
            "merged local component must match the oracle local support intersection");
    require(component.touchesModified,
            "merged external-connectivity component must remain recomposed");
    require(vectorContainsIndex(component.oracleThresholdOwners, 103),
            "merged component must record the shared oracle threshold owner");

    const auto recompositionPlan = noninc::buildLocalRecompositionPlan(oracle, partitionPlan);
    require(recompositionPlan.nodes.size() == 2,
            "external-connectivity fixture must create the level-13 node and its level-14 parent");
    const auto level13It = std::find_if(
        recompositionPlan.nodes.begin(),
        recompositionPlan.nodes.end(),
        [](const auto &node) { return node.level == 13; });
    const auto level14It = std::find_if(
        recompositionPlan.nodes.begin(),
        recompositionPlan.nodes.end(),
        [](const auto &node) { return node.level == 14; });
    require(level13It != recompositionPlan.nodes.end(),
            "external-connectivity recomposition plan must include level 13");
    require(level14It != recompositionPlan.nodes.end(),
            "external-connectivity recomposition plan must include level 14");
    require(level13It->oracleNodeId == 103 && level14It->oracleNodeId == 104,
            "external-connectivity plan nodes must match the oracle hierarchy");
    require(vectorContainsPixel(level13It->directProperParts, 9),
            "externally connected boundary pixel must become a direct part of the matched oracle node");
    require(level13It->parentPlanIndex ==
                static_cast<int>(std::distance(recompositionPlan.nodes.begin(), level14It)),
            "level-13 external-connectivity node must attach under the local level-14 oracle parent");
}

void testNonOracleLocalPartitionMatchesExternalConnectivityFixture() {
    const auto oracle = makeExternalConnectivityPartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto oraclePlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto nonOraclePlan = noninc::buildLocalPartitionPlanWithoutOracle(input, *adjacency);

    require(partitionPlansMatchIgnoringOracleLabels(nonOraclePlan, oraclePlan),
            "non-oracle partition must match external-connectivity topology via pre-snapshot support");
    require(nonOraclePlan.levels.front().components.size() == 1,
            "non-oracle external-connectivity fixture must merge the level-13 pieces");

    const auto report = noninc::auditLocalPlanOracleDependencies(&nonOraclePlan);
    require(report.isOracleFree(),
            "non-oracle external-connectivity partition must not use oracle labels");
}

void testNonOracleLocalRecompositionMatchesExternalConnectivityFixture() {
    const auto oracle = makeExternalConnectivityPartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto oraclePartition = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto oracleRecomposition = noninc::buildLocalRecompositionPlan(oracle, oraclePartition);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto strippedBuild = noninc::buildStrippedNonOracleLocalRecomposition(
        input,
        *adjacency);

    require(recompositionPlansMatchIgnoringOracleLabels(strippedBuild.recompositionPlan, oracleRecomposition),
            "non-oracle recomposition must match external-connectivity topology");
    require(strippedBuild.recompositionPlan.nodes.size() == 2,
            "non-oracle external-connectivity recomposition must create level-13 and level-14 nodes");
    require(strippedBuild.recompositionPlan.nodes.front().parentPlanIndex == 1,
            "non-oracle external-connectivity child must attach under the same-support level-14 parent");

    require(strippedBuild.isOracleFree(),
            "non-oracle external-connectivity recomposition must not use oracle labels");
}

noninc::LocalMergeOracle makeAncestorSplitPreservationPartitionOracle() {
    auto updated = ImageUInt8::create(1, 5);
    (*updated)[0] = 80;
    (*updated)[1] = 80;
    (*updated)[2] = 120;
    (*updated)[3] = 80;
    (*updated)[4] = 80;

    noninc::LocalMergeOracle oracle;
    oracle.preImage = updated->clone();
    oracle.updatedImage = updated;
    oracle.contraction.a = 92;
    oracle.contraction.b = 80;
    oracle.contraction.dualIsMaxtree = false;
    oracle.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    oracle.contraction.modifiedProperPart = {3};
    oracle.contraction.modifiedComponents = {{3}};
    oracle.boundaryPixels = {0, 1, 4};
    oracle.localPixels = {0, 1, 3, 4};

    constexpr NodeId preLeftChild = 701;
    constexpr NodeId preRightChild = 702;
    constexpr NodeId preAncestor = 703;
    constexpr NodeId oracleLeftChild = 801;
    constexpr NodeId oracleRightChild = 802;
    constexpr NodeId oracleAncestor = 803;

    auto addOwner = [&](PixelId pixelId,
                        noninc::LocalMergePixelRole role,
                        NodeId preOwner,
                        NodeId oracleOwner,
                        int preOwnerLevel,
                        int oracleOwnerLevel) {
        noninc::LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = role;
        owner.preDualOwner = preOwner;
        owner.oracleDualOwner = oracleOwner;
        owner.preDualOwnerLevel = preOwnerLevel;
        owner.oracleDualOwnerLevel = oracleOwnerLevel;
        oracle.pixelOwners.push_back(owner);
    };

    addOwner(0, noninc::LocalMergePixelRole::Boundary, preLeftChild, oracleLeftChild, 80, 80);
    addOwner(1, noninc::LocalMergePixelRole::Boundary, preLeftChild, oracleLeftChild, 80, 80);
    addOwner(3, noninc::LocalMergePixelRole::Modified, preAncestor, oracleRightChild, 92, 80);
    addOwner(4, noninc::LocalMergePixelRole::Boundary, preRightChild, oracleRightChild, 80, 80);

    noninc::LocalDualNodeSnapshot preLeft;
    preLeft.nodeId = preLeftChild;
    preLeft.parentId = preAncestor;
    preLeft.parentInSnapshot = true;
    preLeft.altitude = 80;
    preLeft.directLocalProperParts = {0, 1};
    preLeft.supportLocalPixels = {0, 1};
    preLeft.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preRight;
    preRight.nodeId = preRightChild;
    preRight.parentId = preAncestor;
    preRight.parentInSnapshot = true;
    preRight.altitude = 80;
    preRight.directLocalProperParts = {4};
    preRight.supportLocalPixels = {4};
    preRight.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preParent;
    preParent.nodeId = preAncestor;
    preParent.parentId = InvalidNode;
    preParent.parentInSnapshot = false;
    preParent.altitude = 92;
    preParent.directLocalProperParts = {3};
    preParent.supportLocalPixels = {0, 1, 3, 4};
    preParent.childrenInSnapshot = {preLeftChild, preRightChild};
    preParent.supportTouchesModified = true;
    preParent.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot oracleLeft = preLeft;
    oracleLeft.nodeId = oracleLeftChild;
    oracleLeft.parentId = oracleAncestor;

    noninc::LocalDualNodeSnapshot oracleRight;
    oracleRight.nodeId = oracleRightChild;
    oracleRight.parentId = oracleAncestor;
    oracleRight.parentInSnapshot = true;
    oracleRight.altitude = 80;
    oracleRight.directLocalProperParts = {3, 4};
    oracleRight.supportLocalPixels = {3, 4};
    oracleRight.supportTouchesModified = true;
    oracleRight.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot oracleParent;
    oracleParent.nodeId = oracleAncestor;
    oracleParent.parentId = InvalidNode;
    oracleParent.parentInSnapshot = false;
    oracleParent.altitude = 92;
    oracleParent.supportLocalPixels = {0, 1, 3, 4};
    oracleParent.childrenInSnapshot = {oracleLeftChild, oracleRightChild};
    oracleParent.supportTouchesModified = true;
    oracleParent.supportTouchesBoundary = true;

    oracle.preDualNodes = {preLeft, preRight, preParent};
    oracle.oracleDualNodes = {oracleLeft, oracleRight, oracleParent};
    oracle.preDualMatchesPrimal = true;
    oracle.oracleDualMatchesUpdatedImage = true;
    return oracle;
}

void testNonOracleLocalPartitionPreservesProperThresholdChildrenUnderAncestor() {
    const auto oracle = makeAncestorSplitPreservationPartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto oraclePartition = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto nonOraclePartition =
        noninc::buildLocalPartitionPlanWithoutOracle(input, *adjacency);

    require(oraclePartition.levels.size() == 12 &&
            nonOraclePartition.levels.size() == 12,
            "ancestor split fixture must span the expected lower interval");
    require(oraclePartition.levels.front().components.size() == 2,
            "oracle fixture must keep the two level-b siblings separate");
    require(nonOraclePartition.levels.front().components.size() == 2,
            "non-oracle partition must not merge siblings through their ancestor");
    require(partitionPlansMatchIgnoringOracleLabels(nonOraclePartition, oraclePartition),
            "non-oracle partition must preserve proper threshold children under an ancestor");

    const auto oracleRecomposition =
        noninc::buildLocalRecompositionPlan(oracle, oraclePartition);
    const auto nonOracleRecomposition =
        noninc::buildLocalRecompositionPlanWithoutOracle(input, nonOraclePartition);

    require(recompositionPlansMatchIgnoringOracleLabels(nonOracleRecomposition, oracleRecomposition),
            "non-oracle recomposition must match ancestor split fixture topology");
    require(nonOracleRecomposition.nodes.size() == 1,
            "ancestor split fixture must materialize only the modified level-b sibling");
    require(nonOracleRecomposition.nodes.front().supportPixels == std::vector<PixelId>({3, 4}),
            "ancestor split fixture recomposition node must exclude the stable sibling");
}

noninc::LocalMergeOracle makeOracleSnapshotCoveragePartitionOracle() {
    auto updated = ImageUInt8::create(3, 5);
    for (int pixelId = 0; pixelId < updated->getSize(); ++pixelId) {
        (*updated)[pixelId] = 30;
    }

    const std::vector<PixelId> separatedLevel13Pixels = {1, 3, 11};
    for (PixelId pixelId : separatedLevel13Pixels) {
        (*updated)[pixelId] = 13;
    }

    noninc::LocalMergeOracle oracle;
    oracle.preImage = updated->clone();
    oracle.updatedImage = updated;
    oracle.contraction.a = 14;
    oracle.contraction.b = 13;
    oracle.contraction.dualIsMaxtree = false;
    oracle.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    oracle.contraction.modifiedProperPart = {1};
    oracle.contraction.modifiedComponents = {{1}};
    oracle.boundaryPixels = {3, 11};
    oracle.localPixels = separatedLevel13Pixels;

    constexpr NodeId oracleTarget = 500;
    constexpr NodeId preIslandA = 501;
    constexpr NodeId preIslandB = 502;
    constexpr NodeId preIslandC = 503;
    constexpr NodeId oracleParent = 504;

    auto addOwner = [&](PixelId pixelId,
                        noninc::LocalMergePixelRole role,
                        NodeId ownerId) {
        noninc::LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = role;
        owner.preDualOwner = ownerId;
        owner.oracleDualOwner = InvalidNode;
        owner.preDualOwnerLevel = 13;
        owner.oracleDualOwnerLevel = -1;
        oracle.pixelOwners.push_back(owner);
    };

    addOwner(1, noninc::LocalMergePixelRole::Modified, preIslandA);
    addOwner(3, noninc::LocalMergePixelRole::Boundary, preIslandB);
    addOwner(11, noninc::LocalMergePixelRole::Boundary, preIslandC);

    auto makeIslandSnapshot = [](NodeId nodeId, PixelId pixelId) {
        noninc::LocalDualNodeSnapshot snapshot;
        snapshot.nodeId = nodeId;
        snapshot.parentId = InvalidNode;
        snapshot.parentInSnapshot = false;
        snapshot.altitude = 13;
        snapshot.directLocalProperParts = {pixelId};
        snapshot.supportLocalPixels = {pixelId};
        return snapshot;
    };

    noninc::LocalDualNodeSnapshot islandA = makeIslandSnapshot(preIslandA, 1);
    islandA.supportTouchesModified = true;
    noninc::LocalDualNodeSnapshot islandB = makeIslandSnapshot(preIslandB, 3);
    islandB.supportTouchesBoundary = true;
    noninc::LocalDualNodeSnapshot islandC = makeIslandSnapshot(preIslandC, 11);
    islandC.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot target;
    target.nodeId = oracleTarget;
    target.parentId = oracleParent;
    target.parentInSnapshot = true;
    target.altitude = 13;
    target.directLocalProperParts = separatedLevel13Pixels;
    target.supportLocalPixels = separatedLevel13Pixels;
    target.supportTouchesModified = true;
    target.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot parent;
    parent.nodeId = oracleParent;
    parent.parentId = InvalidNode;
    parent.parentInSnapshot = false;
    parent.altitude = 14;
    parent.supportLocalPixels = separatedLevel13Pixels;
    parent.childrenInSnapshot = {oracleTarget};
    parent.supportTouchesModified = true;
    parent.supportTouchesBoundary = true;

    oracle.preDualNodes = {islandA, islandB, islandC, parent};
    oracle.oracleDualNodes = {target, parent};
    oracle.preDualMatchesPrimal = true;
    oracle.oracleDualMatchesUpdatedImage = true;
    return oracle;
}

void testOracleSnapshotCoverageMergesDisjointLocalThresholdComponents() {
    const auto oracle = makeOracleSnapshotCoveragePartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto partitionPlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    require(partitionPlan.levels.size() == 1,
            "snapshot-coverage fixture must have exactly one affected level");
    const auto &level13Partition = partitionPlan.levels.front();
    require(level13Partition.level == 13,
            "snapshot-coverage fixture must affect level 13");
    require(level13Partition.components.size() == 1,
            "oracle snapshot coverage must merge disjoint local pieces");

    const auto &component = level13Partition.components.front();
    require((component.pixels == std::vector<PixelId>{1, 3, 11}),
            "snapshot-coverage merged component must match the oracle support");
    require(component.touchesModified && component.touchesBoundary,
            "snapshot-coverage merged component must retain local roles");
    require(component.oracleThresholdOwners.empty(),
            "snapshot-coverage merge must not depend on a shared threshold owner");

    const auto recompositionPlan = noninc::buildLocalRecompositionPlan(oracle, partitionPlan);
    require(recompositionPlan.nodes.size() == 1,
            "snapshot-coverage fixture must create only the affected target node");
    const auto &targetNode = recompositionPlan.nodes.front();
    require(targetNode.oracleNodeId == 500 && targetNode.oracleParentId == 504,
            "snapshot-coverage recomposition plan must preserve the oracle target and external parent");
    require(targetNode.parentPlanIndex == -1,
            "snapshot-coverage target parent must remain external to the affected levels");
}

void testNonOracleLocalPartitionMergesSnapshotCoverageByPreSnapshotSupport() {
    const auto oracle = makeOracleSnapshotCoveragePartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto oraclePlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto nonOraclePlan = noninc::buildLocalPartitionPlanWithoutOracle(input, *adjacency);

    require(oraclePlan.levels.size() == 1 && nonOraclePlan.levels.size() == 1,
            "snapshot-coverage fixture must keep one affected level");
    require(oraclePlan.levels.front().components.size() == 1,
            "oracle-guided snapshot-coverage partition must merge the disjoint pieces");
    require(nonOraclePlan.levels.front().components.size() == 1,
            "non-oracle partition must merge snapshot coverage through pre-update support");
    require(nonOraclePlan.recomposedComponents.size() == 1,
            "non-oracle builder must identify the merged modified local component");
    require(partitionPlansMatchIgnoringOracleLabels(nonOraclePlan, oraclePlan),
            "non-oracle snapshot-coverage partition must match the oracle-guided topology");

    const auto nonOracleReport = noninc::auditLocalPlanOracleDependencies(&nonOraclePlan);
    require(nonOracleReport.isOracleFree(),
            "snapshot-coverage non-oracle partition must not use hidden oracle labels");
    require(nonOracleReport.partitionComponents == 1,
            "snapshot-coverage report must count the merged local component");
}

void testNonOracleLocalRecompositionMatchesSnapshotCoverageFixture() {
    const auto oracle = makeOracleSnapshotCoveragePartitionOracle();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        oracle.updatedImage->getNumRows(),
        oracle.updatedImage->getNumCols(),
        1.0);

    const auto oraclePartition = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto oracleRecomposition = noninc::buildLocalRecompositionPlan(oracle, oraclePartition);
    const auto input = noninc::makeLocalPartitionInputFromOraclePreState(oracle);
    const auto strippedBuild = noninc::buildStrippedNonOracleLocalRecomposition(
        input,
        *adjacency);

    require(recompositionPlansMatchIgnoringOracleLabels(strippedBuild.recompositionPlan, oracleRecomposition),
            "non-oracle recomposition must match snapshot-coverage topology");
    require(strippedBuild.recompositionPlan.nodes.size() == 1,
            "non-oracle snapshot-coverage recomposition must create one target node");
    require(strippedBuild.recompositionPlan.nodes.front().parentPlanIndex == -1,
            "non-oracle snapshot-coverage target parent must remain external");

    require(strippedBuild.isOracleFree(),
            "non-oracle snapshot-coverage recomposition must not use oracle labels");
}

void testNonOracleLocalRecompositionRecordsMovedStableChildren() {
    auto image = make_split_child_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId firstStableChild = minTree.getSmallestComponent(0);
    const NodeId oldSourceOwner = minTree.getSmallestComponent(2);
    const NodeId secondStableChild = minTree.createDetachedNode(minTree.getAltitude(firstStableChild));
    require(firstStableChild != InvalidNode &&
            oldSourceOwner != InvalidNode &&
            secondStableChild != InvalidNode,
            "non-oracle stable-child setup must create valid nodes");
    minTree.moveProperPart(secondStableChild, firstStableChild, 1);
    minTree.attachNode(oldSourceOwner, secondStableChild);

    noninc::LocalDualNodeSnapshot firstChildSnapshot;
    firstChildSnapshot.nodeId = firstStableChild;
    firstChildSnapshot.parentId = oldSourceOwner;
    firstChildSnapshot.parentInSnapshot = true;
    firstChildSnapshot.altitude = 10;
    firstChildSnapshot.directLocalProperParts = {0};
    firstChildSnapshot.supportLocalPixels = {0};
    firstChildSnapshot.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot secondChildSnapshot = firstChildSnapshot;
    secondChildSnapshot.nodeId = secondStableChild;
    secondChildSnapshot.directLocalProperParts = {1};
    secondChildSnapshot.supportLocalPixels = {1};

    noninc::LocalDualNodeSnapshot sourceSnapshot;
    sourceSnapshot.nodeId = oldSourceOwner;
    sourceSnapshot.parentId = InvalidNode;
    sourceSnapshot.parentInSnapshot = false;
    sourceSnapshot.altitude = 30;
    sourceSnapshot.directLocalProperParts = {2};
    sourceSnapshot.supportLocalPixels = {0, 1, 2};
    sourceSnapshot.childrenInSnapshot = {firstStableChild, secondStableChild};
    sourceSnapshot.supportTouchesModified = true;
    sourceSnapshot.supportTouchesPreserved = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 10, 20, 40});
    input.contraction.dualIsMaxtree = false;
    input.contraction.a = 30;
    input.contraction.b = 20;
    input.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    input.localPixels = {0, 1, 2};
    input.preDualNodes = {firstChildSnapshot, secondChildSnapshot, sourceSnapshot};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Preserved, firstStableChild},
        {1, noninc::LocalMergePixelRole::Preserved, secondStableChild},
        {2, noninc::LocalMergePixelRole::Modified, oldSourceOwner},
    };

    const auto strippedBuild = noninc::buildStrippedNonOracleLocalRecomposition(
        input,
        *adjacency);
    const auto &recompositionPlan = strippedBuild.recompositionPlan;

    require(recompositionPlan.nodes.size() == 1,
            "non-oracle stable-child fixture must create one recomposition node");
    const auto &node = recompositionPlan.nodes.front();
    require(node.level == 20 && node.directProperParts == std::vector<PixelId>{2},
            "non-oracle stable-child recomposition node must own the modified pixel at level 20");
    require((node.stableChildPreNodes == std::vector<NodeId>{firstStableChild, secondStableChild}),
            "non-oracle recomposition must record both stable children to move");
    require(strippedBuild.isOracleFree(),
            "non-oracle stable-child recomposition helper must keep the plan oracle-free");
}

void testStrippedNonOraclePipelineMaterializesStableChildrenEndToEnd() {
    auto image = make_split_child_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId firstStableChild = minTree.getSmallestComponent(0);
    const NodeId oldSourceOwner = minTree.getSmallestComponent(2);
    const NodeId secondStableChild = minTree.createDetachedNode(minTree.getAltitude(firstStableChild));
    require(firstStableChild != InvalidNode &&
            oldSourceOwner != InvalidNode &&
            secondStableChild != InvalidNode,
            "stripped non-oracle setup must create valid source and stable-child nodes");
    minTree.moveProperPart(secondStableChild, firstStableChild, 1);
    minTree.attachNode(oldSourceOwner, secondStableChild);
    requireTreeConsistency(minTree);

    noninc::LocalDualNodeSnapshot firstChildSnapshot;
    firstChildSnapshot.nodeId = firstStableChild;
    firstChildSnapshot.parentId = oldSourceOwner;
    firstChildSnapshot.parentInSnapshot = true;
    firstChildSnapshot.altitude = 10;
    firstChildSnapshot.directLocalProperParts = {0};
    firstChildSnapshot.supportLocalPixels = {0};
    firstChildSnapshot.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot secondChildSnapshot = firstChildSnapshot;
    secondChildSnapshot.nodeId = secondStableChild;
    secondChildSnapshot.directLocalProperParts = {1};
    secondChildSnapshot.supportLocalPixels = {1};

    noninc::LocalDualNodeSnapshot sourceSnapshot;
    sourceSnapshot.nodeId = oldSourceOwner;
    sourceSnapshot.parentId = InvalidNode;
    sourceSnapshot.parentInSnapshot = false;
    sourceSnapshot.altitude = 30;
    sourceSnapshot.directLocalProperParts = {2};
    sourceSnapshot.supportLocalPixels = {0, 1, 2};
    sourceSnapshot.childrenInSnapshot = {firstStableChild, secondStableChild};
    sourceSnapshot.supportTouchesModified = true;
    sourceSnapshot.supportTouchesPreserved = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 10, 20, 40});
    input.contraction.dualIsMaxtree = false;
    input.contraction.a = 30;
    input.contraction.b = 20;
    input.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    input.localPixels = {0, 1, 2};
    input.preDualNodes = {firstChildSnapshot, secondChildSnapshot, sourceSnapshot};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Preserved, firstStableChild},
        {1, noninc::LocalMergePixelRole::Preserved, secondStableChild},
        {2, noninc::LocalMergePixelRole::Modified, oldSourceOwner},
    };

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecomposition(
        minTree,
        input,
        *adjacency);
    const auto &recompositionPlan = strippedResult.build.recompositionPlan;

    require(recompositionPlan.nodes.size() == 1,
            "stripped non-oracle pipeline must create one recomposition node");
    const auto &plannedNode = recompositionPlan.nodes.front();
    require(plannedNode.level == 20,
            "stripped non-oracle recomposition node must use the lowered level");
    require(plannedNode.directProperParts == std::vector<PixelId>{2},
            "stripped non-oracle recomposition node must own the modified pixel");
    require((plannedNode.stableChildPreNodes == std::vector<NodeId>{firstStableChild, secondStableChild}),
            "stripped non-oracle plan must carry stable pre-child nodes into materialization");

    require(strippedResult.build.isOracleFree(),
            "stripped non-oracle partition/recomposition plan must not use oracle labels");

    require(strippedResult.materializationContext.oracleDualNodes.empty() &&
            strippedResult.materializationContext.pixelOwners.empty(),
            "stripped materialization helper must keep oracle snapshots and owners out of the context");
    require(strippedResult.materializationDependencyReport.isOracleFreeMaterializationInput(),
            "stripped materialization context must not contain oracle materialization inputs");
    require(strippedResult.materializationDependencyReport.usesLocalPreStateInputs(),
            "stripped materialization context must rely on local pre-state snapshots");
    require(strippedResult.materializationDependencyReport.oracleDualSnapshots == 0 &&
            strippedResult.materializationDependencyReport.pixelOwnersWithOracleDualOwner == 0,
            "stripped materialization context must keep oracle snapshots and owners empty");
    require(strippedResult.isOracleFree(),
            "stripped non-oracle helper result must be fully oracle-free");

    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 1,
            "stripped non-oracle materialization must create one visible recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "stripped non-oracle materialization must reconstruct the target image");

    const auto &created = materialization.createdNodes.front();
    require(created.level == 20,
            "stripped non-oracle created node must keep the planned level");
    require(created.movedDirectProperParts == std::vector<PixelId>{2},
            "stripped non-oracle materialization must report the moved modified pixel");
    require(created.movedChildNodes.size() == 2,
            "stripped non-oracle materialization must report both stable children");
    require(vectorContainsIndex(created.movedChildNodes, firstStableChild),
            "stripped non-oracle materialization must move the first stable child");
    require(vectorContainsIndex(created.movedChildNodes, secondStableChild),
            "stripped non-oracle materialization must move the second stable child");
    require(minTree.getNodeParent(firstStableChild) == created.nodeId,
            "first stable child must attach below the stripped recomposition node");
    require(minTree.getNodeParent(secondStableChild) == created.nodeId,
            "second stable child must attach below the stripped recomposition node");
    require(minTree.getSmallestComponent(0) == firstStableChild,
            "pixel 0 must remain owned by the first stable child");
    require(minTree.getSmallestComponent(1) == secondStableChild,
            "pixel 1 must remain owned by the second stable child");
    require(minTree.getSmallestComponent(2) == created.nodeId,
            "pixel 2 must move to the stripped recomposition node");
    require(!minTree.isAlive(oldSourceOwner),
            "obsolete source owner must be removed by stripped local materialization/canonicalization");
    require(minTree.reconstructionImage()->isEqual(input.updatedImage),
            "stripped non-oracle materialization must match the target image exactly");

    DynamicComponentTree expectedTree(input.updatedImage->clone(), false, adjacency);
    require(minTree.reconstructionImage()->isEqual(expectedTree.reconstructionImage()),
            "stripped non-oracle materialization must match a full rebuild reconstruction");
    requireTreeConsistency(minTree);
}

void testNonOracleLocalRecompositionRecordsMaximalStableChild() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId stableChild = minTree.getSmallestComponent(1);
    const NodeId sourceOwner = minTree.getSmallestComponent(2);
    require(lowNode != InvalidNode &&
            stableChild != InvalidNode &&
            sourceOwner != InvalidNode,
            "non-oracle maximal-stable setup must have low, stable, and source nodes");

    noninc::LocalDualNodeSnapshot lowSnapshot;
    lowSnapshot.nodeId = lowNode;
    lowSnapshot.parentId = stableChild;
    lowSnapshot.parentInSnapshot = true;
    lowSnapshot.altitude = 10;
    lowSnapshot.directLocalProperParts = {0};
    lowSnapshot.supportLocalPixels = {0};
    lowSnapshot.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot stableSnapshot;
    stableSnapshot.nodeId = stableChild;
    stableSnapshot.parentId = sourceOwner;
    stableSnapshot.parentInSnapshot = true;
    stableSnapshot.altitude = 20;
    stableSnapshot.directLocalProperParts = {1};
    stableSnapshot.supportLocalPixels = {0, 1};
    stableSnapshot.childrenInSnapshot = {lowNode};
    stableSnapshot.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot sourceSnapshot;
    sourceSnapshot.nodeId = sourceOwner;
    sourceSnapshot.parentId = InvalidNode;
    sourceSnapshot.parentInSnapshot = false;
    sourceSnapshot.altitude = 30;
    sourceSnapshot.directLocalProperParts = {2};
    sourceSnapshot.supportLocalPixels = {0, 1, 2};
    sourceSnapshot.childrenInSnapshot = {stableChild};
    sourceSnapshot.supportTouchesModified = true;
    sourceSnapshot.supportTouchesPreserved = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 20, 25, 40});
    input.contraction.dualIsMaxtree = false;
    input.contraction.a = 30;
    input.contraction.b = 25;
    input.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    input.localPixels = {0, 1, 2};
    input.preDualNodes = {lowSnapshot, stableSnapshot, sourceSnapshot};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Preserved, lowNode},
        {1, noninc::LocalMergePixelRole::Preserved, stableChild},
        {2, noninc::LocalMergePixelRole::Modified, sourceOwner},
    };

    const auto strippedBuild = noninc::buildStrippedNonOracleLocalRecomposition(
        input,
        *adjacency);
    const auto &recompositionPlan = strippedBuild.recompositionPlan;

    require(recompositionPlan.nodes.size() == 1,
            "non-oracle maximal-stable fixture must create one recomposition node");
    const auto &node = recompositionPlan.nodes.front();
    require(node.level == 25 && node.directProperParts == std::vector<PixelId>{2},
            "non-oracle maximal-stable recomposition node must own the modified pixel at level 25");
    require((node.stableChildPreNodes == std::vector<NodeId>{stableChild}),
            "non-oracle recomposition must record only the maximal stable child");
    require(strippedBuild.isOracleFree(),
            "non-oracle maximal-stable recomposition helper must keep the plan oracle-free");
}

void testLocalOracleDependencyAuditSeparatesCurrentOracleUses() {
    auto snapshotOracle = makeOracleSnapshotCoveragePartitionOracle();
    auto snapshotAdjacency = std::make_shared<AdjacencyRelation>(
        snapshotOracle.updatedImage->getNumRows(),
        snapshotOracle.updatedImage->getNumCols(),
        1.0);
    const auto snapshotPartition = noninc::buildLocalPartitionPlan(
        snapshotOracle,
        *snapshotAdjacency);
    const auto snapshotRecomposition = noninc::buildLocalRecompositionPlan(
        snapshotOracle,
        snapshotPartition);
    const auto snapshotReport = noninc::auditLocalOracleDependencies(
        snapshotOracle,
        &snapshotPartition,
        &snapshotRecomposition);

    require(snapshotReport.pixelOwners == 3,
            "snapshot dependency audit must count pixel-owner metadata");
    require(snapshotReport.pixelOwnersWithOracleDualOwner == 0,
            "snapshot-coverage fixture must not depend on oracle owner labels");
    require(snapshotReport.partitionComponentsWithOracleThresholdOwners == 0,
            "snapshot-coverage fixture must not use shared threshold owners");
    require(snapshotReport.partitionComponentsMatchingOracleSnapshot == 1,
            "snapshot-coverage fixture must record the oracle support match");
    require(snapshotReport.usesOracleSnapshotSupport(),
            "snapshot-coverage fixture must be classified as snapshot-support guided");
    require(snapshotReport.usesOracleRecompositionGuidance(),
            "snapshot-coverage fixture still uses oracle recomposition ids");
    require(!snapshotReport.isOracleFree(),
            "current snapshot-coverage path must not be classified as oracle-free");

    auto externalOracle = makeExternalConnectivityPartitionOracle();
    auto externalAdjacency = std::make_shared<AdjacencyRelation>(
        externalOracle.updatedImage->getNumRows(),
        externalOracle.updatedImage->getNumCols(),
        1.0);
    const auto externalPartition = noninc::buildLocalPartitionPlan(
        externalOracle,
        *externalAdjacency);
    const auto externalRecomposition = noninc::buildLocalRecompositionPlan(
        externalOracle,
        externalPartition);
    const auto externalReport = noninc::auditLocalOracleDependencies(
        externalOracle,
        &externalPartition,
        &externalRecomposition);

    require(externalReport.pixelOwnersWithOracleDualOwner > 0,
            "external-connectivity fixture must expose oracle owner labels");
    require(externalReport.partitionComponentsWithOracleThresholdOwners > 0,
            "external-connectivity fixture must expose shared threshold-owner guidance");
    require(externalReport.usesOracleDualOwnerLabels(),
            "external-connectivity fixture must be classified as owner-label guided");
    require(externalReport.usesOraclePartitionGuidance(),
            "external-connectivity fixture must require oracle partition guidance today");

    const auto serialized = noninc::debugString(snapshotReport);
    require(serialized.find("LocalOracleDependencyReport") != std::string::npos,
            "dependency audit must provide deterministic debug serialization");
    require(serialized.find("uses_oracle_snapshot_support=true") != std::string::npos,
            "dependency audit serialization must expose snapshot-support usage");
    require(serialized.find("is_oracle_free=false") != std::string::npos,
            "dependency audit serialization must expose oracle-free status");
}

void testLocalMaterializationDependencyAuditSeparatesNonOracleInputs() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId parentNode = minTree.getSmallestComponent(1);
    require(lowNode != InvalidNode && parentNode != InvalidNode,
            "materialization audit setup must have low and parent nodes");

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = parentNode;
    preLow.parentInSnapshot = true;
    preLow.altitude = minTree.getAltitude(lowNode);
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot preParent;
    preParent.nodeId = parentNode;
    preParent.parentId = InvalidNode;
    preParent.parentInSnapshot = false;
    preParent.altitude = minTree.getAltitude(parentNode);
    preParent.directLocalProperParts = {1};
    preParent.supportLocalPixels = {0, 1};
    preParent.childrenInSnapshot = {lowNode};
    preParent.supportTouchesModified = true;
    preParent.supportTouchesPreserved = true;

    noninc::LocalMergeOracle nonOracleContext;
    nonOracleContext.preDualNodes = {preLow, preParent};

    noninc::LocalRecompositionPlan nonOraclePlan;
    nonOraclePlan.nodes.resize(1);
    nonOraclePlan.nodes[0].level = 20;
    nonOraclePlan.nodes[0].supportPixels = {0, 1};
    nonOraclePlan.nodes[0].directProperParts = {1};
    nonOraclePlan.nodes[0].stableChildPreNodes = {lowNode};

    const auto nonOracleReport = noninc::auditLocalMaterializationOracleDependencies(
        nonOracleContext,
        nonOraclePlan);
    require(nonOracleReport.preDualSnapshots == 2,
            "materialization audit must count local pre-state snapshots");
    require(nonOracleReport.oracleDualSnapshots == 0,
            "stripped materialization context must not expose oracle snapshots");
    require(nonOracleReport.planNodesWithOracleNodeId == 0 &&
            nonOracleReport.planNodesWithOracleParent == 0 &&
            nonOracleReport.planNodesWithOracleChildren == 0,
            "non-oracle materialization plan must not expose oracle ids");
    require(nonOracleReport.planNodesWithStablePreChildren == 1 &&
            nonOracleReport.stablePreChildren == 1,
            "materialization audit must count non-oracular stable-child inputs");
    require(nonOracleReport.usesLocalPreStateInputs(),
            "non-oracle materialization path must be classified as local-pre-state driven");
    require(!nonOracleReport.mayUseOracleGuidedMaterialization(),
            "stripped non-oracle materialization must not be classified as oracle-guided");
    require(nonOracleReport.isOracleFreeMaterializationInput(),
            "stripped non-oracle materialization input must be oracle-free");

    constexpr NodeId oracleNodeId = 7001;
    constexpr NodeId oracleParentId = 7002;
    constexpr NodeId oracleChildId = 7003;
    noninc::LocalMergeOracle oracleContext;
    oracleContext.preDualNodes = {preLow, preParent};
    noninc::LocalMergePixelOwner owner;
    owner.pixelId = 0;
    owner.role = noninc::LocalMergePixelRole::Preserved;
    owner.preDualOwner = lowNode;
    owner.oracleDualOwner = oracleChildId;
    oracleContext.pixelOwners.push_back(owner);

    noninc::LocalDualNodeSnapshot oracleChild;
    oracleChild.nodeId = oracleChildId;
    oracleChild.parentId = oracleNodeId;
    oracleChild.parentInSnapshot = true;
    oracleChild.altitude = 10;
    oracleChild.directLocalProperParts = {0};
    oracleChild.supportLocalPixels = {0};
    oracleChild.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot oracleNode;
    oracleNode.nodeId = oracleNodeId;
    oracleNode.parentId = oracleParentId;
    oracleNode.parentInSnapshot = true;
    oracleNode.altitude = 20;
    oracleNode.directLocalProperParts = {1};
    oracleNode.supportLocalPixels = {0, 1};
    oracleNode.childrenInSnapshot = {oracleChildId};
    oracleNode.supportTouchesModified = true;
    oracleNode.supportTouchesPreserved = true;

    oracleContext.oracleDualNodes = {oracleChild, oracleNode};

    noninc::LocalRecompositionPlan oraclePlan;
    oraclePlan.nodes.resize(1);
    oraclePlan.nodes[0].level = 20;
    oraclePlan.nodes[0].supportPixels = {0, 1};
    oraclePlan.nodes[0].directProperParts = {1};
    oraclePlan.nodes[0].oracleNodeId = oracleNodeId;
    oraclePlan.nodes[0].oracleParentId = oracleParentId;
    oraclePlan.nodes[0].oracleChildrenInSnapshot = {oracleChildId};

    const auto oracleReport = noninc::auditLocalMaterializationOracleDependencies(
        oracleContext,
        oraclePlan);
    require(oracleReport.oracleDualSnapshots == 2,
            "oracle-guided materialization audit must count oracle snapshots");
    require(oracleReport.pixelOwnersWithOracleDualOwner == 1,
            "oracle-guided materialization audit must count oracle owner labels");
    require(oracleReport.planNodesWithOracleNodeId == 1 &&
            oracleReport.planNodesWithOracleParent == 1 &&
            oracleReport.planNodesWithOracleChildren == 1,
            "oracle-guided materialization audit must count oracle plan fields");
    require(oracleReport.mayUseOracleGuidedMaterialization(),
            "oracle-guided materialization input must be classified as oracle-dependent");
    require(!oracleReport.isOracleFreeMaterializationInput(),
            "oracle-guided materialization input must not be oracle-free");

    const auto serialized = noninc::debugString(nonOracleReport);
    require(serialized.find("LocalMaterializationOracleDependencyReport") != std::string::npos,
            "materialization dependency audit must provide debug serialization");
    require(serialized.find("is_oracle_free_materialization_input=true") != std::string::npos,
            "materialization dependency audit serialization must expose oracle-free status");
}

void testSyntheticRecompositionParentIndicesRespectPolarity() {
    noninc::LocalRecompositionPlan minPlan;
    minPlan.nodes.resize(3);
    minPlan.nodes[0].level = 2;
    minPlan.nodes[0].supportPixels = {4};
    minPlan.nodes[1].level = 5;
    minPlan.nodes[1].supportPixels = {4, 5};
    minPlan.nodes[2].level = 6;
    minPlan.nodes[2].supportPixels = {9};

    noninc::assignLocalRecompositionParentIndices(minPlan, false);
    require(minPlan.nodes[0].parentPlanIndex == 1,
            "min-tree parent must be the smallest higher-level containing support");
    require(vectorContainsIndex(minPlan.nodes[1].childPlanIndices, 0),
            "min-tree parent must record its local child");
    require(minPlan.nodes[1].parentPlanIndex == -1,
            "min-tree parent without containing local support must remain externally attached");
    require(minPlan.nodes[2].parentPlanIndex == -1,
            "disjoint min-tree support must not receive a local parent");

    noninc::LocalRecompositionPlan maxPlan;
    maxPlan.nodes.resize(2);
    maxPlan.nodes[0].level = 7;
    maxPlan.nodes[0].supportPixels = {1};
    maxPlan.nodes[1].level = 3;
    maxPlan.nodes[1].supportPixels = {1, 2};

    noninc::assignLocalRecompositionParentIndices(maxPlan, true);
    require(maxPlan.nodes[0].parentPlanIndex == 1,
            "max-tree parent must be the smallest lower-level containing support");
    require(vectorContainsIndex(maxPlan.nodes[1].childPlanIndices, 0),
            "max-tree parent must record its local child");
}

void testReducedFixtureMaterializesEffectiveRecompositionNode() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "fixture should not fail before materialization test setup");
    }

    const NodeId nodeId = findReducedFixtureMergeNode(runner.maxTree());
    require(nodeId != InvalidNode,
            "reduced fixture merge node must be present before materialization");

    DynamicComponentTree materializationDualTree(
        runner.maxTree().reconstructionImage()->clone(),
        false,
        adjacency);

    const auto oracle = noninc::buildLocalMergeOracle(
        runner.maxTree(),
        materializationDualTree,
        adjacency,
        nodeId);
    const auto partitionPlan = noninc::buildLocalPartitionPlan(oracle, *adjacency);
    const auto recompositionPlan = noninc::buildLocalRecompositionPlan(oracle, partitionPlan);
    const auto materialization = noninc::materializeLocalRecompositionNodes(
        materializationDualTree,
        oracle,
        recompositionPlan);

    require(materialization.numCreatedNodes() == 1,
            "fixture materialization must create one effective recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "materialized local recomposition must reconstruct the oracle target image");
    require(materializationDualTree.reconstructionImage()->isEqual(oracle.updatedImage),
            "dual tree reconstruction must match the local oracle target image");

    const auto &created = materialization.createdNodes.front();
    const PixelId modifiedPixel = oracle.contraction.modifiedProperPart.front();
    require(materializationDualTree.isAlive(created.nodeId),
            "created recomposition node must be alive");
    require(materializationDualTree.getAltitude(created.nodeId) == oracle.contraction.b,
            "created recomposition node must have level b");
    require(materializationDualTree.getNodeParent(created.nodeId) == materializationDualTree.getRoot(),
            "first materialization step attaches top-level local nodes to the root boundary");
    require(materializationDualTree.getSmallestComponent(modifiedPixel) == created.nodeId,
            "modified pixel must be directly owned by the created recomposition node");
    require(vectorContainsPixel(created.movedDirectProperParts, modifiedPixel),
            "materialization result must report the moved modified pixel");
    require(!created.sourceOwnersBefore.empty(),
            "materialization result must record source owners before moving proper parts");
    require(!created.movedChildNodes.empty(),
            "materialization splice must move the stable oracle child under the created node");
    for (NodeId childNodeId : created.movedChildNodes) {
        require(materializationDualTree.isAlive(childNodeId),
                "moved stable child must remain alive");
        require(materializationDualTree.getNodeParent(childNodeId) == created.nodeId,
                "moved stable child must be attached under the created recomposition node");
    }
    require(!materialization.emptiedSourceOwners.empty(),
            "materialization must record the source owner emptied by moving C");
    require(!materialization.releasedSourceOwners.empty(),
            "materialization must release obsolete empty source owners after moving stable children");
    for (NodeId releasedOwner : materialization.releasedSourceOwners) {
        require(!materializationDualTree.isAlive(releasedOwner),
                "released source owner must no longer be alive");
    }
    requireTreeConsistency(materializationDualTree);

    const auto serialized = noninc::debugString(materialization);
    require(serialized.find("LocalRecompositionMaterializationResult") != std::string::npos,
            "materialization result must provide a deterministic debug serialization");
    require(serialized.find("reconstruction_matches_target=true") != std::string::npos,
            "materialization serialization must expose target reconstruction status");
    require(serialized.find("moved_children=") != std::string::npos,
            "materialization serialization must expose moved child nodes");
    require(serialized.find("released_source_owners=") != std::string::npos,
            "materialization serialization must expose released source owners");
}

void testReducedFixtureMaterializesStrippedEffectiveRecompositionNode() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "fixture should not fail before stripped materialization setup");
    }

    const NodeId nodeId = findReducedFixtureMergeNode(runner.maxTree());
    require(nodeId != InvalidNode,
            "reduced fixture merge node must be present before stripped materialization");

    DynamicComponentTree materializationDualTree(
        runner.maxTree().reconstructionImage()->clone(),
        false,
        adjacency);

    const auto oracle = noninc::buildLocalMergeOracle(
        runner.maxTree(),
        materializationDualTree,
        adjacency,
        nodeId);
    const auto input = noninc::makeLocalPartitionInputFromLiveTrees(
        runner.maxTree(),
        materializationDualTree,
        adjacency,
        nodeId);
    require(input.updatedImage->isEqual(oracle.updatedImage),
            "live-tree local partition input must produce the same target image as the checker oracle");
    require(input.localPixels == oracle.localPixels,
            "live-tree local partition input must identify the same local support as the checker oracle");
    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecomposition(
        materializationDualTree,
        input,
        *adjacency);
    const auto &materialization = strippedResult.materialization;

    require(strippedResult.isOracleFree(),
            "reduced stripped materialization result must remain oracle-free");
    require(strippedResult.build.recompositionPlan.nodes.size() == 1,
            "reduced stripped materialization must keep one effective recomposition node");
    const auto plannedNodeCount = std::count_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.kind == noninc::LocalCreatedRecompositionNodeKind::Planned;
        });
    require(plannedNodeCount == 1,
            "reduced stripped materialization must create one planned recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "reduced stripped materialization must reconstruct the local target image");
    require(materializationDualTree.reconstructionImage()->isEqual(oracle.updatedImage),
            "reduced stripped dual reconstruction must match the oracle target image");

    const auto plannedIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.kind == noninc::LocalCreatedRecompositionNodeKind::Planned;
        });
    require(plannedIt != materialization.createdNodes.end(),
            "reduced stripped materialization must expose the planned recomposition node");
    const auto &created = *plannedIt;
    const PixelId modifiedPixel = oracle.contraction.modifiedProperPart.front();
    require(materializationDualTree.isAlive(created.nodeId),
            "reduced stripped recomposition node must remain alive");
    require(materializationDualTree.getAltitude(created.nodeId) == oracle.contraction.b,
            "reduced stripped recomposition node must have level b");
    require(materializationDualTree.getSmallestComponent(modifiedPixel) == created.nodeId,
            "reduced stripped modified pixel must be directly owned by the created node");
    require(!created.movedChildNodes.empty(),
            "reduced stripped splice must move stable children from local pre-state snapshots");

    const auto postLocalSnapshot = noninc::buildPostUpdateLocalDualSnapshot(
        materializationDualTree,
        oracle);
    const bool shapesMatch = noninc::localDualNodeShapesMatch(
        std::span<const noninc::LocalDualNodeSnapshot>(postLocalSnapshot),
        std::span<const noninc::LocalDualNodeSnapshot>(oracle.oracleDualNodes));
    require(shapesMatch,
            noninc::localDualNodeShapeMismatchSummary(
                std::span<const noninc::LocalDualNodeSnapshot>(postLocalSnapshot),
                std::span<const noninc::LocalDualNodeSnapshot>(oracle.oracleDualNodes)));
    requireTreeConsistency(materializationDualTree);
}

void testReducedFixtureAppliesStrippedNonOracleLocalMergeUpdate() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "fixture should not fail before stripped update setup");
    }

    const NodeId nodeId = findReducedFixtureMergeNode(runner.maxTree());
    require(nodeId != InvalidNode,
            "reduced fixture merge node must be present before stripped update");

    DynamicComponentTree primalTree = runner.maxTree();
    DynamicComponentTree dualTree(
        primalTree.reconstructionImage()->clone(),
        false,
        adjacency);

    const NodeId localNodeId = findReducedFixtureMergeNode(primalTree);
    require(localNodeId != InvalidNode,
            "reduced fixture merge node must remain present in the local primal copy");

    const auto checkerOracle = noninc::buildLocalMergeOracle(
        primalTree,
        dualTree,
        adjacency,
        localNodeId);
    const auto result = noninc::applyStrippedNonOracleLocalMergeUpdate(
        primalTree,
        dualTree,
        adjacency,
        localNodeId);

    require(result.applied,
            "stripped local update must apply the reduced fixture merge");
    require(!result.skippedNoImageChange,
            "reduced fixture stripped update must perform a real local recomposition");
    require(result.inputBuilt,
            "stripped local update must build the live-tree local input");
    require(result.localPlanOracleFree && result.materializationInputOracleFree,
            "stripped local update must keep plan and materialization inputs oracle-free");
    require(result.localRecompositionApplied && result.localRecompositionMatchesTarget,
            "stripped local update must materialize the local target");
    require(result.recompositionNodeCount == 1 && result.createdPlannedNodeCount == 1,
            "reduced stripped update must create one planned recomposition node");
    require(result.createdNodeCount >= result.createdPlannedNodeCount,
            "reduced stripped update may create deterministic connector nodes");
    require(result.dualMatchesPrimalAfterMerge,
            "stripped local update must keep primal and dual reconstructions aligned");
    require(primalTree.reconstructionImage()->isEqual(checkerOracle.updatedImage),
            "stripped updated primal tree must match the checker oracle target");
    require(dualTree.reconstructionImage()->isEqual(checkerOracle.updatedImage),
            "stripped updated dual tree must match the checker oracle target");

    const auto postLocalSnapshot = noninc::buildPostUpdateLocalDualSnapshot(
        dualTree,
        checkerOracle);
    const bool shapesMatch = noninc::localDualNodeShapesMatch(
        std::span<const noninc::LocalDualNodeSnapshot>(postLocalSnapshot),
        std::span<const noninc::LocalDualNodeSnapshot>(checkerOracle.oracleDualNodes));
    require(shapesMatch,
            noninc::localDualNodeShapeMismatchSummary(
                std::span<const noninc::LocalDualNodeSnapshot>(postLocalSnapshot),
                std::span<const noninc::LocalDualNodeSnapshot>(checkerOracle.oracleDualNodes)));
    requireTreeConsistency(primalTree);
    requireTreeConsistency(dualTree);
}

void testSyntheticMaterializationCreatesLocalParentChildChain() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {5, 15, 30, 40});
    input.contraction.dualIsMaxtree = false;

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {5, 15};
    plan.nodes.resize(2);
    plan.nodes[0].level = 5;
    plan.nodes[0].supportPixels = {0};
    plan.nodes[0].directProperParts = {0};
    plan.nodes[0].parentPlanIndex = 1;
    plan.nodes[1].level = 15;
    plan.nodes[1].supportPixels = {0, 1};
    plan.nodes[1].directProperParts = {1};
    plan.nodes[1].parentPlanIndex = -1;
    plan.nodes[1].childPlanIndices = {0};

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "synthetic local chain helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 2,
            "synthetic chain must create two effective recomposition nodes");
    require(materialization.reconstructionMatchesTarget,
            "synthetic chain must reconstruct the target image");
    require(minTree.reconstructionImage()->isEqual(input.updatedImage),
            "synthetic chain reconstruction must match the stripped target");

    const auto &child = materialization.createdNodes[0];
    const auto &parent = materialization.createdNodes[1];
    require(minTree.isAlive(child.nodeId) && minTree.isAlive(parent.nodeId),
            "created chain nodes must remain alive");
    require(minTree.getAltitude(child.nodeId) == 5,
            "synthetic child node must have level 5");
    require(minTree.getAltitude(parent.nodeId) == 15,
            "synthetic parent node must have level 15");
    require(minTree.getNodeParent(child.nodeId) == parent.nodeId,
            "synthetic local child must attach to the created local parent");
    require(minTree.getSmallestComponent(0) == child.nodeId,
            "pixel 0 must move to the created child node");
    require(minTree.getSmallestComponent(1) == parent.nodeId,
            "pixel 1 must move to the created parent node");
    require(!materialization.releasedSourceOwners.empty(),
            "empty source owners in the synthetic chain should be released after stabilization");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationKeepsEmptyOwnerWithResidualChildUnreleased() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId residualChild = minTree.getSmallestComponent(0);
    const NodeId sourceOwner = minTree.getSmallestComponent(1);
    require(residualChild != InvalidNode &&
            sourceOwner != InvalidNode &&
            !minTree.isRoot(sourceOwner),
            "empty-owner with child setup must have a residual child and non-root source owner");
    require(minTree.getNodeParent(residualChild) == sourceOwner,
            "residual child must start below the source owner");
    const NodeId promotedParent = minTree.getNodeParent(sourceOwner);
    require(promotedParent != InvalidNode &&
            promotedParent != sourceOwner &&
            minTree.isAlive(promotedParent),
            "empty-owner with child setup must have a live promoted parent");

    noninc::LocalPartitionInput input;
    input.updatedImage = minTree.reconstructionImage()->clone();
    input.contraction.dualIsMaxtree = false;
    input.localPixels = {0, 1, 2};

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {0, 1};
    plan.nodes[0].directProperParts = {1};
    plan.nodes[0].preDualOwners = {promotedParent};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "empty-owner with residual child materialization must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 1,
            "empty-owner with residual child case must create one recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "empty-owner with residual child case must preserve reconstruction");
    require(vectorContainsIndex(materialization.emptiedSourceOwners, sourceOwner),
            "source owner must be recorded as emptied");
    require(!vectorContainsIndex(materialization.releasedSourceOwners, sourceOwner),
            "source owner with residual child must not be released by the simple release step");
    require(minTree.isAlive(sourceOwner),
            "source owner with residual child must remain alive outside the promoted subdomain");

    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    require(minTree.getNodeParent(residualChild) == sourceOwner,
            "residual child must remain below the unreleased source owner");
    require(minTree.getNodeParent(recomposedNode) == promotedParent,
            "recomposition node must attach to the promoted parent");
    require(minTree.getSmallestComponent(0) == residualChild,
            "residual child pixel ownership must be preserved");
    require(minTree.getSmallestComponent(1) == recomposedNode,
            "moved proper part must be owned by the recomposition node");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationReleasesPromotableEmptyOwnerWithResidualChildWhenAuthorized() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId residualChild = minTree.getSmallestComponent(0);
    const NodeId sourceOwner = minTree.getSmallestComponent(1);
    require(residualChild != InvalidNode &&
            sourceOwner != InvalidNode &&
            !minTree.isRoot(sourceOwner),
            "authorized residual release setup must have a residual child and non-root source owner");
    require(minTree.getNodeParent(residualChild) == sourceOwner,
            "authorized residual release child must start below the source owner");
    const NodeId promotedParent = minTree.getNodeParent(sourceOwner);
    require(promotedParent != InvalidNode &&
            promotedParent != sourceOwner &&
            minTree.isAlive(promotedParent),
            "authorized residual release setup must have a live promoted parent");

    noninc::LocalPartitionInput input;
    input.updatedImage = minTree.reconstructionImage()->clone();
    input.contraction.dualIsMaxtree = false;
    input.localPixels = {0, 1, 2};

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {1};
    plan.nodes[0].directProperParts = {1};
    plan.nodes[0].preDualOwners = {promotedParent};
    plan.nodes[0].parentPlanIndex = -1;

    const auto oracle = noninc::makeStrippedNonOracleMaterializationContext(
        minTree,
        input);
    noninc::LocalRecompositionMaterializationOptions options;
    options.releaseEmptySourceOwnersWithResidualChildren = true;
    options.checkNoEmptyLiveNodesAfterCanonicalization = true;
    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan,
        nullptr,
        options);

    require(materialization.reconstructionMatchesTarget,
            "authorized residual release must preserve reconstruction");
    require(vectorContainsIndex(materialization.emptiedSourceOwners, sourceOwner),
            "authorized residual release must record the emptied source owner");
    require(vectorContainsIndex(materialization.releasedSourceOwners, sourceOwner),
            "authorized residual release must release the emptied source owner");
    require(!minTree.isAlive(sourceOwner),
            "authorized residual release must remove the source owner");
    require(minTree.getNodeParent(residualChild) == promotedParent,
            "authorized residual release must promote the residual child to the source parent");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationUsesMappedNonRootExternalAttachment() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId boundaryParent = minTree.getSmallestComponent(1);
    require(boundaryParent != InvalidNode,
            "synthetic boundary parent must exist");
    require(!minTree.isRoot(boundaryParent),
            "synthetic boundary parent must be a non-root node");

    constexpr NodeId oracleParentId = 777;
    noninc::LocalMergeOracle oracle;
    oracle.contraction.dualIsMaxtree = false;
    oracle.preImage = minTree.reconstructionImage();
    oracle.updatedImage = clone_with_values(oracle.preImage, {5, 20, 30, 40});
    noninc::LocalMergePixelOwner boundaryOwner;
    boundaryOwner.pixelId = 1;
    boundaryOwner.role = noninc::LocalMergePixelRole::Boundary;
    boundaryOwner.preDualOwner = boundaryParent;
    boundaryOwner.oracleDualOwner = oracleParentId;
    boundaryOwner.preDualOwnerLevel = minTree.getAltitude(boundaryParent);
    boundaryOwner.oracleDualOwnerLevel = minTree.getAltitude(boundaryParent);
    oracle.pixelOwners.push_back(boundaryOwner);

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {5};
    plan.nodes.resize(1);
    plan.nodes[0].level = 5;
    plan.nodes[0].supportPixels = {0};
    plan.nodes[0].directProperParts = {0};
    plan.nodes[0].oracleParentId = oracleParentId;
    plan.nodes[0].parentPlanIndex = -1;

    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan);

    require(materialization.numCreatedNodes() == 1,
            "non-root external attachment case must create one node");
    require(materialization.reconstructionMatchesTarget,
            "non-root external attachment case must reconstruct the target image");
    const auto &created = materialization.createdNodes.front();
    require(created.parentNodeId == boundaryParent,
            "created node must record the mapped non-root external parent");
    require(minTree.getNodeParent(created.nodeId) == boundaryParent,
            "created node must be attached to the mapped non-root external parent");
    require(minTree.getSmallestComponent(0) == created.nodeId,
            "pixel 0 must move to the created non-root-attached node");
    require(minTree.getSmallestComponent(1) == boundaryParent,
            "boundary pixel owner must remain unchanged");
    require(!materialization.releasedSourceOwners.empty(),
            "old leaf owner should be released after moving the only direct proper part");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationCreatesLocalChainUnderMappedBoundary() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId boundaryParent = minTree.getSmallestComponent(2);
    require(boundaryParent != InvalidNode,
            "synthetic mixed boundary parent must exist");
    require(!minTree.isRoot(boundaryParent),
            "synthetic mixed boundary parent must be non-root");

    constexpr NodeId oracleParentId = 888;
    noninc::LocalMergeOracle oracle;
    oracle.contraction.dualIsMaxtree = false;
    oracle.preImage = minTree.reconstructionImage();
    oracle.updatedImage = clone_with_values(oracle.preImage, {5, 15, 30, 40});

    noninc::LocalMergePixelOwner boundaryOwner;
    boundaryOwner.pixelId = 2;
    boundaryOwner.role = noninc::LocalMergePixelRole::Boundary;
    boundaryOwner.preDualOwner = boundaryParent;
    boundaryOwner.oracleDualOwner = oracleParentId;
    boundaryOwner.preDualOwnerLevel = minTree.getAltitude(boundaryParent);
    boundaryOwner.oracleDualOwnerLevel = minTree.getAltitude(boundaryParent);
    oracle.pixelOwners.push_back(boundaryOwner);

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {5, 15};
    plan.nodes.resize(2);
    plan.nodes[0].level = 5;
    plan.nodes[0].supportPixels = {0};
    plan.nodes[0].directProperParts = {0};
    plan.nodes[0].parentPlanIndex = 1;
    plan.nodes[1].level = 15;
    plan.nodes[1].supportPixels = {0, 1};
    plan.nodes[1].directProperParts = {1};
    plan.nodes[1].oracleParentId = oracleParentId;
    plan.nodes[1].parentPlanIndex = -1;
    plan.nodes[1].childPlanIndices = {0};

    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan);

    require(materialization.numCreatedNodes() == 2,
            "mixed local/external splice must create two local nodes");
    require(materialization.reconstructionMatchesTarget,
            "mixed local/external splice must reconstruct the target image");

    const auto &child = materialization.createdNodes[0];
    const auto &parent = materialization.createdNodes[1];
    require(minTree.getNodeParent(child.nodeId) == parent.nodeId,
            "mixed splice child must attach to the created local parent");
    require(minTree.getNodeParent(parent.nodeId) == boundaryParent,
            "mixed splice parent must attach to the mapped external boundary");
    require(parent.parentNodeId == boundaryParent,
            "materialization result must report the mapped external boundary");
    require(minTree.getSmallestComponent(0) == child.nodeId,
            "pixel 0 must move to the mixed splice child");
    require(minTree.getSmallestComponent(1) == parent.nodeId,
            "pixel 1 must move to the mixed splice parent");
    require(minTree.getSmallestComponent(2) == boundaryParent,
            "boundary pixel must remain owned by the mapped external parent");
    require(!materialization.releasedSourceOwners.empty(),
            "mixed splice should release obsolete source owners");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationMovesMultipleStableChildren() {
    auto image = make_split_child_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId firstStableChild = minTree.getSmallestComponent(0);
    const NodeId sharedLowOwner = minTree.getSmallestComponent(1);
    require(firstStableChild != InvalidNode,
            "synthetic stable child owner must exist");
    require(sharedLowOwner == firstStableChild,
            "test setup expects the first two pixels to start in the same low component");
    require(minTree.getAltitude(firstStableChild) == 10,
            "synthetic low component must have altitude 10");

    const NodeId oldSourceOwner = minTree.getSmallestComponent(2);
    require(oldSourceOwner != InvalidNode && !minTree.isRoot(oldSourceOwner),
            "synthetic source owner must be a live non-root node");
    require(minTree.getNodeParent(firstStableChild) == oldSourceOwner,
            "synthetic low component must initially be attached under the source owner");

    const NodeId secondStableChild = minTree.createDetachedNode(minTree.getAltitude(firstStableChild));
    require(secondStableChild != InvalidNode,
            "synthetic second stable child must be creatable");
    minTree.moveProperPart(secondStableChild, firstStableChild, 1);
    minTree.attachNode(oldSourceOwner, secondStableChild);

    require(minTree.getSmallestComponent(0) == firstStableChild,
            "pixel 0 must remain owned by the first stable child");
    require(minTree.getSmallestComponent(1) == secondStableChild,
            "pixel 1 must be owned by the second stable child after setup split");
    require(minTree.getNodeParent(firstStableChild) == oldSourceOwner,
            "first stable child must remain under the source owner before recomposition");
    require(minTree.getNodeParent(secondStableChild) == oldSourceOwner,
            "second stable child must attach under the source owner before recomposition");
    requireTreeConsistency(minTree);

    constexpr NodeId oracleFirstChildId = 501;
    constexpr NodeId oracleSecondChildId = 502;
    constexpr NodeId oracleParentId = 900;

    noninc::LocalMergeOracle oracle;
    oracle.contraction.dualIsMaxtree = false;
    oracle.preImage = minTree.reconstructionImage();
    oracle.updatedImage = clone_with_values(oracle.preImage, {10, 10, 20, 40});

    auto addOwner = [&](PixelId pixelId,
                        noninc::LocalMergePixelRole role,
                        NodeId preOwner,
                        NodeId oracleOwner) {
        noninc::LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = role;
        owner.preDualOwner = preOwner;
        owner.oracleDualOwner = oracleOwner;
        owner.preDualOwnerLevel =
            preOwner == InvalidNode || !minTree.isAlive(preOwner) ? -1 : minTree.getAltitude(preOwner);
        owner.oracleDualOwnerLevel = owner.preDualOwnerLevel;
        oracle.pixelOwners.push_back(owner);
    };

    addOwner(0, noninc::LocalMergePixelRole::Preserved, firstStableChild, oracleFirstChildId);
    addOwner(1, noninc::LocalMergePixelRole::Preserved, secondStableChild, oracleSecondChildId);
    addOwner(3, noninc::LocalMergePixelRole::Boundary, minTree.getRoot(), oracleParentId);

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {0, 1, 2};
    plan.nodes[0].directProperParts = {2};
    plan.nodes[0].oracleParentId = oracleParentId;
    plan.nodes[0].oracleChildrenInSnapshot = {oracleFirstChildId, oracleSecondChildId};
    plan.nodes[0].parentPlanIndex = -1;

    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan);

    require(materialization.numCreatedNodes() == 1,
            "multiple-child splice must create one recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "multiple-child splice must reconstruct the target image");

    const auto &created = materialization.createdNodes.front();
    require(minTree.isAlive(created.nodeId),
            "created multiple-child recomposition node must remain alive");
    require(minTree.getAltitude(created.nodeId) == 20,
            "created multiple-child recomposition node must have level 20");
    require(created.movedChildNodes.size() == 2,
            "materialization must report both moved stable children");
    require(vectorContainsIndex(created.movedChildNodes, firstStableChild),
            "first stable child must be reported as moved");
    require(vectorContainsIndex(created.movedChildNodes, secondStableChild),
            "second stable child must be reported as moved");
    require(minTree.getNodeParent(firstStableChild) == created.nodeId,
            "first stable child must be spliced under the created recomposition node");
    require(minTree.getNodeParent(secondStableChild) == created.nodeId,
            "second stable child must be spliced under the created recomposition node");
    require(minTree.getSmallestComponent(0) == firstStableChild,
            "pixel 0 owner must remain the first stable child");
    require(minTree.getSmallestComponent(1) == secondStableChild,
            "pixel 1 owner must remain the second stable child");
    require(minTree.getSmallestComponent(2) == created.nodeId,
            "pixel 2 must move to the created recomposition node");
    require(vectorContainsIndex(materialization.releasedSourceOwners, oldSourceOwner),
            "obsolete source owner must be released after both stable children are moved");
    require(!minTree.isAlive(oldSourceOwner),
            "released source owner must no longer be alive");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationMovesStableChildrenFromNonOraclePlan() {
    auto image = make_split_child_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId firstStableChild = minTree.getSmallestComponent(0);
    const NodeId oldSourceOwner = minTree.getSmallestComponent(2);
    require(firstStableChild != InvalidNode &&
            oldSourceOwner != InvalidNode &&
            !minTree.isRoot(oldSourceOwner),
            "non-oracle materialization setup must have a stable child and source owner");

    const NodeId secondStableChild = minTree.createDetachedNode(minTree.getAltitude(firstStableChild));
    require(secondStableChild != InvalidNode,
            "non-oracle materialization setup must create a second stable child");
    minTree.moveProperPart(secondStableChild, firstStableChild, 1);
    minTree.attachNode(oldSourceOwner, secondStableChild);

    require(minTree.getSmallestComponent(0) == firstStableChild,
            "pixel 0 must be owned by the first stable child before non-oracle materialization");
    require(minTree.getSmallestComponent(1) == secondStableChild,
            "pixel 1 must be owned by the second stable child before non-oracle materialization");
    requireTreeConsistency(minTree);

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 10, 20, 40});
    input.contraction.dualIsMaxtree = false;

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {0, 1, 2};
    plan.nodes[0].directProperParts = {2};
    plan.nodes[0].stableChildPreNodes = {firstStableChild, secondStableChild};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "non-oracle stable-child explicit helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 1,
            "non-oracle stable-child materialization must create one recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle stable-child materialization must reconstruct the target image");

    const auto &created = materialization.createdNodes.front();
    require(created.movedChildNodes.size() == 2,
            "non-oracle materialization must report both stable children");
    require(vectorContainsIndex(created.movedChildNodes, firstStableChild),
            "non-oracle materialization must report the first stable child");
    require(vectorContainsIndex(created.movedChildNodes, secondStableChild),
            "non-oracle materialization must report the second stable child");
    require(minTree.getNodeParent(firstStableChild) == created.nodeId,
            "first stable child must be moved under the non-oracle recomposition node");
    require(minTree.getNodeParent(secondStableChild) == created.nodeId,
            "second stable child must be moved under the non-oracle recomposition node");
    require(minTree.getSmallestComponent(0) == firstStableChild,
            "pixel 0 owner must stay on the first stable child");
    require(minTree.getSmallestComponent(1) == secondStableChild,
            "pixel 1 owner must stay on the second stable child");
    require(minTree.getSmallestComponent(2) == created.nodeId,
            "pixel 2 must move to the non-oracle recomposition node");
    require(vectorContainsIndex(materialization.releasedSourceOwners, oldSourceOwner),
            "obsolete source owner must be released after non-oracle stable-child moves");
    require(!minTree.isAlive(oldSourceOwner),
            "released source owner must no longer be alive after non-oracle materialization");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationCreatesNonOracleBridgeBelowMovedStableChild() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId stableChild = minTree.getSmallestComponent(1);
    const NodeId sourceOwner = minTree.getSmallestComponent(2);
    require(lowNode != InvalidNode &&
            stableChild != InvalidNode &&
            sourceOwner != InvalidNode,
            "non-oracle bridge setup must have low, stable, and source nodes");
    require(minTree.getAltitude(lowNode) == 10 &&
            minTree.getAltitude(stableChild) == 20 &&
            minTree.getAltitude(sourceOwner) == 30,
            "non-oracle bridge setup expects levels 10, 20, and 30");
    require(minTree.getNodeParent(lowNode) == stableChild &&
            minTree.getNodeParent(stableChild) == sourceOwner,
            "non-oracle bridge setup expects a nested low child under the stable child");

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = stableChild;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot preStable;
    preStable.nodeId = stableChild;
    preStable.parentId = sourceOwner;
    preStable.parentInSnapshot = false;
    preStable.altitude = 20;
    preStable.directLocalProperParts = {1};
    preStable.supportLocalPixels = {0, 1};
    preStable.childrenInSnapshot = {lowNode};
    preStable.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot preSource;
    preSource.nodeId = sourceOwner;
    preSource.parentId = InvalidNode;
    preSource.parentInSnapshot = false;
    preSource.altitude = 30;
    preSource.directLocalProperParts = {2};
    preSource.supportLocalPixels = {0, 1, 2};
    preSource.childrenInSnapshot = {stableChild};
    preSource.supportTouchesModified = true;
    preSource.supportTouchesPreserved = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 20, 25, 40});
    input.contraction.dualIsMaxtree = false;
    input.contraction.a = 30;
    input.contraction.b = 25;
    input.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    input.localPixels = {0, 1, 2};
    input.preDualNodes = {preLow, preStable, preSource};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Preserved, lowNode},
        {1, noninc::LocalMergePixelRole::Preserved, stableChild},
        {2, noninc::LocalMergePixelRole::Modified, sourceOwner},
    };

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecomposition(
        minTree,
        input,
        *adjacency);
    require(strippedResult.isOracleFree(),
            "non-oracle bridge helper result must remain oracle-free");
    const auto &plan = strippedResult.build.recompositionPlan;
    require(plan.nodes.size() == 1 &&
            plan.nodes.front().stableChildPreNodes == std::vector<NodeId>{stableChild},
            "non-oracle bridge helper must derive the moved stable child from pre-state snapshots");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 1,
            "non-oracle bridge case must create only the recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle bridge case must reconstruct the target image");
    require(materialization.canonicalizationBreakdown.emptyNodeReleases == 0 &&
                materialization.canonicalizationBreakdown.sameLevelParentChildCollapses == 0 &&
                materialization.canonicalizationBreakdown.sameLevelSiblingMerges == 0 &&
                materialization.canonicalizationBreakdown.comparableSiblingMoves == 0,
            "non-oracle bridge case must not require corrective canonicalization");

    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    const auto bridgeIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.planIndex == -1 && created.level == 15;
        });
    require(bridgeIt == materialization.createdNodes.end(),
            "non-oracle bridge case must omit the empty intermediate bridge");

    require(minTree.getNodeParent(stableChild) == recomposedNode,
            "stable child must move under the non-oracle recomposed node");
    require(minTree.getNodeParent(lowNode) == stableChild,
            "low descendant must remain directly below the moved stable child");
    require(minTree.getSmallestComponent(0) == lowNode,
            "pixel 0 must remain owned by the low descendant");
    require(minTree.getSmallestComponent(1) == stableChild,
            "pixel 1 must remain owned by the moved stable child");
    require(minTree.getSmallestComponent(2) == recomposedNode,
            "pixel 2 must move to the recomposed node");
    require(noninc::supportContainsPixel(minTree, recomposedNode, 0),
            "non-oracle recomposed support must include the bridged descendant");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "non-oracle bridge canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationRepairsNonOracleVisiblePreParentLink() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId stableParent = minTree.getSmallestComponent(1);
    const NodeId sourceOwner = minTree.getSmallestComponent(2);
    require(lowNode != InvalidNode &&
            stableParent != InvalidNode &&
            sourceOwner != InvalidNode,
            "non-oracle parent-link setup must have low, stable parent, and source nodes");
    require(minTree.getAltitude(lowNode) == 10 &&
            minTree.getAltitude(stableParent) == 20 &&
            minTree.getAltitude(sourceOwner) == 30,
            "non-oracle parent-link setup expects levels 10, 20, and 30");
    require(minTree.getNodeParent(lowNode) == stableParent &&
            minTree.getNodeParent(stableParent) == sourceOwner,
            "non-oracle parent-link setup expects low under the stable parent");

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = stableParent;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot preStable;
    preStable.nodeId = stableParent;
    preStable.parentId = sourceOwner;
    preStable.parentInSnapshot = true;
    preStable.altitude = 20;
    preStable.directLocalProperParts = {1};
    preStable.supportLocalPixels = {0, 1};
    preStable.childrenInSnapshot = {lowNode};
    preStable.supportTouchesPreserved = true;

    noninc::LocalDualNodeSnapshot preSource;
    preSource.nodeId = sourceOwner;
    preSource.parentId = InvalidNode;
    preSource.parentInSnapshot = false;
    preSource.altitude = 30;
    preSource.directLocalProperParts = {2};
    preSource.supportLocalPixels = {0, 1, 2};
    preSource.childrenInSnapshot = {stableParent};
    preSource.supportTouchesModified = true;
    preSource.supportTouchesPreserved = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 20, 25, 40});
    input.contraction.dualIsMaxtree = false;
    input.contraction.a = 30;
    input.contraction.b = 25;
    input.contraction.direction = ::debug::region::RegionUpdateDirection::Lower;
    input.localPixels = {0, 1, 2};
    input.preDualNodes = {preLow, preStable, preSource};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Preserved, lowNode},
        {1, noninc::LocalMergePixelRole::Preserved, stableParent},
        {2, noninc::LocalMergePixelRole::Modified, sourceOwner},
    };

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {25};
    plan.nodes.resize(1);
    plan.nodes[0].level = 25;
    plan.nodes[0].supportPixels = {0, 1, 2};
    plan.nodes[0].directProperParts = {2};
    plan.nodes[0].stableChildPreNodes = {lowNode, stableParent};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "non-oracle parent-link helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 1,
            "non-oracle parent-link case must only create the recomposition node");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle parent-link repair must preserve reconstruction");
    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    require(minTree.getNodeParent(stableParent) == recomposedNode,
            "stable parent must remain under the recomposed node");
    require(minTree.getNodeParent(lowNode) == stableParent,
            "local pre-state enforcement must restore the visible low-parent link");
    require(minTree.getSmallestComponent(0) == lowNode,
            "pixel 0 must remain owned by the low child");
    require(minTree.getSmallestComponent(1) == stableParent,
            "pixel 1 must remain owned by the stable parent");
    require(minTree.getSmallestComponent(2) == recomposedNode,
            "pixel 2 must move to the recomposed node");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationCreatesNonOracleMissingParentBridgeAboveLiveChild() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId sourceOwner = minTree.getSmallestComponent(1);
    const NodeId contextParent = minTree.getSmallestComponent(2);
    require(lowNode != InvalidNode &&
            sourceOwner != InvalidNode &&
            contextParent != InvalidNode,
            "non-oracle missing-parent setup must have a low child, source owner, and context parent");
    require(minTree.getAltitude(lowNode) == 10 &&
            minTree.getAltitude(sourceOwner) == 20 &&
            minTree.getAltitude(contextParent) == 30,
            "non-oracle missing-parent setup expects levels 10, 20, and 30");
    require(minTree.getNodeParent(lowNode) == sourceOwner &&
            minTree.getNodeParent(sourceOwner) == contextParent,
            "non-oracle missing-parent setup expects the bridge to be absent initially");

    constexpr NodeId missingParentId = 12015;

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = missingParentId;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot missingParent;
    missingParent.nodeId = missingParentId;
    missingParent.parentId = contextParent;
    missingParent.parentInSnapshot = true;
    missingParent.altitude = 15;
    missingParent.supportLocalPixels = {0};
    missingParent.childrenInSnapshot = {lowNode};
    missingParent.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preContext;
    preContext.nodeId = contextParent;
    preContext.parentId = InvalidNode;
    preContext.parentInSnapshot = false;
    preContext.altitude = 30;
    preContext.directLocalProperParts = {2};
    preContext.supportLocalPixels = {0, 1, 2};
    preContext.childrenInSnapshot = {missingParentId};
    preContext.supportTouchesModified = true;
    preContext.supportTouchesBoundary = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {10, 25, 30, 40});
    input.contraction.dualIsMaxtree = false;
    input.localPixels = {0, 1, 2};
    input.preDualNodes = {preLow, missingParent, preContext};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Boundary, lowNode},
        {1, noninc::LocalMergePixelRole::Modified, sourceOwner},
        {2, noninc::LocalMergePixelRole::Boundary, contextParent},
    };

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {25};
    plan.nodes.resize(1);
    plan.nodes[0].level = 25;
    plan.nodes[0].supportPixels = {0, 1, 2};
    plan.nodes[0].directProperParts = {1};
    plan.nodes[0].preDualOwners = {contextParent};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "non-oracle missing-parent helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 2,
            "non-oracle missing-parent case must create the recomposition node and one bridge");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle missing-parent bridge must preserve reconstruction");
    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    const auto bridgeIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.planIndex == -1 && created.level == 15;
        });
    require(bridgeIt != materialization.createdNodes.end(),
            "non-oracle missing parent must be materialized at its pre-snapshot level");
    const NodeId bridgeNode = bridgeIt->nodeId;
    require(!minTree.isAlive(bridgeNode),
            "empty non-oracle missing-parent bridge must be released by local canonicalization");
    require(minTree.getNodeParent(lowNode) == contextParent,
            "live child must be promoted after the empty non-oracle missing-parent bridge is released");
    require(minTree.getNodeParent(recomposedNode) == contextParent,
            "recomposition node must stay attached under the context parent");
    require(minTree.getSmallestComponent(0) == lowNode,
            "boundary pixel must remain owned by the low child");
    require(minTree.getSmallestComponent(1) == recomposedNode,
            "modified context pixel must move to the recomposition node");
    require(minTree.getSmallestComponent(2) == contextParent,
            "context proper part must remain owned by the context parent");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "non-oracle missing-parent canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationCreatesNonOracleMissingChildBridgeBelowLiveParent() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId oldParent = minTree.getSmallestComponent(1);
    require(lowNode != InvalidNode && oldParent != InvalidNode,
            "non-oracle missing-child setup must have a low child and old parent");
    require(minTree.getAltitude(lowNode) == 10 &&
            minTree.getAltitude(oldParent) == 20,
            "non-oracle missing-child setup expects levels 10 and 20");
    require(minTree.getNodeParent(lowNode) == oldParent,
            "non-oracle missing-child setup expects the child bridge to be absent initially");

    constexpr NodeId missingChildId = 13015;

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = missingChildId;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot missingChild;
    missingChild.nodeId = missingChildId;
    missingChild.parentId = oldParent;
    missingChild.parentInSnapshot = true;
    missingChild.altitude = 15;
    missingChild.supportLocalPixels = {0};
    missingChild.childrenInSnapshot = {lowNode};
    missingChild.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preParent;
    preParent.nodeId = oldParent;
    preParent.parentId = InvalidNode;
    preParent.parentInSnapshot = false;
    preParent.altitude = 20;
    preParent.directLocalProperParts = {1};
    preParent.supportLocalPixels = {0, 1};
    preParent.childrenInSnapshot = {missingChildId};
    preParent.supportTouchesModified = true;
    preParent.supportTouchesBoundary = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = minTree.reconstructionImage()->clone();
    input.contraction.dualIsMaxtree = false;
    input.localPixels = {0, 1};
    input.preDualNodes = {preLow, missingChild, preParent};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Boundary, lowNode},
        {1, noninc::LocalMergePixelRole::Modified, oldParent},
    };

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {0, 1};
    plan.nodes[0].directProperParts = {1};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "non-oracle missing-child helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 2,
            "non-oracle missing-child case must create the recomposition node and one child bridge");
    require(materialization.nonOracleChildLinkBridgeAttemptCount > 0,
            "non-oracle missing-child case must attempt a child-link bridge");
    require(materialization.nonOracleChildLinkBridgeMaterializableCount == 1,
            "non-oracle missing-child case must materialize exactly one child-link bridge");
    require(materialization.nonOracleChildLinkBridgeInvalidLevelCount == 0,
            "non-oracle missing-child child-link bridge must have valid levels");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle missing-child bridge must preserve reconstruction");
    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    const auto bridgeIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.planIndex == -1 &&
                   created.level == 15 &&
                   created.kind == noninc::LocalCreatedRecompositionNodeKind::NonOracleChildLinkBridge;
        });
    require(bridgeIt != materialization.createdNodes.end(),
            "non-oracle missing child must be materialized at its pre-snapshot level");
    const NodeId bridgeNode = bridgeIt->nodeId;
    require(!minTree.isAlive(bridgeNode),
            "empty non-oracle missing-child bridge must be released by local canonicalization");
    require(minTree.getNodeParent(lowNode) == recomposedNode,
            "live descendant must be promoted after the empty non-oracle child bridge is released");
    require(minTree.getSmallestComponent(0) == lowNode,
            "boundary pixel must remain owned by the low child");
    require(minTree.getSmallestComponent(1) == recomposedNode,
            "parent pixel must move to the recomposition node");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "non-oracle missing-child canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationCreatesBridgeBelowMovedStableChild() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId stableChild = minTree.getSmallestComponent(1);
    const NodeId sourceOwner = minTree.getSmallestComponent(2);
    require(lowNode != InvalidNode &&
            stableChild != InvalidNode &&
            sourceOwner != InvalidNode,
            "bridge-below-stable setup must have low, stable, and source nodes");
    require(minTree.getAltitude(lowNode) == 10 &&
            minTree.getAltitude(stableChild) == 20 &&
            minTree.getAltitude(sourceOwner) == 30,
            "bridge-below-stable setup expects levels 10, 20, and 30");
    require(minTree.getNodeParent(lowNode) == stableChild &&
            minTree.getNodeParent(stableChild) == sourceOwner,
            "bridge-below-stable setup expects a nested low child under the stable child");

    constexpr NodeId oracleLowId = 1010;
    constexpr NodeId oracleBridgeId = 1015;
    constexpr NodeId oracleStableId = 1020;
    constexpr NodeId oracleParentId = 1025;

    noninc::LocalMergeOracle oracle;
    oracle.contraction.dualIsMaxtree = false;
    oracle.preImage = minTree.reconstructionImage();
    oracle.updatedImage = clone_with_values(oracle.preImage, {10, 20, 25, 40});

    auto addOwner = [&](PixelId pixelId,
                        noninc::LocalMergePixelRole role,
                        NodeId preOwner,
                        NodeId oracleOwner,
                        int oracleOwnerLevel) {
        noninc::LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = role;
        owner.preDualOwner = preOwner;
        owner.oracleDualOwner = oracleOwner;
        owner.preDualOwnerLevel = minTree.getAltitude(preOwner);
        owner.oracleDualOwnerLevel = oracleOwnerLevel;
        oracle.pixelOwners.push_back(owner);
    };

    addOwner(0, noninc::LocalMergePixelRole::Preserved, lowNode, oracleLowId, 10);
    addOwner(1, noninc::LocalMergePixelRole::Preserved, stableChild, oracleStableId, 20);
    addOwner(2, noninc::LocalMergePixelRole::Modified, sourceOwner, oracleParentId, 25);

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = stableChild;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};

    noninc::LocalDualNodeSnapshot preStable;
    preStable.nodeId = stableChild;
    preStable.parentId = sourceOwner;
    preStable.parentInSnapshot = false;
    preStable.altitude = 20;
    preStable.directLocalProperParts = {1};
    preStable.supportLocalPixels = {0, 1};
    preStable.childrenInSnapshot = {lowNode};

    noninc::LocalDualNodeSnapshot oracleLow = preLow;
    oracleLow.nodeId = oracleLowId;
    oracleLow.parentId = oracleBridgeId;

    noninc::LocalDualNodeSnapshot oracleBridge;
    oracleBridge.nodeId = oracleBridgeId;
    oracleBridge.parentId = oracleStableId;
    oracleBridge.parentInSnapshot = true;
    oracleBridge.altitude = 15;
    oracleBridge.supportLocalPixels = {0};
    oracleBridge.childrenInSnapshot = {oracleLowId};

    noninc::LocalDualNodeSnapshot oracleStable = preStable;
    oracleStable.nodeId = oracleStableId;
    oracleStable.parentId = oracleParentId;
    oracleStable.parentInSnapshot = true;
    oracleStable.childrenInSnapshot = {oracleBridgeId};

    noninc::LocalDualNodeSnapshot oracleParent;
    oracleParent.nodeId = oracleParentId;
    oracleParent.parentId = InvalidNode;
    oracleParent.parentInSnapshot = false;
    oracleParent.altitude = 25;
    oracleParent.directLocalProperParts = {2};
    oracleParent.supportLocalPixels = {0, 1, 2};
    oracleParent.childrenInSnapshot = {oracleStableId};
    oracleParent.supportTouchesModified = true;

    oracle.preDualNodes = {preLow, preStable};
    oracle.oracleDualNodes = {oracleLow, oracleBridge, oracleStable, oracleParent};

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {25};
    plan.nodes.resize(1);
    plan.nodes[0].level = 25;
    plan.nodes[0].supportPixels = {0, 1, 2};
    plan.nodes[0].directProperParts = {2};
    plan.nodes[0].oracleNodeId = oracleParentId;
    plan.nodes[0].oracleParentId = InvalidNode;
    plan.nodes[0].oracleChildrenInSnapshot = {oracleStableId};
    plan.nodes[0].parentPlanIndex = -1;

    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan);

    require(materialization.numCreatedNodes() == 2,
            "bridge-below-stable case must create the recomposition node and one oracle bridge");
    require(materialization.reconstructionMatchesTarget,
            "bridge-below-stable case must reconstruct the oracle target");
    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    const auto bridgeIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.planIndex == -1 && created.level == 15;
        });
    require(bridgeIt != materialization.createdNodes.end(),
            "missing oracle descendant must be materialized as a level-15 bridge");
    const NodeId bridgeNode = bridgeIt->nodeId;

    require(minTree.getNodeParent(stableChild) == recomposedNode,
            "stable child must move under the recomposed node");
    require(!minTree.isAlive(bridgeNode),
            "empty missing oracle bridge must be released by local canonicalization");
    require(minTree.getNodeParent(lowNode) == stableChild,
            "low descendant must be promoted after the empty missing oracle bridge is released");
    require(minTree.getSmallestComponent(0) == lowNode,
            "pixel 0 must remain owned by the stable low descendant");
    require(minTree.getSmallestComponent(1) == stableChild,
            "pixel 1 must remain owned by the moved stable child");
    require(minTree.getSmallestComponent(2) == recomposedNode,
            "pixel 2 must move to the recomposed node");
    require(noninc::supportContainsPixel(minTree, recomposedNode, 0),
            "recomposed support must include the bridged stable descendant");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "oracle bridge canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationUsesHiddenConnectorForDetachedSupportBranch() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId oldTarget = minTree.getSmallestComponent(1);
    require(lowNode != InvalidNode && oldTarget != InvalidNode,
            "hidden-connector setup must have low and target owners");
    require(minTree.getAltitude(lowNode) == 10 && minTree.getAltitude(oldTarget) == 20,
            "hidden-connector setup expects levels 10 and 20");

    const NodeId detachedBranch = minTree.createDetachedNode(15);
    require(detachedBranch != InvalidNode,
            "hidden-connector setup must create the detached branch root");
    minTree.attachNode(oldTarget, detachedBranch);
    minTree.moveNode(lowNode, detachedBranch);
    require(minTree.getNodeParent(lowNode) == detachedBranch,
            "low owner must be under the detached branch before recomposition");
    require(minTree.getNodeParent(detachedBranch) == oldTarget,
            "detached branch initially has a visible target parent");
    requireTreeConsistency(minTree);

    constexpr NodeId oracleLowId = 100;
    constexpr NodeId oracleDetachedId = 150;
    constexpr NodeId oracleTargetId = 200;

    noninc::LocalMergeOracle oracle;
    oracle.contraction.dualIsMaxtree = false;
    oracle.preImage = minTree.reconstructionImage();
    oracle.updatedImage = oracle.preImage->clone();

    auto addOwner = [&](PixelId pixelId,
                        noninc::LocalMergePixelRole role,
                        NodeId preOwner,
                        NodeId oracleOwner) {
        noninc::LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = role;
        owner.preDualOwner = preOwner;
        owner.oracleDualOwner = oracleOwner;
        owner.preDualOwnerLevel = minTree.getAltitude(preOwner);
        owner.oracleDualOwnerLevel =
            oracleOwner == oracleLowId ? 10 :
            oracleOwner == oracleDetachedId ? 15 : 20;
        oracle.pixelOwners.push_back(owner);
    };

    addOwner(0, noninc::LocalMergePixelRole::Boundary, lowNode, oracleLowId);
    addOwner(1, noninc::LocalMergePixelRole::Modified, oldTarget, oracleTargetId);

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = detachedBranch;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preDetached;
    preDetached.nodeId = detachedBranch;
    preDetached.parentId = oldTarget;
    preDetached.parentInSnapshot = false;
    preDetached.altitude = 15;
    preDetached.supportLocalPixels = {0};
    preDetached.childrenInSnapshot = {lowNode};
    preDetached.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot oracleLow = preLow;
    oracleLow.nodeId = oracleLowId;
    oracleLow.parentId = oracleDetachedId;

    noninc::LocalDualNodeSnapshot oracleDetached = preDetached;
    oracleDetached.nodeId = oracleDetachedId;
    oracleDetached.parentId = InvalidNode;
    oracleDetached.childrenInSnapshot = {oracleLowId};

    noninc::LocalDualNodeSnapshot oracleTarget;
    oracleTarget.nodeId = oracleTargetId;
    oracleTarget.parentId = InvalidNode;
    oracleTarget.parentInSnapshot = false;
    oracleTarget.altitude = 20;
    oracleTarget.directLocalProperParts = {1};
    oracleTarget.supportLocalPixels = {0, 1};
    oracleTarget.supportTouchesModified = true;
    oracleTarget.supportTouchesBoundary = true;

    oracle.preDualNodes = {preLow, preDetached};
    oracle.oracleDualNodes = {oracleLow, oracleDetached, oracleTarget};

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {0, 1};
    plan.nodes[0].directProperParts = {1};
    plan.nodes[0].oracleNodeId = oracleTargetId;
    plan.nodes[0].oracleParentId = InvalidNode;
    plan.nodes[0].parentPlanIndex = -1;

    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan);

    require(materialization.reconstructionMatchesTarget,
            "hidden connector materialization must preserve reconstruction");
    const NodeId createdTarget = materialization.createdNodes.front().nodeId;
    require(noninc::supportContainsPixel(minTree, createdTarget, 0),
            "created target must contain the detached boundary branch in its support");
    require(!minTree.isAlive(detachedBranch),
            "empty detached branch must be released by local canonicalization");
    require(minTree.getNodeParent(lowNode) == createdTarget,
            "low owner must be promoted under the recomposed target after empty connector release");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "hidden connector canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationUsesNonOracleHiddenSupportConnectorForDetachedBranch() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId lowNode = minTree.getSmallestComponent(0);
    const NodeId oldTarget = minTree.getSmallestComponent(1);
    require(lowNode != InvalidNode && oldTarget != InvalidNode,
            "non-oracle hidden-support setup must have low and target owners");
    require(minTree.getAltitude(lowNode) == 10 && minTree.getAltitude(oldTarget) == 20,
            "non-oracle hidden-support setup expects levels 10 and 20");

    const NodeId detachedBranch = minTree.createDetachedNode(15);
    require(detachedBranch != InvalidNode,
            "non-oracle hidden-support setup must create the detached branch root");
    minTree.attachNode(oldTarget, detachedBranch);
    minTree.moveNode(lowNode, detachedBranch);
    require(minTree.getNodeParent(lowNode) == detachedBranch,
            "low owner must be under the detached branch before non-oracle recomposition");
    require(minTree.getNodeParent(detachedBranch) == oldTarget,
            "detached branch initially has the visible target as parent");
    requireTreeConsistency(minTree);

    noninc::LocalDualNodeSnapshot preLow;
    preLow.nodeId = lowNode;
    preLow.parentId = detachedBranch;
    preLow.parentInSnapshot = true;
    preLow.altitude = 10;
    preLow.directLocalProperParts = {0};
    preLow.supportLocalPixels = {0};
    preLow.supportTouchesBoundary = true;

    noninc::LocalDualNodeSnapshot preDetached;
    preDetached.nodeId = detachedBranch;
    preDetached.parentId = oldTarget;
    preDetached.parentInSnapshot = false;
    preDetached.altitude = 15;
    preDetached.supportLocalPixels = {0};
    preDetached.childrenInSnapshot = {lowNode};
    preDetached.supportTouchesBoundary = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = minTree.reconstructionImage()->clone();
    input.contraction.dualIsMaxtree = false;
    input.localPixels = {0, 1};
    input.preDualNodes = {preLow, preDetached};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Boundary, lowNode},
        {1, noninc::LocalMergePixelRole::Modified, oldTarget},
    };

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {20};
    plan.nodes.resize(1);
    plan.nodes[0].level = 20;
    plan.nodes[0].supportPixels = {0, 1};
    plan.nodes[0].directProperParts = {1};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "non-oracle hidden-support helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 2,
            "non-oracle hidden-support case must create the target and one connector");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle hidden-support connector must preserve reconstruction");
    const NodeId createdTarget = materialization.createdNodes.front().nodeId;
    require(noninc::supportContainsPixel(minTree, createdTarget, 0),
            "non-oracle created target must contain the detached branch in its support");
    require(!minTree.isAlive(detachedBranch),
            "empty non-oracle detached branch must be released by local canonicalization");
    require(minTree.getNodeParent(lowNode) == createdTarget,
            "boundary descendant must be promoted under the non-oracle recomposition node");
    require(minTree.getSmallestComponent(0) == lowNode,
            "boundary pixel must remain owned by the low detached descendant");
    require(minTree.getSmallestComponent(1) == createdTarget,
            "target pixel must move to the non-oracle recomposition node");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "non-oracle hidden-support canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationUsesHiddenExternalParentForSnapshotExternalOracleParent() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId sourceOwner = minTree.getSmallestComponent(0);
    const NodeId visibleExternalParent = minTree.getSmallestComponent(1);
    require(sourceOwner != InvalidNode && visibleExternalParent != InvalidNode,
            "hidden external-parent setup must have source and boundary owners");
    require(minTree.getAltitude(sourceOwner) == 10 &&
            minTree.getAltitude(visibleExternalParent) == 20,
            "hidden external-parent setup expects levels 10 and 20");

    constexpr NodeId oracleNodeId = 700;
    constexpr NodeId oracleParentOutsideSnapshot = 900;

    noninc::LocalMergeOracle oracle;
    oracle.contraction.dualIsMaxtree = false;
    oracle.preImage = minTree.reconstructionImage();
    oracle.updatedImage = clone_with_values(oracle.preImage, {15, 20, 30, 40});

    noninc::LocalMergePixelOwner modifiedOwner;
    modifiedOwner.pixelId = 0;
    modifiedOwner.role = noninc::LocalMergePixelRole::Modified;
    modifiedOwner.preDualOwner = sourceOwner;
    modifiedOwner.oracleDualOwner = oracleNodeId;
    modifiedOwner.preDualOwnerLevel = minTree.getAltitude(sourceOwner);
    modifiedOwner.oracleDualOwnerLevel = 15;
    oracle.pixelOwners.push_back(modifiedOwner);

    noninc::LocalMergePixelOwner boundaryOwner;
    boundaryOwner.pixelId = 1;
    boundaryOwner.role = noninc::LocalMergePixelRole::Boundary;
    boundaryOwner.preDualOwner = visibleExternalParent;
    boundaryOwner.oracleDualOwner = oracleParentOutsideSnapshot;
    boundaryOwner.preDualOwnerLevel = minTree.getAltitude(visibleExternalParent);
    boundaryOwner.oracleDualOwnerLevel = minTree.getAltitude(visibleExternalParent);
    oracle.pixelOwners.push_back(boundaryOwner);

    noninc::LocalDualNodeSnapshot oracleNode;
    oracleNode.nodeId = oracleNodeId;
    oracleNode.parentId = oracleParentOutsideSnapshot;
    oracleNode.parentInSnapshot = false;
    oracleNode.altitude = 15;
    oracleNode.directLocalProperParts = {0};
    oracleNode.supportLocalPixels = {0};
    oracleNode.supportTouchesModified = true;

    oracle.oracleDualNodes = {oracleNode};

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {15};
    plan.nodes.resize(1);
    plan.nodes[0].level = 15;
    plan.nodes[0].supportPixels = {0};
    plan.nodes[0].directProperParts = {0};
    plan.nodes[0].oracleNodeId = oracleNodeId;
    plan.nodes[0].oracleParentId = oracleParentOutsideSnapshot;
    plan.nodes[0].oracleDualOwners = {oracleNodeId};
    plan.nodes[0].preDualOwners = {sourceOwner, visibleExternalParent};
    plan.nodes[0].parentPlanIndex = -1;

    const auto materialization = noninc::materializeLocalRecompositionNodes(
        minTree,
        oracle,
        plan);

    require(materialization.numCreatedNodes() == 2,
            "snapshot-external oracle parent should create the recomposition node and one hidden parent");
    require(materialization.reconstructionMatchesTarget,
            "hidden external-parent connector must preserve reconstruction");
    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    const auto hiddenParentIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.planIndex == -1 && created.level == 16;
        });
    require(hiddenParentIt != materialization.createdNodes.end(),
            "hidden external parent must still be materialized before canonicalization");
    require(!minTree.isAlive(hiddenParentIt->nodeId),
            "empty hidden external parent must be released by local canonicalization");
    require(minTree.getNodeParent(recomposedNode) == visibleExternalParent,
            "recomposed node must attach directly after the empty hidden parent is released");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "hidden external-parent canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

void testSyntheticMaterializationUsesNonOracleHiddenExternalParentConnector() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const NodeId sourceOwner = minTree.getSmallestComponent(0);
    const NodeId visibleExternalParent = minTree.getSmallestComponent(1);
    require(sourceOwner != InvalidNode && visibleExternalParent != InvalidNode,
            "non-oracle hidden external-parent setup must have source and boundary owners");
    require(minTree.getAltitude(sourceOwner) == 10 &&
            minTree.getAltitude(visibleExternalParent) == 20,
            "non-oracle hidden external-parent setup expects levels 10 and 20");

    noninc::LocalDualNodeSnapshot preSource;
    preSource.nodeId = sourceOwner;
    preSource.parentId = visibleExternalParent;
    preSource.parentInSnapshot = true;
    preSource.altitude = 10;
    preSource.directLocalProperParts = {0};
    preSource.supportLocalPixels = {0};
    preSource.supportTouchesModified = true;

    noninc::LocalDualNodeSnapshot preExternalParent;
    preExternalParent.nodeId = visibleExternalParent;
    preExternalParent.parentId = InvalidNode;
    preExternalParent.parentInSnapshot = false;
    preExternalParent.altitude = 20;
    preExternalParent.directLocalProperParts = {1};
    preExternalParent.supportLocalPixels = {0, 1};
    preExternalParent.childrenInSnapshot = {sourceOwner};
    preExternalParent.supportTouchesModified = true;
    preExternalParent.supportTouchesBoundary = true;

    noninc::LocalPartitionInput input;
    input.updatedImage = clone_with_values(minTree.reconstructionImage(), {15, 20, 30, 40});
    input.contraction.dualIsMaxtree = false;
    input.localPixels = {0, 1};
    input.preDualNodes = {preSource, preExternalParent};
    input.pixelMetadata = {
        {0, noninc::LocalMergePixelRole::Modified, sourceOwner},
        {1, noninc::LocalMergePixelRole::Boundary, visibleExternalParent},
    };

    noninc::LocalRecompositionPlan plan;
    plan.affectedLevels = {15};
    plan.nodes.resize(1);
    plan.nodes[0].level = 15;
    plan.nodes[0].supportPixels = {0};
    plan.nodes[0].directProperParts = {0};
    plan.nodes[0].preDualOwners = {sourceOwner, visibleExternalParent};
    plan.nodes[0].parentPlanIndex = -1;

    const auto strippedResult = noninc::materializeStrippedNonOracleLocalRecompositionPlan(
        minTree,
        input,
        plan);
    require(strippedResult.isOracleFree(),
            "non-oracle hidden external-parent helper result must remain oracle-free");
    const auto &materialization = strippedResult.materialization;

    require(materialization.numCreatedNodes() == 2,
            "non-oracle external context parent should create one hidden parent connector");
    require(materialization.reconstructionMatchesTarget,
            "non-oracle hidden external-parent connector must preserve reconstruction");
    const NodeId recomposedNode = materialization.createdNodes.front().nodeId;
    const auto hiddenParentIt = std::find_if(
        materialization.createdNodes.begin(),
        materialization.createdNodes.end(),
        [](const auto &created) {
            return created.planIndex == -1 && created.level == 16;
        });
    require(hiddenParentIt != materialization.createdNodes.end(),
            "non-oracle hidden external parent must still be materialized before canonicalization");
    require(!minTree.isAlive(hiddenParentIt->nodeId),
            "empty non-oracle hidden external parent must be released by local canonicalization");
    require(minTree.getNodeParent(recomposedNode) == visibleExternalParent,
            "non-oracle recomposed node must attach directly after the empty hidden parent is released");
    require(minTree.getSmallestComponent(0) == recomposedNode,
            "modified pixel must move to the non-oracle recomposition node");
    require(minTree.getSmallestComponent(1) == visibleExternalParent,
            "context pixel must remain owned by the external parent");
    require(materialization.emptyCreatedNodesAfterCanonicalization.empty(),
            "non-oracle hidden external-parent canonicalization must not leave created empty nodes");
    requireTreeConsistency(minTree);
}

} // namespace

int main() {
    try {
        testAffectedLevelsFollowDualInsertionIntervals();
        testDualThresholdPredicateUsesTreePolarity();
        testReducedFixtureBuildsLocalThresholdPartitions();
        testNonOracleLocalPartitionMatchesReducedFixture();
        testNonOracleLocalRecompositionMatchesReducedFixture();
        testReducedFixtureBuildsEffectiveRecompositionPlan();
        testExternalConnectivityMergesLocalThresholdComponentsByOracleOwner();
        testNonOracleLocalPartitionMatchesExternalConnectivityFixture();
        testNonOracleLocalRecompositionMatchesExternalConnectivityFixture();
        testNonOracleLocalPartitionPreservesProperThresholdChildrenUnderAncestor();
        testOracleSnapshotCoverageMergesDisjointLocalThresholdComponents();
        testNonOracleLocalPartitionMergesSnapshotCoverageByPreSnapshotSupport();
        testNonOracleLocalRecompositionMatchesSnapshotCoverageFixture();
        testNonOracleLocalRecompositionRecordsMovedStableChildren();
        testStrippedNonOraclePipelineMaterializesStableChildrenEndToEnd();
        testNonOracleLocalRecompositionRecordsMaximalStableChild();
        testLocalOracleDependencyAuditSeparatesCurrentOracleUses();
        testLocalMaterializationDependencyAuditSeparatesNonOracleInputs();
        testSyntheticRecompositionParentIndicesRespectPolarity();
        testReducedFixtureMaterializesEffectiveRecompositionNode();
        testReducedFixtureMaterializesStrippedEffectiveRecompositionNode();
        testReducedFixtureAppliesStrippedNonOracleLocalMergeUpdate();
        testSyntheticMaterializationCreatesLocalParentChildChain();
        testSyntheticMaterializationKeepsEmptyOwnerWithResidualChildUnreleased();
        testSyntheticMaterializationReleasesPromotableEmptyOwnerWithResidualChildWhenAuthorized();
        testSyntheticMaterializationUsesMappedNonRootExternalAttachment();
        testSyntheticMaterializationCreatesLocalChainUnderMappedBoundary();
        testSyntheticMaterializationMovesMultipleStableChildren();
        testSyntheticMaterializationMovesStableChildrenFromNonOraclePlan();
        testSyntheticMaterializationCreatesNonOracleBridgeBelowMovedStableChild();
        testSyntheticMaterializationRepairsNonOracleVisiblePreParentLink();
        testSyntheticMaterializationCreatesNonOracleMissingParentBridgeAboveLiveChild();
        testSyntheticMaterializationCreatesNonOracleMissingChildBridgeBelowLiveParent();
        testSyntheticMaterializationCreatesBridgeBelowMovedStableChild();
        testSyntheticMaterializationUsesHiddenConnectorForDetachedSupportBranch();
        testSyntheticMaterializationUsesNonOracleHiddenSupportConnectorForDetachedBranch();
        testSyntheticMaterializationUsesHiddenExternalParentForSnapshotExternalOracleParent();
        testSyntheticMaterializationUsesNonOracleHiddenExternalParentConnector();
    } catch (const std::exception &e) {
        std::cerr << "details_non_increasing_local_partition_unit_tests: FAIL\n" << e.what() << "\n";
        return 1;
    }

    std::cout << "details_non_increasing_local_partition_unit_tests: OK\n";
    return 0;
}
