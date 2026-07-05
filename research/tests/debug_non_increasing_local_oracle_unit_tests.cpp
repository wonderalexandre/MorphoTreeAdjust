#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/NonIncreasingLocalOracle.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingCasfKernel.hpp"

namespace noninc = debug::non_increasing;
namespace noninc_kernel = morphotree::details::non_increasing;

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

bool anyNodeTouchesModified(const std::vector<noninc::LocalDualNodeSnapshot> &nodes) {
    return std::any_of(nodes.begin(), nodes.end(), [](const auto &node) {
        return node.supportTouchesModified;
    });
}

bool anyNodeTouchesPreserved(const std::vector<noninc::LocalDualNodeSnapshot> &nodes) {
    return std::any_of(nodes.begin(), nodes.end(), [](const auto &node) {
        return node.supportTouchesPreserved;
    });
}

bool anyNodeTouchesBoundary(const std::vector<noninc::LocalDualNodeSnapshot> &nodes) {
    return std::any_of(nodes.begin(), nodes.end(), [](const auto &node) {
        return node.supportTouchesBoundary;
    });
}

void requirePixelOwnerRolesMatchTargetImage(const noninc::LocalMergeOracle &oracle) {
    int modifiedCount = 0;
    int preservedCount = 0;
    int boundaryCount = 0;

    for (const auto &owner : oracle.pixelOwners) {
        require(owner.preDualOwner != InvalidNode, "pre-oracle owner must be valid");
        require(owner.oracleDualOwner != InvalidNode, "oracle owner must be valid");

        const auto expectedLevel = static_cast<int>((*oracle.updatedImage)[owner.pixelId]);
        require(owner.oracleDualOwnerLevel == expectedLevel,
                "oracle owner level must reconstruct the target image value");

        if (owner.role == noninc::LocalMergePixelRole::Modified) {
            ++modifiedCount;
            require(owner.oracleDualOwnerLevel == oracle.contraction.b,
                    "modified pixels must be assigned to the parent level b");
        } else if (owner.role == noninc::LocalMergePixelRole::Preserved) {
            ++preservedCount;
            require((*oracle.updatedImage)[owner.pixelId] == (*oracle.preImage)[owner.pixelId],
                    "preserved pixels must keep their original gray level");
        } else {
            ++boundaryCount;
            require((*oracle.updatedImage)[owner.pixelId] == (*oracle.preImage)[owner.pixelId],
                    "boundary pixels must keep their original gray level");
        }
    }

    require(modifiedCount == static_cast<int>(oracle.contraction.modifiedProperPart.size()),
            "pixel owner snapshot must include every modified pixel exactly once");
    require(preservedCount == static_cast<int>(oracle.preservedPixels.size()),
            "pixel owner snapshot must include every preserved pixel exactly once");
    require(boundaryCount == static_cast<int>(oracle.boundaryPixels.size()),
            "pixel owner snapshot must include every boundary pixel exactly once");
}

void testReducedFixtureLocalOracleCapturesSplitCase() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    noninc_kernel::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc_kernel::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "fixture should not fail before the diagnostic threshold");
    }

    DynamicComponentTree &maxTree = runner.maxTree();
    DynamicComponentTree &minTree = runner.minTree();
    const NodeId nodeId = findReducedFixtureMergeNode(maxTree);
    require(nodeId != InvalidNode, "reduced fixture merge node must be present before threshold 10");

    const auto oracle = noninc::buildLocalMergeOracle(maxTree, minTree, adjacency, nodeId);

    require(oracle.preDualMatchesPrimal, "pre dual tree must reconstruct the same image as the primal tree");
    require(oracle.oracleDualMatchesUpdatedImage, "oracle dual rebuild must reconstruct the target image");
    require(oracle.contraction.primalIsMaxtree, "fixture failure is in the max-tree phase");
    require(!oracle.contraction.dualIsMaxtree, "fixture oracle must rebuild the min-tree");
    require(oracle.contraction.a == 92, "descriptor must preserve level a");
    require(oracle.contraction.b == 80, "descriptor must preserve parent level b");
    require(oracle.contraction.direction == ::debug::region::RegionUpdateDirection::Lower,
            "fixture merge must be a lower update");
    require(oracle.contraction.modifiedProperPart.size() == 1,
            "fixture modified region must contain one direct proper part");
    require(oracle.contraction.modifiedComponents.size() == 1,
            "fixture modified region must have one connected component");
    require(oracle.contraction.preservedChildren.size() == 1,
            "fixture merge must preserve one child");
    require(!oracle.preservedPixels.empty(), "preserved child support must be captured");
    require(!oracle.boundaryPixels.empty(), "local oracle must capture an external boundary");
    require(!oracle.preDualNodes.empty(), "pre dual local snapshot must not be empty");
    require(!oracle.oracleDualNodes.empty(), "oracle dual local snapshot must not be empty");

    requirePixelOwnerRolesMatchTargetImage(oracle);
    require(anyNodeTouchesModified(oracle.oracleDualNodes),
            "oracle snapshot must include a node touching the modified region");
    require(anyNodeTouchesPreserved(oracle.oracleDualNodes),
            "oracle snapshot must include a node touching the preserved child context");
    require(anyNodeTouchesBoundary(oracle.oracleDualNodes),
            "oracle snapshot must include a node touching the external boundary");

    DynamicComponentTree primalCopy = maxTree;
    primalCopy.mergeNodeIntoParent(nodeId);
    require(primalCopy.reconstructionImage()->isEqual(oracle.updatedImage),
            "descriptor target image must match the real primal merge");

    const auto serialized = noninc::debugString(oracle);
    require(serialized.find("LocalMergeOracle") != std::string::npos,
            "local oracle must provide a deterministic debug serialization");
    require(serialized.find("a=92") != std::string::npos &&
                serialized.find("b=80") != std::string::npos,
            "local oracle serialization must preserve the failing level transition");
    require(serialized.find("role=modified") != std::string::npos &&
                serialized.find("role=preserved") != std::string::npos &&
                serialized.find("role=boundary") != std::string::npos,
            "local oracle serialization must expose all local pixel roles");
}

void testLocalOracleRejectsInvalidNode() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    try {
        (void) noninc::buildLocalMergeOracle(maxTree, minTree, adjacency, maxTree.getRoot());
    } catch (const std::runtime_error &) {
        return;
    }
    throw std::runtime_error("local oracle must reject root merges");
}

void testRunnerFirstMismatchOracleCollectionIsIdleOnCleanRun() {
    auto image = make_reduced_split_fixture_image();
    noninc_kernel::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc_kernel::NonIncreasingDualMaintenance::MaterializeThenExplicit);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setFirstInvariantMismatchOracleCollectionEnabled(true);

    for (float threshold : {0.0f, 1.0f, 5.0f, 10.0f}) {
        runner.applyThreshold(threshold);
    }

    require(!runner.hasFirstInvariantMismatch(),
            "clean diagnostic run must preserve the primal/dual reconstruction invariant");
    require(!runner.hasFirstInvariantMismatchOracle(),
            "first-mismatch oracle collection must not retain an oracle when no mismatch occurred");
}

void testRunnerCapturesTargetedLocalOracle() {
    auto image = make_reduced_split_fixture_image();
    noninc_kernel::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc_kernel::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
    }

    const NodeId nodeId = findReducedFixtureMergeNode(runner.maxTree());
    require(nodeId != InvalidNode, "reduced fixture merge node must be present before targeted capture");

    runner.setLocalOracleCaptureTarget(
        10.0f,
        noninc_kernel::NonIncreasingPhase::MaxToMin,
        nodeId);
    runner.applyThreshold(10.0f);

    require(runner.hasCapturedLocalOracle(),
            "runner should capture a targeted local oracle without per-merge mismatch diagnostics");
    const auto &oracle = runner.capturedLocalOracle();
    require(oracle.contraction.selectedNode == nodeId,
            "targeted oracle must refer to the requested node");
    require(oracle.contraction.a == 92 && oracle.contraction.b == 80,
            "targeted oracle must preserve the requested fixture transition");

    runner.clearLocalOracleCaptureTarget();
    require(!runner.hasCapturedLocalOracle(),
            "clearing the local oracle target should clear the captured state");
}

} // namespace

int main() {
    try {
        testReducedFixtureLocalOracleCapturesSplitCase();
        testLocalOracleRejectsInvalidNode();
        testRunnerFirstMismatchOracleCollectionIsIdleOnCleanRun();
        testRunnerCapturesTargetedLocalOracle();
    } catch (const std::exception &e) {
        std::cerr << "debug_non_increasing_local_oracle_unit_tests: FAIL\n" << e.what() << "\n";
        return 1;
    }

    std::cout << "debug_non_increasing_local_oracle_unit_tests: OK\n";
    return 0;
}
