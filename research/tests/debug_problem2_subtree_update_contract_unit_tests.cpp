#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/AttributeComputer.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DualMinMaxTreeIncrementalFilter.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Problem2.hpp"

namespace {

using namespace debug::problem2;

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
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

std::vector<NodeId> getNodesToPrune(DynamicComponentTree &tree,
                                    const std::vector<float> &attribute,
                                    int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree.getRoot());
    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (attribute[static_cast<size_t>(nodeId)] > threshold) {
            for (NodeId childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else if (!tree.isRoot(nodeId)) {
            out.push_back(nodeId);
        }
    }
    return out;
}

std::vector<PixelId> collectSubtreeFootprint(const DynamicComponentTree &tree,
                                             NodeId subtreeRoot) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(pixelId);
        }
    }
    return pixels;
}

struct AcceptedProblem2Fixture {
    NodeId subtreeRoot = InvalidNode;
    Problem2SubtreeUpdateContract<AltitudeType> contract;
};

AcceptedProblem2Fixture findAcceptedFixture(DynamicComponentTree &primalTree,
                                            DynamicComponentTree &dualTree,
                                            const std::vector<float> &primalArea,
                                            int threshold) {
    const auto nodesToPrune = getNodesToPrune(primalTree, primalArea, threshold);
    for (NodeId subtreeRoot : nodesToPrune) {
        auto contract =
            analyzeProblem2SubtreeUpdateContract<AltitudeType>(
                &primalTree,
                &dualTree,
                subtreeRoot);
        if (contract.accepted()) {
            return AcceptedProblem2Fixture{subtreeRoot, std::move(contract)};
        }
    }
    return {};
}

void setAreaComputers(DualMinMaxTreeIncrementalFilter<AltitudeType> &adjust,
                      DynamicAreaComputer &minAreaComputer,
                      DynamicAreaComputer &maxAreaComputer,
                      std::vector<float> &minArea,
                      std::vector<float> &maxArea) {
    minArea = minAreaComputer.compute();
    maxArea = maxAreaComputer.compute();
    adjust.setAttributeComputer(
        minAreaComputer,
        maxAreaComputer,
        std::span<float>(minArea),
        std::span<float>(maxArea));
}

NodeId releaseFirstNonRootDualOwnerForSubtree(DynamicComponentTree &primalTree,
                                              DynamicComponentTree &dualTree,
                                              NodeId subtreeRoot) {
    for (NodeId subtreeNodeId : primalTree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : primalTree.getProperParts(subtreeNodeId)) {
            const NodeId owner = dualTree.getSmallestComponent(pixelId);
            if (owner == InvalidNode ||
                dualTree.isRoot(owner) ||
                dualTree.getNodeParent(owner) == InvalidNode ||
                dualTree.getNodeParent(owner) == owner) {
                continue;
            }

            const NodeId parent = dualTree.getNodeParent(owner);
            dualTree.removeChild(parent, owner, true);
            return owner;
        }
    }
    return InvalidNode;
}

void requireAcceptedContractForCasfFixture(bool updateDualMaxTree) {
    constexpr int kThreshold = 4;
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree &primalTree = updateDualMaxTree ? minTree : maxTree;
    DynamicComponentTree &dualTree = updateDualMaxTree ? maxTree : minTree;

    DynamicAreaComputer primalAreaComputer(&primalTree);
    std::vector<float> primalArea = primalAreaComputer.compute();

    const auto fixture =
        findAcceptedFixture(primalTree, dualTree, primalArea, kThreshold);
    const std::string direction = updateDualMaxTree ? "min-to-max" : "max-to-min";
    require(fixture.subtreeRoot != InvalidNode,
            direction + ": expected an accepted Problem 2 subtree contract");

    const auto &contract = fixture.contract;
    require(contract.accepted(),
            direction + ": contract must be accepted");
    require(contract.readyForRunUpdateTree(),
            direction + ": accepted contract must be ready for runUpdateTree");
    require(contract.isMaxtree == dualTree.isMaxtree(),
            direction + ": contract must mirror the dual tree kind");
    require(contract.subtreeParent == primalTree.getNodeParent(fixture.subtreeRoot),
            direction + ": contract must expose the primal parent of subtreeRoot");
    require(contract.b == static_cast<AltitudeType>(primalTree.getAltitude(contract.subtreeParent)),
            direction + ": contract b must be the primal parent altitude");
    require(contract.nodeCa != InvalidNode && dualTree.isAlive(contract.nodeCa),
            direction + ": contract nodeCa must be live in the dual tree");

    const auto footprint = collectSubtreeFootprint(primalTree, fixture.subtreeRoot);
    require(contract.pixelsC == footprint,
            direction + ": contract C must be exactly the subtree footprint");
    require(contract.regionSize == static_cast<int>(footprint.size()),
            direction + ": contract regionSize must match C");

    if (dualTree.isMaxtree()) {
        require(contract.b > contract.maxC,
                direction + ": max-tree dual contract must satisfy b > max I(C)");
        require(contract.a == contract.minC,
                direction + ": max-tree dual nodeCa level must be min I(C)");
    } else {
        require(contract.b < contract.minC,
                direction + ": min-tree dual contract must satisfy b < min I(C)");
        require(contract.a == contract.maxC,
                direction + ": min-tree dual nodeCa level must be max I(C)");
    }
}

void requireContractMatchesExplicitSweepResult(bool updateDualMaxTree) {
    constexpr int kThreshold = 4;
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DynamicComponentTree &primalTree = updateDualMaxTree ? minTree : maxTree;
    DynamicComponentTree &dualTree = updateDualMaxTree ? maxTree : minTree;

    DynamicAreaComputer maxAreaComputer(&maxTree);
    DynamicAreaComputer minAreaComputer(&minTree);
    std::vector<float> maxArea = maxAreaComputer.compute();
    std::vector<float> minArea = minAreaComputer.compute();
    std::vector<float> &primalArea = updateDualMaxTree ? minArea : maxArea;

    const auto fixture =
        findAcceptedFixture(primalTree, dualTree, primalArea, kThreshold);
    const std::string direction = updateDualMaxTree ? "min-to-max" : "max-to-min";
    require(fixture.subtreeRoot != InvalidNode,
            direction + ": expected an accepted Problem 2 subtree contract");

    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adjacency);
    adjust.setAttributeComputer(
        minAreaComputer,
        maxAreaComputer,
        std::span<float>(minArea),
        std::span<float>(maxArea));

    std::vector<NodeId> nodesToPrune = {fixture.subtreeRoot};
    const int appliedCount = adjust.pruneAndUpdate(&primalTree, &dualTree, nodesToPrune);

    require(appliedCount == 1,
            direction + ": explicit-sweep path must accept the fixture");

    DynamicComponentTree diagnosticMaxTree(image->clone(), true, adjacency);
    DynamicComponentTree diagnosticMinTree(image->clone(), false, adjacency);
    DynamicComponentTree &diagnosticPrimalTree =
        updateDualMaxTree ? diagnosticMinTree : diagnosticMaxTree;
    DynamicComponentTree &diagnosticDualTree =
        updateDualMaxTree ? diagnosticMaxTree : diagnosticMinTree;
    DynamicAreaComputer diagnosticMaxAreaComputer(&diagnosticMaxTree);
    DynamicAreaComputer diagnosticMinAreaComputer(&diagnosticMinTree);
    std::vector<float> diagnosticMaxArea = diagnosticMaxAreaComputer.compute();
    std::vector<float> diagnosticMinArea = diagnosticMinAreaComputer.compute();
    std::vector<float> &diagnosticPrimalArea =
        updateDualMaxTree ? diagnosticMinArea : diagnosticMaxArea;
    const auto diagnosticFixture =
        findAcceptedFixture(
            diagnosticPrimalTree,
            diagnosticDualTree,
            diagnosticPrimalArea,
            kThreshold);
    require(diagnosticFixture.subtreeRoot != InvalidNode,
            direction + ": expected an accepted diagnostic Problem 2 subtree contract");

    DualMinMaxTreeIncrementalFilter<AltitudeType> diagnosticAdjust(
        &diagnosticMinTree,
        &diagnosticMaxTree,
        *adjacency);
    diagnosticAdjust.setAttributeComputer(
        diagnosticMinAreaComputer,
        diagnosticMaxAreaComputer,
        std::span<float>(diagnosticMinArea),
        std::span<float>(diagnosticMaxArea));

    std::vector<NodeId> diagnosticNodesToPrune = {diagnosticFixture.subtreeRoot};
    const auto diagnosticResult =
        debug::problem2::pruneAndUpdateWithDiagnostics(
            diagnosticAdjust,
            &diagnosticPrimalTree,
            &diagnosticDualTree,
            diagnosticNodesToPrune,
            true);
    require(diagnosticResult.inputValid && diagnosticResult.fastGateReady,
            direction + ": diagnostic explicit-sweep path must accept the fixture");
    require(!diagnosticResult.trustedFastGateUsed,
            direction + ": diagnostic explicit-sweep path must use explicit gate bounds");
    require(diagnosticResult.fullGateBoundsComputed,
            direction + ": diagnostic explicit-sweep path must compute full gate bounds");
    require(diagnosticResult.minC == diagnosticFixture.contract.minC,
            direction + ": diagnostic explicit-sweep minC must match the contract");
    require(diagnosticResult.maxC == diagnosticFixture.contract.maxC,
            direction + ": diagnostic explicit-sweep maxC must match the contract");

    DynamicComponentTree trustedMaxTree(image->clone(), true, adjacency);
    DynamicComponentTree trustedMinTree(image->clone(), false, adjacency);
    DynamicComponentTree &trustedPrimalTree =
        updateDualMaxTree ? trustedMinTree : trustedMaxTree;
    DynamicComponentTree &trustedDualTree =
        updateDualMaxTree ? trustedMaxTree : trustedMinTree;
    DynamicAreaComputer trustedMaxAreaComputer(&trustedMaxTree);
    DynamicAreaComputer trustedMinAreaComputer(&trustedMinTree);
    std::vector<float> trustedMaxArea = trustedMaxAreaComputer.compute();
    std::vector<float> trustedMinArea = trustedMinAreaComputer.compute();
    std::vector<float> &trustedPrimalArea =
        updateDualMaxTree ? trustedMinArea : trustedMaxArea;
    const auto trustedFixture =
        findAcceptedFixture(
            trustedPrimalTree,
            trustedDualTree,
            trustedPrimalArea,
            kThreshold);
    require(trustedFixture.subtreeRoot != InvalidNode,
            direction + ": expected an accepted trusted-gate Problem 2 subtree contract");

    DualMinMaxTreeIncrementalFilter<AltitudeType> trustedAdjust(
        &trustedMinTree,
        &trustedMaxTree,
        *adjacency);
    trustedAdjust.setAttributeComputer(
        trustedMinAreaComputer,
        trustedMaxAreaComputer,
        std::span<float>(trustedMinArea),
        std::span<float>(trustedMaxArea));

    std::vector<NodeId> trustedNodesToPrune = {trustedFixture.subtreeRoot};
    const auto trustedResult =
        debug::problem2::pruneAndUpdateWithDiagnostics(
            trustedAdjust,
            &trustedPrimalTree,
            &trustedDualTree,
            trustedNodesToPrune,
            false);
    require(trustedResult.inputValid && trustedResult.fastGateReady,
            direction + ": trusted diagnostic explicit-sweep path must accept the fixture");
    require(trustedResult.regionSize == trustedFixture.contract.regionSize,
            direction + ": trusted diagnostic region size must match the contract");
    require(trustedResult.b == trustedFixture.contract.b,
            direction + ": trusted diagnostic b must match the contract");
    require(trustedResult.trustedFastGateUsed,
            direction + ": trusted diagnostic path must use the trusted fast gate");
    require(!trustedResult.fullGateBoundsComputed,
            direction + ": trusted diagnostic path must not recompute full gate bounds");
    if (trustedDualTree.isMaxtree()) {
        require(trustedResult.minC == trustedFixture.contract.minC,
                direction + ": trusted diagnostic min representative must match the contract");
    } else {
        require(trustedResult.maxC == trustedFixture.contract.maxC,
                direction + ": trusted diagnostic max representative must match the contract");
    }
}

void testProblem2SubtreeContractAcceptsCasfFixtures() {
    requireAcceptedContractForCasfFixture(false);
    requireAcceptedContractForCasfFixture(true);
}

void testProblem2SubtreeContractMatchesExplicitSweep() {
    requireContractMatchesExplicitSweepResult(false);
    requireContractMatchesExplicitSweepResult(true);
}

void testProblem2SubtreeContractRejectsInvalidInputs() {
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const auto nullPrimal =
        analyzeProblem2SubtreeUpdateContract<AltitudeType>(
            nullptr,
            &maxTree,
            maxTree.getRoot());
    require(nullPrimal.status == Problem2SubtreeUpdateStatus::NullPrimalTree,
            "null primal tree must be reported as a formal rejection");

    const auto sameKind =
        analyzeProblem2SubtreeUpdateContract<AltitudeType>(
            &maxTree,
            &maxTree,
            maxTree.getRoot());
    require(sameKind.status == Problem2SubtreeUpdateStatus::MismatchedTreeKinds,
            "same-kind primal/dual trees must be rejected");

    auto smallerImage = ImageUInt8::create(4, 4);
    for (int i = 0; i < smallerImage->getSize(); ++i) {
        (*smallerImage)[i] = static_cast<uint8_t>((i % 3) + 1);
    }
    auto smallerAdjacency = std::make_shared<AdjacencyRelation>(
        smallerImage->getNumRows(),
        smallerImage->getNumCols(),
        1.0);
    DynamicComponentTree smallerMinTree(smallerImage->clone(), false, smallerAdjacency);
    const auto mismatchedDomain =
        analyzeProblem2SubtreeUpdateContract<AltitudeType>(
            &maxTree,
            &smallerMinTree,
            maxTree.getRoot());
    require(mismatchedDomain.status == Problem2SubtreeUpdateStatus::MismatchedTreeDomains,
            "opposite-kind trees over different image domains must be rejected");

    const auto root =
        analyzeProblem2SubtreeUpdateContract<AltitudeType>(
            &minTree,
            &maxTree,
            minTree.getRoot());
    require(root.status == Problem2SubtreeUpdateStatus::InvalidSubtreeRoot,
            "root subtree must be rejected before mutation");
    require(!root.accepted(),
            "invalid root rejection must not be accepted");
}

void testProblem2SubtreeBatchContractAcceptsAndSkipsInvalidRoots() {
    constexpr int kThreshold = 4;
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    DynamicAreaComputer maxAreaComputer(&maxTree);
    std::vector<float> maxArea = maxAreaComputer.compute();
    const auto fixture = findAcceptedFixture(maxTree, minTree, maxArea, kThreshold);
    require(fixture.subtreeRoot != InvalidNode,
            "expected an accepted max-to-min fixture for batch contract");

    const std::vector<NodeId> roots = {
        InvalidNode,
        fixture.subtreeRoot,
    };
    const auto batch =
        analyzeProblem2SubtreeBatchUpdateContract<AltitudeType>(
            &maxTree,
            &minTree,
            roots);

    require(batch.accepted(),
            "batch contract with one skipped root and one accepted event must be accepted");
    require(batch.inputRootCount == 2,
            "batch contract must report all input roots");
    require(batch.skippedNodeCount == 1 && batch.skippedInputIndices.size() == 1,
            "batch contract must skip the invalid root");
    require(batch.skippedInputIndices[0] == 0,
            "batch contract must preserve skipped input index");
    require(batch.attemptedEventCount == 1 && batch.acceptedEventCount == 1,
            "batch contract must attempt and accept exactly one event");
    require(batch.acceptedEvents.size() == 1,
            "batch contract must expose the accepted event descriptor");
    require(batch.acceptedEvents[0].subtreeRoot == fixture.subtreeRoot,
            "batch contract accepted event must match the requested subtree root");
    require(batch.regionSize == fixture.contract.regionSize,
            "batch region size must aggregate accepted event region sizes");
}

void testProblem2SubtreeBatchContractReportsNoValidRoots() {
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const std::vector<NodeId> roots = {
        InvalidNode,
        maxTree.getRoot(),
    };
    const auto batch =
        analyzeProblem2SubtreeBatchUpdateContract<AltitudeType>(
            &maxTree,
            &minTree,
            roots);

    require(batch.status == Problem2SubtreeBatchUpdateStatus::NoValidSubtreeRoots,
            "batch contract with only invalid roots must report no valid roots");
    require(batch.skippedNodeCount == 2,
            "batch contract must count every invalid skipped root");
    require(batch.attemptedEventCount == 0 && batch.acceptedEventCount == 0,
            "batch contract must not attempt events for invalid roots");
}

void testProblem2SubtreeBatchContractStopsAtRejectedEvent() {
    constexpr int kThreshold = 4;
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    DynamicAreaComputer maxAreaComputer(&maxTree);
    std::vector<float> maxArea = maxAreaComputer.compute();
    const auto accepted = findAcceptedFixture(maxTree, minTree, maxArea, kThreshold);
    require(accepted.subtreeRoot != InvalidNode,
            "expected an accepted fixture after the rejected batch event");

    NodeId ownerToRelease = InvalidNode;
    for (PixelId pixelId : accepted.contract.pixelsC) {
        const NodeId owner = minTree.getSmallestComponent(pixelId);
        if (owner != InvalidNode &&
            !minTree.isRoot(owner) &&
            minTree.getNodeParent(owner) != InvalidNode &&
            minTree.getNodeParent(owner) != owner) {
            ownerToRelease = owner;
            break;
        }
    }
    require(ownerToRelease != InvalidNode,
            "expected a releasable dual owner to force a formal rejected event");
    minTree.removeChild(minTree.getNodeParent(ownerToRelease), ownerToRelease, true);

    const auto rejected =
        analyzeProblem2SubtreeUpdateContract<AltitudeType>(
            &maxTree,
            &minTree,
            accepted.subtreeRoot);
    require(rejected.status == Problem2SubtreeUpdateStatus::MissingLiveDualOwner,
            "released dual owner must force a formal MissingLiveDualOwner rejection");

    const std::vector<NodeId> roots = {
        InvalidNode,
        accepted.subtreeRoot,
        accepted.subtreeRoot,
    };
    const auto batch =
        analyzeProblem2SubtreeBatchUpdateContract<AltitudeType>(
            &maxTree,
            &minTree,
            roots);

    require(batch.stoppedOnRejectedEvent(),
            "batch contract must stop at the first real rejected event");
    require(batch.skippedNodeCount == 1,
            "batch contract must skip invalid roots before the rejection");
    require(batch.attemptedEventCount == 1 && batch.acceptedEventCount == 0,
            "batch contract must not analyze events after the first rejection");
    require(batch.stoppedInputIndex == 1,
            "batch contract must report the input index of the rejected event");
    require(batch.stoppedSubtreeRoot == accepted.subtreeRoot,
            "batch contract must report the rejected subtree root");
    require(batch.rejectedEvent.status == rejected.status,
            "batch contract must preserve the formal rejection cause");
    require(batch.acceptedEvents.empty(),
            "batch contract must not expose events after the first rejection");
}

void testProblem2SubtreeBatchContractMatchesPruneAndUpdateAcceptedBatch() {
    constexpr int kThreshold = 4;
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree contractMaxTree(image->clone(), true, adjacency);
    DynamicComponentTree contractMinTree(image->clone(), false, adjacency);
    DynamicAreaComputer contractMaxAreaComputer(&contractMaxTree);
    std::vector<float> contractMaxArea = contractMaxAreaComputer.compute();
    const auto fixture =
        findAcceptedFixture(contractMaxTree, contractMinTree, contractMaxArea, kThreshold);
    require(fixture.subtreeRoot != InvalidNode,
            "expected an accepted fixture for accepted batch/updater comparison");

    const std::vector<NodeId> roots = {
        InvalidNode,
        fixture.subtreeRoot,
    };
    const auto batch =
        analyzeProblem2SubtreeBatchUpdateContract<AltitudeType>(
            &contractMaxTree,
            &contractMinTree,
            roots);
    require(batch.accepted(),
            "contract batch fixture must be accepted");

    DynamicComponentTree updateMaxTree(image->clone(), true, adjacency);
    DynamicComponentTree updateMinTree(image->clone(), false, adjacency);
    DynamicAreaComputer updateMaxAreaComputer(&updateMaxTree);
    DynamicAreaComputer updateMinAreaComputer(&updateMinTree);
    std::vector<float> updateMaxArea;
    std::vector<float> updateMinArea;
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(
        &updateMinTree,
        &updateMaxTree,
        *adjacency);
    setAreaComputers(
        adjust,
        updateMinAreaComputer,
        updateMaxAreaComputer,
        updateMinArea,
        updateMaxArea);

    std::vector<NodeId> mutableRoots = roots;
    const auto result =
        debug::problem2::pruneAndUpdateWithDiagnostics(
            adjust,
            &updateMaxTree,
            &updateMinTree,
            mutableRoots,
            false);

    require(result.inputValid && result.fastGateReady,
            "accepted contract batch must be accepted by pruneAndUpdate");
    require(result.localUpdateApplied && result.primalPruneApplied,
            "accepted contract batch must be applied by pruneAndUpdate");
    require(result.eventCount == batch.attemptedEventCount,
            "pruneAndUpdate event count must match batch contract attempts");
    require(result.skippedNodeCount == batch.skippedNodeCount,
            "pruneAndUpdate skipped count must match batch contract");
    require(result.regionSize == batch.regionSize,
            "pruneAndUpdate region size must match batch contract");
    require(result.prunedNodeCount == batch.acceptedEventCount,
            "pruneAndUpdate pruned count must match accepted batch events");
}

void testProblem2SubtreeBatchContractMatchesPruneAndUpdateRejectedBatch() {
    constexpr int kThreshold = 4;
    auto image = make_large_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    DynamicComponentTree contractMaxTree(image->clone(), true, adjacency);
    DynamicComponentTree contractMinTree(image->clone(), false, adjacency);
    DynamicAreaComputer contractMaxAreaComputer(&contractMaxTree);
    std::vector<float> contractMaxArea = contractMaxAreaComputer.compute();
    const auto fixture =
        findAcceptedFixture(contractMaxTree, contractMinTree, contractMaxArea, kThreshold);
    require(fixture.subtreeRoot != InvalidNode,
            "expected an accepted fixture before forcing rejected batch");
    require(releaseFirstNonRootDualOwnerForSubtree(
                contractMaxTree,
                contractMinTree,
                fixture.subtreeRoot) != InvalidNode,
            "expected a releasable owner in contract dual tree");

    const std::vector<NodeId> roots = {
        InvalidNode,
        fixture.subtreeRoot,
        fixture.subtreeRoot,
    };
    const auto batch =
        analyzeProblem2SubtreeBatchUpdateContract<AltitudeType>(
            &contractMaxTree,
            &contractMinTree,
            roots);
    require(batch.stoppedOnRejectedEvent(),
            "contract batch fixture must stop at a rejected event");
    require(batch.rejectedEvent.status == Problem2SubtreeUpdateStatus::MissingLiveDualOwner,
            "contract batch rejection must be MissingLiveDualOwner");

    DynamicComponentTree updateMaxTree(image->clone(), true, adjacency);
    DynamicComponentTree updateMinTree(image->clone(), false, adjacency);
    require(releaseFirstNonRootDualOwnerForSubtree(
                updateMaxTree,
                updateMinTree,
                fixture.subtreeRoot) != InvalidNode,
            "expected a releasable owner in updater dual tree");
    DynamicAreaComputer updateMaxAreaComputer(&updateMaxTree);
    DynamicAreaComputer updateMinAreaComputer(&updateMinTree);
    std::vector<float> updateMaxArea;
    std::vector<float> updateMinArea;
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(
        &updateMinTree,
        &updateMaxTree,
        *adjacency);
    setAreaComputers(
        adjust,
        updateMinAreaComputer,
        updateMaxAreaComputer,
        updateMinArea,
        updateMaxArea);

    std::vector<NodeId> mutableRoots = roots;
    const auto result =
        debug::problem2::pruneAndUpdateWithDiagnostics(
            adjust,
            &updateMaxTree,
            &updateMinTree,
            mutableRoots,
            false);

    require(!result.inputValid,
            "pruneAndUpdate must reject the missing-owner batch event");
    require(result.eventCount == batch.attemptedEventCount,
            "rejected pruneAndUpdate event count must match batch contract attempts");
    require(result.skippedNodeCount == batch.skippedNodeCount,
            "rejected pruneAndUpdate skipped count must match batch contract");
    require(result.regionSize == batch.regionSize,
            "rejected pruneAndUpdate region size must match batch contract");
    require(result.prunedNodeCount == batch.acceptedEventCount,
            "rejected pruneAndUpdate pruned count must match accepted batch events");
    require(result.message.find("without a live dual owner") != std::string::npos,
            "pruneAndUpdate rejection message must match MissingLiveDualOwner cause");
}

} // namespace

int main() {
    testProblem2SubtreeContractAcceptsCasfFixtures();
    testProblem2SubtreeContractMatchesExplicitSweep();
    testProblem2SubtreeContractRejectsInvalidInputs();
    testProblem2SubtreeBatchContractAcceptsAndSkipsInvalidRoots();
    testProblem2SubtreeBatchContractReportsNoValidRoots();
    testProblem2SubtreeBatchContractStopsAtRejectedEvent();
    testProblem2SubtreeBatchContractMatchesPruneAndUpdateAcceptedBatch();
    testProblem2SubtreeBatchContractMatchesPruneAndUpdateRejectedBatch();
    return 0;
}
