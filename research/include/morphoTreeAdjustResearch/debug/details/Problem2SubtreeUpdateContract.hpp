#pragma once

#include <span>
#include <string>
#include <utility>
#include <vector>

#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem2 {

enum class Problem2SubtreeUpdateStatus {
    Accepted,
    NullPrimalTree,
    NullDualTree,
    MismatchedTreeKinds,
    MismatchedTreeDomains,
    InvalidSubtreeRoot,
    RootSubtreeHasNoParent,
    EmptyImplicitRegion,
    OutOfBoundsPixel,
    MissingLiveDualOwner,
    TargetLevelNotOutsideRegionRange,
};

inline const char *problem2SubtreeUpdateStatusName(Problem2SubtreeUpdateStatus status) {
    switch (status) {
    case Problem2SubtreeUpdateStatus::Accepted:
        return "accepted";
    case Problem2SubtreeUpdateStatus::NullPrimalTree:
        return "null_primal_tree";
    case Problem2SubtreeUpdateStatus::NullDualTree:
        return "null_dual_tree";
    case Problem2SubtreeUpdateStatus::MismatchedTreeKinds:
        return "mismatched_tree_kinds";
    case Problem2SubtreeUpdateStatus::MismatchedTreeDomains:
        return "mismatched_tree_domains";
    case Problem2SubtreeUpdateStatus::InvalidSubtreeRoot:
        return "invalid_subtree_root";
    case Problem2SubtreeUpdateStatus::RootSubtreeHasNoParent:
        return "root_subtree_has_no_parent";
    case Problem2SubtreeUpdateStatus::EmptyImplicitRegion:
        return "empty_implicit_region";
    case Problem2SubtreeUpdateStatus::OutOfBoundsPixel:
        return "out_of_bounds_pixel";
    case Problem2SubtreeUpdateStatus::MissingLiveDualOwner:
        return "missing_live_dual_owner";
    case Problem2SubtreeUpdateStatus::TargetLevelNotOutsideRegionRange:
        return "target_level_not_outside_region_range";
    }
    return "unknown";
}

enum class Problem2SubtreeBatchUpdateStatus {
    Accepted,
    NullPrimalTree,
    NullDualTree,
    MismatchedTreeKinds,
    MismatchedTreeDomains,
    NoValidSubtreeRoots,
    RejectedEvent,
};

inline const char *problem2SubtreeBatchUpdateStatusName(Problem2SubtreeBatchUpdateStatus status) {
    switch (status) {
    case Problem2SubtreeBatchUpdateStatus::Accepted:
        return "accepted";
    case Problem2SubtreeBatchUpdateStatus::NullPrimalTree:
        return "null_primal_tree";
    case Problem2SubtreeBatchUpdateStatus::NullDualTree:
        return "null_dual_tree";
    case Problem2SubtreeBatchUpdateStatus::MismatchedTreeKinds:
        return "mismatched_tree_kinds";
    case Problem2SubtreeBatchUpdateStatus::MismatchedTreeDomains:
        return "mismatched_tree_domains";
    case Problem2SubtreeBatchUpdateStatus::NoValidSubtreeRoots:
        return "no_valid_subtree_roots";
    case Problem2SubtreeBatchUpdateStatus::RejectedEvent:
        return "rejected_event";
    }
    return "unknown";
}

template<typename PixelType = AltitudeType>
struct Problem2SubtreeUpdateContract {
    Problem2SubtreeUpdateStatus status = Problem2SubtreeUpdateStatus::InvalidSubtreeRoot;
    bool isMaxtree = false;
    NodeId subtreeRoot = InvalidNode;
    NodeId subtreeParent = InvalidNode;
    NodeId nodeCa = InvalidNode;
    PixelType a{};
    PixelType minC{};
    PixelType maxC{};
    PixelType b{};
    int subtreeNodeCount = 0;
    int regionSize = 0;
    bool fastGateReady = false;
    std::vector<PixelId> pixelsC;
    std::string message;

    bool accepted() const {
        return status == Problem2SubtreeUpdateStatus::Accepted;
    }

    bool readyForRunUpdateTree() const {
        return accepted() && fastGateReady && nodeCa != InvalidNode && regionSize > 0;
    }
};

template<typename PixelType = AltitudeType>
struct Problem2SubtreeBatchUpdateContract {
    Problem2SubtreeBatchUpdateStatus status = Problem2SubtreeBatchUpdateStatus::NoValidSubtreeRoots;
    int inputRootCount = 0;
    int skippedNodeCount = 0;
    int attemptedEventCount = 0;
    int acceptedEventCount = 0;
    int regionSize = 0;
    int stoppedInputIndex = -1;
    NodeId stoppedSubtreeRoot = InvalidNode;
    std::vector<int> skippedInputIndices;
    std::vector<NodeId> skippedSubtreeRoots;
    std::vector<Problem2SubtreeUpdateContract<PixelType>> acceptedEvents;
    Problem2SubtreeUpdateContract<PixelType> rejectedEvent;
    std::string message;

    bool accepted() const {
        return status == Problem2SubtreeBatchUpdateStatus::Accepted;
    }

    bool stoppedOnRejectedEvent() const {
        return status == Problem2SubtreeBatchUpdateStatus::RejectedEvent;
    }
};

namespace detail {

template<typename PixelType = AltitudeType>
Problem2SubtreeUpdateContract<PixelType> rejectSubtreeUpdateContract(
    Problem2SubtreeUpdateStatus status,
    NodeId subtreeRoot,
    const std::string &message) {
    Problem2SubtreeUpdateContract<PixelType> out;
    out.status = status;
    out.subtreeRoot = subtreeRoot;
    out.message = message;
    return out;
}

template<typename PixelType = AltitudeType>
Problem2SubtreeBatchUpdateContract<PixelType> rejectSubtreeBatchUpdateContract(
    Problem2SubtreeBatchUpdateStatus status,
    const std::string &message) {
    Problem2SubtreeBatchUpdateContract<PixelType> out;
    out.status = status;
    out.message = message;
    return out;
}

} // namespace detail

template<typename PixelType = AltitudeType>
Problem2SubtreeUpdateContract<PixelType> analyzeProblem2SubtreeUpdateContract(
    const DynamicComponentTree *primalTree,
    const DynamicComponentTree *dualTree,
    NodeId subtreeRoot) {
    if (primalTree == nullptr) {
        return detail::rejectSubtreeUpdateContract<PixelType>(
            Problem2SubtreeUpdateStatus::NullPrimalTree,
            subtreeRoot,
            "Problem 2 subtree update requires a valid primal tree");
    }
    if (dualTree == nullptr) {
        return detail::rejectSubtreeUpdateContract<PixelType>(
            Problem2SubtreeUpdateStatus::NullDualTree,
            subtreeRoot,
            "Problem 2 subtree update requires a valid dual tree");
    }

    Problem2SubtreeUpdateContract<PixelType> out;
    out.subtreeRoot = subtreeRoot;
    out.isMaxtree = dualTree->isMaxtree();

    if (primalTree->isMaxtree() == dualTree->isMaxtree()) {
        out.status = Problem2SubtreeUpdateStatus::MismatchedTreeKinds;
        out.message = "Problem 2 subtree update requires opposite primal and dual tree kinds";
        return out;
    }
    if (primalTree->getNumRowsOfImage() != dualTree->getNumRowsOfImage() ||
        primalTree->getNumColsOfImage() != dualTree->getNumColsOfImage() ||
        primalTree->getNumTotalProperParts() != dualTree->getNumTotalProperParts()) {
        out.status = Problem2SubtreeUpdateStatus::MismatchedTreeDomains;
        out.message = "Problem 2 subtree update requires primal and dual trees over the same image domain";
        return out;
    }

    if (subtreeRoot == InvalidNode ||
        subtreeRoot < 0 ||
        subtreeRoot >= primalTree->getNumInternalNodeSlots() ||
        !primalTree->isNode(subtreeRoot) ||
        !primalTree->isAlive(subtreeRoot) ||
        primalTree->isRoot(subtreeRoot)) {
        out.status = Problem2SubtreeUpdateStatus::InvalidSubtreeRoot;
        out.message = "Problem 2 subtree update requires a live non-root primal subtree root";
        return out;
    }

    out.subtreeParent = primalTree->getNodeParent(subtreeRoot);
    if (out.subtreeParent == InvalidNode ||
        out.subtreeParent == subtreeRoot ||
        !primalTree->isAlive(out.subtreeParent)) {
        out.status = Problem2SubtreeUpdateStatus::RootSubtreeHasNoParent;
        out.message = "Problem 2 subtree update cannot locate a live parent for subtreeRoot";
        return out;
    }
    out.b = static_cast<PixelType>(primalTree->getAltitude(out.subtreeParent));

    bool firstOwner = true;
    for (NodeId subtreeNodeId : primalTree->getNodeSubtree(subtreeRoot)) {
        ++out.subtreeNodeCount;
        for (PixelId pixelId : primalTree->getProperParts(subtreeNodeId)) {
            if (pixelId < 0 ||
                pixelId >= primalTree->getNumTotalProperParts() ||
                pixelId >= dualTree->getNumTotalProperParts()) {
                out.status = Problem2SubtreeUpdateStatus::OutOfBoundsPixel;
                out.message = "Problem 2 implicit region C contains an out-of-bounds pixel";
                return out;
            }

            out.pixelsC.push_back(pixelId);

            const NodeId ownerNode = dualTree->getSmallestComponent(pixelId);
            if (ownerNode == InvalidNode || !dualTree->isAlive(ownerNode)) {
                out.status = Problem2SubtreeUpdateStatus::MissingLiveDualOwner;
                out.message = "Problem 2 implicit region C contains a pixel without a live dual owner";
                return out;
            }

            const PixelType level = static_cast<PixelType>(dualTree->getAltitude(ownerNode));
            if (firstOwner) {
                out.minC = level;
                out.maxC = level;
                firstOwner = false;
            } else {
                if (level < out.minC) {
                    out.minC = level;
                }
                if (level > out.maxC) {
                    out.maxC = level;
                }
            }

            if (out.nodeCa == InvalidNode ||
                ((out.isMaxtree && level < out.a) ||
                 (!out.isMaxtree && level > out.a))) {
                out.nodeCa = ownerNode;
                out.a = level;
            }
        }
    }

    out.regionSize = static_cast<int>(out.pixelsC.size());
    if (out.regionSize == 0 || out.nodeCa == InvalidNode || firstOwner) {
        out.status = Problem2SubtreeUpdateStatus::EmptyImplicitRegion;
        out.message = "Problem 2 subtree update received an empty implicit region C";
        return out;
    }

    out.fastGateReady = out.isMaxtree ? (out.b > out.maxC) : (out.b < out.minC);
    if (!out.fastGateReady) {
        out.status = Problem2SubtreeUpdateStatus::TargetLevelNotOutsideRegionRange;
        out.message = out.isMaxtree
            ? "Problem 2 runUpdateTree path requires b > max I(C) for a max-tree dual"
            : "Problem 2 runUpdateTree path requires b < min I(C) for a min-tree dual";
        return out;
    }

    out.status = Problem2SubtreeUpdateStatus::Accepted;
    out.message = "Problem 2 subtree update contract accepted";
    return out;
}

template<typename PixelType = AltitudeType>
Problem2SubtreeBatchUpdateContract<PixelType> analyzeProblem2SubtreeBatchUpdateContract(
    const DynamicComponentTree *primalTree,
    const DynamicComponentTree *dualTree,
    std::span<const NodeId> subtreeRoots) {
    if (primalTree == nullptr) {
        return detail::rejectSubtreeBatchUpdateContract<PixelType>(
            Problem2SubtreeBatchUpdateStatus::NullPrimalTree,
            "Problem 2 subtree batch update requires a valid primal tree");
    }
    if (dualTree == nullptr) {
        return detail::rejectSubtreeBatchUpdateContract<PixelType>(
            Problem2SubtreeBatchUpdateStatus::NullDualTree,
            "Problem 2 subtree batch update requires a valid dual tree");
    }

    Problem2SubtreeBatchUpdateContract<PixelType> out;
    out.inputRootCount = static_cast<int>(subtreeRoots.size());

    if (primalTree->isMaxtree() == dualTree->isMaxtree()) {
        out.status = Problem2SubtreeBatchUpdateStatus::MismatchedTreeKinds;
        out.message = "Problem 2 subtree batch update requires opposite primal and dual tree kinds";
        return out;
    }
    if (primalTree->getNumRowsOfImage() != dualTree->getNumRowsOfImage() ||
        primalTree->getNumColsOfImage() != dualTree->getNumColsOfImage() ||
        primalTree->getNumTotalProperParts() != dualTree->getNumTotalProperParts()) {
        out.status = Problem2SubtreeBatchUpdateStatus::MismatchedTreeDomains;
        out.message = "Problem 2 subtree batch update requires primal and dual trees over the same image domain";
        return out;
    }

    std::vector<unsigned char> logicallyPruned(
        static_cast<std::size_t>(primalTree->getNumInternalNodeSlots()),
        0);

    for (int inputIndex = 0; inputIndex < static_cast<int>(subtreeRoots.size()); ++inputIndex) {
        const NodeId subtreeRoot = subtreeRoots[static_cast<std::size_t>(inputIndex)];
        if (subtreeRoot == InvalidNode ||
            subtreeRoot < 0 ||
            subtreeRoot >= primalTree->getNumInternalNodeSlots() ||
            primalTree->isRoot(subtreeRoot) ||
            !primalTree->isNode(subtreeRoot) ||
            !primalTree->isAlive(subtreeRoot) ||
            logicallyPruned[static_cast<std::size_t>(subtreeRoot)] != 0) {
            ++out.skippedNodeCount;
            out.skippedInputIndices.push_back(inputIndex);
            out.skippedSubtreeRoots.push_back(subtreeRoot);
            continue;
        }

        auto event = analyzeProblem2SubtreeUpdateContract<PixelType>(
            primalTree,
            dualTree,
            subtreeRoot);
        ++out.attemptedEventCount;
        if (!event.accepted()) {
            out.status = Problem2SubtreeBatchUpdateStatus::RejectedEvent;
            out.stoppedInputIndex = inputIndex;
            out.stoppedSubtreeRoot = subtreeRoot;
            out.rejectedEvent = std::move(event);
            out.message = "Problem 2 subtree batch update stopped at the first rejected event";
            return out;
        }

        out.regionSize += event.regionSize;
        ++out.acceptedEventCount;
        for (NodeId subtreeNodeId : primalTree->getNodeSubtree(subtreeRoot)) {
            if (subtreeNodeId >= 0 &&
                subtreeNodeId < static_cast<NodeId>(logicallyPruned.size())) {
                logicallyPruned[static_cast<std::size_t>(subtreeNodeId)] = 1;
            }
        }
        out.acceptedEvents.push_back(std::move(event));
    }

    if (out.acceptedEventCount == 0) {
        out.status = Problem2SubtreeBatchUpdateStatus::NoValidSubtreeRoots;
        out.message = "Problem 2 subtree batch update found no valid subtree roots";
        return out;
    }

    out.status = Problem2SubtreeBatchUpdateStatus::Accepted;
    out.message = "Problem 2 subtree batch update contract accepted";
    return out;
}

template<typename PixelType = AltitudeType>
Problem2SubtreeBatchUpdateContract<PixelType> analyzeProblem2SubtreeBatchUpdateContract(
    const DynamicComponentTree *primalTree,
    const DynamicComponentTree *dualTree,
    const std::vector<NodeId> &subtreeRoots) {
    return analyzeProblem2SubtreeBatchUpdateContract<PixelType>(
        primalTree,
        dualTree,
        std::span<const NodeId>(subtreeRoots.data(), subtreeRoots.size()));
}

} // namespace debug::problem2
