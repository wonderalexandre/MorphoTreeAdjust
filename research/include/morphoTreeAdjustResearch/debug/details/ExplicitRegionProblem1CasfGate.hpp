#pragma once

#include <algorithm>
#include <span>
#include <string>

#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem1 {

template<typename PixelType = AltitudeType>
struct ExplicitRegionProblem1CasfFastGate {
    bool inputValid = false;
    bool acceptedByPublishedProblem = false;
    bool readyUnderCasfAssumptions = false;
    bool isMaxtree = false;
    PixelType minC{};
    PixelType maxC{};
    PixelType b{};
    int regionSize = 0;
    std::string rejectionReason;
};

template<typename PixelType = AltitudeType>
ExplicitRegionProblem1CasfFastGate<PixelType> rejectCasfProblem1FastGate(
    const std::string &reason) {
    ExplicitRegionProblem1CasfFastGate<PixelType> out;
    out.rejectionReason = reason;
    return out;
}

/**
 * @brief Cheap CASF-specific Problem 1 direction gate.
 *
 * This helper is intentionally weaker than the full H1-H8 preflight. It only
 * checks that the explicit region extracted from a valid CASF prune event is
 * moving in the published Problem 1 direction in the current dual tree:
 *
 * - max-tree update: b > max I-(C);
 * - min-tree update: b < min I-(C).
 *
 * The caller must provide the CASF-side assumptions: C comes from a live
 * primal subtree, is non-empty and connected, and the explicit updater is
 * allowed to materialize partial dual owners and create the target level b.
 */
template<typename PixelType = AltitudeType>
ExplicitRegionProblem1CasfFastGate<PixelType> checkCasfSubtreeProblem1FastGate(
    const DynamicComponentTree &dualTree,
    std::span<const PixelId> pixelsC,
    PixelType b) {
    if (pixelsC.empty()) {
        return rejectCasfProblem1FastGate<PixelType>(
            "CASF Problem 1 gate requires a non-empty region C");
    }

    ExplicitRegionProblem1CasfFastGate<PixelType> out;
    out.isMaxtree = dualTree.isMaxtree();
    out.b = b;
    out.regionSize = static_cast<int>(pixelsC.size());

    const int domainSize = dualTree.getNumTotalProperParts();
    bool firstPixel = true;
    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            return rejectCasfProblem1FastGate<PixelType>(
                "CASF Problem 1 gate found an out-of-bounds pixel");
        }

        const NodeId ownerNode = dualTree.getSmallestComponent(pixelId);
        if (ownerNode == InvalidNode || !dualTree.isAlive(ownerNode)) {
            return rejectCasfProblem1FastGate<PixelType>(
                "CASF Problem 1 gate found a pixel without a live dual owner");
        }

        const PixelType level =
            static_cast<PixelType>(dualTree.getAltitude(ownerNode));
        if (firstPixel) {
            out.minC = level;
            out.maxC = level;
            firstPixel = false;
        } else {
            out.minC = std::min(out.minC, level);
            out.maxC = std::max(out.maxC, level);
        }
    }

    out.inputValid = true;
    out.acceptedByPublishedProblem =
        dualTree.isMaxtree() ? (b > out.maxC) : (b < out.minC);
    out.readyUnderCasfAssumptions = out.acceptedByPublishedProblem;
    if (!out.acceptedByPublishedProblem) {
        out.rejectionReason = dualTree.isMaxtree()
            ? "CASF max-tree update requires b > max I-(C)"
            : "CASF min-tree update requires b < min I-(C)";
    }
    return out;
}

} // namespace debug::problem1
