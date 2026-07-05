#pragma once

#include <span>
#include <vector>

#include "Algorithm2FlatRegionMergeUpdater.hpp"
#include "AttributeDirtyRecomputer.hpp"
#include "../../AttributeComputer.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::problem2 {

using ::debug::region::AttributeDirtyRecomputeResult;
using ::debug::region::AttributeDirtyRecomputer;

using Algorithm2AreaAttributeUpdateResult = AttributeDirtyRecomputeResult;

class Algorithm2FlatRaiseAreaAttributeUpdater {
public:
    static Algorithm2AreaAttributeUpdateResult recomputeArea(
        DynamicComponentTree &tree,
        const Algorithm2FlatRaiseApplyResult &structuralResult,
        std::vector<float> &areaBuffer) {
        std::vector<NodeId> dirtySeeds;
        auto appendSeed = [&](NodeId nodeId) {
            if (nodeId != InvalidNode) {
                dirtySeeds.push_back(nodeId);
            }
        };

        for (NodeId nodeId : structuralResult.attributeDirtySeeds) {
            appendSeed(nodeId);
        }
        appendSeed(structuralResult.nodeQb);
        appendSeed(structuralResult.recompositionRoot);
        appendSeed(structuralResult.sourceNodeQa);
        for (NodeId nodeId : structuralResult.createdNodes) {
            appendSeed(nodeId);
        }
        for (NodeId nodeId : structuralResult.xiCircleNodes) {
            appendSeed(nodeId);
        }

        DynamicAreaComputer computer(&tree);
        return AttributeDirtyRecomputer::recomputeAncestorClosure(
            tree,
            computer,
            std::span<const NodeId>(dirtySeeds.data(), dirtySeeds.size()),
            std::span<const NodeId>(structuralResult.releasedNodes.data(),
                                    structuralResult.releasedNodes.size()),
            areaBuffer);
    }
};

} // namespace debug::problem2
