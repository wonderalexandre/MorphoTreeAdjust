#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "../../AttributeComputer.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::region {

struct AttributeDirtyRecomputeResult {
    std::vector<NodeId> dirtySeeds;
    std::vector<NodeId> dirtyClosure;
    std::vector<NodeId> recomputedNodes;

    int numRecomputedNodes() const {
        return static_cast<int>(recomputedNodes.size());
    }
};

class AttributeDirtyRecomputer {
private:
    static void prepareBufferAndReleasedNodes(
        const DynamicComponentTree &tree,
        DynamicAttributeComputer &computer,
        std::span<const NodeId> releasedNodes,
        std::vector<float> &buffer) {
        const auto slots = static_cast<size_t>(tree.getNumInternalNodeSlots());
        if (buffer.size() < slots) {
            buffer.resize(slots, 0.0f);
        }

        for (NodeId nodeId : releasedNodes) {
            if (nodeId >= 0 && static_cast<size_t>(nodeId) < buffer.size()) {
                buffer[static_cast<size_t>(nodeId)] = 0.0f;
            }
            computer.onNodeRemoved(nodeId);
        }
    }

    static void recomputeBottomUp(
        DynamicComponentTree &tree,
        DynamicAttributeComputer &computer,
        const std::vector<NodeId> &nodesTopDown,
        std::vector<float> &buffer,
        std::vector<NodeId> &recomputedNodes,
        bool rebuildLocalFromTree) {
        std::span<float> bufferSpan(buffer.data(), buffer.size());
        for (auto it = nodesTopDown.rbegin();
             it != nodesTopDown.rend();
             ++it) {
            const NodeId nodeId = *it;
            if (rebuildLocalFromTree) {
                computer.preProcessingFromTree(nodeId, bufferSpan);
            } else {
                computer.preProcessing(nodeId, bufferSpan);
            }
            for (NodeId childId : tree.getChildren(nodeId)) {
                computer.mergeProcessing(nodeId, childId, bufferSpan);
            }
            computer.postProcessing(nodeId, bufferSpan);
            recomputedNodes.push_back(nodeId);
        }
    }

public:
    static std::vector<NodeId> buildAncestorClosure(
        const DynamicComponentTree &tree,
        std::span<const NodeId> dirtySeeds) {
        const auto slots = static_cast<size_t>(tree.getNumInternalNodeSlots());
        std::vector<unsigned char> dirty(slots, 0);

        for (NodeId seed : dirtySeeds) {
            if (seed == InvalidNode ||
                seed < 0 ||
                static_cast<size_t>(seed) >= slots ||
                !tree.isAlive(seed)) {
                continue;
            }

            NodeId nodeId = seed;
            while (nodeId != InvalidNode &&
                   nodeId >= 0 &&
                   static_cast<size_t>(nodeId) < slots &&
                   tree.isAlive(nodeId)) {
                dirty[static_cast<size_t>(nodeId)] = 1;
                const NodeId parentId = tree.getNodeParent(nodeId);
                if (parentId == InvalidNode || parentId == nodeId) {
                    break;
                }
                nodeId = parentId;
            }
        }

        std::vector<NodeId> closure;
        const NodeId root = tree.getRoot();
        if (root == InvalidNode || !tree.isAlive(root)) {
            return closure;
        }

        for (NodeId nodeId : tree.getNodeSubtree(root)) {
            if (nodeId >= 0 &&
                static_cast<size_t>(nodeId) < dirty.size() &&
                dirty[static_cast<size_t>(nodeId)] != 0) {
                closure.push_back(nodeId);
            }
        }
        return closure;
    }

    static AttributeDirtyRecomputeResult recomputeAncestorClosure(
        DynamicComponentTree &tree,
        DynamicAttributeComputer &computer,
        std::span<const NodeId> dirtySeeds,
        std::span<const NodeId> releasedNodes,
        std::vector<float> &buffer,
        bool rebuildLocalFromTree = false) {
        AttributeDirtyRecomputeResult result;
        result.dirtySeeds.assign(dirtySeeds.begin(), dirtySeeds.end());

        prepareBufferAndReleasedNodes(tree, computer, releasedNodes, buffer);

        result.dirtyClosure = buildAncestorClosure(tree, dirtySeeds);
        recomputeBottomUp(tree,
                          computer,
                          result.dirtyClosure,
                          buffer,
                          result.recomputedNodes,
                          rebuildLocalFromTree);

        return result;
    }

    static AttributeDirtyRecomputeResult recomputeFullTree(
        DynamicComponentTree &tree,
        DynamicAttributeComputer &computer,
        std::span<const NodeId> releasedNodes,
        std::vector<float> &buffer,
        bool rebuildLocalFromTree = false) {
        AttributeDirtyRecomputeResult result;
        prepareBufferAndReleasedNodes(tree, computer, releasedNodes, buffer);

        const NodeId root = tree.getRoot();
        if (root == InvalidNode || !tree.isAlive(root)) {
            return result;
        }

        result.dirtySeeds.push_back(root);
        for (NodeId nodeId : tree.getNodeSubtree(root)) {
            result.dirtyClosure.push_back(nodeId);
        }
        recomputeBottomUp(tree,
                          computer,
                          result.dirtyClosure,
                          buffer,
                          result.recomputedNodes,
                          rebuildLocalFromTree);

        return result;
    }
};

} // namespace debug::region
