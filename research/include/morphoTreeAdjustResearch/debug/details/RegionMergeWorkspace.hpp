#pragma once

#include <algorithm>
#include <map>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "RegionLevelChange.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::region {

template<typename PixelType = AltitudeType>
struct RegionMergeLevelWorkspace {
    PixelType level{};
    std::vector<NodeId> mergeNodesVMinus;
    NodeId mergedNodeVPlusCandidate = InvalidNode;
};

template<typename PixelType = AltitudeType>
struct RegionMergeWorkspace {
    std::vector<PixelId> pixelsC;
    PixelType requestedA{};
    PixelType a{};
    PixelType b{};
    RegionUpdateDirection direction = RegionUpdateDirection::Unchanged;
    bool isMaxtree = false;
    NodeId nodeCa = InvalidNode;

    std::vector<NodeId> mergeNodesCMinus;
    std::vector<NodeId> mergeNodesIntervalMinus;
    std::vector<RegionMergeLevelWorkspace<PixelType>> mergeNodesByLevel;
    std::vector<NodeId> frontierNodesAboveB;
    std::vector<PixelType> affectedThresholds;

    const std::vector<NodeId> *mergeNodesVMinus(PixelType level) const {
        for (const auto &bucket : mergeNodesByLevel) {
            if (bucket.level == level) {
                return &bucket.mergeNodesVMinus;
            }
        }
        return nullptr;
    }

    NodeId mergedNodeVPlusCandidate(PixelType level) const {
        for (const auto &bucket : mergeNodesByLevel) {
            if (bucket.level == level) {
                return bucket.mergedNodeVPlusCandidate;
            }
        }
        return InvalidNode;
    }

    std::string debugString() const {
        std::ostringstream out;

        auto appendNodes = [&out](const std::vector<NodeId> &nodes) {
            out << "{";
            for (size_t i = 0; i < nodes.size(); ++i) {
                if (i != 0) {
                    out << ",";
                }
                out << nodes[i];
            }
            out << "}";
        };

        out << "C={";
        for (size_t i = 0; i < pixelsC.size(); ++i) {
            if (i != 0) {
                out << ",";
            }
            out << pixelsC[i];
        }
        out << "}\n";
        out << "requestedA=" << static_cast<int>(requestedA)
            << ", a=" << static_cast<int>(a)
            << ", b=" << static_cast<int>(b)
            << ", nodeCa=" << nodeCa << "\n";

        out << "MergeNodes(C, -)=";
        appendNodes(mergeNodesCMinus);
        out << "\n";

        out << "MergeNodes([a,b], -)=";
        appendNodes(mergeNodesIntervalMinus);
        out << "\n";

        for (const auto &bucket : mergeNodesByLevel) {
            out << "MergeNodes(v, -)[v=" << static_cast<int>(bucket.level) << "]=";
            appendNodes(bucket.mergeNodesVMinus);
            out << "\n";
            out << "MergedNode(v, +)[v=" << static_cast<int>(bucket.level) << "]="
                << bucket.mergedNodeVPlusCandidate << "\n";
        }

        out << "frontierNodesAboveB=";
        appendNodes(frontierNodesAboveB);
        out << "\n";

        out << "affectedThresholds={";
        for (size_t i = 0; i < affectedThresholds.size(); ++i) {
            if (i != 0) {
                out << ",";
            }
            out << static_cast<int>(affectedThresholds[i]);
        }
        out << "}\n";
        return out.str();
    }
};

template<typename PixelType = AltitudeType>
RegionMergeWorkspace<PixelType> buildRegionMergeWorkspace(const DynamicComponentTree &dualTree,
                                                          AdjacencyRelation &graph,
                                                          const RegionLevelChange<PixelType> &change) {
    if (change.pixelsC.empty()) {
        throw std::runtime_error("buildRegionMergeWorkspace requires a non-empty region C");
    }

    const int domainSize = dualTree.getNumTotalProperParts();
    std::vector<unsigned char> pixelsInC(static_cast<size_t>(domainSize), 0);
    for (PixelId pixelId : change.pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("buildRegionMergeWorkspace received an out-of-bounds pixel in C");
        }
        if (pixelsInC[static_cast<size_t>(pixelId)] != 0) {
            throw std::runtime_error("buildRegionMergeWorkspace received duplicated pixels in C");
        }
        pixelsInC[static_cast<size_t>(pixelId)] = 1;
    }
    if (!change.assumeConnected &&
        !isConnectedRegion(graph, std::span<const PixelId>(change.pixelsC))) {
        throw std::runtime_error("buildRegionMergeWorkspace requires connected region C");
    }

    auto [nodeCa, altitudeCa] = findNodeCa<PixelType>(dualTree, std::span<const PixelId>(change.pixelsC));

    RegionMergeWorkspace<PixelType> workspace;
    workspace.pixelsC = change.pixelsC;
    workspace.requestedA = change.a;
    workspace.a = altitudeCa;
    workspace.b = change.b;
    workspace.direction = directionFromLevels(altitudeCa, change.b);
    workspace.isMaxtree = dualTree.isMaxtree();
    workspace.nodeCa = nodeCa;

    const int nodeSlots = dualTree.getNumInternalNodeSlots();
    std::vector<unsigned char> adjacentSeedMarks(static_cast<size_t>(nodeSlots), 0);
    std::vector<unsigned char> climbedNodeMarks(static_cast<size_t>(nodeSlots), 0);
    std::vector<unsigned char> collectedNodeMarks(static_cast<size_t>(nodeSlots), 0);
    std::map<int, std::vector<NodeId>> bucketsByLevel;

    auto markNode = [nodeSlots](std::vector<unsigned char> &marks, NodeId nodeId) -> bool {
        if (nodeId < 0 || nodeId >= nodeSlots) {
            return false;
        }
        const size_t index = static_cast<size_t>(nodeId);
        if (marks[index] != 0) {
            return false;
        }
        marks[index] = 1;
        return true;
    };

    auto isMarked = [nodeSlots](const std::vector<unsigned char> &marks, NodeId nodeId) -> bool {
        return nodeId >= 0 &&
               nodeId < nodeSlots &&
               marks[static_cast<size_t>(nodeId)] != 0;
    };

    auto addMergeNode = [&](NodeId nodeId) {
        if (!markNode(collectedNodeMarks, nodeId)) {
            return;
        }
        const int level = dualTree.getAltitude(nodeId);
        bucketsByLevel[level].push_back(nodeId);
        workspace.mergeNodesIntervalMinus.push_back(nodeId);
    };

    auto addFrontierNode = [&](NodeId nodeId) {
        if (markNode(collectedNodeMarks, nodeId)) {
            workspace.frontierNodesAboveB.push_back(nodeId);
        }
    };

    for (PixelId pixelInC : change.pixelsC) {
        for (PixelId neighbor : graph.getNeighborPixels(pixelInC)) {
            if (neighbor < 0 || neighbor >= domainSize) {
                continue;
            }
            if (pixelsInC[static_cast<size_t>(neighbor)] != 0) {
                continue;
            }

            const NodeId seed = dualTree.getSmallestComponent(neighbor);
            if (seed == InvalidNode || !dualTree.isAlive(seed)) {
                continue;
            }

            const PixelType seedAltitude = static_cast<PixelType>(dualTree.getAltitude(seed));
            const bool validSeed = (workspace.isMaxtree && seedAltitude >= altitudeCa) ||
                                   (!workspace.isMaxtree && seedAltitude <= altitudeCa);
            if (!validSeed) {
                continue;
            }

            if (!markNode(adjacentSeedMarks, seed)) {
                continue;
            }
            workspace.mergeNodesCMinus.push_back(seed);

            NodeId current = seed;
            while (current != InvalidNode &&
                   dualTree.isAlive(current) &&
                   !isMarked(climbedNodeMarks, current)) {
                const PixelType currentLevel = static_cast<PixelType>(dualTree.getAltitude(current));
                const bool insideNodeCaInterval = (workspace.isMaxtree && currentLevel >= altitudeCa) ||
                                                  (!workspace.isMaxtree && currentLevel <= altitudeCa);
                if (!insideNodeCaInterval) {
                    break;
                }

                markNode(climbedNodeMarks, current);

                const bool insideMergeInterval = (workspace.isMaxtree && currentLevel <= change.b) ||
                                                 (!workspace.isMaxtree && currentLevel >= change.b);
                if (insideMergeInterval) {
                    addMergeNode(current);
                } else {
                    NodeId parentId = dualTree.getNodeParent(current);
                    if (parentId == current) {
                        parentId = InvalidNode;
                    }
                    const bool parentStillAboveB =
                        parentId != InvalidNode &&
                        ((workspace.isMaxtree && dualTree.getAltitude(parentId) > change.b) ||
                         (!workspace.isMaxtree && dualTree.getAltitude(parentId) < change.b));
                    if (!parentStillAboveB) {
                        addFrontierNode(current);
                    }
                }

                const NodeId parentId = dualTree.getNodeParent(current);
                if (parentId == current) {
                    break;
                }
                current = parentId;
            }
        }
    }

    std::vector<int> levels;
    levels.reserve(bucketsByLevel.size());
    for (const auto &[level, nodes] : bucketsByLevel) {
        static_cast<void>(nodes);
        levels.push_back(level);
    }
    if (workspace.isMaxtree) {
        std::reverse(levels.begin(), levels.end());
    }

    for (int level : levels) {
        RegionMergeLevelWorkspace<PixelType> bucket;
        bucket.level = static_cast<PixelType>(level);
        bucket.mergeNodesVMinus = bucketsByLevel[level];
        if (!bucket.mergeNodesVMinus.empty()) {
            bucket.mergedNodeVPlusCandidate = bucket.mergeNodesVMinus.front();
        }
        workspace.affectedThresholds.push_back(bucket.level);
        workspace.mergeNodesByLevel.push_back(std::move(bucket));
    }

    return workspace;
}

} // namespace debug::region
