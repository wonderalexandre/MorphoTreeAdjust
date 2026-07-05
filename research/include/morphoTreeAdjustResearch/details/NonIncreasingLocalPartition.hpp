#pragma once

#include <algorithm>
#include <cstdlib>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../AdjacencyRelation.hpp"
#include "../Common.hpp"
#include "NonIncreasingLocalOracle.hpp"

namespace morphotree::details::non_increasing {

struct LocalThresholdComponent {
    AltitudeType level = 0;
    std::vector<PixelId> pixels;
    bool touchesModified = false;
    bool touchesPreserved = false;
    bool touchesBoundary = false;
    std::vector<NodeId> preDualOwners;
    std::vector<NodeId> oracleDualOwners;
    std::vector<NodeId> preThresholdOwners;
    std::vector<NodeId> oracleThresholdOwners;

    bool isRecomposedComponent() const {
        return touchesModified;
    }
};

struct LocalLevelPartition {
    AltitudeType level = 0;
    std::vector<LocalThresholdComponent> components;
};

struct LocalPartitionPlan {
    std::vector<AltitudeType> affectedLevels;
    std::vector<LocalLevelPartition> levels;
    std::vector<LocalThresholdComponent> recomposedComponents;
};

struct LocalRecompositionNodePlan {
    AltitudeType level = 0;
    std::vector<PixelId> supportPixels;
    std::vector<PixelId> directProperParts;
    bool touchesModified = false;
    bool touchesPreserved = false;
    bool touchesBoundary = false;
    std::vector<NodeId> preDualOwners;
    std::vector<NodeId> oracleDualOwners;
    NodeId oracleNodeId = InvalidNode;
    NodeId oracleParentId = InvalidNode;
    std::vector<NodeId> oracleChildrenInSnapshot;
    std::vector<NodeId> stableChildPreNodes;
    int parentPlanIndex = -1;
    std::vector<int> childPlanIndices;

    bool matchesOracleNode() const {
        return oracleNodeId != InvalidNode;
    }
};

struct LocalRecompositionPlan {
    std::vector<AltitudeType> affectedLevels;
    std::vector<LocalRecompositionNodePlan> nodes;

    bool allNodesMatchOracle() const {
        return std::all_of(nodes.begin(), nodes.end(), [](const auto &node) {
            return node.matchesOracleNode();
        });
    }
};

inline std::vector<AltitudeType> inclusiveLevelRange(int first, int last) {
    std::vector<AltitudeType> levels;
    if (first > last) {
        return levels;
    }
    levels.reserve(static_cast<std::size_t>(last - first + 1));
    for (int level = first; level <= last; ++level) {
        levels.push_back(static_cast<AltitudeType>(level));
    }
    return levels;
}

inline std::vector<AltitudeType> affectedLevelsForLocalMerge(
    const LocalMergeNodeContractionDescriptor &descriptor) {
    const int a = static_cast<int>(descriptor.a);
    const int b = static_cast<int>(descriptor.b);

    if (descriptor.direction == ::debug::region::RegionUpdateDirection::Unchanged) {
        return {};
    }

    if (descriptor.dualIsMaxtree) {
        if (descriptor.direction == ::debug::region::RegionUpdateDirection::Raise) {
            return inclusiveLevelRange(a + 1, b);
        }
        return inclusiveLevelRange(b + 1, a);
    }

    if (descriptor.direction == ::debug::region::RegionUpdateDirection::Lower) {
        return inclusiveLevelRange(b, a - 1);
    }
    return inclusiveLevelRange(a, b - 1);
}

inline bool pixelPassesDualThreshold(int value, AltitudeType level, bool dualIsMaxtree) {
    if (dualIsMaxtree) {
        return value >= level;
    }
    return value <= level;
}

inline std::vector<NodeId> sortedUniqueNodes(std::vector<NodeId> nodes) {
    nodes.erase(std::remove(nodes.begin(), nodes.end(), InvalidNode), nodes.end());
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
    return nodes;
}

inline bool equalPixelSets(std::span<const PixelId> lhs, std::span<const PixelId> rhs) {
    return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

inline bool containsAllPixels(std::span<const PixelId> container, std::span<const PixelId> contained) {
    return std::includes(container.begin(), container.end(), contained.begin(), contained.end());
}

inline const LocalDualNodeSnapshot *findPartitionSnapshotByNodeId(
    std::span<const LocalDualNodeSnapshot> snapshots,
    NodeId nodeId) {
    for (const auto &snapshot : snapshots) {
        if (snapshot.nodeId == nodeId) {
            return &snapshot;
        }
    }
    return nullptr;
}

inline NodeId thresholdOwnerInLocalSnapshot(std::span<const LocalDualNodeSnapshot> snapshots,
                                            NodeId smallestOwner,
                                            AltitudeType level,
                                            bool dualIsMaxtree) {
    const auto *current = findPartitionSnapshotByNodeId(snapshots, smallestOwner);
    if (current == nullptr ||
        !pixelPassesDualThreshold(current->altitude, level, dualIsMaxtree)) {
        return InvalidNode;
    }

    NodeId thresholdOwner = current->nodeId;
    while (current->parentInSnapshot) {
        const auto *parent = findPartitionSnapshotByNodeId(snapshots, current->parentId);
        if (parent == nullptr ||
            !pixelPassesDualThreshold(parent->altitude, level, dualIsMaxtree)) {
            break;
        }
        thresholdOwner = parent->nodeId;
        current = parent;
    }
    return thresholdOwner;
}

inline LocalThresholdComponent buildLocalThresholdComponent(
    const LocalMergeOracle &oracle,
    AltitudeType level,
    std::vector<PixelId> pixels,
    const std::vector<unsigned char> &roleKnown,
    const std::vector<LocalMergePixelRole> &roleByPixel,
    const std::vector<NodeId> &preOwnerByPixel,
    const std::vector<NodeId> &oracleOwnerByPixel) {
    LocalThresholdComponent component;
    component.level = level;
    component.pixels = sortedUniquePixels(std::move(pixels));

    for (PixelId pixelId : component.pixels) {
        const auto index = static_cast<std::size_t>(pixelId);
        if (index >= roleKnown.size() || roleKnown[index] == 0) {
            throw std::runtime_error("local partition component contains a pixel without local owner metadata");
        }

        const auto role = roleByPixel[index];
        component.touchesModified = component.touchesModified || role == LocalMergePixelRole::Modified;
        component.touchesPreserved = component.touchesPreserved || role == LocalMergePixelRole::Preserved;
        component.touchesBoundary = component.touchesBoundary || role == LocalMergePixelRole::Boundary;
        const NodeId preOwner = preOwnerByPixel[index];
        const NodeId oracleOwner = oracleOwnerByPixel[index];
        component.preDualOwners.push_back(preOwner);
        component.oracleDualOwners.push_back(oracleOwner);
        component.preThresholdOwners.push_back(thresholdOwnerInLocalSnapshot(
            std::span<const LocalDualNodeSnapshot>(oracle.preDualNodes),
            preOwner,
            level,
            oracle.contraction.dualIsMaxtree));
        component.oracleThresholdOwners.push_back(thresholdOwnerInLocalSnapshot(
            std::span<const LocalDualNodeSnapshot>(oracle.oracleDualNodes),
            oracleOwner,
            level,
            oracle.contraction.dualIsMaxtree));
    }

    component.preDualOwners = sortedUniqueNodes(std::move(component.preDualOwners));
    component.oracleDualOwners = sortedUniqueNodes(std::move(component.oracleDualOwners));
    component.preThresholdOwners = sortedUniqueNodes(std::move(component.preThresholdOwners));
    component.oracleThresholdOwners = sortedUniqueNodes(std::move(component.oracleThresholdOwners));
    return component;
}

inline bool shareAnyNode(std::span<const NodeId> lhs, std::span<const NodeId> rhs) {
    auto lhsIt = lhs.begin();
    auto rhsIt = rhs.begin();
    while (lhsIt != lhs.end() && rhsIt != rhs.end()) {
        if (*lhsIt == *rhsIt) {
            return true;
        }
        if (*lhsIt < *rhsIt) {
            ++lhsIt;
        } else {
            ++rhsIt;
        }
    }
    return false;
}

inline void appendNodes(std::vector<NodeId> &target, std::span<const NodeId> source) {
    target.insert(target.end(), source.begin(), source.end());
}

inline LocalThresholdComponent mergeLocalThresholdComponents(LocalThresholdComponent lhs,
                                                             const LocalThresholdComponent &rhs) {
    if (lhs.level != rhs.level) {
        throw std::runtime_error("local threshold component merge requires equal levels");
    }
    appendPixels(lhs.pixels, std::span<const PixelId>(rhs.pixels));
    lhs.pixels = sortedUniquePixels(std::move(lhs.pixels));
    lhs.touchesModified = lhs.touchesModified || rhs.touchesModified;
    lhs.touchesPreserved = lhs.touchesPreserved || rhs.touchesPreserved;
    lhs.touchesBoundary = lhs.touchesBoundary || rhs.touchesBoundary;
    appendNodes(lhs.preDualOwners, std::span<const NodeId>(rhs.preDualOwners));
    appendNodes(lhs.oracleDualOwners, std::span<const NodeId>(rhs.oracleDualOwners));
    appendNodes(lhs.preThresholdOwners, std::span<const NodeId>(rhs.preThresholdOwners));
    appendNodes(lhs.oracleThresholdOwners, std::span<const NodeId>(rhs.oracleThresholdOwners));
    lhs.preDualOwners = sortedUniqueNodes(std::move(lhs.preDualOwners));
    lhs.oracleDualOwners = sortedUniqueNodes(std::move(lhs.oracleDualOwners));
    lhs.preThresholdOwners = sortedUniqueNodes(std::move(lhs.preThresholdOwners));
    lhs.oracleThresholdOwners = sortedUniqueNodes(std::move(lhs.oracleThresholdOwners));
    return lhs;
}

inline bool componentMatchesOracleSnapshot(const LocalMergeOracle &oracle,
                                           const LocalThresholdComponent &component) {
    for (const auto &node : oracle.oracleDualNodes) {
        if (node.altitude != static_cast<int>(component.level)) {
            continue;
        }
        if (equalPixelSets(std::span<const PixelId>(node.supportLocalPixels),
                           std::span<const PixelId>(component.pixels))) {
            return true;
        }
    }
    return false;
}

inline bool mergedComponentMatchesOracleSnapshot(const LocalMergeOracle &oracle,
                                                 const LocalThresholdComponent &lhs,
                                                 const LocalThresholdComponent &rhs) {
    std::vector<PixelId> mergedPixels = lhs.pixels;
    appendPixels(mergedPixels, std::span<const PixelId>(rhs.pixels));
    mergedPixels = sortedUniquePixels(std::move(mergedPixels));

    for (const auto &node : oracle.oracleDualNodes) {
        if (node.altitude != static_cast<int>(lhs.level) ||
            !node.supportTouchesModified) {
            continue;
        }
        if (equalPixelSets(std::span<const PixelId>(node.supportLocalPixels),
                           std::span<const PixelId>(mergedPixels))) {
            return true;
        }
    }
    return false;
}

inline bool mergeFirstPairwiseOracleMatchedComponentSet(
    const LocalMergeOracle &oracle,
    std::vector<LocalThresholdComponent> &components) {
    for (std::size_t i = 0; i < components.size(); ++i) {
        for (std::size_t j = i + 1; j < components.size(); ++j) {
            if (!shareAnyNode(
                    std::span<const NodeId>(components[i].oracleThresholdOwners),
                    std::span<const NodeId>(components[j].oracleThresholdOwners))) {
                continue;
            }
            if (!mergedComponentMatchesOracleSnapshot(oracle, components[i], components[j])) {
                continue;
            }
            components[i] = mergeLocalThresholdComponents(
                std::move(components[i]),
                components[j]);
            components.erase(components.begin() + static_cast<std::ptrdiff_t>(j));
            return true;
        }
    }
    return false;
}

inline bool mergeFirstOracleSnapshotCoveredComponentSet(
    const LocalMergeOracle &oracle,
    std::vector<LocalThresholdComponent> &components) {
    for (const auto &node : oracle.oracleDualNodes) {
        if (!node.supportTouchesModified) {
            continue;
        }

        std::vector<std::size_t> coveredIndices;
        std::vector<PixelId> mergedPixels;
        bool touchesModified = false;
        bool modifiedComponentAlreadyMatchesOracle = false;
        for (std::size_t index = 0; index < components.size(); ++index) {
            const auto &component = components[index];
            if (node.altitude != static_cast<int>(component.level) ||
                !containsAllPixels(std::span<const PixelId>(node.supportLocalPixels),
                                   std::span<const PixelId>(component.pixels))) {
                continue;
            }
            coveredIndices.push_back(index);
            appendPixels(mergedPixels, std::span<const PixelId>(component.pixels));
            touchesModified = touchesModified || component.touchesModified;
            modifiedComponentAlreadyMatchesOracle =
                modifiedComponentAlreadyMatchesOracle ||
                (component.touchesModified && componentMatchesOracleSnapshot(oracle, component));
        }

        if (coveredIndices.size() < 2 ||
            !touchesModified ||
            modifiedComponentAlreadyMatchesOracle) {
            continue;
        }

        mergedPixels = sortedUniquePixels(std::move(mergedPixels));
        if (!equalPixelSets(std::span<const PixelId>(node.supportLocalPixels),
                            std::span<const PixelId>(mergedPixels))) {
            continue;
        }

        LocalThresholdComponent merged = components[coveredIndices.front()];
        for (std::size_t i = 1; i < coveredIndices.size(); ++i) {
            merged = mergeLocalThresholdComponents(
                std::move(merged),
                components[coveredIndices[i]]);
        }

        components[coveredIndices.front()] = std::move(merged);
        for (auto it = coveredIndices.rbegin(); it != coveredIndices.rend(); ++it) {
            if (*it == coveredIndices.front()) {
                continue;
            }
            components.erase(components.begin() + static_cast<std::ptrdiff_t>(*it));
        }
        return true;
    }
    return false;
}

inline void mergeExternallyConnectedLocalThresholdComponents(
    const LocalMergeOracle &oracle,
    std::vector<LocalThresholdComponent> &components) {
    while (mergeFirstPairwiseOracleMatchedComponentSet(oracle, components)) {
    }
    while (mergeFirstOracleSnapshotCoveredComponentSet(oracle, components)) {
    }
}

inline std::vector<LocalThresholdComponent> collectLocalThresholdComponents(
    const LocalMergeOracle &oracle,
    AdjacencyRelation &adjacency,
    AltitudeType level,
    const std::vector<unsigned char> &localMask,
    const std::vector<unsigned char> &roleKnown,
    const std::vector<LocalMergePixelRole> &roleByPixel,
    const std::vector<NodeId> &preOwnerByPixel,
    const std::vector<NodeId> &oracleOwnerByPixel) {
    const int domainSize = oracle.updatedImage->getSize();
    std::vector<unsigned char> allowed(static_cast<std::size_t>(domainSize), 0);
    std::vector<unsigned char> visited(static_cast<std::size_t>(domainSize), 0);
    std::vector<PixelId> seeds;
    seeds.reserve(oracle.localPixels.size());

    for (PixelId pixelId : oracle.localPixels) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("local partition received an out-of-bounds local pixel");
        }
        if (pixelPassesDualThreshold((*oracle.updatedImage)[pixelId],
                                     level,
                                     oracle.contraction.dualIsMaxtree)) {
            allowed[static_cast<std::size_t>(pixelId)] = 1;
            seeds.push_back(pixelId);
        }
    }

    std::vector<LocalThresholdComponent> components;
    std::vector<PixelId> stack;
    stack.reserve(seeds.size());

    for (PixelId seed : seeds) {
        const auto seedIndex = static_cast<std::size_t>(seed);
        if (visited[seedIndex] != 0) {
            continue;
        }

        std::vector<PixelId> pixels;
        visited[seedIndex] = 1;
        stack.push_back(seed);

        while (!stack.empty()) {
            const PixelId current = stack.back();
            stack.pop_back();
            pixels.push_back(current);

            for (PixelId neighbor : adjacency.getNeighborPixels(current)) {
                if (neighbor < 0 || neighbor >= domainSize) {
                    continue;
                }
                const auto neighborIndex = static_cast<std::size_t>(neighbor);
                if (localMask[neighborIndex] == 0 ||
                    allowed[neighborIndex] == 0 ||
                    visited[neighborIndex] != 0) {
                    continue;
                }
                visited[neighborIndex] = 1;
                stack.push_back(neighbor);
            }
        }

        components.push_back(buildLocalThresholdComponent(
            oracle,
            level,
            std::move(pixels),
            roleKnown,
            roleByPixel,
            preOwnerByPixel,
            oracleOwnerByPixel));
    }

    mergeExternallyConnectedLocalThresholdComponents(oracle, components);
    std::sort(components.begin(), components.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.pixels < rhs.pixels;
    });
    return components;
}

inline LocalPartitionPlan buildLocalPartitionPlan(const LocalMergeOracle &oracle,
                                                  AdjacencyRelation &adjacency) {
    if (!oracle.updatedImage) {
        throw std::runtime_error("buildLocalPartitionPlan requires a valid updated image");
    }
    if (oracle.localPixels.empty()) {
        throw std::runtime_error("buildLocalPartitionPlan requires non-empty local pixels");
    }

    const int domainSize = oracle.updatedImage->getSize();
    const auto localMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(oracle.localPixels),
        "local partition");

    std::vector<unsigned char> roleKnown(static_cast<std::size_t>(domainSize), 0);
    std::vector<LocalMergePixelRole> roleByPixel(
        static_cast<std::size_t>(domainSize),
        LocalMergePixelRole::Boundary);
    std::vector<NodeId> preOwnerByPixel(static_cast<std::size_t>(domainSize), InvalidNode);
    std::vector<NodeId> oracleOwnerByPixel(static_cast<std::size_t>(domainSize), InvalidNode);

    for (const auto &owner : oracle.pixelOwners) {
        if (owner.pixelId < 0 || owner.pixelId >= domainSize) {
            throw std::runtime_error("buildLocalPartitionPlan received an out-of-bounds pixel owner");
        }
        const auto index = static_cast<std::size_t>(owner.pixelId);
        if (roleKnown[index] != 0) {
            throw std::runtime_error("buildLocalPartitionPlan received duplicated pixel owner metadata");
        }
        roleKnown[index] = 1;
        roleByPixel[index] = owner.role;
        preOwnerByPixel[index] = owner.preDualOwner;
        oracleOwnerByPixel[index] = owner.oracleDualOwner;
    }

    for (PixelId pixelId : oracle.localPixels) {
        const auto index = static_cast<std::size_t>(pixelId);
        if (roleKnown[index] == 0) {
            throw std::runtime_error("buildLocalPartitionPlan requires owner metadata for every local pixel");
        }
    }

    LocalPartitionPlan plan;
    plan.affectedLevels = affectedLevelsForLocalMerge(oracle.contraction);
    plan.levels.reserve(plan.affectedLevels.size());

    for (AltitudeType level : plan.affectedLevels) {
        LocalLevelPartition partition;
        partition.level = level;
        partition.components = collectLocalThresholdComponents(
            oracle,
            adjacency,
            level,
            localMask,
            roleKnown,
            roleByPixel,
            preOwnerByPixel,
            oracleOwnerByPixel);

        for (const auto &component : partition.components) {
            if (component.isRecomposedComponent()) {
                plan.recomposedComponents.push_back(component);
            }
        }
        plan.levels.push_back(std::move(partition));
    }

    return plan;
}

inline const LocalDualNodeSnapshot *findOracleSnapshotForComponent(
    const LocalMergeOracle &oracle,
    const LocalThresholdComponent &component) {
    for (const auto &node : oracle.oracleDualNodes) {
        if (node.altitude != static_cast<int>(component.level)) {
            continue;
        }
        if (equalPixelSets(std::span<const PixelId>(node.supportLocalPixels),
                           std::span<const PixelId>(component.pixels))) {
            return &node;
        }
    }
    return nullptr;
}

inline std::vector<PixelId> collectDirectPixelsAtLevel(const LocalMergeOracle &oracle,
                                                       const LocalThresholdComponent &component) {
    std::vector<PixelId> directPixels;
    directPixels.reserve(component.pixels.size());
    for (PixelId pixelId : component.pixels) {
        if ((*oracle.updatedImage)[pixelId] == component.level) {
            directPixels.push_back(pixelId);
        }
    }
    return sortedUniquePixels(std::move(directPixels));
}

inline bool validLocalParentLevel(AltitudeType childLevel,
                                  AltitudeType parentLevel,
                                  bool dualIsMaxtree) {
    if (dualIsMaxtree) {
        return parentLevel < childLevel;
    }
    return parentLevel > childLevel;
}

inline void assignLocalRecompositionParentIndices(LocalRecompositionPlan &plan,
                                                  bool dualIsMaxtree) {
    for (auto &node : plan.nodes) {
        node.parentPlanIndex = -1;
        node.childPlanIndices.clear();
    }

    for (std::size_t childIndex = 0; childIndex < plan.nodes.size(); ++childIndex) {
        const auto &child = plan.nodes[childIndex];
        int bestParentIndex = -1;

        if (child.oracleParentId != InvalidNode) {
            for (std::size_t parentIndex = 0; parentIndex < plan.nodes.size(); ++parentIndex) {
                if (parentIndex == childIndex) {
                    continue;
                }
                const auto &candidate = plan.nodes[parentIndex];
                if (candidate.oracleNodeId != child.oracleParentId) {
                    continue;
                }
                if (!validLocalParentLevel(child.level, candidate.level, dualIsMaxtree)) {
                    continue;
                }
                if (!containsAllPixels(std::span<const PixelId>(candidate.supportPixels),
                                       std::span<const PixelId>(child.supportPixels))) {
                    continue;
                }
                bestParentIndex = static_cast<int>(parentIndex);
                break;
            }
        }

        if (bestParentIndex != -1) {
            plan.nodes[childIndex].parentPlanIndex = bestParentIndex;
            continue;
        }

        std::size_t bestParentSupportSize = 0;
        int bestParentLevelDistance = 0;

        for (std::size_t parentIndex = 0; parentIndex < plan.nodes.size(); ++parentIndex) {
            if (parentIndex == childIndex) {
                continue;
            }
            const auto &candidate = plan.nodes[parentIndex];
            if (!validLocalParentLevel(child.level, candidate.level, dualIsMaxtree)) {
                continue;
            }
            if (!containsAllPixels(std::span<const PixelId>(candidate.supportPixels),
                                   std::span<const PixelId>(child.supportPixels))) {
                continue;
            }
            if (candidate.supportPixels.size() == child.supportPixels.size()) {
                continue;
            }

            const int levelDistance = std::abs(static_cast<int>(candidate.level) - static_cast<int>(child.level));
            if (bestParentIndex == -1 ||
                candidate.supportPixels.size() < bestParentSupportSize ||
                (candidate.supportPixels.size() == bestParentSupportSize &&
                 levelDistance < bestParentLevelDistance)) {
                bestParentIndex = static_cast<int>(parentIndex);
                bestParentSupportSize = candidate.supportPixels.size();
                bestParentLevelDistance = levelDistance;
            }
        }

        plan.nodes[childIndex].parentPlanIndex = bestParentIndex;
    }

    for (std::size_t childIndex = 0; childIndex < plan.nodes.size(); ++childIndex) {
        const int parentIndex = plan.nodes[childIndex].parentPlanIndex;
        if (parentIndex >= 0) {
            plan.nodes[static_cast<std::size_t>(parentIndex)].childPlanIndices.push_back(
                static_cast<int>(childIndex));
        }
    }
}

inline LocalRecompositionPlan buildLocalRecompositionPlan(const LocalMergeOracle &oracle,
                                                          const LocalPartitionPlan &partitionPlan) {
    if (!oracle.updatedImage) {
        throw std::runtime_error("buildLocalRecompositionPlan requires a valid updated image");
    }

    LocalRecompositionPlan plan;
    plan.affectedLevels = partitionPlan.affectedLevels;

    for (const auto &partition : partitionPlan.levels) {
        for (const auto &component : partition.components) {
            if (!component.isRecomposedComponent()) {
                continue;
            }

            const LocalDualNodeSnapshot *oracleNode = findOracleSnapshotForComponent(oracle, component);
            std::vector<PixelId> directPixels = collectDirectPixelsAtLevel(oracle, component);
            if (oracleNode != nullptr && !oracleNode->directLocalProperParts.empty()) {
                directPixels = oracleNode->directLocalProperParts;
            }
            directPixels = sortedUniquePixels(std::move(directPixels));

            if (directPixels.empty() && oracleNode == nullptr) {
                continue;
            }

            LocalRecompositionNodePlan node;
            node.level = component.level;
            node.supportPixels = component.pixels;
            node.directProperParts = std::move(directPixels);
            node.touchesModified = component.touchesModified;
            node.touchesPreserved = component.touchesPreserved;
            node.touchesBoundary = component.touchesBoundary;
            node.preDualOwners = component.preDualOwners;
            node.oracleDualOwners = component.oracleDualOwners;

            if (oracleNode != nullptr) {
                node.oracleNodeId = oracleNode->nodeId;
                node.oracleParentId = oracleNode->parentId;
                node.oracleChildrenInSnapshot = oracleNode->childrenInSnapshot;
            }

            plan.nodes.push_back(std::move(node));
        }
    }

    std::sort(plan.nodes.begin(), plan.nodes.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.level != rhs.level) {
            return lhs.level < rhs.level;
        }
        return lhs.supportPixels < rhs.supportPixels;
    });
    assignLocalRecompositionParentIndices(plan, oracle.contraction.dualIsMaxtree);
    return plan;
}

inline std::string debugString(const LocalPartitionPlan &plan) {
    std::ostringstream out;
    out << "LocalPartitionPlan\n";
    out << "affected_levels=[";
    for (std::size_t i = 0; i < plan.affectedLevels.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        out << static_cast<int>(plan.affectedLevels[i]);
    }
    out << "]\n";
    out << "recomposed_components=" << plan.recomposedComponents.size() << "\n";

    for (const auto &partition : plan.levels) {
        out << "level=" << static_cast<int>(partition.level)
            << " components=" << partition.components.size()
            << "\n";
        for (const auto &component : partition.components) {
            out << "  pixels=";
            appendPixelVector(out, std::span<const PixelId>(component.pixels));
            out << " touches_modified=" << (component.touchesModified ? "true" : "false")
                << " touches_preserved=" << (component.touchesPreserved ? "true" : "false")
                << " touches_boundary=" << (component.touchesBoundary ? "true" : "false")
                << " pre_owners=";
            appendNodeVector(out, std::span<const NodeId>(component.preDualOwners));
            out << " oracle_owners=";
            appendNodeVector(out, std::span<const NodeId>(component.oracleDualOwners));
            out << " pre_threshold_owners=";
            appendNodeVector(out, std::span<const NodeId>(component.preThresholdOwners));
            out << " oracle_threshold_owners=";
            appendNodeVector(out, std::span<const NodeId>(component.oracleThresholdOwners));
            out << "\n";
        }
    }

    return out.str();
}

inline std::string debugString(const LocalRecompositionPlan &plan) {
    std::ostringstream out;
    out << "LocalRecompositionPlan\n";
    out << "affected_levels=[";
    for (std::size_t i = 0; i < plan.affectedLevels.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        out << static_cast<int>(plan.affectedLevels[i]);
    }
    out << "]\n";
    out << "nodes=" << plan.nodes.size()
        << " all_nodes_match_oracle=" << (plan.allNodesMatchOracle() ? "true" : "false")
        << "\n";

    for (std::size_t i = 0; i < plan.nodes.size(); ++i) {
        const auto &node = plan.nodes[i];
        out << "  plan_node=" << i
            << " level=" << static_cast<int>(node.level)
            << " parent_plan=" << node.parentPlanIndex
            << " oracle_node=" << node.oracleNodeId
            << " oracle_parent=" << node.oracleParentId
            << " support=";
        appendPixelVector(out, std::span<const PixelId>(node.supportPixels));
        out << " direct=";
        appendPixelVector(out, std::span<const PixelId>(node.directProperParts));
        out << " touches_modified=" << (node.touchesModified ? "true" : "false")
            << " touches_preserved=" << (node.touchesPreserved ? "true" : "false")
            << " touches_boundary=" << (node.touchesBoundary ? "true" : "false")
            << " pre_owners=";
        appendNodeVector(out, std::span<const NodeId>(node.preDualOwners));
        out << " oracle_owners=";
        appendNodeVector(out, std::span<const NodeId>(node.oracleDualOwners));
        out << " oracle_children=";
        appendNodeVector(out, std::span<const NodeId>(node.oracleChildrenInSnapshot));
        out << " stable_pre_children=";
        appendNodeVector(out, std::span<const NodeId>(node.stableChildPreNodes));
        out << " child_plan_indices=[";
        for (std::size_t childIndex = 0; childIndex < node.childPlanIndices.size(); ++childIndex) {
            if (childIndex != 0) {
                out << ",";
            }
            out << node.childPlanIndices[childIndex];
        }
        out << "]\n";
    }

    return out.str();
}

} // namespace morphotree::details::non_increasing
