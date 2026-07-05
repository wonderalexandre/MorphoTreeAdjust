#pragma once

#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::region {

struct RegionOwnerMaterialization {
    NodeId nodeId = InvalidNode;
    int altitude = 0;
    std::vector<PixelId> insideDirectProperParts;
    std::vector<PixelId> outsideDirectProperParts;
    std::vector<std::vector<PixelId>> insideComponents;
    std::vector<std::vector<PixelId>> outsideComponents;

    bool coversAllDirectProperParts() const {
        return !insideDirectProperParts.empty() && outsideDirectProperParts.empty();
    }

    bool cutsDirectProperParts() const {
        return !insideDirectProperParts.empty() && !outsideDirectProperParts.empty();
    }
};

struct RegionMaterializationPlan {
    std::vector<RegionOwnerMaterialization> affectedOwners;

    bool requiresLocalMaterialization() const {
        for (const auto &owner : affectedOwners) {
            if (owner.cutsDirectProperParts()) {
                return true;
            }
        }
        return false;
    }

    int numPartialOwners() const {
        int count = 0;
        for (const auto &owner : affectedOwners) {
            if (owner.cutsDirectProperParts()) {
                ++count;
            }
        }
        return count;
    }

    int numFullyCoveredOwners() const {
        int count = 0;
        for (const auto &owner : affectedOwners) {
            if (owner.coversAllDirectProperParts()) {
                ++count;
            }
        }
        return count;
    }
};

inline std::vector<std::vector<PixelId>> collectLocalConnectedComponents(
    AdjacencyRelation &graph,
    std::span<const PixelId> pixels,
    int domainSize) {
    std::vector<std::vector<PixelId>> components;
    if (pixels.empty()) {
        return components;
    }

    std::vector<unsigned char> allowed(static_cast<size_t>(domainSize), 0);
    std::vector<unsigned char> visited(static_cast<size_t>(domainSize), 0);
    for (PixelId pixelId : pixels) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("collectLocalConnectedComponents received an out-of-bounds pixel");
        }
        allowed[static_cast<size_t>(pixelId)] = 1;
    }

    std::vector<PixelId> stack;
    stack.reserve(pixels.size());
    for (PixelId seed : pixels) {
        if (visited[static_cast<size_t>(seed)] != 0) {
            continue;
        }

        std::vector<PixelId> component;
        stack.push_back(seed);
        visited[static_cast<size_t>(seed)] = 1;

        while (!stack.empty()) {
            const PixelId current = stack.back();
            stack.pop_back();
            component.push_back(current);

            for (PixelId neighbor : graph.getNeighborPixels(current)) {
                if (neighbor < 0 || neighbor >= domainSize) {
                    continue;
                }
                const size_t neighborIndex = static_cast<size_t>(neighbor);
                if (allowed[neighborIndex] == 0 || visited[neighborIndex] != 0) {
                    continue;
                }
                visited[neighborIndex] = 1;
                stack.push_back(neighbor);
            }
        }

        components.push_back(std::move(component));
    }

    return components;
}

inline RegionMaterializationPlan analyzeRegionMaterialization(
    const DynamicComponentTree &tree,
    AdjacencyRelation &graph,
    std::span<const PixelId> pixelsC) {
    RegionMaterializationPlan plan;
    if (pixelsC.empty()) {
        return plan;
    }

    const int domainSize = tree.getNumTotalProperParts();
    std::vector<unsigned char> inRegion(static_cast<size_t>(domainSize), 0);
    std::unordered_map<NodeId, size_t> ownerToEntry;
    ownerToEntry.reserve(pixelsC.size());

    for (PixelId pixelId : pixelsC) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("analyzeRegionMaterialization received an out-of-bounds pixel");
        }
        const size_t pixelIndex = static_cast<size_t>(pixelId);
        if (inRegion[pixelIndex] != 0) {
            throw std::runtime_error("analyzeRegionMaterialization received duplicated pixels in C");
        }
        inRegion[pixelIndex] = 1;

        const NodeId ownerNodeId = tree.getSmallestComponent(pixelId);
        if (ownerNodeId == InvalidNode || !tree.isAlive(ownerNodeId)) {
            throw std::runtime_error("analyzeRegionMaterialization cannot locate a live owner for a pixel in C");
        }

        if (ownerToEntry.find(ownerNodeId) == ownerToEntry.end()) {
            ownerToEntry.emplace(ownerNodeId, plan.affectedOwners.size());
            RegionOwnerMaterialization owner;
            owner.nodeId = ownerNodeId;
            owner.altitude = tree.getAltitude(ownerNodeId);
            plan.affectedOwners.push_back(std::move(owner));
        }
    }

    for (auto &owner : plan.affectedOwners) {
        for (PixelId pixelId : tree.getProperParts(owner.nodeId)) {
            if (inRegion[static_cast<size_t>(pixelId)] != 0) {
                owner.insideDirectProperParts.push_back(pixelId);
            } else {
                owner.outsideDirectProperParts.push_back(pixelId);
            }
        }

        owner.insideComponents = collectLocalConnectedComponents(
            graph,
            std::span<const PixelId>(owner.insideDirectProperParts),
            domainSize);
        owner.outsideComponents = collectLocalConnectedComponents(
            graph,
            std::span<const PixelId>(owner.outsideDirectProperParts),
            domainSize);
    }

    return plan;
}

} // namespace debug::region
