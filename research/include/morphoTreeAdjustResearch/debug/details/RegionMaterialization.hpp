#pragma once

#include <span>
#include <stdexcept>
#include <vector>

#include "RegionLevelChange.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"
#include "RegionMaterializationPreflight.hpp"

namespace debug::region {

struct LocalMaterializationResult {
    NodeId originalOwner = InvalidNode;
    std::vector<NodeId> createdInsideNodes;
    std::vector<NodeId> createdOutsideNodes;

    int numCreatedNodes() const {
        return static_cast<int>(createdInsideNodes.size() + createdOutsideNodes.size());
    }
};

struct RegionMaterializationResult {
    std::vector<LocalMaterializationResult> partialOwnerResults;

    int numMaterializedOwners() const {
        return static_cast<int>(partialOwnerResults.size());
    }

    int numCreatedNodes() const {
        int count = 0;
        for (const auto &result : partialOwnerResults) {
            count += result.numCreatedNodes();
        }
        return count;
    }
};

template<typename PixelType = AltitudeType>
RegionMaterializationPlan analyzeRegionMaterialization(const DynamicComponentTree &tree,
                                                       AdjacencyRelation &graph,
                                                       const RegionLevelChange<PixelType> &change) {
    return ::debug::region::analyzeRegionMaterialization(
        tree,
        graph,
        std::span<const PixelId>(change.pixelsC));
}

inline NodeId createMaterializedChildFromComponent(DynamicComponentTree &tree,
                                                   const RegionOwnerMaterialization &owner,
                                                   std::span<const PixelId> component) {
    if (component.empty()) {
        throw std::runtime_error("createMaterializedChildFromComponent requires a non-empty component");
    }
    if (owner.nodeId == InvalidNode || !tree.isAlive(owner.nodeId)) {
        throw std::runtime_error("createMaterializedChildFromComponent requires a live owner");
    }
    if (tree.getAltitude(owner.nodeId) != owner.altitude) {
        throw std::runtime_error("createMaterializedChildFromComponent received a stale materialization plan");
    }

    for (PixelId pixelId : component) {
        if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
            throw std::runtime_error("createMaterializedChildFromComponent received an out-of-bounds pixel");
        }
        if (tree.getSmallestComponent(pixelId) != owner.nodeId) {
            throw std::runtime_error("createMaterializedChildFromComponent can only move direct proper parts of the owner");
        }
    }

    const NodeId newNodeId = tree.createDetachedNode(owner.altitude);
    if (newNodeId == InvalidNode) {
        throw std::runtime_error("createMaterializedChildFromComponent could not allocate a new node");
    }

    for (PixelId pixelId : component) {
        tree.moveProperPart(newNodeId, owner.nodeId, pixelId);
    }
    tree.attachNode(owner.nodeId, newNodeId);
    return newNodeId;
}

inline LocalMaterializationResult materializeSinglePartialOwner(DynamicComponentTree &tree,
                                                                const RegionOwnerMaterialization &owner) {
    if (!owner.cutsDirectProperParts()) {
        throw std::runtime_error("materializeSinglePartialOwner requires a partial owner");
    }
    if (owner.insideComponents.empty()) {
        throw std::runtime_error("materializeSinglePartialOwner requires at least one inside component");
    }
    if (owner.outsideComponents.empty()) {
        throw std::runtime_error("materializeSinglePartialOwner requires at least one outside component");
    }

    LocalMaterializationResult result;
    result.originalOwner = owner.nodeId;

    for (const auto &insideComponent : owner.insideComponents) {
        result.createdInsideNodes.push_back(createMaterializedChildFromComponent(
            tree,
            owner,
            std::span<const PixelId>(insideComponent)));
    }

    for (size_t i = 1; i < owner.outsideComponents.size(); ++i) {
        result.createdOutsideNodes.push_back(createMaterializedChildFromComponent(
            tree,
            owner,
            std::span<const PixelId>(owner.outsideComponents[i])));
    }

    return result;
}

inline LocalMaterializationResult materializeSinglePartialOwnerInsideOnly(
    DynamicComponentTree &tree,
    const RegionOwnerMaterialization &owner) {
    if (!owner.cutsDirectProperParts()) {
        throw std::runtime_error("materializeSinglePartialOwnerInsideOnly requires a partial owner");
    }
    if (owner.insideComponents.empty()) {
        throw std::runtime_error("materializeSinglePartialOwnerInsideOnly requires at least one inside component");
    }
    if (owner.outsideComponents.empty()) {
        throw std::runtime_error("materializeSinglePartialOwnerInsideOnly requires at least one outside component");
    }

    LocalMaterializationResult result;
    result.originalOwner = owner.nodeId;

    for (const auto &insideComponent : owner.insideComponents) {
        result.createdInsideNodes.push_back(createMaterializedChildFromComponent(
            tree,
            owner,
            std::span<const PixelId>(insideComponent)));
    }

    return result;
}

inline LocalMaterializationResult materializeSinglePartialOwner(DynamicComponentTree &tree,
                                                                const RegionMaterializationPlan &plan) {
    if (plan.numPartialOwners() != 1) {
        throw std::runtime_error("materializeSinglePartialOwner requires exactly one partial owner");
    }

    for (const auto &owner : plan.affectedOwners) {
        if (owner.cutsDirectProperParts()) {
            return materializeSinglePartialOwner(tree, owner);
        }
    }

    throw std::runtime_error("materializeSinglePartialOwner could not locate the partial owner");
}

inline LocalMaterializationResult materializeSinglePartialOwnerInsideOnly(
    DynamicComponentTree &tree,
    const RegionMaterializationPlan &plan) {
    if (plan.numPartialOwners() != 1) {
        throw std::runtime_error("materializeSinglePartialOwnerInsideOnly requires exactly one partial owner");
    }

    for (const auto &owner : plan.affectedOwners) {
        if (owner.cutsDirectProperParts()) {
            return materializeSinglePartialOwnerInsideOnly(tree, owner);
        }
    }

    throw std::runtime_error("materializeSinglePartialOwnerInsideOnly could not locate the partial owner");
}

inline RegionMaterializationResult materializePartialOwners(DynamicComponentTree &tree,
                                                            const RegionMaterializationPlan &plan) {
    if (plan.numPartialOwners() == 0) {
        throw std::runtime_error("materializePartialOwners requires at least one partial owner");
    }

    RegionMaterializationResult result;
    result.partialOwnerResults.reserve(static_cast<size_t>(plan.numPartialOwners()));
    for (const auto &owner : plan.affectedOwners) {
        if (owner.cutsDirectProperParts()) {
            result.partialOwnerResults.push_back(materializeSinglePartialOwner(tree, owner));
        }
    }
    return result;
}

inline RegionMaterializationResult materializePartialOwnersInsideOnly(
    DynamicComponentTree &tree,
    const RegionMaterializationPlan &plan) {
    if (plan.numPartialOwners() == 0) {
        throw std::runtime_error("materializePartialOwnersInsideOnly requires at least one partial owner");
    }

    RegionMaterializationResult result;
    result.partialOwnerResults.reserve(static_cast<size_t>(plan.numPartialOwners()));
    for (const auto &owner : plan.affectedOwners) {
        if (owner.cutsDirectProperParts()) {
            result.partialOwnerResults.push_back(
                materializeSinglePartialOwnerInsideOnly(tree, owner));
        }
    }
    return result;
}

} // namespace debug::region
