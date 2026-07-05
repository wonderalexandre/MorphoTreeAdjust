#pragma once

#include <algorithm>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../AdjacencyRelation.hpp"
#include "../Common.hpp"
#include "NonIncreasingLocalPartition.hpp"

namespace morphotree::details::non_increasing {

inline constexpr std::size_t CompactLocalPartitionPixelLimit = 16384;

struct LocalPartitionPixelMetadata {
    PixelId pixelId = -1;
    LocalMergePixelRole role = LocalMergePixelRole::Boundary;
    NodeId preDualOwner = InvalidNode;
};

struct LocalPartitionInput {
    ImageUInt8Ptr updatedImage;
    int domainSize = 0;
    LocalMergeNodeContractionDescriptor contraction;
    std::vector<PixelId> localPixels;
    std::vector<uint8_t> updatedLocalPixelValues;
    std::vector<LocalPartitionPixelMetadata> pixelMetadata;
    std::vector<LocalDualNodeSnapshot> preDualNodes;
};

struct LocalPartitionInputOptions {
    bool checkPreReconstruction = true;
    bool buildFullUpdatedImage = true;
};

inline unsigned char localRoleCode(LocalMergePixelRole role) {
    switch (role) {
    case LocalMergePixelRole::Modified:
        return 1;
    case LocalMergePixelRole::Preserved:
        return 2;
    case LocalMergePixelRole::Boundary:
        return 3;
    }
    return 0;
}

inline LocalMergePixelRole localRoleFromCode(unsigned char code) {
    if (code == localRoleCode(LocalMergePixelRole::Modified)) {
        return LocalMergePixelRole::Modified;
    }
    if (code == localRoleCode(LocalMergePixelRole::Preserved)) {
        return LocalMergePixelRole::Preserved;
    }
    return LocalMergePixelRole::Boundary;
}

inline LocalPartitionInput makeLocalPartitionInputFromOraclePreState(
    const LocalMergeOracle &oracle) {
    LocalPartitionInput input;
    input.updatedImage = oracle.updatedImage;
    input.domainSize = localMergeOracleDomainSize(oracle);
    input.contraction = oracle.contraction;
    input.localPixels = oracle.localPixels;
    input.updatedLocalPixelValues = oracle.updatedLocalPixelValues;
    input.preDualNodes = oracle.preDualNodes;
    input.pixelMetadata.reserve(oracle.pixelOwners.size());
    for (const auto &owner : oracle.pixelOwners) {
        LocalPartitionPixelMetadata metadata;
        metadata.pixelId = owner.pixelId;
        metadata.role = owner.role;
        metadata.preDualOwner = owner.preDualOwner;
        input.pixelMetadata.push_back(metadata);
    }
    return input;
}

inline uint8_t reconstructedPixelValueFromTree(const DynamicComponentTree &tree,
                                               PixelId pixelId) {
    if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
        throw std::runtime_error("reconstructedPixelValueFromTree received an out-of-bounds pixel");
    }
    const NodeId ownerId = tree.getSmallestComponent(pixelId);
    if (ownerId == InvalidNode || !tree.isNode(ownerId) || !tree.isAlive(ownerId)) {
        return 0;
    }
    return static_cast<uint8_t>(tree.getAltitude(ownerId));
}

inline ImageUInt8Ptr buildUpdatedImageForLocalPartition(
    const DynamicComponentTree &primalTree,
    const ImageUInt8Ptr &preImage,
    const LocalMergeNodeContractionDescriptor &contraction,
    std::span<const PixelId> localPixels,
    bool buildFullImage) {
    ImageUInt8Ptr updatedImage;
    if (buildFullImage) {
        if (!preImage) {
            throw std::runtime_error("full local partition update image requires a pre-update image");
        }
        updatedImage = applyLocalMergeDescriptorToImage(preImage, contraction);
        return updatedImage;
    }

    updatedImage = ImageUInt8::create(
        primalTree.getNumRowsOfImage(),
        primalTree.getNumColsOfImage());
    for (PixelId pixelId : localPixels) {
        (*updatedImage)[pixelId] = reconstructedPixelValueFromTree(primalTree, pixelId);
    }
    for (PixelId pixelId : contraction.modifiedProperPart) {
        (*updatedImage)[pixelId] = static_cast<uint8_t>(contraction.b);
    }
    return updatedImage;
}

inline std::vector<uint8_t> buildUpdatedLocalPixelValuesForLocalPartition(
    const DynamicComponentTree &primalTree,
    const LocalMergeNodeContractionDescriptor &contraction,
    std::span<const PixelId> localPixels) {
    auto modifiedProperPart = sortedUniquePixels(contraction.modifiedProperPart);
    std::vector<uint8_t> values;
    values.reserve(localPixels.size());
    for (PixelId pixelId : localPixels) {
        if (std::binary_search(modifiedProperPart.begin(),
                               modifiedProperPart.end(),
                               pixelId)) {
            values.push_back(static_cast<uint8_t>(contraction.b));
        } else {
            values.push_back(reconstructedPixelValueFromTree(primalTree, pixelId));
        }
    }
    return values;
}

inline int localPartitionDomainSize(const LocalPartitionInput &input) {
    if (input.updatedImage) {
        return input.updatedImage->getSize();
    }
    return input.domainSize;
}

inline bool localPartitionHasSparseUpdatedTarget(const LocalPartitionInput &input) {
    return input.domainSize > 0 &&
           input.updatedLocalPixelValues.size() == input.localPixels.size();
}

inline bool localPartitionHasUpdatedTarget(const LocalPartitionInput &input) {
    return input.updatedImage != nullptr ||
           localPartitionHasSparseUpdatedTarget(input);
}

inline uint8_t localPartitionUpdatedPixelValue(const LocalPartitionInput &input,
                                               PixelId pixelId) {
    if (input.updatedImage) {
        return (*input.updatedImage)[pixelId];
    }
    const auto localIndex = findSortedPixelIndex(
        std::span<const PixelId>(input.localPixels),
        pixelId);
    if (localIndex == input.localPixels.size() ||
        localIndex >= input.updatedLocalPixelValues.size()) {
        throw std::runtime_error("local partition has no sparse target value for the requested pixel");
    }
    return input.updatedLocalPixelValues[localIndex];
}

inline std::vector<PixelId> collectBoundaryPixelsAndMarkLocalRoles(
    AdjacencyRelation &adjacency,
    std::span<const PixelId> corePixels,
    std::vector<unsigned char> &roleByPixel) {
    std::vector<PixelId> boundary;
    for (PixelId pixelId : corePixels) {
        for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
            if (neighbor < 0 || static_cast<std::size_t>(neighbor) >= roleByPixel.size()) {
                continue;
            }
            const auto neighborIndex = static_cast<std::size_t>(neighbor);
            if (roleByPixel[neighborIndex] != 0) {
                continue;
            }
            roleByPixel[neighborIndex] = localRoleCode(LocalMergePixelRole::Boundary);
            boundary.push_back(neighbor);
        }
    }
    return sortedUniquePixels(std::move(boundary));
}

inline std::vector<PixelId> collectBoundaryPixelsAndMarkCompactLocalRoles(
    AdjacencyRelation &adjacency,
    std::span<const PixelId> corePixels,
    int domainSize,
    std::unordered_map<PixelId, unsigned char> &roleByPixel) {
    std::vector<PixelId> boundary;
    for (PixelId pixelId : corePixels) {
        for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
            if (neighbor < 0 || neighbor >= domainSize) {
                continue;
            }
            if (roleByPixel.find(neighbor) != roleByPixel.end()) {
                continue;
            }
            roleByPixel.emplace(neighbor, localRoleCode(LocalMergePixelRole::Boundary));
            boundary.push_back(neighbor);
        }
    }
    return sortedUniquePixels(std::move(boundary));
}

inline std::vector<LocalDualNodeSnapshot> buildLocalDualNodeSnapshotFromLocalRoles(
    const DynamicComponentTree &tree,
    std::span<const PixelId> localPixels,
    const std::vector<unsigned char> &roleByPixel) {
    const int nodeSlots = tree.getNumInternalNodeSlots();
    std::vector<unsigned char> nodeInSnapshot(static_cast<std::size_t>(nodeSlots), 0);
    std::vector<NodeId> snapshotNodeIds;
    snapshotNodeIds.reserve(localPixels.size() * 2);

    auto addNode = [&](NodeId nodeId) {
        if (nodeId != InvalidNode && tree.isNode(nodeId) && tree.isAlive(nodeId)) {
            const auto nodeIndex = static_cast<std::size_t>(nodeId);
            if (nodeInSnapshot[nodeIndex] == 0) {
                nodeInSnapshot[nodeIndex] = 1;
                snapshotNodeIds.push_back(nodeId);
            }
        }
    };

    std::unordered_map<NodeId, std::vector<PixelId>> directByNode;
    directByNode.reserve(localPixels.size());
    std::vector<NodeId> ownerByLocalIndex;
    ownerByLocalIndex.reserve(localPixels.size());
    for (PixelId pixelId : localPixels) {
        const NodeId ownerId = tree.getSmallestComponent(pixelId);
        ownerByLocalIndex.push_back(ownerId);
        addNode(ownerId);
        if (ownerId != InvalidNode && tree.isAlive(ownerId)) {
            directByNode[ownerId].push_back(pixelId);
            const NodeId parentId = tree.getNodeParent(ownerId);
            if (parentId != ownerId) {
                addNode(parentId);
            }
        }
    }

    std::unordered_map<NodeId, std::vector<PixelId>> supportByNode;
    supportByNode.reserve(localPixels.size());
    for (std::size_t localIndex = 0; localIndex < localPixels.size(); ++localIndex) {
        const PixelId pixelId = localPixels[localIndex];
        NodeId current = ownerByLocalIndex[localIndex];
        while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
            if (nodeInSnapshot[static_cast<std::size_t>(current)] != 0) {
                supportByNode[current].push_back(pixelId);
            }
            const NodeId parentId = tree.getNodeParent(current);
            if (parentId == current) {
                break;
            }
            current = parentId;
        }
    }

    std::vector<LocalDualNodeSnapshot> snapshots;
    snapshots.reserve(snapshotNodeIds.size());
    for (NodeId nodeId : snapshotNodeIds) {
        LocalDualNodeSnapshot snapshot;
        snapshot.nodeId = nodeId;
        snapshot.parentId = tree.getNodeParent(nodeId);
        snapshot.parentInSnapshot =
            snapshot.parentId != InvalidNode &&
            snapshot.parentId != nodeId &&
            tree.isNode(snapshot.parentId) &&
            nodeInSnapshot[static_cast<std::size_t>(snapshot.parentId)] != 0;
        snapshot.altitude = tree.getAltitude(nodeId);
        snapshot.numProperParts = tree.getNumProperParts(nodeId);
        snapshot.numChildren = tree.getNumChildren(nodeId);

        auto directIt = directByNode.find(nodeId);
        if (directIt != directByNode.end()) {
            snapshot.directLocalProperParts = std::move(directIt->second);
        }
        snapshot.directLocalProperParts = sortedUniquePixels(std::move(snapshot.directLocalProperParts));

        auto supportIt = supportByNode.find(nodeId);
        if (supportIt != supportByNode.end()) {
            snapshot.supportLocalPixels = sortedUniquePixels(std::move(supportIt->second));
        }
        for (PixelId pixelId : snapshot.supportLocalPixels) {
            const auto pixelIndex = static_cast<std::size_t>(pixelId);
            const unsigned char role = pixelIndex < roleByPixel.size() ? roleByPixel[pixelIndex] : 0;
            snapshot.supportTouchesModified =
                snapshot.supportTouchesModified ||
                role == localRoleCode(LocalMergePixelRole::Modified);
            snapshot.supportTouchesPreserved =
                snapshot.supportTouchesPreserved ||
                role == localRoleCode(LocalMergePixelRole::Preserved);
            snapshot.supportTouchesBoundary =
                snapshot.supportTouchesBoundary ||
                role == localRoleCode(LocalMergePixelRole::Boundary);
        }

        for (NodeId childId : tree.getChildren(nodeId)) {
            if (childId >= 0 &&
                childId < nodeSlots &&
                nodeInSnapshot[static_cast<std::size_t>(childId)] != 0) {
                snapshot.childrenInSnapshot.push_back(childId);
            }
        }
        std::sort(snapshot.childrenInSnapshot.begin(), snapshot.childrenInSnapshot.end());

        snapshots.push_back(std::move(snapshot));
    }

    std::sort(snapshots.begin(), snapshots.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.altitude != rhs.altitude) {
            return lhs.altitude < rhs.altitude;
        }
        return lhs.nodeId < rhs.nodeId;
    });
    return snapshots;
}

inline std::vector<LocalDualNodeSnapshot> buildLocalDualNodeSnapshotFromCompactLocalRoles(
    const DynamicComponentTree &tree,
    std::span<const PixelId> localPixels,
    const std::unordered_map<PixelId, unsigned char> &roleByPixel) {
    std::unordered_set<NodeId> nodeInSnapshot;
    nodeInSnapshot.reserve(localPixels.size() * 2);
    std::vector<NodeId> snapshotNodeIds;
    snapshotNodeIds.reserve(localPixels.size() * 2);

    auto addNode = [&](NodeId nodeId) {
        if (nodeId != InvalidNode && tree.isNode(nodeId) && tree.isAlive(nodeId)) {
            const auto [_, inserted] = nodeInSnapshot.insert(nodeId);
            if (inserted) {
                snapshotNodeIds.push_back(nodeId);
            }
        }
    };

    std::unordered_map<NodeId, std::vector<PixelId>> directByNode;
    directByNode.reserve(localPixels.size());
    std::vector<NodeId> ownerByLocalIndex;
    ownerByLocalIndex.reserve(localPixels.size());
    for (PixelId pixelId : localPixels) {
        const NodeId ownerId = tree.getSmallestComponent(pixelId);
        ownerByLocalIndex.push_back(ownerId);
        addNode(ownerId);
        if (ownerId != InvalidNode && tree.isAlive(ownerId)) {
            directByNode[ownerId].push_back(pixelId);
            const NodeId parentId = tree.getNodeParent(ownerId);
            if (parentId != ownerId) {
                addNode(parentId);
            }
        }
    }

    std::unordered_map<NodeId, std::vector<PixelId>> supportByNode;
    supportByNode.reserve(localPixels.size());
    for (std::size_t localIndex = 0; localIndex < localPixels.size(); ++localIndex) {
        const PixelId pixelId = localPixels[localIndex];
        NodeId current = ownerByLocalIndex[localIndex];
        while (current != InvalidNode && tree.isNode(current) && tree.isAlive(current)) {
            if (nodeInSnapshot.find(current) != nodeInSnapshot.end()) {
                supportByNode[current].push_back(pixelId);
            }
            const NodeId parentId = tree.getNodeParent(current);
            if (parentId == current) {
                break;
            }
            current = parentId;
        }
    }

    std::vector<LocalDualNodeSnapshot> snapshots;
    snapshots.reserve(snapshotNodeIds.size());
    for (NodeId nodeId : snapshotNodeIds) {
        LocalDualNodeSnapshot snapshot;
        snapshot.nodeId = nodeId;
        snapshot.parentId = tree.getNodeParent(nodeId);
        snapshot.parentInSnapshot =
            snapshot.parentId != InvalidNode &&
            snapshot.parentId != nodeId &&
            tree.isNode(snapshot.parentId) &&
            nodeInSnapshot.find(snapshot.parentId) != nodeInSnapshot.end();
        snapshot.altitude = tree.getAltitude(nodeId);
        snapshot.numProperParts = tree.getNumProperParts(nodeId);
        snapshot.numChildren = tree.getNumChildren(nodeId);

        auto directIt = directByNode.find(nodeId);
        if (directIt != directByNode.end()) {
            snapshot.directLocalProperParts = std::move(directIt->second);
        }
        snapshot.directLocalProperParts = sortedUniquePixels(std::move(snapshot.directLocalProperParts));

        auto supportIt = supportByNode.find(nodeId);
        if (supportIt != supportByNode.end()) {
            snapshot.supportLocalPixels = sortedUniquePixels(std::move(supportIt->second));
        }
        for (PixelId pixelId : snapshot.supportLocalPixels) {
            const auto roleIt = roleByPixel.find(pixelId);
            const unsigned char role =
                roleIt == roleByPixel.end() ? 0 : roleIt->second;
            snapshot.supportTouchesModified =
                snapshot.supportTouchesModified ||
                role == localRoleCode(LocalMergePixelRole::Modified);
            snapshot.supportTouchesPreserved =
                snapshot.supportTouchesPreserved ||
                role == localRoleCode(LocalMergePixelRole::Preserved);
            snapshot.supportTouchesBoundary =
                snapshot.supportTouchesBoundary ||
                role == localRoleCode(LocalMergePixelRole::Boundary);
        }

        for (NodeId childId : tree.getChildren(nodeId)) {
            if (childId != InvalidNode &&
                nodeInSnapshot.find(childId) != nodeInSnapshot.end()) {
                snapshot.childrenInSnapshot.push_back(childId);
            }
        }
        std::sort(snapshot.childrenInSnapshot.begin(), snapshot.childrenInSnapshot.end());

        snapshots.push_back(std::move(snapshot));
    }

    std::sort(snapshots.begin(), snapshots.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.altitude != rhs.altitude) {
            return lhs.altitude < rhs.altitude;
        }
        return lhs.nodeId < rhs.nodeId;
    });
    return snapshots;
}

inline LocalPartitionInput makeLocalPartitionInputFromLiveTrees(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId selectedNode,
    const LocalPartitionInputOptions &options = {}) {
    if (!adjacency) {
        throw std::runtime_error("makeLocalPartitionInputFromLiveTrees requires a valid adjacency relation");
    }

    ImageUInt8Ptr preImage;
    if (options.checkPreReconstruction || options.buildFullUpdatedImage) {
        preImage = primalTree.reconstructionImage();
    }
    if (options.checkPreReconstruction &&
        !dualTree.reconstructionImage()->isEqual(preImage)) {
        throw std::runtime_error("makeLocalPartitionInputFromLiveTrees requires matching primal/dual reconstructions");
    }

    LocalPartitionInput input;
    input.contraction = buildLocalMergeNodeContractionDescriptor(
        primalTree,
        dualTree,
        adjacency,
        selectedNode);

    std::vector<PixelId> preservedPixels;
    for (const auto &support : input.contraction.preservedChildSupports) {
        appendPixels(preservedPixels, std::span<const PixelId>(support));
    }
    preservedPixels = sortedUniquePixels(std::move(preservedPixels));

    std::vector<PixelId> corePixels = input.contraction.modifiedProperPart;
    appendPixels(corePixels, std::span<const PixelId>(preservedPixels));
    corePixels = sortedUniquePixels(std::move(corePixels));

    const int domainSize = primalTree.getNumTotalProperParts();
    if (dualTree.getNumTotalProperParts() != domainSize) {
        throw std::runtime_error("makeLocalPartitionInputFromLiveTrees requires matching tree domains");
    }
    input.domainSize = domainSize;

    if (!options.checkPreReconstruction && !options.buildFullUpdatedImage) {
        std::unordered_map<PixelId, unsigned char> compactRoleByPixel;
        compactRoleByPixel.reserve(corePixels.size() * 2);
        for (PixelId pixelId : input.contraction.modifiedProperPart) {
            if (pixelId < 0 || pixelId >= domainSize) {
                throw std::runtime_error("live local partition input modified region contains an out-of-bounds pixel");
            }
            compactRoleByPixel[pixelId] = localRoleCode(LocalMergePixelRole::Modified);
        }
        for (PixelId pixelId : preservedPixels) {
            if (pixelId < 0 || pixelId >= domainSize) {
                throw std::runtime_error("live local partition input preserved region contains an out-of-bounds pixel");
            }
            if (compactRoleByPixel.find(pixelId) == compactRoleByPixel.end()) {
                compactRoleByPixel[pixelId] =
                    localRoleCode(LocalMergePixelRole::Preserved);
            }
        }

        auto boundaryPixels = collectBoundaryPixelsAndMarkCompactLocalRoles(
            *adjacency,
            std::span<const PixelId>(corePixels),
            domainSize,
            compactRoleByPixel);

        input.localPixels = corePixels;
        appendPixels(input.localPixels, std::span<const PixelId>(boundaryPixels));
        input.localPixels = sortedUniquePixels(std::move(input.localPixels));

        input.updatedLocalPixelValues = buildUpdatedLocalPixelValuesForLocalPartition(
            primalTree,
            input.contraction,
            std::span<const PixelId>(input.localPixels));

        input.pixelMetadata.reserve(input.localPixels.size());
        for (PixelId pixelId : input.localPixels) {
            const auto roleIt = compactRoleByPixel.find(pixelId);
            if (roleIt == compactRoleByPixel.end()) {
                throw std::runtime_error("compact live local partition input has a pixel without metadata");
            }
            LocalPartitionPixelMetadata metadata;
            metadata.pixelId = pixelId;
            metadata.role = localRoleFromCode(roleIt->second);
            metadata.preDualOwner = dualTree.getSmallestComponent(pixelId);
            input.pixelMetadata.push_back(metadata);
        }

        input.preDualNodes = buildLocalDualNodeSnapshotFromCompactLocalRoles(
            dualTree,
            std::span<const PixelId>(input.localPixels),
            compactRoleByPixel);
        return input;
    }

    std::vector<unsigned char> roleByPixel(static_cast<std::size_t>(domainSize), 0);
    for (PixelId pixelId : input.contraction.modifiedProperPart) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("live local partition input modified region contains an out-of-bounds pixel");
        }
        roleByPixel[static_cast<std::size_t>(pixelId)] =
            localRoleCode(LocalMergePixelRole::Modified);
    }
    for (PixelId pixelId : preservedPixels) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("live local partition input preserved region contains an out-of-bounds pixel");
        }
        const auto pixelIndex = static_cast<std::size_t>(pixelId);
        if (roleByPixel[pixelIndex] == 0) {
            roleByPixel[pixelIndex] = localRoleCode(LocalMergePixelRole::Preserved);
        }
    }

    auto boundaryPixels = collectBoundaryPixelsAndMarkLocalRoles(
        *adjacency,
        std::span<const PixelId>(corePixels),
        roleByPixel);

    input.localPixels = corePixels;
    appendPixels(input.localPixels, std::span<const PixelId>(boundaryPixels));
    input.localPixels = sortedUniquePixels(std::move(input.localPixels));

    if (options.buildFullUpdatedImage) {
        input.updatedImage = buildUpdatedImageForLocalPartition(
            primalTree,
            preImage,
            input.contraction,
            std::span<const PixelId>(input.localPixels),
            true);
    } else {
        input.updatedLocalPixelValues = buildUpdatedLocalPixelValuesForLocalPartition(
            primalTree,
            input.contraction,
            std::span<const PixelId>(input.localPixels));
    }

    input.pixelMetadata.reserve(input.localPixels.size());
    for (PixelId pixelId : input.localPixels) {
        LocalPartitionPixelMetadata metadata;
        metadata.pixelId = pixelId;
        metadata.role = localRoleFromCode(roleByPixel[static_cast<std::size_t>(pixelId)]);
        metadata.preDualOwner = dualTree.getSmallestComponent(pixelId);
        input.pixelMetadata.push_back(metadata);
    }

    input.preDualNodes = buildLocalDualNodeSnapshotFromLocalRoles(
        dualTree,
        std::span<const PixelId>(input.localPixels),
        roleByPixel);
    return input;
}

inline std::vector<const LocalDualNodeSnapshot *> buildSnapshotIndexByNodeId(
    std::span<const LocalDualNodeSnapshot> snapshots) {
    NodeId maxNodeId = InvalidNode;
    for (const auto &snapshot : snapshots) {
        if (snapshot.nodeId > maxNodeId) {
            maxNodeId = snapshot.nodeId;
        }
        if (snapshot.parentInSnapshot && snapshot.parentId > maxNodeId) {
            maxNodeId = snapshot.parentId;
        }
    }

    std::vector<const LocalDualNodeSnapshot *> index(
        static_cast<std::size_t>(std::max<NodeId>(maxNodeId + 1, 0)),
        nullptr);
    for (const auto &snapshot : snapshots) {
        if (snapshot.nodeId != InvalidNode &&
            snapshot.nodeId >= 0 &&
            static_cast<std::size_t>(snapshot.nodeId) < index.size()) {
            index[static_cast<std::size_t>(snapshot.nodeId)] = &snapshot;
        }
    }
    return index;
}

inline const LocalDualNodeSnapshot *findIndexedSnapshotByNodeId(
    std::span<const LocalDualNodeSnapshot *const> snapshotByNodeId,
    NodeId nodeId) {
    if (nodeId == InvalidNode ||
        nodeId < 0 ||
        static_cast<std::size_t>(nodeId) >= snapshotByNodeId.size()) {
        return nullptr;
    }
    return snapshotByNodeId[static_cast<std::size_t>(nodeId)];
}

inline NodeId thresholdOwnerInIndexedLocalSnapshot(
    std::span<const LocalDualNodeSnapshot *const> snapshotByNodeId,
    NodeId smallestOwner,
    AltitudeType level,
    bool dualIsMaxtree) {
    const auto *current = findIndexedSnapshotByNodeId(snapshotByNodeId, smallestOwner);
    if (current == nullptr ||
        !pixelPassesDualThreshold(current->altitude, level, dualIsMaxtree)) {
        return InvalidNode;
    }

    NodeId thresholdOwner = current->nodeId;
    while (current->parentInSnapshot) {
        const auto *parent = findIndexedSnapshotByNodeId(snapshotByNodeId, current->parentId);
        if (parent == nullptr ||
            !pixelPassesDualThreshold(parent->altitude, level, dualIsMaxtree)) {
            break;
        }
        thresholdOwner = parent->nodeId;
        current = parent;
    }
    return thresholdOwner;
}

inline LocalThresholdComponent buildLocalThresholdComponentWithoutOracle(
    const LocalPartitionInput &input,
    AltitudeType level,
    std::vector<PixelId> pixels,
    const std::vector<unsigned char> &roleKnown,
    const std::vector<LocalMergePixelRole> &roleByPixel,
    const std::vector<NodeId> &preOwnerByPixel,
    std::span<const LocalDualNodeSnapshot *const> preSnapshotByNodeId) {
    LocalThresholdComponent component;
    component.level = level;
    component.pixels = sortedUniquePixels(std::move(pixels));

    for (PixelId pixelId : component.pixels) {
        const auto index = static_cast<std::size_t>(pixelId);
        if (index >= roleKnown.size() || roleKnown[index] == 0) {
            throw std::runtime_error("non-oracle local partition component contains a pixel without metadata");
        }

        const auto role = roleByPixel[index];
        component.touchesModified = component.touchesModified || role == LocalMergePixelRole::Modified;
        component.touchesPreserved = component.touchesPreserved || role == LocalMergePixelRole::Preserved;
        component.touchesBoundary = component.touchesBoundary || role == LocalMergePixelRole::Boundary;

        const NodeId preOwner = preOwnerByPixel[index];
        component.preDualOwners.push_back(preOwner);
        component.preThresholdOwners.push_back(thresholdOwnerInIndexedLocalSnapshot(
            preSnapshotByNodeId,
            preOwner,
            level,
            input.contraction.dualIsMaxtree));
    }

    component.preDualOwners = sortedUniqueNodes(std::move(component.preDualOwners));
    component.preThresholdOwners = sortedUniqueNodes(std::move(component.preThresholdOwners));
    return component;
}

struct CompactLocalPartitionMetadata {
    std::vector<PixelId> localPixels;
    std::vector<LocalMergePixelRole> roleByLocalIndex;
    std::vector<NodeId> preOwnerByLocalIndex;
};

inline std::size_t findLocalPixelIndex(std::span<const PixelId> localPixels,
                                       PixelId pixelId) {
    const auto it = std::lower_bound(localPixels.begin(), localPixels.end(), pixelId);
    if (it == localPixels.end() || *it != pixelId) {
        return localPixels.size();
    }
    return static_cast<std::size_t>(std::distance(localPixels.begin(), it));
}

inline CompactLocalPartitionMetadata buildCompactLocalPartitionMetadata(
    const LocalPartitionInput &input,
    int domainSize) {
    CompactLocalPartitionMetadata metadata;
    metadata.localPixels = input.localPixels;
    metadata.roleByLocalIndex.assign(
        metadata.localPixels.size(),
        LocalMergePixelRole::Boundary);
    metadata.preOwnerByLocalIndex.assign(metadata.localPixels.size(), InvalidNode);
    std::vector<unsigned char> known(metadata.localPixels.size(), 0);

    for (const auto &pixelMetadata : input.pixelMetadata) {
        if (pixelMetadata.pixelId < 0 || pixelMetadata.pixelId >= domainSize) {
            throw std::runtime_error("compact local partition received an out-of-bounds pixel");
        }
        const auto localIndex = findLocalPixelIndex(
            std::span<const PixelId>(metadata.localPixels),
            pixelMetadata.pixelId);
        if (localIndex == metadata.localPixels.size()) {
            throw std::runtime_error("compact local partition received metadata outside the local pixel set");
        }
        if (known[localIndex] != 0) {
            throw std::runtime_error("compact local partition received duplicated pixel metadata");
        }
        known[localIndex] = 1;
        metadata.roleByLocalIndex[localIndex] = pixelMetadata.role;
        metadata.preOwnerByLocalIndex[localIndex] = pixelMetadata.preDualOwner;
    }

    for (std::size_t index = 0; index < metadata.localPixels.size(); ++index) {
        if (known[index] == 0) {
            throw std::runtime_error("compact local partition requires metadata for every local pixel");
        }
    }
    return metadata;
}

inline LocalThresholdComponent buildCompactLocalThresholdComponentWithoutOracle(
    const LocalPartitionInput &input,
    const CompactLocalPartitionMetadata &metadata,
    AltitudeType level,
    std::vector<PixelId> pixels,
    std::span<const LocalDualNodeSnapshot *const> preSnapshotByNodeId) {
    LocalThresholdComponent component;
    component.level = level;
    component.pixels = sortedUniquePixels(std::move(pixels));

    for (PixelId pixelId : component.pixels) {
        const auto localIndex = findLocalPixelIndex(
            std::span<const PixelId>(metadata.localPixels),
            pixelId);
        if (localIndex == metadata.localPixels.size()) {
            throw std::runtime_error("compact local partition component contains a pixel without metadata");
        }

        const auto role = metadata.roleByLocalIndex[localIndex];
        component.touchesModified = component.touchesModified || role == LocalMergePixelRole::Modified;
        component.touchesPreserved = component.touchesPreserved || role == LocalMergePixelRole::Preserved;
        component.touchesBoundary = component.touchesBoundary || role == LocalMergePixelRole::Boundary;

        const NodeId preOwner = metadata.preOwnerByLocalIndex[localIndex];
        component.preDualOwners.push_back(preOwner);
        component.preThresholdOwners.push_back(thresholdOwnerInIndexedLocalSnapshot(
            preSnapshotByNodeId,
            preOwner,
            level,
            input.contraction.dualIsMaxtree));
    }

    component.preDualOwners = sortedUniqueNodes(std::move(component.preDualOwners));
    component.preThresholdOwners = sortedUniqueNodes(std::move(component.preThresholdOwners));
    return component;
}

inline bool preSnapshotCanGuideLocalMerge(const LocalPartitionInput &input,
                                          const LocalDualNodeSnapshot &snapshot,
                                          AltitudeType level) {
    if (!snapshot.supportTouchesModified ||
        snapshot.supportLocalPixels.empty()) {
        return false;
    }

    const auto snapshotLevel = static_cast<AltitudeType>(snapshot.altitude);
    if (snapshotLevel == level) {
        return true;
    }

    return !pixelPassesDualThreshold(snapshot.altitude,
                                    level,
                                    input.contraction.dualIsMaxtree) &&
           validLocalParentLevel(level,
                                 snapshotLevel,
                                 input.contraction.dualIsMaxtree);
}

inline bool preSnapshotHasProperThresholdChildAtLevel(
    const LocalPartitionInput &input,
    const LocalDualNodeSnapshot &snapshot,
    AltitudeType level) {
    auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
    for (NodeId childId : snapshot.childrenInSnapshot) {
        const auto *child = findPartitionSnapshotByNodeId(
            std::span<const LocalDualNodeSnapshot>(input.preDualNodes),
            childId);
        auto childSupport = child != nullptr
                                ? sortedUniquePixels(child->supportLocalPixels)
                                : std::vector<PixelId>{};
        if (child != nullptr &&
            !childSupport.empty() &&
            !equalPixelSets(std::span<const PixelId>(snapshotSupport),
                            std::span<const PixelId>(childSupport)) &&
            pixelPassesDualThreshold(child->altitude,
                                     level,
                                     input.contraction.dualIsMaxtree)) {
            return true;
        }
    }
    return false;
}

inline bool componentsSharePreThresholdOwnerAtLevel(
    const LocalPartitionInput &input,
    const LocalThresholdComponent &lhs,
    const LocalThresholdComponent &rhs,
    AltitudeType level) {
    auto lhsIt = lhs.preThresholdOwners.begin();
    auto rhsIt = rhs.preThresholdOwners.begin();
    while (lhsIt != lhs.preThresholdOwners.end() &&
           rhsIt != rhs.preThresholdOwners.end()) {
        if (*lhsIt == *rhsIt) {
            const auto *snapshot = findPartitionSnapshotByNodeId(
                std::span<const LocalDualNodeSnapshot>(input.preDualNodes),
                *lhsIt);
            if (snapshot != nullptr &&
                snapshot->altitude == static_cast<int>(level)) {
                return true;
            }
            ++lhsIt;
            ++rhsIt;
            continue;
        }
        if (*lhsIt < *rhsIt) {
            ++lhsIt;
        } else {
            ++rhsIt;
        }
    }
    return false;
}

inline bool mergeFirstPreThresholdOwnerMatchedComponentSet(
    const LocalPartitionInput &input,
    AltitudeType level,
    std::vector<LocalThresholdComponent> &components) {
    for (std::size_t i = 0; i < components.size(); ++i) {
        for (std::size_t j = i + 1; j < components.size(); ++j) {
            if (!componentsSharePreThresholdOwnerAtLevel(
                    input,
                    components[i],
                    components[j],
                    level)) {
                continue;
            }
            if (!components[i].touchesModified &&
                !components[j].touchesModified) {
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

inline void mergePreThresholdOwnerMatchedLocalThresholdComponents(
    const LocalPartitionInput &input,
    AltitudeType level,
    std::vector<LocalThresholdComponent> &components) {
    while (mergeFirstPreThresholdOwnerMatchedComponentSet(input, level, components)) {
    }
}

inline bool mergeFirstPreSnapshotCoveredComponentSet(
    const LocalPartitionInput &input,
    AltitudeType level,
    std::vector<LocalThresholdComponent> &components) {
    for (const auto &snapshot : input.preDualNodes) {
        if (!preSnapshotCanGuideLocalMerge(input, snapshot, level)) {
            continue;
        }

        auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
        std::vector<std::size_t> coveredIndices;
        std::vector<PixelId> mergedPixels;
        bool touchesModified = false;
        bool modifiedComponentAlreadyMatchesSnapshot = false;
        int coveredModifiedComponentCount = 0;

        for (std::size_t index = 0; index < components.size(); ++index) {
            const auto &component = components[index];
            if (!containsAllPixels(std::span<const PixelId>(snapshotSupport),
                                   std::span<const PixelId>(component.pixels))) {
                continue;
            }
            coveredIndices.push_back(index);
            appendPixels(mergedPixels, std::span<const PixelId>(component.pixels));
            touchesModified = touchesModified || component.touchesModified;
            if (component.touchesModified) {
                ++coveredModifiedComponentCount;
            }
            modifiedComponentAlreadyMatchesSnapshot =
                modifiedComponentAlreadyMatchesSnapshot ||
                (component.touchesModified &&
                 equalPixelSets(std::span<const PixelId>(snapshotSupport),
                                std::span<const PixelId>(component.pixels)));
        }

        if (coveredIndices.size() < 2 ||
            !touchesModified ||
            modifiedComponentAlreadyMatchesSnapshot ||
            (snapshot.altitude != static_cast<int>(level) &&
             preSnapshotHasProperThresholdChildAtLevel(input, snapshot, level)) ||
            (snapshot.altitude != static_cast<int>(level) &&
             coveredModifiedComponentCount > 1)) {
            continue;
        }

        mergedPixels = sortedUniquePixels(std::move(mergedPixels));
        if (!equalPixelSets(std::span<const PixelId>(snapshotSupport),
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

inline void mergePreSnapshotCoveredLocalThresholdComponents(
    const LocalPartitionInput &input,
    AltitudeType level,
    std::vector<LocalThresholdComponent> &components) {
    while (mergeFirstPreSnapshotCoveredComponentSet(input, level, components)) {
    }
}

inline std::vector<LocalThresholdComponent> collectLocalThresholdComponentsWithoutOracle(
    const LocalPartitionInput &input,
    AdjacencyRelation &adjacency,
    AltitudeType level,
    const std::vector<unsigned char> &localMask,
    const std::vector<unsigned char> &roleKnown,
    const std::vector<LocalMergePixelRole> &roleByPixel,
    const std::vector<NodeId> &preOwnerByPixel,
    std::span<const LocalDualNodeSnapshot *const> preSnapshotByNodeId) {
    const int domainSize = localPartitionDomainSize(input);
    std::vector<unsigned char> allowed(static_cast<std::size_t>(domainSize), 0);
    std::vector<unsigned char> visited(static_cast<std::size_t>(domainSize), 0);
    std::vector<PixelId> seeds;
    seeds.reserve(input.localPixels.size());

    for (PixelId pixelId : input.localPixels) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("non-oracle local partition received an out-of-bounds local pixel");
        }
        if (pixelPassesDualThreshold(localPartitionUpdatedPixelValue(input, pixelId),
                                     level,
                                     input.contraction.dualIsMaxtree)) {
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

        components.push_back(buildLocalThresholdComponentWithoutOracle(
            input,
            level,
            std::move(pixels),
            roleKnown,
            roleByPixel,
            preOwnerByPixel,
            preSnapshotByNodeId));
    }

    mergePreThresholdOwnerMatchedLocalThresholdComponents(input, level, components);
    mergePreSnapshotCoveredLocalThresholdComponents(input, level, components);
    std::sort(components.begin(), components.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.pixels < rhs.pixels;
    });
    return components;
}

inline std::vector<LocalThresholdComponent> collectCompactLocalThresholdComponentsWithoutOracle(
    const LocalPartitionInput &input,
    AdjacencyRelation &adjacency,
    AltitudeType level,
    const CompactLocalPartitionMetadata &metadata,
    std::span<const LocalDualNodeSnapshot *const> preSnapshotByNodeId) {
    const int domainSize = localPartitionDomainSize(input);
    std::vector<unsigned char> allowed(metadata.localPixels.size(), 0);
    std::vector<unsigned char> visited(metadata.localPixels.size(), 0);
    std::vector<std::size_t> seeds;
    seeds.reserve(metadata.localPixels.size());

    for (std::size_t localIndex = 0; localIndex < metadata.localPixels.size(); ++localIndex) {
        const PixelId pixelId = metadata.localPixels[localIndex];
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error("compact non-oracle local partition received an out-of-bounds local pixel");
        }
        if (pixelPassesDualThreshold(localPartitionUpdatedPixelValue(input, pixelId),
                                     level,
                                     input.contraction.dualIsMaxtree)) {
            allowed[localIndex] = 1;
            seeds.push_back(localIndex);
        }
    }

    std::vector<LocalThresholdComponent> components;
    std::vector<std::size_t> stack;
    stack.reserve(seeds.size());

    for (std::size_t seedIndex : seeds) {
        if (visited[seedIndex] != 0) {
            continue;
        }

        std::vector<PixelId> pixels;
        visited[seedIndex] = 1;
        stack.push_back(seedIndex);

        while (!stack.empty()) {
            const std::size_t currentIndex = stack.back();
            stack.pop_back();
            const PixelId currentPixel = metadata.localPixels[currentIndex];
            pixels.push_back(currentPixel);

            for (PixelId neighbor : adjacency.getNeighborPixels(currentPixel)) {
                if (neighbor < 0 || neighbor >= domainSize) {
                    continue;
                }
                const auto neighborIndex = findLocalPixelIndex(
                    std::span<const PixelId>(metadata.localPixels),
                    neighbor);
                if (neighborIndex == metadata.localPixels.size() ||
                    allowed[neighborIndex] == 0 ||
                    visited[neighborIndex] != 0) {
                    continue;
                }
                visited[neighborIndex] = 1;
                stack.push_back(neighborIndex);
            }
        }

        components.push_back(buildCompactLocalThresholdComponentWithoutOracle(
            input,
            metadata,
            level,
            std::move(pixels),
            preSnapshotByNodeId));
    }

    mergePreThresholdOwnerMatchedLocalThresholdComponents(input, level, components);
    mergePreSnapshotCoveredLocalThresholdComponents(input, level, components);
    std::sort(components.begin(), components.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.pixels < rhs.pixels;
    });
    return components;
}

inline LocalPartitionPlan buildLocalPartitionPlanWithoutOracle(
    const LocalPartitionInput &input,
    AdjacencyRelation &adjacency) {
    if (!localPartitionHasUpdatedTarget(input)) {
        throw std::runtime_error("buildLocalPartitionPlanWithoutOracle requires target values");
    }
    if (input.localPixels.empty()) {
        throw std::runtime_error("buildLocalPartitionPlanWithoutOracle requires non-empty local pixels");
    }

    const int domainSize = localPartitionDomainSize(input);
    const auto preSnapshotByNodeId = buildSnapshotIndexByNodeId(
        std::span<const LocalDualNodeSnapshot>(input.preDualNodes));

    if (input.localPixels.size() <= CompactLocalPartitionPixelLimit) {
        const auto compactMetadata = buildCompactLocalPartitionMetadata(input, domainSize);
        LocalPartitionPlan plan;
        plan.affectedLevels = affectedLevelsForLocalMerge(input.contraction);
        plan.levels.reserve(plan.affectedLevels.size());

        for (AltitudeType level : plan.affectedLevels) {
            LocalLevelPartition partition;
            partition.level = level;
            partition.components = collectCompactLocalThresholdComponentsWithoutOracle(
                input,
                adjacency,
                level,
                compactMetadata,
                std::span<const LocalDualNodeSnapshot *const>(preSnapshotByNodeId));

            for (const auto &component : partition.components) {
                if (component.isRecomposedComponent()) {
                    plan.recomposedComponents.push_back(component);
                }
            }
            plan.levels.push_back(std::move(partition));
        }

        return plan;
    }

    const auto localMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(input.localPixels),
        "non-oracle local partition");

    std::vector<unsigned char> roleKnown(static_cast<std::size_t>(domainSize), 0);
    std::vector<LocalMergePixelRole> roleByPixel(
        static_cast<std::size_t>(domainSize),
        LocalMergePixelRole::Boundary);
    std::vector<NodeId> preOwnerByPixel(static_cast<std::size_t>(domainSize), InvalidNode);

    for (const auto &metadata : input.pixelMetadata) {
        if (metadata.pixelId < 0 || metadata.pixelId >= domainSize) {
            throw std::runtime_error("buildLocalPartitionPlanWithoutOracle received an out-of-bounds pixel");
        }
        const auto index = static_cast<std::size_t>(metadata.pixelId);
        if (roleKnown[index] != 0) {
            throw std::runtime_error("buildLocalPartitionPlanWithoutOracle received duplicated pixel metadata");
        }
        roleKnown[index] = 1;
        roleByPixel[index] = metadata.role;
        preOwnerByPixel[index] = metadata.preDualOwner;
    }

    for (PixelId pixelId : input.localPixels) {
        const auto index = static_cast<std::size_t>(pixelId);
        if (roleKnown[index] == 0) {
            throw std::runtime_error("buildLocalPartitionPlanWithoutOracle requires metadata for every local pixel");
        }
    }

    LocalPartitionPlan plan;
    plan.affectedLevels = affectedLevelsForLocalMerge(input.contraction);
    plan.levels.reserve(plan.affectedLevels.size());

    for (AltitudeType level : plan.affectedLevels) {
        LocalLevelPartition partition;
        partition.level = level;
        partition.components = collectLocalThresholdComponentsWithoutOracle(
            input,
            adjacency,
            level,
            localMask,
            roleKnown,
            roleByPixel,
            preOwnerByPixel,
            std::span<const LocalDualNodeSnapshot *const>(preSnapshotByNodeId));

        for (const auto &component : partition.components) {
            if (component.isRecomposedComponent()) {
                plan.recomposedComponents.push_back(component);
            }
        }
        plan.levels.push_back(std::move(partition));
    }

    return plan;
}

inline bool componentMatchesPreSnapshotAtSameLevel(const LocalPartitionInput &input,
                                                   const LocalThresholdComponent &component) {
    for (const auto &snapshot : input.preDualNodes) {
        if (snapshot.altitude != static_cast<int>(component.level)) {
            continue;
        }
        auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
        if (equalPixelSets(std::span<const PixelId>(snapshotSupport),
                           std::span<const PixelId>(component.pixels))) {
            return true;
        }
    }
    return false;
}

inline std::vector<PixelId> collectDirectPixelsAtLevelWithoutOracle(
    const LocalPartitionInput &input,
    const LocalThresholdComponent &component) {
    std::vector<PixelId> directPixels;
    directPixels.reserve(component.pixels.size());
    for (PixelId pixelId : component.pixels) {
        if (localPartitionUpdatedPixelValue(input, pixelId) == component.level) {
            directPixels.push_back(pixelId);
        }
    }
    return sortedUniquePixels(std::move(directPixels));
}

inline bool localSupportTouchesRole(const LocalPartitionInput &input,
                                    std::span<const PixelId> support,
                                    LocalMergePixelRole role) {
    for (PixelId pixelId : support) {
        for (const auto &metadata : input.pixelMetadata) {
            if (metadata.pixelId == pixelId && metadata.role == role) {
                return true;
            }
        }
    }
    return false;
}

inline bool preSnapshotIsStableChildCandidate(const LocalPartitionInput &input,
                                              const LocalRecompositionNodePlan &node,
                                              const LocalDualNodeSnapshot &snapshot) {
    if (snapshot.nodeId == InvalidNode ||
        snapshot.supportLocalPixels.empty() ||
        snapshot.supportTouchesModified) {
        return false;
    }

    const bool touchesPreserved = localSupportTouchesRole(
        input,
        std::span<const PixelId>(snapshot.supportLocalPixels),
        LocalMergePixelRole::Preserved);
    const bool visibleBoundaryChild =
        snapshot.supportTouchesBoundary &&
        !snapshot.directLocalProperParts.empty();
    if (!touchesPreserved && !visibleBoundaryChild) {
        return false;
    }

    auto snapshotSupport = sortedUniquePixels(snapshot.supportLocalPixels);
    if (equalPixelSets(std::span<const PixelId>(snapshotSupport),
                       std::span<const PixelId>(node.supportPixels))) {
        return false;
    }
    if (!containsAllPixels(std::span<const PixelId>(node.supportPixels),
                           std::span<const PixelId>(snapshotSupport))) {
        return false;
    }
    return validLocalParentLevel(static_cast<AltitudeType>(snapshot.altitude),
                                 node.level,
                                 input.contraction.dualIsMaxtree);
}

inline std::vector<NodeId> collectStableChildPreNodesWithoutOracle(
    const LocalPartitionInput &input,
    const LocalRecompositionNodePlan &node) {
    std::vector<const LocalDualNodeSnapshot *> candidates;
    for (const auto &snapshot : input.preDualNodes) {
        if (preSnapshotIsStableChildCandidate(input, node, snapshot)) {
            candidates.push_back(&snapshot);
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const auto *lhs, const auto *rhs) {
        if (lhs->supportLocalPixels.size() != rhs->supportLocalPixels.size()) {
            return lhs->supportLocalPixels.size() > rhs->supportLocalPixels.size();
        }
        if (lhs->altitude != rhs->altitude) {
            return lhs->altitude > rhs->altitude;
        }
        return lhs->nodeId < rhs->nodeId;
    });

    std::vector<const LocalDualNodeSnapshot *> selected;
    for (const auto *candidate : candidates) {
        auto candidateSupport = sortedUniquePixels(candidate->supportLocalPixels);
        bool coveredBySelected = false;
        for (const auto *selectedSnapshot : selected) {
            auto selectedSupport = sortedUniquePixels(selectedSnapshot->supportLocalPixels);
            if (containsAllPixels(std::span<const PixelId>(selectedSupport),
                                  std::span<const PixelId>(candidateSupport))) {
                coveredBySelected = true;
                break;
            }
        }
        if (!coveredBySelected) {
            selected.push_back(candidate);
        }
    }

    std::vector<NodeId> nodeIds;
    nodeIds.reserve(selected.size());
    for (const auto *snapshot : selected) {
        nodeIds.push_back(snapshot->nodeId);
    }
    return sortedUniqueNodes(std::move(nodeIds));
}

inline bool canBeNonOracleParentCandidate(const LocalRecompositionNodePlan &child,
                                          const LocalRecompositionNodePlan &candidate,
                                          bool dualIsMaxtree) {
    if (!validLocalParentLevel(child.level, candidate.level, dualIsMaxtree)) {
        return false;
    }
    if (!containsAllPixels(std::span<const PixelId>(candidate.supportPixels),
                           std::span<const PixelId>(child.supportPixels))) {
        return false;
    }
    if (candidate.supportPixels.size() > child.supportPixels.size()) {
        return true;
    }
    return candidate.supportPixels.size() == child.supportPixels.size() &&
           candidate.directProperParts.empty();
}

inline void assignLocalRecompositionParentIndicesWithoutOracle(
    LocalRecompositionPlan &plan,
    bool dualIsMaxtree) {
    for (auto &node : plan.nodes) {
        node.parentPlanIndex = -1;
        node.childPlanIndices.clear();
    }

    for (std::size_t childIndex = 0; childIndex < plan.nodes.size(); ++childIndex) {
        const auto &child = plan.nodes[childIndex];
        int bestParentIndex = -1;
        std::size_t bestParentSupportSize = 0;
        int bestParentLevelDistance = 0;

        for (std::size_t parentIndex = 0; parentIndex < plan.nodes.size(); ++parentIndex) {
            if (parentIndex == childIndex) {
                continue;
            }
            const auto &candidate = plan.nodes[parentIndex];
            if (!canBeNonOracleParentCandidate(child, candidate, dualIsMaxtree)) {
                continue;
            }

            const int levelDistance =
                std::abs(static_cast<int>(candidate.level) - static_cast<int>(child.level));
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

inline LocalRecompositionPlan buildLocalRecompositionPlanWithoutOracle(
    const LocalPartitionInput &input,
    const LocalPartitionPlan &partitionPlan) {
    if (!localPartitionHasUpdatedTarget(input)) {
        throw std::runtime_error("buildLocalRecompositionPlanWithoutOracle requires target values");
    }

    LocalRecompositionPlan plan;
    plan.affectedLevels = partitionPlan.affectedLevels;

    for (const auto &partition : partitionPlan.levels) {
        for (const auto &component : partition.components) {
            if (!component.isRecomposedComponent()) {
                continue;
            }

            auto directPixels = collectDirectPixelsAtLevelWithoutOracle(input, component);
            const bool matchesSameLevelPreSnapshot =
                componentMatchesPreSnapshotAtSameLevel(input, component);
            if (directPixels.empty() && !matchesSameLevelPreSnapshot) {
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
            node.stableChildPreNodes = collectStableChildPreNodesWithoutOracle(input, node);
            plan.nodes.push_back(std::move(node));
        }
    }

    std::sort(plan.nodes.begin(), plan.nodes.end(), [](const auto &lhs, const auto &rhs) {
        if (lhs.level != rhs.level) {
            return lhs.level < rhs.level;
        }
        return lhs.supportPixels < rhs.supportPixels;
    });
    assignLocalRecompositionParentIndicesWithoutOracle(
        plan,
        input.contraction.dualIsMaxtree);
    return plan;
}

} // namespace morphotree::details::non_increasing
