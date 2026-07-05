#pragma once

#include <algorithm>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../AdjacencyRelation.hpp"
#include "../Common.hpp"
#include "../DynamicComponentTree.hpp"
#include "../debug/details/RegionLevelChange.hpp"
#include "../debug/details/RegionMaterialization.hpp"

namespace morphotree::details::non_increasing {

enum class LocalMergePixelRole {
    Modified,
    Preserved,
    Boundary,
};

struct LocalMergeNodeContractionDescriptor {
    NodeId selectedNode = InvalidNode;
    NodeId parentNode = InvalidNode;
    bool primalIsMaxtree = false;
    bool dualIsMaxtree = false;
    AltitudeType a = 0;
    AltitudeType b = 0;
    ::debug::region::RegionUpdateDirection direction = ::debug::region::RegionUpdateDirection::Unchanged;
    std::vector<PixelId> modifiedProperPart;
    std::vector<std::vector<PixelId>> modifiedComponents;
    std::vector<NodeId> preservedChildren;
    std::vector<std::vector<PixelId>> preservedChildSupports;
};

struct LocalMergePixelOwner {
    PixelId pixelId = InvalidNode;
    LocalMergePixelRole role = LocalMergePixelRole::Boundary;
    NodeId preDualOwner = InvalidNode;
    NodeId oracleDualOwner = InvalidNode;
    int preDualOwnerLevel = 0;
    int oracleDualOwnerLevel = 0;
};

struct LocalDualNodeSnapshot {
    NodeId nodeId = InvalidNode;
    NodeId parentId = InvalidNode;
    bool parentInSnapshot = false;
    int altitude = 0;
    int numProperParts = 0;
    int numChildren = 0;
    std::vector<PixelId> directLocalProperParts;
    std::vector<PixelId> supportLocalPixels;
    std::vector<NodeId> childrenInSnapshot;
    bool supportTouchesModified = false;
    bool supportTouchesPreserved = false;
    bool supportTouchesBoundary = false;
};

struct LocalMergeOracle {
    LocalMergeNodeContractionDescriptor contraction;
    ImageUInt8Ptr preImage;
    ImageUInt8Ptr updatedImage;
    int domainSize = 0;
    std::vector<PixelId> preservedPixels;
    std::vector<PixelId> boundaryPixels;
    std::vector<PixelId> localPixels;
    std::vector<uint8_t> updatedLocalPixelValues;
    std::vector<LocalMergePixelOwner> pixelOwners;
    std::vector<LocalDualNodeSnapshot> preDualNodes;
    std::vector<LocalDualNodeSnapshot> oracleDualNodes;
    bool preDualMatchesPrimal = false;
    bool oracleDualMatchesUpdatedImage = false;
};

inline std::vector<PixelId> sortedUniquePixels(std::vector<PixelId> pixels) {
    std::sort(pixels.begin(), pixels.end());
    pixels.erase(std::unique(pixels.begin(), pixels.end()), pixels.end());
    return pixels;
}

inline std::size_t findSortedPixelIndex(std::span<const PixelId> pixels,
                                        PixelId pixelId) {
    const auto it = std::lower_bound(pixels.begin(), pixels.end(), pixelId);
    if (it == pixels.end() || *it != pixelId) {
        return pixels.size();
    }
    return static_cast<std::size_t>(std::distance(pixels.begin(), it));
}

inline int localMergeOracleDomainSize(const LocalMergeOracle &oracle) {
    if (oracle.updatedImage) {
        return oracle.updatedImage->getSize();
    }
    return oracle.domainSize;
}

inline bool localMergeOracleHasSparseUpdatedTarget(const LocalMergeOracle &oracle) {
    return oracle.domainSize > 0 &&
           oracle.updatedLocalPixelValues.size() == oracle.localPixels.size();
}

inline bool localMergeOracleHasUpdatedTarget(const LocalMergeOracle &oracle) {
    return oracle.updatedImage != nullptr ||
           localMergeOracleHasSparseUpdatedTarget(oracle);
}

inline uint8_t localMergeOracleUpdatedPixelValue(const LocalMergeOracle &oracle,
                                                 PixelId pixelId) {
    if (oracle.updatedImage) {
        return (*oracle.updatedImage)[pixelId];
    }
    const auto localIndex = findSortedPixelIndex(
        std::span<const PixelId>(oracle.localPixels),
        pixelId);
    if (localIndex == oracle.localPixels.size() ||
        localIndex >= oracle.updatedLocalPixelValues.size()) {
        throw std::runtime_error("local merge oracle has no sparse target value for the requested pixel");
    }
    return oracle.updatedLocalPixelValues[localIndex];
}

inline void appendPixels(std::vector<PixelId> &target, std::span<const PixelId> source) {
    target.insert(target.end(), source.begin(), source.end());
}

inline std::vector<PixelId> collectDirectProperParts(const DynamicComponentTree &tree,
                                                     NodeId nodeId) {
    std::vector<PixelId> pixels;
    if (nodeId == InvalidNode || !tree.isAlive(nodeId)) {
        return pixels;
    }
    pixels.reserve(static_cast<std::size_t>(tree.getNumProperParts(nodeId)));
    for (PixelId pixelId : tree.getProperParts(nodeId)) {
        pixels.push_back(pixelId);
    }
    return sortedUniquePixels(std::move(pixels));
}

inline std::vector<PixelId> collectNodeSupportPixels(const DynamicComponentTree &tree,
                                                     NodeId nodeId) {
    std::vector<PixelId> pixels;
    if (nodeId == InvalidNode || !tree.isAlive(nodeId)) {
        return pixels;
    }
    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(pixelId);
        }
    }
    return sortedUniquePixels(std::move(pixels));
}

inline std::vector<unsigned char> makePixelMask(int domainSize,
                                                std::span<const PixelId> pixels,
                                                const char *context) {
    std::vector<unsigned char> mask(static_cast<std::size_t>(domainSize), 0);
    for (PixelId pixelId : pixels) {
        if (pixelId < 0 || pixelId >= domainSize) {
            throw std::runtime_error(std::string(context) + " received an out-of-bounds pixel");
        }
        mask[static_cast<std::size_t>(pixelId)] = 1;
    }
    return mask;
}

inline std::vector<std::vector<PixelId>> sortedConnectedComponents(AdjacencyRelation &adjacency,
                                                                   std::span<const PixelId> pixels,
                                                                   int domainSize) {
    auto components = ::debug::region::collectLocalConnectedComponents(adjacency, pixels, domainSize);
    for (auto &component : components) {
        component = sortedUniquePixels(std::move(component));
    }
    std::sort(components.begin(), components.end());
    return components;
}

inline LocalMergeNodeContractionDescriptor buildLocalMergeNodeContractionDescriptor(
    const DynamicComponentTree &primalTree,
    const DynamicComponentTree &dualTree,
    const AdjacencyRelationPtr &adjacency,
    NodeId selectedNode) {
    if (!adjacency) {
        throw std::runtime_error("buildLocalMergeNodeContractionDescriptor requires a valid adjacency relation");
    }
    if (selectedNode == InvalidNode ||
        selectedNode == primalTree.getRoot() ||
        !primalTree.isNode(selectedNode) ||
        !primalTree.isAlive(selectedNode)) {
        throw std::runtime_error("buildLocalMergeNodeContractionDescriptor requires a live non-root selected node");
    }

    const NodeId parentNode = primalTree.getNodeParent(selectedNode);
    if (parentNode == InvalidNode || parentNode == selectedNode || !primalTree.isAlive(parentNode)) {
        throw std::runtime_error("buildLocalMergeNodeContractionDescriptor cannot locate a live parent");
    }

    LocalMergeNodeContractionDescriptor descriptor;
    descriptor.selectedNode = selectedNode;
    descriptor.parentNode = parentNode;
    descriptor.primalIsMaxtree = primalTree.isMaxtree();
    descriptor.dualIsMaxtree = dualTree.isMaxtree();
    descriptor.a = static_cast<AltitudeType>(primalTree.getAltitude(selectedNode));
    descriptor.b = static_cast<AltitudeType>(primalTree.getAltitude(parentNode));
    descriptor.direction = ::debug::region::directionFromLevels(descriptor.a, descriptor.b);
    descriptor.modifiedProperPart = collectDirectProperParts(primalTree, selectedNode);
    descriptor.modifiedComponents = sortedConnectedComponents(
        *adjacency,
        std::span<const PixelId>(descriptor.modifiedProperPart),
        primalTree.getNumTotalProperParts());

    for (NodeId childId : primalTree.getChildren(selectedNode)) {
        descriptor.preservedChildren.push_back(childId);
        descriptor.preservedChildSupports.push_back(collectNodeSupportPixels(primalTree, childId));
    }

    return descriptor;
}

inline ImageUInt8Ptr applyLocalMergeDescriptorToImage(
    const ImageUInt8Ptr &input,
    const LocalMergeNodeContractionDescriptor &descriptor) {
    if (!input) {
        throw std::runtime_error("applyLocalMergeDescriptorToImage requires a valid input image");
    }
    auto output = input->clone();
    for (PixelId pixelId : descriptor.modifiedProperPart) {
        if (pixelId < 0 || pixelId >= output->getSize()) {
            throw std::runtime_error("applyLocalMergeDescriptorToImage received an out-of-bounds pixel");
        }
        (*output)[pixelId] = static_cast<uint8_t>(descriptor.b);
    }
    return output;
}

inline std::vector<PixelId> collectBoundaryPixels(AdjacencyRelation &adjacency,
                                                  std::span<const PixelId> corePixels,
                                                  const std::vector<unsigned char> &coreMask) {
    std::vector<unsigned char> boundaryMask(coreMask.size(), 0);
    std::vector<PixelId> boundary;
    for (PixelId pixelId : corePixels) {
        for (PixelId neighbor : adjacency.getNeighborPixels(pixelId)) {
            if (neighbor < 0 || static_cast<std::size_t>(neighbor) >= coreMask.size()) {
                continue;
            }
            const auto neighborIndex = static_cast<std::size_t>(neighbor);
            if (coreMask[neighborIndex] != 0 || boundaryMask[neighborIndex] != 0) {
                continue;
            }
            boundaryMask[neighborIndex] = 1;
            boundary.push_back(neighbor);
        }
    }
    return sortedUniquePixels(std::move(boundary));
}

inline LocalMergePixelRole classifyLocalPixel(PixelId pixelId,
                                              const std::vector<unsigned char> &modifiedMask,
                                              const std::vector<unsigned char> &preservedMask) {
    const auto index = static_cast<std::size_t>(pixelId);
    if (modifiedMask[index] != 0) {
        return LocalMergePixelRole::Modified;
    }
    if (preservedMask[index] != 0) {
        return LocalMergePixelRole::Preserved;
    }
    return LocalMergePixelRole::Boundary;
}

inline std::vector<LocalDualNodeSnapshot> buildLocalDualNodeSnapshot(
    const DynamicComponentTree &tree,
    std::span<const PixelId> localPixels,
    const std::vector<unsigned char> &modifiedMask,
    const std::vector<unsigned char> &preservedMask,
    const std::vector<unsigned char> &boundaryMask) {
    const int nodeSlots = tree.getNumInternalNodeSlots();
    std::vector<unsigned char> nodeInSnapshot(static_cast<std::size_t>(nodeSlots), 0);

    auto addNode = [&](NodeId nodeId) {
        if (nodeId != InvalidNode && tree.isNode(nodeId) && tree.isAlive(nodeId)) {
            nodeInSnapshot[static_cast<std::size_t>(nodeId)] = 1;
        }
    };

    for (PixelId pixelId : localPixels) {
        const NodeId ownerId = tree.getSmallestComponent(pixelId);
        addNode(ownerId);
        if (ownerId != InvalidNode && tree.isAlive(ownerId)) {
            const NodeId parentId = tree.getNodeParent(ownerId);
            if (parentId != ownerId) {
                addNode(parentId);
            }
        }
    }

    std::unordered_map<NodeId, std::vector<PixelId>> supportByNode;
    supportByNode.reserve(localPixels.size());
    for (PixelId pixelId : localPixels) {
        NodeId current = tree.getSmallestComponent(pixelId);
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
    for (NodeId nodeId = 0; nodeId < nodeSlots; ++nodeId) {
        if (nodeInSnapshot[static_cast<std::size_t>(nodeId)] == 0) {
            continue;
        }

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

        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            const auto pixelIndex = static_cast<std::size_t>(pixelId);
            if (pixelIndex < modifiedMask.size() &&
                (modifiedMask[pixelIndex] != 0 ||
                 preservedMask[pixelIndex] != 0 ||
                 boundaryMask[pixelIndex] != 0)) {
                snapshot.directLocalProperParts.push_back(pixelId);
            }
        }
        snapshot.directLocalProperParts = sortedUniquePixels(std::move(snapshot.directLocalProperParts));

        auto supportIt = supportByNode.find(nodeId);
        if (supportIt != supportByNode.end()) {
            snapshot.supportLocalPixels = sortedUniquePixels(std::move(supportIt->second));
        }
        for (PixelId pixelId : snapshot.supportLocalPixels) {
            const auto pixelIndex = static_cast<std::size_t>(pixelId);
            snapshot.supportTouchesModified = snapshot.supportTouchesModified || modifiedMask[pixelIndex] != 0;
            snapshot.supportTouchesPreserved = snapshot.supportTouchesPreserved || preservedMask[pixelIndex] != 0;
            snapshot.supportTouchesBoundary = snapshot.supportTouchesBoundary || boundaryMask[pixelIndex] != 0;
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

inline LocalMergeOracle buildLocalMergeOracle(const DynamicComponentTree &primalTree,
                                              const DynamicComponentTree &dualTree,
                                              const AdjacencyRelationPtr &adjacency,
                                              NodeId selectedNode) {
    if (!adjacency) {
        throw std::runtime_error("buildLocalMergeOracle requires a valid adjacency relation");
    }

    LocalMergeOracle oracle;
    oracle.contraction = buildLocalMergeNodeContractionDescriptor(primalTree, dualTree, adjacency, selectedNode);
    oracle.preImage = primalTree.reconstructionImage();
    oracle.preDualMatchesPrimal = dualTree.reconstructionImage()->isEqual(oracle.preImage);
    oracle.updatedImage = applyLocalMergeDescriptorToImage(oracle.preImage, oracle.contraction);
    oracle.domainSize = oracle.updatedImage->getSize();

    DynamicComponentTree rebuiltDual(oracle.updatedImage->clone(), dualTree.isMaxtree(), adjacency);
    oracle.oracleDualMatchesUpdatedImage = rebuiltDual.reconstructionImage()->isEqual(oracle.updatedImage);

    for (const auto &support : oracle.contraction.preservedChildSupports) {
        appendPixels(oracle.preservedPixels, std::span<const PixelId>(support));
    }
    oracle.preservedPixels = sortedUniquePixels(std::move(oracle.preservedPixels));

    std::vector<PixelId> corePixels = oracle.contraction.modifiedProperPart;
    appendPixels(corePixels, std::span<const PixelId>(oracle.preservedPixels));
    corePixels = sortedUniquePixels(std::move(corePixels));

    const int domainSize = primalTree.getNumTotalProperParts();
    const auto modifiedMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(oracle.contraction.modifiedProperPart),
        "modified region");
    const auto preservedMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(oracle.preservedPixels),
        "preserved region");
    const auto coreMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(corePixels),
        "local core");

    oracle.boundaryPixels = collectBoundaryPixels(*adjacency, std::span<const PixelId>(corePixels), coreMask);
    auto boundaryMask = makePixelMask(
        domainSize,
        std::span<const PixelId>(oracle.boundaryPixels),
        "boundary region");

    oracle.localPixels = corePixels;
    appendPixels(oracle.localPixels, std::span<const PixelId>(oracle.boundaryPixels));
    oracle.localPixels = sortedUniquePixels(std::move(oracle.localPixels));

    oracle.pixelOwners.reserve(oracle.localPixels.size());
    for (PixelId pixelId : oracle.localPixels) {
        const NodeId preOwner = dualTree.getSmallestComponent(pixelId);
        const NodeId oracleOwner = rebuiltDual.getSmallestComponent(pixelId);
        LocalMergePixelOwner owner;
        owner.pixelId = pixelId;
        owner.role = classifyLocalPixel(pixelId, modifiedMask, preservedMask);
        owner.preDualOwner = preOwner;
        owner.oracleDualOwner = oracleOwner;
        owner.preDualOwnerLevel =
            preOwner == InvalidNode || !dualTree.isAlive(preOwner) ? -1 : dualTree.getAltitude(preOwner);
        owner.oracleDualOwnerLevel =
            oracleOwner == InvalidNode || !rebuiltDual.isAlive(oracleOwner) ? -1 : rebuiltDual.getAltitude(oracleOwner);
        oracle.pixelOwners.push_back(owner);
    }

    oracle.preDualNodes = buildLocalDualNodeSnapshot(
        dualTree,
        std::span<const PixelId>(oracle.localPixels),
        modifiedMask,
        preservedMask,
        boundaryMask);
    oracle.oracleDualNodes = buildLocalDualNodeSnapshot(
        rebuiltDual,
        std::span<const PixelId>(oracle.localPixels),
        modifiedMask,
        preservedMask,
        boundaryMask);

    return oracle;
}

inline const char *localMergePixelRoleName(LocalMergePixelRole role) {
    switch (role) {
    case LocalMergePixelRole::Modified:
        return "modified";
    case LocalMergePixelRole::Preserved:
        return "preserved";
    case LocalMergePixelRole::Boundary:
        return "boundary";
    }
    return "unknown";
}

inline void appendPixelVector(std::ostringstream &out, std::span<const PixelId> pixels) {
    out << "[";
    for (std::size_t i = 0; i < pixels.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        out << pixels[i];
    }
    out << "]";
}

inline void appendNodeVector(std::ostringstream &out, std::span<const NodeId> nodes) {
    out << "[";
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        out << nodes[i];
    }
    out << "]";
}

inline std::string debugString(const LocalMergeOracle &oracle) {
    std::ostringstream out;
    out << "LocalMergeOracle\n";
    out << "selected=" << oracle.contraction.selectedNode
        << " parent=" << oracle.contraction.parentNode
        << " primal_is_max=" << (oracle.contraction.primalIsMaxtree ? "true" : "false")
        << " dual_is_max=" << (oracle.contraction.dualIsMaxtree ? "true" : "false")
        << " a=" << static_cast<int>(oracle.contraction.a)
        << " b=" << static_cast<int>(oracle.contraction.b)
        << "\n";
    out << "modified=";
    appendPixelVector(out, std::span<const PixelId>(oracle.contraction.modifiedProperPart));
    out << "\n";
    out << "components=[";
    for (std::size_t i = 0; i < oracle.contraction.modifiedComponents.size(); ++i) {
        if (i != 0) {
            out << ",";
        }
        appendPixelVector(out, std::span<const PixelId>(oracle.contraction.modifiedComponents[i]));
    }
    out << "]\n";
    out << "preserved_children=";
    appendNodeVector(out, std::span<const NodeId>(oracle.contraction.preservedChildren));
    out << "\n";
    out << "preserved_pixels=";
    appendPixelVector(out, std::span<const PixelId>(oracle.preservedPixels));
    out << "\n";
    out << "boundary_pixels=";
    appendPixelVector(out, std::span<const PixelId>(oracle.boundaryPixels));
    out << "\n";
    out << "local_pixels=";
    appendPixelVector(out, std::span<const PixelId>(oracle.localPixels));
    out << "\n";
    out << "pre_dual_matches_primal=" << (oracle.preDualMatchesPrimal ? "true" : "false")
        << " oracle_dual_matches_updated=" << (oracle.oracleDualMatchesUpdatedImage ? "true" : "false")
        << "\n";

    out << "pixel_owners\n";
    for (const auto &owner : oracle.pixelOwners) {
        out << "  pixel=" << owner.pixelId
            << " role=" << localMergePixelRoleName(owner.role)
            << " pre_owner=" << owner.preDualOwner
            << " pre_level=" << owner.preDualOwnerLevel
            << " oracle_owner=" << owner.oracleDualOwner
            << " oracle_level=" << owner.oracleDualOwnerLevel
            << "\n";
    }

    auto appendSnapshot = [&out](const char *name, const std::vector<LocalDualNodeSnapshot> &nodes) {
        out << name << "\n";
        for (const auto &node : nodes) {
            out << "  node=" << node.nodeId
                << " parent=" << node.parentId
                << " parent_in_snapshot=" << (node.parentInSnapshot ? "true" : "false")
                << " level=" << node.altitude
                << " proper_parts=" << node.numProperParts
                << " children=" << node.numChildren
                << " touches_modified=" << (node.supportTouchesModified ? "true" : "false")
                << " touches_preserved=" << (node.supportTouchesPreserved ? "true" : "false")
                << " touches_boundary=" << (node.supportTouchesBoundary ? "true" : "false")
                << " direct_local=";
            appendPixelVector(out, std::span<const PixelId>(node.directLocalProperParts));
            out << " support_local=";
            appendPixelVector(out, std::span<const PixelId>(node.supportLocalPixels));
            out << " children_in_snapshot=";
            appendNodeVector(out, std::span<const NodeId>(node.childrenInSnapshot));
            out << "\n";
        }
    };

    appendSnapshot("pre_dual_nodes", oracle.preDualNodes);
    appendSnapshot("oracle_dual_nodes", oracle.oracleDualNodes);
    return out.str();
}

} // namespace morphotree::details::non_increasing
