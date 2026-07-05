#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"

namespace region_debug = debug::region;

namespace {

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

ImageUInt8Ptr make_demo_image() {
    auto image = ImageUInt8::create(4, 4);
    const std::vector<uint8_t> raw = {
        4, 4, 2, 1,
        4, 3, 2, 1,
        5, 5, 2, 0,
        5, 6, 6, 0,
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

NodeId firstNonRootNode(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isRoot(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

std::vector<PixelId> sortedPixels(std::vector<PixelId> pixels) {
    std::sort(pixels.begin(), pixels.end());
    return pixels;
}

std::vector<PixelId> collectSubtreePixels(const DynamicComponentTree &tree, NodeId subtreeRoot) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(pixelId);
        }
    }
    return sortedPixels(std::move(pixels));
}

bool matchesAnySubtreeFootprint(const DynamicComponentTree &tree, std::vector<PixelId> pixels) {
    pixels = sortedPixels(std::move(pixels));
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (collectSubtreePixels(tree, nodeId) == pixels) {
            return true;
        }
    }
    return false;
}

template<typename PixelType>
ImageUInt8Ptr applyRegionLevelChangeToImage(const ImageUInt8Ptr &input,
                                            const region_debug::RegionLevelChange<PixelType> &change) {
    auto output = input->clone();
    for (PixelId pixelId : change.pixelsC) {
        require(pixelId >= 0 && pixelId < output->getSize(), "region pixel must be in image bounds");
        (*output)[pixelId] = static_cast<uint8_t>(change.b);
    }
    return output;
}

struct RebuildOracle {
    ImageUInt8Ptr updatedImage;
    DynamicComponentTree maxTree;
    DynamicComponentTree minTree;

    RebuildOracle(ImageUInt8Ptr image, const AdjacencyRelationPtr &adj)
        : updatedImage(std::move(image)),
          maxTree(updatedImage->clone(), true, adj),
          minTree(updatedImage->clone(), false, adj) {}
};

template<typename PixelType>
RebuildOracle rebuildTreesAfterChange(const ImageUInt8Ptr &input,
                                      const AdjacencyRelationPtr &adj,
                                      const region_debug::RegionLevelChange<PixelType> &change) {
    return RebuildOracle(applyRegionLevelChangeToImage(input, change), adj);
}

void requireSameReconstruction(const RebuildOracle &oracle) {
    const auto maxReconstruction = oracle.maxTree.reconstructionImage();
    const auto minReconstruction = oracle.minTree.reconstructionImage();
    require(maxReconstruction->isEqual(oracle.updatedImage),
            "rebuilt max-tree reconstruction must match the globally updated image");
    require(minReconstruction->isEqual(oracle.updatedImage),
            "rebuilt min-tree reconstruction must match the globally updated image");
    require(maxReconstruction->isEqual(minReconstruction),
            "rebuilt max-tree and min-tree reconstructions must match");
}

template<typename PixelType>
void requireOnlyRegionChanged(const ImageUInt8Ptr &input,
                              const ImageUInt8Ptr &updated,
                              const region_debug::RegionLevelChange<PixelType> &change) {
    std::vector<bool> inRegion(static_cast<size_t>(input->getSize()), false);
    for (PixelId pixelId : change.pixelsC) {
        require(pixelId >= 0 && pixelId < input->getSize(), "region pixel must be in input bounds");
        inRegion[static_cast<size_t>(pixelId)] = true;
    }

    for (int pixelId = 0; pixelId < input->getSize(); ++pixelId) {
        if (inRegion[static_cast<size_t>(pixelId)]) {
            require((*updated)[pixelId] == static_cast<uint8_t>(change.b),
                    "region pixels must be assigned to b");
        } else {
            require((*updated)[pixelId] == (*input)[pixelId],
                    "pixels outside C must remain unchanged");
        }
    }
}

void testSubtreeRegionOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input->clone(), true, adj);
    DynamicComponentTree minTree(input->clone(), false, adj);

    const NodeId subtreeRoot = firstNonRootNode(maxTree);
    require(subtreeRoot != InvalidNode, "test image must have a non-root subtree");

    auto anchored = region_debug::buildSubtreeRegionChange<AltitudeType>(maxTree, minTree, subtreeRoot);
    auto oracle = rebuildTreesAfterChange(input, adj, anchored.change);

    requireOnlyRegionChanged(input, oracle.updatedImage, anchored.change);
    requireSameReconstruction(oracle);
}

void testExplicitConnectedNonSubtreeRegionOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);
    DynamicComponentTree maxTree(input->clone(), true, adj);

    const std::vector<PixelId> pixelsC = {0, 1};
    require(region_debug::isConnectedRegion(*adj, std::span<const PixelId>(pixelsC)),
            "explicit region must be connected");
    require(!matchesAnySubtreeFootprint(maxTree, pixelsC),
            "test region should not already be a subtree footprint");

    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(4),
        static_cast<AltitudeType>(7),
        true,
        false);
    auto oracle = rebuildTreesAfterChange(input, adj, change);

    require(change.direction == region_debug::RegionUpdateDirection::Raise,
            "explicit non-subtree region should be a raise update");
    requireOnlyRegionChanged(input, oracle.updatedImage, change);
    requireSameReconstruction(oracle);
}

void testExplicitContractionOracle() {
    auto input = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(input->getNumRows(), input->getNumCols(), 1.5);

    const std::vector<PixelId> pixelsC = {13, 14};
    require(region_debug::isConnectedRegion(*adj, std::span<const PixelId>(pixelsC)),
            "contraction region must be connected");

    auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
        std::span<const PixelId>(pixelsC),
        static_cast<AltitudeType>(6),
        static_cast<AltitudeType>(4),
        true,
        false);
    auto oracle = rebuildTreesAfterChange(input, adj, change);

    require(change.direction == region_debug::RegionUpdateDirection::Lower,
            "contraction region should be a lower update");
    requireOnlyRegionChanged(input, oracle.updatedImage, change);
    requireSameReconstruction(oracle);
}

} // namespace

int main() {
    testSubtreeRegionOracle();
    testExplicitConnectedNonSubtreeRegionOracle();
    testExplicitContractionOracle();
    return 0;
}
