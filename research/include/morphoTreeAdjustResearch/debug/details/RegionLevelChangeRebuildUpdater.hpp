#pragma once

#include <span>
#include <stdexcept>

#include "RegionLevelChange.hpp"
#include "../../AdjacencyRelation.hpp"
#include "../../Common.hpp"
#include "../../DynamicComponentTree.hpp"

namespace debug::region {

class RegionLevelChangeRebuildUpdater {
public:
    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr applyToImage(const ImageUInt8Ptr &input,
                                      const RegionLevelChange<PixelType> &change) {
        if (!input) {
            throw std::runtime_error("RegionLevelChangeRebuildUpdater requires a valid input image");
        }

        auto output = input->clone();
        for (PixelId pixelId : change.pixelsC) {
            if (pixelId < 0 || pixelId >= output->getSize()) {
                throw std::runtime_error("RegionLevelChangeRebuildUpdater received an out-of-bounds region pixel");
            }
            (*output)[pixelId] = static_cast<uint8_t>(change.b);
        }
        return output;
    }

    template<typename PixelType = AltitudeType>
    static ImageUInt8Ptr updatePair(DynamicComponentTree &minTree,
                                    DynamicComponentTree &maxTree,
                                    const AdjacencyRelationPtr &adjacency,
                                    const RegionLevelChange<PixelType> &change) {
        if (!adjacency) {
            throw std::runtime_error("RegionLevelChangeRebuildUpdater requires a valid adjacency relation");
        }
        if (change.pixelsC.empty()) {
            throw std::runtime_error("RegionLevelChangeRebuildUpdater requires a non-empty region C");
        }
        if (!change.assumeConnected &&
            !isConnectedRegion(*adjacency, std::span<const PixelId>(change.pixelsC))) {
            throw std::runtime_error("RegionLevelChangeRebuildUpdater requires connected region C");
        }

        auto updatedImage = applyToImage(maxTree.reconstructionImage(), change);
        maxTree.build(updatedImage->clone(), true, adjacency);
        minTree.build(updatedImage->clone(), false, adjacency);
        return updatedImage;
    }
};

} // namespace debug::region
