#pragma once

#include "../../Common.hpp"

namespace debug::region {

enum class RegionUpdateDirection {
    Unchanged,
    Raise,
    Lower,
};

template<typename PixelType = AltitudeType>
RegionUpdateDirection directionFromLevels(PixelType a, PixelType b) {
    if (b > a) {
        return RegionUpdateDirection::Raise;
    }
    if (b < a) {
        return RegionUpdateDirection::Lower;
    }
    return RegionUpdateDirection::Unchanged;
}

} // namespace debug::region
