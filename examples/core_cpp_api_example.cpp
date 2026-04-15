#include <MorphoTreeAdjust/MorphoTreeAdjust.hpp>

#include <memory>
#include <vector>

int main() {
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

    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);
    DynamicComponentTree maxTree(image->clone(), true, adj);
    DynamicComponentTree minTree(image->clone(), false, adj);

    DynamicComponentTreeAdjustment<AltitudeType> adjust(&minTree, &maxTree, *adj);
    DynamicAreaComputer maxAreaComputer(&maxTree);
    DynamicAreaComputer minAreaComputer(&minTree);
    auto maxArea = maxAreaComputer.compute();
    auto minArea = minAreaComputer.compute();
    adjust.setAttributeComputer(minAreaComputer, maxAreaComputer,
                                std::span<float>(minArea), std::span<float>(maxArea));

    ComponentTreeCasf<AltitudeType> casf(image->clone(), adj, AREA);
    auto filtered = casf.filter(std::vector<int>{1, 2}, ComponentTreeCasf<AltitudeType>::Mode::Hybrid);

    return filtered == nullptr ? 1 : 0;
}
