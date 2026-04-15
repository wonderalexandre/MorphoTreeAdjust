#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../morphoTreeAdjust/include/MorphoTreeAdjust.hpp"

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

} // namespace

int main() {
    auto image = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.5);

    DynamicComponentTree maxTree(image->clone(), true, adj);
    DynamicComponentTree minTree(image->clone(), false, adj);
    require(maxTree.getRoot() != InvalidNode, "max-tree root must be valid");
    require(minTree.getRoot() != InvalidNode, "min-tree root must be valid");

    DynamicAreaComputer maxAreaComputer(&maxTree);
    const auto maxArea = maxAreaComputer.compute();
    require(maxArea.size() == static_cast<size_t>(maxTree.getGlobalIdSpaceSize()),
            "area buffer size must match the tree global id space");

    DynamicComponentTreeAdjustment<AltitudeType> subtreeAdjust(&minTree, &maxTree, *adj);
    DynamicComponentTreeAdjustmentLeaf<AltitudeType> leafAdjust(&minTree, &maxTree, *adj);
    static_cast<void>(subtreeAdjust);
    static_cast<void>(leafAdjust);

    ComponentTreeCasf<AltitudeType> casf(image->clone(), adj, AREA);
    auto filtered = casf.filter(std::vector<int>{1, 2}, ComponentTreeCasf<AltitudeType>::Mode::Hybrid);
    require(filtered != nullptr, "ComponentTreeCasf must return a valid filtered image");
    require(filtered->getNumRows() == image->getNumRows(), "filtered image must preserve the number of rows");
    require(filtered->getNumCols() == image->getNumCols(), "filtered image must preserve the number of cols");

    return 0;
}
