#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"

namespace {

using Offset = AdjacencyRelation::Offset;

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::vector<int> collect_range(AdjacencyRelation &adjacency) {
    std::vector<int> values;
    for (int value : adjacency) {
        values.push_back(value);
    }
    return values;
}

std::vector<int> expected_rectangular_neighbors(int numCols, int row, int col, int halfHeight, int halfWidth, bool forwardOnly) {
    std::vector<int> values;
    for (int dy = -halfHeight; dy <= halfHeight; ++dy) {
        for (int dx = -halfWidth; dx <= halfWidth; ++dx) {
            if (dy == 0 && dx == 0) {
                continue;
            }
            if (forwardOnly && !(dy > 0 || (dy == 0 && dx > 0))) {
                continue;
            }
            values.push_back((row + dy) * numCols + (col + dx));
        }
    }
    std::sort(values.begin(), values.end());
    return values;
}

void test_circular_constructor_preserves_existing_behavior() {
    AdjacencyRelation adjacency(5, 5, 1.5);

    require(adjacency.hasRadius(), "circular adjacency should retain radius metadata");
    require(adjacency.getSize() == 9, "radius 1.5 should produce 8-neighborhood plus origin");
    require(adjacency.isAdjacent(12, 18), "diagonal pixels should be adjacent for radius 1.5");
    require(!adjacency.isAdjacent(12, 24), "far pixels should not be adjacent");
}

void test_rectangular_factory_builds_expected_offsets_and_forward_half() {
    auto adjacency = AdjacencyRelation::rectangular(7, 8, 1, 2);

    require(!adjacency.hasRadius(), "rectangular adjacency should not advertise circular radius");
    require(adjacency.getSize() == 15, "1x2 rectangle should contain 15 offsets including origin");
    require(adjacency.isAdjacent(3 * 8 + 3, 4 * 8 + 5), "offset (+1,+2) should belong to rectangular adjacency");
    require(!adjacency.isAdjacent(3 * 8 + 3, 5 * 8 + 3), "offset (+2,0) should be outside rectangular adjacency");

    auto neighbors = collect_range(adjacency.getNeighborPixels(3, 3));
    auto expected = expected_rectangular_neighbors(8, 3, 3, 1, 2, false);
    std::sort(neighbors.begin(), neighbors.end());
    require(neighbors == expected, "rectangular neighborhood should match the full stencil");

    auto forwardNeighbors = collect_range(adjacency.getNeighborPixelsForward(3, 3));
    auto expectedForward = expected_rectangular_neighbors(8, 3, 3, 1, 2, true);
    std::sort(forwardNeighbors.begin(), forwardNeighbors.end());
    require(forwardNeighbors == expectedForward, "forward rectangular neighborhood should keep one representative per symmetric pair");
}

void test_generic_constructor_accepts_custom_symmetric_offsets() {
    const std::vector<Offset> offsets = {
        {0, 0},
        {-1, 0},
        {1, 0},
        {0, -2},
        {0, 2}
    };

    AdjacencyRelation adjacency(6, 7, offsets);
    require(!adjacency.hasRadius(), "custom offset adjacency should not advertise circular radius");
    require(adjacency.getSize() == 5, "custom symmetric adjacency should preserve unique offsets");
    require(adjacency.isAdjacent(2 * 7 + 3, 2 * 7 + 5), "custom horizontal jump should be adjacent");
    require(!adjacency.isAdjacent(2 * 7 + 3, 3 * 7 + 4), "offset absent from the stencil should not be adjacent");

    auto neighbors = collect_range(adjacency.getNeighborPixels(2, 3));
    std::sort(neighbors.begin(), neighbors.end());
    require(neighbors == std::vector<int>({10, 15, 19, 24}), "custom adjacency neighbors should match the declared offsets");
}

void test_generic_constructor_rejects_non_symmetric_offsets() {
    bool threw = false;
    try {
        AdjacencyRelation adjacency(5, 5, std::vector<Offset>{{0, 0}, {1, 0}});
        (void) adjacency;
    } catch (const std::invalid_argument &) {
        threw = true;
    }

    require(threw, "adjacency should reject offset sets without central symmetry");
}

} // namespace

int main() {
    try {
        test_circular_constructor_preserves_existing_behavior();
        test_rectangular_factory_builds_expected_offsets_and_forward_half();
        test_generic_constructor_accepts_custom_symmetric_offsets();
        test_generic_constructor_rejects_non_symmetric_offsets();
    } catch (const std::exception &ex) {
        std::cerr << "[adjacency_relation_unit_tests] " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}
