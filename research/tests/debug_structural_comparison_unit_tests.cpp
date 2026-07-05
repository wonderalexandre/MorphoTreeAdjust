#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/StructuralComparison.hpp"

namespace structural_debug = debug::structural;

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

DynamicComponentTree make_tree(const ImageUInt8Ptr &image, bool isMaxtree) {
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.5);
    return DynamicComponentTree(image->clone(), isMaxtree, adjacency);
}

NodeId firstNonRootGrandchild(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (nodeId == tree.getRoot()) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId != InvalidNode &&
            parentId != tree.getRoot() &&
            tree.getAltitude(tree.getRoot()) < tree.getAltitude(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

void testEquivalentRebuiltTreesHaveSameStructure() {
    auto image = make_demo_image();
    auto actual = make_tree(image, true);
    auto expected = make_tree(image, true);

    const auto comparison =
        structural_debug::compareDynamicComponentTreeStructure(actual, expected);

    require(comparison.equivalent,
            "rebuilt trees from the same image should be structurally equivalent: " +
                comparison.message);
    require(comparison.actualLiveNodes == comparison.expectedLiveNodes,
            "equivalent trees should expose the same number of live nodes");
}

void testEmptyConnectorIsRejectedEvenWhenReconstructionMatches() {
    auto image = make_demo_image();
    auto actual = make_tree(image, true);
    auto expected = make_tree(image, true);

    const NodeId connector = actual.createDetachedNode(actual.getAltitude(actual.getRoot()) + 1);
    actual.attachNode(actual.getRoot(), connector);

    require(actual.reconstructionImage()->isEqual(expected.reconstructionImage()),
            "empty connector fixture should preserve reconstructed image");

    const auto comparison =
        structural_debug::compareDynamicComponentTreeStructure(actual, expected);

    require(!comparison.equivalent,
            "structural comparison should reject an empty live connector");
    require(!comparison.actualValid,
            "tree with an empty connector should be structurally invalid");
}

void testChangedParentRelationIsDetectedEvenWhenReconstructionMatches() {
    auto image = make_demo_image();
    auto actual = make_tree(image, true);
    auto expected = make_tree(image, true);

    const NodeId movedNode = firstNonRootGrandchild(actual);
    require(movedNode != InvalidNode,
            "fixture should contain a max-tree grandchild");
    actual.moveNode(movedNode, actual.getRoot());

    require(actual.reconstructionImage()->isEqual(expected.reconstructionImage()),
            "parent-move fixture should preserve reconstructed image");

    const auto comparison =
        structural_debug::compareDynamicComponentTreeStructure(actual, expected);

    require(!comparison.equivalent,
            "structural comparison should detect changed parent/support relations");
    require(comparison.actualValid,
            "parent-move fixture should remain internally valid");
    require(comparison.expectedValid,
            "rebuilt expected tree should be internally valid");
}

} // namespace

int main() {
    testEquivalentRebuiltTreesHaveSameStructure();
    testEmptyConnectorIsRejectedEvenWhenReconstructionMatches();
    testChangedParentRelationIsDetectedEvenWhenReconstructionMatches();
    return 0;
}
