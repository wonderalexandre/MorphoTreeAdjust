#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/DynamicComponentTree.hpp"

namespace {

using tree_t = DynamicComponentTree;

/**
 * @brief Reference image used to validate the initial topology of the dynamic tree.
 *
 * The expected topology of this max-tree serves as the baseline for checking
 * ids, children, proper-part ownership, and structural invariants after local
 * mutations.
 */
ImageUInt8Ptr make_demo_image() {
    auto image = ImageUInt8::create(12, 12);
    const std::vector<uint8_t> raw = {
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
        2, 5, 5, 5, 3, 3, 3, 2, 2, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 5, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 6, 2, 3, 1,
        2, 5, 6, 5, 3, 3, 3, 2, 5, 2, 3, 1,
        2, 5, 5, 5, 3, 3, 3, 2, 2, 2, 3, 1,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1,
        2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 1,
        6, 6, 6, 6, 6, 2, 4, 4, 4, 4, 4, 1,
        7, 3, 6, 4, 6, 2, 4, 8, 4, 8, 4, 1,
        7, 8, 6, 6, 6, 2, 4, 4, 4, 4, 4, 1,
        2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1
    };

    for (size_t i = 0; i < raw.size(); ++i) {
        (*image)[(int) i] = raw[i];
    }
    return image;
}

tree_t make_component_tree_from_demo(bool isMaxtree) {
    auto image = make_demo_image();
    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    return tree_t(image, isMaxtree, adj);
}

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template<typename RangeT>
std::vector<int> collect_range(const RangeT &range) {
    std::vector<int> values;
    for (int value : range) {
        values.push_back(value);
    }
    return values;
}

std::vector<int> collect_proper_parts(const tree_t &tree, int nodeId) {
    return collect_range(tree.getProperParts(nodeId));
}

bool vector_equal(const std::vector<int> &lhs, const std::vector<int> &rhs) {
    return lhs == rhs;
}

bool vector_same(std::vector<int> lhs, std::vector<int> rhs) {
    std::sort(lhs.begin(), lhs.end());
    std::sort(rhs.begin(), rhs.end());
    return lhs == rhs;
}

int nid(const tree_t &tree, int oldNodeId) {
    return oldNodeId - tree.getNumTotalProperParts();
}

int count_alive_nodes(const tree_t &tree) {
    int count = 0;
    for (int nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId)) {
            ++count;
        }
    }
    return count;
}

bool has_unreachable_alive_nodes(const tree_t &tree) {
    std::vector<bool> reachable((size_t) tree.getNumInternalNodeSlots(), false);
    if (tree.getRoot() != InvalidNode && tree.isAlive(tree.getRoot())) {
        for (int nodeId : tree.getNodeSubtree(tree.getRoot())) {
            reachable[(size_t) nodeId] = true;
        }
    }

    for (int nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId) && !reachable[(size_t) nodeId]) {
            return true;
        }
    }
    return false;
}

/**
 * @brief Checks global invariants of the dynamic hierarchy.
 *
 * The goal is to check the minimum conditions that must remain true after any
 * structural mutation: valid root, parent/child symmetry, reachability of live
 * nodes, and complete partitioning of the proper parts.
 */
void require_tree_consistency(const tree_t &tree) {
    require(tree.getRoot() != InvalidNode, "root must be valid");
    require(tree.isAlive(tree.getRoot()), "root must be alive");
    require(tree.getNodeParent(tree.getRoot()) == tree.getRoot(), "root parent must be itself");
    require(!has_unreachable_alive_nodes(tree), "all alive nodes must be reachable from root");

    int totalProperParts = 0;
    std::vector<bool> seen((size_t) tree.getNumTotalProperParts(), false);

    for (int nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }

        int countedChildren = 0;
        for (int childId : tree.getChildren(nodeId)) {
            require(tree.isAlive(childId), "child must be alive");
            require(tree.getNodeParent(childId) == nodeId, "child parent must match");
            require(tree.hasChild(nodeId, childId), "parent-child relation must be symmetric");
            ++countedChildren;
        }
        require(countedChildren == tree.getNumChildren(nodeId), "child count mismatch");

        int countedProperParts = 0;
        for (int pixelId : tree.getProperParts(nodeId)) {
            require(pixelId >= 0 && pixelId < tree.getNumTotalProperParts(), "proper part must be a leaf id");
            require(tree.getSmallestComponent(pixelId) == nodeId, "proper part owner mismatch");
            require(!seen[(size_t) pixelId], "proper part visited twice");
            seen[(size_t) pixelId] = true;
            ++countedProperParts;
        }
        require(countedProperParts == tree.getNumProperParts(nodeId), "proper part count mismatch");
        totalProperParts += countedProperParts;
    }

    require(totalProperParts == tree.getNumTotalProperParts(), "global proper part count mismatch");
    require(std::all_of(seen.begin(), seen.end(), [](bool v) { return v; }), "proper part partition is incomplete");
}

void test_initial_maxtree_matches_reference() {
    // The initial construction must reproduce the expected reference topology.
    auto tree = make_component_tree_from_demo(true);

    require(tree.getNumTotalProperParts() == 144, "unexpected number of proper parts");
    require(tree.getNumInternalNodeSlots() == 15, "unexpected node domain size");
    require(count_alive_nodes(tree) == 15, "unexpected number of alive internal nodes");
    require(tree.getRoot() == nid(tree, 144), "unexpected root id");
    require(tree.getAltitude(tree.getRoot()) == 1, "unexpected root altitude");
    require(tree.getSmallestComponent(32) == nid(tree, 151), "unexpected owner for pixel 32");
    require(tree.getSmallestComponent(44) == nid(tree, 153), "unexpected owner for pixel 44");
    require(tree.getSmallestComponent(109) == nid(tree, 147), "unexpected owner for pixel 109");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 145))), {nid(tree, 146), nid(tree, 147), nid(tree, 151)}), "unexpected children for node 145");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 146))), {nid(tree, 148), nid(tree, 150)}), "unexpected children for node 146");
    require(vector_same(collect_proper_parts(tree, nid(tree, 151)), {32, 56}), "unexpected proper parts for node 151");
    require(vector_same(collect_proper_parts(tree, nid(tree, 153)), {44}), "unexpected proper parts for node 153");
    require_tree_consistency(tree);
}

void test_topology_primitives_preserve_contracts() {
    // detach/attach/move/removeChild must preserve parentage links and the local sibling order.
    auto tree = make_component_tree_from_demo(true);
    require_tree_consistency(tree);

    tree.detachNode(nid(tree, 150));
    require(tree.getNodeParent(nid(tree, 150)) == nid(tree, 150), "detached node must become its own parent");
    require(!tree.hasChild(nid(tree, 146), nid(tree, 150)), "detached node should not remain a child");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 146))), {nid(tree, 148)}), "unexpected children after detach");

    tree.attachNode(nid(tree, 147), nid(tree, 150));
    require(tree.getNodeParent(nid(tree, 150)) == nid(tree, 147), "attachNode must update parent");
    require(tree.hasChild(nid(tree, 147), nid(tree, 150)), "attachNode must update child list");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 147))), {nid(tree, 149), nid(tree, 150)}), "unexpected children after attach");

    tree.moveNode(nid(tree, 148), nid(tree, 147));
    require(tree.getNodeParent(nid(tree, 148)) == nid(tree, 147), "moveNode must update parent");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 147))), {nid(tree, 149), nid(tree, 150), nid(tree, 148)}), "unexpected children after moveNode");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 146))), {}), "source node should be empty after moveNode");

    tree.removeChild(nid(tree, 147), nid(tree, 150), false);
    require(tree.getNodeParent(nid(tree, 150)) == nid(tree, 150), "removeChild(false) must detach node");
    require(!tree.hasChild(nid(tree, 147), nid(tree, 150)), "removed child must disappear from parent");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 147))), {nid(tree, 149), nid(tree, 148)}), "unexpected children after removeChild");

    tree.attachNode(nid(tree, 146), nid(tree, 150));
    require(tree.getNodeParent(nid(tree, 150)) == nid(tree, 146), "reattach must update parent");
    require(vector_equal(collect_range(tree.getChildren(nid(tree, 146))), {nid(tree, 150)}), "unexpected children after reattach");
    require_tree_consistency(tree);
}

void test_move_proper_parts_operations() {
    // Single and bulk proper-part moves must keep ownership and direct lists consistent.
    {
        auto tree = make_component_tree_from_demo(true);
        require_tree_consistency(tree);

        const int sourceNode = tree.getSmallestComponent(32);
        const int targetNode = tree.getSmallestComponent(44);
        tree.moveProperPart(targetNode, sourceNode, 32);

        require(tree.getSmallestComponent(32) == targetNode, "moveProperPart must update owner");
        require(tree.getNumProperParts(sourceNode) == 1, "source node should lose one proper part");
        require(tree.getNumProperParts(targetNode) == 2, "target node should gain one proper part");
        require(vector_same(collect_proper_parts(tree, sourceNode), {56}), "unexpected source proper parts after moveProperPart");
        require(vector_same(collect_proper_parts(tree, targetNode), {44, 32}), "unexpected target proper parts after moveProperPart");
        require_tree_consistency(tree);
    }

    {
        auto tree = make_component_tree_from_demo(true);
        require_tree_consistency(tree);

        tree.moveProperParts(nid(tree, 153), nid(tree, 151));

        require(tree.getNumProperParts(nid(tree, 151)) == 0, "source node should be emptied by moveProperParts");
        require(tree.getNumProperParts(nid(tree, 153)) == 3, "target node should receive all direct proper parts");
        require(vector_same(collect_proper_parts(tree, nid(tree, 153)), {44, 32, 56}), "unexpected target proper parts after moveProperParts");
        require(tree.getSmallestComponent(32) == nid(tree, 153), "owner of 32 must be updated");
        require(tree.getSmallestComponent(56) == nid(tree, 153), "owner of 56 must be updated");
        require_tree_consistency(tree);
    }
}

void test_move_children_merge_and_prune() {
    // moveChildren, mergeNodeIntoParent, and pruneNode are the core mutations of incremental adjustment.
    {
        auto tree = make_component_tree_from_demo(true);
        require_tree_consistency(tree);

        tree.moveChildren(nid(tree, 147), nid(tree, 146));
        require(tree.getNumChildren(nid(tree, 146)) == 0, "source must lose all children after moveChildren");
        require(tree.getNumChildren(nid(tree, 147)) == 3, "target must gain all children after moveChildren");
        require(vector_equal(collect_range(tree.getChildren(nid(tree, 147))), {nid(tree, 149), nid(tree, 148), nid(tree, 150)}), "child order must be preserved");
        require(tree.getNodeParent(nid(tree, 148)) == nid(tree, 147), "moved child 148 has wrong parent");
        require(tree.getNodeParent(nid(tree, 150)) == nid(tree, 147), "moved child 150 has wrong parent");
        require_tree_consistency(tree);
    }

    {
        auto tree = make_component_tree_from_demo(true);
        require_tree_consistency(tree);

        const int parentNode = tree.getNodeParent(nid(tree, 150));
        const int oldAlive = count_alive_nodes(tree);
        const int oldParentProperParts = tree.getNumProperParts(parentNode);
        tree.mergeNodeIntoParent(nid(tree, 150));

        require(!tree.isAlive(nid(tree, 150)), "merged node must be released");
        require(count_alive_nodes(tree) == oldAlive - 1, "mergeNodeIntoParent must remove exactly one internal node");
        require(tree.getNumProperParts(parentNode) == oldParentProperParts + 12, "parent must inherit proper parts of merged node");
        require(vector_equal(collect_range(tree.getChildren(parentNode)), {nid(tree, 148), nid(tree, 152)}), "children must be promoted in order");
        require(tree.getNodeParent(nid(tree, 152)) == parentNode, "promoted child must point to parent");
        require(tree.getSmallestComponent(13) == parentNode, "proper parts must be re-owned by parent");
        require_tree_consistency(tree);
    }

    {
        auto tree = make_component_tree_from_demo(true);
        require_tree_consistency(tree);

        const int parentNode = tree.getNodeParent(nid(tree, 150));
        const int oldAlive = count_alive_nodes(tree);
        const int oldParentProperParts = tree.getNumProperParts(parentNode);
        tree.pruneNode(nid(tree, 150));

        require(!tree.isAlive(nid(tree, 150)), "pruned node must be released");
        require(count_alive_nodes(tree) == oldAlive - 2, "pruneNode must remove subtree internal nodes");
        require(tree.getNumProperParts(parentNode) == oldParentProperParts + 15, "parent must inherit pruned proper parts");
        require(vector_equal(collect_range(tree.getChildren(parentNode)), {nid(tree, 148)}), "remaining child list is incorrect after prune");
        require(tree.getSmallestComponent(13) == parentNode, "proper parts from pruned subtree must move to parent");
        require(tree.getSmallestComponent(26) == parentNode, "descendant proper parts must move to parent");
        require(tree.getSmallestComponent(44) == nid(tree, 153), "unrelated subtree must remain untouched");
        require_tree_consistency(tree);
    }
}

void test_surviving_node_ids_stay_stable_under_topology_updates() {
    // Purely topological operations must not renumber surviving internal nodes.
    auto tree = make_component_tree_from_demo(true);
    require_tree_consistency(tree);

    const std::vector<int> trackedNodes = {
        nid(tree, 148), nid(tree, 149), nid(tree, 150), nid(tree, 151), nid(tree, 152),
        nid(tree, 153), nid(tree, 154), nid(tree, 155), nid(tree, 156), nid(tree, 157),
        nid(tree, 158)
    };

    std::vector<std::vector<int>> trackedProperParts;
    trackedProperParts.reserve(trackedNodes.size());
    for (int nodeId : trackedNodes) {
        require(tree.isAlive(nodeId), "tracked node must start alive");
        trackedProperParts.push_back(collect_proper_parts(tree, nodeId));
    }

    tree.detachNode(nid(tree, 150));
    tree.attachNode(nid(tree, 147), nid(tree, 150));
    tree.moveNode(nid(tree, 148), nid(tree, 147));
    tree.moveChildren(nid(tree, 147), nid(tree, 146));

    for (size_t i = 0; i < trackedNodes.size(); ++i) {
        const int nodeId = trackedNodes[i];
        require(tree.isAlive(nodeId), "surviving node must stay alive");
        require(vector_same(collect_proper_parts(tree, nodeId), trackedProperParts[i]),
                "pure topology updates must keep direct proper parts attached to the same surviving node");
        for (int pixelId : trackedProperParts[i]) {
            require(tree.getSmallestComponent(pixelId) == nodeId,
                    "proper part owner must stay stable under pure topology updates");
        }
    }

    require(tree.getRoot() == nid(tree, 144), "root id should stay unchanged");
    require(tree.getNodeParent(nid(tree, 145)) == nid(tree, 144), "node 145 must remain under the root");
    require(tree.getNodeParent(nid(tree, 146)) == nid(tree, 145), "node 146 must remain under 145");
    require(tree.getNodeParent(nid(tree, 147)) == nid(tree, 145), "node 147 must remain under 145");
    require(tree.getNodeParent(nid(tree, 148)) == nid(tree, 147), "node 148 must move under 147");
    require(tree.getNodeParent(nid(tree, 149)) == nid(tree, 147), "node 149 must remain under 147");
    require(tree.getNodeParent(nid(tree, 150)) == nid(tree, 147), "node 150 must move under 147");
    require(tree.getNodeParent(nid(tree, 151)) == nid(tree, 145), "node 151 must remain under 145");
    require(tree.getNodeParent(nid(tree, 152)) == nid(tree, 150), "node 152 must remain under 150");
    require(tree.getNodeParent(nid(tree, 153)) == nid(tree, 151), "node 153 must remain under 151");
    require_tree_consistency(tree);
}

void test_subtree_traversal_tracks_mutated_topology() {
    // Depth-first subtree traversal must reflect the topology after moveChildren.
    auto tree = make_component_tree_from_demo(true);
    require_tree_consistency(tree);

    tree.moveChildren(nid(tree, 147), nid(tree, 146));

    require(vector_equal(collect_range(tree.getChildren(nid(tree, 147))),
                         {nid(tree, 149), nid(tree, 148), nid(tree, 150)}),
            "children order should reflect appended moved children");
    require(vector_equal(collect_range(tree.getNodeSubtree(nid(tree, 147))),
                         {nid(tree, 147), nid(tree, 149), nid(tree, 154), nid(tree, 155), nid(tree, 158),
                          nid(tree, 148), nid(tree, 156), nid(tree, 157), nid(tree, 150), nid(tree, 152)}),
            "subtree traversal should reflect the updated topology");
    require(tree.countDescendants(nid(tree, 147)) == 9, "descendant count should match subtree size minus root");
    require(tree.countProperPartsInDescendants(nid(tree, 147)) == 44,
            "descendant proper part count should match the moved subtree content");
    require_tree_consistency(tree);
}

void test_root_and_node_slot_reuse() {
    // setRoot, releaseNode, and allocateNode exercise the lifecycle of internal ids.
    {
        auto tree = make_component_tree_from_demo(true);
        const int oldRoot = tree.getRoot();
        tree.setRoot(nid(tree, 145));

        require(tree.getRoot() == nid(tree, 145), "setRoot must update root id");
        require(tree.getNodeParent(nid(tree, 145)) == nid(tree, 145), "new root must parent itself");
        require(tree.getNodeParent(oldRoot) == oldRoot, "old root becomes detached");
        require(!tree.hasChild(oldRoot, nid(tree, 145)), "new root must not remain child of old root");
        require(has_unreachable_alive_nodes(tree), "detached old rooted component should become unreachable");
        require(collect_range(tree.getNodeSubtree(nid(tree, 145))).front() == nid(tree, 145), "subtree traversal must start at new root");
    }

    {
        auto tree = make_component_tree_from_demo(true);
        require_tree_consistency(tree);

        const int reusableNode = nid(tree, 153);
        const int parentNode = tree.getNodeParent(reusableNode);
        tree.moveProperPart(nid(tree, 151), reusableNode, 44);
        tree.removeChild(parentNode, reusableNode, false);
        tree.releaseNode(reusableNode);

        require(!tree.isAlive(reusableNode), "released node must not stay alive");
        require(vector_equal(collect_range(tree.getChildren(parentNode)), {}), "released node must disappear from children");
        require(vector_same(collect_proper_parts(tree, nid(tree, 151)), {32, 56, 44}), "proper parts must stay in sibling owner");

        const int newNode = tree.allocateNode();
        require(newNode == reusableNode, "allocateNode should reuse the freed id first");
        require(tree.isAlive(newNode), "reallocated node must be alive");
        require(tree.getNodeParent(newNode) == newNode, "fresh node must initially parent itself");
        require(tree.getNumChildren(newNode) == 0, "fresh node must start childless");
        require(tree.getNumProperParts(newNode) == 0, "fresh node must start with no proper parts");
        tree.attachNode(parentNode, newNode);
        require(tree.getNodeParent(newNode) == parentNode, "reattached node must point to parent");
        require(tree.hasChild(parentNode, newNode), "parent must contain reattached node");
        require_tree_consistency(tree);
    }
}

void test_remove_child_with_release_node_removes_empty_leaf_slot() {
    // removeChild(..., true) must detach and release a leaf node with no proper parts.
    auto tree = make_component_tree_from_demo(true);
    require_tree_consistency(tree);

    const int parentNode = nid(tree, 151);
    const int releasedNode = nid(tree, 153);
    const int oldAlive = count_alive_nodes(tree);

    tree.moveProperPart(parentNode, releasedNode, 44);
    tree.removeChild(parentNode, releasedNode, true);

    require(!tree.isAlive(releasedNode), "released child must no longer be alive");
    require(count_alive_nodes(tree) == oldAlive - 1, "removeChild(true) must release exactly one node");
    require(vector_equal(collect_range(tree.getChildren(parentNode)), {}), "released child must disappear from children");
    require(vector_same(collect_proper_parts(tree, parentNode), {32, 56, 44}),
            "released child proper parts must stay with the chosen surviving owner");

    const int reusedNode = tree.allocateNode();
    require(reusedNode == releasedNode, "released node id should be reusable first");
    require(tree.isAlive(reusedNode), "reallocated node must become alive again");
    require(tree.getNodeParent(reusedNode) == reusedNode, "freshly reallocated node must start detached");
    tree.attachNode(parentNode, reusedNode);
    require(tree.getNodeParent(reusedNode) == parentNode, "reattached node must point to the original parent");
    require_tree_consistency(tree);
}

} // namespace

int main() {
    try {
        test_initial_maxtree_matches_reference();
        test_topology_primitives_preserve_contracts();
        test_move_proper_parts_operations();
        test_move_children_merge_and_prune();
        test_surviving_node_ids_stay_stable_under_topology_updates();
        test_subtree_traversal_tracks_mutated_topology();
        test_root_and_node_slot_reuse();
        test_remove_child_with_release_node_removes_empty_leaf_slot();
    } catch (const std::exception &e) {
        std::cerr << "dynamic_component_tree_unit_tests: FAIL\n" << e.what() << "\n";
        return 1;
    }

    std::cout << "dynamic_component_tree_unit_tests: OK\n";
    return 0;
}
