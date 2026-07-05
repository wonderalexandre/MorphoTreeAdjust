#pragma once

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../../DynamicComponentTree.hpp"

namespace debug::structural {

struct DynamicComponentTreeStructuralComparisonOptions {
    bool requireNoEmptyLiveNodes = true;
    bool requireStrictAltitudeOrder = true;
};

struct DynamicComponentTreeStructuralComparisonResult {
    bool equivalent = false;
    bool actualValid = false;
    bool expectedValid = false;
    int actualLiveNodes = 0;
    int expectedLiveNodes = 0;
    int actualRecords = 0;
    int expectedRecords = 0;
    std::string message;
};

namespace detail {

struct DynamicComponentTreeNodeKey {
    int altitude = 0;
    std::vector<PixelId> componentPixels;
    std::vector<PixelId> directProperParts;

    bool operator<(const DynamicComponentTreeNodeKey &other) const {
        if (altitude != other.altitude) {
            return altitude < other.altitude;
        }
        if (componentPixels != other.componentPixels) {
            return componentPixels < other.componentPixels;
        }
        return directProperParts < other.directProperParts;
    }

    bool operator==(const DynamicComponentTreeNodeKey &other) const {
        return altitude == other.altitude &&
               componentPixels == other.componentPixels &&
               directProperParts == other.directProperParts;
    }
};

struct DynamicComponentTreeNodeRecord {
    DynamicComponentTreeNodeKey node;
    bool hasParent = false;
    DynamicComponentTreeNodeKey parent;

    bool operator<(const DynamicComponentTreeNodeRecord &other) const {
        if (!(node == other.node)) {
            return node < other.node;
        }
        if (hasParent != other.hasParent) {
            return hasParent < other.hasParent;
        }
        if (!hasParent) {
            return false;
        }
        return parent < other.parent;
    }

    bool operator==(const DynamicComponentTreeNodeRecord &other) const {
        return node == other.node &&
               hasParent == other.hasParent &&
               (!hasParent || parent == other.parent);
    }
};

struct DynamicComponentTreeStructuralSnapshot {
    bool valid = false;
    int liveNodes = 0;
    std::vector<DynamicComponentTreeNodeRecord> records;
    std::string message;
};

inline std::string vectorPreview(const std::vector<PixelId> &values) {
    std::ostringstream out;
    out << "{";
    const size_t limit = std::min<size_t>(values.size(), 8);
    for (size_t i = 0; i < limit; ++i) {
        if (i != 0) {
            out << ",";
        }
        out << values[i];
    }
    if (values.size() > limit) {
        out << ",...";
    }
    out << "}";
    return out.str();
}

inline std::string nodeKeyToString(const DynamicComponentTreeNodeKey &key) {
    std::ostringstream out;
    out << "alt=" << key.altitude
        << " support_size=" << key.componentPixels.size()
        << " support=" << vectorPreview(key.componentPixels)
        << " direct_size=" << key.directProperParts.size()
        << " direct=" << vectorPreview(key.directProperParts);
    return out.str();
}

inline std::string nodeRecordToString(const DynamicComponentTreeNodeRecord &record) {
    std::ostringstream out;
    out << "node(" << nodeKeyToString(record.node) << ")";
    if (record.hasParent) {
        out << " parent(" << nodeKeyToString(record.parent) << ")";
    } else {
        out << " parent(ROOT)";
    }
    return out.str();
}

inline DynamicComponentTreeNodeKey makeNodeKey(const DynamicComponentTree &tree,
                                               NodeId nodeId) {
    DynamicComponentTreeNodeKey key;
    key.altitude = tree.getAltitude(nodeId);
    key.componentPixels = tree.getPixelsOfCC(nodeId);
    key.directProperParts.reserve(static_cast<size_t>(tree.getNumProperParts(nodeId)));

    for (PixelId pixelId : tree.getProperParts(nodeId)) {
        key.directProperParts.push_back(pixelId);
    }

    std::sort(key.componentPixels.begin(), key.componentPixels.end());
    std::sort(key.directProperParts.begin(), key.directProperParts.end());
    return key;
}

inline DynamicComponentTreeStructuralSnapshot makeStructuralSnapshot(
    const DynamicComponentTree &tree,
    const DynamicComponentTreeStructuralComparisonOptions &options,
    const std::string &label) {
    DynamicComponentTreeStructuralSnapshot snapshot;

    auto fail = [&](const std::string &message) {
        snapshot.valid = false;
        snapshot.message = label + ": " + message;
        return snapshot;
    };

    if (tree.getRoot() == InvalidNode) {
        return fail("root is invalid");
    }
    if (!tree.isAlive(tree.getRoot())) {
        return fail("root is not alive");
    }
    if (tree.getNodeParent(tree.getRoot()) != tree.getRoot()) {
        return fail("root parent is not itself");
    }

    std::vector<unsigned char> reachable(
        static_cast<size_t>(tree.getNumInternalNodeSlots()),
        0);
    std::vector<NodeId> nodes;
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            return fail("subtree traversal reached an invalid or dead node");
        }
        if (reachable[static_cast<size_t>(nodeId)] != 0) {
            return fail("subtree traversal visited a node twice");
        }
        reachable[static_cast<size_t>(nodeId)] = 1;
        nodes.push_back(nodeId);
    }

    int liveNodes = 0;
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }
        ++liveNodes;
        if (reachable[static_cast<size_t>(nodeId)] == 0) {
            return fail("alive node is unreachable from the root");
        }
    }
    snapshot.liveNodes = liveNodes;

    if (liveNodes != tree.getNumNodes()) {
        return fail("live node count differs from tree.getNumNodes()");
    }

    std::vector<unsigned char> seenProperParts(
        static_cast<size_t>(tree.getNumTotalProperParts()),
        0);
    int totalProperParts = 0;

    for (NodeId nodeId : nodes) {
        if (options.requireNoEmptyLiveNodes &&
            tree.getNumProperParts(nodeId) == 0) {
            return fail("live node without direct proper parts");
        }

        int countedChildren = 0;
        for (NodeId childId : tree.getChildren(nodeId)) {
            if (!tree.isNode(childId) || !tree.isAlive(childId)) {
                return fail("child is invalid or dead");
            }
            if (tree.getNodeParent(childId) != nodeId) {
                return fail("child parent does not point back to parent");
            }
            if (options.requireStrictAltitudeOrder) {
                if (tree.isMaxtree() &&
                    tree.getAltitude(nodeId) >= tree.getAltitude(childId)) {
                    std::ostringstream out;
                    out << "max-tree edge is not strictly increasing in altitude"
                        << " parent=" << nodeId
                        << " parent_alt=" << tree.getAltitude(nodeId)
                        << " child=" << childId
                        << " child_alt=" << tree.getAltitude(childId);
                    return fail(out.str());
                }
                if (!tree.isMaxtree() &&
                    tree.getAltitude(nodeId) <= tree.getAltitude(childId)) {
                    std::ostringstream out;
                    out << "min-tree edge is not strictly decreasing in altitude"
                        << " parent=" << nodeId
                        << " parent_alt=" << tree.getAltitude(nodeId)
                        << " child=" << childId
                        << " child_alt=" << tree.getAltitude(childId);
                    return fail(out.str());
                }
            }
            ++countedChildren;
        }
        if (countedChildren != tree.getNumChildren(nodeId)) {
            return fail("child count mismatch");
        }

        int countedProperParts = 0;
        for (PixelId pixelId : tree.getProperParts(nodeId)) {
            if (pixelId < 0 || pixelId >= tree.getNumTotalProperParts()) {
                return fail("direct proper part is out of bounds");
            }
            if (tree.getSmallestComponent(pixelId) != nodeId) {
                return fail("direct proper part owner mismatch");
            }
            if (seenProperParts[static_cast<size_t>(pixelId)] != 0) {
                return fail("direct proper part appears in more than one live node");
            }
            seenProperParts[static_cast<size_t>(pixelId)] = 1;
            ++countedProperParts;
        }
        if (countedProperParts != tree.getNumProperParts(nodeId)) {
            return fail("direct proper part count mismatch");
        }
        totalProperParts += countedProperParts;
    }

    if (totalProperParts != tree.getNumTotalProperParts()) {
        return fail("global direct proper part count mismatch");
    }
    for (PixelId pixelId = 0; pixelId < tree.getNumTotalProperParts(); ++pixelId) {
        if (seenProperParts[static_cast<size_t>(pixelId)] == 0) {
            return fail("some pixel has no direct proper-part owner");
        }
    }

    std::map<NodeId, DynamicComponentTreeNodeKey> keysByNode;
    std::map<DynamicComponentTreeNodeKey, int> keyCounts;
    for (NodeId nodeId : nodes) {
        auto key = makeNodeKey(tree, nodeId);
        keyCounts[key]++;
        keysByNode.emplace(nodeId, std::move(key));
    }
    for (const auto &[key, count] : keyCounts) {
        if (count > 1) {
            return fail("duplicate structural node key: " + nodeKeyToString(key));
        }
    }

    snapshot.records.reserve(nodes.size());
    for (NodeId nodeId : nodes) {
        DynamicComponentTreeNodeRecord record;
        record.node = keysByNode.at(nodeId);
        const NodeId parentId = tree.getNodeParent(nodeId);
        record.hasParent = parentId != nodeId;
        if (record.hasParent) {
            record.parent = keysByNode.at(parentId);
        }
        snapshot.records.push_back(std::move(record));
    }
    std::sort(snapshot.records.begin(), snapshot.records.end());

    snapshot.valid = true;
    snapshot.message = label + ": valid";
    return snapshot;
}

} // namespace detail

inline DynamicComponentTreeStructuralComparisonResult compareDynamicComponentTreeStructure(
    const DynamicComponentTree &actual,
    const DynamicComponentTree &expected,
    const DynamicComponentTreeStructuralComparisonOptions &options = {}) {
    DynamicComponentTreeStructuralComparisonResult result;

    if (actual.isMaxtree() != expected.isMaxtree()) {
        result.message = "tree type mismatch";
        return result;
    }
    if (actual.getNumRowsOfImage() != expected.getNumRowsOfImage() ||
        actual.getNumColsOfImage() != expected.getNumColsOfImage()) {
        result.message = "image domain mismatch";
        return result;
    }
    if (actual.getNumTotalProperParts() != expected.getNumTotalProperParts()) {
        result.message = "proper-part domain size mismatch";
        return result;
    }

    const auto actualSnapshot =
        detail::makeStructuralSnapshot(actual, options, "actual");
    const auto expectedSnapshot =
        detail::makeStructuralSnapshot(expected, options, "expected");

    result.actualValid = actualSnapshot.valid;
    result.expectedValid = expectedSnapshot.valid;
    result.actualLiveNodes = actualSnapshot.liveNodes;
    result.expectedLiveNodes = expectedSnapshot.liveNodes;
    result.actualRecords = static_cast<int>(actualSnapshot.records.size());
    result.expectedRecords = static_cast<int>(expectedSnapshot.records.size());

    if (!actualSnapshot.valid) {
        result.message = actualSnapshot.message;
        return result;
    }
    if (!expectedSnapshot.valid) {
        result.message = expectedSnapshot.message;
        return result;
    }

    if (actualSnapshot.records.size() != expectedSnapshot.records.size()) {
        std::ostringstream out;
        out << "structural node count mismatch: actual="
            << actualSnapshot.records.size()
            << " expected=" << expectedSnapshot.records.size();
        result.message = out.str();
        return result;
    }

    for (size_t index = 0; index < actualSnapshot.records.size(); ++index) {
        if (!(actualSnapshot.records[index] == expectedSnapshot.records[index])) {
            std::ostringstream out;
            out << "structural record mismatch at sorted index " << index
                << "; actual="
                << detail::nodeRecordToString(actualSnapshot.records[index])
                << "; expected="
                << detail::nodeRecordToString(expectedSnapshot.records[index]);
            result.message = out.str();
            return result;
        }
    }

    result.equivalent = true;
    result.message = "dynamic component trees are structurally equivalent";
    return result;
}

inline bool dynamicComponentTreesAreStructurallyEquivalent(
    const DynamicComponentTree &actual,
    const DynamicComponentTree &expected,
    const DynamicComponentTreeStructuralComparisonOptions &options = {}) {
    return compareDynamicComponentTreeStructure(actual, expected, options).equivalent;
}

} // namespace debug::structural
