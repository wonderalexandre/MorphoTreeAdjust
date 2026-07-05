#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/NonIncreasingMcpcDirectPreflight.hpp"
#include "../include/morphoTreeAdjustResearch/details/NonIncreasingCasfKernel.hpp"

#include "../../dev-tools/external/stb/stb_image.h"

#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#elif defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

namespace noninc = morphotree::details::non_increasing;

namespace {

template<typename Adjust>
concept HasBooleanMergeNodeOverload = requires(
    Adjust adjust,
    DynamicComponentTree *primalTree,
    DynamicComponentTree *dualTree,
    NodeId nodeId) {
    adjust.mergeNodeIntoParentAndUpdate(primalTree, dualTree, nodeId, true);
};

static_assert(
    !HasBooleanMergeNodeOverload<DualMinMaxTreeIncrementalFilter<AltitudeType>>,
    "mergeNodeIntoParentAndUpdate must not expose a boolean target-creation overload");

void require(bool condition, const std::string &message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_tree_has_no_empty_live_nodes(const DynamicComponentTree &tree,
                                          const std::string &treeName) {
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        require(tree.getNumProperParts(nodeId) > 0,
                treeName + ": live node without proper parts at id " + std::to_string(nodeId));
    }
}

ImageUInt8Ptr load_gray_image(const std::filesystem::path &path) {
    int numCols = 0;
    int numRows = 0;
    int numChannels = 0;
    uint8_t *data = stbi_load(path.string().c_str(), &numCols, &numRows, &numChannels, 1);
    if (data == nullptr) {
        throw std::runtime_error("failed to load image: " + path.string());
    }

    auto image = ImageUInt8::create(numRows, numCols);
    std::copy(data, data + (numRows * numCols), image->rawData());
    stbi_image_free(data);
    return image;
}

ImageUInt8Ptr crop_top_left(const ImageUInt8Ptr &image, int rows, int cols) {
    require(image != nullptr, "cannot crop a null image");
    require(rows > 0 && cols > 0, "crop dimensions must be positive");
    require(rows <= image->getNumRows() && cols <= image->getNumCols(),
            "crop dimensions must fit inside the image");

    auto crop = ImageUInt8::create(rows, cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            (*crop)[row * cols + col] =
                (*image)[row * image->getNumCols() + col];
        }
    }
    return crop;
}

ImageUInt8Ptr load_dat_crop(const std::string &relativePath, int rows, int cols) {
    const auto repoRoot = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path();
    return crop_top_left(load_gray_image(repoRoot / relativePath), rows, cols);
}

ImageUInt8Ptr make_monotone_line_image() {
    auto image = ImageUInt8::create(1, 4);
    (*image)[0] = 10;
    (*image)[1] = 20;
    (*image)[2] = 30;
    (*image)[3] = 40;
    return image;
}

ImageUInt8Ptr make_demo_image() {
    auto image = ImageUInt8::create(10, 10);
    for (int r = 0; r < image->getNumRows(); ++r) {
        for (int c = 0; c < image->getNumCols(); ++c) {
            const int p = r * image->getNumCols() + c;
            int value = 20;
            if (r >= 2 && r <= 7 && c >= 2 && c <= 7) {
                value = 80;
            }
            if (r >= 4 && r <= 5 && c >= 4 && c <= 5) {
                value = 130;
            }
            if ((r == 1 && c == 8) || (r == 8 && c == 1)) {
                value = 160;
            }
            (*image)[p] = static_cast<uint8_t>(value);
        }
    }
    return image;
}

ImageUInt8Ptr make_partial_f_order_independence_fixture_image() {
    auto image = ImageUInt8::create(5, 5);
    const std::vector<uint8_t> raw = {
        0, 0, 0, 0, 0,
        0, 20, 0, 10, 0,
        0, 10, 0, 0, 0,
        0, 0, 0, 0, 0,
        0, 0, 0, 0, 0,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_multiple_nodeca_p1_restricted_fixture_image() {
    auto image = ImageUInt8::create(2, 2);
    const std::vector<uint8_t> raw = {
        12, 10,
        11, 10,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_reduced_split_fixture_image() {
    auto image = ImageUInt8::create(6, 6);
    const std::vector<uint8_t> raw = {
        20, 44, 68, 44, 68, 44,
        56, 140, 104, 116, 104, 116,
        68, 140, 104, 116, 140, 128,
        56, 104, 56, 20, 32, 104,
        80, 104, 152, 104, 56, 128,
        68, 116, 140, 68, 20, 92,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_p1_mcpc_target_b_fixture_image() {
    auto image = ImageUInt8::create(6, 6);
    const std::vector<uint8_t> raw = {
        44, 116, 104, 68, 104, 44,
        104, 104, 128, 116, 92, 152,
        116, 80, 140, 92, 152, 92,
        56, 56, 44, 20, 104, 128,
        68, 92, 104, 68, 104, 128,
        80, 80, 128, 128, 116, 152,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_p1_mcpc_connector_fixture_image() {
    auto image = ImageUInt8::create(9, 9);
    const std::vector<uint8_t> raw = {
        80, 116, 80, 80, 44, 56, 68, 56, 68,
        80, 80, 92, 104, 68, 92, 80, 140, 92,
        128, 20, 104, 128, 20, 116, 140, 20, 152,
        92, 44, 44, 44, 32, 128, 44, 92, 56,
        56, 68, 68, 56, 68, 152, 80, 56, 68,
        140, 128, 128, 68, 140, 128, 104, 152, 128,
        104, 116, 116, 20, 104, 44, 44, 56, 104,
        80, 68, 104, 92, 140, 104, 116, 92, 152,
        128, 128, 152, 116, 140, 104, 140, 92, 152,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_p1_mcpc_structure_mismatch_fixture_image() {
    auto image = ImageUInt8::create(9, 9);
    const std::vector<uint8_t> raw = {
        104, 68, 68, 80, 92, 68, 32, 32, 116,
        104, 56, 152, 56, 116, 68, 128, 92, 140,
        116, 152, 140, 116, 20, 140, 80, 80, 104,
        56, 44, 104, 68, 32, 116, 56, 32, 56,
        56, 104, 128, 128, 80, 104, 128, 68, 140,
        104, 116, 116, 152, 152, 116, 92, 80, 152,
        56, 56, 128, 56, 116, 92, 56, 56, 92,
        44, 80, 68, 44, 116, 92, 44, 68, 104,
        92, 80, 116, 68, 104, 152, 80, 104, 20,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_icdar_connector_crop_fixture_image() {
    auto image = ImageUInt8::create(10, 8);
    const std::vector<uint8_t> raw = {
        53, 32, 6, 6, 6, 5, 5, 5,
        53, 37, 7, 6, 5, 5, 5, 5,
        53, 36, 5, 7, 6, 6, 6, 6,
        53, 34, 5, 7, 8, 7, 7, 7,
        55, 34, 5, 8, 8, 7, 7, 7,
        55, 29, 4, 8, 8, 7, 7, 7,
        56, 27, 2, 8, 8, 7, 7, 7,
        59, 28, 2, 7, 7, 7, 7, 7,
        64, 30, 2, 7, 7, 7, 7, 7,
        68, 30, 2, 7, 7, 7, 7, 7,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

ImageUInt8Ptr make_icdar_phase_candidate_crop_fixture_image() {
    auto image = ImageUInt8::create(7, 7);
    const std::vector<uint8_t> raw = {
        148, 182, 206, 200, 182, 202, 205,
        173, 223, 227, 227, 220, 202, 215,
        191, 227, 227, 227, 229, 220, 216,
        204, 220, 227, 227, 229, 229, 217,
        232, 220, 223, 229, 229, 229, 206,
        239, 233, 209, 209, 228, 212, 167,
        239, 239, 230, 183, 167, 149, 148,
    };
    for (std::size_t i = 0; i < raw.size(); ++i) {
        (*image)[static_cast<int>(i)] = raw[i];
    }
    return image;
}

double direct_image_variance(ImageUInt8Ptr image) {
    const int n = image->getNumRows() * image->getNumCols();
    double sum = 0.0;
    double sumSquares = 0.0;
    for (int p = 0; p < n; ++p) {
        const double value = static_cast<double>((*image)[p]);
        sum += value;
        sumSquares += value * value;
    }
    const double mean = sum / static_cast<double>(n);
    return (sumSquares / static_cast<double>(n)) - (mean * mean);
}

void test_root_gray_variance_matches_image_variance() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image, true, adjacency);

    const auto variance = noninc::computeGrayLevelVariance(maxTree);
    const double expected = direct_image_variance(image);
    const double observed = static_cast<double>(variance[static_cast<std::size_t>(maxTree.getRoot())]);

    require(std::abs(expected - observed) < 1.0e-3,
            "root gray-level variance must match direct image variance");
}

void test_root_box_area_perimeter_ratio_matches_direct_geometry() {
    auto image = ImageUInt8::create(2, 3);
    for (int p = 0; p < image->getNumRows() * image->getNumCols(); ++p) {
        (*image)[p] = 42;
    }

    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image, true, adjacency);

    const auto ratio = noninc::computeBoxAreaPerimeterRatio(maxTree);
    const double expected = 6.0 / (2.0 * 3.0 + 2.0 * 2.0);
    const double observed = static_cast<double>(ratio[static_cast<std::size_t>(maxTree.getRoot())]);

    require(std::abs(expected - observed) < 1.0e-6,
            "root area-per-box-perimeter ratio must match direct bounding-box geometry");
}

void test_node_storage_reallocation_counters_track_capacity_growth() {
    auto image = ImageUInt8::create(1, 1);
    (*image)[0] = 42;
    auto adjacency = std::make_shared<AdjacencyRelation>(1, 1, 1.0);
    DynamicComponentTree tree(image, true, adjacency);

    require(tree.getNodeStorageCapacityGrowthEvents() == 0,
            "freshly built tree should not report dynamic capacity growth");
    require(tree.getNodeStorageVectorReallocationEvents() == 0,
            "freshly built tree should not report node-vector reallocations");

    for (int i = 0; i < 1024 && tree.getNodeStorageCapacityGrowthEvents() == 0; ++i) {
        const NodeId nodeId = tree.createDetachedNode(42);
        require(nodeId != InvalidNode,
                "detached node creation should succeed in the synthetic reallocation test");
    }

    require(tree.getNodeStorageCapacityGrowthEvents() > 0,
            "capacity growth counter must increase when dynamic node slots exceed capacity");
    require(tree.getNodeStorageVectorReallocationEvents() > 0,
            "vector reallocation counter must increase when node storage capacity grows");
    require(tree.getNodeStorageReallocatedElements() > 0,
            "reallocated element estimate must be positive after capacity growth");
    require(tree.getNodeStorageReallocatedBytes() > 0,
            "reallocated byte estimate must be positive after capacity growth");

    tree.resetNodeStorageReallocationCounters();
    require(tree.getNodeStorageCapacityGrowthEvents() == 0 &&
                tree.getNodeStorageVectorReallocationEvents() == 0 &&
                tree.getNodeStorageReallocatedElements() == 0 &&
                tree.getNodeStorageReallocatedBytes() == 0,
            "reset must clear all node-storage reallocation counters");
}

void test_box_area_perimeter_clean_update_matches_naive_sequence() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 64, 64);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    const std::vector<float> thresholds = {0.25f, 0.5f, 1.0f, 2.0f, 4.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCleanUnchecked,
        noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio);

    auto updated = runner.filter(thresholds);
    auto naive = noninc::applyNaiveAttributeMergeSequence(
        image,
        adjacency,
        thresholds,
        noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio);

    require(updated->isEqual(naive),
            "area-per-box-perimeter clean local update must match the naive sequence");
}

void require_live_box_area_perimeter_attributes_match_full_recompute(
    DynamicComponentTree &tree,
    const std::vector<float> &actual,
    const std::string &treeName) {
    const auto expected = noninc::computeBoxAreaPerimeterRatio(tree);
    require(actual.size() >= expected.size(),
            treeName + ": local attribute buffer must cover the current node slots");
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        const auto index = static_cast<std::size_t>(nodeId);
        const float delta = std::abs(actual[index] - expected[index]);
        require(delta < 1.0e-5f,
                treeName + ": local bbox attribute diverged at live node " +
                    std::to_string(nodeId));
    }
}

void test_box_area_perimeter_clean_update_keeps_local_attribute_buffers_exact() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 64, 64);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    const std::vector<float> thresholds = {0.25f, 0.5f, 1.0f, 2.0f, 4.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCleanUnchecked,
        noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio);

    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require_live_box_area_perimeter_attributes_match_full_recompute(
            runner.maxTree(),
            runner.maxAttribute(),
            "max-tree");
        require_live_box_area_perimeter_attributes_match_full_recompute(
            runner.minTree(),
            runner.minAttribute(),
            "min-tree");
    }

    const auto naive = noninc::applyNaiveAttributeMergeSequence(
        image,
        adjacency,
        thresholds,
        noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio);
    require(runner.reconstructImage()->isEqual(naive),
            "bbox local attribute regression must still match the naive sequence");
}

struct FullSupportCandidateOrderEntry {
    NodeId nodeId = InvalidNode;
    int altitude = 0;
    std::vector<PixelId> supportPixels;
};

std::vector<NodeId> collect_reference_candidates_by_full_support_order(
    const DynamicComponentTree &tree,
    const std::vector<float> &attribute,
    float threshold) {
    std::vector<NodeId> queue;
    std::vector<FullSupportCandidateOrderEntry> ordered;
    queue.push_back(tree.getRoot());
    std::size_t head = 0;
    while (head < queue.size()) {
        const NodeId nodeId = queue[head++];
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        if (nodeId != tree.getRoot() &&
            static_cast<std::size_t>(nodeId) < attribute.size() &&
            attribute[static_cast<std::size_t>(nodeId)] <= threshold) {
            ordered.push_back({nodeId, tree.getAltitude(nodeId), {}});
            noninc::collectNodeSupportPixelsInto(
                tree,
                nodeId,
                ordered.back().supportPixels);
        }
        for (NodeId childId : tree.getChildren(nodeId)) {
            queue.push_back(childId);
        }
    }
    std::sort(
        ordered.begin(),
        ordered.end(),
        [](const FullSupportCandidateOrderEntry &lhs,
           const FullSupportCandidateOrderEntry &rhs) {
            return std::tie(lhs.altitude, lhs.supportPixels, lhs.nodeId) <
                   std::tie(rhs.altitude, rhs.supportPixels, rhs.nodeId);
        });

    std::vector<NodeId> output;
    output.reserve(ordered.size());
    for (const auto &entry : ordered) {
        output.push_back(entry.nodeId);
    }
    return output;
}

void require_candidate_order_matches_full_support_reference(
    const DynamicComponentTree &tree,
    const std::vector<float> &attribute,
    float threshold,
    const std::string &treeName) {
    const auto observed =
        noninc::collectNonIncreasingMergeCandidates(tree, attribute, threshold);
    const auto expected =
        collect_reference_candidates_by_full_support_order(tree, attribute, threshold);
    require(observed == expected,
            treeName + ": cheap candidate key must preserve full-support order");
}

void test_candidate_key_matches_full_support_order_on_clean_bbox_path() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 64, 64);
    const std::vector<float> thresholds = {0.25f, 0.5f, 1.0f, 2.0f, 4.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCleanUnchecked,
        noninc::NonIncreasingAttributeKind::BoxAreaPerimeterRatio);

    for (float threshold : thresholds) {
        require_candidate_order_matches_full_support_reference(
            runner.maxTree(),
            runner.maxAttribute(),
            threshold,
            "max-tree");
        require_candidate_order_matches_full_support_reference(
            runner.minTree(),
            runner.minAttribute(),
            threshold,
            "min-tree");
        runner.applyThreshold(threshold);
    }
}

void test_partial_f_default_records_formal_acceptance_without_fallback() {
    auto image = make_demo_image();
    const std::vector<float> thresholds = {0.0f, 10.0f, 200.0f};

    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    int attempts = 0;
    int accepted = 0;
    int rejected = 0;
    int fallbackEvents = 0;
    int copiedEvents = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        for (const auto &event : runner.lastMergeEvents()) {
            if (!event.partialFAttempted) {
                continue;
            }
            ++attempts;
            if (event.partialFAccepted) {
                ++accepted;
                require(event.applied,
                        "accepted partial F events must be applied");
                require(event.partialFMutationAttempted,
                        "accepted partial F events must attempt the direct mutation");
                require(event.mergeNodeLocalDualMatchesPrimalAfterMerge,
                        "accepted partial F events must preserve primal/dual reconstruction equality");
            }
            if (event.partialFRejected) {
                ++rejected;
                require(!event.applied,
                        "rejected partial F events must not be marked applied");
                require(!event.partialFMutationAttempted,
                        "rejected partial F events must not attempt mutation");
                require(!event.partialFRejectReasonCode.empty(),
                        "rejected partial F events must expose a formal reason");
            }
            if (event.mergeNodeLocalFallbackToPhaseRebuild ||
                event.strippedLocalFallbackToOracle ||
                event.strippedLocalFallbackToRebuild ||
                event.strippedLocalFallbackToPhaseRebuild) {
                ++fallbackEvents;
            }
            if (event.mergeNodeLocalCopyMicros > 0) {
                ++copiedEvents;
            }
        }
    }

    require(attempts > 0,
            "default non-increasing CASF must exercise partial F attempts");
    require(accepted > 0,
            "default partial F run must accept at least one event on the demo fixture");
    require(attempts == accepted + rejected,
            "each partial F attempt must be accepted or rejected");
    require(fallbackEvents == 0,
            "partial F must not use fallback paths");
    require(copiedEvents == 0,
            "partial F must not use copy-based decision paths");
}

void test_partial_f_rejects_invalid_node_without_mutation() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adjacency);

    const auto beforeMax = maxTree.reconstructionImage();
    const auto beforeMin = minTree.reconstructionImage();
    const int beforeMaxAlive = maxTree.getNumNodes();
    const int beforeMinAlive = minTree.getNumNodes();

    const auto result = noninc::applyNonIncreasingCasfPartialF(
        maxTree,
        minTree,
        adjust,
        *adjacency,
        maxTree.getRoot());

    require(result.rejected,
            "partial F must reject root merge events");
    require(result.rejectReasonCode == "invalid_event",
            "partial F root rejection must use invalid_event");
    require(!result.mutationAttempted,
            "partial F invalid-event rejection must not attempt mutation");
    require(maxTree.reconstructionImage()->isEqual(beforeMax),
            "partial F invalid-event rejection must preserve the primal reconstruction");
    require(minTree.reconstructionImage()->isEqual(beforeMin),
            "partial F invalid-event rejection must preserve the dual reconstruction");
    require(maxTree.getNumNodes() == beforeMaxAlive,
            "partial F invalid-event rejection must preserve primal live node count");
    require(minTree.getNumNodes() == beforeMinAlive,
            "partial F invalid-event rejection must preserve dual live node count");
}

void test_partial_f_mcpc_formal_accounting_without_fallback() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    bool foundMcpcRejection = false;
    int partialFAttempts = 0;
    int partialFAcceptances = 0;
    int partialFRejections = 0;
    int partialFMutationAttempts = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "partial F MCPC rejection fixture must preserve primal/dual reconstruction equality");
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.partialFAttempted) {
                ++partialFAttempts;
            }
            if (event.partialFAccepted) {
                ++partialFAcceptances;
            }
            if (event.partialFRejected) {
                ++partialFRejections;
                require(!event.partialFMutationAttempted,
                        "partial F rejection must not attempt mutation");
            }
            if (event.partialFMutationAttempted) {
                ++partialFMutationAttempts;
            }
            if (event.partialFRejected &&
                event.partialFRejectReasonCode.rfind("d12_", 0) == 0 &&
                event.mergeNodeLocalPartialOwnerComponents > 0) {
                foundMcpcRejection = true;
                require(event.partialFAttempted,
                        "partial F MCPC rejection must be recorded as a partial F attempt");
                require(!event.partialFAccepted,
                        "partial F MCPC rejection must not be recorded as accepted");
                require(event.mergeNodeLocalAttempted,
                        "partial F MCPC rejection must still expose merge-node diagnostics");
                require(!event.mergeNodeLocalAccepted,
                        "partial F MCPC rejection must not accept the merge-node update");
                require(!event.mergeNodeLocalApplied,
                        "partial F MCPC rejection must not apply the merge-node update");
                require(event.mergeNodeLocalRejectedComponents > 0,
                        "partial F MCPC rejection must expose rejected component evidence");
                require(!event.partialFMutationAttempted,
                        "partial F must reject MCPC events before mutation");
                require(!event.applied,
                        "partial F MCPC rejection must leave the event unapplied");
                require(!event.mergeNodeLocalFallbackToPhaseRebuild,
                        "partial F MCPC rejection must not fallback to phase rebuild");
                require(!event.strippedLocalAttempted,
                        "partial F MCPC rejection must not try stripped copied recomposition");
                require(!event.mcpcDirectPreflightAttempted,
                        "partial F MCPC rejection must not enter MCPC direct preflight");
                require(!event.mcpcNoCopyCandidateEvaluated,
                        "partial F MCPC rejection must not evaluate no-copy candidate diagnostics");
                require(!event.mcpcDirectCommitCandidateEvaluated,
                        "partial F MCPC rejection must not evaluate direct-commit diagnostics");
                require(!event.mcpcNoCopyDirectCommitAttempted,
                        "partial F MCPC rejection must not attempt no-copy direct commit");
                require(!event.mcpcValidatedDirectCommitAttempted,
                        "partial F MCPC rejection must not attempt validated direct commit");
                break;
            }
        }
        if (foundMcpcRejection) {
            break;
        }
    }

    require(partialFAttempts > 0,
            "5.2.10 fixture must exercise partial F attempts");
    require(partialFAttempts == partialFAcceptances + partialFRejections,
            "partial F MCPC accounting must partition attempts into acceptances and rejections");
    require(partialFMutationAttempts == partialFAcceptances,
            "partial F mutations must coincide exactly with accepted events");
    (void)foundMcpcRejection;
}

struct PartialFMcpcLiveFixture {
    bool found = false;
    std::size_t eventIndex = 0;
    noninc::NonIncreasingMergeEvent event;
    noninc::NonIncreasingMcpcLiveMutationDryRun dryRun;
};

PartialFMcpcLiveFixture find_ruler128_first_partial_f_mcpc_live_fixture() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 128, 128);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.applyThreshold(0.0f);

    auto &maxTree = runner.maxTree();
    auto &minTree = runner.minTree();
    auto maxAttribute = noninc::computeGrayLevelVariance(maxTree);
    auto maxCandidates = noninc::collectNonIncreasingMergeCandidates(
        maxTree,
        maxAttribute,
        1.0f);
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&minTree, &maxTree, *adjacency);

    PartialFMcpcLiveFixture out;
    std::size_t eventIndex = 0;
    for (NodeId nodeId : maxCandidates) {
        if (nodeId == InvalidNode ||
            !maxTree.isNode(nodeId) ||
            !maxTree.isAlive(nodeId) ||
            static_cast<std::size_t>(nodeId) >= maxAttribute.size()) {
            continue;
        }

        auto dryRun = noninc::NonIncreasingMcpcLiveMutationDryRun{};
        const auto partialFPreflight = noninc::preflightNonIncreasingCasfPartialF(
            maxTree,
            minTree,
            *adjacency,
            nodeId);
        if (partialFPreflight.partialOwnerComponents > 0) {
            const auto artifacts = noninc::buildMcpcDirectPreflightArtifacts(
                maxTree,
                minTree,
                adjacency,
                nodeId);
            dryRun = noninc::diagnoseMcpcLiveMutationDryRun(minTree, artifacts);
        }

        noninc::NonIncreasingMergeEvent event;
        event.threshold = 1.0f;
        event.phase = noninc::NonIncreasingPhase::MaxToMin;
        event.attributeValue = maxAttribute[static_cast<std::size_t>(nodeId)];
        noninc::mergeNodeIntoParentWithExplicitDualUpdate(
            maxTree,
            minTree,
            adjust,
            adjacency,
            noninc::NonIncreasingDualMaintenance::PartialF,
            nodeId,
            &event);

        if (dryRun.accepted) {
            require(event.partialFMcpcLiveUpdate,
                    "accepted direct dry-run must correspond to a live MCPC partial F event");
            out.found = true;
            out.eventIndex = eventIndex;
            out.event = event;
            out.dryRun = dryRun;
            return out;
        }
        ++eventIndex;
    }

    return out;
}

void test_partial_f_accepts_mcpc_live_subdomain_without_fallback() {
    const auto fixture = find_ruler128_first_partial_f_mcpc_live_fixture();
    require(fixture.found,
            "ruler crop must expose at least one accepted MCPC live subdomain event");

    const auto &event = fixture.event;
    require(event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.threshold == 1.0f,
            "ruler MCPC live fixture must remain in the max->min threshold-1 phase");
    require(event.componentCount > 0 && event.totalComponentPixels > 0,
            "ruler MCPC live fixture must remain a non-empty event");
    require(event.directProperParts > 0,
            "ruler MCPC live fixture must preserve direct proper parts");
    require(event.partialFAccepted && !event.partialFRejected,
            "MCPC live event must be accepted by partial F only");
    require(event.partialFMutationAttempted,
            "MCPC live event must mutate only after accepted preflight");
    require(event.applied && event.mergeNodeLocalAccepted && event.mergeNodeLocalApplied,
            "MCPC live event must apply the formal merge-node update");
    require(event.mergeNodeLocalPartialOwnerComponents > 0,
            "MCPC live event must exercise the partial-owner MCPC path");
    require(event.mergeNodeLocalComponentCount > 0 &&
                event.mergeNodeLocalPixelCount > 0,
            "MCPC live event must report a non-empty local update");
    require(event.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "MCPC live event must preserve primal/dual reconstruction equality");
    require(event.partialFRejectReasonCode.empty(),
            "MCPC live acceptance must not carry a reject reason");
    require(!event.mergeNodeLocalFallbackToPhaseRebuild &&
                !event.strippedLocalFallbackToOracle &&
                !event.strippedLocalFallbackToRebuild &&
                !event.strippedLocalFallbackToPhaseRebuild,
            "MCPC live acceptance must not use fallback paths");
    require(!event.strippedLocalAttempted,
            "MCPC live acceptance must not use stripped copied recomposition");
    require(event.mergeNodeLocalCopyMicros == 0,
            "MCPC live acceptance must not use copy-based decision paths");

    const auto &dryRun = fixture.dryRun;
    require(dryRun.evaluated && dryRun.accepted,
            "direct dry-run must accept the ruler MCPC live fixture");
    require(dryRun.rejectReasonCode.empty(),
            "accepted direct dry-run must not expose a reject reason");
    require(dryRun.predictedAuxiliaryNodeCount ==
                dryRun.predictedHiddenExternalParentConnectorNodeCount +
                    dryRun.predictedStableChildBridgeNodeCount,
            "accepted direct dry-run must account for all classified auxiliary nodes");
    require(dryRun.plannedStableChildRefCount == dryRun.stableChildToMoveCount,
            "accepted direct dry-run must account for stable-child moves");
}

void test_partial_f_accepts_parent_link_bridge_live_subdomain_without_fallback() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 128, 128);
    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    runner.applyThreshold(0.0f);
    runner.applyThreshold(1.0f);

    const auto &events = runner.lastMergeEvents();
    const noninc::NonIncreasingMergeEvent *selected = nullptr;
    std::size_t eventIndex = 0;
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcParentLinkBridgeLiveUpdate) {
            selected = &events[i];
            eventIndex = i;
            break;
        }
    }

    require(selected != nullptr,
            "ruler crop must expose a parent-link bridge MCPC live partial F event");

    const auto &event = *selected;
    require(event.threshold == 1.0f,
            "parent-link bridge fixture must remain in the threshold-1 phase");
    require(eventIndex < events.size(),
            "parent-link bridge fixture must report a valid event index");
    require(event.componentCount > 0 && event.totalComponentPixels > 0,
            "parent-link bridge fixture must remain a non-empty event");
    require(event.directProperParts > 0,
            "parent-link bridge fixture must preserve direct proper parts");
    require(event.partialFAccepted && !event.partialFRejected,
            "parent-link bridge event must be accepted by partial F only");
    require(event.partialFMutationAttempted &&
                event.partialFMcpcLiveUpdate &&
                event.partialFMcpcParentLinkBridgeLiveUpdate,
            "parent-link bridge event must be a live MCPC partial F mutation");
    require(event.applied && event.mergeNodeLocalAccepted && event.mergeNodeLocalApplied,
            "parent-link bridge event must apply the formal merge-node update");
    require(event.mergeNodeLocalPartialOwnerComponents > 0,
            "parent-link bridge event must exercise the partial-owner MCPC path");
    require(event.mergeNodeLocalDualMatchesPrimalBeforeMerge &&
                event.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "parent-link bridge event must preserve primal/dual reconstruction equality");
    require(event.partialFRejectReasonCode.empty(),
            "parent-link bridge acceptance must not carry a reject reason");
    require(!event.mergeNodeLocalFallbackToPhaseRebuild &&
                !event.strippedLocalFallbackToOracle &&
                !event.strippedLocalFallbackToRebuild &&
                !event.strippedLocalFallbackToPhaseRebuild,
            "parent-link bridge acceptance must not use fallback paths");
    require(!event.strippedLocalAttempted,
            "parent-link bridge acceptance must not use stripped copied recomposition");
    require(event.mergeNodeLocalCopyMicros == 0,
            "parent-link bridge acceptance must not use copy-based decision paths");
}

void test_partial_f_accepts_hidden_support_redundant_live_subdomain_without_fallback() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 128, 128);
    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    runner.applyThreshold(0.0f);
    runner.applyThreshold(1.0f);

    const auto &events = runner.lastMergeEvents();
    const noninc::NonIncreasingMergeEvent *selected = nullptr;
    for (const auto &event : events) {
        if (event.partialFMcpcHiddenSupportRedundantLiveUpdate) {
            selected = &event;
            break;
        }
    }

    require(selected != nullptr,
            "ruler crop must expose a redundant hidden-support MCPC live partial F event");

    const auto &event = *selected;
    require(event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.threshold == 1.0f,
            "hidden-support fixture must remain in the max->min threshold-1 phase");
    require(event.componentCount > 0 && event.totalComponentPixels > 0,
            "hidden-support fixture must remain a non-empty event");
    require(event.directProperParts > 0,
            "hidden-support fixture must preserve direct proper parts");
    require(event.partialFAccepted && !event.partialFRejected,
            "hidden-support event must be accepted by partial F only");
    require(event.partialFMutationAttempted &&
                event.partialFMcpcLiveUpdate &&
                event.partialFMcpcHiddenSupportRedundantLiveUpdate,
            "hidden-support event must be a live MCPC partial F mutation");
    require(event.applied && event.mergeNodeLocalAccepted && event.mergeNodeLocalApplied,
            "hidden-support event must apply the formal merge-node update");
    require(event.mergeNodeLocalPartialOwnerComponents > 0,
            "hidden-support event must exercise the partial-owner MCPC path");
    require(event.mergeNodeLocalDualMatchesPrimalBeforeMerge &&
                event.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "hidden-support event must preserve primal/dual reconstruction equality");
    require(event.partialFRejectReasonCode.empty(),
            "hidden-support acceptance must not carry a reject reason");
    require(!event.mergeNodeLocalFallbackToPhaseRebuild &&
                !event.strippedLocalFallbackToOracle &&
                !event.strippedLocalFallbackToRebuild &&
                !event.strippedLocalFallbackToPhaseRebuild,
            "hidden-support acceptance must not use fallback paths");
    require(!event.strippedLocalAttempted,
            "hidden-support acceptance must not use stripped copied recomposition");
    require(event.mergeNodeLocalCopyMicros == 0,
            "hidden-support acceptance must not use copy-based decision paths");
}

void test_partial_f_omits_stable_child_bridge_before_live_mutation() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 128, 128);
    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    runner.applyThreshold(0.0f);
    runner.applyThreshold(1.0f);

    const auto &events = runner.lastMergeEvents();
    const noninc::NonIncreasingMergeEvent *selected = nullptr;
    bool rejectedByAuxiliaryBridge = false;
    for (const auto &event : events) {
        if (event.partialFRejected &&
            event.partialFRejectReasonCode == "d12_auxiliary_nodes_required" &&
            event.mergeNodeLocalPartialOwnerComponents > 0) {
            rejectedByAuxiliaryBridge = true;
        }
        if (event.partialFAccepted &&
            event.partialFMcpcLiveUpdate &&
            event.partialFMcpcStableChildMoveLiveUpdate &&
            event.mergeNodeLocalPartialOwnerComponents > 0) {
            selected = &event;
        }
    }

    require(!rejectedByAuxiliaryBridge,
            "stable-child bridge cases must no longer reject as auxiliary-node events");
    require(selected != nullptr,
            "ruler crop must expose a stable-child MCPC event accepted without bridge materialization");

    const auto &event = *selected;
    require(event.threshold == 1.0f,
            "stable-child bridge fixture must remain in the threshold-1 phase");
    require(event.componentCount > 0 && event.totalComponentPixels > 0,
            "stable-child bridge fixture must remain a non-empty event");
    require(event.directProperParts > 0,
            "stable-child bridge fixture must preserve direct proper parts");
    require(event.partialFAccepted && !event.partialFRejected,
            "stable-child bridge event must be accepted by partial F after canonical construction");
    require(event.partialFMutationAttempted &&
                event.partialFMcpcLiveUpdate &&
                event.partialFMcpcStableChildMoveLiveUpdate &&
                !event.partialFMcpcStableChildBridgeLiveUpdate,
            "stable-child bridge event must be represented as a canonical stable-child move");
    require(event.applied && event.mergeNodeLocalAccepted && event.mergeNodeLocalApplied,
            "stable-child bridge event must apply the formal merge-node update");
    require(event.mergeNodeLocalPartialOwnerComponents > 0,
            "stable-child bridge event must exercise the partial-owner MCPC path");
    require(event.mergeNodeLocalDualMatchesPrimalBeforeMerge &&
                event.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "stable-child bridge event must preserve primal/dual reconstruction equality");
    require(event.partialFRejectReasonCode.empty(),
            "stable-child bridge acceptance must not carry a rejection code");
    require(!event.mergeNodeLocalFallbackToPhaseRebuild &&
                !event.strippedLocalFallbackToOracle &&
                !event.strippedLocalFallbackToRebuild &&
                !event.strippedLocalFallbackToPhaseRebuild,
            "stable-child bridge acceptance must not use fallback paths");
    require(!event.strippedLocalAttempted,
            "stable-child bridge acceptance must not use stripped copied recomposition");
    require(event.mergeNodeLocalCopyMicros == 0,
            "stable-child bridge acceptance must not use copy-based decision paths");
}

void test_partial_f_accepts_stable_child_move_live_subdomain_without_fallback() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 128, 128);
    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    runner.applyThreshold(0.0f);
    runner.applyThreshold(1.0f);

    const auto &events = runner.lastMergeEvents();
    const noninc::NonIncreasingMergeEvent *selected = nullptr;
    for (std::size_t i = 0; i < events.size(); ++i) {
        if (events[i].partialFMcpcStableChildMoveLiveUpdate) {
            selected = &events[i];
            break;
        }
    }

    require(selected != nullptr,
            "ruler crop must expose a stable-child move MCPC live partial F event");

    const auto &event = *selected;
    require(event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.threshold == 1.0f,
            "stable-child move fixture must remain in the max->min threshold-1 phase");
    require(event.nodeId != InvalidNode && event.parentId != InvalidNode,
            "stable-child move fixture must expose live node ids");
    require(event.componentCount > 0 && event.totalComponentPixels > 0,
            "stable-child move fixture must remain a non-empty event");
    require(event.directProperParts > 0,
            "stable-child move fixture must preserve direct proper parts");
    require(event.partialFAccepted && !event.partialFRejected,
            "stable-child move event must be accepted by partial F only");
    require(event.partialFMutationAttempted &&
                event.partialFMcpcLiveUpdate &&
                event.partialFMcpcStableChildMoveLiveUpdate,
            "stable-child move event must be a live MCPC partial F mutation");
    require(event.applied && event.mergeNodeLocalAccepted && event.mergeNodeLocalApplied,
            "stable-child move event must apply the formal merge-node update");
    require(event.mergeNodeLocalPartialOwnerComponents > 0,
            "stable-child move event must exercise the partial-owner MCPC path");
    require(event.mergeNodeLocalDualMatchesPrimalBeforeMerge &&
                event.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "stable-child move event must preserve primal/dual reconstruction equality");
    require(event.partialFRejectReasonCode.empty(),
            "stable-child move acceptance must not carry a reject reason");
    require(!event.mergeNodeLocalFallbackToPhaseRebuild &&
                !event.strippedLocalFallbackToOracle &&
                !event.strippedLocalFallbackToRebuild &&
                !event.strippedLocalFallbackToPhaseRebuild,
            "stable-child move acceptance must not use fallback paths");
    require(!event.strippedLocalAttempted,
            "stable-child move acceptance must not use stripped copied recomposition");
    require(event.mergeNodeLocalCopyMicros == 0,
            "stable-child move acceptance must not use copy-based decision paths");
}

void test_partial_f_accepts_empty_owner_release_live_subdomain_without_fallback() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 128, 128);
    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    runner.applyThreshold(0.0f);
    runner.applyThreshold(1.0f);

    const auto &events = runner.lastMergeEvents();
    const noninc::NonIncreasingMergeEvent *selected = nullptr;
    for (const auto &event : events) {
        if (event.partialFMcpcEmptyOwnerReleaseLiveUpdate) {
            selected = &event;
            break;
        }
    }

    require(selected != nullptr,
            "ruler crop must expose an empty-owner release MCPC live partial F event");

    const auto &event = *selected;
    require(event.threshold == 1.0f,
            "empty-owner release fixture must remain in the threshold-1 phase");
    require(event.componentCount > 0 && event.totalComponentPixels > 0,
            "empty-owner release fixture must remain a non-empty event");
    require(event.directProperParts > 0,
            "empty-owner release fixture must preserve direct proper parts");
    require(event.partialFAccepted && !event.partialFRejected,
            "empty-owner release event must be accepted by partial F only");
    require(event.partialFMutationAttempted &&
                event.partialFMcpcLiveUpdate &&
                event.partialFMcpcEmptyOwnerReleaseLiveUpdate,
            "empty-owner release event must be a live MCPC partial F mutation");
    require(event.applied && event.mergeNodeLocalAccepted && event.mergeNodeLocalApplied,
            "empty-owner release event must apply the formal merge-node update");
    require(event.mergeNodeLocalPartialOwnerComponents > 0,
            "empty-owner release event must exercise the partial-owner MCPC path");
    require(event.mergeNodeLocalDualMatchesPrimalBeforeMerge &&
                event.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "empty-owner release event must preserve primal/dual reconstruction equality");
    require(event.partialFRejectReasonCode.empty(),
            "empty-owner release acceptance must not carry a reject reason");
    require(!event.mergeNodeLocalFallbackToPhaseRebuild &&
                !event.strippedLocalFallbackToOracle &&
                !event.strippedLocalFallbackToRebuild &&
                !event.strippedLocalFallbackToPhaseRebuild,
            "empty-owner release acceptance must not use fallback paths");
    require(!event.strippedLocalAttempted,
            "empty-owner release acceptance must not use stripped copied recomposition");
    require(event.mergeNodeLocalCopyMicros == 0,
            "empty-owner release acceptance must not use copy-based decision paths");
}

void test_partial_f_accepts_materializable_missing_target_b() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::RebuildAfterPhase);
    runner.applyThreshold(0.0f);

    DynamicComponentTree &primalTree = runner.maxTree();
    DynamicComponentTree &dualTree = runner.minTree();
    const NodeId nodeId = 48;
    require(primalTree.isNode(nodeId) && primalTree.isAlive(nodeId),
            "ruler fixture must expose the recorded missing target-b node before threshold 1");

    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(&dualTree, &primalTree, *adjacency);
    const auto preflight = noninc::preflightNonIncreasingCasfPartialF(
        primalTree,
        dualTree,
        *adjacency,
        nodeId);
    require(preflight.accepted,
            "partial F preflight must accept the materializable missing target-b event");
    require(preflight.missingTargetComponents > 0,
            "partial F preflight must record missing target-b materialization");

    const auto result = noninc::applyNonIncreasingCasfPartialF(
        primalTree,
        dualTree,
        adjust,
        *adjacency,
        nodeId);
    require(result.accepted,
            "partial F must accept materializable missing target-b events");
    require(result.mutationAttempted,
            "partial F accepted missing-target events must mutate directly");
    require(result.localUpdateApplied && result.primalMergeApplied,
            "partial F accepted missing-target events must update dual and primal");
    require(dualTree.reconstructionImage()->isEqual(primalTree.reconstructionImage()),
            "partial F accepted missing-target events must preserve reconstruction equality");
}

NodeId find_order_independence_partial_f_node(const DynamicComponentTree &primalTree,
                                              const DynamicComponentTree &dualTree,
                                              AdjacencyRelation &adjacency) {
    for (NodeId nodeId = 0; nodeId < primalTree.getNumInternalNodeSlots(); ++nodeId) {
        if (!primalTree.isNode(nodeId) ||
            !primalTree.isAlive(nodeId) ||
            primalTree.isRoot(nodeId)) {
            continue;
        }
        const auto preflight = noninc::preflightNonIncreasingCasfPartialF(
            primalTree,
            dualTree,
            adjacency,
            nodeId);
        if (preflight.accepted && preflight.componentCount > 1) {
            return nodeId;
        }
    }
    return InvalidNode;
}

struct OrderedComponentApplyResult {
    DynamicComponentTree primalTree;
    DynamicComponentTree dualTree;
    int componentCount = 0;
};

struct PartialFOrderIndependenceTarget {
    bool primalIsMaxTree = true;
    NodeId nodeId = InvalidNode;
};

PartialFOrderIndependenceTarget find_order_independence_partial_f_target(
    const DynamicComponentTree &maxTree,
    const DynamicComponentTree &minTree,
    AdjacencyRelation &adjacency) {
    const auto findInPair = [&](const DynamicComponentTree &primalTree,
                                const DynamicComponentTree &dualTree) {
        return find_order_independence_partial_f_node(
            primalTree,
            dualTree,
            adjacency);
    };

    if (const NodeId maxNode = findInPair(maxTree, minTree); maxNode != InvalidNode) {
        return {true, maxNode};
    }
    if (const NodeId minNode = findInPair(minTree, maxTree); minNode != InvalidNode) {
        return {false, minNode};
    }
    return {};
}

OrderedComponentApplyResult apply_partial_f_components_in_order(
    const ImageUInt8Ptr &image,
    const AdjacencyRelationPtr &adjacency,
    bool reverseComponentOrder) {
    DynamicComponentTree maxTree(image->clone(), true, adjacency);
    DynamicComponentTree minTree(image->clone(), false, adjacency);
    const auto target = find_order_independence_partial_f_target(
        maxTree,
        minTree,
        *adjacency);
    require(target.nodeId != InvalidNode,
            "order-independence fixture must expose an accepted multi-component partial F node");

    OrderedComponentApplyResult out{
        target.primalIsMaxTree ? std::move(maxTree) : std::move(minTree),
        target.primalIsMaxTree ? std::move(minTree) : std::move(maxTree),
        0};
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust(
        out.primalTree.isMaxtree() ? &out.dualTree : &out.primalTree,
        out.primalTree.isMaxtree() ? &out.primalTree : &out.dualTree,
        *adjacency);
    const NodeId nodeId = target.nodeId;

    auto components = noninc::directProperPartComponents(
        out.primalTree,
        *adjacency,
        nodeId);
    require(components.size() > 1,
            "order-independence fixture must have multiple connected components");
    out.componentCount = static_cast<int>(components.size());
    if (reverseComponentOrder) {
        std::reverse(components.begin(), components.end());
    }

    const NodeId parentId = out.primalTree.getNodeParent(nodeId);
    const auto b = static_cast<AltitudeType>(out.primalTree.getAltitude(parentId));
    for (const auto &component : components) {
        adjust.updateTreeWithExplicitRegion(
            &out.dualTree,
            std::span<const PixelId>(component),
            b,
            true);
    }
    out.primalTree.mergeNodeIntoParent(nodeId);
    return out;
}

void test_partial_f_component_order_independence_modulo_isomorphism() {
    auto image = make_partial_f_order_independence_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);

    auto forward = apply_partial_f_components_in_order(image, adjacency, false);
    auto reverse = apply_partial_f_components_in_order(image, adjacency, true);

    require(forward.componentCount == reverse.componentCount && forward.componentCount > 1,
            "forward and reverse component schedules must process the same family CC_G(C)");
    require(forward.dualTree.reconstructionImage()->isEqual(forward.primalTree.reconstructionImage()),
            "forward component order must preserve primal/dual reconstruction equality");
    require(reverse.dualTree.reconstructionImage()->isEqual(reverse.primalTree.reconstructionImage()),
            "reverse component order must preserve primal/dual reconstruction equality");
    require(forward.primalTree.reconstructionImage()->isEqual(reverse.primalTree.reconstructionImage()),
            "component order must preserve the primal output image");
    require(forward.dualTree.reconstructionImage()->isEqual(reverse.dualTree.reconstructionImage()),
            "component order must preserve the dual output image");
    require(noninc::collectNonIncreasingTreeStructureSignatures(forward.primalTree) ==
                noninc::collectNonIncreasingTreeStructureSignatures(reverse.primalTree),
            "component order must preserve the primal tree modulo ids");
    require(noninc::collectNonIncreasingTreeStructureSignatures(forward.dualTree) ==
                noninc::collectNonIncreasingTreeStructureSignatures(reverse.dualTree),
            "component order must preserve the dual tree modulo ids");
}

struct PartialFReplayEvent {
    float threshold = 0.0f;
    noninc::NonIncreasingPhase phase = noninc::NonIncreasingPhase::MaxToMin;
    NodeId nodeId = InvalidNode;
    NodeId parentId = InvalidNode;
    int nodeLevel = 0;
    int parentLevel = 0;
    int componentCount = 0;
    int totalComponentPixels = 0;
    bool attempted = false;
    bool accepted = false;
    bool rejected = false;
    bool mutationAttempted = false;
    bool mcpcLiveUpdate = false;
    bool mcpcHiddenSupportRedundantLiveUpdate = false;
    bool mcpcParentLinkBridgeLiveUpdate = false;
    bool mcpcStableChildBridgeLiveUpdate = false;
    bool mcpcChildLinkBridgeLiveUpdate = false;
    bool mcpcEmptyOwnerReleaseLiveUpdate = false;
    bool mcpcStableChildMoveLiveUpdate = false;
    bool applied = false;
    std::string rejectReasonCode;

    bool operator==(const PartialFReplayEvent &other) const {
        return threshold == other.threshold &&
               phase == other.phase &&
               nodeId == other.nodeId &&
               parentId == other.parentId &&
               nodeLevel == other.nodeLevel &&
               parentLevel == other.parentLevel &&
               componentCount == other.componentCount &&
               totalComponentPixels == other.totalComponentPixels &&
               attempted == other.attempted &&
               accepted == other.accepted &&
               rejected == other.rejected &&
               mutationAttempted == other.mutationAttempted &&
               mcpcLiveUpdate == other.mcpcLiveUpdate &&
               mcpcHiddenSupportRedundantLiveUpdate ==
                   other.mcpcHiddenSupportRedundantLiveUpdate &&
               mcpcParentLinkBridgeLiveUpdate ==
                   other.mcpcParentLinkBridgeLiveUpdate &&
               mcpcStableChildBridgeLiveUpdate ==
                   other.mcpcStableChildBridgeLiveUpdate &&
               mcpcChildLinkBridgeLiveUpdate ==
                   other.mcpcChildLinkBridgeLiveUpdate &&
               mcpcEmptyOwnerReleaseLiveUpdate ==
                   other.mcpcEmptyOwnerReleaseLiveUpdate &&
               mcpcStableChildMoveLiveUpdate ==
                   other.mcpcStableChildMoveLiveUpdate &&
               applied == other.applied &&
               rejectReasonCode == other.rejectReasonCode;
    }
};

std::vector<PartialFReplayEvent> collect_partial_f_replay_trace(
    const ImageUInt8Ptr &image,
    const std::vector<float> &thresholds) {
    noninc::NonIncreasingCasfKernel runner(image, 1.0);
    runner.setMergeEventCollectionEnabled(true);

    std::vector<PartialFReplayEvent> trace;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        for (const auto &event : runner.lastMergeEvents()) {
            if (!event.partialFAttempted) {
                continue;
            }
            trace.push_back({
                threshold,
                event.phase,
                event.nodeId,
                event.parentId,
                event.nodeLevel,
                event.parentLevel,
                event.componentCount,
                event.totalComponentPixels,
                event.partialFAttempted,
                event.partialFAccepted,
                event.partialFRejected,
                event.partialFMutationAttempted,
                event.partialFMcpcLiveUpdate,
                event.partialFMcpcHiddenSupportRedundantLiveUpdate,
                event.partialFMcpcParentLinkBridgeLiveUpdate,
                event.partialFMcpcStableChildBridgeLiveUpdate,
                event.partialFMcpcChildLinkBridgeLiveUpdate,
                event.partialFMcpcEmptyOwnerReleaseLiveUpdate,
                event.partialFMcpcStableChildMoveLiveUpdate,
                event.applied,
                event.partialFRejectReasonCode});
        }
    }
    return trace;
}

void test_partial_f_replay_is_deterministic() {
    auto image = make_demo_image();
    const std::vector<float> thresholds = {0.0f, 10.0f, 200.0f};

    const auto run1 = collect_partial_f_replay_trace(image, thresholds);
    const auto run2 = collect_partial_f_replay_trace(image, thresholds);

    require(!run1.empty(),
            "partial F deterministic replay must observe at least one event");
    require(run1.size() == run2.size(),
            "partial F deterministic replay must preserve event count");
    require(run1 == run2,
            "partial F deterministic replay must preserve decisions and reasons");
}

void test_phase_rebuild_mode_matches_naive_baseline() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 10.0f, 200.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::RebuildAfterPhase);
    const auto explicitOutput = runner.filter(thresholds);
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(explicitOutput->isEqual(naiveOutput),
            "phase-rebuild non-increasing CASF must match naive rebuild baseline");
}

void test_materialized_explicit_mode_matches_naive_baseline() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 10.0f, 200.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MaterializeThenExplicit);
    const auto explicitOutput = runner.filter(thresholds);
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(explicitOutput->isEqual(naiveOutput),
            "materialize-then-explicit non-increasing CASF must match naive rebuild baseline");
}

void test_explicit_only_reduced_fixture_captures_missing_split_case() {
    auto image = make_reduced_split_fixture_image();
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    for (float threshold : {0.0f, 1.0f, 5.0f, 10.0f}) {
        runner.applyThreshold(threshold);
        if (runner.hasFirstInvariantMismatch()) {
            break;
        }
    }

    require(runner.hasFirstInvariantMismatch(),
            "reduced fixture must expose the current explicit-only split limitation");

    const auto &event = runner.firstInvariantMismatchEvent();
    require(event.phase == noninc::NonIncreasingPhase::MaxToMin,
            "reduced split fixture should fail in the max->min phase");
    require(event.nodeLevel == 92 && event.parentLevel == 80,
            "reduced split fixture should preserve the failing level transition");
    require(event.numChildren == 1,
            "reduced split fixture should fail on a node merge that preserves children");
    require(event.directProperParts == 1 && event.componentCount == 1,
            "reduced split fixture should show that even one connected direct proper part can require a split");

    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    noninc::NonIncreasingCasfKernel phaseRebuildRunner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::RebuildAfterPhase);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};
    const auto phaseRebuildOutput = phaseRebuildRunner.filter(thresholds);
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(phaseRebuildOutput->isEqual(naiveOutput),
            "phase-rebuild guard must still match the naive baseline on the reduced split fixture");

    noninc::NonIncreasingCasfKernel materializedRunner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MaterializeThenExplicit);
    materializedRunner.setStopOnFirstInvariantMismatchEnabled(true);
    const auto materializedOutput = materializedRunner.filter(thresholds);

    require(!materializedRunner.hasFirstInvariantMismatch(),
            "materialize-then-explicit mode should preserve the reconstruction invariant on the reduced split fixture");
    require(materializedOutput->isEqual(naiveOutput),
            "materialize-then-explicit mode must match the naive baseline on the reduced split fixture");
}

void test_explicit_only_marks_multiple_nodeca_preflight_failure() {
    auto image = make_multiple_nodeca_p1_restricted_fixture_image();
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);
    runner.setMergeEventCollectionEnabled(true);

    bool foundEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold = 0.0f; threshold <= 12.0f; threshold += 1.0f) {
        runner.applyThreshold(threshold);
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.mergeNodeLocalAttempted &&
                event.mergeNodeLocalProblem1FailureMode == "multiple_nodeCa") {
                targetEvent = event;
                foundEvent = true;
                break;
            }
        }
        if (foundEvent) {
            break;
        }
    }

    require(foundEvent,
            "multiple-nodeCa synthetic fixture should trigger explicit-only Problem 1 preflight rejection");
    require(!targetEvent.mergeNodeLocalAccepted,
            "multiple-nodeCa preflight must not be accepted as local explicit-only update");
    require(targetEvent.mergeNodeLocalRejectedComponents > 0,
            "multiple-nodeCa preflight should count rejected components");
    require(targetEvent.mergeNodeLocalProblem1PreflightRejectedByPublishedProblem > 0,
            "multiple-nodeCa preflight should count publication-level rejection");
    require(targetEvent.mergeNodeLocalFailureSummary.find("single nodeCa") != std::string::npos,
            "multiple-nodeCa rejection should expose the single-nodeCa diagnostic message");
}

NodeId find_reduced_fixture_merge_node(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (nodeId == tree.getRoot() || !tree.isAlive(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
            continue;
        }
        if (tree.getAltitude(nodeId) == 92 &&
            tree.getAltitude(parentId) == 80 &&
            tree.getNumChildren(nodeId) == 1 &&
            tree.getNumProperParts(nodeId) == 1) {
            return nodeId;
        }
    }
    return InvalidNode;
}

NodeId find_p1_mcpc_fixture_merge_node(const DynamicComponentTree &tree,
                                       int expectedDirectProperParts) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (nodeId == tree.getRoot() || !tree.isAlive(nodeId)) {
            continue;
        }
        const NodeId parentId = tree.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId || !tree.isAlive(parentId)) {
            continue;
        }
        if (tree.getAltitude(nodeId) == 92 &&
            tree.getAltitude(parentId) == 80 &&
            tree.getNumChildren(nodeId) == 1 &&
            tree.getNumProperParts(nodeId) == expectedDirectProperParts) {
            return nodeId;
        }
    }
    return InvalidNode;
}

std::vector<int> copied_materialization_non_connector_levels(
    const noninc::NonIncreasingStrippedLocalUpdateResult &result) {
    std::vector<int> levels;
    for (std::size_t index = 0;
         index < result.createdNodeLevels.size() &&
         index < result.createdNodePlanIndices.size();
         ++index) {
        if (result.createdNodePlanIndices[index] >= 0) {
            levels.push_back(result.createdNodeLevels[index]);
        }
    }
    std::sort(levels.begin(), levels.end());
    return levels;
}

void test_oracle_guided_local_update_repairs_reduced_split_event() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::ExplicitOnly);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "reduced fixture should not fail before the local update target event");
    }

    const NodeId nodeId = find_reduced_fixture_merge_node(runner.maxTree());
    require(nodeId != InvalidNode,
            "reduced fixture target merge node must exist before threshold 10");

    const auto result = noninc::applyOracleGuidedLocalMergeUpdate(
        runner.maxTree(),
        runner.minTree(),
        adjacency,
        nodeId);

    require(result.applied,
            "oracle-guided local update must apply the reduced fixture target merge");
    require(result.oracleBuilt,
            "oracle-guided local update must build the local oracle");
    require(result.localRecompositionApplied,
            "oracle-guided local update must materialize local recomposition nodes");
    require(result.localRecompositionMatchesTarget,
            "oracle-guided local update must reconstruct the local target image");
    require(result.dualMatchesPrimalAfterMerge,
            "oracle-guided local update must leave primal and dual reconstructions equal");
    require(result.localStructureMatchesOracle,
            "oracle-guided local update must match the reduced fixture local oracle structure");
    require(result.recompositionNodeCount >= 1 && result.createdNodeCount >= 1,
            "oracle-guided local update must create the reduced fixture recomposition node");
}

void test_stripped_p1_mcpc_fixture_creates_target_b_node() {
    auto image = make_p1_mcpc_target_b_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);
    runner.setMergeEventCollectionEnabled(true);

    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.checkNoEmptyLiveNodesAfterCanonicalization = true;
    options.compareLocalStructureAgainstRebuild = true;
    runner.setStrippedLocalUpdateOptions(options);

    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
    }

    const noninc::NonIncreasingMergeEvent *targetEvent = nullptr;
    for (const auto &event : runner.lastMergeEvents()) {
        if (event.threshold == 10.0f &&
            event.phase == noninc::NonIncreasingPhase::MaxToMin &&
            event.nodeLevel == 92 &&
            event.parentLevel == 80 &&
            event.numChildren == 1 &&
            event.directProperParts == 2 &&
            event.componentCount == 2 &&
            event.strippedLocalCreatedNodes == 1) {
            targetEvent = &event;
            break;
        }
    }

    require(targetEvent != nullptr,
            "P1-MCPC fixture must expose a multi-component preserved-context merge");
    require(targetEvent->strippedLocalAccepted,
            "P1-MCPC fixture local stripped update must be accepted");
    require(targetEvent->strippedLocalMatchesTarget,
            "P1-MCPC fixture local recomposition must match the target image");
    require(targetEvent->strippedLocalStructureChecked,
            "P1-MCPC fixture must compare local structure against rebuild");
    require(targetEvent->strippedLocalStructureMatchesRebuild,
            "P1-MCPC fixture local structure must match rebuild");
    require(targetEvent->strippedLocalDualMatchesPrimalAfterMerge,
            "P1-MCPC fixture must keep primal and dual reconstructions equal");
    require(!targetEvent->strippedLocalFallbackToRebuild,
            "P1-MCPC fixture should not need rebuild fallback");
    require(targetEvent->strippedCreatedNodeLevels.size() == 1 &&
                targetEvent->strippedCreatedNodeLevels[0] == targetEvent->parentLevel,
            "P1-MCPC fixture must create the target-b recomposition node");
    require(targetEvent->strippedEmptyCreatedNodesAfterCanonicalization == 0,
            "P1-MCPC fixture must not leave empty created nodes");

    const auto naiveOutput =
        noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(runner.reconstructImage()->isEqual(naiveOutput),
            "P1-MCPC fixture stripped sequence must match the naive baseline");
}

void test_stripped_p1_mcpc_fixture_records_auxiliary_connector() {
    auto image = make_p1_mcpc_connector_fixture_image();
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);
    runner.setMergeEventCollectionEnabled(true);

    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.checkNoEmptyLiveNodesAfterCanonicalization = true;
    options.compareLocalStructureAgainstRebuild = true;
    runner.setStrippedLocalUpdateOptions(options);

    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
    }

    const noninc::NonIncreasingMergeEvent *targetEvent = nullptr;
    for (const auto &event : runner.lastMergeEvents()) {
        if (event.threshold == 10.0f &&
            event.phase == noninc::NonIncreasingPhase::MaxToMin &&
            event.nodeLevel == 92 &&
            event.parentLevel == 80 &&
            event.numChildren == 1 &&
            event.directProperParts == 2 &&
            event.componentCount == 2) {
            targetEvent = &event;
            break;
        }
    }

    require(targetEvent != nullptr,
            "P1-MCPC connector fixture must expose the target event");
    require(targetEvent->strippedCreatedConnectorNodes == 0,
            "P1-MCPC connector fixture must no longer create an auxiliary connector");
    require(targetEvent->strippedCreatedStableChildBridgeNodes == 0,
            "P1-MCPC connector fixture must omit stable-child bridges by construction");
    require(targetEvent->strippedLocalStructureMatchesRebuild,
            "P1-MCPC connector fixture must still match rebuild structure");
    require(targetEvent->strippedLocalCreatedNodes == 1,
            "P1-MCPC connector fixture must create only the target node");
    require(targetEvent->strippedCreatedNodeLevels.size() == 1 &&
                targetEvent->strippedCreatedNodeLevels[0] == 80,
            "P1-MCPC connector fixture must preserve the target node level");
    require(targetEvent->strippedCreatedNodePlanIndices.size() == 1 &&
                targetEvent->strippedCreatedNodePlanIndices[0] == 0,
            "P1-MCPC connector fixture must keep only the planned target node");
    require(targetEvent->strippedEmptyCreatedNodesAfterCanonicalization == 0,
            "P1-MCPC connector fixture must not leave empty created nodes");
    require(targetEvent->strippedPostCanonicalizationChangesRequired == 0,
            "P1-MCPC connector fixture must not require corrective post-canonicalization");
}

void test_stripped_p1_mcpc_fixture_preserves_multicomponent_level_b_structure() {
    auto image = make_p1_mcpc_structure_mismatch_fixture_image();
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);
    runner.setMergeEventCollectionEnabled(true);

    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.checkNoEmptyLiveNodesAfterCanonicalization = true;
    options.compareLocalStructureAgainstRebuild = true;
    runner.setStrippedLocalUpdateOptions(options);

    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
    }

    const noninc::NonIncreasingMergeEvent *targetEvent = nullptr;
    for (const auto &event : runner.lastMergeEvents()) {
        if (event.threshold == 10.0f &&
            event.phase == noninc::NonIncreasingPhase::MaxToMin &&
            event.nodeLevel == 92 &&
            event.parentLevel == 80 &&
            event.numChildren == 1 &&
            event.directProperParts == 3 &&
            event.componentCount == 2 &&
            event.strippedLocalCreatedNodes > 0 &&
            event.strippedLocalStructureChecked &&
            event.strippedLocalStructureMatchesRebuild) {
            targetEvent = &event;
            break;
        }
    }

    require(targetEvent != nullptr,
            "P1-MCPC structure fixture must preserve local/rebuild structure");
    require(targetEvent->strippedLocalMatchesTarget,
            "P1-MCPC multicomponent fixture must match target reconstruction");
    require(targetEvent->strippedLocalDualMatchesPrimalAfterMerge,
            "P1-MCPC multicomponent fixture must preserve primal/dual reconstruction equality");
    require(std::find(targetEvent->strippedCreatedNodeLevels.begin(),
                      targetEvent->strippedCreatedNodeLevels.end(),
                      targetEvent->parentLevel) !=
                targetEvent->strippedCreatedNodeLevels.end(),
            "P1-MCPC multicomponent fixture must create target-b node");
    require(targetEvent->strippedLocalRecompositionNodes == 1,
            "P1-MCPC multicomponent fixture must keep one canonical recomposition node");
    require(targetEvent->strippedEmptyConnectorNodesAfterCanonicalization == 0,
            "P1-MCPC multicomponent fixture must not leave empty auxiliary connectors");
    require(targetEvent->strippedLocalStructureMismatchSummary.empty(),
            "P1-MCPC multicomponent fixture must not record a structure mismatch");
}

void test_mcpc_direct_preflight_matches_copied_target_b_materialization_plan() {
    auto image = make_p1_mcpc_target_b_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
    }

    const NodeId nodeId = find_p1_mcpc_fixture_merge_node(runner.maxTree(), 2);
    require(nodeId != InvalidNode,
            "P1-MCPC direct preflight target-b fixture must expose the merge before threshold 10");

    const auto preflight = noninc::preflightMcpcDirectMergeNodeUpdate(
        runner.maxTree(),
        runner.minTree(),
        adjacency,
        nodeId);
    require(preflight.checked && preflight.inputBuilt,
            "P1-MCPC direct preflight must build local input");
    require(preflight.partitionPlanBuilt && preflight.recompositionPlanBuilt,
            "P1-MCPC direct preflight must build partition and recomposition plans");
    require(preflight.planOracleAuditChecked && preflight.planOracleFree,
            "P1-MCPC direct preflight plan must be oracle-free");
    require(preflight.materializationOracleAuditChecked &&
                preflight.materializationOracleFree,
            "P1-MCPC direct preflight materialization context must be oracle-free");
    require(preflight.materializationInputsValid && preflight.rollbackFreeSummaryBuilt,
            "P1-MCPC direct preflight must validate materialization inputs and summarize mutation");
    require(preflight.passed,
            "P1-MCPC direct preflight must accept the target-b fixture before direct mutation");
    require(preflight.directV1PrecheckEvaluated && preflight.directV1PrecheckAccepted,
            "P1-MCPC direct-v1 precheck must accept the target-b fixture before copy materialization");

    DynamicComponentTree primalTrial = runner.maxTree();
    DynamicComponentTree dualTrial = runner.minTree();
    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.checkNoEmptyLiveNodesAfterCanonicalization = true;
    options.compareLocalStructureAgainstRebuild = true;
    const auto copied = noninc::applyStrippedNonOracleLocalMergeUpdate(
        primalTrial,
        dualTrial,
        adjacency,
        nodeId,
        options);

    require(copied.applied && copied.localRecompositionMatchesTarget,
            "copied P1-MCPC target-b fixture must remain accepted");
    require(preflight.summary.plannedNodeCount == copied.recompositionNodeCount,
            "P1-MCPC direct preflight must predict the target-b recomposition node count");
    require(preflight.summary.plannedDirectProperPartCount > 0,
            "P1-MCPC direct preflight must predict moved direct proper parts");
    require(preflight.summary.plannedNodeLevels ==
                copied_materialization_non_connector_levels(copied),
            "P1-MCPC direct preflight must predict non-connector created-node levels");
    require(preflight.summary.invalidStableChildRefCount == 0 &&
                preflight.summary.stableChildCycleRiskCount == 0,
            "P1-MCPC direct preflight must not detect rollback-free risks on target-b fixture");
}

void test_mcpc_direct_preflight_separates_connector_from_planned_nodes() {
    auto image = make_p1_mcpc_connector_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
    }

    const NodeId nodeId = find_p1_mcpc_fixture_merge_node(runner.maxTree(), 2);
    require(nodeId != InvalidNode,
            "P1-MCPC direct preflight connector fixture must expose the merge before threshold 10");

    const auto preflight = noninc::preflightMcpcDirectMergeNodeUpdate(
        runner.maxTree(),
        runner.minTree(),
        adjacency,
        nodeId);
    require(preflight.passed,
            "P1-MCPC direct preflight must accept the connector fixture as a plan candidate");
    require(preflight.directV1PrecheckEvaluated && preflight.directV1PrecheckAccepted,
            "P1-MCPC direct-v1 precheck must accept the connector fixture after canonical construction");
    require(preflight.directV1PrecheckRejectionSummary.empty(),
            "P1-MCPC direct-v1 connector fixture must not report a canonical rejection");

    DynamicComponentTree primalTrial = runner.maxTree();
    DynamicComponentTree dualTrial = runner.minTree();
    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.checkNoEmptyLiveNodesAfterCanonicalization = true;
    options.compareLocalStructureAgainstRebuild = true;
    const auto copied = noninc::applyStrippedNonOracleLocalMergeUpdate(
        primalTrial,
        dualTrial,
        adjacency,
        nodeId,
        options);

    require(copied.applied && copied.localRecompositionMatchesTarget,
            "copied P1-MCPC connector fixture must remain accepted");
    require(copied.createdConnectorNodeCount == 0,
            "copied P1-MCPC connector fixture must not materialize an auxiliary connector");
    require(preflight.summary.predictedAuxiliaryNodeCount == copied.createdConnectorNodeCount,
            "P1-MCPC direct preflight must predict the connector fixture auxiliary count");
    require(preflight.summary.predictedAuxiliaryNodeCount == 0,
            "P1-MCPC direct preflight must expose canonical zero-auxiliary construction");
    require(preflight.summary.predictedStableChildBridgeNodeCount ==
                copied.createdStableChildBridgeNodeCount,
            "P1-MCPC direct preflight must predict omitted stable-child bridges");
    require(preflight.summary.predictedHiddenExternalParentConnectorNodeCount ==
                copied.createdHiddenExternalParentConnectorNodeCount,
            "P1-MCPC direct preflight must not predict hidden external connectors for the connector fixture");
    require(preflight.summary.plannedNodeCount == copied.recompositionNodeCount,
            "P1-MCPC direct preflight must count only canonical planned recomposition nodes");
    require(copied.createdNodeCount == preflight.summary.plannedNodeCount,
            "P1-MCPC copied materialization must create only planned canonical nodes");
    require(copied.postCanonicalizationChangesRequired == 0,
            "P1-MCPC connector fixture must not require corrective post-canonicalization");
    require(preflight.summary.plannedNodeLevels ==
                copied_materialization_non_connector_levels(copied),
            "P1-MCPC direct preflight connector fixture must predict non-connector node levels");
}

void test_mcpc_direct_v1_acceptance_diagnostics_enforce_nonoracle_guards() {
    noninc::NonIncreasingMcpcDirectPreflightResult preflight;
    preflight.passed = true;

    noninc::NonIncreasingStrippedLocalUpdateResult strippedUpdate;
    strippedUpdate.applied = true;
    strippedUpdate.localRecompositionTargetChecked = true;
    strippedUpdate.localRecompositionMatchesTarget = true;
    strippedUpdate.finalReconstructionChecked = true;
    strippedUpdate.dualMatchesPrimalAfterMerge = true;

    auto diagnostics =
        noninc::diagnoseMcpcDirectPreflightV1Acceptance(
            preflight,
            strippedUpdate);
    require(diagnostics.evaluated && diagnostics.accepted,
            "direct-v1 diagnostics should accept a clean stripped update");

    strippedUpdate.nonOracleHiddenSupportCandidateCount = 1;
    diagnostics =
        noninc::diagnoseMcpcDirectPreflightV1Acceptance(
            preflight,
            strippedUpdate);
    require(!diagnostics.accepted &&
                diagnostics.rejectionSummary.find("hidden-support candidate") !=
                    std::string::npos,
            "direct-v1 diagnostics must reject non-redundant hidden-support candidates");

    strippedUpdate.nonOracleHiddenSupportAlreadyRepresentedCount = 1;
    strippedUpdate.nonOracleParentLinkTouchesPlanCount = 2;
    strippedUpdate.nonOracleParentLinkBothRepresentedCount = 1;
    diagnostics =
        noninc::diagnoseMcpcDirectPreflightV1Acceptance(
            preflight,
            strippedUpdate);
    require(!diagnostics.accepted &&
                diagnostics.rejectionSummary.find("without both endpoints represented") !=
                    std::string::npos,
            "direct-v1 diagnostics must reject parent links without both endpoints represented");

    strippedUpdate.nonOracleParentLinkBothRepresentedCount = 2;
    strippedUpdate.nonOracleParentLinkBridgeMaterializableCount = 1;
    diagnostics =
        noninc::diagnoseMcpcDirectPreflightV1Acceptance(
            preflight,
            strippedUpdate);
    require(!diagnostics.accepted &&
                diagnostics.rejectionSummary.find("parent-link bridge") !=
                    std::string::npos,
            "direct-v1 diagnostics must reject materializable parent-link bridges");
}

void test_mcpc_direct_v1_precheck_rejects_synthetic_nonoracle_failures() {
    noninc::NonIncreasingMcpcDirectPreflightResult preflight;
    preflight.passed = true;

    std::string rejectionSummary;
    require(noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary),
            "direct-v1 precheck should accept a clean synthetic preflight");
    require(rejectionSummary.empty(),
            "clean direct-v1 precheck should not report a rejection reason");

    preflight.directV1Precheck.predictedAuxiliaryNodeCount = 1;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("auxiliary nodes") != std::string::npos,
            "direct-v1 precheck must reject predicted auxiliary connectors");

    preflight.directV1Precheck = {};
    preflight.directV1Precheck.sameLevelPlannedParentChildCount = 1;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("same-level planned parent-child") !=
                    std::string::npos,
            "direct-v1 precheck must reject planned same-level parent-child relations");

    preflight.directV1Precheck = {};
    preflight.directV1Precheck.sameLevelExternalAttachmentCount = 1;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("same-level external attachment") !=
                    std::string::npos,
            "direct-v1 precheck must reject same-level external attachments");

    preflight.directV1Precheck = {};
    preflight.directV1Precheck.nonOracleHiddenSupportCandidateCount = 1;
    preflight.directV1Precheck.nonOracleHiddenSupportNotRedundantCount = 1;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("hidden-support candidate") !=
                    std::string::npos,
            "direct-v1 precheck must reject non-redundant hidden-support candidates");

    preflight.directV1Precheck = {};
    preflight.directV1Precheck.nonOracleParentLinkTouchesPlanCount = 2;
    preflight.directV1Precheck.nonOracleParentLinkBothRepresentedCount = 1;
    preflight.directV1Precheck.nonOracleParentLinkNotBothRepresentedCount = 1;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("without both endpoints represented") !=
                    std::string::npos,
            "direct-v1 precheck must reject parent links without both endpoints represented");

    preflight.directV1Precheck = {};
    preflight.directV1Precheck.nonOracleParentLinkTouchesPlanCount = 1;
    preflight.directV1Precheck.nonOracleParentLinkBothRepresentedCount = 1;
    preflight.directV1Precheck.nonOracleParentLinkInvalidLevelCount = 1;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("invalid non-oracle parent-link level") !=
                    std::string::npos,
            "direct-v1 precheck must reject invalid parent-link levels");

    preflight.directV1Precheck = {};
    preflight.passed = false;
    require(!noninc::evaluateMcpcDirectV1PrecheckAcceptance(
                preflight,
                &rejectionSummary) &&
                rejectionSummary.find("base preflight did not pass") !=
                    std::string::npos,
            "direct-v1 precheck must reject when the base preflight fails");
}

void test_oracle_guided_singleton_merge_preserves_equal_support_parent_chain() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image, true, adjacency);
    DynamicComponentTree minTree(image, false, adjacency);

    const PixelId singletonPeak = 1 * image->getNumCols() + 8;
    const NodeId nodeId = maxTree.getSmallestComponent(singletonPeak);
    require(nodeId != InvalidNode,
            "controlled singleton peak must have a max-tree node");
    require(maxTree.getAltitude(nodeId) == 160,
            "controlled singleton peak must be the 160-level max-tree node");

    const auto result = noninc::applyOracleGuidedLocalMergeUpdate(
        maxTree,
        minTree,
        adjacency,
        nodeId);

    require(result.applied,
            "oracle-guided local update must apply the singleton peak merge");
    require(result.localRecompositionMatchesTarget,
            "singleton peak local update must reconstruct the target image");
    require(result.dualMatchesPrimalAfterMerge,
            "singleton peak local update must keep primal and dual reconstructions equal");
    require(result.localStructureMatchesOracle,
            "singleton peak local update must preserve equal-support local parent chains");
}

void test_canonicalization_reparents_adjacent_comparable_siblings() {
    auto image = ImageUInt8::create(1, 3);
    (*image)[0] = 5;
    (*image)[1] = 9;
    (*image)[2] = 0;
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image, true, adjacency);

    const NodeId midNode = maxTree.getSmallestComponent(0);
    const NodeId highNode = maxTree.getSmallestComponent(1);
    require(midNode != InvalidNode && highNode != InvalidNode && midNode != highNode,
            "adjacent comparable sibling fixture must have distinct max-tree nodes");
    require(maxTree.getNodeParent(highNode) == midNode,
            "fixture must initially place the high node under the mid-level node");

    maxTree.moveNode(highNode, maxTree.getRoot());
    require(maxTree.getNodeParent(highNode) == maxTree.getRoot(),
            "fixture setup must expose adjacent comparable siblings");

    const int changes = noninc::canonicalizeAfterLocalRecomposition(maxTree, adjacency.get());
    require(changes > 0,
            "canonicalization must modify the stale adjacent comparable sibling topology");
    require(maxTree.getNodeParent(highNode) == midNode,
            "canonicalization must reparent the more extreme adjacent sibling under the less extreme sibling");
    require(maxTree.reconstructionImage()->isEqual(image),
            "canonicalization must preserve the reconstructed image");
}

void test_canonicalization_removes_empty_focused_connector() {
    auto image = ImageUInt8::create(1, 2);
    (*image)[0] = 0;
    (*image)[1] = 10;
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree maxTree(image, true, adjacency);

    const NodeId childNode = maxTree.getSmallestComponent(1);
    require(childNode != InvalidNode && childNode != maxTree.getRoot(),
            "empty connector fixture must have a non-root high-level child");
    const NodeId connectorNode = maxTree.createDetachedNode(5);
    require(connectorNode != InvalidNode,
            "empty connector fixture must create an intermediate connector");

    maxTree.attachNode(maxTree.getRoot(), connectorNode);
    maxTree.moveNode(childNode, connectorNode);
    require(maxTree.getNumProperParts(connectorNode) == 0,
            "fixture connector must have no direct proper parts before canonicalization");

    const std::vector<NodeId> focusNodes = {connectorNode, childNode};
    const int changes = noninc::canonicalizeAfterLocalRecomposition(
        maxTree,
        adjacency.get(),
        std::span<const NodeId>(focusNodes),
        std::span<const NodeId>{});
    require(changes > 0,
            "canonicalization must remove the focused empty connector");
    require(!maxTree.isAlive(connectorNode),
            "focused connector without direct proper parts must not survive canonicalization");
    require(maxTree.getNodeParent(childNode) == maxTree.getRoot(),
            "connector child must be promoted to the connector parent");
    require(maxTree.reconstructionImage()->isEqual(image),
            "removing an empty connector must preserve reconstruction");
    require_tree_has_no_empty_live_nodes(maxTree, "max-tree after empty connector canonicalization");
}

void test_phase_structure_signature_detects_candidate_equivalent_connector() {
    auto image = ImageUInt8::create(1, 3);
    (*image)[0] = 5;
    (*image)[1] = 9;
    (*image)[2] = 0;
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree localTree(image, true, adjacency);
    DynamicComponentTree rebuiltTree(image, true, adjacency);

    const NodeId midNode = localTree.getSmallestComponent(0);
    const NodeId highNode = localTree.getSmallestComponent(1);
    require(midNode != InvalidNode && highNode != InvalidNode && midNode != highNode,
            "candidate-equivalent connector fixture must expose a parent/child chain");
    require(localTree.getNodeParent(highNode) == midNode,
            "candidate-equivalent connector fixture must initially place the high node under the mid node");

    const NodeId connectorNode = localTree.createDetachedNode(7);
    require(connectorNode != InvalidNode,
            "candidate-equivalent connector fixture must create a connector");
    localTree.attachNode(midNode, connectorNode);
    localTree.moveNode(highNode, connectorNode);

    require(localTree.reconstructionImage()->isEqual(image),
            "inserted connector must preserve reconstruction");

    auto comparison = noninc::compareNonIncreasingPhaseAgainstRebuild(
        localTree,
        rebuiltTree,
        -1.0f,
        noninc::NonIncreasingPhase::MinToMax,
        false,
        adjacency.get());
    require(comparison.candidateSetMatchesRebuild,
            "negative threshold should make candidate signatures equal despite the connector");
    require(comparison.candidateScanOrderMatchesRebuild,
            "empty candidate stream must also match in scan order");
    require(comparison.releaseH1ReconstructionMatches,
            "release diagnostic H1 must record reconstruction equality");
    require(comparison.releaseH2CandidateSetMatches,
            "release diagnostic H2 must record candidate-set equality");
    require(comparison.releaseH5CandidateScanOrderMatches,
            "release diagnostic H5-order must record candidate stream equality");
    require(!comparison.releasePredicateEvaluated,
            "release predicate must not be evaluated without a phase focus");
    require(!comparison.boundaryStructureCompared,
            "candidate-equivalent connector comparison should not need mismatch-driven boundary diagnostics");
    require(!noninc::nonIncreasingTreeStructureSignaturesMatch(localTree, rebuiltTree),
            "structural signature must detect the candidate-equivalent connector");
    const std::vector<NodeId> focusNodes = {connectorNode, highNode};
    const auto localPixels = noninc::collectNodeSupportPixels(localTree, connectorNode);
    require(!noninc::nonIncreasingStructureEnvelopeSignaturesMatch(
                localTree,
                rebuiltTree,
                std::span<const NodeId>(focusNodes),
                std::span<const PixelId>(localPixels)),
            "focused structural envelope must detect the candidate-equivalent connector");
    const auto localEnvelopeNodes = noninc::collectNonIncreasingStructureEnvelopeNodes(
        localTree,
        std::span<const NodeId>(focusNodes),
        std::span<const PixelId>(localPixels),
        0);
    const auto rebuiltEnvelopeNodes = noninc::collectNonIncreasingStructureEnvelopeNodes(
        rebuiltTree,
        std::span<const NodeId>{},
        std::span<const PixelId>(localPixels),
        0);
    const auto localBoundaryNodes = noninc::collectNonIncreasingStructureBoundaryNodes(
        localTree,
        std::span<const NodeId>(localEnvelopeNodes),
        adjacency.get());
    const auto rebuiltBoundaryNodes = noninc::collectNonIncreasingStructureBoundaryNodes(
        rebuiltTree,
        std::span<const NodeId>(rebuiltEnvelopeNodes),
        adjacency.get());
    const auto localBoundaryIncidences = noninc::collectNonIncreasingBoundaryIncidenceSignatures(
        localTree,
        std::span<const NodeId>(localEnvelopeNodes),
        std::span<const NodeId>(localBoundaryNodes),
        adjacency.get());
    const auto rebuiltBoundaryIncidences = noninc::collectNonIncreasingBoundaryIncidenceSignatures(
        rebuiltTree,
        std::span<const NodeId>(rebuiltEnvelopeNodes),
        std::span<const NodeId>(rebuiltBoundaryNodes),
        adjacency.get());
    require(!localBoundaryIncidences.empty() || !rebuiltBoundaryIncidences.empty(),
            "focused structural boundary must expose local/rebuild incidences around the connector fixture");

    const int changes = noninc::canonicalizeAfterLocalRecomposition(
        localTree,
        adjacency.get(),
        std::span<const NodeId>(focusNodes),
        std::span<const NodeId>{});
    require(changes > 0,
            "canonicalization must remove the candidate-equivalent connector");

    DynamicComponentTree rebuiltAfter(image, true, adjacency);
    require(noninc::nonIncreasingTreeStructureSignaturesMatch(localTree, rebuiltAfter),
            "structural signature must match rebuild after removing the connector");
    require(noninc::nonIncreasingStructureEnvelopeSignaturesMatch(
                localTree,
                rebuiltAfter,
                std::span<const NodeId>(focusNodes),
                std::span<const PixelId>(localPixels)),
            "focused structural envelope must match rebuild after removing the connector");
}

void test_oracle_guided_sequence_matches_naive_on_controlled_demo() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 10.0f, 200.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::OracleGuidedLocalRecomposition);
    runner.setMergeEventCollectionEnabled(true);
    runner.setPhaseOracleComparisonEnabled(true);
    runner.setLeafExplicitFastPathEnabled(false);
    auto naiveOutput = image->clone();
    int oracleLocalEvents = 0;
    int oracleLocalTargetMatches = 0;
    int oracleLocalStructureMatches = 0;
    int oracleLocalEventsWithAffectedRegion = 0;
    int phaseOracleComparisons = 0;
    int phaseOracleCandidateMismatches = 0;
    int phaseOracleReconstructionMismatches = 0;
    noninc::NonIncreasingPhaseOracleMismatchFixture firstMismatchFixture;
    noninc::NonIncreasingPhaseOracleMismatchFixture smallestMismatchFixture;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        naiveOutput = noninc::applyNaiveGrayVarianceMergeThreshold(naiveOutput, adjacency, threshold);
        const auto comparisonSpan =
            std::span<const noninc::NonIncreasingPhaseOracleComparison>(runner.lastPhaseOracleComparisons());
        if (!firstMismatchFixture.found) {
            firstMismatchFixture = noninc::buildFirstPhaseOracleMismatchFixture(comparisonSpan, 1);
        }
        auto thresholdSmallestFixture =
            noninc::buildSmallestPhaseOracleMismatchFixture(comparisonSpan, 0);
        if (thresholdSmallestFixture.found &&
            (!smallestMismatchFixture.found ||
             thresholdSmallestFixture.croppedImage->getSize() < smallestMismatchFixture.croppedImage->getSize())) {
            smallestMismatchFixture = std::move(thresholdSmallestFixture);
        }
        for (const auto &comparison : runner.lastPhaseOracleComparisons()) {
            ++phaseOracleComparisons;
            if (!comparison.localReconstructionMatchesRebuild) {
                ++phaseOracleReconstructionMismatches;
            }
            if (!comparison.candidateSetMatchesRebuild) {
                ++phaseOracleCandidateMismatches;
            }
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (!event.oracleGuidedLocalApplied) {
                continue;
            }
            ++oracleLocalEvents;
            if (event.totalComponentPixels > 0 &&
                event.minRow >= 0 &&
                event.minCol >= 0 &&
                event.maxRow >= event.minRow &&
                event.maxCol >= event.minCol &&
                !event.affectedPixelsPreview.empty()) {
                ++oracleLocalEventsWithAffectedRegion;
            }
            if (event.oracleGuidedLocalMatchesTarget) {
                ++oracleLocalTargetMatches;
            }
            if (event.oracleGuidedLocalStructureMatchesOracle) {
                ++oracleLocalStructureMatches;
            }
        }
    }
    const auto oracleLocalOutput = runner.reconstructImage();

    require(oracleLocalEvents > 0,
            "oracle-guided sequence diagnostic must record local recomposition events");
    require(oracleLocalTargetMatches == oracleLocalEvents,
            "oracle-guided sequence must keep each local recomposition image-correct");
    require(oracleLocalEventsWithAffectedRegion > 0,
            "oracle-guided diagnostics must record affected region bbox and pixel previews");
    require(oracleLocalStructureMatches == oracleLocalEvents,
            "oracle-guided controlled sequence must match each local oracle structure");
    require(phaseOracleComparisons > 0,
            "oracle-guided sequence diagnostic must compare local trees against rebuilds by phase");
    require(phaseOracleReconstructionMismatches == 0,
            "oracle-guided phase comparisons should remain image-correct in the controlled sequence");
    require(phaseOracleCandidateMismatches == 0,
            "oracle-guided controlled sequence must keep phase candidate sets aligned with rebuild");
    require(!firstMismatchFixture.found,
            "phase-oracle diagnostics should not find a first mismatch in the controlled sequence");
    require(!smallestMismatchFixture.found,
            "phase-oracle diagnostics should not find a smallest mismatch in the controlled sequence");
    require(oracleLocalOutput->isEqual(naiveOutput),
            "oracle-guided controlled sequence must match the naive rebuild baseline");
}

void test_stripped_then_oracle_sequence_matches_naive_on_reduced_fixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenOracleLocalRecomposition);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    int strippedAttempts = 0;
    int strippedAccepted = 0;
    int strippedFallbacks = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "stripped-then-oracle mode must preserve reconstruction equality inside each threshold");
        for (const auto &event : runner.lastMergeEvents()) {
            if (!event.strippedLocalAttempted) {
                continue;
            }
            ++strippedAttempts;
            if (event.strippedLocalAccepted) {
                ++strippedAccepted;
            }
            if (event.strippedLocalFallbackToOracle) {
                ++strippedFallbacks;
                require(event.oracleGuidedLocalApplied,
                        "stripped fallback must complete through the oracle-guided local update");
            }
        }
    }

    const auto strippedOutput = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(strippedAttempts > 0,
            "stripped-then-oracle diagnostics must record stripped local attempts");
    require(strippedAccepted > 0,
            "reduced fixture must exercise at least one accepted stripped local update");
    require(strippedAccepted + strippedFallbacks == strippedAttempts,
            "each stripped attempt must either be accepted or fall back to the oracle-guided update");
    require(strippedOutput->isEqual(naiveOutput),
            "stripped-then-oracle reduced fixture sequence must match the naive rebuild baseline");
}

void test_stripped_then_rebuild_fallback_matches_naive_on_reduced_fixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);
    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.buildFullUpdatedImage = false;
    options.checkMaterializationTargetReconstruction = true;
    runner.setStrippedLocalUpdateOptions(options);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    int strippedAttempts = 0;
    int rebuildFallbacks = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "stripped-then-rebuild mode must preserve reconstruction equality after fallback");
        for (const auto &event : runner.lastMergeEvents()) {
            if (!event.strippedLocalAttempted) {
                continue;
            }
            ++strippedAttempts;
            if (event.strippedLocalFallbackToRebuild) {
                ++rebuildFallbacks;
                require(event.dualRebuiltAfterLocalMismatch,
                        "stripped rebuild fallback must mark the dual rebuild event");
                require(event.strippedLocalDualMatchesPrimalAfterMerge,
                        "stripped rebuild fallback must leave dual and primal reconstructions equal");
            }
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(strippedAttempts > 0,
            "stripped-then-rebuild diagnostics must record stripped local attempts");
    require(rebuildFallbacks > 0,
            "forced invalid stripped options must exercise the rebuild fallback path");
    require(output->isEqual(naiveOutput),
            "stripped-then-rebuild reduced fixture sequence must match the naive rebuild baseline");
}

void test_stripped_then_phase_rebuild_fallback_matches_naive_on_reduced_fixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition);
    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.buildFullUpdatedImage = false;
    options.checkMaterializationTargetReconstruction = true;
    runner.setStrippedLocalUpdateOptions(options);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    int phaseFallbacks = 0;
    int phaseSkips = 0;
    int phaseRebuilds = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "stripped-then-phase-rebuild mode must preserve reconstruction equality after phase fallback");
        phaseRebuilds += runner.lastPhaseDualRebuilds();
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.strippedLocalFallbackToPhaseRebuild) {
                ++phaseFallbacks;
            }
            if (event.strippedLocalSkippedByPhaseRebuild) {
                ++phaseSkips;
            }
            require(!event.strippedLocalFallbackToRebuild,
                    "phase-rebuild mode must not use per-event rebuild fallback");
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(phaseFallbacks > 0,
            "forced invalid stripped options must exercise the phase rebuild fallback path");
    require(phaseSkips > 0,
            "after the first phase fallback, later events in the dirty phase must skip dual local updates");
    require(phaseRebuilds > 0,
            "phase fallback must rebuild the stale dual tree at phase boundaries");
    require(output->isEqual(naiveOutput),
            "stripped-then-phase-rebuild reduced fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_phase_rebuild_mode_matches_naive_on_reduced_fixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeUpdateThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    int attempts = 0;
    int accepted = 0;
    int phaseFallbacks = 0;
    int phaseSkips = 0;
    int phaseRebuilds = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "merge-node phase-rebuild mode must preserve reconstruction equality after phase fallback");
        phaseRebuilds += runner.lastPhaseDualRebuilds();
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.mergeNodeLocalAttempted) {
                ++attempts;
            }
            if (event.mergeNodeLocalAccepted) {
                ++accepted;
            }
            if (event.mergeNodeLocalFallbackToPhaseRebuild) {
                ++phaseFallbacks;
            }
            if (event.mergeNodeLocalSkippedByPhaseRebuild) {
                ++phaseSkips;
            }
            require(!event.strippedLocalAttempted,
                    "merge-node phase-rebuild mode must not use stripped local recomposition");
            require(!event.strippedLocalFallbackToRebuild,
                    "merge-node phase-rebuild mode must not use per-event rebuild fallback");
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(attempts > 0,
            "merge-node phase-rebuild mode must exercise direct merge-node update attempts");
    require(accepted > 0,
            "merge-node phase-rebuild mode must accept at least one copied direct merge-node update");
    require(phaseFallbacks > 0,
            "reduced fixture must exercise the merge-node phase rebuild fallback path");
    require(phaseSkips > 0,
            "after the first merge-node phase fallback, later events in the dirty phase must skip dual local updates");
    require(phaseRebuilds > 0,
            "merge-node phase fallback must rebuild the stale dual tree at phase boundaries");
    require(output->isEqual(naiveOutput),
            "merge-node phase-rebuild reduced fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_certified_phase_rebuild_mode_matches_naive_on_reduced_fixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    int attempts = 0;
    int certificationChecks = 0;
    int certificationPassed = 0;
    int directCommits = 0;
    int phaseFallbacks = 0;
    int phaseSkips = 0;
    int phaseRebuilds = 0;
    int copyTimedEvents = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "certified merge-node phase-rebuild mode must preserve reconstruction equality after phase fallback");
        phaseRebuilds += runner.lastPhaseDualRebuilds();
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.mergeNodeLocalAttempted) {
                ++attempts;
            }
            if (event.mergeNodeLocalCertificationChecked) {
                ++certificationChecks;
            }
            if (event.mergeNodeLocalCertificationPassed) {
                ++certificationPassed;
            }
            if (event.mergeNodeLocalDirectCommit) {
                ++directCommits;
            }
            if (event.mergeNodeLocalFallbackToPhaseRebuild) {
                ++phaseFallbacks;
            }
            if (event.mergeNodeLocalSkippedByPhaseRebuild) {
                ++phaseSkips;
            }
            if (event.mergeNodeLocalCopyMicros > 0) {
                ++copyTimedEvents;
            }
            require(!event.strippedLocalAttempted,
                    "certified merge-node phase-rebuild mode must not use stripped local recomposition");
            require(!event.strippedLocalFallbackToRebuild,
                    "certified merge-node phase-rebuild mode must not use per-event rebuild fallback");
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(attempts > 0,
            "certified merge-node phase-rebuild mode must exercise merge-node attempts");
    require(certificationChecks > 0,
            "certified merge-node phase-rebuild mode must run the Problem 1 certification preflight");
    require(certificationPassed > 0,
            "reduced fixture must expose at least one certified direct merge-node event");
    require(directCommits > 0,
            "certified merge-node phase-rebuild mode must commit at least one event without copying");
    require(phaseFallbacks > 0,
            "reduced fixture must still exercise certified merge-node phase fallback");
    require(phaseSkips > 0,
            "after certified merge-node phase fallback, later events in the dirty phase must skip dual local updates");
    require(phaseRebuilds > 0,
            "certified merge-node phase fallback must rebuild the stale dual tree at phase boundaries");
    require(copyTimedEvents == 0,
            "certified merge-node mode must not run the copy-based acceptance path");
    require(output->isEqual(naiveOutput),
            "certified merge-node phase-rebuild reduced fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_certified_create_b_phase_rebuild_mode_matches_naive_on_reduced_fixture() {
    auto image = make_reduced_split_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    int attempts = 0;
    int certificationChecks = 0;
    int certificationPassed = 0;
    int targetCreationAllowed = 0;
    int targetCreationRequired = 0;
    int rawMissingTargetBComponents = 0;
    int certifiedMissingTargetBComponents = 0;
    int partialOwnerAndMissingTargetBComponents = 0;
    int accepted = 0;
    int applied = 0;
    int directCommits = 0;
    int copiedAcceptedEvents = 0;
    int phaseFallbacks = 0;
    int phaseSkips = 0;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "certified create-b merge-node phase-rebuild mode must preserve reconstruction equality");
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.mergeNodeLocalAttempted) {
                ++attempts;
            }
            if (event.mergeNodeLocalCertificationChecked) {
                ++certificationChecks;
            }
            if (event.mergeNodeLocalCertificationPassed) {
                ++certificationPassed;
            }
            if (event.mergeNodeLocalTargetCreationAllowed) {
                ++targetCreationAllowed;
            }
            if (event.mergeNodeLocalTargetCreationRequired) {
                ++targetCreationRequired;
            }
            rawMissingTargetBComponents += event.mergeNodeLocalRawMissingTargetComponents;
            certifiedMissingTargetBComponents += event.mergeNodeLocalMissingTargetComponents;
            partialOwnerAndMissingTargetBComponents +=
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents;
            if (event.mergeNodeLocalAccepted) {
                ++accepted;
            }
            if (event.mergeNodeLocalApplied) {
                ++applied;
            }
            if (event.mergeNodeLocalDirectCommit) {
                ++directCommits;
            }
            if (event.mergeNodeLocalAccepted &&
                event.mergeNodeLocalTargetCreationAllowed &&
                !event.mergeNodeLocalDirectCommit) {
                ++copiedAcceptedEvents;
            }
            if (event.mergeNodeLocalFallbackToPhaseRebuild) {
                ++phaseFallbacks;
            }
            if (event.mergeNodeLocalSkippedByPhaseRebuild) {
                ++phaseSkips;
            }
            require(!event.strippedLocalAttempted,
                    "certified create-b merge-node mode must not use stripped local recomposition");
            require(!event.strippedLocalFallbackToRebuild,
                    "certified create-b merge-node mode must not use per-event rebuild fallback");
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);

    require(attempts > 0,
            "certified create-b merge-node mode must exercise merge-node attempts");
    require(certificationChecks > 0,
            "certified create-b merge-node mode must run certification preflight");
    require(certificationPassed > 0,
            "reduced fixture must expose at least one certification accepted with target-b creation allowed");
    require(targetCreationAllowed > 0,
            "certified create-b merge-node mode must record that target-b creation was allowed");
    require(targetCreationRequired > 0,
            "reduced fixture must exercise a missing target-b case accepted for trial creation");
    require(rawMissingTargetBComponents >= certifiedMissingTargetBComponents,
            "raw missing target-b count must cover the certified missing target-b count");
    require(certifiedMissingTargetBComponents > 0,
            "reduced fixture must expose certified missing target-b components");
    require(partialOwnerAndMissingTargetBComponents <= rawMissingTargetBComponents,
            "partial-owner and missing target-b count must be contained in the raw missing target-b count");
    require(accepted > 0,
            "certified create-b merge-node mode must accept at least one copied trial update");
    require(applied > 0,
            "certified create-b merge-node mode must apply at least one copied trial update");
    require(directCommits == 0,
            "certified create-b merge-node mode must not mutate the live trees through the direct commit path");
    require(copiedAcceptedEvents > 0,
            "certified create-b merge-node mode must accept through the guarded copy path");
    require(output->isEqual(naiveOutput),
            "certified create-b merge-node reduced fixture sequence must match the naive rebuild baseline");
    (void) phaseFallbacks;
    (void) phaseSkips;
}

void test_merge_node_create_b_real_ruler_fixture_accepts_certified_missing_target_b() {
    auto image = load_dat_crop("dat/misc256/ruler.512.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "ruler create-b fixture must preserve reconstruction equality");
        if (threshold != 1.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 224 &&
                event.parentLevel == 157 &&
                event.mergeNodeLocalRawMissingTargetComponents == 1 &&
                event.mergeNodeLocalMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "ruler crop must expose a real certified missing target-b merge-node event");
    require(targetEvent.nodeId == 48 && targetEvent.parentId == 5,
            "ruler missing target-b fixture must preserve the recorded node ids");
    require(targetEvent.componentCount == 1 && targetEvent.totalComponentPixels == 11,
            "ruler missing target-b fixture must remain a connected Problem 1 event");
    require(targetEvent.mergeNodeLocalPartialOwnerComponents == 0,
            "ruler missing target-b fixture must be free of partial owners");
    require(targetEvent.mergeNodeLocalTargetCreationRequired,
            "ruler fixture must require target-b creation");
    require(targetEvent.mergeNodeLocalAccepted,
            "ruler fixture must accept the guarded create-b trial update");
    require(!targetEvent.mergeNodeLocalDirectCommit,
            "ruler create-b fixture must not use the live direct commit path");
    require(!targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "ruler create-b fixture must not fallback once the copied trial is accepted");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "ruler create-b fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_create_b_real_partial_owner_missing_b_fixture_falls_back() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "5.2.10 partial-owner missing-b fixture must preserve reconstruction equality");
        if (threshold != 10.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 110 &&
                event.parentLevel == 109 &&
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "5.2.10 crop must expose a real partial-owner and missing target-b event");
    require(targetEvent.nodeId == 14 && targetEvent.parentId == 12,
            "5.2.10 partial-owner missing-b fixture must preserve the recorded node ids");
    require(targetEvent.componentCount == 3 && targetEvent.totalComponentPixels == 3,
            "5.2.10 missing-b fixture must remain a multi-component event");
    require(targetEvent.mergeNodeLocalRawMissingTargetComponents == 1,
            "5.2.10 fixture must record raw missing target-b");
    require(targetEvent.mergeNodeLocalMissingTargetComponents == 0,
            "5.2.10 fixture must not classify the event as certified missing target-b");
    require(targetEvent.mergeNodeLocalPartialOwnerComponents == 3,
            "5.2.10 fixture must retain the partial-owner rejection");
    require(!targetEvent.mergeNodeLocalTargetCreationRequired,
            "partial-owner missing-b fixture must not enter the simple create-b acceptance class");
    require(!targetEvent.mergeNodeLocalAccepted,
            "partial-owner missing-b fixture must not accept the simple create-b update");
    require(targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "partial-owner missing-b fixture must use the conservative phase fallback");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "partial-owner missing-b fixture sequence must match the naive rebuild baseline");
}

noninc::NonIncreasingMcpcDirectPreflightResult make_clean_mcpc_d12_preflight() {
    noninc::NonIncreasingMcpcDirectPreflightResult preflight;
    preflight.checked = true;
    preflight.inputBuilt = true;
    preflight.partitionPlanBuilt = true;
    preflight.recompositionPlanBuilt = true;
    preflight.planOracleAuditChecked = true;
    preflight.planOracleFree = true;
    preflight.materializationContextBuilt = true;
    preflight.materializationOracleAuditChecked = true;
    preflight.materializationOracleFree = true;
    preflight.materializationInputsValid = true;
    preflight.rollbackFreeSummaryBuilt = true;
    preflight.passed = true;
    preflight.directV1PrecheckEvaluated = true;
    preflight.directV1PrecheckAccepted = true;
    return preflight;
}

void require_mcpc_d12_reject_reason(
    const noninc::NonIncreasingMcpcDirectPreflightResult &preflight,
    const std::string &reasonCode,
    const std::string &message) {
    const auto d12 = noninc::diagnoseMcpcD12Domain(preflight);
    require(d12.evaluated,
            message + ": D12 diagnostics must mark evaluation");
    require(!d12.accepted,
            message + ": D12 must reject until the live mutation proof is complete");
    require(d12.rejectReasonCode == reasonCode,
            message + ": unexpected D12 reason " + d12.rejectReasonCode);
}

void test_mcpc_d12_domain_diagnostics_reject_formal_predicates() {
    const auto clean = make_clean_mcpc_d12_preflight();
    require_mcpc_d12_reject_reason(
        clean,
        "d12_exception_free_unproven",
        "clean read-only P1-MCPC preflight");

    auto preflight = clean;
    preflight.checked = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_base_preflight_failed",
        "unchecked P1-MCPC preflight");

    preflight = clean;
    preflight.partitionPlanBuilt = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_plan_incomplete",
        "incomplete P1-MCPC plan");

    preflight = clean;
    preflight.planOracleFree = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_plan_oracle_not_free",
        "non oracle-free P1-MCPC plan");

    preflight = clean;
    preflight.materializationContextBuilt = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_materialization_context_missing",
        "missing P1-MCPC materialization context");

    preflight = clean;
    preflight.materializationOracleFree = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_materialization_oracle_not_free",
        "non oracle-free P1-MCPC materialization input");

    preflight = clean;
    preflight.materializationInputsValid = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_materialization_inputs_invalid",
        "invalid P1-MCPC materialization inputs");

    preflight = clean;
    preflight.rollbackFreeSummaryBuilt = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_rollback_summary_missing",
        "missing P1-MCPC rollback-free summary");

    preflight = clean;
    preflight.summary.invalidStableChildRefCount = 1;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_invalid_stable_child_ref",
        "invalid stable-child refs");

    preflight = clean;
    preflight.summary.stableChildCycleRiskCount = 1;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_stable_child_cycle_risk",
        "stable-child cycle risk");

    preflight = clean;
    preflight.directV1PrecheckEvaluated = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_direct_v1_precheck_missing",
        "missing direct-v1 precheck");

    preflight = clean;
    preflight.directV1Precheck.nonOracleHiddenSupportNotRedundantCount = 1;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_hidden_support_not_redundant",
        "non-redundant hidden support");

    preflight = clean;
    preflight.directV1Precheck.nonOracleParentLinkNotBothRepresentedCount = 1;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_parent_link_not_represented",
        "parent link without both endpoints represented");

    preflight = clean;
    preflight.directV1Precheck.nonOracleParentLinkInvalidLevelCount = 1;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_parent_link_invalid_level",
        "invalid parent-link level");

    preflight = clean;
    preflight.directV1PrecheckAccepted = false;
    require_mcpc_d12_reject_reason(
        preflight,
        "d12_direct_v1_precheck_rejected",
        "direct-v1 precheck rejected without finer counters");
}

noninc::NonIncreasingMcpcDirectPreflightArtifacts
make_empty_owner_with_child_dry_run_artifacts(DynamicComponentTree &minTree,
                                             bool moveResidualChild,
                                             bool externalThroughSourceOwner = false,
                                             std::vector<PixelId> planSupport = {0, 1},
                                             bool includeResidualChildSnapshot = true,
                                             bool residualTouchesModified = false) {
    const NodeId residualChild = minTree.getSmallestComponent(0);
    const NodeId sourceOwner = minTree.getSmallestComponent(1);
    require(residualChild != InvalidNode &&
            sourceOwner != InvalidNode &&
            !minTree.isRoot(sourceOwner),
            "empty-owner dry-run fixture must have residual child and source owner");
    require(minTree.getNodeParent(residualChild) == sourceOwner,
            "empty-owner dry-run residual child must start below source owner");
    const NodeId promotedParent = minTree.getNodeParent(sourceOwner);
    require(promotedParent != InvalidNode &&
            promotedParent != sourceOwner &&
            minTree.isAlive(promotedParent),
            "empty-owner dry-run fixture must have a promoted parent");

    noninc::NonIncreasingMcpcDirectPreflightArtifacts artifacts;
    artifacts.preflight = make_clean_mcpc_d12_preflight();
    artifacts.preflight.summary.plannedNodeCount = 1;
    artifacts.preflight.summary.plannedDirectProperPartCount = 1;
    artifacts.preflight.summary.plannedExternalAttachmentCount = 1;
    artifacts.preflight.summary.sourceOwnerTouchedCount = 1;
    artifacts.preflight.summary.sourceOwnerToBeEmptiedCount = 1;
    artifacts.preflight.summary.plannedNodeLevels = {20};
    if (moveResidualChild) {
        artifacts.preflight.summary.plannedStableChildRefCount = 1;
        artifacts.preflight.summary.stableChildToMoveCount = 1;
    }

    artifacts.input.updatedImage = minTree.reconstructionImage()->clone();
    artifacts.input.contraction.dualIsMaxtree = false;
    artifacts.input.localPixels = {0, 1, 2};

    if (includeResidualChildSnapshot) {
        noninc::LocalDualNodeSnapshot residualSnapshot;
        residualSnapshot.nodeId = residualChild;
        residualSnapshot.parentId = sourceOwner;
        residualSnapshot.parentInSnapshot = true;
        residualSnapshot.altitude = minTree.getAltitude(residualChild);
        residualSnapshot.numProperParts = minTree.getNumProperParts(residualChild);
        residualSnapshot.numChildren = minTree.getNumChildren(residualChild);
        residualSnapshot.directLocalProperParts = {0};
        residualSnapshot.supportLocalPixels = {0};
        residualSnapshot.supportTouchesModified = residualTouchesModified;
        residualSnapshot.supportTouchesPreserved = !residualTouchesModified;
        artifacts.input.preDualNodes = {residualSnapshot};
    }

    artifacts.materializationContext.contraction = artifacts.input.contraction;
    artifacts.materializationContext.preImage = minTree.reconstructionImage();
    artifacts.materializationContext.updatedImage = artifacts.input.updatedImage;
    artifacts.materializationContext.localPixels = artifacts.input.localPixels;
    artifacts.materializationContext.preDualNodes = artifacts.input.preDualNodes;

    artifacts.recompositionPlan.affectedLevels = {20};
    artifacts.recompositionPlan.nodes.resize(1);
    artifacts.recompositionPlan.nodes[0].level = 20;
    artifacts.recompositionPlan.nodes[0].supportPixels = std::move(planSupport);
    artifacts.recompositionPlan.nodes[0].directProperParts = {1};
    if (moveResidualChild) {
        artifacts.recompositionPlan.nodes[0].stableChildPreNodes = {residualChild};
    }
    artifacts.recompositionPlan.nodes[0].preDualOwners =
        {externalThroughSourceOwner ? sourceOwner : promotedParent};
    artifacts.recompositionPlan.nodes[0].parentPlanIndex = -1;
    return artifacts;
}

void test_mcpc_live_dry_run_accepts_empty_owner_with_all_children_moved() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const auto artifacts = make_empty_owner_with_child_dry_run_artifacts(
        minTree,
        true);
    const auto dryRun = noninc::diagnoseMcpcLiveMutationDryRun(
        minTree,
        artifacts);

    require(dryRun.evaluated,
            "empty-owner with child dry-run must mark evaluation");
    require(dryRun.accepted,
            "empty-owner whose children are all moved must enter the live D12 subdomain: " +
                dryRun.rejectReasonCode);
    require(dryRun.sourceOwnerToBeEmptiedCount == 1,
            "empty-owner with child dry-run must expose one emptied source owner");
    require(dryRun.rejectReasonCode.empty(),
            "accepted empty-owner with child dry-run must not report rejection");
}

void test_mcpc_live_dry_run_rejects_empty_owner_with_residual_child() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const auto artifacts = make_empty_owner_with_child_dry_run_artifacts(
        minTree,
        false);
    const auto dryRun = noninc::diagnoseMcpcLiveMutationDryRun(
        minTree,
        artifacts);

    require(dryRun.evaluated,
            "invalid empty-owner with child dry-run must mark evaluation");
    require(!dryRun.accepted,
            "empty-owner with a residual child must stay outside Dom(F)");
    require(dryRun.rejectReasonCode ==
                "d12_dry_run_empty_owner_with_children_residual_child_support_overlaps_plan",
            "residual child overlapping a planned node must expose the formal subcause");
}

void test_mcpc_live_dry_run_rejects_empty_owner_with_residual_child_missing_snapshot() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const auto artifacts = make_empty_owner_with_child_dry_run_artifacts(
        minTree,
        false,
        false,
        {1},
        false);
    const auto dryRun = noninc::diagnoseMcpcLiveMutationDryRun(
        minTree,
        artifacts);

    require(dryRun.evaluated,
            "missing residual snapshot dry-run must mark evaluation");
    require(!dryRun.accepted,
            "empty-owner residual promotion without a snapshot must stay outside Dom(F)");
    require(dryRun.rejectReasonCode ==
                "d12_dry_run_empty_owner_with_children_residual_child_snapshot_missing",
            "missing residual snapshot must expose the formal subcause");
}

void test_mcpc_live_dry_run_accepts_empty_owner_with_promotable_residual_child() {
    auto image = make_monotone_line_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    DynamicComponentTree minTree(image->clone(), false, adjacency);

    const auto artifacts = make_empty_owner_with_child_dry_run_artifacts(
        minTree,
        false,
        false,
        {1});
    const auto dryRun = noninc::diagnoseMcpcLiveMutationDryRun(
        minTree,
        artifacts);

    require(dryRun.evaluated,
            "promotable residual child dry-run must mark evaluation");
    require(dryRun.accepted,
            "empty-owner residual child disjoint from planned supports must enter Dom(F): " +
                dryRun.rejectReasonCode);
    require(dryRun.rejectReasonCode.empty(),
            "accepted promotable residual child dry-run must not report rejection");
}

void test_mcpc_live_mutation_dry_run_rejects_target_b_fixture() {
    auto image = make_p1_mcpc_target_b_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        1.0);
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenRebuildLocalRecomposition);

    for (float threshold : {0.0f, 1.0f, 5.0f}) {
        runner.applyThreshold(threshold);
    }

    const NodeId nodeId = find_p1_mcpc_fixture_merge_node(runner.maxTree(), 2);
    require(nodeId != InvalidNode,
            "P1-MCPC dry-run fixture must expose the merge before threshold 10");

    const auto artifacts = noninc::buildMcpcDirectPreflightArtifacts(
        runner.maxTree(),
        runner.minTree(),
        adjacency,
        nodeId);
    const auto d12 = noninc::diagnoseMcpcD12Domain(artifacts.preflight);
    require(d12.rejectReasonCode == "d12_exception_free_unproven",
            "P1-MCPC dry-run fixture must pass read-only D12 predicates");

    const auto dryRun =
        noninc::diagnoseMcpcLiveMutationDryRun(runner.minTree(), artifacts);
    require(dryRun.evaluated,
            "P1-MCPC live dry-run must mark evaluation");
    require(!dryRun.accepted,
            "target-b P1-MCPC fixture must remain outside the promoted live subdomain");
    require(!dryRun.rejectReasonCode.empty(),
            "rejected P1-MCPC live dry-run must expose a formal reason");
}

void test_mcpc_no_copy_candidate_diagnostics_require_read_only_preflight_conditions() {
    noninc::NonIncreasingMcpcDirectPreflightResult acceptedPreflight;
    acceptedPreflight.checked = true;
    acceptedPreflight.inputBuilt = true;
    acceptedPreflight.partitionPlanBuilt = true;
    acceptedPreflight.recompositionPlanBuilt = true;
    acceptedPreflight.planOracleAuditChecked = true;
    acceptedPreflight.planOracleFree = true;
    acceptedPreflight.materializationContextBuilt = true;
    acceptedPreflight.materializationOracleAuditChecked = true;
    acceptedPreflight.materializationOracleFree = true;
    acceptedPreflight.materializationInputsValid = true;
    acceptedPreflight.rollbackFreeSummaryBuilt = true;
    acceptedPreflight.passed = true;
    acceptedPreflight.directV1PrecheckEvaluated = true;
    acceptedPreflight.directV1PrecheckAccepted = true;

    const auto accepted =
        noninc::diagnoseMcpcNoCopyCandidate(acceptedPreflight);
    require(accepted.evaluated,
            "P1-MCPC no-copy candidate diagnostics must mark evaluation");
    require(accepted.accepted,
            "P1-MCPC no-copy candidate diagnostics must accept a clean read-only preflight");
    require(accepted.rejectionSummary.empty(),
            "accepted P1-MCPC no-copy candidate must not report rejection reasons");

    auto rejectedPreflight = acceptedPreflight;
    rejectedPreflight.directV1PrecheckAccepted = false;
    const auto rejectedByDirectV1 =
        noninc::diagnoseMcpcNoCopyCandidate(rejectedPreflight);
    require(rejectedByDirectV1.evaluated && !rejectedByDirectV1.accepted,
            "P1-MCPC no-copy candidate diagnostics must reject failed direct-v1 precheck");
    require(rejectedByDirectV1.rejectionSummary.find("direct-v1 precheck did not accept") !=
                std::string::npos,
            "P1-MCPC no-copy direct-v1 rejection must be reported");

    rejectedPreflight = acceptedPreflight;
    rejectedPreflight.directV1Precheck.predictedAuxiliaryNodeCount = 1;
    rejectedPreflight.directV1PrecheckAccepted =
        noninc::evaluateMcpcDirectV1PrecheckAcceptance(rejectedPreflight);
    const auto rejectedByAuxiliary =
        noninc::diagnoseMcpcNoCopyCandidate(rejectedPreflight);
    require(rejectedByAuxiliary.evaluated && !rejectedByAuxiliary.accepted,
            "P1-MCPC no-copy candidate diagnostics must reject predicted auxiliary connectors");
    require(rejectedByAuxiliary.rejectionSummary.find("direct-v1 precheck did not accept") !=
                std::string::npos,
            "P1-MCPC no-copy auxiliary-connector rejection must be reported through direct-v1");

    rejectedPreflight = acceptedPreflight;
    rejectedPreflight.materializationInputsValid = false;
    const auto rejectedByInputs =
        noninc::diagnoseMcpcNoCopyCandidate(rejectedPreflight);
    require(rejectedByInputs.evaluated && !rejectedByInputs.accepted,
            "P1-MCPC no-copy candidate diagnostics must reject invalid materialization inputs");
    require(rejectedByInputs.rejectionSummary.find("materialization inputs were not validated") !=
                std::string::npos,
            "P1-MCPC no-copy materialization-input rejection must be reported");

    rejectedPreflight = acceptedPreflight;
    rejectedPreflight.summary.stableChildCycleRiskCount = 1;
    const auto rejectedByCycle =
        noninc::diagnoseMcpcNoCopyCandidate(rejectedPreflight);
    require(rejectedByCycle.evaluated && !rejectedByCycle.accepted,
            "P1-MCPC no-copy candidate diagnostics must reject stable-child cycle risks");
    require(rejectedByCycle.rejectionSummary.find("stable child cycle risks") !=
                std::string::npos,
            "P1-MCPC no-copy cycle-risk rejection must be reported");
}

void test_merge_node_create_b_mcpc_real_partial_owner_missing_b_fixture_is_accepted() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setMcpcDirectPreflightDiagnosticsEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "5.2.10 P1-MCPC mode must preserve reconstruction equality");
        if (threshold != 10.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 110 &&
                event.parentLevel == 109 &&
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "5.2.10 crop must expose the P1-MCPC missing target-b event");
    require(targetEvent.mergeNodeLocalRawMissingTargetComponents == 1,
            "P1-MCPC mode must preserve the raw missing target-b diagnostic");
    require(targetEvent.mergeNodeLocalMissingTargetComponents == 0,
            "P1-MCPC mode must preserve the non-certified missing target-b diagnostic");
    require(targetEvent.mergeNodeLocalPartialOwnerComponents == 3,
            "P1-MCPC mode must preserve partial-owner diagnostics");
    require(targetEvent.mcpcDirectPreflightAttempted,
            "P1-MCPC mode must run direct preflight diagnostics when enabled");
    require(targetEvent.mcpcDirectPreflightPassed,
            "P1-MCPC direct preflight diagnostics must accept the recorded fixture");
    require(targetEvent.mcpcDirectPreflightPlanOracleFree &&
                targetEvent.mcpcDirectPreflightMaterializationOracleFree,
            "P1-MCPC direct preflight diagnostics must remain oracle-free");
    require(targetEvent.mcpcDirectPreflightMaterializationInputsValid,
            "P1-MCPC direct preflight diagnostics must validate materialization inputs");
    require(targetEvent.mcpcDirectPreflightSharedWithStrippedUpdate,
            "P1-MCPC direct preflight artifacts must be reused by the stripped copied trial");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckEvaluated,
            "P1-MCPC accepted trial must evaluate direct-v1 precheck diagnostics");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckAccepted,
            "P1-MCPC accepted fixture should satisfy the conservative direct-v1 precheck");
    require(targetEvent.mcpcDirectPreflightDirectV1Evaluated,
            "P1-MCPC accepted trial must evaluate direct-v1 acceptance diagnostics");
    require(targetEvent.mcpcDirectPreflightDirectV1Accepted,
            "P1-MCPC accepted fixture should satisfy the conservative direct-v1 criteria");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckAccepted ==
                targetEvent.mcpcDirectPreflightDirectV1Accepted,
            "P1-MCPC direct-v1 precheck and copied diagnostics should agree on the fixture");
    require(targetEvent.mcpcDirectPreflightSummaryComparedWithStrippedUpdate,
            "P1-MCPC accepted trial must compare direct summary against copied materialization");
    require(targetEvent.mcpcDirectPreflightSummaryMatchesStrippedUpdate,
            "P1-MCPC direct summary must match copied materialization on the recorded fixture");
    require(targetEvent.mcpcDirectPreflightSummaryMismatchCount == 0,
            "P1-MCPC direct summary comparison must report no mismatch on the recorded fixture");
    require(targetEvent.mcpcNoCopyCandidateEvaluated,
            "P1-MCPC accepted fixture must evaluate no-copy candidacy from read-only preflight");
    require(targetEvent.mcpcNoCopyCandidateAccepted,
            "P1-MCPC accepted fixture should be classified as a no-copy candidate");
    require(targetEvent.mcpcDirectCommitCandidateEvaluated,
            "P1-MCPC accepted fixture must evaluate direct commit candidacy diagnostics");
    require(targetEvent.mcpcDirectCommitCandidateAccepted,
            "P1-MCPC accepted fixture should be classified as a direct commit candidate");
    require(!targetEvent.mcpcDirectPreflightDirectV1StrictGateEnabled &&
                !targetEvent.mcpcDirectPreflightDirectV1StrictGateRejected,
            "P1-MCPC direct-v1 strict gate must stay disabled by default");
    require(targetEvent.mcpcDirectPreflightPlannedNodes ==
                targetEvent.strippedLocalRecompositionNodes,
            "P1-MCPC direct preflight must predict canonical recomposition node count");
    require(targetEvent.mcpcDirectPreflightPlannedDirectProperParts > 0,
            "P1-MCPC direct preflight must predict moved direct proper parts");
    require(targetEvent.mcpcDirectPreflightPredictedAuxiliaryNodes ==
                targetEvent.strippedCreatedConnectorNodes,
            "P1-MCPC direct preflight must predict observed auxiliary nodes on the recorded fixture");
    require(targetEvent.mcpcDirectPreflightPredictedStableChildBridgeNodes ==
                targetEvent.strippedCreatedStableChildBridgeNodes,
            "P1-MCPC direct preflight must predict observed stable-child bridges on the recorded fixture");
    require(targetEvent.mcpcDirectPreflightPredictedHiddenExternalParentConnectorNodes ==
                targetEvent.strippedCreatedHiddenExternalParentConnectorNodes,
            "P1-MCPC direct preflight must predict observed hidden external connectors on the recorded fixture");
    require(targetEvent.strippedLocalAttempted,
            "P1-MCPC mode must try stripped local recomposition after certification rejection");
    require(targetEvent.strippedLocalAccepted,
            "P1-MCPC stripped local recomposition must be accepted on the recorded fixture");
    require(targetEvent.mergeNodeLocalAccepted,
            "P1-MCPC mode must count the accepted stripped trial as a local merge-node acceptance");
    require(targetEvent.mergeNodeLocalDualMatchesPrimalAfterMerge,
            "P1-MCPC accepted trial must leave primal and dual reconstructions equal");
    require(!targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "P1-MCPC accepted trial must avoid the phase fallback");
    require(targetEvent.strippedLocalCreatedNodes > 0,
            "P1-MCPC accepted trial must materialize local recomposition nodes");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "P1-MCPC accepted fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_create_b_mcpc_strict_precheck_accepts_real_fixture() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setMcpcDirectV1StrictPrecheckEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "strict P1-MCPC precheck must preserve reconstruction equality");
        if (threshold != 10.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 110 &&
                event.parentLevel == 109 &&
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "strict P1-MCPC test must expose the recorded missing target-b event");
    require(targetEvent.mcpcDirectPreflightDirectV1StrictGateEnabled,
            "strict P1-MCPC mode must mark the direct-v1 gate as enabled");
    require(!targetEvent.mcpcDirectPreflightDirectV1StrictGateRejected,
            "strict P1-MCPC mode must not reject an accepted real fixture");
    require(targetEvent.mcpcDirectPreflightAttempted,
            "strict P1-MCPC mode must build preflight artifacts without diagnostic mode");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckEvaluated,
            "strict P1-MCPC mode must evaluate the direct-v1 precheck");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckAccepted,
            "strict P1-MCPC real fixture must pass the direct-v1 precheck");
    require(targetEvent.strippedLocalAttempted && targetEvent.strippedLocalAccepted,
            "strict P1-MCPC accepted precheck must continue through copied local recomposition");
    require(targetEvent.mcpcDirectPreflightSharedWithStrippedUpdate,
            "strict P1-MCPC accepted precheck must reuse preflight artifacts in the copied trial");
    require(targetEvent.mcpcDirectPreflightSummaryComparedWithStrippedUpdate &&
                targetEvent.mcpcDirectPreflightSummaryMatchesStrippedUpdate,
            "strict P1-MCPC accepted precheck must preserve direct-summary/copy agreement");
    require(targetEvent.mcpcDirectPreflightSummaryMismatchCount == 0,
            "strict P1-MCPC accepted fixture must have no direct-summary mismatch");
    require(targetEvent.mcpcNoCopyCandidateEvaluated &&
                targetEvent.mcpcNoCopyCandidateAccepted,
            "strict P1-MCPC accepted fixture should be classified as a no-copy candidate");
    require(targetEvent.mcpcDirectCommitCandidateEvaluated &&
                targetEvent.mcpcDirectCommitCandidateAccepted,
            "strict P1-MCPC accepted fixture should be classified as a direct commit candidate");
    require(!targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "strict P1-MCPC accepted fixture must avoid phase fallback");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "strict P1-MCPC fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_create_b_mcpc_validated_direct_commit_accepts_real_fixture() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setMcpcValidatedDirectCommitEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "validated direct P1-MCPC mode must preserve reconstruction equality");
        if (threshold != 10.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 110 &&
                event.parentLevel == 109 &&
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "validated direct P1-MCPC test must expose the recorded missing target-b event");
    require(targetEvent.mcpcValidatedDirectCommitEnabled,
            "validated direct P1-MCPC mode must mark direct commit as enabled");
    require(targetEvent.mcpcValidatedDirectCommitAttempted,
            "validated direct P1-MCPC mode must attempt live direct commit after copy validation");
    require(targetEvent.mcpcValidatedDirectCommitAccepted,
            "validated direct P1-MCPC mode must accept the live direct commit on the fixture");
    require(!targetEvent.mcpcValidatedDirectCommitFellBackToCopy,
            "validated direct P1-MCPC mode must not fall back to the copied result on the fixture");
    require(targetEvent.mcpcValidatedDirectCommitDualMatchesPrimal,
            "validated direct P1-MCPC live commit must leave dual/primal reconstructions equal");
    require(targetEvent.mcpcValidatedDirectCommitMatchesCopiedOutput,
            "validated direct P1-MCPC live commit must match the copied validation output");
    require(targetEvent.mergeNodeLocalDirectCommit && targetEvent.strippedDirectCommit,
            "validated direct P1-MCPC accepted fixture must record the live direct commit path");
    require(targetEvent.mcpcDirectCommitCandidateEvaluated &&
                targetEvent.mcpcDirectCommitCandidateAccepted,
            "validated direct P1-MCPC fixture must pass candidate diagnostics before live commit");
    require(targetEvent.mcpcNoCopyCandidateEvaluated &&
                targetEvent.mcpcNoCopyCandidateAccepted,
            "validated direct P1-MCPC fixture should be classified as a no-copy candidate");
    require(!targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "validated direct P1-MCPC fixture must avoid phase fallback");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "validated direct P1-MCPC fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_create_b_mcpc_no_copy_direct_commit_accepts_real_fixture() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setMcpcNoCopyDirectCommitEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "no-copy direct P1-MCPC mode must preserve reconstruction equality");
        if (threshold != 10.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 110 &&
                event.parentLevel == 109 &&
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "no-copy direct P1-MCPC test must expose the recorded missing target-b event");
    require(targetEvent.mcpcNoCopyCandidateEvaluated &&
                targetEvent.mcpcNoCopyCandidateAccepted,
            "no-copy direct P1-MCPC fixture must pass read-only no-copy candidacy");
    require(targetEvent.mcpcNoCopyDirectCommitEnabled,
            "no-copy direct P1-MCPC mode must mark no-copy direct commit as enabled");
    require(targetEvent.mcpcNoCopyDirectCommitAttempted,
            "no-copy direct P1-MCPC mode must attempt live direct commit before copy");
    require(targetEvent.mcpcNoCopyDirectCommitAccepted,
            "no-copy direct P1-MCPC mode must accept the live direct commit on the fixture");
    require(!targetEvent.mcpcNoCopyDirectCommitFellBackToPhaseRebuild,
            "no-copy direct P1-MCPC fixture must not fall back to phase rebuild");
    require(targetEvent.mcpcNoCopyDirectCommitDualMatchesPrimal,
            "no-copy direct P1-MCPC commit must leave dual/primal reconstructions equal");
    require(targetEvent.mergeNodeLocalDirectCommit && targetEvent.strippedDirectCommit,
            "no-copy direct P1-MCPC accepted fixture must record direct commit flags");
    require(targetEvent.strippedLocalAttempted && targetEvent.strippedLocalAccepted,
            "no-copy direct P1-MCPC accepted fixture must record accepted stripped recomposition");
    require(targetEvent.strippedLocalCopyMicros == 0,
            "no-copy direct P1-MCPC accepted fixture must avoid the copied trial");
    require(!targetEvent.mcpcDirectCommitCandidateEvaluated,
            "no-copy direct P1-MCPC fixture should not require post-copy commit candidacy");
    require(!targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "no-copy direct P1-MCPC fixture must avoid phase fallback");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "no-copy direct P1-MCPC fixture sequence must match the naive rebuild baseline");
}

void test_merge_node_create_b_mcpc_no_copy_direct_commit_falls_back_on_thrown_failure() {
    auto image = load_dat_crop("dat/misc256/5.2.10.png", 32, 32);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 5.0f, 10.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setMcpcNoCopyDirectCommitEnabled(true);

    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.testThrowAfterMaterializationBeforePrimalMerge = true;
    runner.setStrippedLocalUpdateOptions(options);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "no-copy direct fallback P1-MCPC mode must preserve reconstruction equality");
        if (threshold != 10.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeLevel == 110 &&
                event.parentLevel == 109 &&
                event.mergeNodeLocalPartialOwnerAndMissingTargetComponents == 1) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "no-copy direct fallback P1-MCPC test must expose the recorded missing target-b event");
    require(targetEvent.mcpcNoCopyCandidateEvaluated,
            "fallback P1-MCPC fixture must classify no-copy candidacy");
    require(targetEvent.mcpcNoCopyCandidateAccepted,
            "fallback P1-MCPC fixture must still classify as no-copy candidate");
    require(targetEvent.mcpcNoCopyDirectCommitEnabled,
            "fallback P1-MCPC mode must mark no-copy direct commit as enabled");
    require(targetEvent.mcpcNoCopyDirectCommitAttempted,
            "fallback P1-MCPC mode must attempt live direct commit before copy");
    require(!targetEvent.mcpcNoCopyDirectCommitAccepted,
            "fallback P1-MCPC mode must reject the live direct commit path");
    require(targetEvent.mcpcNoCopyDirectCommitFailureSummary.find(
                "forced test failure") != std::string::npos,
            "fallback P1-MCPC mode must expose the injected test failure reason");
    require(targetEvent.mcpcNoCopyDirectCommitFellBackToPhaseRebuild,
            "fallback P1-MCPC mode must record phase rebuild fallback");
    require(targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "fallback P1-MCPC mode must execute merge-node phase fallback");
    require(!targetEvent.mcpcNoCopyDirectCommitDualMatchesPrimal,
            "fallback P1-MCPC test should mark dual/primal mismatch after failed direct path");
    require(!targetEvent.mergeNodeLocalDirectCommit && !targetEvent.strippedDirectCommit,
            "fallback P1-MCPC test must not report a direct commit success");
    require(!targetEvent.mergeNodeLocalAccepted,
            "fallback P1-MCPC test must not accept merge-node local recomposition");
    require(!targetEvent.strippedLocalAccepted,
            "fallback P1-MCPC test must not accept stripped local recomposition");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "no-copy direct fallback P1-MCPC sequence must match the naive rebuild baseline");
}

void test_merge_node_create_b_mcpc_real_hidden_external_fixture_predicts_auxiliary_connectors() {
    auto image = load_dat_crop("dat/misc256/7.1.10.png", 128, 128);
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    const std::vector<float> thresholds = {0.0f, 1.0f, 2.0f, 5.0f, 10.0f, 20.0f};

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::MergeNodeCertifiedCreateBMcpcThenPhaseRebuild);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);
    runner.setMcpcDirectPreflightDiagnosticsEnabled(true);

    bool foundTargetEvent = false;
    noninc::NonIncreasingMergeEvent targetEvent;
    for (float threshold : thresholds) {
        runner.applyThreshold(threshold);
        require(!runner.hasFirstInvariantMismatch(),
                "7.1.10 hidden-external P1-MCPC fixture must preserve reconstruction equality");
        if (threshold != 20.0f) {
            continue;
        }
        for (const auto &event : runner.lastMergeEvents()) {
            if (event.phase == noninc::NonIncreasingPhase::MaxToMin &&
                event.nodeId == 74 &&
                event.parentId == 70 &&
                event.nodeLevel == 110 &&
                event.parentLevel == 108) {
                targetEvent = event;
                foundTargetEvent = true;
                break;
            }
        }
    }

    require(foundTargetEvent,
            "7.1.10 crop must expose the recorded hidden-external P1-MCPC event");
    require(targetEvent.componentCount == 47 && targetEvent.totalComponentPixels == 130,
            "7.1.10 hidden-external fixture must preserve the recorded multicomponent event");
    require(targetEvent.mcpcDirectPreflightAttempted && targetEvent.mcpcDirectPreflightPassed,
            "hidden-external fixture must pass P1-MCPC direct preflight diagnostics");
    require(targetEvent.mcpcDirectPreflightSharedWithStrippedUpdate,
            "hidden-external fixture must reuse P1-MCPC direct preflight artifacts");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckEvaluated,
            "hidden-external fixture must evaluate direct-v1 precheck diagnostics");
    require(!targetEvent.mcpcDirectPreflightDirectV1PrecheckAccepted,
            "hidden-external fixture must fail direct-v1 precheck because it needs auxiliary connectors");
    require(targetEvent.mcpcDirectPreflightDirectV1Evaluated,
            "hidden-external fixture must evaluate direct-v1 acceptance diagnostics");
    require(!targetEvent.mcpcDirectPreflightDirectV1Accepted,
            "hidden-external fixture must fail direct-v1 copied diagnostics because it creates auxiliary connectors");
    require(targetEvent.mcpcDirectPreflightDirectV1PrecheckAccepted ==
                targetEvent.mcpcDirectPreflightDirectV1Accepted,
            "hidden-external direct-v1 precheck and copied diagnostics should agree");
    require(targetEvent.mcpcDirectPreflightSummaryComparedWithStrippedUpdate,
            "hidden-external fixture must compare direct summary against copied materialization");
    require(targetEvent.mcpcDirectPreflightSummaryMatchesStrippedUpdate,
            "hidden-external direct summary must match copied materialization");
    require(targetEvent.mcpcDirectPreflightSummaryMismatchCount == 0,
            "hidden-external direct summary comparison must report no mismatch");
    require(targetEvent.mcpcNoCopyCandidateEvaluated,
            "hidden-external fixture must evaluate no-copy candidacy from read-only preflight");
    require(!targetEvent.mcpcNoCopyCandidateAccepted,
            "hidden-external fixture must not be classified as a no-copy candidate");
    require(targetEvent.mcpcDirectCommitCandidateEvaluated,
            "hidden-external fixture must evaluate direct commit candidacy diagnostics");
    require(!targetEvent.mcpcDirectCommitCandidateAccepted,
            "hidden-external fixture must not be classified as a direct commit candidate");
    require(targetEvent.strippedLocalAccepted && targetEvent.mergeNodeLocalAccepted,
            "hidden-external fixture must accept the stripped copied trial");
    require(!targetEvent.mergeNodeLocalFallbackToPhaseRebuild,
            "hidden-external fixture must avoid phase fallback after stripped acceptance");
    require(targetEvent.strippedCreatedHiddenExternalParentConnectorNodes == 3,
            "hidden-external fixture must materialize the recorded hidden external connectors");
    require(targetEvent.strippedCreatedStableChildBridgeNodes == 0,
            "hidden-external fixture must isolate hidden external connectors from stable-child bridges");
    require(targetEvent.mcpcDirectPreflightPredictedAuxiliaryNodes ==
                targetEvent.strippedCreatedConnectorNodes,
            "P1-MCPC direct preflight must predict all auxiliary connectors on hidden-external fixture");
    require(targetEvent.mcpcDirectPreflightPredictedHiddenExternalParentConnectorNodes ==
                targetEvent.strippedCreatedHiddenExternalParentConnectorNodes,
            "P1-MCPC direct preflight must predict hidden external connectors on the real fixture");
    require(targetEvent.mcpcDirectPreflightPredictedStableChildBridgeNodes == 0,
            "P1-MCPC direct preflight must not predict stable-child bridges on the hidden-external fixture");

    const auto output = runner.reconstructImage();
    const auto naiveOutput = noninc::applyNaiveGrayVarianceMergeSequence(image, adjacency, thresholds);
    require(output->isEqual(naiveOutput),
            "hidden-external P1-MCPC fixture sequence must match the naive rebuild baseline");
}

void test_stripped_local_connector_crop_releases_empty_created_bridges() {
    auto image = make_icdar_connector_crop_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition);
    noninc::NonIncreasingStrippedLocalUpdateOptions options;
    options.checkInputPreReconstruction = false;
    options.buildFullUpdatedImage = false;
    options.checkPreMaterializationReconstruction = false;
    options.checkMaterializationTargetReconstruction = false;
    options.checkNoEmptyLiveNodesAfterCanonicalization = false;
    options.checkFinalDualPrimalReconstruction = false;
    options.releaseEmptyNodesAfterPrimalMerge = false;
    runner.setStrippedLocalUpdateOptions(options);
    runner.setStrippedCommitMode(noninc::NonIncreasingStrippedCommitMode::Direct);
    runner.setMergeEventCollectionEnabled(true);
    runner.setStopOnFirstInvariantMismatchEnabled(true);

    runner.applyThreshold(10.0f);
    require(!runner.hasFirstInvariantMismatch(),
            "connector crop stripped local update must preserve reconstruction equality");

    int strippedAttempts = 0;
    int hiddenCreatedNodes = 0;
    int emptyCreatedNodes = 0;
    int phaseFallbacks = 0;
    for (const auto &event : runner.lastMergeEvents()) {
        if (!event.strippedLocalAttempted) {
            continue;
        }
        ++strippedAttempts;
        hiddenCreatedNodes +=
            std::max(0, event.strippedLocalCreatedNodes - event.strippedLocalRecompositionNodes);
        emptyCreatedNodes += event.strippedEmptyCreatedNodesAfterCanonicalization;
        if (event.strippedLocalFallbackToPhaseRebuild) {
            ++phaseFallbacks;
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput =
        noninc::applyNaiveGrayVarianceMergeThreshold(image, adjacency, 10.0f);

    require(strippedAttempts > 0,
            "connector crop must exercise stripped local updates");
    require(hiddenCreatedNodes > 0,
            "connector crop must exercise auxiliary bridge/connector creation");
    require(emptyCreatedNodes == 0,
            "local canonicalization must release created bridge nodes without proper parts");
    require(phaseFallbacks == 0,
            "connector crop must not need phase fallback after empty bridge release");
    require(output->isEqual(naiveOutput),
            "connector crop stripped local update must match the naive rebuild baseline");
}

NodeId find_first_non_root_node(const DynamicComponentTree &tree) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (nodeId != tree.getRoot() && tree.isAlive(nodeId)) {
            return nodeId;
        }
    }
    return InvalidNode;
}

NodeId find_node_by_level_and_support_size(const DynamicComponentTree &tree,
                                           int level,
                                           int supportSize) {
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (nodeId == tree.getRoot() || !tree.isAlive(nodeId)) {
            continue;
        }
        if (tree.getAltitude(nodeId) != level) {
            continue;
        }
        if (static_cast<int>(noninc::collectNodeSupportPixels(tree, nodeId).size()) == supportSize) {
            return nodeId;
        }
    }
    return InvalidNode;
}

struct PhaseCandidateCropNonCanonicalNodes {
    NodeId localOnly = InvalidNode;
    NodeId highChild = InvalidNode;
};

PhaseCandidateCropNonCanonicalNodes insert_phase_candidate_crop_local_only_intermediate(
    DynamicComponentTree &tree,
    ImageUInt8Ptr image) {
    const NodeId highChild = find_node_by_level_and_support_size(tree, 223, 16);
    require(highChild != InvalidNode,
            "phase-candidate crop must expose the level-223 support-16 child");

    NodeId attachParent = tree.getNodeParent(highChild);
    if (attachParent != InvalidNode &&
        attachParent != highChild &&
        tree.isAlive(attachParent) &&
        tree.getAltitude(attachParent) == 220 &&
        !tree.isRoot(attachParent)) {
        attachParent = tree.getNodeParent(attachParent);
    }
    require(attachParent != InvalidNode &&
                attachParent != highChild &&
                tree.isAlive(attachParent) &&
                tree.getAltitude(attachParent) < 220,
            "phase-candidate crop must have an external parent below level 220");

    const NodeId localOnly = tree.createDetachedNode(220);
    tree.attachNode(attachParent, localOnly);
    tree.moveNode(highChild, localOnly);

    const std::vector<PixelId> directLevel220Pixels = {11, 19};
    for (PixelId pixelId : directLevel220Pixels) {
        require(static_cast<int>((*image)[pixelId]) == 220,
                "phase-candidate crop setup must move only level-220 pixels");
        const NodeId sourceOwner = tree.getSmallestComponent(pixelId);
        require(sourceOwner != InvalidNode &&
                    sourceOwner != localOnly &&
                    tree.isAlive(sourceOwner) &&
                    tree.getAltitude(sourceOwner) == 220,
                "phase-candidate crop level-220 pixel must be directly owned by a level-220 source");
        tree.moveProperPart(localOnly, sourceOwner, pixelId);
    }

    return {localOnly, highChild};
}

void insert_redundant_same_level_bridge(DynamicComponentTree &tree, NodeId childNode) {
    const NodeId parentNode = tree.getNodeParent(childNode);
    require(parentNode != InvalidNode && parentNode != childNode,
            "test bridge target must have a real parent");
    const NodeId bridgeNode = tree.createDetachedNode(tree.getAltitude(childNode));
    tree.attachNode(parentNode, bridgeNode);
    tree.moveNode(childNode, bridgeNode);
}

void test_phase_candidate_crop_local_only_intermediate_is_canonicalized() {
    auto image = make_icdar_phase_candidate_crop_fixture_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);
    DynamicComponentTree localTree(image, true, adjacency);
    DynamicComponentTree rebuiltTree(image, true, adjacency);

    const auto nodes =
        insert_phase_candidate_crop_local_only_intermediate(localTree, image);

    require(localTree.reconstructionImage()->isEqual(image),
            "inserted local-only intermediate must preserve reconstruction");

    auto before = noninc::compareNonIncreasingPhaseAgainstRebuild(
        localTree,
        rebuiltTree,
        10.0f,
        noninc::NonIncreasingPhase::MinToMax);
    require(!before.candidateSetMatchesRebuild,
            "inserted intermediate must create a local-only phase candidate");
    require(!before.candidateScanOrderMatchesRebuild,
            "local-only phase candidate must also break candidate scan order");
    require(before.releasePredicateEvaluated,
            "phase-candidate mismatch must evaluate the release predicate on the focus");
    require(before.releaseFocusAvailable,
            "phase-candidate mismatch must provide a release focus");
    require(!before.releasePredicateSatisfied,
            "release predicate must reject the focused local-only intermediate");
    require(before.releaseH1ReconstructionMatches,
            "release H1 should still pass for the reconstruction-equivalent intermediate");
    require(!before.releaseH2CandidateSetMatches,
            "release H2 must fail for the local-only candidate");
    require(!before.releaseH3EnvelopeStructureMatches,
            "release H3 must fail for the focused local-only intermediate");
    require(before.firstLocalOnlyContext.altitude == 220,
            "phase-candidate crop local-only node must be the level-220 intermediate");
    require(before.envelopeStructureCompared,
            "phase-candidate mismatch must compare the focused structural envelope");
    require(!before.envelopeStructureMatchesRebuild,
            "phase-candidate crop local-only intermediate must differ structurally from rebuild");
    require(before.releaseFocusPixelCount > 0,
            "phase release diagnostic must expose the focus pixel set F");
    require(before.releaseFocusPixelCount ==
                static_cast<int>(before.releaseFocusPixels.size()),
            "phase release focus pixel count must match the stored focus pixels");
    require(before.releaseClosureNodeCount > 0,
            "phase release diagnostic must expose the local closure node set");
    require(before.releaseClosurePixelCount ==
                static_cast<int>(before.releaseClosurePixels.size()),
            "phase release closure pixel count must match the stored closure pixels");
    require(before.releaseClosurePixelCount >= before.releaseFocusPixelCount,
            "the incident closure support must cover the focus pixels");

    const std::vector<NodeId> focusNodes = {nodes.localOnly, nodes.highChild};
    const int changes = noninc::canonicalizeAfterLocalRecomposition(
        localTree,
        adjacency.get(),
        std::span<const NodeId>(focusNodes),
        std::span<const NodeId>{});
    require(changes > 0,
            "canonicalization must change the local-only intermediate topology");

    DynamicComponentTree rebuiltAfter(image, true, adjacency);
    auto after = noninc::compareNonIncreasingPhaseAgainstRebuild(
        localTree,
        rebuiltAfter,
        10.0f,
        noninc::NonIncreasingPhase::MinToMax);
    require(after.candidateSetMatchesRebuild,
            "canonicalization must remove the crop local-only candidate");
    require(after.candidateScanOrderMatchesRebuild,
            "canonicalization must restore the candidate scan order");
    require(!after.releasePredicateEvaluated,
            "release predicate remains unevaluated when no phase focus is needed");
}

void test_phase_candidate_guard_canonicalizes_reconstruction_equivalent_noncanonical_tree() {
    auto image = make_demo_image();
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), 1.0);

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition);
    runner.setPhaseCandidateEquivalenceGuardEnabled(true);
    runner.setPhaseOracleComparisonEnabled(true);

    const NodeId targetNode = find_first_non_root_node(runner.maxTree());
    require(targetNode != InvalidNode,
            "demo max-tree must have a non-root node for the phase-candidate guard test");
    insert_redundant_same_level_bridge(runner.maxTree(), targetNode);

    require(runner.maxTree().reconstructionImage()->isEqual(runner.minTree().reconstructionImage()),
            "redundant same-level bridge must preserve image reconstruction before the guard");

    runner.applyThreshold(100000.0f);

    bool sawInitialCandidateMismatch = false;
    for (const auto &comparison : runner.lastPhaseOracleComparisons()) {
        if (comparison.beforeMaxToMin && !comparison.candidateSetMatchesRebuild) {
            sawInitialCandidateMismatch = true;
            break;
        }
    }

    const auto output = runner.reconstructImage();
    const auto naiveOutput =
        noninc::applyNaiveGrayVarianceMergeThreshold(image, adjacency, 100000.0f);

    require(sawInitialCandidateMismatch,
            "phase oracle must expose the reconstruction-equivalent candidate mismatch");
    require(runner.lastPhaseCandidateGuardCanonicalizationChanges() > 0,
            "phase candidate guard must canonicalize the noncanonical phase tree before rebuilding");
    require(runner.lastPhaseCandidateGuardRebuilds() == 0,
            "phase candidate guard should avoid rebuild when canonicalization restores candidate equivalence");
    require(output->isEqual(naiveOutput),
            "phase candidate guard must preserve the naive threshold result after canonicalizing");
}

void test_phase_candidate_guard_canonicalizes_crop_local_only_intermediate() {
    auto image = make_icdar_phase_candidate_crop_fixture_image();

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition);
    runner.setPhaseCandidateEquivalenceGuardEnabled(true);
    runner.setPhaseOracleComparisonEnabled(true);
    runner.setMaxAppliedMergeEventsPerThreshold(1);

    insert_phase_candidate_crop_local_only_intermediate(runner.maxTree(), image);

    require(runner.maxTree().reconstructionImage()->isEqual(runner.minTree().reconstructionImage()),
            "phase-candidate crop intermediate must preserve reconstruction before the guard");

    runner.applyThreshold(10.0f);

    bool sawInitialCandidateMismatch = false;
    for (const auto &comparison : runner.lastPhaseOracleComparisons()) {
        if (comparison.beforeMaxToMin &&
            !comparison.candidateSetMatchesRebuild &&
            comparison.firstLocalOnlyContext.altitude == 220) {
            sawInitialCandidateMismatch = true;
            break;
        }
    }

    require(sawInitialCandidateMismatch,
            "phase oracle must expose the injected crop local-only candidate before repair");
    require(runner.lastPhaseCandidateGuardCanonicalizationChanges() > 0,
            "phase candidate guard must canonicalize the injected crop local-only intermediate");
    require(runner.lastPhaseCandidateGuardRebuilds() == 0,
            "phase candidate guard should not rebuild when crop canonicalization restores candidate equivalence");
    require(runner.maxTree().reconstructionImage()->isEqual(runner.minTree().reconstructionImage()),
            "phase candidate guard crop repair must preserve primal/dual reconstruction equality");
}

void test_phase_release_rejection_stop_captures_incident_closure() {
    auto image = make_icdar_phase_candidate_crop_fixture_image();

    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::StrippedThenPhaseRebuildLocalRecomposition);
    runner.setPhaseOracleComparisonEnabled(true);
    runner.setStopOnFirstPhaseReleaseRejectionEnabled(true);

    insert_phase_candidate_crop_local_only_intermediate(runner.maxTree(), image);

    runner.applyThreshold(10.0f);

    require(runner.lastPhaseReleaseRejectionStopReached(),
            "diagnostic runner should stop at the first focused Release_t rejection");
    require(!runner.lastPhaseOracleComparisons().empty(),
            "diagnostic stop must retain the phase-oracle comparison that triggered it");

    const auto &comparison = runner.lastPhaseOracleComparisons().front();
    require(comparison.releasePredicateEvaluated,
            "incident closure diagnostic must evaluate the release predicate");
    require(comparison.releaseFocusAvailable,
            "incident closure diagnostic must expose a focused comparison");
    require(!comparison.releasePredicateSatisfied,
            "injected local-only candidate should be rejected by Release_t");
    require(comparison.releaseFocusPixelCount > 0,
            "incident closure diagnostic must retain focus pixels");
    require(comparison.releaseClosureNodeCount > 0,
            "incident closure diagnostic must retain closure nodes");
    require(comparison.releaseClosurePixelCount >= comparison.releaseFocusPixelCount,
            "incident closure support must cover the focus");
}

void test_oracle_guided_event_limit_stops_threshold() {
    auto image = make_demo_image();
    noninc::NonIncreasingCasfKernel runner(
        image,
        1.0,
        noninc::NonIncreasingDualMaintenance::OracleGuidedLocalRecomposition);
    runner.setMergeEventCollectionEnabled(true);
    runner.setMaxAppliedMergeEventsPerThreshold(2);

    runner.applyThreshold(200.0f);

    require(runner.lastMergeEventLimitReached(),
            "oracle-guided runner must report when the merge-event limit stops a threshold");
    require(runner.lastMergeEvents().size() == 2,
            "oracle-guided runner must stop after the configured number of applied merge events");
}

void test_empty_threshold_sequence_is_noop() {
    auto image = make_demo_image();
    noninc::NonIncreasingCasfKernel runner(image, 1.0);

    require(runner.filter({})->isEqual(image),
            "empty non-increasing threshold sequence must preserve the input image");
}

} // namespace

int main() {
    try {
        test_root_gray_variance_matches_image_variance();
        test_root_box_area_perimeter_ratio_matches_direct_geometry();
        test_node_storage_reallocation_counters_track_capacity_growth();
        test_box_area_perimeter_clean_update_matches_naive_sequence();
        test_box_area_perimeter_clean_update_keeps_local_attribute_buffers_exact();
        test_candidate_key_matches_full_support_order_on_clean_bbox_path();
        test_partial_f_default_records_formal_acceptance_without_fallback();
        test_partial_f_rejects_invalid_node_without_mutation();
        test_partial_f_mcpc_formal_accounting_without_fallback();
        test_partial_f_accepts_mcpc_live_subdomain_without_fallback();
        test_partial_f_accepts_parent_link_bridge_live_subdomain_without_fallback();
        test_partial_f_accepts_hidden_support_redundant_live_subdomain_without_fallback();
        test_partial_f_omits_stable_child_bridge_before_live_mutation();
        test_partial_f_accepts_stable_child_move_live_subdomain_without_fallback();
        test_partial_f_accepts_materializable_missing_target_b();
        test_partial_f_component_order_independence_modulo_isomorphism();
        test_partial_f_replay_is_deterministic();
        test_mcpc_d12_domain_diagnostics_reject_formal_predicates();
        test_mcpc_live_dry_run_accepts_empty_owner_with_all_children_moved();
        test_mcpc_live_dry_run_rejects_empty_owner_with_residual_child();
        test_mcpc_live_dry_run_rejects_empty_owner_with_residual_child_missing_snapshot();
        test_mcpc_live_dry_run_accepts_empty_owner_with_promotable_residual_child();
        test_mcpc_live_mutation_dry_run_rejects_target_b_fixture();
        test_mcpc_direct_v1_precheck_rejects_synthetic_nonoracle_failures();
        test_mcpc_no_copy_candidate_diagnostics_require_read_only_preflight_conditions();
        test_merge_node_create_b_mcpc_real_partial_owner_missing_b_fixture_is_accepted();
        test_empty_threshold_sequence_is_noop();
    } catch (const std::exception &e) {
        std::cerr << "details_non_increasing_casf_kernel_unit_tests: FAIL\n" << e.what() << "\n";
        return 1;
    }

    std::cout << "details_non_increasing_casf_kernel_unit_tests: OK\n";
    return 0;
}
