#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "dynamic_casf_apply_common.hpp"

#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1CasfGate.hpp"
#include "../include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1Updater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Problem2.hpp"

namespace common = dynamic_casf_apply_common;
namespace problem1_debug = debug::problem1;
namespace problem2_debug = debug::problem2;

namespace {

enum class BackendKind {
    UpdateTree,
    ExplicitRegionDirect,
    Problem1FastGateExplicit,
    Problem2ExplicitSweep,
};

struct BackendRunResult {
    ImageUInt8Ptr output;
    long long totalUs = 0;
    long long updateUs = 0;
    long long phase1UpdateUs = 0;
    long long phase2UpdateUs = 0;
    long long attributeRefreshUs = 0;
    long long problem2BatchWallUs = 0;
    long long problem2EventTotalUs = 0;
    long long problem2EventSetupUs = 0;
    long long collectGateUs = 0;
    long long problem2RunUpdateTreeWallUs = 0;
    long long bucketBuildUs = 0;
    long long targetNodeCreationUs = 0;
    long long metadataUs = 0;
    long long levelPrepUs = 0;
    long long sweepUs = 0;
    long long finalizeUs = 0;
    long long finalReconnectUs = 0;
    long long releaseEmptyUs = 0;
    long long sameLevelCollapseUs = 0;
    long long problem2PruneUs = 0;
    long long problem2EventUnaccountedUs = 0;
    long long problem2KernelUnaccountedUs = 0;
    long long problem2BatchOverheadUs = 0;
    long long regionPixels = 0;
    long long mergeNodes = 0;
    long long mergeLevels = 0;
    long long frontierNodes = 0;
    long long sameLevelEdgesCollapsed = 0;
    int events = 0;
    int skippedNodes = 0;
    int fastGateReady = 0;
    int fastGateRejected = 0;
    int localApplied = 0;
    int localRejected = 0;
    int fallbackToUpdateTree = 0;
    int attributeBufferRebinds = 0;
    int maxRegionPixels = 0;
    int maxBucketSize = 0;
    int targetNodeCreatedEvents = 0;
};

std::vector<int> parseThresholdCsv(const std::string &text) {
    std::vector<int> thresholds;
    std::stringstream stream(text);
    std::string token;
    while (std::getline(stream, token, ',')) {
        if (token.empty()) {
            continue;
        }
        thresholds.push_back(std::stoi(token));
    }
    return thresholds;
}

std::vector<PixelId> collectSubtreePixels(const DynamicComponentTree &tree,
                                          NodeId subtreeRoot) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.getNodeSubtree(subtreeRoot)) {
        for (PixelId pixelId : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(pixelId);
        }
    }
    return pixels;
}

bool isValidPruneRoot(const DynamicComponentTree &tree, NodeId rootSubtree) {
    return rootSubtree != InvalidNode &&
           rootSubtree != tree.getRoot() &&
           tree.isNode(rootSubtree) &&
           tree.isAlive(rootSubtree);
}

int countValidPruneRoots(const DynamicComponentTree &tree,
                         const std::vector<NodeId> &nodesToPrune) {
    int count = 0;
    for (NodeId rootSubtree : nodesToPrune) {
        if (isValidPruneRoot(tree, rootSubtree)) {
            ++count;
        }
    }
    return count;
}

class UpdateBackendRunner {
private:
    ImageUInt8Ptr inputImage_;
    AdjacencyRelationPtr adj_;
    DynamicComponentTree maxTree_;
    DynamicComponentTree minTree_;
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust_;
    std::unique_ptr<DynamicAttributeComputer> maxAttributeComputer_;
    std::unique_ptr<DynamicAttributeComputer> minAttributeComputer_;
    Attribute attribute_ = AREA;
    std::vector<float> maxArea_;
    std::vector<float> minArea_;
    BackendKind backend_ = BackendKind::UpdateTree;
    bool profileProblem2ExplicitSweep_ = false;

    void rebindAttributeBuffers() {
        adjust_.setAttributeComputer(*minAttributeComputer_,
                                     *maxAttributeComputer_,
                                     std::span<float>(minArea_),
                                     std::span<float>(maxArea_));
    }

    long long refreshAllAttributes() {
        const auto start = std::chrono::steady_clock::now();
        maxArea_ = common::computeAttributeVector(&maxTree_, attribute_);
        minArea_ = common::computeAttributeVector(&minTree_, attribute_);
        rebindAttributeBuffers();
        return common::elapsedMicros(start, std::chrono::steady_clock::now());
    }

    long long ensureAreaBuffersHaveSpareNodeSlot(BackendRunResult &result) {
        const auto start = std::chrono::steady_clock::now();
        bool resized = false;
        constexpr std::size_t spareSlots = 16;

        const std::size_t requiredMax =
            static_cast<std::size_t>(maxTree_.getNumInternalNodeSlots()) + 1;
        if (maxArea_.size() < requiredMax) {
            maxArea_.resize(requiredMax + spareSlots, 0.0f);
            resized = true;
        }

        const std::size_t requiredMin =
            static_cast<std::size_t>(minTree_.getNumInternalNodeSlots()) + 1;
        if (minArea_.size() < requiredMin) {
            minArea_.resize(requiredMin + spareSlots, 0.0f);
            resized = true;
        }

        if (resized) {
            rebindAttributeBuffers();
            ++result.attributeBufferRebinds;
        }

        return common::elapsedMicros(start, std::chrono::steady_clock::now());
    }

    long long pruneMaxTreeAndUpdateMinTreeExplicit(std::vector<NodeId> &nodesToPrune) {
        long long updateUs = 0;
        for (NodeId rootSubtree : nodesToPrune) {
            if (!isValidPruneRoot(maxTree_, rootSubtree)) {
                continue;
            }

            const auto start = std::chrono::steady_clock::now();
            const NodeId parentId = maxTree_.getNodeParent(rootSubtree);
            const auto b = static_cast<AltitudeType>(maxTree_.getAltitude(parentId));
            const auto pixelsC = collectSubtreePixels(maxTree_, rootSubtree);
            adjust_.updateTreeWithExplicitRegion(&minTree_, std::span<const PixelId>(pixelsC), b);
            maxTree_.pruneNode(rootSubtree);
            updateUs += common::elapsedMicros(start, std::chrono::steady_clock::now());
        }
        return updateUs;
    }

    long long pruneMinTreeAndUpdateMaxTreeExplicit(std::vector<NodeId> &nodesToPrune) {
        long long updateUs = 0;
        for (NodeId rootSubtree : nodesToPrune) {
            if (!isValidPruneRoot(minTree_, rootSubtree)) {
                continue;
            }

            const auto start = std::chrono::steady_clock::now();
            const NodeId parentId = minTree_.getNodeParent(rootSubtree);
            const auto b = static_cast<AltitudeType>(minTree_.getAltitude(parentId));
            const auto pixelsC = collectSubtreePixels(minTree_, rootSubtree);
            adjust_.updateTreeWithExplicitRegion(&maxTree_, std::span<const PixelId>(pixelsC), b);
            minTree_.pruneNode(rootSubtree);
            updateUs += common::elapsedMicros(start, std::chrono::steady_clock::now());
        }
        return updateUs;
    }

    long long applyProblem1FastGateEvent(DynamicComponentTree &primalTree,
                                         DynamicComponentTree &dualTree,
                                         NodeId rootSubtree,
                                         bool updateMaxTree,
                                         BackendRunResult &result) {
        if (!isValidPruneRoot(primalTree, rootSubtree)) {
            return 0;
        }

        ++result.events;
        const auto start = std::chrono::steady_clock::now();
        const NodeId parentId = primalTree.getNodeParent(rootSubtree);
        const auto b = static_cast<AltitudeType>(primalTree.getAltitude(parentId));
        const auto pixelsC = collectSubtreePixels(primalTree, rootSubtree);
        const auto gate =
            problem1_debug::checkCasfSubtreeProblem1FastGate<AltitudeType>(
                dualTree,
                std::span<const PixelId>(pixelsC),
                b);

        if (!gate.readyUnderCasfAssumptions) {
            ++result.fastGateRejected;
            ++result.fallbackToUpdateTree;
            adjust_.updateTree(&dualTree, rootSubtree);
        } else {
            ++result.fastGateReady;
            problem1_debug::ExplicitRegionProblem1UpdateOptions options;
            options.assumeConnected = true;
            options.allowLocalMaterialization = true;
            options.allowTargetNodeCreation = true;
            options.enforceLocalCorrectnessValidation = false;
            options.validateAgainstRebuild = false;
            options.validateStructureAgainstRebuild = false;
            options.runtimePostConditionValidation = false;

            const auto updateResult =
                problem1_debug::updateExplicitRegionProblem1WithDiagnostics<AltitudeType>(
                    minTree_,
                    maxTree_,
                    adj_,
                    updateMaxTree,
                    std::span<const PixelId>(pixelsC),
                    b,
                    options);
            if (!updateResult.appliedLocalUpdate()) {
                ++result.localRejected;
                throw std::runtime_error(
                    "problem1 fast-gate explicit update did not apply locally: " +
                    updateResult.message);
            }
            ++result.localApplied;
        }

        primalTree.pruneNode(rootSubtree);
        return common::elapsedMicros(start, std::chrono::steady_clock::now());
    }

    long long applyProblem2ExplicitSweepBatch(DynamicComponentTree &primalTree,
                                              DynamicComponentTree &dualTree,
                                              std::vector<NodeId> &nodesToPrune,
                                              BackendRunResult &result) {
        const auto start = std::chrono::steady_clock::now();
        if (!profileProblem2ExplicitSweep_) {
            const int appliedCount =
                adjust_.pruneAndUpdate(&primalTree, &dualTree, nodesToPrune);
            const long long batchWallUs =
                common::elapsedMicros(start, std::chrono::steady_clock::now());

            result.problem2BatchWallUs += batchWallUs;
            result.events += appliedCount;
            result.fastGateReady += appliedCount;
            result.localApplied += appliedCount;
            return batchWallUs;
        }

        const auto updateResult =
            problem2_debug::pruneAndUpdateWithDiagnostics(
                adjust_,
                &primalTree,
                &dualTree,
                nodesToPrune,
                true);
        const long long batchWallUs =
            common::elapsedMicros(start, std::chrono::steady_clock::now());

        const long long eventAccountedUs =
            updateResult.eventSetupUs +
            updateResult.collectGateUs +
            updateResult.runUpdateTreeWallUs +
            updateResult.pruneUs;
        const long long kernelAccountedUs =
            updateResult.bucketBuildUs +
            updateResult.targetNodeCreationUs +
            updateResult.metadataUs +
            updateResult.levelPrepUs +
            updateResult.sweepUs +
            updateResult.finalizeUs +
            updateResult.finalReconnectUs +
            updateResult.releaseEmptyUs +
            updateResult.sameLevelCollapseUs;
        result.problem2BatchWallUs += batchWallUs;
        result.problem2EventTotalUs += updateResult.eventTotalUs;
        result.problem2EventSetupUs += updateResult.eventSetupUs;
        result.problem2RunUpdateTreeWallUs += updateResult.runUpdateTreeWallUs;
        result.problem2PruneUs += updateResult.pruneUs;
        if (updateResult.eventTotalUs > 0) {
            result.problem2EventUnaccountedUs +=
                std::max<long long>(0, updateResult.eventTotalUs - eventAccountedUs);
            result.problem2BatchOverheadUs +=
                std::max<long long>(0, batchWallUs - updateResult.eventTotalUs);
        }
        if (updateResult.runUpdateTreeWallUs > 0) {
            result.problem2KernelUnaccountedUs +=
                std::max<long long>(0, updateResult.runUpdateTreeWallUs - kernelAccountedUs);
        }
        if (!updateResult.inputValid) {
            result.collectGateUs += updateResult.collectGateUs;
            ++result.localRejected;
            throw std::runtime_error(
                "problem2 explicit-sweep prune/update received an invalid event: " +
                updateResult.message);
        }

        result.events += updateResult.eventCount;
        result.skippedNodes += updateResult.skippedNodeCount;
        result.collectGateUs += updateResult.collectGateUs;
        result.bucketBuildUs += updateResult.bucketBuildUs;
        result.targetNodeCreationUs += updateResult.targetNodeCreationUs;
        result.metadataUs += updateResult.metadataUs;
        result.levelPrepUs += updateResult.levelPrepUs;
        result.sweepUs += updateResult.sweepUs;
        result.finalizeUs += updateResult.finalizeUs;
        result.finalReconnectUs += updateResult.finalReconnectUs;
        result.releaseEmptyUs += updateResult.releaseEmptyUs;
        result.sameLevelCollapseUs += updateResult.sameLevelCollapseUs;
        result.regionPixels += updateResult.regionSize;
        result.mergeNodes += updateResult.mergeNodeCount;
        result.mergeLevels += updateResult.mergeLevelCount;
        result.frontierNodes += updateResult.frontierNodeCount;
        result.sameLevelEdgesCollapsed += updateResult.sameLevelEdgesCollapsed;
        result.maxRegionPixels = std::max(result.maxRegionPixels, updateResult.regionSize);
        result.maxBucketSize = std::max(result.maxBucketSize, updateResult.maxBucketSize);
        if (updateResult.targetNodeCreated) {
            ++result.targetNodeCreatedEvents;
        }
        if (updateResult.eventCount == 0) {
            return batchWallUs;
        }

        if (updateResult.fastGateReady && updateResult.localUpdateApplied && updateResult.primalPruneApplied) {
            result.fastGateReady += updateResult.eventCount;
            result.localApplied += updateResult.prunedNodeCount;
        } else if (!updateResult.fastGateReady) {
            result.fastGateRejected += updateResult.eventCount;
            ++result.localRejected;
            throw std::runtime_error(
                "problem2 explicit-sweep prune/update rejected the event without fallback: " +
                updateResult.message);
        } else {
            ++result.localRejected;
            throw std::runtime_error(
                "problem2 explicit-sweep prune/update did not complete locally: " +
                updateResult.message);
        }

        return batchWallUs;
    }

    long long pruneMaxTreeAndUpdateMinTreeProblem1FastGate(
        std::vector<NodeId> &nodesToPrune,
        BackendRunResult &result) {
        long long updateUs = 0;
        for (NodeId rootSubtree : nodesToPrune) {
            updateUs += applyProblem1FastGateEvent(
                maxTree_,
                minTree_,
                rootSubtree,
                false,
                result);
        }
        return updateUs;
    }

    long long pruneMinTreeAndUpdateMaxTreeProblem1FastGate(
        std::vector<NodeId> &nodesToPrune,
        BackendRunResult &result) {
        long long updateUs = 0;
        for (NodeId rootSubtree : nodesToPrune) {
            updateUs += applyProblem1FastGateEvent(
                minTree_,
                maxTree_,
                rootSubtree,
                true,
                result);
        }
        return updateUs;
    }

    long long pruneMaxTreeAndUpdateMinTreeProblem2ExplicitSweep(
        std::vector<NodeId> &nodesToPrune,
        BackendRunResult &result) {
        return applyProblem2ExplicitSweepBatch(
            maxTree_,
            minTree_,
            nodesToPrune,
            result);
    }

    long long pruneMinTreeAndUpdateMaxTreeProblem2ExplicitSweep(
        std::vector<NodeId> &nodesToPrune,
        BackendRunResult &result) {
        return applyProblem2ExplicitSweepBatch(
            minTree_,
            maxTree_,
            nodesToPrune,
            result);
    }

public:
    UpdateBackendRunner(ImageUInt8Ptr image,
                        double radioAdj,
                        Attribute attribute,
                        BackendKind backend,
                        bool profileProblem2ExplicitSweep)
        : inputImage_(image->clone()),
          adj_(std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radioAdj)),
          maxTree_(inputImage_, true, adj_),
          minTree_(inputImage_, false, adj_),
          adjust_(&minTree_, &maxTree_, *adj_),
          maxAttributeComputer_(common::makeAttributeComputer(&maxTree_, attribute)),
          minAttributeComputer_(common::makeAttributeComputer(&minTree_, attribute)),
          attribute_(attribute),
          maxArea_(maxTree_.getNumInternalNodeSlots(), 0.0f),
          minArea_(minTree_.getNumInternalNodeSlots(), 0.0f),
          backend_(backend),
          profileProblem2ExplicitSweep_(profileProblem2ExplicitSweep) {
        maxArea_ = common::computeAttributeVector(&maxTree_, attribute_);
        minArea_ = common::computeAttributeVector(&minTree_, attribute_);
        rebindAttributeBuffers();
        adjust_.setRuntimePostConditionValidationEnabled(false);
    }

    BackendRunResult run(const std::vector<int> &thresholds) {
        if ((backend_ == BackendKind::Problem1FastGateExplicit ||
             backend_ == BackendKind::Problem2ExplicitSweep) &&
            attribute_ != AREA) {
            throw std::runtime_error(
                "Experimental explicit-sweep comparison backends currently require --attribute area");
        }

        BackendRunResult result;
        const auto totalStart = std::chrono::steady_clock::now();

        for (int threshold : thresholds) {
            auto nodesToPrune = common::getDynamicNodesThreshold(maxTree_, maxArea_, threshold);
            if (backend_ == BackendKind::ExplicitRegionDirect) {
                result.events += countValidPruneRoots(maxTree_, nodesToPrune);
                result.phase1UpdateUs += pruneMaxTreeAndUpdateMinTreeExplicit(nodesToPrune);
            } else if (backend_ == BackendKind::Problem1FastGateExplicit) {
                result.phase1UpdateUs +=
                    pruneMaxTreeAndUpdateMinTreeProblem1FastGate(nodesToPrune, result);
                result.attributeRefreshUs += refreshAllAttributes();
            } else if (backend_ == BackendKind::Problem2ExplicitSweep) {
                result.phase1UpdateUs +=
                    pruneMaxTreeAndUpdateMinTreeProblem2ExplicitSweep(nodesToPrune, result);
            } else {
                result.events += countValidPruneRoots(maxTree_, nodesToPrune);
                const auto start = std::chrono::steady_clock::now();
                adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
                result.phase1UpdateUs += common::elapsedMicros(start, std::chrono::steady_clock::now());
            }

            nodesToPrune = common::getDynamicNodesThreshold(minTree_, minArea_, threshold);
            if (backend_ == BackendKind::ExplicitRegionDirect) {
                result.events += countValidPruneRoots(minTree_, nodesToPrune);
                result.phase2UpdateUs += pruneMinTreeAndUpdateMaxTreeExplicit(nodesToPrune);
            } else if (backend_ == BackendKind::Problem1FastGateExplicit) {
                result.phase2UpdateUs +=
                    pruneMinTreeAndUpdateMaxTreeProblem1FastGate(nodesToPrune, result);
                result.attributeRefreshUs += refreshAllAttributes();
            } else if (backend_ == BackendKind::Problem2ExplicitSweep) {
                result.phase2UpdateUs +=
                    pruneMinTreeAndUpdateMaxTreeProblem2ExplicitSweep(nodesToPrune, result);
            } else {
                result.events += countValidPruneRoots(minTree_, nodesToPrune);
                const auto start = std::chrono::steady_clock::now();
                adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
                result.phase2UpdateUs += common::elapsedMicros(start, std::chrono::steady_clock::now());
            }
        }

        result.output = minTree_.reconstructionImage();
        result.totalUs = common::elapsedMicros(totalStart, std::chrono::steady_clock::now());
        result.updateUs =
            result.phase1UpdateUs + result.phase2UpdateUs + result.attributeRefreshUs;
        return result;
    }
};

BackendRunResult runBackend(ImageUInt8Ptr image,
                            double radioAdj,
                            Attribute attribute,
                            const std::vector<int> &thresholds,
                            BackendKind backend,
                            bool profileProblem2ExplicitSweep) {
    UpdateBackendRunner runner(
        image,
        radioAdj,
        attribute,
        backend,
        profileProblem2ExplicitSweep);
    return runner.run(thresholds);
}

std::string basename(const std::string &path) {
    return std::filesystem::path(path).filename().string();
}

void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0
              << " [--radio-adj <radius>]"
              << " [--attribute area|bbox_width|bbox_height|bbox_diagonal]"
              << " [--thresholds <csv>]"
              << " [--include-problem1-fast-gate]"
              << " [--include-problem2-explicit-sweep]"
              << " [--only-problem2-explicit-sweep]"
              << " [--profile-problem2-explicit-sweep]"
              << " <image1> [image2 ...]\n";
}

} // namespace

int main(int argc, char **argv) {
    try {
        double radioAdj = 1.5;
        Attribute attribute = AREA;
        bool includeProblem1FastGate = false;
        bool includeProblem2ExplicitSweep = false;
        bool onlyProblem2ExplicitSweep = false;
        bool profileProblem2ExplicitSweep = false;
        std::string thresholdsCsv =
            "10,20,30,40,49,59,69,79,89,99,109,119,128,138,148,158,168,178,188,198,"
            "207,217,227,237,247,257,267,277,286,296,306,316,326,336,346,356,365,375,"
            "385,395,405,415,425,435,444,454,464,474,484,494";

        int argi = 1;
        while (argi < argc && std::string(argv[argi]).rfind("--", 0) == 0) {
            const std::string arg = argv[argi];
            if (arg == "--radio-adj") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                radioAdj = std::stod(argv[argi + 1]);
                argi += 2;
                continue;
            }
            if (arg == "--attribute") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                attribute = common::parseAttributeName(argv[argi + 1]);
                argi += 2;
                continue;
            }
            if (arg == "--thresholds") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                thresholdsCsv = argv[argi + 1];
                argi += 2;
                continue;
            }
            if (arg == "--include-problem1-fast-gate") {
                includeProblem1FastGate = true;
                ++argi;
                continue;
            }
            if (arg == "--include-problem2-explicit-sweep") {
                includeProblem2ExplicitSweep = true;
                ++argi;
                continue;
            }
            if (arg == "--only-problem2-explicit-sweep") {
                includeProblem2ExplicitSweep = true;
                onlyProblem2ExplicitSweep = true;
                ++argi;
                continue;
            }
            if (arg == "--profile-problem2-explicit-sweep") {
                profileProblem2ExplicitSweep = true;
                ++argi;
                continue;
            }
            throw std::runtime_error("Unknown option: " + arg);
        }

        if (argi >= argc) {
            printUsage(argv[0]);
            return 1;
        }

        const auto thresholds = parseThresholdCsv(thresholdsCsv);
        std::cout
            << "image,rows,cols,thresholds,backend,total_ms,update_ms,"
            << "phase1_update_ms,phase2_update_ms,attribute_refresh_ms,"
            << "problem2_batch_wall_ms,problem2_event_total_ms,"
            << "problem2_event_setup_ms,"
            << "collect_gate_ms,problem2_run_update_tree_wall_ms,"
            << "bucket_build_ms,target_creation_ms,metadata_ms,level_prep_ms,"
            << "sweep_ms,finalize_ms,final_reconnect_ms,"
            << "release_empty_ms,same_level_collapse_ms,problem2_prune_ms,"
            << "problem2_event_unaccounted_ms,problem2_kernel_unaccounted_ms,"
            << "problem2_batch_overhead_ms,"
            << "events,skipped_nodes,fast_gate_ready,fast_gate_rejected,local_applied,"
            << "local_rejected,fallback_to_update_tree,attribute_buffer_rebinds,"
            << "region_pixels,max_region_pixels,merge_nodes,merge_levels,"
            << "frontier_nodes,max_bucket_size,same_level_edges_collapsed,"
            << "target_node_created_events,"
            << "matches_update_tree\n";

        for (; argi < argc; ++argi) {
            const std::string imagePath = argv[argi];
            const auto image = common::loadGrayImage(imagePath);

            const auto updateTree =
                runBackend(
                    image,
                    radioAdj,
                    attribute,
                    thresholds,
                    BackendKind::UpdateTree,
                    profileProblem2ExplicitSweep);

            auto printRow = [&](const char *backend, const BackendRunResult &result) {
                const bool matches = updateTree.output->isEqual(result.output);
                std::cout << basename(imagePath) << ","
                          << image->getNumRows() << ","
                          << image->getNumCols() << ","
                          << thresholds.size() << ","
                          << backend << ","
                          << common::microsToMillis(result.totalUs) << ","
                          << common::microsToMillis(result.updateUs) << ","
                          << common::microsToMillis(result.phase1UpdateUs) << ","
                          << common::microsToMillis(result.phase2UpdateUs) << ","
                          << common::microsToMillis(result.attributeRefreshUs) << ","
                          << common::microsToMillis(result.problem2BatchWallUs) << ","
                          << common::microsToMillis(result.problem2EventTotalUs) << ","
                          << common::microsToMillis(result.problem2EventSetupUs) << ","
                          << common::microsToMillis(result.collectGateUs) << ","
                          << common::microsToMillis(result.problem2RunUpdateTreeWallUs) << ","
                          << common::microsToMillis(result.bucketBuildUs) << ","
                          << common::microsToMillis(result.targetNodeCreationUs) << ","
                          << common::microsToMillis(result.metadataUs) << ","
                          << common::microsToMillis(result.levelPrepUs) << ","
                          << common::microsToMillis(result.sweepUs) << ","
                          << common::microsToMillis(result.finalizeUs) << ","
                          << common::microsToMillis(result.finalReconnectUs) << ","
                          << common::microsToMillis(result.releaseEmptyUs) << ","
                          << common::microsToMillis(result.sameLevelCollapseUs) << ","
                          << common::microsToMillis(result.problem2PruneUs) << ","
                          << common::microsToMillis(result.problem2EventUnaccountedUs) << ","
                          << common::microsToMillis(result.problem2KernelUnaccountedUs) << ","
                          << common::microsToMillis(result.problem2BatchOverheadUs) << ","
                          << result.events << ","
                          << result.skippedNodes << ","
                          << result.fastGateReady << ","
                          << result.fastGateRejected << ","
                          << result.localApplied << ","
                          << result.localRejected << ","
                          << result.fallbackToUpdateTree << ","
                          << result.attributeBufferRebinds << ","
                          << result.regionPixels << ","
                          << result.maxRegionPixels << ","
                          << result.mergeNodes << ","
                          << result.mergeLevels << ","
                          << result.frontierNodes << ","
                          << result.maxBucketSize << ","
                          << result.sameLevelEdgesCollapsed << ","
                          << result.targetNodeCreatedEvents << ","
                          << (matches ? "true" : "false")
                          << "\n";
            };

            printRow("updateTree", updateTree);
            if (!onlyProblem2ExplicitSweep) {
                const auto explicitRegion =
                    runBackend(
                        image,
                        radioAdj,
                        attribute,
                        thresholds,
                        BackendKind::ExplicitRegionDirect,
                        profileProblem2ExplicitSweep);
                printRow("updateTreeWithExplicitRegion", explicitRegion);
            }
            if (includeProblem1FastGate) {
                const auto problem1FastGate =
                    runBackend(image,
                               radioAdj,
                               attribute,
                               thresholds,
                               BackendKind::Problem1FastGateExplicit,
                               profileProblem2ExplicitSweep);
                printRow("problem1FastGateExplicit", problem1FastGate);
            }
            if (includeProblem2ExplicitSweep) {
                const auto problem2ExplicitSweep =
                    runBackend(image,
                               radioAdj,
                               attribute,
                               thresholds,
                               BackendKind::Problem2ExplicitSweep,
                               profileProblem2ExplicitSweep);
                printRow("problem2ExplicitSweep", problem2ExplicitSweep);
            }
        }
    } catch (const std::exception &e) {
        std::cerr << "compare_update_backends: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
