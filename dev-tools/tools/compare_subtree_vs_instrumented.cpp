#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "../../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../../morphoTreeAdjust/include/Common.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTreeAdjustment.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTreeAdjustmentInstrumented.hpp"

#include "../external/stb/stb_image.h"

struct TimerHooksContext {
    Stopwatch *sw = nullptr;
};

static double elapsedMs(const Stopwatch &sw) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(sw.elapsed()).count();
}

static void pauseTimerHooks(void *ctx) {
    auto *hooks = static_cast<TimerHooksContext *>(ctx);
    if (hooks != nullptr && hooks->sw != nullptr && hooks->sw->running()) {
        hooks->sw->pause();
    }
}

static void resumeTimerHooks(void *ctx) {
    auto *hooks = static_cast<TimerHooksContext *>(ctx);
    if (hooks != nullptr && hooks->sw != nullptr && !hooks->sw->running()) {
        hooks->sw->resume();
    }
}

static ImageUInt8Ptr loadGrayImage(const std::string &path) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char *data = stbi_load(path.c_str(), &width, &height, &channels, 1);
    if (data == nullptr) {
        throw std::runtime_error("failed to load image: " + path);
    }

    ImageUInt8Ptr image = ImageUInt8::create(height, width);
    std::copy(data,
              data + static_cast<std::size_t>(width) * static_cast<std::size_t>(height),
              image->rawData());
    stbi_image_free(data);
    return image;
}

static std::vector<float> computeAreaAttribute(DynamicComponentTree *tree) {
    DynamicAreaComputer computer(tree);
    return computer.compute();
}

static std::vector<NodeId> getNodesThreshold(DynamicComponentTree *tree, const std::vector<float> &attribute, int threshold) {
    std::vector<NodeId> nodesToPrune;
    nodesToPrune.reserve(static_cast<std::size_t>(std::max(0, tree->getGlobalIdSpaceSize() / 8)));
    for (NodeId nodeId = 0; nodeId < tree->getGlobalIdSpaceSize(); ++nodeId) {
        if (!tree->isNode(nodeId) || !tree->isAlive(nodeId) || tree->isRoot(nodeId)) {
            continue;
        }
        if (static_cast<std::size_t>(nodeId) >= attribute.size()) {
            continue;
        }
        if (attribute[static_cast<std::size_t>(nodeId)] < static_cast<float>(threshold)) {
            nodesToPrune.push_back(nodeId);
        }
    }
    return nodesToPrune;
}

static std::vector<int> parseThresholds(const std::string &text) {
    std::vector<int> out;
    std::stringstream ss(text);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (!token.empty()) {
            out.push_back(std::stoi(token));
        }
    }
    return out;
}

static bool sameImage(const ImageUInt8Ptr &a, const ImageUInt8Ptr &b) {
    if (a == nullptr || b == nullptr) {
        return false;
    }
    if (a->getNumRows() != b->getNumRows() || a->getNumCols() != b->getNumCols()) {
        return false;
    }
    return a->isEqual(b);
}

struct RunResult {
    double phase1UpdateMs = 0.0;
    double phase2UpdateMs = 0.0;
    ImageUInt8Ptr output;
};

static RunResult runBase(ImageUInt8Ptr image, double radioAdj, const std::vector<int> &thresholds) {
    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radioAdj);
    DynamicComponentTree maxTree(image, true, adj);
    DynamicComponentTree minTree(image, false, adj);
    DynamicComponentTreeAdjustment<AltitudeType> adjust(&minTree, &maxTree, *adj);
    adjust.setRuntimePostConditionValidationEnabled(true);

    Stopwatch phase1UpdateSw;
    Stopwatch phase2UpdateSw;

    for (int threshold : thresholds) {
        const auto attributeMax = computeAreaAttribute(&maxTree);
        const auto nodesToPruneMax = getNodesThreshold(&maxTree, attributeMax, threshold);
        for (NodeId subtreeRoot : nodesToPruneMax) {
            if (subtreeRoot == InvalidNode || subtreeRoot == maxTree.getRoot() || !maxTree.isAlive(subtreeRoot)) {
                continue;
            }
            phase1UpdateSw.resume();
            adjust.updateTree(&minTree, subtreeRoot);
            phase1UpdateSw.pause();
            maxTree.pruneNode(subtreeRoot);
        }

        const auto attributeMin = computeAreaAttribute(&minTree);
        const auto nodesToPruneMin = getNodesThreshold(&minTree, attributeMin, threshold);
        for (NodeId subtreeRoot : nodesToPruneMin) {
            if (subtreeRoot == InvalidNode || subtreeRoot == minTree.getRoot() || !minTree.isAlive(subtreeRoot)) {
                continue;
            }
            phase2UpdateSw.resume();
            adjust.updateTree(&maxTree, subtreeRoot);
            phase2UpdateSw.pause();
            minTree.pruneNode(subtreeRoot);
        }
    }

    return {elapsedMs(phase1UpdateSw), elapsedMs(phase2UpdateSw), minTree.reconstructionImage()};
}

static RunResult runInstrumented(ImageUInt8Ptr image, double radioAdj, const std::vector<int> &thresholds, bool enableMetrics) {
    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radioAdj);
    DynamicComponentTree maxTree(image, true, adj);
    DynamicComponentTree minTree(image, false, adj);
    DynamicComponentTreeAdjustmentInstrumented<AltitudeType> adjust(&minTree, &maxTree, *adj);
    adjust.setRuntimePostConditionValidationEnabled(true);

    Stopwatch phase1UpdateSw;
    Stopwatch phase2UpdateSw;

    for (int threshold : thresholds) {
        const auto attributeMax = computeAreaAttribute(&maxTree);
        std::vector<float> targetAreaPhase1;
        if (enableMetrics) {
            DynamicAreaComputer targetAreaComputer(&minTree);
            targetAreaPhase1 = targetAreaComputer.compute();
        }
        const auto nodesToPruneMax = getNodesThreshold(&maxTree, attributeMax, threshold);
        TimerHooksContext phase1Hooks{&phase1UpdateSw};
        for (NodeId subtreeRoot : nodesToPruneMax) {
            if (subtreeRoot == InvalidNode || subtreeRoot == maxTree.getRoot() || !maxTree.isAlive(subtreeRoot)) {
                continue;
            }
            DynamicSubtreeMetrics metrics;
            if (enableMetrics) {
                adjust.setMetricsAreaBuffer(&minTree, std::span<const float>(targetAreaPhase1));
                adjust.setMetrics(&metrics, &phase1Hooks, pauseTimerHooks, resumeTimerHooks);
            } else {
                adjust.setMetricsAreaBuffer(nullptr, {});
                adjust.setMetrics(nullptr);
            }
            phase1UpdateSw.resume();
            adjust.updateTree(&minTree, subtreeRoot);
            phase1UpdateSw.pause();
            maxTree.pruneNode(subtreeRoot);
        }
        adjust.setMetrics(nullptr);

        const auto attributeMin = computeAreaAttribute(&minTree);
        std::vector<float> targetAreaPhase2;
        if (enableMetrics) {
            DynamicAreaComputer targetAreaComputer(&maxTree);
            targetAreaPhase2 = targetAreaComputer.compute();
        }
        const auto nodesToPruneMin = getNodesThreshold(&minTree, attributeMin, threshold);
        TimerHooksContext phase2Hooks{&phase2UpdateSw};
        for (NodeId subtreeRoot : nodesToPruneMin) {
            if (subtreeRoot == InvalidNode || subtreeRoot == minTree.getRoot() || !minTree.isAlive(subtreeRoot)) {
                continue;
            }
            DynamicSubtreeMetrics metrics;
            if (enableMetrics) {
                adjust.setMetricsAreaBuffer(&maxTree, std::span<const float>(targetAreaPhase2));
                adjust.setMetrics(&metrics, &phase2Hooks, pauseTimerHooks, resumeTimerHooks);
            } else {
                adjust.setMetricsAreaBuffer(nullptr, {});
                adjust.setMetrics(nullptr);
            }
            phase2UpdateSw.resume();
            adjust.updateTree(&maxTree, subtreeRoot);
            phase2UpdateSw.pause();
            minTree.pruneNode(subtreeRoot);
        }
        adjust.setMetrics(nullptr);
    }

    return {elapsedMs(phase1UpdateSw), elapsedMs(phase2UpdateSw), minTree.reconstructionImage()};
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <image> [radio-adj] [thresholds-csv]\n";
        return 1;
    }

    const std::string imagePath = argv[1];
    const double radioAdj = (argc >= 3) ? std::stod(argv[2]) : 1.5;
    const std::string thresholdsArg = (argc >= 4) ? argv[3] : "50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500";
    const std::vector<int> thresholds = parseThresholds(thresholdsArg);

    const ImageUInt8Ptr image = loadGrayImage(imagePath);

    const RunResult base = runBase(image, radioAdj, thresholds);
    const RunResult instrumentedNoMetrics = runInstrumented(image, radioAdj, thresholds, false);
    const RunResult instrumentedWithMetrics = runInstrumented(image, radioAdj, thresholds, true);

    std::cout << "Input: " << imagePath << "\n";
    std::cout << "Image size: " << image->getNumRows() << "x" << image->getNumCols() << "\n";
    std::cout << "Threshold count: " << thresholds.size() << "\n";
    std::cout << "Adjacency radius: " << radioAdj << "\n\n";

    std::cout << "Update-only timing (ms)\n";
    std::cout << "  base phase1=" << base.phase1UpdateMs << " phase2=" << base.phase2UpdateMs << " total=" << (base.phase1UpdateMs + base.phase2UpdateMs) << "\n";
    std::cout << "  instrumented(no-metrics) phase1=" << instrumentedNoMetrics.phase1UpdateMs << " phase2=" << instrumentedNoMetrics.phase2UpdateMs << " total=" << (instrumentedNoMetrics.phase1UpdateMs + instrumentedNoMetrics.phase2UpdateMs) << "\n";
    std::cout << "  instrumented(with-metrics-paused) phase1=" << instrumentedWithMetrics.phase1UpdateMs << " phase2=" << instrumentedWithMetrics.phase2UpdateMs << " total=" << (instrumentedWithMetrics.phase1UpdateMs + instrumentedWithMetrics.phase2UpdateMs) << "\n\n";

    std::cout << "Final image comparison\n";
    std::cout << "  base == instrumented(no-metrics): " << (sameImage(base.output, instrumentedNoMetrics.output) ? "true" : "false") << "\n";
    std::cout << "  base == instrumented(with-metrics-paused): " << (sameImage(base.output, instrumentedWithMetrics.output) ? "true" : "false") << "\n";

    return 0;
}
