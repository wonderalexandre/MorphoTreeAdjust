
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <memory>
#include <sstream>
#include <span>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>

#include "../../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../../morphoTreeAdjust/include/Common.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilterInstrumented.hpp"

namespace {

enum class BenchmarkImageModel {
    structured,
    natural_like,
    piecewise_scene
};

const char *benchmark_image_model_name(BenchmarkImageModel model) {
    switch (model) {
        case BenchmarkImageModel::structured:
            return "structured";
        case BenchmarkImageModel::natural_like:
            return "natural_like";
        case BenchmarkImageModel::piecewise_scene:
            return "piecewise_scene";
    }
    return "unknown";
}

std::string benchmark_case_label(BenchmarkImageModel model, int numThresholds) {
    return std::string(benchmark_image_model_name(model)) + ", thresholds=" + std::to_string(numThresholds);
}

uint8_t clamp_to_u8(double value) {
    return static_cast<uint8_t>(std::clamp<int>(static_cast<int>(std::lround(value)), 0, 255));
}

// Structured synthetic pattern with gradients, periodic bands, and local texture.
ImageUInt8Ptr make_structured_benchmark_image(int numRows, int numCols) {
    auto image = ImageUInt8::create(numRows, numCols);
    uint8_t *data = image->rawData();
    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            const int p = r * numCols + c;
            const int value = (37 * r + 19 * c + 11 * ((r / 4) % 7) + 23 * ((c / 8) % 5) + ((r * c) % 29)) % 256;
            data[p] = static_cast<uint8_t>(value);
        }
    }
    return image;
}

// Natural-like synthetic pattern built from smooth multi-scale value noise and illumination.
ImageUInt8Ptr make_natural_like_benchmark_image(int numRows, int numCols) {
    auto image = ImageUInt8::create(numRows, numCols);
    uint8_t *data = image->rawData();

    const auto hash01 = [](uint32_t y, uint32_t x, uint32_t seed) -> double {
        uint32_t h = seed;
        h ^= y + 0x9e3779b9u + (h << 6) + (h >> 2);
        h ^= x + 0x9e3779b9u + (h << 6) + (h >> 2);
        h *= 0x85ebca6bu;
        h ^= h >> 13;
        h *= 0xc2b2ae35u;
        h ^= h >> 16;
        return static_cast<double>(h & 0x00ffffffu) / static_cast<double>(0x01000000u);
    };

    const auto smoothstep = [](double t) -> double {
        return t * t * (3.0 - 2.0 * t);
    };

    const auto value_noise = [&](double y, double x, double scale, uint32_t seed) -> double {
        const double ys = y / scale;
        const double xs = x / scale;
        const auto y0 = static_cast<uint32_t>(std::floor(ys));
        const auto x0 = static_cast<uint32_t>(std::floor(xs));
        const auto y1 = y0 + 1u;
        const auto x1 = x0 + 1u;
        const double fy = smoothstep(ys - std::floor(ys));
        const double fx = smoothstep(xs - std::floor(xs));

        const double v00 = hash01(y0, x0, seed);
        const double v01 = hash01(y0, x1, seed);
        const double v10 = hash01(y1, x0, seed);
        const double v11 = hash01(y1, x1, seed);

        const double vx0 = v00 + fx * (v01 - v00);
        const double vx1 = v10 + fx * (v11 - v10);
        return vx0 + fy * (vx1 - vx0);
    };

    const auto fbm = [&](double y, double x) -> double {
        double sum = 0.0;
        double amplitude = 1.0;
        double norm = 0.0;
        std::array<double, 4> scales = {48.0, 24.0, 12.0, 6.0};
        std::array<uint32_t, 4> seeds = {0x12345u, 0x23456u, 0x34567u, 0x45678u};
        for (std::size_t i = 0; i < scales.size(); ++i) {
            sum += amplitude * value_noise(y, x, scales[i], seeds[i]);
            norm += amplitude;
            amplitude *= 0.5;
        }
        return sum / norm;
    };

    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            const int p = r * numCols + c;
            const double y = static_cast<double>(r);
            const double x = static_cast<double>(c);

            const double illumination =
                    58.0
                    + 34.0 * std::sin(0.010 * x + 0.014 * y)
                    + 27.0 * std::cos(0.008 * x - 0.012 * y);

            const double textureLarge = 110.0 * fbm(0.85 * y, 0.85 * x);
            const double textureSmall = 26.0 * value_noise(y, x, 3.0, 0x56789u);
            const double oriented =
                    11.0 * std::sin(0.041 * x + 0.063 * y)
                    + 8.0 * std::cos(0.057 * x - 0.034 * y);

            const double edgeLike =
                    ((r / 19 + c / 23) % 2 == 0) ? 9.0 : -9.0;

            data[p] = clamp_to_u8(illumination + textureLarge + textureSmall + oriented + edgeLike - 40.0);
        }
    }

    return image;
}

// Scene-like synthetic pattern with smooth objects, background illumination, and mild texture.
ImageUInt8Ptr make_piecewise_scene_benchmark_image(int numRows, int numCols) {
    auto image = ImageUInt8::create(numRows, numCols);
    uint8_t *data = image->rawData();

    const auto smoothstep = [](double edge0, double edge1, double x) -> double {
        const double t = std::clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
        return t * t * (3.0 - 2.0 * t);
    };

    const auto soft_disk = [&](double x, double y, double cx, double cy, double radius, double feather) -> double {
        const double d = std::sqrt((x - cx) * (x - cx) + (y - cy) * (y - cy));
        return 1.0 - smoothstep(radius - feather, radius + feather, d);
    };

    const auto soft_ellipse = [&](double x, double y, double cx, double cy, double rx, double ry, double feather) -> double {
        const double dx = (x - cx) / rx;
        const double dy = (y - cy) / ry;
        const double d = std::sqrt(dx * dx + dy * dy);
        const double normalizedFeather = feather / std::max(rx, ry);
        return 1.0 - smoothstep(1.0 - normalizedFeather, 1.0 + normalizedFeather, d);
    };

    const auto soft_rect = [&](double x, double y, double x0, double y0, double x1, double y1, double feather) -> double {
        const double left = smoothstep(x0 - feather, x0 + feather, x);
        const double right = 1.0 - smoothstep(x1 - feather, x1 + feather, x);
        const double top = smoothstep(y0 - feather, y0 + feather, y);
        const double bottom = 1.0 - smoothstep(y1 - feather, y1 + feather, y);
        return left * right * top * bottom;
    };

    const double cx1 = 0.26 * static_cast<double>(numCols);
    const double cy1 = 0.34 * static_cast<double>(numRows);
    const double rx1 = 0.18 * static_cast<double>(numCols);
    const double ry1 = 0.24 * static_cast<double>(numRows);

    const double cx2 = 0.71 * static_cast<double>(numCols);
    const double cy2 = 0.60 * static_cast<double>(numRows);
    const double radius2 = 0.17 * static_cast<double>(std::min(numRows, numCols));

    const double rectX0 = 0.58 * static_cast<double>(numCols);
    const double rectY0 = 0.14 * static_cast<double>(numRows);
    const double rectX1 = 0.90 * static_cast<double>(numCols);
    const double rectY1 = 0.39 * static_cast<double>(numRows);

    for (int r = 0; r < numRows; ++r) {
        for (int c = 0; c < numCols; ++c) {
            const int p = r * numCols + c;
            const double y = static_cast<double>(r);
            const double x = static_cast<double>(c);

            double value = 34.0 + 0.26 * y + 0.14 * x;
            value += 18.0 * std::sin(0.009 * x + 0.013 * y);
            value += 12.0 * std::cos(0.012 * x - 0.010 * y);
            value += 6.0 * std::sin(0.047 * x) * std::cos(0.041 * y);

            const double ellipseMask = soft_ellipse(x, y, cx1, cy1, rx1, ry1, 4.0);
            const double diskMask = soft_disk(x, y, cx2, cy2, radius2, 3.0);
            const double rectMask = soft_rect(x, y, rectX0, rectY0, rectX1, rectY1, 3.5);

            value += ellipseMask * (78.0 + 8.0 * std::sin(0.08 * x + 0.03 * y));
            value += diskMask * (122.0 + 10.0 * std::cos(0.06 * x - 0.04 * y));
            value += rectMask * (54.0 + 5.0 * std::sin(0.05 * y));
            value += ((r + 2 * c) % 9);
            value += 3.5 * std::sin(0.19 * x + 0.11 * y);

            data[p] = clamp_to_u8(value);
        }
    }

    return image;
}

ImageUInt8Ptr make_benchmark_image(BenchmarkImageModel model, int numRows, int numCols) {
    switch (model) {
        case BenchmarkImageModel::structured:
            return make_structured_benchmark_image(numRows, numCols);
        case BenchmarkImageModel::natural_like:
            return make_natural_like_benchmark_image(numRows, numCols);
        case BenchmarkImageModel::piecewise_scene:
            return make_piecewise_scene_benchmark_image(numRows, numCols);
    }
    return make_structured_benchmark_image(numRows, numCols);
}

std::vector<int> make_area_thresholds(int numPixels, int numThresholds) {
    assert(numThresholds > 0);
    std::vector<int> thresholds;
    thresholds.reserve(static_cast<std::size_t>(numThresholds));
    for (int step = 1; step <= numThresholds; ++step) {
        const double ratio = static_cast<double>(step) / static_cast<double>(numThresholds + 1);
        const int threshold = std::max(1, static_cast<int>(std::llround(ratio * static_cast<double>(numPixels))));
        if (thresholds.empty() || threshold > thresholds.back()) {
            thresholds.push_back(threshold);
        }
    }
    return thresholds;
}

bool same_image(const ImageUInt8Ptr &lhs, const ImageUInt8Ptr &rhs) {
    if (lhs == nullptr || rhs == nullptr) {
        return lhs == rhs;
    }
    return lhs->isEqual(rhs);
}

std::vector<NodeId> get_nodes_threshold(DynamicComponentTree *tree,
                                        const std::vector<float> &attribute,
                                        int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree->getRoot());

    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (attribute[static_cast<std::size_t>(nodeId)] > threshold) {
            for (NodeId childId : tree->getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }

    return out;
}

ImageUInt8Ptr apply_naive_threshold_by_area(ImageUInt8Ptr image, double radius, int threshold) {
    auto adjacency = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radius);

    auto maxTreePtr = std::make_shared<DynamicComponentTree>(image, true, adjacency);
    auto *maxTree = maxTreePtr.get();
    DynamicAreaComputer maxAreaComputer(maxTree);
    std::vector<float> maxArea = maxAreaComputer.compute();
    for (NodeId nodeId : get_nodes_threshold(maxTree, maxArea, threshold)) {
        maxTree->pruneNode(nodeId);
    }
    ImageUInt8Ptr filtered = maxTree->reconstructionImage();

    auto minTreePtr = std::make_shared<DynamicComponentTree>(filtered, false, adjacency);
    auto *minTree = minTreePtr.get();
    DynamicAreaComputer minAreaComputer(minTree);
    std::vector<float> minArea = minAreaComputer.compute();
    for (NodeId nodeId : get_nodes_threshold(minTree, minArea, threshold)) {
        minTree->pruneNode(nodeId);
    }
    return minTree->reconstructionImage();
}

ImageUInt8Ptr run_naive_area_sequence(ImageUInt8Ptr image,
                                      double radius,
                                      const std::vector<int> &thresholds) {
    ImageUInt8Ptr current = image->clone();
    for (int threshold : thresholds) {
        current = apply_naive_threshold_by_area(current, radius, threshold);
    }
    return current;
}

template<typename Fn>
double time_seconds(Fn &&fn) {
    const auto start = std::chrono::steady_clock::now();
    fn();
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

struct AreaBenchmarkTiming {
    double casfSeconds = 0.0;
    double naiveSeconds = 0.0;
};

constexpr std::array<int, 4> kThresholdRounds = {2, 8, 16, 32};

std::vector<std::vector<int64_t>> make_benchmark_argument_product() {
    std::vector<int64_t> sizes = benchmark::CreateRange(32, 1024, 2);
    std::vector<int64_t> thresholds(kThresholdRounds.begin(), kThresholdRounds.end());
    return {std::move(sizes), std::move(thresholds)};
}

std::vector<int> get_dynamic_nodes_threshold(const DynamicComponentTree &tree,
                                             const std::vector<float> &area,
                                             int threshold) {
    std::vector<int> out;
    FastQueue<int> queue;
    queue.push(tree.getRoot());

    while (!queue.empty()) {
        const int nodeId = queue.pop();
        if (area[static_cast<std::size_t>(nodeId)] > threshold) {
            for (int childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }

    return out;
}

bool validate_area_buffer(const DynamicComponentTree &tree,
                          DynamicAreaComputer &areaComputer,
                          const std::vector<float> &incrementalArea,
                          const char *treeName,
                          int threshold,
                          std::string &error) {
    std::vector<float> fullArea(tree.getNumInternalNodeSlots(), 0.0f);
    areaComputer.compute(std::span<float>(fullArea));

    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }

        const auto index = static_cast<std::size_t>(nodeId);
        if (incrementalArea[index] != fullArea[index]) {
            std::ostringstream stream;
            stream << treeName
                   << " area mismatch after threshold=" << threshold
                   << " at node=" << nodeId
                   << " incremental=" << incrementalArea[index]
                   << " full=" << fullArea[index]
                   << " parent=" << tree.getNodeParent(nodeId)
                   << " children=" << tree.getNumChildren(nodeId)
                   << " proper_parts=" << tree.getNumProperParts(nodeId);
            error = stream.str();
            return false;
        }
    }

    return true;
}

class ComponentTreeCasfSubtree {
private:
    ImageUInt8Ptr input_;
    AdjacencyRelationPtr adjacency_;
    DynamicComponentTree maxTree_;
    DynamicComponentTree minTree_;
    DualMinMaxTreeIncrementalFilterInstrumented<AltitudeType> adjust_;
    DynamicAreaComputer maxAreaComputer_;
    DynamicAreaComputer minAreaComputer_;
    std::vector<float> maxArea_;
    std::vector<float> minArea_;

public:
    ComponentTreeCasfSubtree(double radius, ImageUInt8Ptr image)
        : input_(image->clone()),
          adjacency_(std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radius)),
          maxTree_(input_, true, adjacency_),
          minTree_(input_, false, adjacency_),
          adjust_(&minTree_, &maxTree_, *adjacency_),
          maxAreaComputer_(&maxTree_),
          minAreaComputer_(&minTree_),
          maxArea_(maxTree_.getNumInternalNodeSlots(), 0.0f),
          minArea_(minTree_.getNumInternalNodeSlots(), 0.0f) {
        maxAreaComputer_.compute(std::span<float>(maxArea_));
        minAreaComputer_.compute(std::span<float>(minArea_));
        adjust_.setAttributeComputer(minAreaComputer_,
                                     maxAreaComputer_,
                                     std::span<float>(minArea_),
                                     std::span<float>(maxArea_));
        adjust_.setRuntimePostConditionValidationEnabled(false);
    }

    ImageUInt8Ptr filter(const std::vector<int> &thresholds,
                         bool validateAttributes = false,
                         std::string *error = nullptr) {
        auto validate_all = [&](int threshold) -> bool {
            if (!validateAttributes) {
                return true;
            }

            std::string localError;
            if (!validate_area_buffer(maxTree_, maxAreaComputer_, maxArea_, "maxTree", threshold, localError)) {
                if (error != nullptr) {
                    *error = localError;
                }
                return false;
            }

            if (!validate_area_buffer(minTree_, minAreaComputer_, minArea_, "minTree", threshold, localError)) {
                if (error != nullptr) {
                    *error = localError;
                }
                return false;
            }

            return true;
        };

        if (!validate_all(-1)) {
            return nullptr;
        }

        for (int threshold : thresholds) {
            auto maxNodes = get_dynamic_nodes_threshold(maxTree_, maxArea_, threshold);
            adjust_.pruneMaxTreeAndUpdateMinTree(maxNodes);
            if (!validate_all(threshold)) {
                return nullptr;
            }

            auto minNodes = get_dynamic_nodes_threshold(minTree_, minArea_, threshold);
            adjust_.pruneMinTreeAndUpdateMaxTree(minNodes);
            if (!validate_all(threshold)) {
                return nullptr;
            }
        }
        return minTree_.reconstructionImage();
    }
};

AreaBenchmarkTiming measure_reference_timings(ImageUInt8Ptr image,
                                              double radius,
                                              const std::vector<int> &thresholds) {
    AreaBenchmarkTiming timing;
    timing.casfSeconds = time_seconds([&]() {
        ComponentTreeCasfSubtree runner(radius, image);
        const auto output = runner.filter(thresholds);
        benchmark::DoNotOptimize(output->rawData());
    });
    timing.naiveSeconds = time_seconds([&]() {
        const auto output = run_naive_area_sequence(image, radius, thresholds);
        benchmark::DoNotOptimize(output->rawData());
    });
    return timing;
}

bool validate_against_naive(ImageUInt8Ptr image,
                            double radius,
                            const std::vector<int> &thresholds,
                            std::string &error) {
    const auto expected = run_naive_area_sequence(image, radius, thresholds);
    ComponentTreeCasfSubtree runner(radius, image);
    const auto actual = runner.filter(thresholds, false, &error);
    if (actual == nullptr) {
        return false;
    }
    if (!same_image(expected, actual)) {
        error = "CASF subtree and naive outputs must match in the benchmark fixture.";
        return false;
    }
    return true;
}

static void BM_component_tree_casf_area(benchmark::State &state, BenchmarkImageModel model) {
    const int size = static_cast<int>(state.range(0));
    const int numThresholds = static_cast<int>(state.range(1));
    const int numPixels = size * size;
    const double radius = 1.0;
    const auto image = make_benchmark_image(model, size, size);
    const auto thresholds = make_area_thresholds(numPixels, numThresholds);
    const auto timing = measure_reference_timings(image, radius, thresholds);

    std::string error;
    if (!validate_against_naive(image, radius, thresholds, error)) {
        state.SkipWithError(error.c_str());
        return;
    }

    state.SetLabel(benchmark_case_label(model, numThresholds));
    state.counters["num_thresholds"] = benchmark::Counter(static_cast<double>(numThresholds), benchmark::Counter::kDefaults, benchmark::Counter::OneK::kIs1000);
    state.counters["speedup_vs_naive"] = benchmark::Counter(timing.naiveSeconds / timing.casfSeconds);

    for (auto _ : state) {
        ComponentTreeCasfSubtree runner(radius, image);
        const auto output = runner.filter(thresholds);
        benchmark::DoNotOptimize(output->rawData());
        benchmark::ClobberMemory();
    }
}

static void BM_component_tree_casf_naive_area(benchmark::State &state, BenchmarkImageModel model) {
    const int size = static_cast<int>(state.range(0));
    const int numThresholds = static_cast<int>(state.range(1));
    const int numPixels = size * size;
    const double radius = 1.0;
    const auto image = make_benchmark_image(model, size, size);
    const auto thresholds = make_area_thresholds(numPixels, numThresholds);
    const auto timing = measure_reference_timings(image, radius, thresholds);

    std::string error;
    if (!validate_against_naive(image, radius, thresholds, error)) {
        state.SkipWithError(error.c_str());
        return;
    }

    state.SetLabel(benchmark_case_label(model, numThresholds));
    state.counters["num_thresholds"] = benchmark::Counter(static_cast<double>(numThresholds), benchmark::Counter::kDefaults, benchmark::Counter::OneK::kIs1000);
    state.counters["speedup_vs_naive"] = benchmark::Counter(timing.naiveSeconds / timing.naiveSeconds);

    for (auto _ : state) {
        const auto output = run_naive_area_sequence(image, radius, thresholds);
        benchmark::DoNotOptimize(output->rawData());
        benchmark::ClobberMemory();
    }
}

// Threshold rounds span a short sequence, a moderate CASF stack, and deeper
// runs where the incremental strategy tends to amortize its setup cost.
BENCHMARK_CAPTURE(BM_component_tree_casf_area, structured, BenchmarkImageModel::structured)->ArgsProduct(make_benchmark_argument_product());
BENCHMARK_CAPTURE(BM_component_tree_casf_area, natural_like, BenchmarkImageModel::natural_like)->ArgsProduct(make_benchmark_argument_product());
BENCHMARK_CAPTURE(BM_component_tree_casf_area, piecewise_scene, BenchmarkImageModel::piecewise_scene)->ArgsProduct(make_benchmark_argument_product());

BENCHMARK_CAPTURE(BM_component_tree_casf_naive_area, structured, BenchmarkImageModel::structured)->ArgsProduct(make_benchmark_argument_product());
BENCHMARK_CAPTURE(BM_component_tree_casf_naive_area, natural_like, BenchmarkImageModel::natural_like)->ArgsProduct(make_benchmark_argument_product());
BENCHMARK_CAPTURE(BM_component_tree_casf_naive_area, piecewise_scene, BenchmarkImageModel::piecewise_scene)->ArgsProduct(make_benchmark_argument_product());

} // namespace
