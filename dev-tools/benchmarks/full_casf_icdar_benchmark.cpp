#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "../tools/dynamic_casf_apply_common.hpp"
#include "../../morphoTreeAdjust/include/ComponentTreeCasfNonIncreasing.hpp"

namespace fs = std::filesystem;
namespace common = dynamic_casf_apply_common;

using Clock = std::chrono::steady_clock;

struct Options {
    fs::path datasetRoot;
    fs::path outputDir = "results/full_casf_icdar_benchmark";
    fs::path rectangularityThresholdFile =
        "results/non_increasing_bbox_attributes_benchmark/thresholds/rectangularity_50_thresholds.txt";
    fs::path boxAspectThresholdFile =
        "results/non_increasing_bbox_attributes_benchmark/thresholds/box_aspect_ratio_50_thresholds.txt";
    std::vector<std::string> resolutions = {"icdar_480p", "icdar_720p", "icdar_1080p"};
    std::vector<std::string> methods = {
        "classic_old_legacy",
        "classic_new_public",
        "non_increasing_rectangularity",
        "non_increasing_box_aspect_ratio"};
    std::string imageGlob = "*.png";
    std::vector<int> classicThresholds;
    std::map<std::string, std::vector<int>> classicThresholdsByResolution;
    std::vector<float> rectangularityThresholds;
    std::vector<float> boxAspectThresholds;
    fs::path rawCsvPath;
    std::size_t imageLimit = 0;
    int repeat = 1;
    int warmup = 0;
    double radioAdj = 1.0;
    bool resume = true;
};

struct TimedImage {
    ImageUInt8Ptr image;
    double ms = 0.0;
};

struct Timing {
    double initMs = 0.0;
    double filterMs = 0.0;
    double reconstructMs = 0.0;
    double totalMs = 0.0;
    double wallMs = 0.0;
    ImageUInt8Ptr output;
    std::string error;
};

struct ImageStats {
    long long changedPixels = 0;
    double meanAbs = 0.0;
    int uniqueGrayLevels = 0;
    std::uint64_t checksum = 1469598103934665603ull;
};

static double elapsedMs(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

template<class Fn>
static double timeCall(Fn &&fn) {
    const auto start = Clock::now();
    fn();
    return elapsedMs(start, Clock::now());
}

static std::vector<std::string> split(std::string_view text, char separator) {
    std::vector<std::string> out;
    std::string current;
    for (char ch : text) {
        if (ch == separator) {
            if (!current.empty()) {
                out.push_back(current);
            }
            current.clear();
        } else if (!std::isspace(static_cast<unsigned char>(ch))) {
            current.push_back(ch);
        }
    }
    if (!current.empty()) {
        out.push_back(current);
    }
    return out;
}

static std::string trimCopy(std::string_view text) {
    std::size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start]))) {
        ++start;
    }
    std::size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    return std::string(text.substr(start, end - start));
}

static std::string toLowerCopy(std::string_view text) {
    std::string out(text);
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return out;
}

static std::string normalizeListSeparators(std::string text) {
    std::replace(text.begin(), text.end(), ';', ',');
    return text;
}

static bool parseBool(std::string_view text) {
    const std::string value = toLowerCopy(trimCopy(text));
    if (value == "true" || value == "1" || value == "yes") {
        return true;
    }
    if (value == "false" || value == "0" || value == "no") {
        return false;
    }
    throw std::runtime_error("invalid boolean value: " + std::string(text));
}

static bool isKnownMethod(const std::string &method) {
    return method == "classic_old_legacy" ||
           method == "classic_new_public" ||
           method == "non_increasing_rectangularity" ||
           method == "non_increasing_box_aspect_ratio";
}

static bool methodEnabled(const Options &options, const std::string &method) {
    return std::find(options.methods.begin(), options.methods.end(), method) != options.methods.end();
}

static std::vector<int> defaultClassicThresholds() {
    std::vector<int> thresholds;
    thresholds.reserve(50);
    for (int threshold = 50; threshold <= 2500; threshold += 50) {
        thresholds.push_back(threshold);
    }
    return thresholds;
}

static std::vector<int> scaleJmivAreaThresholds(double factor) {
    std::vector<int> thresholds;
    thresholds.reserve(50);
    for (int threshold = 50; threshold <= 2500; threshold += 50) {
        thresholds.push_back(static_cast<int>(std::lround(static_cast<double>(threshold) * factor)));
    }
    return thresholds;
}

static std::map<std::string, std::vector<int>> defaultJmivThresholdsByResolution() {
    return {
        {"icdar_1080p", defaultClassicThresholds()},
        {"icdar_720p", scaleJmivAreaThresholds(4.0 / 9.0)},
        {"icdar_480p", scaleJmivAreaThresholds(16.0 / 81.0)}};
}

static std::vector<int> parseIntCsv(const std::string &text) {
    std::vector<int> values;
    for (const auto &part : split(normalizeListSeparators(text), ',')) {
        values.push_back(std::stoi(part));
    }
    if (values.empty()) {
        throw std::runtime_error("empty integer threshold list");
    }
    return values;
}

static std::vector<float> parseFloatCsv(const std::string &text) {
    std::vector<float> values;
    for (const auto &part : split(normalizeListSeparators(text), ',')) {
        values.push_back(std::stof(part));
    }
    if (values.empty()) {
        throw std::runtime_error("empty float threshold list");
    }
    return values;
}

static std::string readTextFile(const fs::path &path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot read threshold file: " + path.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

static std::vector<float> rectangularityFallbackThresholds() {
    std::vector<float> thresholds;
    thresholds.reserve(50);
    for (int i = 1; i <= 50; ++i) {
        thresholds.push_back(static_cast<float>(0.004 * i));
    }
    return thresholds;
}

static std::vector<float> boxAspectFallbackThresholds() {
    std::vector<float> thresholds;
    thresholds.reserve(50);
    for (int i = 1; i <= 50; ++i) {
        thresholds.push_back(static_cast<float>(1.0 / (1.0 - 0.002 * i)));
    }
    return thresholds;
}

static std::vector<float> loadThresholdsOrFallback(const fs::path &path,
                                                   const std::vector<float> &fallback) {
    if (!fs::exists(path)) {
        return fallback;
    }
    return parseFloatCsv(readTextFile(path));
}

static std::string joinInts(const std::vector<int> &values) {
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            out << ';';
        }
        out << values[i];
    }
    return out.str();
}

static std::string joinFloats(const std::vector<float> &values) {
    std::ostringstream out;
    out << std::setprecision(7);
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            out << ';';
        }
        out << values[i];
    }
    return out.str();
}

static bool globMatch(std::string_view pattern, std::string_view text) {
    std::size_t p = 0;
    std::size_t t = 0;
    std::size_t star = std::string::npos;
    std::size_t match = 0;
    while (t < text.size()) {
        if (p < pattern.size() && (pattern[p] == '?' || pattern[p] == text[t])) {
            ++p;
            ++t;
            continue;
        }
        if (p < pattern.size() && pattern[p] == '*') {
            star = p++;
            match = t;
            continue;
        }
        if (star != std::string::npos) {
            p = star + 1;
            t = ++match;
            continue;
        }
        return false;
    }
    while (p < pattern.size() && pattern[p] == '*') {
        ++p;
    }
    return p == pattern.size();
}

static std::vector<fs::path> collectImages(const fs::path &dir,
                                           const std::string &imageGlob,
                                           std::size_t limit) {
    std::vector<fs::path> images;
    for (const auto &entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const std::string extension = entry.path().extension().string();
        if ((extension == ".png" || extension == ".jpg" || extension == ".jpeg" ||
             extension == ".tif" || extension == ".tiff") &&
            globMatch(imageGlob, entry.path().filename().string())) {
            images.push_back(entry.path());
        }
    }
    std::sort(images.begin(), images.end());
    if (limit > 0 && images.size() > limit) {
        images.resize(limit);
    }
    return images;
}

static std::size_t bufferCapacity(const DynamicComponentTree &tree) {
    const auto slots = static_cast<std::size_t>(
        std::max(0, tree.getNumInternalNodeSlots()));
    const auto pixels = static_cast<std::size_t>(
        std::max(1, tree.getNumTotalProperParts()));
    return slots + pixels;
}

static std::vector<float> computeAreaBuffer(DynamicComponentTree &tree,
                                            DynamicAreaComputer &computer) {
    std::vector<float> buffer(bufferCapacity(tree), 0.0f);
    computer.compute(std::span<float>(buffer));
    return buffer;
}

static std::vector<NodeId> getNodesToPrune(const DynamicComponentTree &tree,
                                           const std::vector<float> &area,
                                           int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree.getRoot());
    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (nodeId == InvalidNode || !tree.isNode(nodeId) || !tree.isAlive(nodeId)) {
            continue;
        }
        if (area[static_cast<std::size_t>(nodeId)] > threshold) {
            for (NodeId childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else if (nodeId != tree.getRoot()) {
            out.push_back(nodeId);
        }
    }
    return out;
}

static std::vector<NodeId> currentPruneCandidates(const DynamicComponentTree &tree,
                                                  const std::vector<NodeId> &nodes) {
    std::vector<NodeId> out;
    out.reserve(nodes.size());
    for (NodeId nodeId : nodes) {
        if (nodeId != InvalidNode &&
            nodeId != tree.getRoot() &&
            tree.isNode(nodeId) &&
            tree.isAlive(nodeId)) {
            out.push_back(nodeId);
        }
    }
    return out;
}

class ClassicCasfRunner {
private:
    bool usePublicPath_ = false;
    ImageUInt8Ptr input_;
    AdjacencyRelationPtr adjacency_;
    std::unique_ptr<DynamicComponentTree> maxTree_;
    std::unique_ptr<DynamicComponentTree> minTree_;
    std::unique_ptr<DualMinMaxTreeIncrementalFilter<AltitudeType>> adjust_;
    std::unique_ptr<DynamicAreaComputer> maxAreaComputer_;
    std::unique_ptr<DynamicAreaComputer> minAreaComputer_;
    std::vector<float> maxArea_;
    std::vector<float> minArea_;

public:
    ClassicCasfRunner(ImageUInt8Ptr image, double radius, bool usePublicPath)
        : usePublicPath_(usePublicPath),
          input_(image->clone()),
          adjacency_(std::make_shared<AdjacencyRelation>(
              image->getNumRows(),
              image->getNumCols(),
              radius)),
          maxTree_(std::make_unique<DynamicComponentTree>(input_->clone(), true, adjacency_)),
          minTree_(std::make_unique<DynamicComponentTree>(input_->clone(), false, adjacency_)),
          adjust_(std::make_unique<DualMinMaxTreeIncrementalFilter<AltitudeType>>(
              minTree_.get(),
              maxTree_.get(),
              *adjacency_)),
          maxAreaComputer_(std::make_unique<DynamicAreaComputer>(maxTree_.get())),
          minAreaComputer_(std::make_unique<DynamicAreaComputer>(minTree_.get())),
          maxArea_(computeAreaBuffer(*maxTree_, *maxAreaComputer_)),
          minArea_(computeAreaBuffer(*minTree_, *minAreaComputer_)) {
        adjust_->setAttributeComputer(
            *minAreaComputer_,
            *maxAreaComputer_,
            std::span<float>(minArea_),
            std::span<float>(maxArea_));
    }

    void applyThreshold(int threshold) {
        auto maxNodes = getNodesToPrune(*maxTree_, maxArea_, threshold);
        if (usePublicPath_) {
            maxNodes = currentPruneCandidates(*maxTree_, maxNodes);
            adjust_->pruneSubtreesAndUpdateDual(maxTree_.get(), minTree_.get(), std::span<const NodeId>(maxNodes));
        } else {
            adjust_->pruneMaxTreeAndUpdateMinTree(maxNodes);
        }

        auto minNodes = getNodesToPrune(*minTree_, minArea_, threshold);
        if (usePublicPath_) {
            minNodes = currentPruneCandidates(*minTree_, minNodes);
            adjust_->pruneSubtreesAndUpdateDual(minTree_.get(), maxTree_.get(), std::span<const NodeId>(minNodes));
        } else {
            adjust_->pruneMinTreeAndUpdateMaxTree(minNodes);
        }
    }

    ImageUInt8Ptr reconstructImage() const {
        return minTree_->reconstructionImage();
    }
};

static Timing runClassic(ImageUInt8Ptr image,
                         double radius,
                         const std::vector<int> &thresholds,
                         bool usePublicPath) {
    Timing timing;
    const auto wallStart = Clock::now();
    try {
        std::unique_ptr<ClassicCasfRunner> runner;
        timing.initMs = timeCall([&]() {
            runner = std::make_unique<ClassicCasfRunner>(image, radius, usePublicPath);
        });
        timing.filterMs = timeCall([&]() {
            for (int threshold : thresholds) {
                runner->applyThreshold(threshold);
            }
        });
        timing.reconstructMs = timeCall([&]() {
            timing.output = runner->reconstructImage();
        });
        timing.totalMs = timing.initMs + timing.filterMs + timing.reconstructMs;
    } catch (const std::exception &e) {
        timing.error = e.what();
    }
    timing.wallMs = elapsedMs(wallStart, Clock::now());
    return timing;
}

static Timing runNonIncreasing(ImageUInt8Ptr image,
                               double radius,
                               Attribute attribute,
                               const std::vector<float> &thresholds) {
    Timing timing;
    const auto wallStart = Clock::now();
    try {
        std::unique_ptr<ComponentTreeCasfNonIncreasing<AltitudeType>> runner;
        timing.initMs = timeCall([&]() {
            runner = std::make_unique<ComponentTreeCasfNonIncreasing<AltitudeType>>(
                image,
                radius,
                attribute);
        });
        timing.filterMs = timeCall([&]() {
            timing.output = runner->filter(std::span<const float>(thresholds));
        });
        timing.reconstructMs = 0.0;
        timing.totalMs = timing.initMs + timing.filterMs;
    } catch (const std::exception &e) {
        timing.error = e.what();
    }
    timing.wallMs = elapsedMs(wallStart, Clock::now());
    return timing;
}

static ImageStats computeImageStats(const ImageUInt8Ptr &input, const ImageUInt8Ptr &output) {
    ImageStats stats;
    if (input == nullptr || output == nullptr ||
        input->getNumRows() != output->getNumRows() ||
        input->getNumCols() != output->getNumCols()) {
        return stats;
    }

    bool seen[256] = {};
    const int numPixels = input->getNumRows() * input->getNumCols();
    const auto *in = input->rawData();
    const auto *out = output->rawData();
    long long absSum = 0;
    for (int i = 0; i < numPixels; ++i) {
        const int delta = static_cast<int>(out[i]) - static_cast<int>(in[i]);
        if (delta != 0) {
            ++stats.changedPixels;
        }
        absSum += std::llabs(static_cast<long long>(delta));
        seen[out[i]] = true;
        stats.checksum ^= static_cast<std::uint64_t>(out[i]);
        stats.checksum *= 1099511628211ull;
    }
    for (bool value : seen) {
        if (value) {
            ++stats.uniqueGrayLevels;
        }
    }
    stats.meanAbs = numPixels > 0 ? static_cast<double>(absSum) / static_cast<double>(numPixels) : 0.0;
    return stats;
}

static std::string csvEscape(const std::string &text) {
    if (text.find_first_of(",\"\n") == std::string::npos) {
        return text;
    }
    std::string escaped = "\"";
    for (char ch : text) {
        if (ch == '"') {
            escaped += "\"\"";
        } else {
            escaped += ch;
        }
    }
    escaped += '"';
    return escaped;
}

static void writeRawHeader(std::ofstream &out) {
    out << "group,method,attribute,resolution,image_index,image,image_path,rows,cols,"
           "threshold_count,thresholds,radio_adj,init_ms,filter_ms,reconstruction_ms,total_ms,wall_ms,"
           "success,error,classic_matches,changed_pixels,mean_abs_from_original,unique_gray_levels,checksum\n";
}

static void writeRawRow(std::ofstream &out,
                        const std::string &group,
                        const std::string &method,
                        const std::string &attribute,
                        const std::string &resolution,
                        std::size_t imageIndex,
                        const fs::path &imagePath,
                        int rows,
                        int cols,
                        std::size_t thresholdCount,
                        const std::string &thresholds,
                        double radioAdj,
                        const Timing &timing,
                        const ImageStats &stats,
                        const std::string &classicMatches) {
    out << group << ','
        << method << ','
        << attribute << ','
        << resolution << ','
        << imageIndex << ','
        << imagePath.filename().string() << ','
        << csvEscape(imagePath.string()) << ','
        << rows << ','
        << cols << ','
        << thresholdCount << ','
        << csvEscape(thresholds) << ','
        << radioAdj << ','
        << timing.initMs << ','
        << timing.filterMs << ','
        << timing.reconstructMs << ','
        << timing.totalMs << ','
        << timing.wallMs << ','
        << (timing.error.empty() ? "true" : "false") << ','
        << csvEscape(timing.error) << ','
        << classicMatches << ','
        << stats.changedPixels << ','
        << stats.meanAbs << ','
        << stats.uniqueGrayLevels << ','
        << stats.checksum
        << '\n';
    out.flush();
}

static std::set<std::string> readCompletedKeys(const fs::path &rawPath) {
    std::set<std::string> keys;
    std::ifstream input(rawPath);
    if (!input) {
        return keys;
    }
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> columns = split(line, ',');
        if (columns.size() < 6) {
            continue;
        }
        keys.insert(columns[1] + "|" + columns[3] + "|" + columns[4]);
    }
    return keys;
}

static bool hasCompleted(const std::set<std::string> &keys,
                         const std::string &method,
                         const std::string &resolution,
                         std::size_t imageIndex) {
    return keys.count(method + "|" + resolution + "|" + std::to_string(imageIndex)) != 0;
}

static void parseMethodList(const std::string &value, Options &options) {
    options.methods.clear();
    for (const auto &method : split(value, ',')) {
        if (!isKnownMethod(method)) {
            throw std::runtime_error("unknown method in config: " + method);
        }
        options.methods.push_back(method);
    }
    if (options.methods.empty()) {
        throw std::runtime_error("empty methods list");
    }
}

static void parseConfigFile(const fs::path &path, Options &options) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open config: " + path.string());
    }

    std::string line;
    int lineNo = 0;
    while (std::getline(input, line)) {
        ++lineNo;
        const std::size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line.erase(commentPos);
        }
        const std::string trimmed = trimCopy(line);
        if (trimmed.empty()) {
            continue;
        }

        const std::size_t separatorPos = trimmed.find_first_of("=:");
        if (separatorPos == std::string::npos) {
            throw std::runtime_error("invalid config line " + std::to_string(lineNo));
        }

        const std::string key = toLowerCopy(trimCopy(trimmed.substr(0, separatorPos)));
        const std::string value = trimCopy(trimmed.substr(separatorPos + 1));

        if (key == "dataset_root") {
            options.datasetRoot = value;
        } else if (key == "output_dir") {
            options.outputDir = value;
        } else if (key == "raw_csv") {
            options.rawCsvPath = value;
        } else if (key == "resolutions") {
            options.resolutions = split(value, ',');
        } else if (key == "image_glob") {
            options.imageGlob = value;
        } else if (key == "image_limit") {
            options.imageLimit = static_cast<std::size_t>(std::stoll(value));
        } else if (key == "repeat") {
            options.repeat = std::stoi(value);
        } else if (key == "warmup") {
            options.warmup = std::stoi(value);
        } else if (key == "radio_adj") {
            options.radioAdj = std::stod(value);
        } else if (key == "resume") {
            options.resume = parseBool(value);
        } else if (key == "methods") {
            parseMethodList(value, options);
        } else if (key == "method") {
            if (!isKnownMethod(value)) {
                throw std::runtime_error("unknown method in config: " + value);
            }
            options.methods.push_back(value);
        } else if (key == "classic_attribute") {
            if (toLowerCopy(value) != "area") {
                throw std::runtime_error("full_casf_icdar_benchmark currently supports only classic_attribute=area");
            }
        } else if (key == "classic_threshold_policy") {
            // Documentational key. Explicit per-resolution lists below define the effective policy.
        } else if (key == "classic_thresholds") {
            options.classicThresholds = parseIntCsv(value);
            options.classicThresholdsByResolution.clear();
        } else if (key.rfind("classic_thresholds_", 0) == 0) {
            const std::string resolution = key.substr(std::string("classic_thresholds_").size());
            options.classicThresholdsByResolution[resolution] = parseIntCsv(value);
        } else if (key == "non_increasing_attributes") {
            // The selected attributes are represented by the method list.
        } else if (key == "rectangularity_thresholds") {
            options.rectangularityThresholds = parseFloatCsv(value);
        } else if (key == "box_aspect_ratio_thresholds") {
            options.boxAspectThresholds = parseFloatCsv(value);
        } else if (key == "rectangularity_threshold_file") {
            options.rectangularityThresholdFile = value;
        } else if (key == "box_aspect_threshold_file" || key == "box_aspect_ratio_threshold_file") {
            options.boxAspectThresholdFile = value;
        } else {
            throw std::runtime_error("unknown config key on line " + std::to_string(lineNo) + ": " + key);
        }
    }
}

static const std::vector<int> &classicThresholdsForResolution(const Options &options,
                                                              const std::string &resolution) {
    const auto it = options.classicThresholdsByResolution.find(resolution);
    if (it != options.classicThresholdsByResolution.end()) {
        return it->second;
    }
    return options.classicThresholds;
}

static Options parseOptions(int argc, char **argv) {
    Options options;
    options.classicThresholds = defaultClassicThresholds();
    options.classicThresholdsByResolution = defaultJmivThresholdsByResolution();
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto requireValue = [&](const char *name) -> std::string {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("missing value for ") + name);
            }
            return argv[++i];
        };
        if (arg == "--config") {
            parseConfigFile(requireValue("--config"), options);
        } else if (arg == "--dataset-root") {
            options.datasetRoot = requireValue("--dataset-root");
        } else if (arg == "--output-dir") {
            options.outputDir = requireValue("--output-dir");
        } else if (arg == "--raw-csv") {
            options.rawCsvPath = requireValue("--raw-csv");
        } else if (arg == "--image-limit") {
            options.imageLimit = static_cast<std::size_t>(std::stoll(requireValue("--image-limit")));
        } else if (arg == "--radio-adj") {
            options.radioAdj = std::stod(requireValue("--radio-adj"));
        } else if (arg == "--resolutions") {
            options.resolutions = split(requireValue("--resolutions"), ',');
        } else if (arg == "--methods") {
            parseMethodList(requireValue("--methods"), options);
        } else if (arg == "--classic-thresholds") {
            options.classicThresholds = parseIntCsv(requireValue("--classic-thresholds"));
            options.classicThresholdsByResolution.clear();
        } else if (arg == "--classic-thresholds-for-resolution") {
            const auto parts = split(requireValue("--classic-thresholds-for-resolution"), ':');
            if (parts.size() != 2) {
                throw std::runtime_error("--classic-thresholds-for-resolution expects resolution:t1,t2,...");
            }
            options.classicThresholdsByResolution[parts[0]] = parseIntCsv(parts[1]);
        } else if (arg == "--rectangularity-thresholds") {
            options.rectangularityThresholds = parseFloatCsv(requireValue("--rectangularity-thresholds"));
        } else if (arg == "--box-aspect-thresholds") {
            options.boxAspectThresholds = parseFloatCsv(requireValue("--box-aspect-thresholds"));
        } else if (arg == "--rectangularity-threshold-file") {
            options.rectangularityThresholdFile = requireValue("--rectangularity-threshold-file");
        } else if (arg == "--box-aspect-threshold-file") {
            options.boxAspectThresholdFile = requireValue("--box-aspect-threshold-file");
        } else if (arg == "--no-resume") {
            options.resume = false;
        } else {
            throw std::runtime_error("unknown option: " + arg);
        }
    }
    return options;
}

int main(int argc, char **argv) {
    try {
        Options options = parseOptions(argc, argv);
        if (options.datasetRoot.empty()) {
            throw std::runtime_error(
                "dataset root is required; pass --dataset-root or set dataset_root in the config file");
        }
        const auto rectangularityThresholds = options.rectangularityThresholds.empty()
            ? loadThresholdsOrFallback(options.rectangularityThresholdFile, rectangularityFallbackThresholds())
            : options.rectangularityThresholds;
        const auto boxAspectThresholds = options.boxAspectThresholds.empty()
            ? loadThresholdsOrFallback(options.boxAspectThresholdFile, boxAspectFallbackThresholds())
            : options.boxAspectThresholds;

        fs::create_directories(options.outputDir);
        const fs::path rawPath = options.rawCsvPath.empty()
            ? options.outputDir / "full_casf_icdar_raw.csv"
            : options.rawCsvPath;
        if (!rawPath.parent_path().empty()) {
            fs::create_directories(rawPath.parent_path());
        }
        const bool append = options.resume && fs::exists(rawPath);
        std::set<std::string> completed = append ? readCompletedKeys(rawPath) : std::set<std::string>{};

        std::ofstream raw(rawPath, append ? std::ios::app : std::ios::trunc);
        if (!raw) {
            throw std::runtime_error("cannot open output csv: " + rawPath.string());
        }
        raw << std::fixed << std::setprecision(6);
        if (!append) {
            writeRawHeader(raw);
        }

        std::ofstream manifest(options.outputDir / "manifest.txt", std::ios::trunc);
        manifest << "dataset_root=" << options.datasetRoot << '\n'
                 << "resolutions=";
        for (std::size_t i = 0; i < options.resolutions.size(); ++i) {
            if (i != 0) {
                manifest << ';';
            }
            manifest << options.resolutions[i];
        }
        manifest << '\n'
                 << "image_limit=" << options.imageLimit << '\n'
                 << "image_glob=" << options.imageGlob << '\n'
                 << "repeat=" << options.repeat << '\n'
                 << "warmup=" << options.warmup << '\n'
                 << "radio_adj=" << options.radioAdj << '\n'
                 << "raw_csv=" << rawPath << '\n'
                 << "rectangularity_thresholds=" << joinFloats(rectangularityThresholds) << '\n'
                 << "box_aspect_ratio_thresholds=" << joinFloats(boxAspectThresholds) << '\n';
        for (const auto &resolution : options.resolutions) {
            manifest << "classic_thresholds_" << resolution << '='
                     << joinInts(classicThresholdsForResolution(options, resolution)) << '\n';
        }

        std::size_t totalImages = 0;
        for (const auto &resolution : options.resolutions) {
            totalImages += collectImages(options.datasetRoot / resolution,
                                         options.imageGlob,
                                         options.imageLimit).size();
        }

        std::size_t globalImageCounter = 0;
        for (const auto &resolution : options.resolutions) {
            const auto images = collectImages(options.datasetRoot / resolution,
                                              options.imageGlob,
                                              options.imageLimit);
            const auto &classicThresholds = classicThresholdsForResolution(options, resolution);
            const std::string classicThresholdText = joinInts(classicThresholds);
            std::cerr << "resolution=" << resolution << " images=" << images.size() << '\n';
            for (std::size_t imageIndex = 0; imageIndex < images.size(); ++imageIndex) {
                const auto &imagePath = images[imageIndex];
                ++globalImageCounter;
                std::cerr << "[" << globalImageCounter << "/" << totalImages << "] "
                          << resolution << " " << imagePath.filename().string() << std::flush;

                ImageUInt8Ptr input = common::loadGrayImage(imagePath.string());
                const int rows = input->getNumRows();
                const int cols = input->getNumCols();

                Timing classicOld;
                Timing classicNew;
                if (methodEnabled(options, "classic_old_legacy") &&
                    !hasCompleted(completed, "classic_old_legacy", resolution, imageIndex)) {
                    classicOld = runClassic(input, options.radioAdj, classicThresholds, false);
                }
                if (methodEnabled(options, "classic_new_public") &&
                    !hasCompleted(completed, "classic_new_public", resolution, imageIndex)) {
                    classicNew = runClassic(input, options.radioAdj, classicThresholds, true);
                }

                std::string classicMatches;
                if (classicOld.output != nullptr && classicNew.output != nullptr) {
                    classicMatches = classicOld.output->isEqual(classicNew.output) ? "true" : "false";
                }

                if (methodEnabled(options, "classic_old_legacy") &&
                    !hasCompleted(completed, "classic_old_legacy", resolution, imageIndex)) {
                    writeRawRow(raw,
                                "classic",
                                "classic_old_legacy",
                                "area",
                                resolution,
                                imageIndex,
                                imagePath,
                                rows,
                                cols,
                                classicThresholds.size(),
                                classicThresholdText,
                                options.radioAdj,
                                classicOld,
                                computeImageStats(input, classicOld.output),
                                classicMatches);
                }

                if (methodEnabled(options, "classic_new_public") &&
                    !hasCompleted(completed, "classic_new_public", resolution, imageIndex)) {
                    writeRawRow(raw,
                                "classic",
                                "classic_new_public",
                                "area",
                                resolution,
                                imageIndex,
                                imagePath,
                                rows,
                                cols,
                                classicThresholds.size(),
                                classicThresholdText,
                                options.radioAdj,
                                classicNew,
                                computeImageStats(input, classicNew.output),
                                classicMatches);
                }

                if (methodEnabled(options, "non_increasing_rectangularity") &&
                    !hasCompleted(completed, "non_increasing_rectangularity", resolution, imageIndex)) {
                    Timing rectangularity = runNonIncreasing(
                        input,
                        options.radioAdj,
                        RECTANGULARITY,
                        rectangularityThresholds);
                    writeRawRow(raw,
                                "non_increasing",
                                "non_increasing_rectangularity",
                                "rectangularity",
                                resolution,
                                imageIndex,
                                imagePath,
                                rows,
                                cols,
                                rectangularityThresholds.size(),
                                joinFloats(rectangularityThresholds),
                                options.radioAdj,
                                rectangularity,
                                computeImageStats(input, rectangularity.output),
                                "");
                }

                if (methodEnabled(options, "non_increasing_box_aspect_ratio") &&
                    !hasCompleted(completed, "non_increasing_box_aspect_ratio", resolution, imageIndex)) {
                    Timing aspectRatio = runNonIncreasing(
                        input,
                        options.radioAdj,
                        BOX_ASPECT_RATIO,
                        boxAspectThresholds);
                    writeRawRow(raw,
                                "non_increasing",
                                "non_increasing_box_aspect_ratio",
                                "box_aspect_ratio",
                                resolution,
                                imageIndex,
                                imagePath,
                                rows,
                                cols,
                                boxAspectThresholds.size(),
                                joinFloats(boxAspectThresholds),
                                options.radioAdj,
                                aspectRatio,
                                computeImageStats(input, aspectRatio.output),
                                "");
                }

                std::cerr << " done\n";
            }
        }

        std::cerr << "wrote " << rawPath << '\n';
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
