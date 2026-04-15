#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "../../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../../morphoTreeAdjust/include/Common.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTreeAdjustmentInstrumented.hpp"

#include "../external/stb/stb_image.h"

namespace fs = std::filesystem;

enum class AttributeMode {
    Area,
    BBoxWidth,
    BBoxHeight,
    BBoxDiagonal
};

struct BenchOptions {
    std::vector<std::string> images;
    std::vector<std::string> methods;
    std::vector<int> thresholds;
    AttributeMode attributeMode = AttributeMode::Area;
    int repeat = 1;
    int warmup = 0;
    double radioAdj = 1.5;
    bool quiet = false;
    bool json = false;
    bool noValidate = false;
};

struct StatsSummary {
    double mean = 0.0;
    double median = 0.0;
    double stddev = 0.0;
    double min = 0.0;
    double max = 0.0;
};

struct ThresholdMetrics {
    int numDiscardedNodes = 0;
    int flatzones = 0;
    int inputTreeNodes = 0;
    long long numPixelBorders = 0;
    long long area = 0;
    int numSubtrees = 0;
    bool valid = false;
};

struct IterationSample {
    int threshold = 0;
    double phase1Ms = 0.0;
    double phase2Ms = 0.0;
    double totalMs = 0.0;
    double phase1UpdateMs = 0.0;
    double phase1PruningMs = 0.0;
    double phase2UpdateMs = 0.0;
    double phase2PruningMs = 0.0;
    ThresholdMetrics phase1InputMetrics;
    ThresholdMetrics phase2InputMetrics;
    DynamicSubtreeMetrics phase1AlgMetrics;
    DynamicSubtreeMetrics phase2AlgMetrics;
};

struct IterationAggregate {
    int index = 0;
    int threshold = 0;
    StatsSummary phase1;
    StatsSummary phase2;
    StatsSummary total;
    StatsSummary phase1Update;
    StatsSummary phase1Pruning;
    StatsSummary phase2Update;
    StatsSummary phase2Pruning;
    ThresholdMetrics phase1InputMetrics;
    ThresholdMetrics phase2InputMetrics;
    DynamicSubtreeMetrics phase1AlgMetrics;
    DynamicSubtreeMetrics phase2AlgMetrics;
    std::size_t count = 0;
};

struct MethodResult {
    std::string name;
    std::vector<double> totalSamplesMs;
    std::vector<double> buildSamplesMs;
    std::vector<double> buildMaxTreeSamplesMs;
    std::vector<double> buildMinTreeSamplesMs;
    std::vector<IterationAggregate> iterations;
    ImageUInt8Ptr output;
};

struct ImageMetrics {
    int rows = 0;
    int cols = 0;
    int numFlatzonesImage = 0;
    long long numPixelBordersImage = 0;
    int numNodesMaxTree = 0;
    int numNodesMinTree = 0;
};

struct MethodImageResult {
    std::string imagePath;
    ImageMetrics imageMetrics;
    MethodResult result;
};

struct MethodBenchmarkResult {
    std::string name;
    std::vector<MethodImageResult> images;
};

struct FlatzonePartition {
    std::vector<int> labelByPixel;
    int count = 0;
};

struct TimerHooksContext {
    Stopwatch *phase = nullptr;
    Stopwatch *iter = nullptr;
    Stopwatch *total = nullptr;
};

static inline double elapsedMs(const Stopwatch &sw) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(sw.elapsed()).count();
}

static inline double elapsedUsToMsRounded(const Stopwatch &sw) {
    const auto us = std::chrono::duration_cast<std::chrono::microseconds>(sw.elapsed()).count();
    return static_cast<double>(us) / 1000.0;
}

static void pauseTimerHooks(void *ctx) {
    auto *timers = static_cast<TimerHooksContext *>(ctx);
    if (timers == nullptr) {
        return;
    }
    if (timers->phase != nullptr && timers->phase->running()) {
        timers->phase->pause();
    }
    if (timers->iter != nullptr && timers->iter->running()) {
        timers->iter->pause();
    }
    if (timers->total != nullptr && timers->total->running()) {
        timers->total->pause();
    }
}

static void resumeTimerHooks(void *ctx) {
    auto *timers = static_cast<TimerHooksContext *>(ctx);
    if (timers == nullptr) {
        return;
    }
    if (timers->phase != nullptr && !timers->phase->running()) {
        timers->phase->resume();
    }
    if (timers->iter != nullptr && !timers->iter->running()) {
        timers->iter->resume();
    }
    if (timers->total != nullptr && !timers->total->running()) {
        timers->total->resume();
    }
}

static void printUsage(const char *argv0) {
    std::cout << "Usage: " << argv0 << " [options] <image>...\n"
              << "Options:\n"
              << "  -i, --image <path>   Add image path (can repeat)\n"
              << "  -r, --repeat <N>     Repetitions per image (default 1)\n"
              << "  -w, --warmup <N>     Warmup runs per image (default 0)\n"
              << "  -q, --quiet          Reduce console output\n"
              << "  --json               Emit hierarchical JSON summary\n"
              << "  --config <path>      Load options from file\n"
              << "  --attribute <name>   Attribute: area|bbox_width|bbox_height|bbox_diagonal\n"
              << "  --method <name>      Method: naive|our_subtree\n"
              << "  --thresholds <list>  Comma-separated thresholds\n"
              << "  --no-validate        Skip output validation\n"
              << "  -h, --help           Show this help\n";
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
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return out;
}

static bool parseNonNegativeInt(const char *text, int *out) {
    if (text == nullptr || *text == '\0') {
        return false;
    }
    char *end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (end == nullptr || *end != '\0' || value < 0 || value > std::numeric_limits<int>::max()) {
        return false;
    }
    *out = static_cast<int>(value);
    return true;
}

static bool parseBool(std::string_view text, bool *out) {
    const std::string value = toLowerCopy(trimCopy(text));
    if (value == "true" || value == "1" || value == "yes") {
        *out = true;
        return true;
    }
    if (value == "false" || value == "0" || value == "no") {
        *out = false;
        return true;
    }
    return false;
}

static bool parseAttributeMode(std::string_view text, AttributeMode *out) {
    const std::string value = toLowerCopy(trimCopy(text));
    if (value == "area") {
        *out = AttributeMode::Area;
        return true;
    }
    if (value == "bbox_width" || value == "box_width") {
        *out = AttributeMode::BBoxWidth;
        return true;
    }
    if (value == "bbox_height" || value == "box_height") {
        *out = AttributeMode::BBoxHeight;
        return true;
    }
    if (value == "bbox_diagonal" || value == "bbox" || value == "box_diagonal") {
        *out = AttributeMode::BBoxDiagonal;
        return true;
    }
    return false;
}

static bool isKnownMethod(std::string_view method) {
    return method == "naive" || method == "our_subtree";
}

static std::vector<std::string> splitCommaSeparated(std::string_view text) {
    std::vector<std::string> values;
    std::stringstream ss{std::string(text)};
    std::string item;
    while (std::getline(ss, item, ',')) {
        const std::string trimmed = trimCopy(item);
        if (!trimmed.empty()) {
            values.push_back(trimmed);
        }
    }
    return values;
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

static bool hasGlobChars(std::string_view value) {
    return value.find_first_of("*?") != std::string_view::npos;
}

static bool expandGlobPath(const std::string &pattern, std::vector<std::string> &out) {
    if (!hasGlobChars(pattern)) {
        out.push_back(pattern);
        return true;
    }

    fs::path path(pattern);
    fs::path dir = path.parent_path();
    if (dir.empty()) {
        dir = ".";
    }

    std::error_code ec;
    if (!fs::exists(dir, ec)) {
        return false;
    }

    const std::string filenamePattern = path.filename().string();
    for (const auto &entry : fs::directory_iterator(dir, ec)) {
        if (ec) {
            return false;
        }
        if (globMatch(filenamePattern, entry.path().filename().string())) {
            out.push_back(entry.path().string());
        }
    }
    std::sort(out.begin(), out.end());
    return true;
}

static bool parseThresholdList(std::string_view text, std::vector<int> *out) {
    out->clear();
    for (const auto &item : splitCommaSeparated(text)) {
        int value = 0;
        if (!parseNonNegativeInt(item.c_str(), &value)) {
            return false;
        }
        out->push_back(value);
    }
    return !out->empty();
}

static bool parseConfigFile(const std::string &path, BenchOptions &options) {
    std::ifstream in(path);
    if (!in) {
        std::cerr << "Error: could not open config " << path << "\n";
        return false;
    }

    const fs::path baseDir = fs::absolute(fs::path(path)).parent_path();
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        const std::size_t hashPos = line.find('#');
        if (hashPos != std::string::npos) {
            line.erase(hashPos);
        }
        const std::string trimmed = trimCopy(line);
        if (trimmed.empty()) {
            continue;
        }

        const std::size_t sepPos = trimmed.find_first_of("=:");
        if (sepPos == std::string::npos) {
            std::cerr << "Error: invalid line in config (" << lineNo << ")\n";
            return false;
        }

        const std::string key = toLowerCopy(trimCopy(trimmed.substr(0, sepPos)));
        const std::string value = trimCopy(trimmed.substr(sepPos + 1));

        if (key == "image") {
            std::vector<std::string> expanded;
            const fs::path resolved = baseDir / value;
            if (!expandGlobPath(resolved.string(), expanded)) {
                std::cerr << "Error: could not expand image entry in config (" << lineNo << ")\n";
                return false;
            }
            options.images.insert(options.images.end(), expanded.begin(), expanded.end());
            continue;
        }
        if (key == "method") {
            if (!isKnownMethod(value)) {
                std::cerr << "Erro: metodo desconhecido no config (" << lineNo << ")\n";
                return false;
            }
            options.methods.push_back(value);
            continue;
        }
        if (key == "methods") {
            for (const auto &method : splitCommaSeparated(value)) {
                if (!isKnownMethod(method)) {
                    std::cerr << "Erro: metodo desconhecido no config (" << lineNo << ")\n";
                    return false;
                }
                options.methods.push_back(method);
            }
            continue;
        }
        if (key == "thresholds") {
            if (!parseThresholdList(value, &options.thresholds)) {
                std::cerr << "Erro: thresholds invalido no config (" << lineNo << ")\n";
                return false;
            }
            continue;
        }
        if (key == "repeat") {
            if (!parseNonNegativeInt(value.c_str(), &options.repeat)) {
                std::cerr << "Erro: repeat invalido no config (" << lineNo << ")\n";
                return false;
            }
            continue;
        }
        if (key == "warmup") {
            if (!parseNonNegativeInt(value.c_str(), &options.warmup)) {
                std::cerr << "Erro: warmup invalido no config (" << lineNo << ")\n";
                return false;
            }
            continue;
        }
        if (key == "radio_adj") {
            char *end = nullptr;
            options.radioAdj = std::strtod(value.c_str(), &end);
            if (end == nullptr || *end != '\0') {
                std::cerr << "Erro: radio_adj invalido no config (" << lineNo << ")\n";
                return false;
            }
            continue;
        }
        if (key == "attribute") {
            if (!parseAttributeMode(value, &options.attributeMode)) {
                std::cerr << "Erro: attribute invalido no config (" << lineNo << ")\n";
                return false;
            }
            continue;
        }
        if (key == "quiet") {
            if (!parseBool(value, &options.quiet)) {
                return false;
            }
            continue;
        }
        if (key == "json") {
            if (!parseBool(value, &options.json)) {
                return false;
            }
            continue;
        }
        if (key == "no_validate") {
            if (!parseBool(value, &options.noValidate)) {
                return false;
            }
            continue;
        }
        if (key == "validate_only") {
            bool ignored = false;
            if (!parseBool(value, &ignored)) {
                return false;
            }
            continue;
        }

        std::cerr << "Erro: chave desconhecida no config (" << lineNo << "): " << key << "\n";
        return false;
    }

    return true;
}

static bool parseCommandLine(int argc, char **argv, BenchOptions &options) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return false;
        }
        if (arg == "-i" || arg == "--image") {
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para " << arg << "\n";
                return false;
            }
            options.images.push_back(argv[++i]);
            continue;
        }
        if (arg == "-r" || arg == "--repeat") {
            if (i + 1 >= argc || !parseNonNegativeInt(argv[++i], &options.repeat)) {
                std::cerr << "Erro: repeat invalido\n";
                return false;
            }
            continue;
        }
        if (arg == "-w" || arg == "--warmup") {
            if (i + 1 >= argc || !parseNonNegativeInt(argv[++i], &options.warmup)) {
                std::cerr << "Erro: warmup invalido\n";
                return false;
            }
            continue;
        }
        if (arg == "-q" || arg == "--quiet") {
            options.quiet = true;
            continue;
        }
        if (arg == "--json") {
            options.json = true;
            continue;
        }
        if (arg == "--config") {
            if (i + 1 >= argc || !parseConfigFile(argv[++i], options)) {
                return false;
            }
            continue;
        }
        if (arg == "--attribute") {
            if (i + 1 >= argc || !parseAttributeMode(argv[++i], &options.attributeMode)) {
                std::cerr << "Erro: attribute invalido\n";
                return false;
            }
            continue;
        }
        if (arg == "--thresholds") {
            if (i + 1 >= argc || !parseThresholdList(argv[++i], &options.thresholds)) {
                std::cerr << "Erro: thresholds invalido\n";
                return false;
            }
            continue;
        }
        if (arg == "--method") {
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para --method\n";
                return false;
            }
            const std::string method = argv[++i];
            if (!isKnownMethod(method)) {
                std::cerr << "Erro: metodo desconhecido " << method << "\n";
                return false;
            }
            options.methods.push_back(method);
            continue;
        }
        if (arg == "--no-validate") {
            options.noValidate = true;
            continue;
        }
        if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Erro: opcao desconhecida " << arg << "\n";
            return false;
        }
        options.images.push_back(arg);
    }
    return true;
}

static ImageUInt8Ptr loadGrayImage(const std::string &path) {
    int numCols = 0;
    int numRows = 0;
    int numChannels = 0;
    uint8_t *data = stbi_load(path.c_str(), &numCols, &numRows, &numChannels, 1);
    if (data == nullptr) {
        throw std::runtime_error("Failed to load image: " + path);
    }

    ImageUInt8Ptr image = ImageUInt8::create(numRows, numCols);
    std::copy(data, data + (numRows * numCols), image->rawData());
    stbi_image_free(data);
    return image;
}

static Attribute toBoundingBoxAttribute(AttributeMode mode) {
    switch (mode) {
        case AttributeMode::BBoxWidth: return BOX_WIDTH;
        case AttributeMode::BBoxHeight: return BOX_HEIGHT;
        case AttributeMode::BBoxDiagonal: return DIAGONAL_LENGTH;
        case AttributeMode::Area: return AREA;
    }
    return AREA;
}

static std::unique_ptr<DynamicAttributeComputer> makeIncrementalAttributeComputer(DynamicComponentTree *tree, AttributeMode mode) {
    switch (mode) {
        case AttributeMode::Area:
            return std::make_unique<DynamicAreaComputer>(tree);
        case AttributeMode::BBoxWidth:
        case AttributeMode::BBoxHeight:
        case AttributeMode::BBoxDiagonal:
            return std::make_unique<DynamicBoundingBoxComputer>(tree, toBoundingBoxAttribute(mode));
    }
    return nullptr;
}

static std::vector<float> computeAttributeVector(DynamicComponentTree *tree, AttributeMode mode) {
    switch (mode) {
        case AttributeMode::Area: {
            DynamicAreaComputer computer(tree);
            return computer.compute();
        }
        case AttributeMode::BBoxWidth:
        case AttributeMode::BBoxHeight:
        case AttributeMode::BBoxDiagonal: {
            DynamicBoundingBoxComputer computer(tree, toBoundingBoxAttribute(mode));
            return computer.compute();
        }
    }
    return {};
}

static StatsSummary summarize(const std::vector<double> &samples) {
    StatsSummary stats;
    if (samples.empty()) {
        return stats;
    }
    std::vector<double> sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    stats.min = sorted.front();
    stats.max = sorted.back();
    stats.mean = std::accumulate(sorted.begin(), sorted.end(), 0.0) / static_cast<double>(sorted.size());
    if (sorted.size() % 2 == 0) {
        const std::size_t i = sorted.size() / 2;
        stats.median = (sorted[i - 1] + sorted[i]) * 0.5;
    } else {
        stats.median = sorted[sorted.size() / 2];
    }
    double variance = 0.0;
    for (double sample : sorted) {
        const double delta = sample - stats.mean;
        variance += delta * delta;
    }
    variance /= static_cast<double>(sorted.size());
    stats.stddev = std::sqrt(variance);
    return stats;
}

static FlatzonePartition computeFlatzonePartition(ImageUInt8Ptr image, const AdjacencyRelationPtr &adj) {
    FlatzonePartition partition;
    partition.labelByPixel.assign(static_cast<std::size_t>(image->getSize()), -1);

    const auto *data = image->rawData();
    FastQueue<int> queue(std::max(1, image->getSize() / 4));
    int nextLabel = 0;
    for (int p = 0; p < image->getSize(); ++p) {
        if (partition.labelByPixel[static_cast<std::size_t>(p)] >= 0) {
            continue;
        }
        partition.labelByPixel[static_cast<std::size_t>(p)] = nextLabel;
        queue.push(p);
        while (!queue.empty()) {
            const int q = queue.pop();
            for (int n : adj->getNeighborPixels(q)) {
                if (partition.labelByPixel[static_cast<std::size_t>(n)] < 0 && data[n] == data[q]) {
                    partition.labelByPixel[static_cast<std::size_t>(n)] = nextLabel;
                    queue.push(n);
                }
            }
        }
        ++nextLabel;
    }

    partition.count = nextLabel;
    return partition;
}

static long long computePixelBordersFromMask(ImageUInt8Ptr image,
                                             const AdjacencyRelationPtr &adj,
                                             const std::vector<uint8_t> &mask,
                                             const std::vector<PixelId> &pixels) {
    const auto *data = image->rawData();
    long long count = 0;
    for (PixelId p : pixels) {
        bool isBorder = false;
        for (PixelId q : adj->getNeighborPixels(p)) {
            if (mask[static_cast<std::size_t>(q)] == 0 || data[q] != data[p]) {
                isBorder = true;
                break;
            }
        }
        if (isBorder) {
            ++count;
        }
    }
    return count;
}

static ImageMetrics computeImageMetrics(ImageUInt8Ptr image, double radioAdj) {
    ImageMetrics metrics;
    metrics.rows = image->getNumRows();
    metrics.cols = image->getNumCols();

    auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radioAdj);
    const FlatzonePartition partition = computeFlatzonePartition(image, adj);
    DynamicComponentTree maxtree(image, true, adj);
    DynamicComponentTree mintree(image, false, adj);

    std::vector<uint8_t> fullMask(static_cast<std::size_t>(image->getSize()), 1);
    std::vector<PixelId> allPixels(static_cast<std::size_t>(image->getSize()));
    std::iota(allPixels.begin(), allPixels.end(), 0);

    metrics.numFlatzonesImage = partition.count;
    metrics.numPixelBordersImage = computePixelBordersFromMask(image, adj, fullMask, allPixels);
    metrics.numNodesMaxTree = maxtree.getNumNodes();
    metrics.numNodesMinTree = mintree.getNumNodes();
    return metrics;
}

static void collectSupportPixels(DynamicComponentTree *tree,
                                 const std::vector<NodeId> &nodes,
                                 std::vector<uint8_t> &mask,
                                 std::vector<PixelId> &pixels,
                                 int &discardedNodes) {
    pixels.clear();
    std::fill(mask.begin(), mask.end(), 0);
    discardedNodes = 0;
    for (NodeId rootId : nodes) {
        for (NodeId subtreeNodeId : tree->getNodeSubtree(rootId)) {
            ++discardedNodes;
            for (PixelId pixelId : tree->getProperParts(subtreeNodeId)) {
                if (mask[static_cast<std::size_t>(pixelId)] == 0) {
                    mask[static_cast<std::size_t>(pixelId)] = 1;
                    pixels.push_back(pixelId);
                }
            }
        }
    }
}

static ThresholdMetrics collectThresholdMetrics(DynamicComponentTree *tree,
                                                const std::vector<NodeId> &nodes,
                                                ImageUInt8Ptr image,
                                                const AdjacencyRelationPtr &adj,
                                                const FlatzonePartition &partition) {
    ThresholdMetrics metrics;
    if (tree == nullptr || image == nullptr || adj == nullptr) {
        return metrics;
    }

    std::vector<uint8_t> mask(static_cast<std::size_t>(image->getSize()), 0);
    std::vector<PixelId> pixels;
    int discardedNodes = 0;
    collectSupportPixels(tree, nodes, mask, pixels, discardedNodes);

    std::vector<uint8_t> seenFlatzones(static_cast<std::size_t>(std::max(0, partition.count)), 0);
    int flatzones = 0;
    for (PixelId pixelId : pixels) {
        const int label = partition.labelByPixel[static_cast<std::size_t>(pixelId)];
        if (label >= 0 && seenFlatzones[static_cast<std::size_t>(label)] == 0) {
            seenFlatzones[static_cast<std::size_t>(label)] = 1;
            ++flatzones;
        }
    }

    metrics.numDiscardedNodes = discardedNodes;
    metrics.flatzones = flatzones;
    metrics.inputTreeNodes = tree->getNumNodes();
    metrics.numPixelBorders = computePixelBordersFromMask(image, adj, mask, pixels);
    metrics.area = static_cast<long long>(pixels.size());
    metrics.numSubtrees = static_cast<int>(nodes.size());
    metrics.valid = true;
    return metrics;
}

static std::vector<NodeId> getNodesThreshold(DynamicComponentTree *tree, const std::vector<float> &attribute, int threshold) {
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

static void addSubtreeMetrics(DynamicSubtreeMetrics &dst, const DynamicSubtreeMetrics &src) {
    if (!src.valid) {
        return;
    }
    if (!dst.valid) {
        dst = src;
        dst.valid = true;
        return;
    }
    dst.area_subtree += src.area_subtree;
    dst.count_proper_parts += src.count_proper_parts;
    dst.count_nodes_Fb += src.count_nodes_Fb;
    dst.count_total_nodes_merged += src.count_total_nodes_merged;
    dst.count_total_nodes_and_children_merged += src.count_total_nodes_and_children_merged;
    dst.count_total_nodes_removed += src.count_total_nodes_removed;
    dst.count_adjacent_nodes += src.count_adjacent_nodes;
    dst.count_attribute_recomputations += src.count_attribute_recomputations;
    dst.count_unique_nodes_attribute_recomputed += src.count_unique_nodes_attribute_recomputed;
    dst.count_children_scanned_in_attribute_recomputations += src.count_children_scanned_in_attribute_recomputations;
    dst.count_attribute_recomputations_sweep += src.count_attribute_recomputations_sweep;
    dst.count_attribute_recomputations_finalize += src.count_attribute_recomputations_finalize;
    dst.count_attribute_recomputations_absorb += src.count_attribute_recomputations_absorb;
    dst.count_attribute_recomputations_skipped_preserved += src.count_attribute_recomputations_skipped_preserved;
    dst.count_finalize_duplicates_final_union_node += src.count_finalize_duplicates_final_union_node;
    dst.count_finalize_duplicates_node_ca += src.count_finalize_duplicates_node_ca;
    dst.count_finalize_duplicates_node_ca_parent += src.count_finalize_duplicates_node_ca_parent;
    dst.count_finalize_duplicates_candidate_root += src.count_finalize_duplicates_candidate_root;
    dst.count_absorb_duplicates_root += src.count_absorb_duplicates_root;
    dst.count_absorb_duplicates_non_root += src.count_absorb_duplicates_non_root;
    dst.area_tau_star += src.area_tau_star;
    dst.loop_iterations += src.loop_iterations;
    dst.avg_area_nodes_merged += src.avg_area_nodes_merged;
    dst.max_area_nodes_merged = std::max(dst.max_area_nodes_merged, src.max_area_nodes_merged);
}

static DynamicSubtreeMetrics finalizeAggregatedSubtreeMetrics(DynamicSubtreeMetrics metrics, int numContributions) {
    if (metrics.valid && numContributions > 1) {
        metrics.avg_area_nodes_merged /= static_cast<double>(numContributions);
    }
    return metrics;
}

static IterationAggregate aggregateIterationSamples(std::size_t index, int threshold, const std::vector<IterationSample> &samples) {
    IterationAggregate aggregate;
    aggregate.index = static_cast<int>(index) + 1;
    aggregate.threshold = threshold;
    aggregate.count = samples.size();

    std::vector<double> phase1;
    std::vector<double> phase2;
    std::vector<double> total;
    std::vector<double> phase1Update;
    std::vector<double> phase1Pruning;
    std::vector<double> phase2Update;
    std::vector<double> phase2Pruning;
    int phase1AlgCount = 0;
    int phase2AlgCount = 0;

    for (const IterationSample &sample : samples) {
        phase1.push_back(sample.phase1Ms);
        phase2.push_back(sample.phase2Ms);
        total.push_back(sample.totalMs);
        phase1Update.push_back(sample.phase1UpdateMs);
        phase1Pruning.push_back(sample.phase1PruningMs);
        phase2Update.push_back(sample.phase2UpdateMs);
        phase2Pruning.push_back(sample.phase2PruningMs);

        if (!aggregate.phase1InputMetrics.valid && sample.phase1InputMetrics.valid) {
            aggregate.phase1InputMetrics = sample.phase1InputMetrics;
        }
        if (!aggregate.phase2InputMetrics.valid && sample.phase2InputMetrics.valid) {
            aggregate.phase2InputMetrics = sample.phase2InputMetrics;
        }
        if (sample.phase1AlgMetrics.valid) {
            addSubtreeMetrics(aggregate.phase1AlgMetrics, sample.phase1AlgMetrics);
            ++phase1AlgCount;
        }
        if (sample.phase2AlgMetrics.valid) {
            addSubtreeMetrics(aggregate.phase2AlgMetrics, sample.phase2AlgMetrics);
            ++phase2AlgCount;
        }
    }

    aggregate.phase1 = summarize(phase1);
    aggregate.phase2 = summarize(phase2);
    aggregate.total = summarize(total);
    aggregate.phase1Update = summarize(phase1Update);
    aggregate.phase1Pruning = summarize(phase1Pruning);
    aggregate.phase2Update = summarize(phase2Update);
    aggregate.phase2Pruning = summarize(phase2Pruning);
    aggregate.phase1AlgMetrics = finalizeAggregatedSubtreeMetrics(aggregate.phase1AlgMetrics, phase1AlgCount);
    aggregate.phase2AlgMetrics = finalizeAggregatedSubtreeMetrics(aggregate.phase2AlgMetrics, phase2AlgCount);
    return aggregate;
}

static MethodResult runNaiveMethod(ImageUInt8Ptr image, const BenchOptions &options) {
    auto runner = [&](bool collectMetrics, std::vector<IterationSample> *iterationSamples) {
        auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), options.radioAdj);
        ImageUInt8Ptr current = image->clone();
        Stopwatch totalSw;
        totalSw.start();

        for (std::size_t thresholdIndex = 0; thresholdIndex < options.thresholds.size(); ++thresholdIndex) {
            const int threshold = options.thresholds[thresholdIndex];
            IterationSample iter;
            iter.threshold = threshold;
            Stopwatch iterSw;
            iterSw.start();

            Stopwatch phase1Sw;
            phase1Sw.start();
            DynamicComponentTree maxTree(current, true, adj);
            const auto attributeMax = computeAttributeVector(&maxTree, options.attributeMode);
            const auto nodesToPruneMax = getNodesThreshold(&maxTree, attributeMax, threshold);
            phase1Sw.pause();
            iterSw.pause();
            totalSw.pause();
            if (collectMetrics) {
                const FlatzonePartition partition = computeFlatzonePartition(current, adj);
                iter.phase1InputMetrics = collectThresholdMetrics(&maxTree, nodesToPruneMax, current, adj, partition);
            }
            phase1Sw.resume();
            iterSw.resume();
            totalSw.resume();
            Stopwatch phase1PruneSw;
            phase1PruneSw.start();
            for (NodeId nodeId : nodesToPruneMax) {
                if (nodeId != maxTree.getRoot()) {
                    maxTree.pruneNode(nodeId);
                }
            }
            phase1PruneSw.pause();
            current = maxTree.reconstructionImage();
            phase1Sw.pause();
            iter.phase1Ms = elapsedMs(phase1Sw);
            iter.phase1PruningMs = elapsedMs(phase1PruneSw);
            iter.phase1UpdateMs = std::max(0.0, iter.phase1Ms - iter.phase1PruningMs);

            Stopwatch phase2Sw;
            phase2Sw.start();
            DynamicComponentTree minTree(current, false, adj);
            const auto attributeMin = computeAttributeVector(&minTree, options.attributeMode);
            const auto nodesToPruneMin = getNodesThreshold(&minTree, attributeMin, threshold);
            phase2Sw.pause();
            iterSw.pause();
            totalSw.pause();
            if (collectMetrics) {
                const FlatzonePartition partition = computeFlatzonePartition(current, adj);
                iter.phase2InputMetrics = collectThresholdMetrics(&minTree, nodesToPruneMin, current, adj, partition);
            }
            phase2Sw.resume();
            iterSw.resume();
            totalSw.resume();
            Stopwatch phase2PruneSw;
            phase2PruneSw.start();
            for (NodeId nodeId : nodesToPruneMin) {
                if (nodeId != minTree.getRoot()) {
                    minTree.pruneNode(nodeId);
                }
            }
            phase2PruneSw.pause();
            current = minTree.reconstructionImage();
            phase2Sw.pause();
            iter.phase2Ms = elapsedMs(phase2Sw);
            iter.phase2PruningMs = elapsedMs(phase2PruneSw);
            iter.phase2UpdateMs = std::max(0.0, iter.phase2Ms - iter.phase2PruningMs);

            iterSw.pause();
            iter.totalMs = elapsedMs(iterSw);
            if (iterationSamples != nullptr) {
                iterationSamples->push_back(iter);
            }
        }

        totalSw.pause();
        return std::make_pair(current, elapsedMs(totalSw));
    };

    for (int i = 0; i < options.warmup; ++i) {
        runner(false, nullptr);
    }

    MethodResult result;
    result.name = "naive";
    std::vector<std::vector<IterationSample>> allIterationSamples;
    for (int i = 0; i < options.repeat; ++i) {
        std::vector<IterationSample> iterSamples;
        const bool collectMetrics = (i == 0);
        auto [output, totalMs] = runner(collectMetrics, &iterSamples);
        result.output = output;
        result.totalSamplesMs.push_back(totalMs);
        result.buildSamplesMs.push_back(0.0);
        result.buildMaxTreeSamplesMs.push_back(0.0);
        result.buildMinTreeSamplesMs.push_back(0.0);
        allIterationSamples.push_back(std::move(iterSamples));
    }

    result.iterations.reserve(options.thresholds.size());
    for (std::size_t thresholdIndex = 0; thresholdIndex < options.thresholds.size(); ++thresholdIndex) {
        std::vector<IterationSample> samples;
        samples.reserve(allIterationSamples.size());
        for (const auto &runSamples : allIterationSamples) {
            samples.push_back(runSamples[thresholdIndex]);
        }
        result.iterations.push_back(aggregateIterationSamples(thresholdIndex, options.thresholds[thresholdIndex], samples));
    }
    return result;
}

static MethodResult runOurSubtreeMethod(ImageUInt8Ptr image, const BenchOptions &options) {
    auto runner = [&](bool collectMetrics, std::vector<IterationSample> *iterationSamples) {
        Stopwatch totalSw;
        totalSw.start();

        Stopwatch buildSw;
        buildSw.start();
        auto adj = std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), options.radioAdj);
        Stopwatch buildMaxSw;
        buildMaxSw.start();
        DynamicComponentTree maxTree(image, true, adj);
        buildMaxSw.pause();
        Stopwatch buildMinSw;
        buildMinSw.start();
        DynamicComponentTree minTree(image, false, adj);
        buildMinSw.pause();
        buildSw.pause();

        DynamicComponentTreeAdjustmentInstrumented<AltitudeType> adjust(&minTree, &maxTree, *adj);
        auto minAttributeComputer = makeIncrementalAttributeComputer(&minTree, options.attributeMode);
        auto maxAttributeComputer = makeIncrementalAttributeComputer(&maxTree, options.attributeMode);
        std::vector<float> minAreaBuffer = computeAttributeVector(&minTree, options.attributeMode);
        std::vector<float> maxAreaBuffer = computeAttributeVector(&maxTree, options.attributeMode);
        adjust.setAttributeComputer(*minAttributeComputer,
                                    *maxAttributeComputer,
                                    std::span<float>(minAreaBuffer),
                                    std::span<float>(maxAreaBuffer));

        for (std::size_t thresholdIndex = 0; thresholdIndex < options.thresholds.size(); ++thresholdIndex) {
            const int threshold = options.thresholds[thresholdIndex];
            IterationSample iter;
            iter.threshold = threshold;
            Stopwatch iterSw;
            iterSw.start();

            Stopwatch phase1Sw;
            phase1Sw.start();
            const auto attributeMax = computeAttributeVector(&maxTree, options.attributeMode);
            std::vector<float> targetAreaPhase1;
            if (collectMetrics) {
                DynamicAreaComputer targetAreaComputer(&minTree);
                targetAreaPhase1 = targetAreaComputer.compute();
            }
            const auto nodesToPruneMax = getNodesThreshold(&maxTree, attributeMax, threshold);
            phase1Sw.pause();
            iterSw.pause();
            totalSw.pause();
            if (collectMetrics) {
                const ImageUInt8Ptr currentImage = maxTree.reconstructionImage();
                const FlatzonePartition partition = computeFlatzonePartition(currentImage, adj);
                iter.phase1InputMetrics = collectThresholdMetrics(&maxTree, nodesToPruneMax, currentImage, adj, partition);
            }
            phase1Sw.resume();
            iterSw.resume();
            totalSw.resume();

            Stopwatch phase1UpdateSw;
            Stopwatch phase1PruneSw;
            TimerHooksContext hooks{&phase1Sw, &iterSw, &totalSw};
            DynamicSubtreeMetrics phase1MetricsAcc;
            int phase1MetricsCount = 0;
            for (NodeId subtreeRoot : nodesToPruneMax) {
                if (subtreeRoot == InvalidNode || subtreeRoot == maxTree.getRoot() || !maxTree.isAlive(subtreeRoot)) {
                    continue;
                }
                DynamicSubtreeMetrics metrics;
                if (collectMetrics) {
                    adjust.setMetricsAreaBuffer(&minTree, std::span<const float>(targetAreaPhase1));
                    adjust.setMetrics(&metrics, &hooks, pauseTimerHooks, resumeTimerHooks);
                } else {
                    adjust.setMetricsAreaBuffer(nullptr, {});
                    adjust.setMetrics(nullptr);
                }
                phase1UpdateSw.resume();
                adjust.updateTree(&minTree, subtreeRoot);
                phase1UpdateSw.pause();
                if (metrics.valid) {
                    addSubtreeMetrics(phase1MetricsAcc, metrics);
                    ++phase1MetricsCount;
                }
                phase1PruneSw.resume();
                maxTree.pruneNode(subtreeRoot);
                phase1PruneSw.pause();
            }
            adjust.setMetrics(nullptr);
            phase1Sw.pause();
            iter.phase1Ms = elapsedMs(phase1Sw);
            iter.phase1UpdateMs = elapsedUsToMsRounded(phase1UpdateSw);
            iter.phase1PruningMs = elapsedUsToMsRounded(phase1PruneSw);
            iter.phase1AlgMetrics = finalizeAggregatedSubtreeMetrics(phase1MetricsAcc, phase1MetricsCount);

            Stopwatch phase2Sw;
            phase2Sw.start();
            const auto attributeMin = computeAttributeVector(&minTree, options.attributeMode);
            std::vector<float> targetAreaPhase2;
            if (collectMetrics) {
                DynamicAreaComputer targetAreaComputer(&maxTree);
                targetAreaPhase2 = targetAreaComputer.compute();
            }
            const auto nodesToPruneMin = getNodesThreshold(&minTree, attributeMin, threshold);
            phase2Sw.pause();
            iterSw.pause();
            totalSw.pause();
            if (collectMetrics) {
                const ImageUInt8Ptr currentImage = minTree.reconstructionImage();
                const FlatzonePartition partition = computeFlatzonePartition(currentImage, adj);
                iter.phase2InputMetrics = collectThresholdMetrics(&minTree, nodesToPruneMin, currentImage, adj, partition);
            }
            phase2Sw.resume();
            iterSw.resume();
            totalSw.resume();

            Stopwatch phase2UpdateSw;
            Stopwatch phase2PruneSw;
            hooks = TimerHooksContext{&phase2Sw, &iterSw, &totalSw};
            DynamicSubtreeMetrics phase2MetricsAcc;
            int phase2MetricsCount = 0;
            for (NodeId subtreeRoot : nodesToPruneMin) {
                if (subtreeRoot == InvalidNode || subtreeRoot == minTree.getRoot() || !minTree.isAlive(subtreeRoot)) {
                    continue;
                }
                DynamicSubtreeMetrics metrics;
                if (collectMetrics) {
                    adjust.setMetricsAreaBuffer(&maxTree, std::span<const float>(targetAreaPhase2));
                    adjust.setMetrics(&metrics, &hooks, pauseTimerHooks, resumeTimerHooks);
                } else {
                    adjust.setMetricsAreaBuffer(nullptr, {});
                    adjust.setMetrics(nullptr);
                }
                phase2UpdateSw.resume();
                adjust.updateTree(&maxTree, subtreeRoot);
                phase2UpdateSw.pause();
                if (metrics.valid) {
                    addSubtreeMetrics(phase2MetricsAcc, metrics);
                    ++phase2MetricsCount;
                }
                phase2PruneSw.resume();
                minTree.pruneNode(subtreeRoot);
                phase2PruneSw.pause();
            }
            adjust.setMetrics(nullptr);
            phase2Sw.pause();
            iter.phase2Ms = elapsedMs(phase2Sw);
            iter.phase2UpdateMs = elapsedUsToMsRounded(phase2UpdateSw);
            iter.phase2PruningMs = elapsedUsToMsRounded(phase2PruneSw);
            iter.phase2AlgMetrics = finalizeAggregatedSubtreeMetrics(phase2MetricsAcc, phase2MetricsCount);

            iterSw.pause();
            iter.totalMs = elapsedMs(iterSw);
            if (iterationSamples != nullptr) {
                iterationSamples->push_back(iter);
            }
        }

        totalSw.pause();
        return std::make_tuple(minTree.reconstructionImage(),
                               elapsedMs(totalSw),
                               elapsedMs(buildSw),
                               elapsedMs(buildMaxSw),
                               elapsedMs(buildMinSw));
    };

    for (int i = 0; i < options.warmup; ++i) {
        runner(false, nullptr);
    }

    MethodResult result;
    result.name = "our_subtree";
    std::vector<std::vector<IterationSample>> allIterationSamples;
    for (int i = 0; i < options.repeat; ++i) {
        std::vector<IterationSample> iterSamples;
        const bool collectMetrics = (i == 0);
        auto [output, totalMs, buildMs, buildMaxMs, buildMinMs] = runner(collectMetrics, &iterSamples);
        result.output = output;
        result.totalSamplesMs.push_back(totalMs);
        result.buildSamplesMs.push_back(buildMs);
        result.buildMaxTreeSamplesMs.push_back(buildMaxMs);
        result.buildMinTreeSamplesMs.push_back(buildMinMs);
        allIterationSamples.push_back(std::move(iterSamples));
    }

    result.iterations.reserve(options.thresholds.size());
    for (std::size_t thresholdIndex = 0; thresholdIndex < options.thresholds.size(); ++thresholdIndex) {
        std::vector<IterationSample> samples;
        samples.reserve(allIterationSamples.size());
        for (const auto &runSamples : allIterationSamples) {
            samples.push_back(runSamples[thresholdIndex]);
        }
        result.iterations.push_back(aggregateIterationSamples(thresholdIndex, options.thresholds[thresholdIndex], samples));
    }
    return result;
}

static std::string attributeModeToString(AttributeMode mode) {
    switch (mode) {
        case AttributeMode::Area: return "area";
        case AttributeMode::BBoxWidth: return "bbox_width";
        case AttributeMode::BBoxHeight: return "bbox_height";
        case AttributeMode::BBoxDiagonal: return "bbox_diagonal";
    }
    return "area";
}

static std::string jsonEscape(std::string_view text) {
    std::string out;
    for (char c : text) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            default: out += c; break;
        }
    }
    return out;
}

static void printJsonStatsObject(const StatsSummary &stats, std::size_t n) {
    std::cout << "{"
              << "\"mean\":" << stats.mean << ","
              << "\"median\":" << stats.median << ","
              << "\"stddev\":" << stats.stddev << ","
              << "\"min\":" << stats.min << ","
              << "\"max\":" << stats.max << ","
              << "\"n\":" << n
              << "}";
}

static void printJsonInputMetricsObject(const ThresholdMetrics &metrics) {
    std::cout << "{"
              << "\"num_discarded_nodes\":" << metrics.numDiscardedNodes << ","
              << "\"flatzones\":" << metrics.flatzones << ","
              << "\"input_tree_nodes\":" << metrics.inputTreeNodes << ","
              << "\"num_pixel_borders\":" << metrics.numPixelBorders << ","
              << "\"area\":" << metrics.area << ","
              << "\"num_subtrees\":" << metrics.numSubtrees
              << "}";
}

static void printJsonAlgMetricsObject(const DynamicSubtreeMetrics &metrics) {
    std::cout << "{"
              << "\"a\":" << metrics.a << ","
              << "\"b\":" << metrics.b << ","
              << "\"area_substree\":" << metrics.area_subtree << ","
              << "\"count_proper_parts\":" << metrics.count_proper_parts << ","
              << "\"count_nodes_Fb\":" << metrics.count_nodes_Fb << ","
              << "\"count_total_nodes_merged\":" << metrics.count_total_nodes_merged << ","
              << "\"count_total_nodes_and_children_merged\":" << metrics.count_total_nodes_and_children_merged << ","
              << "\"count_total_nodes_removed\":" << metrics.count_total_nodes_removed << ","
              << "\"count_adjacent_nodes\":" << metrics.count_adjacent_nodes << ","
              << "\"count_attribute_recomputations\":" << metrics.count_attribute_recomputations << ","
              << "\"count_unique_nodes_attribute_recomputed\":" << metrics.count_unique_nodes_attribute_recomputed << ","
              << "\"count_children_scanned_in_attribute_recomputations\":" << metrics.count_children_scanned_in_attribute_recomputations << ","
              << "\"count_attribute_recomputations_sweep\":" << metrics.count_attribute_recomputations_sweep << ","
              << "\"count_attribute_recomputations_finalize\":" << metrics.count_attribute_recomputations_finalize << ","
              << "\"count_attribute_recomputations_absorb\":" << metrics.count_attribute_recomputations_absorb << ","
              << "\"count_attribute_recomputations_skipped_preserved\":" << metrics.count_attribute_recomputations_skipped_preserved << ","
              << "\"count_finalize_duplicates_final_union_node\":" << metrics.count_finalize_duplicates_final_union_node << ","
              << "\"count_finalize_duplicates_node_ca\":" << metrics.count_finalize_duplicates_node_ca << ","
              << "\"count_finalize_duplicates_node_ca_parent\":" << metrics.count_finalize_duplicates_node_ca_parent << ","
              << "\"count_finalize_duplicates_candidate_root\":" << metrics.count_finalize_duplicates_candidate_root << ","
              << "\"count_absorb_duplicates_root\":" << metrics.count_absorb_duplicates_root << ","
              << "\"count_absorb_duplicates_non_root\":" << metrics.count_absorb_duplicates_non_root << ","
              << "\"area_tau_star\":" << metrics.area_tau_star << ","
              << "\"loop_iterations\":" << metrics.loop_iterations << ","
              << "\"avg_area_nodes_merged\":" << metrics.avg_area_nodes_merged << ","
              << "\"max_area_nodes_merged\":" << metrics.max_area_nodes_merged
              << "}";
}

static void printJsonIntArray(const std::vector<int> &values) {
    std::cout << "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            std::cout << ",";
        }
        std::cout << values[i];
    }
    std::cout << "]";
}

static void printJsonBenchmarkSummary(const std::vector<MethodBenchmarkResult> &methods,
                                      const BenchOptions &options) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "{\n  \"methods\":[\n";
    for (std::size_t methodIndex = 0; methodIndex < methods.size(); ++methodIndex) {
        const auto &method = methods[methodIndex];
        if (methodIndex > 0) {
            std::cout << ",\n";
        }
        std::cout << "    {\n"
                  << "      \"method\":\"" << jsonEscape(method.name) << "\",\n"
                  << "      \"attribute\":\"" << attributeModeToString(options.attributeMode) << "\",\n"
                  << "      \"thresholds\":";
        printJsonIntArray(options.thresholds);
        std::cout << ",\n"
                  << "      \"radio_adj\":" << options.radioAdj << ",\n"
                  << "      \"images\":[\n";

        for (std::size_t imageIndex = 0; imageIndex < method.images.size(); ++imageIndex) {
            const auto &image = method.images[imageIndex];
            const StatsSummary totalStats = summarize(image.result.totalSamplesMs);
            const StatsSummary buildStats = summarize(image.result.buildSamplesMs);
            const StatsSummary buildMaxStats = summarize(image.result.buildMaxTreeSamplesMs);
            const StatsSummary buildMinStats = summarize(image.result.buildMinTreeSamplesMs);
            if (imageIndex > 0) {
                std::cout << ",\n";
            }
            std::cout << "        {\n"
                      << "          \"image\":\"" << jsonEscape(image.imagePath) << "\",\n"
                      << "          \"rows\":" << image.imageMetrics.rows << ",\n"
                      << "          \"cols\":" << image.imageMetrics.cols << ",\n"
                      << "          \"num_flatzones_image\":" << image.imageMetrics.numFlatzonesImage << ",\n"
                      << "          \"num_pixel_borders_image\":" << image.imageMetrics.numPixelBordersImage << ",\n"
                      << "          \"num_nodes_maxtree\":" << image.imageMetrics.numNodesMaxTree << ",\n"
                      << "          \"num_nodes_mintree\":" << image.imageMetrics.numNodesMinTree << ",\n"
                      << "          \"total_ms\":";
            printJsonStatsObject(totalStats, image.result.totalSamplesMs.size());
            std::cout << ",\n"
                      << "          \"build_ms\":";
            printJsonStatsObject(buildStats, image.result.buildSamplesMs.size());
            std::cout << ",\n"
                      << "          \"build_maxtree_ms\":";
            printJsonStatsObject(buildMaxStats, image.result.buildMaxTreeSamplesMs.size());
            std::cout << ",\n"
                      << "          \"build_mintree_ms\":";
            printJsonStatsObject(buildMinStats, image.result.buildMinTreeSamplesMs.size());
            std::cout << ",\n"
                      << "          \"iterations\":[\n";

            for (std::size_t iterIndex = 0; iterIndex < image.result.iterations.size(); ++iterIndex) {
                const auto &iter = image.result.iterations[iterIndex];
                if (iterIndex > 0) {
                    std::cout << ",\n";
                }
                std::cout << "            {\n"
                          << "              \"index\":" << iter.index << ",\n"
                          << "              \"threshold\":" << iter.threshold << ",\n"
                          << "              \"phase1\":";
                printJsonStatsObject(iter.phase1, iter.count);
                std::cout << ",\n"
                          << "              \"phase1_input_metrics\":";
                printJsonInputMetricsObject(iter.phase1InputMetrics);
                std::cout << ",\n"
                          << "              \"phase1_alg_metrics\":";
                printJsonAlgMetricsObject(iter.phase1AlgMetrics);
                std::cout << ",\n"
                          << "              \"phase1_update_ms\":";
                printJsonStatsObject(iter.phase1Update, iter.count);
                std::cout << ",\n"
                          << "              \"phase1_pruning_ms\":";
                printJsonStatsObject(iter.phase1Pruning, iter.count);
                std::cout << ",\n"
                          << "              \"phase2\":";
                printJsonStatsObject(iter.phase2, iter.count);
                std::cout << ",\n"
                          << "              \"phase2_input_metrics\":";
                printJsonInputMetricsObject(iter.phase2InputMetrics);
                std::cout << ",\n"
                          << "              \"phase2_alg_metrics\":";
                printJsonAlgMetricsObject(iter.phase2AlgMetrics);
                std::cout << ",\n"
                          << "              \"phase2_update_ms\":";
                printJsonStatsObject(iter.phase2Update, iter.count);
                std::cout << ",\n"
                          << "              \"phase2_pruning_ms\":";
                printJsonStatsObject(iter.phase2Pruning, iter.count);
                std::cout << ",\n"
                          << "              \"total\":";
                printJsonStatsObject(iter.total, iter.count);
                std::cout << "\n            }";
            }

            std::cout << "\n          ]\n"
                      << "        }";
        }

        std::cout << "\n      ]\n"
                  << "    }";
    }
    std::cout << "\n  ]\n}\n";
}

static void printConsoleSummary(const MethodImageResult &entry) {
    const StatsSummary totalStats = summarize(entry.result.totalSamplesMs);
    std::cout << std::fixed << std::setprecision(2)
              << "Image: " << entry.imagePath << "\n"
              << "Method: " << entry.result.name << "\n"
              << "Total mean: " << totalStats.mean << " ms"
              << ", median: " << totalStats.median << " ms"
              << ", stddev: " << totalStats.stddev << " ms"
              << ", min: " << totalStats.min << " ms"
              << ", max: " << totalStats.max << " ms"
              << ", n=" << entry.result.totalSamplesMs.size() << "\n";
}

int main(int argc, char **argv) {
    try {
        BenchOptions options;
        if (!parseCommandLine(argc, argv, options)) {
            return 1;
        }

        if (options.images.empty()) {
            printUsage(argv[0]);
            return 1;
        }
        if (options.thresholds.empty()) {
            options.thresholds = {50, 100, 150, 200, 250, 300, 350, 400};
        }
        if (options.methods.empty()) {
            options.methods = {"naive", "our_subtree"};
        }

        for (const auto &method : options.methods) {
            if (!isKnownMethod(method)) {
                throw std::runtime_error("Unknown method: " + method);
            }
        }

        std::vector<MethodBenchmarkResult> benchmarkResults;
        benchmarkResults.reserve(options.methods.size());
        for (const auto &method : options.methods) {
            benchmarkResults.push_back(MethodBenchmarkResult{method, {}});
        }

        const std::size_t numImages = options.images.size();
        for (std::size_t imageIndex = 0; imageIndex < numImages; ++imageIndex) {
            const auto &imagePath = options.images[imageIndex];
            std::cerr << "["
                      << (imageIndex + 1) << "/" << numImages
                      << "] processing image: " << imagePath << "\n";

            ImageUInt8Ptr image = loadGrayImage(imagePath);
            const ImageMetrics imageMetrics = computeImageMetrics(image, options.radioAdj);

            MethodResult naiveResult;
            bool hasNaive = false;
            for (std::size_t methodIndex = 0; methodIndex < options.methods.size(); ++methodIndex) {
                const auto &method = options.methods[methodIndex];
                MethodResult result;
                if (method == "naive") {
                    result = runNaiveMethod(image, options);
                    naiveResult = result;
                    hasNaive = true;
                } else if (method == "our_subtree") {
                    result = runOurSubtreeMethod(image, options);
                    if (!options.noValidate && hasNaive && !naiveResult.output->isEqual(result.output)) {
                        throw std::runtime_error("Validation failed for our_subtree on image: " + imagePath);
                    }
                }

                benchmarkResults[methodIndex].images.push_back(MethodImageResult{imagePath, imageMetrics, std::move(result)});
                if (!options.json && !options.quiet) {
                    printConsoleSummary(benchmarkResults[methodIndex].images.back());
                }
            }
        }

        if (options.json) {
            printJsonBenchmarkSummary(benchmarkResults, options);
        }

        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
