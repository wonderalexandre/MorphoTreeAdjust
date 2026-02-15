
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentBySubtree.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/FlatZonesGraph.hpp"

#include "../morphoTreeAdjust/include/AttributeComputer.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"
#include "../tests/Tests.hpp"

static bool gEnableStats = true;
static bool gCollectMetrics = false;

static inline long long elapsedMs(const Stopwatch& sw) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(sw.elapsed()).count();
}

static inline long long elapsedUs(const Stopwatch& sw) {
    return std::chrono::duration_cast<std::chrono::microseconds>(sw.elapsed()).count();
}

enum class AttributeMode {
    Area,
    BBoxWidth,
    BBoxHeight,
    BBoxDiagonal
};

enum class GraphType {
    Any,
    OnDemand,
    Eager,
    OnPixel
};

struct BenchOptions {
    std::vector<std::string> images;
    std::vector<std::string> methods;
    std::vector<int> thresholds;
    AttributeMode attributeMode = AttributeMode::Area;
    GraphType graphType = GraphType::Any;
    int repeat = 1;
    int warmup = 0;
    double radioAdj = 1.5;
    size_t cutoffHybrid = 2;
    size_t cutoffHybridSub = 1;
    bool quiet = false;
    bool json = false;
    bool noValidate = false;
    bool validateOnly = false;
    bool stats = false;
    bool iterStats = false;
    bool showHelp = false;
};

struct ImageEntry {
    std::string name;
    ImageUInt8Ptr img;
    int rows = 0;
    int cols = 0;
};

struct StatsSummary {
    double mean = 0.0;
    double median = 0.0;
    double stddev = 0.0;
    double min = 0.0;
    double max = 0.0;
};

struct ThresholdMetrics {
    int numNodes = 0;
    int numFlatZones = 0;
    int inputTreeNodes = 0;
    int numSubtrees = 0;
    long long numPixelBorders = 0;
    long long area = 0;
    bool valid = false;
};


struct BuildTimes {
    long long total = 0;
    long long graph = 0;
    long long maxtree = 0;
    long long mintree = 0;
    bool valid = false;
};

struct IterationTimes {
    int threshold = 0;
    long long phase1 = 0;
    long long phase2 = 0;
    long long total = 0;
    long long phase1Update = 0;
    long long phase1Pruning = 0;
    long long phase2Update = 0;
    long long phase2Pruning = 0;
    ThresholdMetrics phase1Metrics;
    ThresholdMetrics phase2Metrics;
    SubtreeMetrics phase1SubtreeMetrics;
    SubtreeMetrics phase2SubtreeMetrics;
};

struct IterationAggregate {
    int threshold = 0;
    StatsSummary phase1;
    StatsSummary phase2;
    StatsSummary total;
    StatsSummary phase1Update;
    StatsSummary phase1Pruning;
    StatsSummary phase2Update;
    StatsSummary phase2Pruning;
    size_t count = 0;
    ThresholdMetrics phase1Metrics;
    ThresholdMetrics phase2Metrics;
    SubtreeMetrics phase1SubtreeMetrics;
    SubtreeMetrics phase2SubtreeMetrics;
};

struct MethodResult {
    std::string name;
    std::vector<long long> samples;
    std::vector<long long> buildSamples;
    std::vector<long long> buildGraphSamples;
    std::vector<long long> buildMaxTreeSamples;
    std::vector<long long> buildMinTreeSamples;
    std::vector<IterationAggregate> iterStats;
    std::string phase1Label;
    std::string phase2Label;
    ImageUInt8Ptr output;
};

struct ImageMetrics {
    int numFlatzones = 0;
    long long numPixelBorders = 0;
    int maxTreeNodes = 0;
    int minTreeNodes = 0;
    bool valid = false;
};

struct ImageBenchmarkResult {
    std::string image;
    int rows = 0;
    int cols = 0;
    ImageMetrics metrics;
    MethodResult result;
};

struct MethodBenchmarkResult {
    std::string name;
    std::vector<ImageBenchmarkResult> images;
};

struct ValidationResult {
    std::string image;
    bool eqLeaf = false;
    bool eqHybrid = false;
    bool eqSubtree = false;
    bool eqHybridSub = false;
};

struct NullBuffer : std::streambuf {
    int overflow(int c) override { return c; }
};

struct CoutSilencer {
    std::streambuf* old = nullptr;
    NullBuffer nullBuffer;
    explicit CoutSilencer(bool active) {
        if (active) {
            old = std::cout.rdbuf(&nullBuffer);
        }
    }
    ~CoutSilencer() {
        if (old) {
            std::cout.rdbuf(old);
        }
    }
};

static void printUsage(const char* argv0) {
    std::cout << "Usage: " << argv0 << " [options] <image>...\n"
              << "Options:\n"
              << "  -i, --image <path>   Add image path (can repeat)\n"
              << "  -r, --repeat <N>     Repetitions per image (default 1)\n"
              << "  -w, --warmup <N>     Warmup runs per image (default 0)\n"
              << "  -q, --quiet          Reduce console output\n"
              << "  --json               Emit JSON lines output\n"
              << "  --config <path>      Load options from file\n"
              << "  --attribute <name>   Attribute: area|bbox_width|bbox_height|bbox_diagonal (bbox=diagonal)\n"
              << "  --method <name>      Run only selected method (can repeat)\n"
              << "  --graph-type <name>  Graph type: any|on_demand|eager|on_pixel (naive/pixel alias)\n"
              << "  --no-validate        Skip output validation\n"
              << "  --validate-only      Only validate outputs (no benchmark)\n"
              << "  --stats              Enable metrics/logs per threshold\n"
              << "  --iter-stats         Print per-iteration time statistics\n"
              << "  -h, --help           Show this help\n";
}

static bool parseNonNegativeInt(const char* text, int* out) {
    if (!text || !*text) return false;
    char* end = nullptr;
    long value = std::strtol(text, &end, 10);
    if (!end || *end != '\0') return false;
    if (value < 0 || value > std::numeric_limits<int>::max()) return false;
    *out = static_cast<int>(value);
    return true;
}

static bool containsString(const std::vector<std::string>& values, std::string_view target) {
    for (const auto& value : values) {
        if (value == target) return true;
    }
    return false;
}

static bool isKnownMethod(std::string_view name) {
    static const std::array<std::string_view, 13> kMethods = {
        "naive",
        "our",
        "our_hybrid",
        "our_subtree",
        "our_hybrid_subtree",
        "our_eager",
        "our_hybrid_eager",
        "our_subtree_eager",
        "our_hybrid_subtree_eager",
        "our_naive",
        "our_hybrid_naive",
        "our_subtree_naive",
        "our_hybrid_subtree_naive"
    };
    for (std::string_view method : kMethods) {
        if (method == name) return true;
    }
    return false;
}

static std::string toLowerCopy(std::string_view text);
static std::string trimCopy(std::string_view text);

static bool parseAttributeMode(std::string_view text, AttributeMode* out) {
    const std::string value = toLowerCopy(trimCopy(text));
    if (value == "area" || value == "areacomputer") {
        *out = AttributeMode::Area;
        return true;
    }
    if (value == "bbox_width" || value == "box_width" || value == "boundingbox_width") {
        *out = AttributeMode::BBoxWidth;
        return true;
    }
    if (value == "bbox_height" || value == "box_height" || value == "boundingbox_height") {
        *out = AttributeMode::BBoxHeight;
        return true;
    }
    if (value == "bbox" || value == "boundingbox" || value == "boundingboxcomputer" ||
        value == "bbox_diag" || value == "bbox_diagonal" || value == "box_diagonal") {
        *out = AttributeMode::BBoxDiagonal;
        return true;
    }
    return false;
}

static const char* attributeModeToString(AttributeMode mode) {
    switch (mode) {
        case AttributeMode::Area: return "area";
        case AttributeMode::BBoxWidth: return "bbox_width";
        case AttributeMode::BBoxHeight: return "bbox_height";
        case AttributeMode::BBoxDiagonal: return "bbox_diagonal";
    }
    return "area";
}

static Attribute toBoundingBoxAttribute(AttributeMode mode) {
    switch (mode) {
        case AttributeMode::BBoxWidth: return BOX_WIDTH;
        case AttributeMode::BBoxHeight: return BOX_HEIGHT;
        case AttributeMode::BBoxDiagonal: return DIAGONAL_LENGTH;
        case AttributeMode::Area: break;
    }
    return DIAGONAL_LENGTH;
}

static const char* graphTypeForMethod(std::string_view method) {
    if (method == "naive") {
        return "on_pixel";
    }
    if (method.find("_eager") != std::string_view::npos) {
        return "eager";
    }
    if (method.find("_naive") != std::string_view::npos) {
        return "on_pixel";
    }
    if (method.rfind("our", 0) == 0) {
        return "on_demand";
    }
    return "unknown";
}

static std::string methodNameWithGraphSuffix(std::string_view method) {
    if (method == "naive") {
        return "naive_on_pixel";
    }
    if (method.ends_with("_naive")) {
        return std::string(method.substr(0, method.size() - 6)) + "_on_pixel";
    }
    if (method.ends_with("_eager")) {
        return std::string(method);
    }
    if (method.rfind("our", 0) == 0) {
        return std::string(method) + "_on_demand";
    }
    return std::string(method);
}

static GraphType graphTypeFromMethod(std::string_view method) {
    if (method == "naive") {
        return GraphType::OnPixel;
    }
    if (method.find("_eager") != std::string_view::npos) {
        return GraphType::Eager;
    }
    if (method.find("_naive") != std::string_view::npos) {
        return GraphType::OnPixel;
    }
    if (method.rfind("our", 0) == 0) {
        return GraphType::OnDemand;
    }
    return GraphType::Any;
}

static const char* graphTypeToString(GraphType type) {
    switch (type) {
        case GraphType::Any: return "any";
        case GraphType::OnDemand: return "on_demand";
        case GraphType::Eager: return "eager";
        case GraphType::OnPixel: return "on_pixel";
    }
    return "any";
}

static bool parseGraphType(std::string_view text, GraphType* out) {
    const std::string value = toLowerCopy(trimCopy(text));
    if (value == "any") {
        *out = GraphType::Any;
        return true;
    }
    if (value == "on_demand" || value == "ondemand" || value == "default") {
        *out = GraphType::OnDemand;
        return true;
    }
    if (value == "eager") {
        *out = GraphType::Eager;
        return true;
    }
    if (value == "on_pixel" || value == "pixel" || value == "naive") {
        *out = GraphType::OnPixel;
        return true;
    }
    return false;
}

static bool isBaseOurMethod(std::string_view method) {
    return method == "our" ||
           method == "our_hybrid" ||
           method == "our_subtree" ||
           method == "our_hybrid_subtree";
}

static bool hasGraphSuffix(std::string_view method) {
    return method.find("_eager") != std::string_view::npos ||
           method.find("_naive") != std::string_view::npos;
}

static std::string applyGraphSuffix(std::string_view method, GraphType type) {
    if (!isBaseOurMethod(method)) {
        return std::string(method);
    }
    switch (type) {
        case GraphType::Eager:
            return std::string(method) + "_eager";
        case GraphType::OnPixel:
            return std::string(method) + "_naive";
        case GraphType::OnDemand:
        case GraphType::Any:
            return std::string(method);
    }
    return std::string(method);
}

static std::vector<std::string> normalizeMethodsForGraphType(const std::vector<std::string>& methods,
                                                             GraphType graphType) {
    if (graphType == GraphType::Any) {
        return methods;
    }
    std::vector<std::string> normalized;
    normalized.reserve(methods.size());
    for (const auto& method : methods) {
        std::string resolved = method;
        if (isBaseOurMethod(method) && !hasGraphSuffix(method)) {
            resolved = applyGraphSuffix(method, graphType);
        }
        if (graphTypeFromMethod(resolved) != graphType) {
            continue;
        }
        if (!containsString(normalized, resolved)) {
            normalized.push_back(resolved);
        }
    }
    return normalized;
}

template <typename CNPsType, typename GraphT>
static std::vector<float> computeAttributeVector(ComponentTree<CNPsType, GraphT>* tree,
                                                 AttributeMode mode) {
    switch (mode) {
        case AttributeMode::Area: {
            AreaComputer<CNPsType, GraphT> computer(tree);
            return computer.compute();
        }
        case AttributeMode::BBoxWidth:
        case AttributeMode::BBoxHeight:
        case AttributeMode::BBoxDiagonal: {
            BoundingBoxComputer<CNPsType, GraphT> computer(tree, toBoundingBoxAttribute(mode));
            return computer.compute();
        }
    }
    AreaComputer<CNPsType, GraphT> computer(tree);
    return computer.compute();
}

static bool globMatch(std::string_view pattern, std::string_view text) {
    size_t p = 0;
    size_t t = 0;
    size_t star = std::string::npos;
    size_t match = 0;
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

static bool expandGlobPath(const std::string& pattern, std::vector<std::string>& out) {
    if (!hasGlobChars(pattern)) {
        out.push_back(pattern);
        return true;
    }
    std::filesystem::path path(pattern);
    std::filesystem::path dir = path.parent_path();
    if (dir.empty()) {
        dir = ".";
    }
    const std::string filenamePattern = path.filename().string();
    std::error_code ec;
    if (!std::filesystem::exists(dir, ec)) {
        return false;
    }
    std::vector<std::string> matches;
    for (const auto& entry : std::filesystem::directory_iterator(dir, ec)) {
        if (ec) break;
        const auto name = entry.path().filename().string();
        if (globMatch(filenamePattern, name)) {
            matches.push_back(entry.path().string());
        }
    }
    if (matches.empty()) return false;
    std::sort(matches.begin(), matches.end());
    out.insert(out.end(), matches.begin(), matches.end());
    return true;
}

static bool addImageFromValue(std::string_view value, std::vector<std::string>& out) {
    std::string trimmed = trimCopy(value);
    if (trimmed.empty()) return false;
    std::vector<std::string> expanded;
    if (!expandGlobPath(trimmed, expanded)) {
        return false;
    }
    out.insert(out.end(), expanded.begin(), expanded.end());
    return true;
}

static std::string toLowerCopy(std::string_view text) {
    std::string out;
    out.reserve(text.size());
    for (char c : text) {
        if (c >= 'A' && c <= 'Z') {
            out.push_back(static_cast<char>(c - 'A' + 'a'));
        } else {
            out.push_back(c);
        }
    }
    return out;
}

static std::string trimCopy(std::string_view text) {
    size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start])) != 0) {
        ++start;
    }
    size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1])) != 0) {
        --end;
    }
    return std::string(text.substr(start, end - start));
}

static bool parseBool(std::string_view text, bool* out) {
    const std::string value = toLowerCopy(trimCopy(text));
    if (value == "1" || value == "true" || value == "yes" || value == "on") {
        *out = true;
        return true;
    }
    if (value == "0" || value == "false" || value == "no" || value == "off") {
        *out = false;
        return true;
    }
    return false;
}

static bool parseDouble(std::string_view text, double* out) {
    std::string value = trimCopy(text);
    if (value.empty()) return false;
    char* end = nullptr;
    double parsed = std::strtod(value.c_str(), &end);
    if (!end || *end != '\0') return false;
    *out = parsed;
    return true;
}

static bool parseIntList(std::string_view text, std::vector<int>& out) {
    std::string value = trimCopy(text);
    if (value.empty()) return false;
    for (char& c : value) {
        if (c == ',') c = ' ';
    }
    std::istringstream iss(value);
    std::vector<int> parsed;
    int v = 0;
    while (iss >> v) {
        parsed.push_back(v);
    }
    if (parsed.empty()) return false;
    out.swap(parsed);
    return true;
}

static void appendStringList(std::string_view text, std::vector<std::string>& out) {
    std::string value = trimCopy(text);
    if (value.empty()) return;
    for (char& c : value) {
        if (c == ',') c = ' ';
    }
    std::istringstream iss(value);
    std::string token;
    while (iss >> token) {
        out.push_back(token);
    }
}

static bool parseConfigFile(const std::string& path, BenchOptions& options) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Erro: nao foi possivel abrir arquivo de config " << path << "\n";
        return false;
    }
    std::string line;
    int lineNo = 0;
    while (std::getline(file, line)) {
        ++lineNo;
        const size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line.erase(commentPos);
        }
        std::string trimmed = trimCopy(line);
        if (trimmed.empty()) continue;
        size_t sepPos = trimmed.find('=');
        if (sepPos == std::string::npos) {
            sepPos = trimmed.find(':');
        }
        if (sepPos == std::string::npos) {
            std::cerr << "Erro: linha invalida no config (" << lineNo << ")\n";
            return false;
        }
        const std::string key = toLowerCopy(trimCopy(trimmed.substr(0, sepPos)));
        const std::string value = trimCopy(trimmed.substr(sepPos + 1));
        if (key == "image") {
            if (!addImageFromValue(value, options.images)) {
                std::cerr << "Erro: image invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "images") {
            std::vector<std::string> tokens;
            appendStringList(value, tokens);
            if (tokens.empty()) {
                std::cerr << "Erro: images invalido no config (" << lineNo << ")\n";
                return false;
            }
            for (const auto& token : tokens) {
                if (!addImageFromValue(token, options.images)) {
                    std::cerr << "Erro: images invalido no config (" << lineNo << ")\n";
                    return false;
                }
            }
        } else if (key == "method") {
            if (!isKnownMethod(value)) {
                std::cerr << "Erro: metodo desconhecido no config (" << lineNo << ")\n";
                return false;
            }
            if (!containsString(options.methods, value)) {
                options.methods.push_back(value);
            }
        } else if (key == "methods") {
            std::vector<std::string> values;
            appendStringList(value, values);
            for (const auto& method : values) {
                if (!isKnownMethod(method)) {
                    std::cerr << "Erro: metodo desconhecido no config (" << lineNo << ")\n";
                    return false;
                }
                if (!containsString(options.methods, method)) {
                    options.methods.push_back(method);
                }
            }
        } else if (key == "thresholds") {
            if (!parseIntList(value, options.thresholds)) {
                std::cerr << "Erro: thresholds invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "attribute" || key == "attribute_type") {
            if (!parseAttributeMode(value, &options.attributeMode)) {
                std::cerr << "Erro: attribute invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "graph_type" || key == "graph") {
            if (!parseGraphType(value, &options.graphType)) {
                std::cerr << "Erro: graph_type invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "repeat") {
            if (!parseNonNegativeInt(value.c_str(), &options.repeat) || options.repeat < 1) {
                std::cerr << "Erro: repeat invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "warmup") {
            if (!parseNonNegativeInt(value.c_str(), &options.warmup)) {
                std::cerr << "Erro: warmup invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "radioadj" || key == "radio_adj" || key == "adj_radius") {
            if (!parseDouble(value, &options.radioAdj)) {
                std::cerr << "Erro: radioAdj invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "cutoffhybrid" || key == "cutoff_hybrid") {
            int parsed = 0;
            if (!parseNonNegativeInt(value.c_str(), &parsed)) {
                std::cerr << "Erro: cutoffHybrid invalido no config (" << lineNo << ")\n";
                return false;
            }
            options.cutoffHybrid = static_cast<size_t>(parsed);
        } else if (key == "cutoffhybridsubtree" || key == "cutoff_hybrid_subtree") {
            int parsed = 0;
            if (!parseNonNegativeInt(value.c_str(), &parsed)) {
                std::cerr << "Erro: cutoffHybridSubtree invalido no config (" << lineNo << ")\n";
                return false;
            }
            options.cutoffHybridSub = static_cast<size_t>(parsed);
        } else if (key == "quiet") {
            if (!parseBool(value, &options.quiet)) {
                std::cerr << "Erro: quiet invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "json") {
            if (!parseBool(value, &options.json)) {
                std::cerr << "Erro: json invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "no_validate") {
            if (!parseBool(value, &options.noValidate)) {
                std::cerr << "Erro: no_validate invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "validate_only") {
            if (!parseBool(value, &options.validateOnly)) {
                std::cerr << "Erro: validate_only invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "stats") {
            if (!parseBool(value, &options.stats)) {
                std::cerr << "Erro: stats invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else if (key == "iter_stats") {
            if (!parseBool(value, &options.iterStats)) {
                std::cerr << "Erro: iter_stats invalido no config (" << lineNo << ")\n";
                return false;
            }
        } else {
            std::cerr << "Erro: chave desconhecida no config (" << lineNo << ")\n";
            return false;
        }
    }
    return true;
}

static ImageEntry loadImageEntry(const std::string& filename) {
    ImageEntry entry;
    entry.name = filename;
    int numCols = 0;
    int numRows = 0;
    int nchannels = 0;
    uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    (void)nchannels;
    if (!data) {
        std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
        entry.img = getCameramanImage();
        entry.cols = entry.img->getNumCols();
        entry.rows = entry.img->getNumRows();
    } else {
        entry.img = ImageUInt8::create(numRows, numCols);
        std::copy(data, data + (numRows * numCols), entry.img->rawData());
        stbi_image_free(data);
        entry.cols = numCols;
        entry.rows = numRows;
    }
    return entry;
}

static bool parseArgs(int argc, char* argv[], BenchOptions& options) {
    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            options.showHelp = true;
            return true;
        }
        if (arg == "--") {
            for (++i; i < argc; ++i) {
                options.images.emplace_back(argv[i]);
            }
            break;
        }
        if (arg == "-i" || arg == "--image") {
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para " << arg << "\n";
                return false;
            }
            options.images.emplace_back(argv[++i]);
            continue;
        }
        if (arg == "-r" || arg == "--repeat") {
            if (i + 1 >= argc || !parseNonNegativeInt(argv[++i], &options.repeat) || options.repeat < 1) {
                std::cerr << "Erro: valor invalido para " << arg << "\n";
                return false;
            }
            continue;
        }
        if (arg == "-w" || arg == "--warmup") {
            if (i + 1 >= argc || !parseNonNegativeInt(argv[++i], &options.warmup)) {
                std::cerr << "Erro: valor invalido para " << arg << "\n";
                return false;
            }
            continue;
        }
        if (arg == "--attribute") {
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para " << arg << "\n";
                return false;
            }
            if (!parseAttributeMode(argv[++i], &options.attributeMode)) {
                std::cerr << "Erro: atributo invalido para " << arg << "\n";
                return false;
            }
            continue;
        }
        if (arg == "--graph-type") {
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para " << arg << "\n";
                return false;
            }
            if (!parseGraphType(argv[++i], &options.graphType)) {
                std::cerr << "Erro: graph_type invalido para " << arg << "\n";
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
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para " << arg << "\n";
                return false;
            }
            if (!parseConfigFile(std::string(argv[++i]), options)) {
                return false;
            }
            continue;
        }
        if (arg == "--method") {
            if (i + 1 >= argc) {
                std::cerr << "Erro: falta argumento para " << arg << "\n";
                return false;
            }
            std::string method(argv[++i]);
            if (!isKnownMethod(method)) {
                std::cerr << "Erro: metodo desconhecido " << method << "\n";
                return false;
            }
            if (!containsString(options.methods, method)) {
                options.methods.push_back(method);
            }
            continue;
        }
        if (arg == "--no-validate") {
            options.noValidate = true;
            continue;
        }
        if (arg == "--validate-only") {
            options.validateOnly = true;
            continue;
        }
        if (arg == "--stats") {
            options.stats = true;
            continue;
        }
        if (arg == "--iter-stats") {
            options.iterStats = true;
            continue;
        }
        if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Erro: opcao desconhecida " << arg << "\n";
            return false;
        }
        options.images.emplace_back(arg);
    }
    return true;
}

static StatsSummary computeStats(std::vector<long long> samples) {
    StatsSummary stats;
    if (samples.empty()) return stats;
    long long minValue = samples[0];
    long long maxValue = samples[0];
    long long sum = 0;
    for (long long v : samples) {
        sum += v;
        minValue = std::min(minValue, v);
        maxValue = std::max(maxValue, v);
    }
    stats.min = static_cast<double>(minValue);
    stats.max = static_cast<double>(maxValue);
    stats.mean = static_cast<double>(sum) / static_cast<double>(samples.size());

    std::sort(samples.begin(), samples.end());
    if (samples.size() % 2 == 0) {
        stats.median = (static_cast<double>(samples[samples.size() / 2 - 1]) +
                        static_cast<double>(samples[samples.size() / 2])) /
                       2.0;
    } else {
        stats.median = static_cast<double>(samples[samples.size() / 2]);
    }

    double variance = 0.0;
    for (long long v : samples) {
        double diff = static_cast<double>(v) - stats.mean;
        variance += diff * diff;
    }
    variance /= static_cast<double>(samples.size());
    stats.stddev = std::sqrt(variance);
    return stats;
}

static std::string jsonEscape(std::string_view value) {
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        switch (c) {
            case '\\': out += "\\\\"; break;
            case '"': out += "\\\""; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default: out += c; break;
        }
    }
    return out;
}

static void printIndent(int level) {
    for (int i = 0; i < level; ++i) {
        std::cout << "  ";
    }
}

static void printJsonStatsObject(const StatsSummary& stats, size_t count);

static void printJsonRows(const std::string& image,
                          const MethodResult& result,
                          int rows,
                          int cols) {
    StatsSummary stats = computeStats(result.samples);
    StatsSummary buildStats = computeStats(result.buildSamples);
    StatsSummary buildGraphStats = computeStats(result.buildGraphSamples);
    StatsSummary buildMaxTreeStats = computeStats(result.buildMaxTreeSamples);
    StatsSummary buildMinTreeStats = computeStats(result.buildMinTreeSamples);
    std::cout << std::fixed << std::setprecision(2);
    const std::string imageEscaped = jsonEscape(image);
    for (size_t i = 0; i < result.samples.size(); ++i) {
        std::cout << "{\n";
        printIndent(1);
        std::cout << "\"image\":\"" << imageEscaped << "\",\n";
        printIndent(1);
        if (result.name == "naive") {
            std::cout << "\"method\":\"naive\",\n";
            printIndent(1);
            std::cout << "\"graph_type\":\"\",\n";
        } else {
            std::cout << "\"method\":\"" << methodNameWithGraphSuffix(result.name) << "\",\n";
            printIndent(1);
            std::cout << "\"graph_type\":\"" << graphTypeForMethod(result.name) << "\",\n";
        }
        printIndent(1);
        std::cout << "\"repeat\":" << (i + 1) << ",\n";
        printIndent(1);
        std::cout << "\"ms\":" << result.samples[i] << ",\n";
        printIndent(1);
        std::cout << "\"mean\":" << stats.mean << ",\n";
        printIndent(1);
        std::cout << "\"median\":" << stats.median << ",\n";
        printIndent(1);
        std::cout << "\"stddev\":" << stats.stddev << ",\n";
        printIndent(1);
        std::cout << "\"min\":" << stats.min << ",\n";
        printIndent(1);
        std::cout << "\"max\":" << stats.max << ",\n";
        printIndent(1);
        std::cout << "\"build_ms\":";
        printJsonStatsObject(buildStats, result.buildSamples.size());
        std::cout << ",\n";
        printIndent(1);
        std::cout << "\"build_graph_ms\":";
        printJsonStatsObject(buildGraphStats, result.buildGraphSamples.size());
        std::cout << ",\n";
        printIndent(1);
        std::cout << "\"build_maxtree_ms\":";
        printJsonStatsObject(buildMaxTreeStats, result.buildMaxTreeSamples.size());
        std::cout << ",\n";
        printIndent(1);
        std::cout << "\"build_mintree_ms\":";
        printJsonStatsObject(buildMinTreeStats, result.buildMinTreeSamples.size());
        std::cout << ",\n";
        printIndent(1);
        std::cout << "\"rows\":" << rows << ",\n";
        printIndent(1);
        std::cout << "\"cols\":" << cols << "\n";
        std::cout << "}\n";
    }
}

static void printJsonValidation(const std::string& image,
                                bool eqLeaf,
                                bool eqHybrid,
                                bool eqSubtree,
                                bool eqHybridSub) {
    std::cout << "{\n";
    printIndent(1);
    std::cout << "\"image\":\"" << jsonEscape(image) << "\",\n";
    printIndent(1);
    std::cout << "\"validation\":{\n";
    printIndent(2);
    std::cout << "\"naive_vs_our\":" << (eqLeaf ? "true" : "false") << ",\n";
    printIndent(2);
    std::cout << "\"naive_vs_our_hybrid\":" << (eqHybrid ? "true" : "false") << ",\n";
    printIndent(2);
    std::cout << "\"naive_vs_our_subtree\":" << (eqSubtree ? "true" : "false") << ",\n";
    printIndent(2);
    std::cout << "\"naive_vs_our_hybrid_subtree\":" << (eqHybridSub ? "true" : "false") << "\n";
    printIndent(1);
    std::cout << "}\n";
    std::cout << "}\n";
}

static void printJsonStatsObject(const StatsSummary& stats, size_t count) {
    std::cout << "{"
              << "\"mean\":" << stats.mean << ","
              << "\"median\":" << stats.median << ","
              << "\"stddev\":" << stats.stddev << ","
              << "\"min\":" << stats.min << ","
              << "\"max\":" << stats.max << ","
              << "\"n\":" << count
              << "}";
}

static void printJsonMetricsObject(const ThresholdMetrics& metrics) {
    std::cout << "{"
              << "\"num_discarded_nodes\":" << metrics.numNodes << ","
              << "\"flatzones\":" << metrics.numFlatZones << ","
              << "\"input_tree_nodes\":" << metrics.inputTreeNodes << ","
              << "\"num_pixel_borders\":" << metrics.numPixelBorders << ","
              << "\"area\":" << metrics.area << ","
              << "\"num_subtrees\":" << metrics.numSubtrees
              << "}";
}

static void printJsonSubtreeMetricsObject(const SubtreeMetrics& metrics) {
    std::cout << "{"
              << "\"a\":" << metrics.a << ","
              << "\"b\":" << metrics.b << ","
              << "\"area_substree\":" << metrics.area_substree << ","
              << "\"count_proper_parts\":" << metrics.count_proper_parts << ","
              << "\"count_nodes_Fb\":" << metrics.count_nodes_Fb << ","
              << "\"count_total_nodes_merged\":" << metrics.count_total_nodes_merged << ","
              << "\"count_total_nodes_and_children_merged\":" << metrics.count_total_nodes_and_children_merged << ","
              << "\"count_total_nodes_removed\":" << metrics.count_total_nodes_removed << ","
              << "\"count_adjacent_nodes\":" << metrics.count_adjacent_nodes << ","
              << "\"area_tau_star\":" << metrics.area_tau_star << ","
              << "\"loop_iterations\":" << metrics.loop_iterations << ","
              << "\"avg_area_nodes_merged\":" << metrics.avg_area_nodes_merged << ","
              << "\"max_area_nodes_merged\":" << metrics.max_area_nodes_merged
              << "}";
}

static void printJsonIntArray(const std::vector<int>& values) {
    std::cout << "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << values[i];
    }
    std::cout << "]";
}

static void printJsonBenchmarkSummary(const std::vector<MethodBenchmarkResult>& methods,
                                      AttributeMode attributeMode,
                                      const std::vector<int>& thresholds,
                                      double radioAdj,
                                      const std::vector<ValidationResult>& validations) {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "{\n";
    printIndent(1);
    std::cout << "\"methods\":[\n";
    for (size_t m = 0; m < methods.size(); ++m) {
        const auto& method = methods[m];
        if (m > 0) std::cout << ",\n";
        printIndent(2);
        std::cout << "{\n";
        printIndent(3);
        if (method.name == "naive") {
            std::cout << "\"method\":\"naive\",\n";
            printIndent(3);
            std::cout << "\"graph_type\":\"\",\n";
        } else {
            std::cout << "\"method\":\"" << methodNameWithGraphSuffix(method.name) << "\",\n";
            printIndent(3);
            std::cout << "\"graph_type\":\"" << graphTypeForMethod(method.name) << "\",\n";
        }
        printIndent(3);
        std::cout << "\"attribute\":\"" << attributeModeToString(attributeMode) << "\",\n";
        printIndent(3);
        std::cout << "\"thresholds\":";
        printJsonIntArray(thresholds);
        std::cout << ",\n";
        printIndent(3);
        std::cout << "\"radio_adj\":" << radioAdj << ",\n";
        printIndent(3);
        std::cout << "\"images\":[\n";
        for (size_t i = 0; i < method.images.size(); ++i) {
            const auto& img = method.images[i];
            if (i > 0) std::cout << ",\n";
            StatsSummary totalStats = computeStats(img.result.samples);
            StatsSummary buildStats = computeStats(img.result.buildSamples);
            StatsSummary buildGraphStats = computeStats(img.result.buildGraphSamples);
            StatsSummary buildMaxTreeStats = computeStats(img.result.buildMaxTreeSamples);
            StatsSummary buildMinTreeStats = computeStats(img.result.buildMinTreeSamples);
            printIndent(4);
            std::cout << "{\n";
            printIndent(5);
            std::cout << "\"image\":\"" << jsonEscape(img.image) << "\",\n";
            printIndent(5);
            std::cout << "\"rows\":" << img.rows << ",\n";
            printIndent(5);
            std::cout << "\"cols\":" << img.cols << ",\n";
            printIndent(5);
            std::cout << "\"num_flatzones_image\":" << img.metrics.numFlatzones << ",\n";
            printIndent(5);
            std::cout << "\"num_pixel_borders_image\":" << img.metrics.numPixelBorders << ",\n";
            printIndent(5);
            std::cout << "\"num_nodes_maxtree\":" << img.metrics.maxTreeNodes << ",\n";
            printIndent(5);
            std::cout << "\"num_nodes_mintree\":" << img.metrics.minTreeNodes << ",\n";
            printIndent(5);
            std::cout << "\"total_ms\":";
            printJsonStatsObject(totalStats, img.result.samples.size());
            std::cout << ",\n";
            printIndent(5);
            std::cout << "\"build_ms\":";
            printJsonStatsObject(buildStats, img.result.buildSamples.size());
            std::cout << ",\n";
            printIndent(5);
            std::cout << "\"build_graph_ms\":";
            printJsonStatsObject(buildGraphStats, img.result.buildGraphSamples.size());
            std::cout << ",\n";
            printIndent(5);
            std::cout << "\"build_maxtree_ms\":";
            printJsonStatsObject(buildMaxTreeStats, img.result.buildMaxTreeSamples.size());
            std::cout << ",\n";
            printIndent(5);
            std::cout << "\"build_mintree_ms\":";
            printJsonStatsObject(buildMinTreeStats, img.result.buildMinTreeSamples.size());
            std::cout << ",\n";
            printIndent(5);
            std::cout << "\"iterations\":[\n";
            for (size_t it = 0; it < img.result.iterStats.size(); ++it) {
                const auto& iter = img.result.iterStats[it];
                if (it > 0) std::cout << ",\n";
                printIndent(6);
                std::cout << "{\n";
                printIndent(7);
                std::cout << "\"index\":" << (it + 1) << ",\n";
                printIndent(7);
                std::cout << "\"threshold\":" << iter.threshold << ",\n";
                printIndent(7);
                std::cout << "\"phase1\":";
                printJsonStatsObject(iter.phase1, iter.count);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase1_input_metrics\":";
                printJsonMetricsObject(iter.phase1Metrics);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase1_alg_metrics\":";
                printJsonSubtreeMetricsObject(iter.phase1SubtreeMetrics);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase1_update_ms\":";
                printJsonStatsObject(iter.phase1Update, iter.count);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase1_pruning_ms\":";
                printJsonStatsObject(iter.phase1Pruning, iter.count);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase2\":";
                printJsonStatsObject(iter.phase2, iter.count);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase2_input_metrics\":";
                printJsonMetricsObject(iter.phase2Metrics);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase2_alg_metrics\":";
                printJsonSubtreeMetricsObject(iter.phase2SubtreeMetrics);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase2_update_ms\":";
                printJsonStatsObject(iter.phase2Update, iter.count);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"phase2_pruning_ms\":";
                printJsonStatsObject(iter.phase2Pruning, iter.count);
                std::cout << ",\n";
                printIndent(7);
                std::cout << "\"total\":";
                printJsonStatsObject(iter.total, iter.count);
                std::cout << "\n";
                printIndent(6);
                std::cout << "}";
            }
            std::cout << "\n";
            printIndent(5);
            std::cout << "]\n";
            printIndent(4);
            std::cout << "}";
        }
        if (!method.images.empty()) {
            std::cout << "\n";
        }
        printIndent(3);
        std::cout << "]\n";
        printIndent(2);
        std::cout << "}";
    }
    if (!methods.empty()) {
        std::cout << "\n";
    }
    printIndent(1);
    std::cout << "]";
    if (!validations.empty()) {
        std::cout << ",\n";
        printIndent(1);
        std::cout << "\"validation\":[\n";
        for (size_t i = 0; i < validations.size(); ++i) {
            const auto& v = validations[i];
            if (i > 0) std::cout << ",\n";
            printIndent(2);
            std::cout << "{\n";
            printIndent(3);
            std::cout << "\"image\":\"" << jsonEscape(v.image) << "\",\n";
            printIndent(3);
            std::cout << "\"naive_vs_our\":" << (v.eqLeaf ? "true" : "false") << ",\n";
            printIndent(3);
            std::cout << "\"naive_vs_our_hybrid\":" << (v.eqHybrid ? "true" : "false") << ",\n";
            printIndent(3);
            std::cout << "\"naive_vs_our_subtree\":" << (v.eqSubtree ? "true" : "false") << ",\n";
            printIndent(3);
            std::cout << "\"naive_vs_our_hybrid_subtree\":" << (v.eqHybridSub ? "true" : "false") << "\n";
            printIndent(2);
            std::cout << "}";
        }
        std::cout << "\n";
        printIndent(1);
        std::cout << "]";
    }
    std::cout << "\n}\n";
}

static void printStatsSummary(const MethodResult& result) {
    StatsSummary stats = computeStats(result.samples);
    std::cout << std::fixed << std::setprecision(2)
              << "Summary " << result.name
              << ": mean=" << stats.mean << " ms"
              << ", median=" << stats.median << " ms"
              << ", stddev=" << stats.stddev << " ms"
              << ", min=" << stats.min << " ms"
              << ", max=" << stats.max << " ms"
              << ", n=" << result.samples.size()
              << std::endl;
}

static void printStatsLine(const std::string& label, const StatsSummary& stats) {
    if (label.empty()) return;
    std::cout << label
              << "mean=" << stats.mean << " ms"
              << ", median=" << stats.median << " ms"
              << ", stddev=" << stats.stddev << " ms"
              << ", min=" << stats.min << " ms"
              << ", max=" << stats.max << " ms"
              << "\n";
}

static void printIterationStats(const MethodResult& result) {
    if (result.iterStats.empty()) return;
    std::cout << "Method: " << result.name << "\n";
    std::cout << std::fixed << std::setprecision(2);
    for (size_t i = 0; i < result.iterStats.size(); ++i) {
        const auto& iter = result.iterStats[i];
        std::cout << "Opening/Closing: " << (i + 1)
                  << " \t\tthreshold:" << iter.threshold << "\n";
        printStatsLine(result.phase1Label, iter.phase1);
        printStatsLine(result.phase2Label, iter.phase2);
        printStatsLine("Time: ", iter.total);
        std::cout << "\n";
    }
}

template <typename Func>
static MethodResult runBenchmarkMethod(const std::string& name,
                                       int warmup,
                                       int repeat,
                                       bool quiet,
                                       bool collectIterStats,
                                       const char* phase1Label,
                                       const char* phase2Label,
                                       Func&& func) {
    MethodResult result;
    result.name = name;
    if (phase1Label) {
        result.phase1Label = phase1Label;
    }
    if (phase2Label) {
        result.phase2Label = phase2Label;
    }
    result.samples.reserve(static_cast<size_t>(repeat));
    result.buildSamples.reserve(static_cast<size_t>(repeat));
    result.buildGraphSamples.reserve(static_cast<size_t>(repeat));
    result.buildMaxTreeSamples.reserve(static_cast<size_t>(repeat));
    result.buildMinTreeSamples.reserve(static_cast<size_t>(repeat));
    std::vector<std::vector<long long>> phase1Samples;
    std::vector<std::vector<long long>> phase2Samples;
    std::vector<std::vector<long long>> totalSamples;
    std::vector<std::vector<long long>> phase1UpdateSamples;
    std::vector<std::vector<long long>> phase1PruningSamples;
    std::vector<std::vector<long long>> phase2UpdateSamples;
    std::vector<std::vector<long long>> phase2PruningSamples;
    std::vector<ThresholdMetrics> phase1Metrics;
    std::vector<ThresholdMetrics> phase2Metrics;
    std::vector<SubtreeMetrics> phase1SubtreeMetrics;
    std::vector<SubtreeMetrics> phase2SubtreeMetrics;
    std::vector<int> thresholds;
    {
        CoutSilencer silencer(quiet);
        for (int i = 0; i < warmup; ++i) {
            func(nullptr, nullptr, nullptr);
        }
        for (int i = 0; i < repeat; ++i) {
            std::vector<IterationTimes> iterTimes;
            Stopwatch sw;
            BuildTimes buildTimes;
            sw.start();
            ImageUInt8Ptr out = func(&sw, collectIterStats ? &iterTimes : nullptr, &buildTimes);
            sw.pause();
            result.samples.push_back(elapsedMs(sw));
            if (buildTimes.valid) {
                result.buildSamples.push_back(buildTimes.total);
                result.buildGraphSamples.push_back(buildTimes.graph);
                result.buildMaxTreeSamples.push_back(buildTimes.maxtree);
                result.buildMinTreeSamples.push_back(buildTimes.mintree);
            }
            if (i == 0) {
                result.output = out;
            }
            if (collectIterStats) {
                if (thresholds.empty()) {
                    thresholds.reserve(iterTimes.size());
                    phase1Samples.resize(iterTimes.size());
                    phase2Samples.resize(iterTimes.size());
                    totalSamples.resize(iterTimes.size());
                    phase1UpdateSamples.resize(iterTimes.size());
                    phase1PruningSamples.resize(iterTimes.size());
                    phase2UpdateSamples.resize(iterTimes.size());
                    phase2PruningSamples.resize(iterTimes.size());
                    phase1Metrics.resize(iterTimes.size());
                    phase2Metrics.resize(iterTimes.size());
                    phase1SubtreeMetrics.resize(iterTimes.size());
                    phase2SubtreeMetrics.resize(iterTimes.size());
                    for (const auto& iter : iterTimes) {
                        thresholds.push_back(iter.threshold);
                    }
                }
                const size_t count = std::min(iterTimes.size(), thresholds.size());
                for (size_t idx = 0; idx < count; ++idx) {
                    phase1Samples[idx].push_back(iterTimes[idx].phase1);
                    phase2Samples[idx].push_back(iterTimes[idx].phase2);
                    totalSamples[idx].push_back(iterTimes[idx].total);
                    phase1UpdateSamples[idx].push_back(iterTimes[idx].phase1Update);
                    phase1PruningSamples[idx].push_back(iterTimes[idx].phase1Pruning);
                    phase2UpdateSamples[idx].push_back(iterTimes[idx].phase2Update);
                    phase2PruningSamples[idx].push_back(iterTimes[idx].phase2Pruning);
                    if (!phase1Metrics[idx].valid && iterTimes[idx].phase1Metrics.valid) {
                        phase1Metrics[idx] = iterTimes[idx].phase1Metrics;
                    }
                    if (!phase2Metrics[idx].valid && iterTimes[idx].phase2Metrics.valid) {
                        phase2Metrics[idx] = iterTimes[idx].phase2Metrics;
                    }
                    if (!phase1SubtreeMetrics[idx].valid && iterTimes[idx].phase1SubtreeMetrics.valid) {
                        phase1SubtreeMetrics[idx] = iterTimes[idx].phase1SubtreeMetrics;
                    }
                    if (!phase2SubtreeMetrics[idx].valid && iterTimes[idx].phase2SubtreeMetrics.valid) {
                        phase2SubtreeMetrics[idx] = iterTimes[idx].phase2SubtreeMetrics;
                    }
                }
            }
        }
    }
    if (collectIterStats && !thresholds.empty()) {
        result.iterStats.resize(thresholds.size());
        for (size_t idx = 0; idx < thresholds.size(); ++idx) {
            result.iterStats[idx].threshold = thresholds[idx];
            result.iterStats[idx].phase1 = computeStats(phase1Samples[idx]);
            result.iterStats[idx].phase2 = computeStats(phase2Samples[idx]);
            result.iterStats[idx].total = computeStats(totalSamples[idx]);
            result.iterStats[idx].phase1Update = computeStats(phase1UpdateSamples[idx]);
            result.iterStats[idx].phase1Pruning = computeStats(phase1PruningSamples[idx]);
            result.iterStats[idx].phase2Update = computeStats(phase2UpdateSamples[idx]);
            result.iterStats[idx].phase2Pruning = computeStats(phase2PruningSamples[idx]);
            result.iterStats[idx].count = phase1Samples[idx].size();
            result.iterStats[idx].phase1Metrics = phase1Metrics[idx];
            result.iterStats[idx].phase2Metrics = phase2Metrics[idx];
            result.iterStats[idx].phase1SubtreeMetrics = phase1SubtreeMetrics[idx];
            result.iterStats[idx].phase2SubtreeMetrics = phase2SubtreeMetrics[idx];
        }
    }
    return result;
}


template <typename CNPsType, typename GraphT>
static inline long long computeNumPixelBorders(ComponentTree<CNPsType, GraphT>* tree,
                                               ImageUInt8Ptr img,
                                               const std::vector<NodeId>& nodes) {
    if (!gEnableStats && !gCollectMetrics) return 0;
    if (!tree || !img) return -1;
    AdjacencyRelationPtr adj = tree->getAdjacencyRelation();
    if (!adj) return -1;

    const int numPixels = img->getSize();
    std::vector<uint8_t> mask(static_cast<size_t>(numPixels), 0);
    std::vector<int> uniquePixels;
    uniquePixels.reserve(static_cast<size_t>(numPixels));
    for (NodeId id : nodes) {
        for (int p : tree->getPixelsOfCCById(id)) {
            if (!mask[p]) {
                mask[p] = 1;
                uniquePixels.push_back(p);
            }
        }
    }

    const auto* data = img->rawData();
    long long count = 0;
    for (int p : uniquePixels) {
        bool isBorder = false;
        for (int q : adj->getNeighborPixels(p)) {
            if (!mask[q] || data[q] != data[p]) {
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

static ImageMetrics computeImageMetricsFZ(const DefaultFlatZonesGraph* graph,
                                         ComponentTreeFZ<>* maxtree,
                                         ComponentTreeFZ<>* mintree,
                                         ImageUInt8Ptr img) {
    ImageMetrics metrics;
    if (!graph || !maxtree || !mintree || !img) return metrics;
    metrics.maxTreeNodes = maxtree->getNumNodes();
    metrics.minTreeNodes = mintree->getNumNodes();
    metrics.numFlatzones = graph->getNumFlatZones();
    metrics.numPixelBorders = computeNumPixelBorders(maxtree, img, {maxtree->getRootById()});
    metrics.valid = true;
    return metrics;
}

static ImageMetrics computeImageMetricsPixel(ComponentTreeP* maxtree,
                                            ComponentTreeP* mintree,
                                            ImageUInt8Ptr img) {
    ImageMetrics metrics;
    if (!maxtree || !mintree || !img) return metrics;
    metrics.maxTreeNodes = maxtree->getNumNodes();
    metrics.minTreeNodes = mintree->getNumNodes();
    metrics.numFlatzones = -1;
    metrics.numPixelBorders = computeNumPixelBorders(maxtree, img, {maxtree->getRootById()});
    metrics.valid = true;
    return metrics;
}

class BenchmarkTimers {
public:
    class PauseScope {
    public:
        explicit PauseScope(BenchmarkTimers* timers)
            : timers_(timers),
              wasRunning_{false, false, false} {
            if (timers_) {
                timers_->pauseAll(wasRunning_);
            }
        }

        ~PauseScope() {
            if (timers_) {
                timers_->resumeAll(wasRunning_);
            }
        }

    private:
        BenchmarkTimers* timers_;
        std::array<bool, 3> wasRunning_;
    };

    explicit BenchmarkTimers(Stopwatch* external = nullptr,
                             std::vector<IterationTimes>* iterTimes = nullptr)
        : external_(external),
          iterTimes_(iterTimes) {}

    void setExternal(Stopwatch* external) { external_ = external; }
    void setIterationStorage(std::vector<IterationTimes>* iterTimes) { iterTimes_ = iterTimes; }

    void startIteration() {
        phase_.start();
        total_.start();
    }

    void startPhase() { phase_.start(); }

    void printIterationHeader(size_t index, int threshold) {
        setIterationInfo(index, threshold);
        PauseScope pause(this);
        std::cout << "Opening/Closing: " << (index + 1) << " \t\tthreshold:" << threshold << std::endl;
    }

    void printPhase(const char* label) {
        phase_.pause();
        storePhase1();
        {
            PauseScope pause(this);
            std::cout << label << elapsedMs(phase_) << " ms\n";
        }
        phase_.start();
    }

    void printPhaseAndTotal(const char* label) {
        phase_.pause();
        total_.pause();
        storePhase2AndTotal();
        {
            PauseScope pause(this);
            std::cout << label << elapsedMs(phase_) << " ms\n";
            std::cout << "Time: " << elapsedMs(total_) << " ms\n\n";
        }
    }

    void printBuild(const char* label) {
        phase_.pause();
        if (!label) return;
        {
            PauseScope pause(this);
            std::cout << label << elapsedMs(phase_) << " ms\n";
        }
    }

    void pauseMetrics() {
        if (metricsPaused_) return;
        pauseAll(metricsWasRunning_);
        metricsPaused_ = true;
    }

    void resumeMetrics() {
        if (!metricsPaused_) return;
        resumeAll(metricsWasRunning_);
        metricsPaused_ = false;
    }

private:
    void pauseAll(std::array<bool, 3>& wasRunning) {
        wasRunning[0] = phase_.running();
        wasRunning[1] = total_.running();
        wasRunning[2] = external_ && external_->running();
        if (wasRunning[0]) {
            phase_.pause();
        }
        if (wasRunning[1]) {
            total_.pause();
        }
        if (wasRunning[2]) {
            external_->pause();
        }
    }

    void resumeAll(const std::array<bool, 3>& wasRunning) {
        if (wasRunning[0]) {
            phase_.resume();
        }
        if (wasRunning[1]) {
            total_.resume();
        }
        if (wasRunning[2]) {
            external_->resume();
        }
    }

    void setIterationInfo(size_t index, int threshold) {
        currentIteration_ = index;
        if (!iterTimes_) return;
        if (iterTimes_->size() <= index) {
            iterTimes_->resize(index + 1);
        }
        auto& entry = (*iterTimes_)[index];
        entry.threshold = threshold;
        entry.phase1 = 0;
        entry.phase2 = 0;
        entry.total = 0;
        entry.phase1Update = 0;
        entry.phase1Pruning = 0;
        entry.phase2Update = 0;
        entry.phase2Pruning = 0;
        entry.phase1Metrics = ThresholdMetrics{};
        entry.phase2Metrics = ThresholdMetrics{};
        entry.phase1SubtreeMetrics.reset();
        entry.phase2SubtreeMetrics.reset();
    }

    IterationTimes* currentIteration() {
        if (!iterTimes_) return nullptr;
        if (currentIteration_ >= iterTimes_->size()) return nullptr;
        return &(*iterTimes_)[currentIteration_];
    }

    void storePhase1() {
        if (auto* iter = currentIteration()) {
            iter->phase1 = elapsedMs(phase_);
        }
    }

    void storePhase2AndTotal() {
        if (auto* iter = currentIteration()) {
            iter->phase2 = elapsedMs(phase_);
            iter->total = elapsedMs(total_);
        }
    }

    Stopwatch phase_;
    Stopwatch total_;
    Stopwatch* external_;
    std::vector<IterationTimes>* iterTimes_ = nullptr;
    size_t currentIteration_ = 0;
    std::array<bool, 3> metricsWasRunning_{false, false, false};
    bool metricsPaused_ = false;
};


template <typename CNPsType, typename GraphT>
std::vector<NodeId> getNodesThreshold(ComponentTree<CNPsType, GraphT>* tree,
                                      std::span<float> attribute,
                                      int threshold,
                                      BenchmarkTimers* timers = nullptr,
                                      ImageUInt8Ptr img = nullptr,
                                      ThresholdMetrics* metrics = nullptr) {
    std::vector<NodeId> lista;
    FastQueue<NodeId> queue;
    queue.push(tree->getRootById());

    int numFlatZones = 0;
    int numNodes = 0;
    long long areaSum = 0;
    while (!queue.empty()) {
        NodeId id = queue.pop();
        if (attribute[id] > threshold) {
            for(NodeId c: tree->getChildrenById(id)) { 
                queue.push(c);
            }
        } else {
            BenchmarkTimers::PauseScope pause(timers);
            numFlatZones += (tree->computerNumFlatzoneDescendants(id) + tree->getNumFlatzoneById(id));
            numNodes     += tree->computerNumDescendants(id) + 1;
            areaSum      += tree->getAreaById(id);
            lista.push_back(id);
        }
    }
    const bool collectMetrics = gEnableStats || gCollectMetrics || metrics != nullptr;
    if (collectMetrics) {
        BenchmarkTimers::PauseScope pause(timers);
        long long numPixelBorders = computeNumPixelBorders(tree, img, lista);
        if (metrics) {
            metrics->numNodes = numNodes;
            metrics->numFlatZones = numFlatZones;
            metrics->inputTreeNodes = tree->getNumNodes();
            metrics->numSubtrees = static_cast<int>(lista.size());
            metrics->numPixelBorders = numPixelBorders;
            metrics->area = areaSum;
            metrics->valid = true;
        }
        if (gEnableStats) {
            std::cout << "	Threshold: " << threshold
                    << ", #Nodes: " << numNodes
                    << ", #FlatZones: " << numFlatZones
                    << ", #InputTreeNodes: " << tree->getNumNodes()
                    << ", #numPixelBorders: " << numPixelBorders
                    << ", #Area: " << areaSum
                    << ", #Subtrees: " << lista.size()
                    << std::endl;
        }
    }
    return lista;
}

static ImageUInt8Ptr runPrunePass(ImageUInt8Ptr imgIn,
                                  AdjacencyRelationPtr adj,
                                  bool isMaxTree,
                                  int threshold,
                                  AttributeMode attributeMode,
                                  BenchmarkTimers* timers,
                                  ThresholdMetrics* metrics) {
    ComponentTreePPtr treePtr = std::make_shared<ComponentTreeP>(imgIn, isMaxTree, adj);
    ComponentTreeP* tree = treePtr.get();
    std::vector<float> attribute = computeAttributeVector(tree, attributeMode);

    for (NodeId node : getNodesThreshold(tree, attribute, threshold, timers, imgIn, metrics)) {
        tree->prunning(node);
    }

    return tree->reconstructionImage();
}

template <typename Adjuster>
class BenchAdjuster : public Adjuster {
public:
    using Adjuster::Adjuster;
    using Adjuster::updateTree;
    using Adjuster::prunning;
};

static void pauseBenchmarkTimers(void* ctx) {
    if (!ctx) return;
    static_cast<BenchmarkTimers*>(ctx)->pauseMetrics();
}

static void resumeBenchmarkTimers(void* ctx) {
    if (!ctx) return;
    static_cast<BenchmarkTimers*>(ctx)->resumeMetrics();
}

static void addSubtreeMetrics(SubtreeMetrics& dst, const SubtreeMetrics& src) {
    if (!src.valid) return;
    if (!dst.valid) {
        dst = src;
        dst.valid = true;
        return;
    }
    dst.area_substree += src.area_substree;
    dst.count_proper_parts += src.count_proper_parts;
    dst.count_nodes_Fb += src.count_nodes_Fb;
    dst.count_total_nodes_merged += src.count_total_nodes_merged;
    dst.count_total_nodes_and_children_merged += src.count_total_nodes_and_children_merged;
    dst.count_total_nodes_removed += src.count_total_nodes_removed;
    dst.count_adjacent_nodes += src.count_adjacent_nodes;
    dst.area_tau_star += src.area_tau_star;
    dst.loop_iterations += src.loop_iterations;
}

template <typename AdjusterType, typename GraphT>
static void adjustMinTreeTimed(AdjusterType& adjust,
                               ComponentTreeFZ<GraphT>* mintree,
                               ComponentTreeFZ<GraphT>* maxtree,
                               std::vector<NodeId>& nodesToPruning,
                               IterationTimes* iter,
                               BenchmarkTimers* timers) {
    (void)timers;
    long long updateUs = 0;
    long long pruningUs = 0;
    Stopwatch sw;
    for (NodeId rSubtreeId : nodesToPruning) {
        assert(rSubtreeId != maxtree->getRootById() && "rSubtree is root");
        if (rSubtreeId == InvalidNode) {
            continue;
        }
        sw.start();
        adjust.updateTree(mintree, rSubtreeId);
        sw.pause();
        updateUs += elapsedUs(sw);

        sw.start();
        adjust.prunning(maxtree, rSubtreeId);
        sw.pause();
        pruningUs += elapsedUs(sw);
    }
    if (iter) {
        // Accumulate in microseconds and convert once to avoid per-call truncation to 0 ms.
        iter->phase1Update += (updateUs + 500) / 1000;
        iter->phase1Pruning += (pruningUs + 500) / 1000;
    }
}

template <typename Computer, typename GraphT>
static void adjustMinTreeTimed(BenchAdjuster<ComponentTreeAdjustmentByLeaf<Computer, GraphT>>& adjust,
                               ComponentTreeFZ<GraphT>* mintree,
                               ComponentTreeFZ<GraphT>* maxtree,
                               std::vector<NodeId>& nodesToPruning,
                               IterationTimes* iter,
                               BenchmarkTimers* timers) {
    (void)timers;
    long long updateUs = 0;
    long long pruningUs = 0;
    Stopwatch sw;
    for (NodeId nodeId : nodesToPruning) {
        for (NodeId leafId : maxtree->getIteratorPostOrderTraversalById(nodeId)) {
            assert(leafId != maxtree->getRootById(leafId) && "Lmax is root");
            assert(maxtree->isLeafById(leafId) && "Lmax não é uma folha");
            if (leafId == InvalidNode) {
                continue;
            }
            sw.start();
            adjust.updateTree(mintree, leafId);
            sw.pause();
            updateUs += elapsedUs(sw);

            sw.start();
            adjust.prunning(maxtree, leafId);
            sw.pause();
            pruningUs += elapsedUs(sw);
        }
    }
    if (iter) {
        iter->phase1Update += (updateUs + 500) / 1000;
        iter->phase1Pruning += (pruningUs + 500) / 1000;
    }
}

template <typename Computer, typename GraphT>
static void adjustMinTreeTimed(BenchAdjuster<ComponentTreeAdjustmentBySubtree<Computer, GraphT>>& adjust,
                               ComponentTreeFZ<GraphT>* mintree,
                               ComponentTreeFZ<GraphT>* maxtree,
                               std::vector<NodeId>& nodesToPruning,
                               IterationTimes* iter,
                               BenchmarkTimers* timers) {
    long long updateUs = 0;
    long long pruningUs = 0;
    Stopwatch sw;
    SubtreeMetrics metrics;
    if (iter && timers) {
        adjust.setMetrics(&metrics, timers, pauseBenchmarkTimers, resumeBenchmarkTimers);
    } else {
        adjust.setMetrics(nullptr);
    }
    for (NodeId rSubtreeId : nodesToPruning) {
        assert(rSubtreeId != maxtree->getRootById() && "rSubtree is root");
        if (rSubtreeId == InvalidNode) {
            continue;
        }
        sw.start();
        adjust.updateTree(mintree, rSubtreeId);
        sw.pause();
        updateUs += elapsedUs(sw);

        if (iter) {
            addSubtreeMetrics(iter->phase1SubtreeMetrics, metrics);
        }

        sw.start();
        adjust.prunning(maxtree, rSubtreeId);
        sw.pause();
        pruningUs += elapsedUs(sw);
    }
    if (iter) {
        iter->phase1Update += (updateUs + 500) / 1000;
        iter->phase1Pruning += (pruningUs + 500) / 1000;
    }
}

template <typename AdjusterType, typename GraphT>
static void adjustMaxTreeTimed(AdjusterType& adjust,
                               ComponentTreeFZ<GraphT>* maxtree,
                               ComponentTreeFZ<GraphT>* mintree,
                               std::vector<NodeId>& nodesToPruning,
                               IterationTimes* iter,
                               BenchmarkTimers* timers) {
    (void)timers;
    long long updateUs = 0;
    long long pruningUs = 0;
    Stopwatch sw;
    for (NodeId rSubtreeId : nodesToPruning) {
        assert(rSubtreeId != mintree->getRootById() && "rSubtree is root");
        if (rSubtreeId == InvalidNode) {
            continue;
        }
        sw.start();
        adjust.updateTree(maxtree, rSubtreeId);
        sw.pause();
        updateUs += elapsedUs(sw);

        sw.start();
        adjust.prunning(mintree, rSubtreeId);
        sw.pause();
        pruningUs += elapsedUs(sw);
    }
    if (iter) {
        iter->phase2Update += (updateUs + 500) / 1000;
        iter->phase2Pruning += (pruningUs + 500) / 1000;
    }
}

template <typename Computer, typename GraphT>
static void adjustMaxTreeTimed(BenchAdjuster<ComponentTreeAdjustmentByLeaf<Computer, GraphT>>& adjust,
                               ComponentTreeFZ<GraphT>* maxtree,
                               ComponentTreeFZ<GraphT>* mintree,
                               std::vector<NodeId>& nodesToPruning,
                               IterationTimes* iter,
                               BenchmarkTimers* timers) {
    (void)timers;
    long long updateUs = 0;
    long long pruningUs = 0;
    Stopwatch sw;
    for (NodeId nodeId : nodesToPruning) {
        for (NodeId leafId : mintree->getIteratorPostOrderTraversalById(nodeId)) {
            assert(leafId != mintree->getRootById(leafId) && "Lmin is root");
            assert(mintree->isLeafById(leafId) && "Lmin não é uma folha");
            if (leafId == InvalidNode) {
                continue;
            }
            sw.start();
            adjust.updateTree(maxtree, leafId);
            sw.pause();
            updateUs += elapsedUs(sw);

            sw.start();
            adjust.prunning(mintree, leafId);
            sw.pause();
            pruningUs += elapsedUs(sw);
        }
    }
    if (iter) {
        iter->phase2Update += (updateUs + 500) / 1000;
        iter->phase2Pruning += (pruningUs + 500) / 1000;
    }
}

template <typename Computer, typename GraphT>
static void adjustMaxTreeTimed(BenchAdjuster<ComponentTreeAdjustmentBySubtree<Computer, GraphT>>& adjust,
                               ComponentTreeFZ<GraphT>* maxtree,
                               ComponentTreeFZ<GraphT>* mintree,
                               std::vector<NodeId>& nodesToPruning,
                               IterationTimes* iter,
                               BenchmarkTimers* timers) {
    long long updateUs = 0;
    long long pruningUs = 0;
    Stopwatch sw;
    SubtreeMetrics metrics;
    if (iter && timers) {
        adjust.setMetrics(&metrics, timers, pauseBenchmarkTimers, resumeBenchmarkTimers);
    } else {
        adjust.setMetrics(nullptr);
    }
    for (NodeId rSubtreeId : nodesToPruning) {
        assert(rSubtreeId != mintree->getRootById() && "rSubtree is root");
        if (rSubtreeId == InvalidNode) {
            continue;
        }
        sw.start();
        adjust.updateTree(maxtree, rSubtreeId);
        sw.pause();
        updateUs += elapsedUs(sw);

        if (iter) {
            addSubtreeMetrics(iter->phase2SubtreeMetrics, metrics);
        }

        sw.start();
        adjust.prunning(mintree, rSubtreeId);
        sw.pause();
        pruningUs += elapsedUs(sw);
    }
    if (iter) {
        iter->phase2Update += (updateUs + 500) / 1000;
        iter->phase2Pruning += (pruningUs + 500) / 1000;
    }
}

template <typename AdjusterType, typename GraphT>
static void runAdjustmentLoop(AdjusterType& adjust,
                              ComponentTreeFZ<GraphT>* maxtree,
                              ComponentTreeFZ<GraphT>* mintree,
                              std::vector<float>& attributeMax,
                              std::vector<float>& attributeMin,
                              const std::vector<int>& thresholds,
                              size_t startIndex,
                              ImageUInt8Ptr img,
                              BenchmarkTimers& timers,
                              const char* labelUpdateMin,
                              const char* labelUpdateMax,
                              std::vector<IterationTimes>* iterTimes) {
    for (size_t i = startIndex; i < thresholds.size(); ++i) {
        int threshold = thresholds[i];
        timers.printIterationHeader(i, threshold);
        timers.startIteration();

        IterationTimes* iter = nullptr;
        if (iterTimes && i < iterTimes->size()) {
            iter = &(*iterTimes)[i];
        }
        auto nodesToPruning = getNodesThreshold(maxtree,
                                                attributeMax,
                                                threshold,
                                                &timers,
                                                img,
                                                iter ? &iter->phase1Metrics : nullptr);
        adjustMinTreeTimed(adjust, mintree, maxtree, nodesToPruning, iter, &timers);
        timers.printPhase(labelUpdateMin);

        nodesToPruning = getNodesThreshold(mintree,
                                           attributeMin,
                                           threshold,
                                           &timers,
                                           img,
                                           iter ? &iter->phase2Metrics : nullptr);
        adjustMaxTreeTimed(adjust, maxtree, mintree, nodesToPruning, iter, &timers);
        timers.printPhaseAndTotal(labelUpdateMax);
    }
}

ImageUInt8Ptr computerCASF_naive(ImageUInt8Ptr img,
                                 double radioAdj,
                                 const std::vector<int>& thresholds,
                                 AttributeMode attributeMode,
                                 Stopwatch* external = nullptr,
                                 std::vector<IterationTimes>* iterTimes = nullptr,
                                 BuildTimes* buildTimes = nullptr) {
    (void)buildTimes;
    BenchmarkTimers timers(external, iterTimes);
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    
    for (size_t i = 0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        timers.printIterationHeader(i, threshold);
        timers.startIteration();

        IterationTimes* iter = nullptr;
        if (iterTimes && i < iterTimes->size()) {
            iter = &(*iterTimes)[i];
        }
        imgOut = runPrunePass(imgOut, adj, true, threshold, attributeMode, &timers,
                              iter ? &iter->phase1Metrics : nullptr);
        timers.printPhase("\t- Time (build/prunning/rec maxtree) naive: ");

        imgOut = runPrunePass(imgOut, adj, false, threshold, attributeMode, &timers,
                              iter ? &iter->phase2Metrics : nullptr);
        timers.printPhaseAndTotal("\t- Time (build/prunning/rec mintree) naive: ");
  	    
	}
    
    return imgOut;
}


template <typename GraphT>
ImageUInt8Ptr computerCASF_hybrid(ImageUInt8Ptr img,
                                  double radioAdj,
                                  const std::vector<int>& thresholds,
                                  size_t cutoffPointHybrid,
                                  AttributeMode attributeMode,
                                  Stopwatch* external = nullptr,
                                  std::vector<IterationTimes>* iterTimes = nullptr,
                                  BuildTimes* buildTimes = nullptr) {
    BenchmarkTimers timers(external, iterTimes);
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();

    for (size_t i = 0; i < cutoffPointHybrid; i++) {
        int threshold = thresholds[i];
        timers.printIterationHeader(i, threshold);
        timers.startIteration();

        IterationTimes* iter = nullptr;
        if (iterTimes && i < iterTimes->size()) {
            iter = &(*iterTimes)[i];
        }
        imgOut = runPrunePass(imgOut, adj, true, threshold, attributeMode, &timers,
                              iter ? &iter->phase1Metrics : nullptr);
        timers.printPhase("\t- Time (update mintree and pruning maxtree) hybrid: ");

        imgOut = runPrunePass(imgOut, adj, false, threshold, attributeMode, &timers,
                              iter ? &iter->phase2Metrics : nullptr);
        timers.printPhaseAndTotal("\t- Time (update maxtree and pruning mintree) hybrid: ");

    }

    timers.startPhase();
    Stopwatch graphSw;
    Stopwatch maxTreeSw;
    Stopwatch minTreeSw;
    if (buildTimes) {
        graphSw.start();
    }
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(imgOut, adj);
    if (buildTimes) {
        graphSw.pause();
        maxTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    if (buildTimes) {
        maxTreeSw.pause();
        minTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    if (buildTimes) {
        minTreeSw.pause();
        buildTimes->graph = elapsedMs(graphSw);
        buildTimes->maxtree = elapsedMs(maxTreeSw);
        buildTimes->mintree = elapsedMs(minTreeSw);
        buildTimes->total = buildTimes->graph + buildTimes->maxtree + buildTimes->mintree;
        buildTimes->valid = true;
    }
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    if (attributeMode == AttributeMode::Area) {
        BenchAdjuster<ComponentTreeAdjustmentByLeaf<DefaultAttributeComputerT<GraphT>, GraphT>> adjust(mintree, maxtree);
        AreaComputerFZT<GraphT> computerAttrMax(maxtree);
        AreaComputerFZT<GraphT> computerAttrMin(mintree);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

        timers.printBuild(nullptr);
        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          cutoffPointHybrid,
                          imgOut,
                          timers,
                          "\t- Time (update mintree and pruning maxtree) hybrid: ",
                          "\t- Time (update maxtree and pruning mintree) hybrid: ",
                          iterTimes);
    } else {
        const Attribute bboxAttr = toBoundingBoxAttribute(attributeMode);
        BenchAdjuster<ComponentTreeAdjustmentByLeaf<BoundingBoxComputerFZT<GraphT>, GraphT>> adjust(mintree, maxtree);
        BoundingBoxComputerFZT<GraphT> computerAttrMax(maxtree, bboxAttr);
        BoundingBoxComputerFZT<GraphT> computerAttrMin(mintree, bboxAttr);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

        timers.printBuild(nullptr);
        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          cutoffPointHybrid,
                          imgOut,
                          timers,
                          "\t- Time (update mintree and pruning maxtree) hybrid: ",
                          "\t- Time (update maxtree and pruning mintree) hybrid: ",
                          iterTimes);
    }
    imgOut = mintree->reconstructionImage();
    return imgOut;
}


template <typename GraphT>
ImageUInt8Ptr computerCASF(ImageUInt8Ptr img,
                           double radioAdj,
                           const std::vector<int>& thresholds,
                           AttributeMode attributeMode,
                           Stopwatch* external = nullptr,
                           std::vector<IterationTimes>* iterTimes = nullptr,
                           BuildTimes* buildTimes = nullptr) {
    BenchmarkTimers timers(external, iterTimes);
    timers.startPhase();
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    Stopwatch graphSw;
    Stopwatch maxTreeSw;
    Stopwatch minTreeSw;
    if (buildTimes) {
        graphSw.start();
    }
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(img, adj);
    if (buildTimes) {
        graphSw.pause();
        maxTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    if (buildTimes) {
        maxTreeSw.pause();
        minTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    if (buildTimes) {
        minTreeSw.pause();
        buildTimes->graph = elapsedMs(graphSw);
        buildTimes->maxtree = elapsedMs(maxTreeSw);
        buildTimes->mintree = elapsedMs(minTreeSw);
        buildTimes->total = buildTimes->graph + buildTimes->maxtree + buildTimes->mintree;
        buildTimes->valid = true;
    }
    
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();
    
    if (attributeMode == AttributeMode::Area) {
        BenchAdjuster<ComponentTreeAdjustmentByLeaf<DefaultAttributeComputerT<GraphT>, GraphT>> adjust(mintree, maxtree);
        AreaComputerFZT<GraphT> computerAttrMax(maxtree);
        AreaComputerFZT<GraphT> computerAttrMin(mintree);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
        timers.printBuild("\n\tTime (build trees): ");

        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          0,
                          img,
                          timers,
                          "\t- Time (update mintree and pruning maxtree): ",
                          "\t- Time (update maxtree and pruning mintree): ",
                          iterTimes);
    } else {
        const Attribute bboxAttr = toBoundingBoxAttribute(attributeMode);
        BenchAdjuster<ComponentTreeAdjustmentByLeaf<BoundingBoxComputerFZT<GraphT>, GraphT>> adjust(mintree, maxtree);
        BoundingBoxComputerFZT<GraphT> computerAttrMax(maxtree, bboxAttr);
        BoundingBoxComputerFZT<GraphT> computerAttrMin(mintree, bboxAttr);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
        timers.printBuild("\n\tTime (build trees): ");

        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          0,
                          img,
                          timers,
                          "\t- Time (update mintree and pruning maxtree): ",
                          "\t- Time (update maxtree and pruning mintree): ",
                          iterTimes);
    }
    auto imgOut = mintree->reconstructionImage();
    return imgOut;
}


template <typename GraphT>
ImageUInt8Ptr computerCASF_hybridSubtree(ImageUInt8Ptr img,
                                         double radioAdj,
                                         const std::vector<int>& thresholds,
                                         size_t cutoffPointHybrid,
                                         AttributeMode attributeMode,
                                         Stopwatch* external = nullptr,
                                         std::vector<IterationTimes>* iterTimes = nullptr,
                                         BuildTimes* buildTimes = nullptr) {
    BenchmarkTimers timers(external, iterTimes);
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();

    for (size_t i = 0; i < cutoffPointHybrid; i++) {
        int threshold = thresholds[i];
        timers.printIterationHeader(i, threshold);
        timers.startIteration();

        IterationTimes* iter = nullptr;
        if (iterTimes && i < iterTimes->size()) {
            iter = &(*iterTimes)[i];
        }
        imgOut = runPrunePass(imgOut, adj, true, threshold, attributeMode, &timers,
                              iter ? &iter->phase1Metrics : nullptr);
        timers.printPhase("\t- Time (update mintree and pruning maxtree) hybrid: ");

        imgOut = runPrunePass(imgOut, adj, false, threshold, attributeMode, &timers,
                              iter ? &iter->phase2Metrics : nullptr);
        timers.printPhaseAndTotal("\t- Time (update maxtree and pruning mintree) hybrid: ");

    }

    timers.startPhase();
    Stopwatch graphSw;
    Stopwatch maxTreeSw;
    Stopwatch minTreeSw;
    if (buildTimes) {
        graphSw.start();
    }
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(imgOut, adj);
    if (buildTimes) {
        graphSw.pause();
        maxTreeSw.start();
    }

    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    if (buildTimes) {
        maxTreeSw.pause();
        minTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    if (buildTimes) {
        minTreeSw.pause();
        buildTimes->graph = elapsedMs(graphSw);
        buildTimes->maxtree = elapsedMs(maxTreeSw);
        buildTimes->mintree = elapsedMs(minTreeSw);
        buildTimes->total = buildTimes->graph + buildTimes->maxtree + buildTimes->mintree;
        buildTimes->valid = true;
    }
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    if (attributeMode == AttributeMode::Area) {
        BenchAdjuster<ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<GraphT>, GraphT>> adjust(mintree, maxtree);
        AreaComputerFZT<GraphT> computerAttrMax(maxtree);
        AreaComputerFZT<GraphT> computerAttrMin(mintree);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
        timers.printBuild(nullptr);

        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          cutoffPointHybrid,
                          imgOut,
                          timers,
                          "\t- Time (update mintree and pruning maxtree) hybrid: ",
                          "\t- Time (update maxtree and pruning mintree) hybrid: ",
                          iterTimes);
    } else {
        const Attribute bboxAttr = toBoundingBoxAttribute(attributeMode);
        BenchAdjuster<ComponentTreeAdjustmentBySubtree<BoundingBoxComputerFZT<GraphT>, GraphT>> adjust(mintree, maxtree);
        BoundingBoxComputerFZT<GraphT> computerAttrMax(maxtree, bboxAttr);
        BoundingBoxComputerFZT<GraphT> computerAttrMin(mintree, bboxAttr);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
        timers.printBuild(nullptr);

        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          cutoffPointHybrid,
                          imgOut,
                          timers,
                          "\t- Time (update mintree and pruning maxtree) hybrid: ",
                          "\t- Time (update maxtree and pruning mintree) hybrid: ",
                          iterTimes);
    }
    imgOut = mintree->reconstructionImage();
    return imgOut;
}



template <typename GraphT>
ImageUInt8Ptr computerCASF_subtree(ImageUInt8Ptr img,
                                   double radioAdj,
                                   const std::vector<int>& thresholds,
                                   AttributeMode attributeMode,
                                   Stopwatch* external = nullptr,
                                   std::vector<IterationTimes>* iterTimes = nullptr,
                                   BuildTimes* buildTimes = nullptr) {
    BenchmarkTimers timers(external, iterTimes);
    timers.startPhase();
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    Stopwatch graphSw;
    Stopwatch maxTreeSw;
    Stopwatch minTreeSw;
    if (buildTimes) {
        graphSw.start();
    }
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(img, adj);
    if (buildTimes) {
        graphSw.pause();
        maxTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    if (buildTimes) {
        maxTreeSw.pause();
        minTreeSw.start();
    }
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    if (buildTimes) {
        minTreeSw.pause();
        buildTimes->graph = elapsedMs(graphSw);
        buildTimes->maxtree = elapsedMs(maxTreeSw);
        buildTimes->mintree = elapsedMs(minTreeSw);
        buildTimes->total = buildTimes->graph + buildTimes->maxtree + buildTimes->mintree;
        buildTimes->valid = true;
    }
    
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    if (attributeMode == AttributeMode::Area) {
        BenchAdjuster<ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<GraphT>, GraphT>> adjust(mintree, maxtree);
        AreaComputerFZT<GraphT> computerAttrMax(maxtree);
        AreaComputerFZT<GraphT> computerAttrMin(mintree);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

        timers.printBuild("\n\tTime (build trees): ");

        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          0,
                          img,
                          timers,
                          "\t- Time (update mintree and pruning maxtree): ",
                          "\t- Time (update maxtree and pruning mintree): ",
                          iterTimes);
    } else {
        const Attribute bboxAttr = toBoundingBoxAttribute(attributeMode);
        BenchAdjuster<ComponentTreeAdjustmentBySubtree<BoundingBoxComputerFZT<GraphT>, GraphT>> adjust(mintree, maxtree);
        BoundingBoxComputerFZT<GraphT> computerAttrMax(maxtree, bboxAttr);
        BoundingBoxComputerFZT<GraphT> computerAttrMin(mintree, bboxAttr);
        std::vector<float> attributeMax = computerAttrMax.compute();
        std::vector<float> attributeMin = computerAttrMin.compute();
        adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

        timers.printBuild("\n\tTime (build trees): ");

        runAdjustmentLoop(adjust,
                          maxtree,
                          mintree,
                          attributeMax,
                          attributeMin,
                          thresholds,
                          0,
                          img,
                          timers,
                          "\t- Time (update mintree and pruning maxtree): ",
                          "\t- Time (update maxtree and pruning mintree): ",
                          iterTimes);
    }
    auto imgOut = mintree->reconstructionImage();
    return imgOut;
}


int main(int argc, char* argv[]) {
    BenchOptions options;
    if (!parseArgs(argc, argv, options)) {
        printUsage(argv[0]);
        return 1;
    }
    if (options.showHelp) {
        printUsage(argv[0]);
        return 0;
    }
    if (options.images.empty()) {
        options.images.emplace_back("Test");
    }
    gEnableStats = options.stats;
    gCollectMetrics = options.stats || options.iterStats;

    std::vector<ImageEntry> loadedImages;
    loadedImages.reserve(options.images.size());
    for (const auto& filename : options.images) {
        loadedImages.push_back(loadImageEntry(filename));
    }

    //std::vector<int> thresholds = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
    //std::vector<int> thresholds = {111, 222, 333, 445, 556, 667, 778, 889, 1000, 1111, 1222, 1333, 1445, 1556, 1667, 1778, 1889, 2000, 2111, 2222, 2333, 2445, 2556, 2667, 2778, 2889, 3000, 3111, 3222, 3333, 3445, 3556, 3667, 3778, 3889, 4000, 4111, 4222, 4333, 4445, 4556, 4667, 4778, 4889, 5000, 5111, 5222, 5333, 5445, 5556, 5667};
    std::vector<int> thresholds = options.thresholds;
    if (thresholds.empty()) {
        thresholds = {50, 100, 150, 200, 250, 300, 350, 400};//, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500};
    }

    const double radioAdj = options.radioAdj;
    const std::array<std::string, 13> kDefaultMethods = {
        "naive",
        "our",
        "our_hybrid",
        "our_subtree",
        "our_hybrid_subtree",
        "our_eager",
        "our_hybrid_eager",
        "our_subtree_eager",
        "our_hybrid_subtree_eager",
        "our_naive",
        "our_hybrid_naive",
        "our_subtree_naive",
        "our_hybrid_subtree_naive"
    };
    std::vector<std::string> methodsToRun;
    if (options.methods.empty()) {
        methodsToRun.assign(kDefaultMethods.begin(), kDefaultMethods.end());
    } else {
        methodsToRun = options.methods;
    }
    methodsToRun = normalizeMethodsForGraphType(methodsToRun, options.graphType);
    if (methodsToRun.empty()) {
        std::cerr << "Erro: nenhum metodo disponivel para graph_type="
                  << graphTypeToString(options.graphType) << "\n";
        return 1;
    }
    const bool quietRuns = options.quiet || options.json || options.iterStats;
    const char* labelNaiveMax = "\t- Time (build/prunning/rec maxtree) naive: ";
    const char* labelNaiveMin = "\t- Time (build/prunning/rec mintree) naive: ";
    const char* labelUpdateMin = "\t- Time (update mintree and pruning maxtree): ";
    const char* labelUpdateMax = "\t- Time (update maxtree and pruning mintree): ";
    const char* labelUpdateMinHybrid = "\t- Time (update mintree and pruning maxtree) hybrid: ";
    const char* labelUpdateMaxHybrid = "\t- Time (update maxtree and pruning mintree) hybrid: ";
    auto shouldRun = [&](std::string_view name) {
        return containsString(methodsToRun, name);
    };
    auto methodIndex = [&](std::string_view name) -> size_t {
        for (size_t i = 0; i < methodsToRun.size(); ++i) {
            if (methodsToRun[i] == name) return i;
        }
        return methodsToRun.size();
    };

    std::vector<MethodBenchmarkResult> jsonMethods;
    if (options.json && options.iterStats && !options.validateOnly) {
        jsonMethods.reserve(methodsToRun.size());
        for (const auto& method : methodsToRun) {
            MethodBenchmarkResult entry;
            entry.name = method;
            jsonMethods.push_back(std::move(entry));
        }
    }
    std::vector<ValidationResult> validationResults;

    const size_t totalImages = loadedImages.size();
    for (size_t imageIndex = 0; imageIndex < totalImages; ++imageIndex) {
        const auto& entry = loadedImages[imageIndex];
        if (!options.quiet && totalImages > 1) {
            const double percent = (static_cast<double>(imageIndex + 1) * 100.0) /
                                   static_cast<double>(totalImages);
            std::ostringstream progress;
            progress << "Progress: " << (imageIndex + 1) << "/" << totalImages
                     << " (" << std::fixed << std::setprecision(1) << percent << "%)"
                     << " - " << entry.name << "\n";
            if (options.json) {
                std::cerr << progress.str();
            } else {
                std::cout << progress.str();
            }
        }
        if (!options.json) {
            std::cout << "Image: " << entry.name << std::endl;
            std::cout << "Resolution: " << entry.cols << "x" << entry.rows << std::endl;
        }

        ImageMetrics imageMetricsFZ;
        ImageMetrics imageMetricsPixel;
        if (options.json && options.iterStats && !options.validateOnly) {
            AdjacencyRelationPtr metricsAdj = std::make_shared<AdjacencyRelation>(entry.rows, entry.cols, radioAdj);
            std::shared_ptr<DefaultFlatZonesGraph> metricsGraph = std::make_shared<DefaultFlatZonesGraph>(entry.img, metricsAdj);
            ComponentTreeFZPtr<> maxTreeFZPtr = std::make_shared<ComponentTreeFZ<>>(metricsGraph, true);
            ComponentTreeFZPtr<> minTreeFZPtr = std::make_shared<ComponentTreeFZ<>>(metricsGraph, false);
            ComponentTreePPtr maxTreePixelPtr = std::make_shared<ComponentTreeP>(entry.img, true, metricsAdj);
            ComponentTreePPtr minTreePixelPtr = std::make_shared<ComponentTreeP>(entry.img, false, metricsAdj);

            imageMetricsFZ = computeImageMetricsFZ(metricsGraph.get(),
                                                   maxTreeFZPtr.get(),
                                                   minTreeFZPtr.get(),
                                                   entry.img);
            imageMetricsPixel = computeImageMetricsPixel(maxTreePixelPtr.get(),
                                                         minTreePixelPtr.get(),
                                                         entry.img);
        }

        auto metricsForMethod = [&](std::string_view methodName) -> ImageMetrics {
            if (methodName == "naive") {
                return imageMetricsPixel;
            }
            return imageMetricsFZ;
        };

        if (options.validateOnly) {
            ImageUInt8Ptr imgOut1;
            ImageUInt8Ptr imgOut2;
            ImageUInt8Ptr imgOut3;
            ImageUInt8Ptr imgOut4;
            ImageUInt8Ptr imgOut5;
            {
                CoutSilencer silencer(options.quiet || options.json || options.iterStats);
                imgOut1 = computerCASF_naive(entry.img, radioAdj, thresholds, options.attributeMode, nullptr, nullptr);
                imgOut2 = computerCASF<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.attributeMode, nullptr, nullptr);
                imgOut3 = computerCASF_hybrid<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.cutoffHybrid, options.attributeMode, nullptr, nullptr);
                imgOut4 = computerCASF_subtree<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.attributeMode, nullptr, nullptr);
                imgOut5 = computerCASF_hybridSubtree<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.cutoffHybridSub, options.attributeMode, nullptr, nullptr);
            }
            const bool eqLeaf = imgOut1->isEqual(imgOut2);
            const bool eqHybrid = imgOut1->isEqual(imgOut3);
            const bool eqSubtree = imgOut1->isEqual(imgOut4);
            const bool eqHybridSub = imgOut1->isEqual(imgOut5);
            if (options.json) {
                printJsonValidation(entry.name, eqLeaf, eqHybrid, eqSubtree, eqHybridSub);
            } else {
                std::cout << "The images (naive/our) are equals: " << (eqLeaf ? "True" : "False") << "\n";
                std::cout << "The images (naive/our_hybrid) are equals: " << (eqHybrid ? "True" : "False") << "\n";
                std::cout << "The images (naive/our_subtree) are equals: " << (eqSubtree ? "True" : "False") << "\n";
                std::cout << "The images (naive/our_hybridSubtree) are equals: " << (eqHybridSub ? "True" : "False") << "\n\n";
            }
            continue;
        }

        MethodResult resNaive;
        if (shouldRun("naive")) {
            resNaive = runBenchmarkMethod(
                "naive",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelNaiveMax,
                labelNaiveMin,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_naive(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resNaive);
                } else {
                    printStatsSummary(resNaive);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resNaive, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("naive");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("naive");
                    imgResult.result = resNaive;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resLeaf;
        if (shouldRun("our")) {
            resLeaf = runBenchmarkMethod(
                "our",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMin,
                labelUpdateMax,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resLeaf);
                } else {
                    printStatsSummary(resLeaf);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resLeaf, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our");
                    imgResult.result = resLeaf;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resHybrid;
        if (shouldRun("our_hybrid")) {
            resHybrid = runBenchmarkMethod(
                "our_hybrid",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMinHybrid,
                labelUpdateMaxHybrid,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_hybrid<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.cutoffHybrid, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resHybrid);
                } else {
                    printStatsSummary(resHybrid);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resHybrid, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_hybrid");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_hybrid");
                    imgResult.result = resHybrid;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resSubtree;
        if (shouldRun("our_subtree")) {
            resSubtree = runBenchmarkMethod(
                "our_subtree",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMin,
                labelUpdateMax,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_subtree<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resSubtree);
                } else {
                    printStatsSummary(resSubtree);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resSubtree, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_subtree");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_subtree");
                    imgResult.result = resSubtree;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resHybridSub;
        if (shouldRun("our_hybrid_subtree")) {
            resHybridSub = runBenchmarkMethod(
                "our_hybrid_subtree",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMinHybrid,
                labelUpdateMaxHybrid,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_hybridSubtree<DefaultFlatZonesGraph>(entry.img, radioAdj, thresholds, options.cutoffHybridSub, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resHybridSub);
                } else {
                    printStatsSummary(resHybridSub);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resHybridSub, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_hybrid_subtree");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_hybrid_subtree");
                    imgResult.result = resHybridSub;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resLeafEager;
        if (shouldRun("our_eager")) {
            resLeafEager = runBenchmarkMethod(
                "our_eager",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMin,
                labelUpdateMax,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF<FlatZonesGraphFullEdges>(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resLeafEager);
                } else {
                    printStatsSummary(resLeafEager);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resLeafEager, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_eager");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_eager");
                    imgResult.result = resLeafEager;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resHybridEager;
        if (shouldRun("our_hybrid_eager")) {
            resHybridEager = runBenchmarkMethod(
                "our_hybrid_eager",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMinHybrid,
                labelUpdateMaxHybrid,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_hybrid<FlatZonesGraphFullEdges>(entry.img, radioAdj, thresholds, options.cutoffHybrid, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resHybridEager);
                } else {
                    printStatsSummary(resHybridEager);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resHybridEager, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_hybrid_eager");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_hybrid_eager");
                    imgResult.result = resHybridEager;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resSubtreeEager;
        if (shouldRun("our_subtree_eager")) {
            resSubtreeEager = runBenchmarkMethod(
                "our_subtree_eager",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMin,
                labelUpdateMax,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_subtree<FlatZonesGraphFullEdges>(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resSubtreeEager);
                } else {
                    printStatsSummary(resSubtreeEager);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resSubtreeEager, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_subtree_eager");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_subtree_eager");
                    imgResult.result = resSubtreeEager;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resHybridSubEager;
        if (shouldRun("our_hybrid_subtree_eager")) {
            resHybridSubEager = runBenchmarkMethod(
                "our_hybrid_subtree_eager",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMinHybrid,
                labelUpdateMaxHybrid,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_hybridSubtree<FlatZonesGraphFullEdges>(entry.img, radioAdj, thresholds, options.cutoffHybridSub, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resHybridSubEager);
                } else {
                    printStatsSummary(resHybridSubEager);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resHybridSubEager, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_hybrid_subtree_eager");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_hybrid_subtree_eager");
                    imgResult.result = resHybridSubEager;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resLeafNaive;
        if (shouldRun("our_naive")) {
            resLeafNaive = runBenchmarkMethod(
                "our_naive",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMin,
                labelUpdateMax,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF<FlatZonesGraphOnDemandEdgesByPixel>(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resLeafNaive);
                } else {
                    printStatsSummary(resLeafNaive);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resLeafNaive, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_naive");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_naive");
                    imgResult.result = resLeafNaive;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resHybridNaive;
        if (shouldRun("our_hybrid_naive")) {
            resHybridNaive = runBenchmarkMethod(
                "our_hybrid_naive",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMinHybrid,
                labelUpdateMaxHybrid,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_hybrid<FlatZonesGraphOnDemandEdgesByPixel>(entry.img, radioAdj, thresholds, options.cutoffHybrid, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resHybridNaive);
                } else {
                    printStatsSummary(resHybridNaive);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resHybridNaive, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_hybrid_naive");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_hybrid_naive");
                    imgResult.result = resHybridNaive;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resSubtreeNaive;
        if (shouldRun("our_subtree_naive")) {
            resSubtreeNaive = runBenchmarkMethod(
                "our_subtree_naive",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMin,
                labelUpdateMax,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_subtree<FlatZonesGraphOnDemandEdgesByPixel>(entry.img, radioAdj, thresholds, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resSubtreeNaive);
                } else {
                    printStatsSummary(resSubtreeNaive);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resSubtreeNaive, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_subtree_naive");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_subtree_naive");
                    imgResult.result = resSubtreeNaive;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        MethodResult resHybridSubNaive;
        if (shouldRun("our_hybrid_subtree_naive")) {
            resHybridSubNaive = runBenchmarkMethod(
                "our_hybrid_subtree_naive",
                options.warmup,
                options.repeat,
                quietRuns,
                options.iterStats,
                labelUpdateMinHybrid,
                labelUpdateMaxHybrid,
                [&](Stopwatch* sw, std::vector<IterationTimes>* iterTimes, BuildTimes* buildTimes) {
                    return computerCASF_hybridSubtree<FlatZonesGraphOnDemandEdgesByPixel>(entry.img, radioAdj, thresholds, options.cutoffHybridSub, options.attributeMode, sw, iterTimes, buildTimes);
                });
            if (!options.json) {
                if (options.iterStats) {
                    printIterationStats(resHybridSubNaive);
                } else {
                    printStatsSummary(resHybridSubNaive);
                }
            }
            if (options.json && !options.iterStats) {
                printJsonRows(entry.name, resHybridSubNaive, entry.rows, entry.cols);
            }
            if (options.json && options.iterStats && !options.validateOnly) {
                const size_t idx = methodIndex("our_hybrid_subtree_naive");
                if (idx < jsonMethods.size()) {
                    ImageBenchmarkResult imgResult;
                    imgResult.image = entry.name;
                    imgResult.rows = entry.rows;
                    imgResult.cols = entry.cols;
                    imgResult.metrics = metricsForMethod("our_hybrid_subtree_naive");
                    imgResult.result = resHybridSubNaive;
                    imgResult.result.output = nullptr;
                    jsonMethods[idx].images.push_back(std::move(imgResult));
                }
            }
        }

        if (!options.noValidate && resNaive.output && resLeaf.output && resHybrid.output && resSubtree.output && resHybridSub.output) {
            const bool eqLeaf = resNaive.output->isEqual(resLeaf.output);
            const bool eqHybrid = resNaive.output->isEqual(resHybrid.output);
            const bool eqSubtree = resNaive.output->isEqual(resSubtree.output);
            const bool eqHybridSub = resNaive.output->isEqual(resHybridSub.output);
            if (options.json && options.iterStats && !options.validateOnly) {
                ValidationResult validation;
                validation.image = entry.name;
                validation.eqLeaf = eqLeaf;
                validation.eqHybrid = eqHybrid;
                validation.eqSubtree = eqSubtree;
                validation.eqHybridSub = eqHybridSub;
                validationResults.push_back(std::move(validation));
            } else if (options.json) {
                printJsonValidation(entry.name, eqLeaf, eqHybrid, eqSubtree, eqHybridSub);
            } else {
                std::cout << "The images (naive/our) are equals: " << (eqLeaf ? "True" : "False") << "\n";
                std::cout << "The images (naive/our_hybrid) are equals: " << (eqHybrid ? "True" : "False") << "\n";
                std::cout << "The images (naive/our_subtree) are equals: " << (eqSubtree ? "True" : "False") << "\n";
                std::cout << "The images (naive/our_hybridSubtree) are equals: " << (eqHybridSub ? "True" : "False") << "\n\n";
            }
        }
    }

    if (options.json && options.iterStats && !options.validateOnly) {
        printJsonBenchmarkSummary(jsonMethods,
                                  options.attributeMode,
                                  thresholds,
                                  radioAdj,
                                  validationResults);
    }

    /*
    data = new unsigned char[n];./
    for (int i = 0; i < numCols * numRows; i++) {
        data[i] = static_cast<int>(imgOut3[i]);  // Converte de `unsigned char` para `int`
    }
    stbi_write_png(("/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/out_our_" + entry.path().filename().string()).c_str(), numCols, numRows, 1, data, 0);
    delete[] data;
    */

    return 0;
}
