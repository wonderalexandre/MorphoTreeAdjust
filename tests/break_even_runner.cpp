#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <thread>
#include <cstdio>
#include <filesystem>
#include <unistd.h>
#include <sstream>
#include <string>
#include <string_view>
#include <span>
#include <vector>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentBySubtree.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/FlatZonesGraph.hpp"
#include "../morphoTreeAdjust/include/AttributeComputer.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"

static inline double elapsedMs(const Stopwatch& sw) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(sw.elapsed()).count();
}

static std::string gSelfPath;

enum class GraphType {
    OnDemand,
    Eager,
    OnPixel
};

static bool parseGraphType(std::string_view value, GraphType* out) {
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

static void printUsage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " --image <path> --thresholds <csv> --naive-cum <csv> [options]\n"
              << "Options:\n"
              << "  --repeat <N>       Timed runs (default 3)\n"
              << "  --warmup <N>       Warmup runs (default 1)\n"
              << "  --radio-adj <val>  Adjacency radius (default 1.5)\n"
              << "  --graph-type <t>  on_demand|eager|on_pixel (default on_pixel)\n"
              << "  --k0-window <N>    Window around m* for k0* (default 5; 0=full)\n"
              << "  --cooldown-ms <N>  Sleep between k0 evaluations (default 1000)\n"
              << "  --no-k0-process    Disable per-k0 subprocess (default: enabled)\n"
              << "  --no-k0            Skip k0* computation (allows early stop on break-even)\n"
              << "  --quiet            Reduce stderr output\n";
}

static std::vector<int> parseCsvInts(const std::string& csv) {
    std::vector<int> out;
    std::stringstream ss(csv);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) continue;
        out.push_back(std::stoi(token));
    }
    return out;
}

static std::vector<double> parseCsvDoubles(const std::string& csv) {
    std::vector<double> out;
    std::stringstream ss(csv);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) continue;
        out.push_back(std::stod(token));
    }
    return out;
}

static double medianOf(std::vector<double> values) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    if (n % 2 == 1) {
        return values[n / 2];
    }
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
}

static int medianOf(std::vector<int> values) {
    if (values.empty()) return -1;
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    if (n % 2 == 1) {
        return values[n / 2];
    }
    return values[n / 2 - 1];
}

struct RunResult {
    int m_star = -1;
    int k0_star = -1;
    int k0_threshold = -1;
    int threshold = -1;
    double naive_cum_ms = 0.0;
    double update_cum_ms = 0.0;
    double pre_ms = 0.0;
    double post_ms = 0.0;
    double r_hybrid_ms = 0.0;
    size_t iterations_executed = 0;
    int k0_window_used = 0;
    int cooldown_ms_used = 0;
};

struct NaiveRun {
    std::vector<double> iter_ms;
    std::vector<ImageUInt8Ptr> images;
    std::vector<std::string> image_paths;
    std::string temp_dir;
};

struct UpdateRun {
    double pre_ms = 0.0;
    double iter_ms = 0.0;
    double post_ms = 0.0;
};

template <typename CNPsType, typename GraphT>
static std::vector<NodeId> getNodesThreshold(ComponentTree<CNPsType, GraphT>* tree,
                                             std::span<float> attribute,
                                             int threshold) {
    std::vector<NodeId> lista;
    FastQueue<NodeId> queue;
    queue.push(tree->getRootById());

    while (!queue.empty()) {
        NodeId id = queue.pop();
        if (attribute[id] > threshold) {
            for (NodeId c : tree->getChildrenById(id)) {
                queue.push(c);
            }
        } else {
            lista.push_back(id);
        }
    }
    return lista;
}

static ImageUInt8Ptr runPrunePass(ImageUInt8Ptr imgIn,
                                  AdjacencyRelationPtr adj,
                                  bool isMaxTree,
                                  int threshold) {
    ComponentTreePPtr treePtr = std::make_shared<ComponentTreeP>(imgIn, isMaxTree, adj);
    ComponentTreeP* tree = treePtr.get();
    AreaComputerP computerAttr(tree);
    std::vector<float> attribute = computerAttr.compute();
    for (NodeId node : getNodesThreshold(tree, attribute, threshold)) {
        tree->prunning(node);
    }
    return tree->reconstructionImage();
}

static std::string shellQuote(const std::string& value) {
    std::string out = "'";
    for (char c : value) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out += c;
        }
    }
    out += "'";
    return out;
}

static UpdateRun parseUpdateRun(const std::string& raw) {
    UpdateRun run;
    if (raw.empty()) return run;
    auto findValue = [&](const char* key) -> double {
        const std::string pattern = std::string("\"") + key + "\":";
        size_t pos = raw.find(pattern);
        if (pos == std::string::npos) return 0.0;
        pos += pattern.size();
        size_t end = pos;
        while (end < raw.size() && (std::isdigit(raw[end]) || raw[end] == '.' || raw[end] == '-')) {
            ++end;
        }
        return std::atof(raw.substr(pos, end - pos).c_str());
    };
    run.pre_ms = findValue("pre_ms");
    run.iter_ms = findValue("iter_ms");
    run.post_ms = findValue("post_ms");
    return run;
}

static UpdateRun runUpdateFromImageChild(const std::string& imagePath,
                                         double radioAdj,
                                         const std::vector<int>& thresholds,
                                         size_t startIndex,
                                         GraphType graphType) {
    std::stringstream thresholdsCsv;
    for (size_t i = 0; i < thresholds.size(); ++i) {
        if (i > 0) thresholdsCsv << ",";
        thresholdsCsv << thresholds[i];
    }
    std::stringstream cmd;
    cmd << shellQuote(gSelfPath)
        << " --eval-update"
        << " --image " << shellQuote(imagePath)
        << " --thresholds " << shellQuote(thresholdsCsv.str())
        << " --start-index " << startIndex
        << " --radio-adj " << radioAdj
        << " --graph-type ";
    switch (graphType) {
        case GraphType::OnDemand: cmd << "on_demand"; break;
        case GraphType::Eager: cmd << "eager"; break;
        case GraphType::OnPixel: cmd << "on_pixel"; break;
    }
    FILE* pipe = popen(cmd.str().c_str(), "r");
    if (!pipe) {
        return UpdateRun{};
    }
    std::string output;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe)) {
        output += buffer;
    }
    pclose(pipe);
    return parseUpdateRun(output);
}

static NaiveRun runNaive(ImageUInt8Ptr img,
                         double radioAdj,
                         const std::vector<int>& thresholds,
                         size_t maxIndex,
                         bool storeImages,
                         bool writeImages) {
    NaiveRun run;
    if (!img) return run;
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    const size_t limit = std::min(maxIndex, thresholds.size());
    run.iter_ms.reserve(limit);
    if (storeImages) {
        run.images.reserve(limit);
    }
    if (writeImages) {
        run.image_paths.reserve(limit);
        char tmpl[] = "/tmp/mta_break_even_XXXXXX";
        char* dir = mkdtemp(tmpl);
        if (dir) {
            run.temp_dir = dir;
        }
    }
    for (size_t i = 0; i < limit; ++i) {
        const int threshold = thresholds[i];
        Stopwatch iterSw;
        iterSw.start();
        imgOut = runPrunePass(imgOut, adj, true, threshold);
        imgOut = runPrunePass(imgOut, adj, false, threshold);
        iterSw.pause();
        run.iter_ms.push_back(elapsedMs(iterSw));
        if (storeImages) {
            run.images.push_back(imgOut);
        }
        if (writeImages && !run.temp_dir.empty()) {
            std::ostringstream path;
            path << run.temp_dir << "/naive_" << i << ".png";
            const std::string outPath = path.str();
            stbi_write_png(outPath.c_str(),
                           imgOut->getNumCols(),
                           imgOut->getNumRows(),
                           1,
                           imgOut->rawData(),
                           imgOut->getNumCols());
            run.image_paths.push_back(outPath);
        }
    }
    return run;
}

template <typename GraphT>
static UpdateRun runUpdateFromImage(ImageUInt8Ptr img,
                                    double radioAdj,
                                    const std::vector<int>& thresholds,
                                    size_t startIndex) {
    UpdateRun run;
    if (!img || startIndex >= thresholds.size()) return run;

    Stopwatch preSw;
    preSw.start();
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(img, adj);
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<GraphT>, GraphT> adjust(mintree, maxtree);
    AreaComputerFZT<GraphT> computerAttrMax(maxtree);
    AreaComputerFZT<GraphT> computerAttrMin(mintree);
    std::vector<float> attributeMax = computerAttrMax.compute();
    std::vector<float> attributeMin = computerAttrMin.compute();
    adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
    preSw.pause();
    run.pre_ms = elapsedMs(preSw);

    double iterSum = 0.0;
    for (size_t i = startIndex; i < thresholds.size(); ++i) {
        const int threshold = thresholds[i];
        Stopwatch iterSw;
        iterSw.start();
        auto nodesToPruning = getNodesThreshold(maxtree, attributeMax, threshold);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        nodesToPruning = getNodesThreshold(mintree, attributeMin, threshold);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning);
        iterSw.pause();
        iterSum += elapsedMs(iterSw);
    }
    run.iter_ms = iterSum;

    Stopwatch postSw;
    postSw.start();
    (void)mintree->reconstructionImage();
    postSw.pause();
    run.post_ms = elapsedMs(postSw);
    return run;
}

static ImageUInt8Ptr loadImage(const std::string& filename) {
    int numCols = 0;
    int numRows = 0;
    int nchannels = 0;
    uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    (void)nchannels;
    if (!data) {
        std::cerr << "Erro: Não foi possível carregar a imagem " << filename << "\n";
        return nullptr;
    }
    ImageUInt8Ptr img = ImageUInt8::create(numRows, numCols);
    std::copy(data, data + (numRows * numCols), img->rawData());
    stbi_image_free(data);
    return img;
}

template <typename GraphT>
static RunResult runOnce(ImageUInt8Ptr img,
                         const std::string& imagePath,
                         double radioAdj,
                         const std::vector<int>& thresholds,
                         const std::vector<double>& naiveCum,
                         bool allowEarlyStop,
                         bool computeK0,
                         GraphType graphType,
                         int k0Window,
                         int cooldownMs,
                         bool useK0Process) {
    RunResult result;
    if (!img) return result;

    Stopwatch preSw;
    preSw.start();
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(img, adj);
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<GraphT>, GraphT> adjust(mintree, maxtree);
    AreaComputerFZT<GraphT> computerAttrMax(maxtree);
    AreaComputerFZT<GraphT> computerAttrMin(mintree);
    std::vector<float> attributeMax = computerAttrMax.compute();
    std::vector<float> attributeMin = computerAttrMin.compute();
    adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
    preSw.pause();
    result.pre_ms = elapsedMs(preSw);

    double iterSum = 0.0;
    double postMs = -1.0;
    std::vector<double> updatePrefix;
    updatePrefix.reserve(thresholds.size());

    for (size_t i = 0; i < thresholds.size(); ++i) {
        const int threshold = thresholds[i];

        Stopwatch iterSw;
        iterSw.start();
        auto nodesToPruning = getNodesThreshold(maxtree, attributeMax, threshold);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);

        nodesToPruning = getNodesThreshold(mintree, attributeMin, threshold);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning);
        iterSw.pause();

        const double iterMs = elapsedMs(iterSw);
        iterSum += iterMs;
        updatePrefix.push_back(iterSum);
        result.iterations_executed = i + 1;

        if (!allowEarlyStop) {
            continue;
        }

        if (i < naiveCum.size() && (result.pre_ms + iterSum) <= naiveCum[i]) {
            if (postMs < 0.0) {
                Stopwatch postSw;
                postSw.start();
                (void)mintree->reconstructionImage();
                postSw.pause();
                postMs = elapsedMs(postSw);
            }
            const double updateCum = result.pre_ms + postMs + iterSum;
            if (updateCum <= naiveCum[i]) {
                result.m_star = static_cast<int>(i + 1);
                result.threshold = threshold;
                result.naive_cum_ms = naiveCum[i];
                result.update_cum_ms = updateCum;
                result.post_ms = postMs;
                return result;
            }
        }
    }

    if (postMs < 0.0) {
        Stopwatch postSw;
        postSw.start();
        (void)mintree->reconstructionImage();
        postSw.pause();
        postMs = elapsedMs(postSw);
    }
    result.post_ms = postMs;
    result.update_cum_ms = result.pre_ms + postMs + iterSum;
    if (!thresholds.empty() && !naiveCum.empty()) {
        const size_t last = std::min(thresholds.size(), naiveCum.size()) - 1;
        result.naive_cum_ms = naiveCum[last];
    }

    if (!allowEarlyStop) {
        for (size_t i = 0; i < updatePrefix.size() && i < naiveCum.size(); ++i) {
            const double updateCum = result.pre_ms + postMs + updatePrefix[i];
            if (updateCum <= naiveCum[i]) {
                result.m_star = static_cast<int>(i + 1);
                result.threshold = thresholds[i];
                result.naive_cum_ms = naiveCum[i];
                result.update_cum_ms = updateCum;
                break;
            }
        }
    }

    if (computeK0) {
        size_t k = thresholds.size();
        size_t startK0 = 0;
        size_t endK0 = k;
        if (k0Window > 0 && result.m_star > 0) {
            const int m = result.m_star - 1;
            startK0 = 0;
            endK0 = static_cast<size_t>(std::min<int>(static_cast<int>(k), m + k0Window + 1));
        }
        if (startK0 > endK0) {
            startK0 = 0;
            endK0 = k;
        }
        const size_t maxK0 = endK0;
        NaiveRun naiveRun = runNaive(img, radioAdj, thresholds, maxK0, !useK0Process, useK0Process);
        if (useK0Process && naiveRun.temp_dir.empty()) {
            useK0Process = false;
            naiveRun = runNaive(img, radioAdj, thresholds, maxK0, true, false);
        }
        if (!naiveRun.iter_ms.empty()) {
            std::vector<double> naiveCumLocal;
            naiveCumLocal.reserve(naiveRun.iter_ms.size());
            double s = 0.0;
            for (double v : naiveRun.iter_ms) {
                s += v;
                naiveCumLocal.push_back(s);
            }
            const size_t kFull = thresholds.size();
            double bestR = std::numeric_limits<double>::infinity();
            int bestK0 = 0;
            const size_t evalStart = std::min(startK0, naiveCumLocal.size());
            const size_t evalEnd = std::min(endK0, naiveCumLocal.size());
            for (size_t k0 = evalStart; k0 <= evalEnd; ++k0) {
                const double naivePrefix = (k0 == 0) ? 0.0 : naiveCumLocal[k0 - 1];
                if (k0 == kFull) {
                    if (naivePrefix < bestR) {
                        bestR = naivePrefix;
                        bestK0 = static_cast<int>(k0);
                    }
                    continue;
                }
                ImageUInt8Ptr startImg = nullptr;
                std::string startPath = imagePath;
                if (k0 > 0) {
                    if (useK0Process) {
                        if (k0 - 1 < naiveRun.image_paths.size()) {
                            startPath = naiveRun.image_paths[k0 - 1];
                        } else {
                            continue;
                        }
                    } else {
                        if (k0 - 1 < naiveRun.images.size()) {
                            startImg = naiveRun.images[k0 - 1];
                        } else {
                            continue;
                        }
                    }
                }
                UpdateRun updateRun;
                if (useK0Process) {
                    updateRun = runUpdateFromImageChild(startPath, radioAdj, thresholds, k0, graphType);
                } else {
                    if (!startImg) startImg = img;
                    updateRun = runUpdateFromImage<GraphT>(startImg, radioAdj, thresholds, k0);
                }
                const double total = naivePrefix + updateRun.pre_ms + updateRun.iter_ms + updateRun.post_ms;
                if (total < bestR) {
                    bestR = total;
                    bestK0 = static_cast<int>(k0);
                }
                if (cooldownMs > 0) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(cooldownMs));
                }
            }
            result.k0_star = bestK0;
            result.r_hybrid_ms = bestR;
        }
        if (useK0Process && !naiveRun.temp_dir.empty()) {
            std::error_code ec;
            std::filesystem::remove_all(naiveRun.temp_dir, ec);
        }
    }
    result.k0_window_used = k0Window;
    result.cooldown_ms_used = cooldownMs;
    return result;
}

static void emitJson(const std::string& image,
                     int mStar,
                     int k0Star,
                     int k0Threshold,
                     int threshold,
                     double naiveCum,
                     double updateCum,
                     double preMs,
                     double postMs,
                     double rHybrid,
                     size_t iterationsExecuted,
                     int k0WindowUsed,
                     int cooldownMsUsed) {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(3);
    std::cout << "{"
              << "\"image\":" << "\"" << image << "\""
              << ",\"m_star\":" << mStar
              << ",\"k0_star\":" << k0Star
              << ",\"k0_threshold\":" << k0Threshold
              << ",\"threshold\":" << threshold
              << ",\"naive_cum_ms\":" << naiveCum
              << ",\"update_cum_ms\":" << updateCum
              << ",\"pre_ms\":" << preMs
              << ",\"post_ms\":" << postMs
              << ",\"r_hybrid_ms\":" << rHybrid
              << ",\"iterations_executed\":" << iterationsExecuted
              << ",\"k0_window_used\":" << k0WindowUsed
              << ",\"cooldown_ms\":" << cooldownMsUsed
              << "}\n";
}

int main(int argc, char* argv[]) {
    gSelfPath = argv[0];
    std::string imagePath;
    std::string thresholdsCsv;
    std::string naiveCumCsv;
    int repeat = 3;
    int warmup = 1;
    double radioAdj = 1.5;
    GraphType graphType = GraphType::OnPixel;
    bool quiet = false;
    bool computeK0 = true;
    bool allowEarlyStop = false;
    bool evalUpdateOnly = false;
    size_t startIndex = 0;
    int k0Window = 5;
    int cooldownMs = 1000;
    bool useK0Process = true;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        if (arg == "--eval-update") {
            evalUpdateOnly = true;
            continue;
        }
        if (arg == "--image" && i + 1 < argc) {
            imagePath = argv[++i];
            continue;
        }
        if (arg == "--thresholds" && i + 1 < argc) {
            thresholdsCsv = argv[++i];
            continue;
        }
        if (arg == "--naive-cum" && i + 1 < argc) {
            naiveCumCsv = argv[++i];
            continue;
        }
        if (arg == "--start-index" && i + 1 < argc) {
            startIndex = static_cast<size_t>(std::max(0, std::atoi(argv[++i])));
            continue;
        }
        if (arg == "--repeat" && i + 1 < argc) {
            repeat = std::max(1, std::atoi(argv[++i]));
            continue;
        }
        if (arg == "--warmup" && i + 1 < argc) {
            warmup = std::max(0, std::atoi(argv[++i]));
            continue;
        }
        if (arg == "--radio-adj" && i + 1 < argc) {
            radioAdj = std::atof(argv[++i]);
            continue;
        }
        if (arg == "--graph-type" && i + 1 < argc) {
            if (!parseGraphType(argv[++i], &graphType)) {
                std::cerr << "Erro: graph-type invalido\n";
                return 1;
            }
            continue;
        }
        if (arg == "--k0-window" && i + 1 < argc) {
            k0Window = std::max(0, std::atoi(argv[++i]));
            continue;
        }
        if (arg == "--cooldown-ms" && i + 1 < argc) {
            cooldownMs = std::max(0, std::atoi(argv[++i]));
            continue;
        }
        if (arg == "--no-k0-process") {
            useK0Process = false;
            continue;
        }
        if (arg == "--no-k0") {
            computeK0 = false;
            allowEarlyStop = true;
            continue;
        }
        if (arg == "--quiet") {
            quiet = true;
            continue;
        }
        std::cerr << "Erro: argumento desconhecido " << arg << "\n";
        printUsage(argv[0]);
        return 1;
    }

    if (imagePath.empty() || thresholdsCsv.empty() || naiveCumCsv.empty()) {
        if (!evalUpdateOnly) {
            printUsage(argv[0]);
            return 1;
        }
    }

    const std::vector<int> thresholds = parseCsvInts(thresholdsCsv);
    const std::vector<double> naiveCum = parseCsvDoubles(naiveCumCsv);
    if (thresholds.empty()) {
        std::cerr << "Erro: thresholds vazios\n";
        return 1;
    }
    if (!evalUpdateOnly) {
        if (naiveCum.empty()) {
            std::cerr << "Erro: naive-cum vazio\n";
            return 1;
        }
        if (thresholds.size() != naiveCum.size()) {
            std::cerr << "Erro: thresholds/naive-cum com tamanhos diferentes\n";
            return 1;
        }
    }

    ImageUInt8Ptr img = loadImage(imagePath);
    if (!img) {
        return 1;
    }

    if (evalUpdateOnly) {
        UpdateRun update;
        switch (graphType) {
            case GraphType::OnDemand:
                update = runUpdateFromImage<DefaultFlatZonesGraph>(img, radioAdj, thresholds, startIndex);
                break;
            case GraphType::Eager:
                update = runUpdateFromImage<FlatZonesGraphFullEdges>(img, radioAdj, thresholds, startIndex);
                break;
            case GraphType::OnPixel:
                update = runUpdateFromImage<FlatZonesGraphOnDemandEdgesByPixel>(img, radioAdj, thresholds, startIndex);
                break;
        }
        std::cout.setf(std::ios::fixed);
        std::cout.precision(3);
        std::cout << "{"
                  << "\"pre_ms\":" << update.pre_ms
                  << ",\"iter_ms\":" << update.iter_ms
                  << ",\"post_ms\":" << update.post_ms
                  << "}\n";
        return 0;
    }

    for (int i = 0; i < warmup; ++i) {
        if (!quiet) {
            std::cerr << "Warmup " << (i + 1) << "...\n";
        }
        switch (graphType) {
            case GraphType::OnDemand:
                (void)runOnce<DefaultFlatZonesGraph>(img, imagePath, radioAdj, thresholds, naiveCum, false, computeK0,
                                                     graphType, k0Window, cooldownMs, useK0Process);
                break;
            case GraphType::Eager:
                (void)runOnce<FlatZonesGraphFullEdges>(img, imagePath, radioAdj, thresholds, naiveCum, false, computeK0,
                                                   graphType, k0Window, cooldownMs, useK0Process);
                break;
            case GraphType::OnPixel:
                (void)runOnce<FlatZonesGraphOnDemandEdgesByPixel>(img, imagePath, radioAdj, thresholds, naiveCum, false, computeK0,
                                                   graphType, k0Window, cooldownMs, useK0Process);
                break;
        }
    }

    std::vector<int> mStars;
    std::vector<double> naiveCums;
    std::vector<double> updateCums;
    std::vector<double> preMs;
    std::vector<double> postMs;
    std::vector<double> rHybrid;
    std::vector<size_t> iterationsExecuted;
    std::vector<int> k0Stars;

    for (int i = 0; i < repeat; ++i) {
        if (!quiet) {
            std::cerr << "Run " << (i + 1) << "/" << repeat << "...\n";
        }
        RunResult run;
        switch (graphType) {
            case GraphType::OnDemand:
                run = runOnce<DefaultFlatZonesGraph>(img, imagePath, radioAdj, thresholds, naiveCum, allowEarlyStop,
                                                     computeK0, graphType, k0Window, cooldownMs, useK0Process);
                break;
            case GraphType::Eager:
                run = runOnce<FlatZonesGraphFullEdges>(img, imagePath, radioAdj, thresholds, naiveCum, allowEarlyStop,
                                                   computeK0, graphType, k0Window, cooldownMs, useK0Process);
                break;
            case GraphType::OnPixel:
                run = runOnce<FlatZonesGraphOnDemandEdgesByPixel>(img, imagePath, radioAdj, thresholds, naiveCum, allowEarlyStop,
                                                   computeK0, graphType, k0Window, cooldownMs, useK0Process);
                break;
        }
        mStars.push_back(run.m_star);
        k0Stars.push_back(run.k0_star);
        naiveCums.push_back(run.naive_cum_ms);
        updateCums.push_back(run.update_cum_ms);
        preMs.push_back(run.pre_ms);
        postMs.push_back(run.post_ms);
        rHybrid.push_back(run.r_hybrid_ms);
        iterationsExecuted.push_back(run.iterations_executed);
    }

    const int mStarMedian = medianOf(mStars);
    const int k0StarMedian = medianOf(k0Stars);
    int threshold = -1;
    int k0Threshold = -1;
    double naiveCumMedian = medianOf(naiveCums);
    double updateCumMedian = medianOf(updateCums);
    double preMsMedian = medianOf(preMs);
    double postMsMedian = medianOf(postMs);
    double rHybridMedian = medianOf(rHybrid);
    size_t iterationsMedian = 0;
    if (!iterationsExecuted.empty()) {
        std::vector<double> tmp;
        tmp.reserve(iterationsExecuted.size());
        for (auto v : iterationsExecuted) {
            tmp.push_back(static_cast<double>(v));
        }
        iterationsMedian = static_cast<size_t>(medianOf(tmp));
    }
    if (mStarMedian > 0 && static_cast<size_t>(mStarMedian) <= thresholds.size()) {
        threshold = thresholds[static_cast<size_t>(mStarMedian - 1)];
    }
    if (k0StarMedian > 0 && static_cast<size_t>(k0StarMedian) <= thresholds.size()) {
        k0Threshold = thresholds[static_cast<size_t>(k0StarMedian - 1)];
    }

    emitJson(imagePath,
             mStarMedian,
             k0StarMedian,
             k0Threshold,
             threshold,
             naiveCumMedian,
             updateCumMedian,
             preMsMedian,
             postMsMedian,
             rHybridMedian,
             iterationsMedian,
             k0Window,
             cooldownMs);
    return 0;
}
