#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"

#include "./external/stb/stb_image.h"

struct ImageEntry {
    std::string name;
    ImageUInt8Ptr img;
    int rows = 0;
    int cols = 0;
};

static std::string trimCopy(std::string_view text) {
    size_t start = 0;
    size_t end = text.size();
    while (start < end && std::isspace(static_cast<unsigned char>(text[start]))) {
        ++start;
    }
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }
    return std::string(text.substr(start, end - start));
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

static ImageEntry loadImageEntry(const std::string& filename) {
    ImageEntry entry;
    entry.name = filename;
    int numCols = 0;
    int numRows = 0;
    int nchannels = 0;
    uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    (void)nchannels;
    entry.img = ImageUInt8::create(numRows, numCols);
    std::copy(data, data + (numRows * numCols), entry.img->rawData());
    stbi_image_free(data);
    entry.cols = numCols;
    entry.rows = numRows;
    
    return entry;
}

static std::vector<long long> histogramAreas(ComponentTreeP* tree, int maxBin) {
    std::vector<long long> counts(static_cast<size_t>(maxBin), 0);
    for (NodeId id : tree->getIteratorBreadthFirstTraversalById()) {
        int area = tree->getAreaById(id);
        int bin = area;
        if (bin < 1) bin = 1;
        if (bin > maxBin) bin = maxBin;
        counts[static_cast<size_t>(bin - 1)]++;
    }
    return counts;
}

static void printJsonIntArray(const std::vector<int>& values) {
    std::cout << "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << values[i];
    }
    std::cout << "]";
}

static void printJsonLongArray(const std::vector<long long>& values) {
    std::cout << "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << values[i];
    }
    std::cout << "]";
}

static void printJsonDoubleArray(const std::vector<double>& values) {
    std::cout << "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << values[i];
    }
    std::cout << "]";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <images...>\n";
        std::cerr << "Example: " << argv[0] << " ./images/*.png\n";
        std::cerr << "Options:\n";
        std::cerr << "  --per-image   Emit JSON per image instead of dataset average\n";
        return 1;
    }

    bool perImage = false;
    std::vector<std::string> inputs;
    for (int i = 1; i < argc; ++i) {
        std::string value = trimCopy(argv[i]);
        if (value.empty()) continue;
        if (value == "--per-image") {
            perImage = true;
            continue;
        }
        if (!expandGlobPath(value, inputs)) {
            std::cerr << "Aviso: sem imagens para o padrao " << value << "\n";
        }
    }

    if (inputs.empty()) {
        std::cerr << "Erro: nenhuma imagem encontrada.\n";
        return 1;
    }

    int maxBin = 0;
    std::vector<long long> sumMax;
    std::vector<long long> sumMin;

    for (const auto& filename : inputs) {
        ImageEntry entry = loadImageEntry(filename);
        const int imageMaxBin = entry.rows * entry.cols;
        if (!perImage && imageMaxBin > maxBin) {
            maxBin = imageMaxBin;
            sumMax.resize(static_cast<size_t>(maxBin), 0);
            sumMin.resize(static_cast<size_t>(maxBin), 0);
        }

        AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(entry.rows, entry.cols, 1.5);
        ComponentTreePPtr maxTreePtr = std::make_shared<ComponentTreeP>(entry.img, true, adj);
        ComponentTreePPtr minTreePtr = std::make_shared<ComponentTreeP>(entry.img, false, adj);

        std::vector<long long> maxHist = histogramAreas(maxTreePtr.get(), imageMaxBin);
        std::vector<long long> minHist = histogramAreas(minTreePtr.get(), imageMaxBin);

        if (perImage) {
            std::vector<int> bins;
            bins.reserve(static_cast<size_t>(imageMaxBin));
            for (int i = 1; i <= imageMaxBin; ++i) {
                bins.push_back(i);
            }

            std::cout << "{\n";
            std::cout << "  \"image\":\"" << jsonEscape(entry.name) << "\",\n";
            std::cout << "  \"rows\":" << entry.rows << ",\n";
            std::cout << "  \"cols\":" << entry.cols << ",\n";
            std::cout << "  \"bins\":";
            printJsonIntArray(bins);
            std::cout << ",\n";
            std::cout << "  \"maxtree_area_hist\":{\"bins\":";
            printJsonIntArray(bins);
            std::cout << ",\"counts\":";
            printJsonLongArray(maxHist);
            std::cout << "},\n";
            std::cout << "  \"mintree_area_hist\":{\"bins\":";
            printJsonIntArray(bins);
            std::cout << ",\"counts\":";
            printJsonLongArray(minHist);
            std::cout << "}\n";
            std::cout << "}\n";
        } else {
            for (size_t i = 0; i < maxHist.size(); ++i) {
                sumMax[i] += maxHist[i];
                sumMin[i] += minHist[i];
            }
        }
    }

    if (perImage) {
        return 0;
    }

    if (maxBin <= 0) {
        std::cerr << "Erro: nenhuma imagem valida para histograma.\n";
        return 1;
    }

    std::vector<int> bins;
    bins.reserve(static_cast<size_t>(maxBin));
    for (int i = 1; i <= maxBin; ++i) {
        bins.push_back(i);
    }

    const double denom = static_cast<double>(inputs.size());
    std::vector<double> avgMax(sumMax.size(), 0.0);
    std::vector<double> avgMin(sumMin.size(), 0.0);
    for (size_t i = 0; i < sumMax.size(); ++i) {
        avgMax[i] = sumMax[i] / denom;
        avgMin[i] = sumMin[i] / denom;
    }

    std::cout << "{\n";
    std::cout << "  \"num_images\":" << inputs.size() << ",\n";
    std::cout << "  \"bins\":";
    printJsonIntArray(bins);
    std::cout << ",\n";
    std::cout << "  \"maxtree_area_hist\":{\"bins\":";
    printJsonIntArray(bins);
    std::cout << ",\"avg_counts\":";
    printJsonDoubleArray(avgMax);
    std::cout << "},\n";
    std::cout << "  \"mintree_area_hist\":{\"bins\":";
    printJsonIntArray(bins);
    std::cout << ",\"avg_counts\":";
    printJsonDoubleArray(avgMin);
    std::cout << "}\n";
    std::cout << "}\n";

    return 0;
}
