#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "dynamic_casf_apply_common.hpp"

#include "../include/morphoTreeAdjustResearch/AdjacencyRelation.hpp"
#include "../include/morphoTreeAdjustResearch/AttributeComputer.hpp"
#include "../include/morphoTreeAdjustResearch/Common.hpp"
#include "../include/morphoTreeAdjustResearch/DynamicComponentTree.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseAreaAttributeUpdater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseEligibility.hpp"
#include "../include/morphoTreeAdjustResearch/debug/Algorithm2FlatRegionMergeUpdater.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChange.hpp"
#include "../include/morphoTreeAdjustResearch/debug/RegionLevelChangeRebuildUpdater.hpp"

namespace common = dynamic_casf_apply_common;
namespace problem2_debug = debug::problem2;
namespace region_debug = debug::region;

namespace {

enum class QaCaseFilter {
    Any,
    Preserve,
    Release,
};

struct Options {
    double radioAdj = 1.5;
    int maxEventsPerImage = 20;
    int maxEventsPerNode = 1;
    int minRegionSize = 1;
    int maxRegionSize = 1;
    QaCaseFilter qaCaseFilter = QaCaseFilter::Any;
    bool requireHigherNeighbor = true;
    int requireAcceptedEvents = 0;
    bool failOnRejected = false;
    bool failOnMismatch = false;
};

struct Candidate {
    PixelId pixel = InvalidNode;
    std::vector<PixelId> pixelsP;
    NodeId sourceNodeQa = InvalidNode;
    int row = 0;
    int col = 0;
    int a = 0;
    int b = 0;
    int regionSize = 0;
    int directProperParts = 0;
    int children = 0;
    int higherNeighbors = 0;
    bool candidateReleasesQa = false;
};

struct RunResult {
    long long localStructUs = 0;
    long long localAreaUs = 0;
    long long rebuildStructUs = 0;
    long long rebuildAreaUs = 0;
    int affectedLevels = 0;
    int createdNodes = 0;
    int absorbedNodes = 0;
    int xiBulletNodes = 0;
    int xiCircleNodes = 0;
    int areaDirtyClosure = 0;
    int areaRecomputedNodes = 0;
    bool preservedQa = false;
    bool releasedQa = false;
    bool reconstructionMatches = false;
    bool areaMatches = false;
};

struct ImageSummary {
    int scannedCandidates = 0;
    int acceptedEvents = 0;
    int rejectedDuringUpdate = 0;
    int mismatchedEvents = 0;
};

void printUsage(const char *argv0) {
    std::cerr << "Usage: " << argv0
              << " [--radio-adj <radius>]"
              << " [--max-events-per-image <n>]"
              << " [--max-events-per-node <n>]"
              << " [--min-region-size <n>]"
              << " [--max-region-size <n>]"
              << " [--qa-case any|preserve|release]"
              << " [--allow-no-higher-neighbor]"
              << " [--require-events <n>]"
              << " [--fail-on-rejected]"
              << " [--fail-on-mismatch]"
              << " <image1> [image2 ...]\n";
}

QaCaseFilter parseQaCaseFilter(const std::string &value) {
    if (value == "any") {
        return QaCaseFilter::Any;
    }
    if (value == "preserve" || value == "preserved") {
        return QaCaseFilter::Preserve;
    }
    if (value == "release" || value == "released") {
        return QaCaseFilter::Release;
    }
    throw std::runtime_error("Invalid --qa-case value: " + value);
}

std::string basename(const std::string &path) {
    return std::filesystem::path(path).filename().string();
}

long long elapsedMicros(std::chrono::steady_clock::time_point start,
                        std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

double microsToMillis(long long micros) {
    return static_cast<double>(micros) / 1000.0;
}

std::vector<float> computeArea(DynamicComponentTree &tree) {
    DynamicAreaComputer computer(&tree);
    return computer.compute();
}

bool areaBufferMatchesFullRecompute(DynamicComponentTree &tree,
                                    const std::vector<float> &areaBuffer) {
    DynamicAreaComputer computer(&tree);
    const auto full = computer.compute();
    for (NodeId nodeId : tree.getNodeSubtree(tree.getRoot())) {
        if (!tree.isAlive(nodeId)) {
            continue;
        }
        const auto index = static_cast<std::size_t>(nodeId);
        if (index >= areaBuffer.size() || index >= full.size()) {
            return false;
        }
        if (areaBuffer[index] != full[index]) {
            return false;
        }
    }
    return true;
}

std::vector<PixelId> growConnectedRegionFromDirectProperParts(
    AdjacencyRelation &adjacency,
    PixelId seed,
    const std::vector<int> &directProperPartStamp,
    int directProperPartMark,
    std::vector<int> &visitedStamp,
    int visitedMark,
    int maxRegionSize) {
    std::vector<PixelId> region;
    std::vector<PixelId> stack;
    region.reserve(static_cast<std::size_t>(maxRegionSize));
    stack.push_back(seed);
    visitedStamp[static_cast<std::size_t>(seed)] = visitedMark;

    while (!stack.empty()) {
        const PixelId pixel = stack.back();
        stack.pop_back();
        if (static_cast<int>(region.size()) < maxRegionSize) {
            region.push_back(pixel);
        }

        for (PixelId neighbor : adjacency.getNeighborPixels(pixel)) {
            const auto index = static_cast<std::size_t>(neighbor);
            if (directProperPartStamp[index] != directProperPartMark ||
                visitedStamp[index] == visitedMark) {
                continue;
            }
            visitedStamp[index] = visitedMark;
            stack.push_back(neighbor);
        }

        if (static_cast<int>(region.size()) >= maxRegionSize) {
            break;
        }
    }

    return region;
}

int chooseTargetLevelFromRegionBoundary(const ImageUInt8Ptr &image,
                                        AdjacencyRelation &adjacency,
                                        std::span<const PixelId> pixelsP,
                                        int a,
                                        std::vector<int> &regionStamp,
                                        int regionMark,
                                        int &higherNeighbors) {
    int b = 256;
    higherNeighbors = 0;

    for (PixelId pixel : pixelsP) {
        regionStamp[static_cast<std::size_t>(pixel)] = regionMark;
    }

    for (PixelId pixel : pixelsP) {
        for (PixelId neighbor : adjacency.getNeighborPixels(pixel)) {
            if (regionStamp[static_cast<std::size_t>(neighbor)] == regionMark) {
                continue;
            }
            const int value = static_cast<int>((*image)[neighbor]);
            if (value > a) {
                ++higherNeighbors;
                b = std::min(b, value);
            }
        }
    }

    if (b == 256) {
        return std::min(255, a + 1);
    }
    return b;
}

std::vector<Candidate> collectFlatRaiseCandidates(const ImageUInt8Ptr &image,
                                                  DynamicComponentTree &maxTree,
                                                  AdjacencyRelation &adjacency,
                                                  const Options &options) {
    std::vector<Candidate> candidates;
    if (maxTree.getRoot() == InvalidNode) {
        return candidates;
    }

    const auto domainSize = static_cast<std::size_t>(image->getSize());
    std::vector<int> directProperPartStamp(domainSize, 0);
    std::vector<int> visitedStamp(domainSize, 0);
    std::vector<int> regionStamp(domainSize, 0);
    int directProperPartMark = 0;
    int visitedMark = 0;
    int regionMark = 0;

    for (NodeId nodeId : maxTree.getNodeSubtree(maxTree.getRoot())) {
        if (!maxTree.isAlive(nodeId) || maxTree.getNumProperParts(nodeId) == 0) {
            continue;
        }

        const int a = maxTree.getAltitude(nodeId);
        if (a >= 255) {
            continue;
        }

        ++directProperPartMark;
        std::vector<PixelId> directProperParts;
        directProperParts.reserve(static_cast<std::size_t>(maxTree.getNumProperParts(nodeId)));
        for (PixelId pixel : maxTree.getProperParts(nodeId)) {
            directProperParts.push_back(pixel);
            directProperPartStamp[static_cast<std::size_t>(pixel)] = directProperPartMark;
        }

        int emittedForNode = 0;
        for (PixelId pixel : directProperParts) {
            if (static_cast<int>((*image)[pixel]) != a) {
                continue;
            }

            ++visitedMark;
            auto pixelsP = growConnectedRegionFromDirectProperParts(
                adjacency,
                pixel,
                directProperPartStamp,
                directProperPartMark,
                visitedStamp,
                visitedMark,
                options.maxRegionSize);
            if (pixelsP.empty()) {
                continue;
            }
            if (static_cast<int>(pixelsP.size()) < options.minRegionSize) {
                continue;
            }

            const int directProperParts = maxTree.getNumProperParts(nodeId);
            const bool candidateReleasesQa =
                static_cast<int>(pixelsP.size()) == directProperParts;
            if (options.qaCaseFilter == QaCaseFilter::Release && !candidateReleasesQa) {
                continue;
            }
            if (options.qaCaseFilter == QaCaseFilter::Preserve && candidateReleasesQa) {
                continue;
            }

            int higherNeighbors = 0;
            ++regionMark;
            const int b = chooseTargetLevelFromRegionBoundary(
                image,
                adjacency,
                std::span<const PixelId>(pixelsP.data(), pixelsP.size()),
                a,
                regionStamp,
                regionMark,
                higherNeighbors);
            if (b <= a) {
                continue;
            }
            if (options.requireHigherNeighbor && higherNeighbors == 0) {
                continue;
            }

            auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
                std::span<const PixelId>(pixelsP.data(), pixelsP.size()),
                static_cast<AltitudeType>(a),
                static_cast<AltitudeType>(b),
                true,
                true);
            const auto eligibility = problem2_debug::checkAlgorithm2FlatRaiseEligibility(
                maxTree,
                adjacency,
                change);
            if (!eligibility.accepted) {
                continue;
            }

            Candidate candidate;
            candidate.pixel = pixel;
            candidate.pixelsP = std::move(pixelsP);
            candidate.sourceNodeQa = eligibility.sourceNodeQa;
            candidate.row = pixel / image->getNumCols();
            candidate.col = pixel % image->getNumCols();
            candidate.a = eligibility.a;
            candidate.b = eligibility.b;
            candidate.regionSize = eligibility.pixelsP;
            candidate.directProperParts = eligibility.directProperPartsOfQa;
            candidate.children = eligibility.childrenOfQa;
            candidate.higherNeighbors = eligibility.higherBoundaryNeighbors;
            candidate.candidateReleasesQa = eligibility.releasesQa;
            candidates.push_back(candidate);

            ++emittedForNode;
            if (emittedForNode >= options.maxEventsPerNode ||
                static_cast<int>(candidates.size()) >= options.maxEventsPerImage) {
                break;
            }
        }

        if (static_cast<int>(candidates.size()) >= options.maxEventsPerImage) {
            break;
        }
    }

    return candidates;
}

RunResult runCandidate(const ImageUInt8Ptr &image,
                       const AdjacencyRelationPtr &adjacency,
                       const DynamicComponentTree &baseMaxTree,
	                       const std::vector<float> &baseArea,
	                       const Candidate &candidate) {
	    RunResult out;
	    const auto change = region_debug::buildExplicitRegionChange<AltitudeType>(
	        std::span<const PixelId>(candidate.pixelsP.data(), candidate.pixelsP.size()),
	        static_cast<AltitudeType>(candidate.a),
	        static_cast<AltitudeType>(candidate.b),
	        false,
        true);

    DynamicComponentTree localMaxTree = baseMaxTree;
    std::vector<float> localArea = baseArea;

    problem2_debug::Algorithm2FlatRaiseApplyResult structuralResult;
    auto start = std::chrono::steady_clock::now();
    structuralResult =
        problem2_debug::Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner(
            localMaxTree,
            *adjacency,
            change);
    out.localStructUs = elapsedMicros(start, std::chrono::steady_clock::now());

    start = std::chrono::steady_clock::now();
    const auto areaResult = problem2_debug::Algorithm2FlatRaiseAreaAttributeUpdater::recomputeArea(
        localMaxTree,
        structuralResult,
        localArea);
    out.localAreaUs = elapsedMicros(start, std::chrono::steady_clock::now());

    const auto expectedImage = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(
        image,
        change);
    out.reconstructionMatches =
        structuralResult.reconstruction->isEqual(expectedImage) &&
        localMaxTree.reconstructionImage()->isEqual(expectedImage);
    out.areaMatches = areaBufferMatchesFullRecompute(localMaxTree, localArea);

    start = std::chrono::steady_clock::now();
    const auto rebuiltImage = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(
        image,
        change);
    DynamicComponentTree rebuiltMaxTree(rebuiltImage->clone(), true, adjacency);
    const auto rebuiltReconstruction = rebuiltMaxTree.reconstructionImage();
    if (!rebuiltReconstruction->isEqual(expectedImage)) {
        throw std::runtime_error("rebuild oracle reconstruction mismatch");
    }
    out.rebuildStructUs = elapsedMicros(start, std::chrono::steady_clock::now());

    start = std::chrono::steady_clock::now();
    const auto rebuiltAreaImage = region_debug::RegionLevelChangeRebuildUpdater::applyToImage(
        image,
        change);
    DynamicComponentTree rebuiltAreaMaxTree(rebuiltAreaImage->clone(), true, adjacency);
    DynamicAreaComputer rebuiltAreaComputer(&rebuiltAreaMaxTree);
    auto rebuiltArea = rebuiltAreaComputer.compute();
    static_cast<void>(rebuiltArea);
    out.rebuildAreaUs = elapsedMicros(start, std::chrono::steady_clock::now());

    out.affectedLevels = static_cast<int>(structuralResult.affectedLevels.size());
    out.createdNodes = static_cast<int>(structuralResult.createdNodes.size());
    out.absorbedNodes = static_cast<int>(structuralResult.absorbedNodes.size());
    out.xiBulletNodes = static_cast<int>(structuralResult.xiBulletNodes.size());
    out.xiCircleNodes = static_cast<int>(structuralResult.xiCircleNodes.size());
    out.areaDirtyClosure = static_cast<int>(areaResult.dirtyClosure.size());
    out.areaRecomputedNodes = static_cast<int>(areaResult.recomputedNodes.size());
    out.preservedQa = structuralResult.preservedSourceNodeQa;
    out.releasedQa = structuralResult.releasedSourceNodeQa;
    return out;
}

	void printCsvHeader() {
	    std::cout << "image,rows,cols,event_index,pixels_p,pixel,row,col,source_node,a,b,"
	              << "qa_direct_proper_parts,qa_children,higher_neighbors,candidate_releases_qa,"
	              << "local_struct_ms,local_area_ms,local_total_ms,"
              << "rebuild_struct_ms,rebuild_area_ms,"
              << "struct_speedup_vs_rebuild,area_speedup_vs_rebuild,"
              << "affected_levels,created_nodes,absorbed_nodes,xi_bullet_nodes,xi_circle_nodes,"
              << "area_dirty_closure,area_recomputed_nodes,preserved_qa,released_qa,"
              << "reconstruction_matches,area_matches\n";
}

void printCsvRow(const std::string &imagePath,
                 const ImageUInt8Ptr &image,
                 int eventIndex,
                 const Candidate &candidate,
                 const RunResult &result) {
    const double localStructMs = microsToMillis(result.localStructUs);
    const double localAreaMs = microsToMillis(result.localAreaUs);
    const double localTotalMs = microsToMillis(result.localStructUs + result.localAreaUs);
    const double rebuildStructMs = microsToMillis(result.rebuildStructUs);
    const double rebuildAreaMs = microsToMillis(result.rebuildAreaUs);
    const double structSpeedup =
        result.localStructUs > 0 ? static_cast<double>(result.rebuildStructUs) / result.localStructUs : 0.0;
    const double areaSpeedup =
        (result.localStructUs + result.localAreaUs) > 0
            ? static_cast<double>(result.rebuildAreaUs) /
                  static_cast<double>(result.localStructUs + result.localAreaUs)
            : 0.0;

    std::cout << basename(imagePath) << ","
	              << image->getNumRows() << ","
	              << image->getNumCols() << ","
	              << eventIndex << ","
	              << candidate.regionSize << ","
	              << candidate.pixel << ","
	              << candidate.row << ","
              << candidate.col << ","
              << candidate.sourceNodeQa << ","
              << candidate.a << ","
              << candidate.b << ","
              << candidate.directProperParts << ","
              << candidate.children << ","
              << candidate.higherNeighbors << ","
              << (candidate.candidateReleasesQa ? "true" : "false") << ","
              << std::fixed << std::setprecision(3)
              << localStructMs << ","
              << localAreaMs << ","
              << localTotalMs << ","
              << rebuildStructMs << ","
              << rebuildAreaMs << ","
              << structSpeedup << ","
              << areaSpeedup << ","
              << result.affectedLevels << ","
              << result.createdNodes << ","
              << result.absorbedNodes << ","
              << result.xiBulletNodes << ","
              << result.xiCircleNodes << ","
              << result.areaDirtyClosure << ","
              << result.areaRecomputedNodes << ","
              << (result.preservedQa ? "true" : "false") << ","
              << (result.releasedQa ? "true" : "false") << ","
              << (result.reconstructionMatches ? "true" : "false") << ","
              << (result.areaMatches ? "true" : "false") << "\n";
}

Options parseOptions(int argc, char **argv, int &firstImageArg) {
    Options options;
    int argi = 1;
    while (argi < argc && std::string(argv[argi]).rfind("--", 0) == 0) {
        const std::string arg = argv[argi];
        if (arg == "--radio-adj") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--radio-adj requires a value");
            }
            options.radioAdj = std::stod(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--max-events-per-image") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--max-events-per-image requires a value");
            }
            options.maxEventsPerImage = std::stoi(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--max-events-per-node") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--max-events-per-node requires a value");
            }
            options.maxEventsPerNode = std::stoi(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--min-region-size") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--min-region-size requires a value");
            }
            options.minRegionSize = std::stoi(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--max-region-size") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--max-region-size requires a value");
            }
            options.maxRegionSize = std::stoi(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--qa-case") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--qa-case requires a value");
            }
            options.qaCaseFilter = parseQaCaseFilter(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--allow-no-higher-neighbor") {
            options.requireHigherNeighbor = false;
            ++argi;
            continue;
        }
        if (arg == "--require-events") {
            if (argi + 1 >= argc) {
                throw std::runtime_error("--require-events requires a value");
            }
            options.requireAcceptedEvents = std::stoi(argv[argi + 1]);
            argi += 2;
            continue;
        }
        if (arg == "--fail-on-rejected") {
            options.failOnRejected = true;
            ++argi;
            continue;
        }
        if (arg == "--fail-on-mismatch") {
            options.failOnMismatch = true;
            ++argi;
            continue;
        }
        throw std::runtime_error("Unknown option: " + arg);
    }

    if (options.radioAdj <= 0.0) {
        throw std::runtime_error("--radio-adj must be positive");
    }
    if (options.maxEventsPerImage <= 0) {
        throw std::runtime_error("--max-events-per-image must be positive");
    }
    if (options.maxEventsPerNode <= 0) {
        throw std::runtime_error("--max-events-per-node must be positive");
    }
    if (options.minRegionSize <= 0) {
        throw std::runtime_error("--min-region-size must be positive");
    }
    if (options.maxRegionSize <= 0) {
        throw std::runtime_error("--max-region-size must be positive");
    }
    if (options.minRegionSize > options.maxRegionSize) {
        throw std::runtime_error("--min-region-size cannot be larger than --max-region-size");
    }
    if (options.requireAcceptedEvents < 0) {
        throw std::runtime_error("--require-events cannot be negative");
    }

    firstImageArg = argi;
    return options;
}

ImageSummary processImage(const std::string &imagePath,
                          const Options &options) {
    ImageSummary summary;
    const auto image = common::loadGrayImage(imagePath);
    auto adjacency = std::make_shared<AdjacencyRelation>(
        image->getNumRows(),
        image->getNumCols(),
        options.radioAdj);
    DynamicComponentTree baseMaxTree(image->clone(), true, adjacency);
    const auto baseArea = computeArea(baseMaxTree);

    const auto candidates = collectFlatRaiseCandidates(
        image,
        baseMaxTree,
        *adjacency,
        options);
    summary.scannedCandidates = static_cast<int>(candidates.size());

    int eventIndex = 0;
    for (const Candidate &candidate : candidates) {
        try {
            const auto result = runCandidate(
                image,
                adjacency,
                baseMaxTree,
                baseArea,
                candidate);
            printCsvRow(imagePath, image, eventIndex, candidate, result);
            ++summary.acceptedEvents;
            if (!result.reconstructionMatches || !result.areaMatches) {
                ++summary.mismatchedEvents;
            }
        } catch (const std::exception &e) {
            ++summary.rejectedDuringUpdate;
            std::cerr << "algorithm2_flat_raise_real_events: rejected "
                      << basename(imagePath)
                      << " event=" << eventIndex
                      << " pixel=" << candidate.pixel
                      << " reason=" << e.what() << "\n";
        }
        ++eventIndex;
    }

    return summary;
}

} // namespace

int main(int argc, char **argv) {
    try {
        int firstImageArg = 1;
        const Options options = parseOptions(argc, argv, firstImageArg);
        if (firstImageArg >= argc) {
            printUsage(argv[0]);
            return 1;
        }

        printCsvHeader();

        ImageSummary total;
        int numImages = 0;
        for (int argi = firstImageArg; argi < argc; ++argi) {
            const auto summary = processImage(argv[argi], options);
            total.scannedCandidates += summary.scannedCandidates;
            total.acceptedEvents += summary.acceptedEvents;
            total.rejectedDuringUpdate += summary.rejectedDuringUpdate;
            total.mismatchedEvents += summary.mismatchedEvents;
            ++numImages;
        }

        std::cerr << "algorithm2_flat_raise_real_events: images=" << numImages
                  << " scanned_candidates=" << total.scannedCandidates
                  << " accepted_events=" << total.acceptedEvents
                  << " rejected_during_update=" << total.rejectedDuringUpdate
                  << " mismatched_events=" << total.mismatchedEvents
                  << "\n";
        if (options.requireAcceptedEvents > 0 &&
            total.acceptedEvents < options.requireAcceptedEvents) {
            std::cerr << "algorithm2_flat_raise_real_events: expected at least "
                      << options.requireAcceptedEvents
                      << " accepted events, got "
                      << total.acceptedEvents
                      << "\n";
            return 2;
        }
        if (options.failOnRejected && total.rejectedDuringUpdate > 0) {
            std::cerr << "algorithm2_flat_raise_real_events: rejected events found in regression mode\n";
            return 3;
        }
        if (options.failOnMismatch && total.mismatchedEvents > 0) {
            std::cerr << "algorithm2_flat_raise_real_events: mismatched events found in regression mode\n";
            return 4;
        }
    } catch (const std::exception &e) {
        std::cerr << "algorithm2_flat_raise_real_events: " << e.what() << "\n";
        printUsage(argv[0]);
        return 1;
    }

    return 0;
}
