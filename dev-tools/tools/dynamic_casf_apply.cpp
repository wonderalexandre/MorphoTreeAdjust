#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "dynamic_casf_apply_common.hpp"

namespace common = dynamic_casf_apply_common;

enum class CasfMode {
    DynamicSubtree,
    DynamicLeaf,
    Naive,
    Compare
};

static void printUsage(const char* argv0) {
    std::cerr << "Usage: " << argv0
              << " [--mode dynamic-subtree|dynamic-leaf|naive|compare]"
              << " [--attribute area|bbox_width|bbox_height|bbox_diagonal]"
              << " [--radio-adj <radius>]"
              << " [--iter-timing]"
              << " [--no-output]"
              << " <input.png> [<output.png>] <threshold1> [threshold2 ...]\n";
}

int main(int argc, char** argv) {
    try {
        if (argc < 4) {
            printUsage(argv[0]);
            return 1;
        }

        int argi = 1;
        CasfMode mode = CasfMode::DynamicSubtree;
        Attribute attribute = AREA;
        double radioAdj = 1.0;
        bool iterTiming = false;
        bool noOutput = false;
        while (argi < argc && std::string(argv[argi]).rfind("--", 0) == 0) {
            if (std::string(argv[argi]) == "--mode") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                const std::string modeValue = argv[argi + 1];
                if (modeValue == "dynamic-subtree") {
                    mode = CasfMode::DynamicSubtree;
                } else if (modeValue == "dynamic-leaf") {
                    mode = CasfMode::DynamicLeaf;
                } else if (modeValue == "naive") {
                    mode = CasfMode::Naive;
                } else if (modeValue == "compare") {
                    mode = CasfMode::Compare;
                } else {
                    throw std::runtime_error("Invalid mode: " + modeValue);
                }
                argi += 2;
                continue;
            }
            if (std::string(argv[argi]) == "--attribute") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                attribute = common::parseAttributeName(argv[argi + 1]);
                argi += 2;
                continue;
            }
            if (std::string(argv[argi]) == "--radio-adj") {
                if (argi + 1 >= argc) {
                    printUsage(argv[0]);
                    return 1;
                }
                char* endPtr = nullptr;
                radioAdj = std::strtod(argv[argi + 1], &endPtr);
                if (endPtr == argv[argi + 1] || *endPtr != '\0' || radioAdj <= 0.0) {
                    throw std::runtime_error("Invalid --radio-adj value. Expected a positive radius.");
                }
                argi += 2;
                continue;
            }
            if (std::string(argv[argi]) == "--iter-timing") {
                iterTiming = true;
                ++argi;
                continue;
            }
            if (std::string(argv[argi]) == "--no-output") {
                noOutput = true;
                ++argi;
                continue;
            }
            throw std::runtime_error("Unknown option: " + std::string(argv[argi]));
        }

        if ((!noOutput && argc - argi < 3) || (noOutput && argc - argi < 2)) {
            printUsage(argv[0]);
            return 1;
        }

        const std::string inputPath = argv[argi++];
        const std::string outputPath = noOutput ? std::string() : argv[argi++];
        auto input = common::loadGrayImage(inputPath);
        std::vector<int> thresholds = common::parseThresholds(argc, argv, argi);
        std::cout << "Input: " << inputPath << "\n";
        if (noOutput) {
            std::cout << "Output: disabled\n";
        } else {
            std::cout << "Output: " << outputPath << "\n";
        }
        std::cout << "Image size: " << input->getNumRows() << "x" << input->getNumCols() << "\n";
        std::cout << "Mode: "
                  << (mode == CasfMode::DynamicSubtree ? "dynamic-subtree"
                      : mode == CasfMode::DynamicLeaf ? "dynamic-leaf"
                      : mode == CasfMode::Naive ? "naive"
                      : "compare")
                  << "\n";
        std::cout << "Attribute: " << common::attributeName(attribute) << "\n";
        std::cout << "Adjacency radius: " << radioAdj << "\n";

        ImageUInt8Ptr output;
        if (mode == CasfMode::DynamicSubtree) {
            common::DynamicComponentTreeCasfRunner runner(input, radioAdj, attribute);
            const auto &initStats = runner.getInitializationStats();
            const long long initUs = common::totalInitializationMicros(initStats);
            long long iterationUs = 0;
            for (int threshold : thresholds) {
                runner.applyThreshold(threshold);
                const auto &stats = runner.getIterationStats();
                iterationUs += stats.phase1_candidates_us +
                               stats.phase1_adjust_us +
                               stats.phase1_refresh_area_us +
                               stats.phase2_candidates_us +
                               stats.phase2_adjust_us +
                               stats.phase2_refresh_area_us;
                if (iterTiming) {
                    common::printDynamicIterationStats(stats, "runner");
                }
            }
            output = runner.reconstructImage();
            const long long reconstructUs = runner.getIterationStats().reconstruct_us;
            const long long totalUs = initUs + iterationUs + reconstructUs;
            std::cout << "Initialization time: " << common::microsToMillis(initUs) << " ms\n";
            std::cout << "Iteration time: " << common::microsToMillis(iterationUs) << " ms\n";
            std::cout << "Reconstruction time: " << common::microsToMillis(reconstructUs) << " ms\n";
            std::cout << "Total time: " << common::microsToMillis(totalUs) << " ms\n";
        } else if (mode == CasfMode::DynamicLeaf) {
            common::DynamicComponentTreeLeafCasfRunner runner(input, radioAdj, attribute);
            const auto &initStats = runner.getInitializationStats();
            const long long initUs = common::totalInitializationMicros(initStats);
            long long iterationUs = 0;
            for (int threshold : thresholds) {
                runner.applyThreshold(threshold);
                const auto &stats = runner.getIterationStats();
                iterationUs += stats.phase1_candidates_us +
                               stats.phase1_adjust_us +
                               stats.phase1_refresh_area_us +
                               stats.phase2_candidates_us +
                               stats.phase2_adjust_us +
                               stats.phase2_refresh_area_us;
                if (iterTiming) {
                    common::printDynamicLeafIterationStats(stats, "runner");
                }
            }
            output = runner.reconstructImage();
            const long long reconstructUs = runner.getIterationStats().reconstruct_us;
            const long long totalUs = initUs + iterationUs + reconstructUs;
            std::cout << "Initialization time: " << common::microsToMillis(initUs) << " ms\n";
            std::cout << "Iteration time: " << common::microsToMillis(iterationUs) << " ms\n";
            std::cout << "Reconstruction time: " << common::microsToMillis(reconstructUs) << " ms\n";
            std::cout << "Total time: " << common::microsToMillis(totalUs) << " ms\n";
        } else if (mode == CasfMode::Naive) {
            long long totalUs = 0;
            output = input->clone();
            for (int threshold : thresholds) {
                const auto iterStart = std::chrono::steady_clock::now();
                output = common::applyNaiveThreshold(output, radioAdj, threshold, attribute);
                const auto iterEnd = std::chrono::steady_clock::now();
                const long long iterUs = common::elapsedMicros(iterStart, iterEnd);
                totalUs += iterUs;
                if (iterTiming) {
                    std::cout << "  naive threshold=" << threshold
                              << ": total=" << common::microsToMillis(iterUs) << " ms\n";
                }
            }
            std::cout << "Total time: " << common::microsToMillis(totalUs) << " ms\n";
        } else if (mode == CasfMode::Compare) {
            struct CompareTimingLine {
                int threshold = 0;
                double dynamicMs = 0.0;
                double leafMs = 0.0;
                double naiveMs = 0.0;
            };

            common::DynamicComponentTreeCasfRunner dynamicRunner(input, radioAdj, attribute);
            common::DynamicComponentTreeLeafCasfRunner leafRunner(input, radioAdj, attribute);
            auto naiveOutput = input->clone();
            std::vector<CompareTimingLine> timingLines;
            timingLines.reserve(thresholds.size());
            const auto &dynamicInitStats = dynamicRunner.getInitializationStats();
            const auto &leafInitStats = leafRunner.getInitializationStats();
            const long long dynamicInitUs = common::totalInitializationMicros(dynamicInitStats);
            const long long leafInitUs = common::totalInitializationMicros(leafInitStats);
            long long dynamicIterationUs = 0;
            long long leafIterationUs = 0;
            long long naiveTotalUs = 0;

            for (int threshold : thresholds) {
                const auto naiveIterStart = std::chrono::steady_clock::now();
                naiveOutput = common::applyNaiveThreshold(naiveOutput, radioAdj, threshold, attribute);
                const auto naiveIterEnd = std::chrono::steady_clock::now();
                const long long naiveIterUs = common::elapsedMicros(naiveIterStart, naiveIterEnd);
                naiveTotalUs += naiveIterUs;

                dynamicRunner.applyThreshold(threshold);
                const auto &dynamicStats = dynamicRunner.getIterationStats();
                dynamicIterationUs += dynamicStats.phase1_candidates_us +
                                      dynamicStats.phase1_adjust_us +
                                      dynamicStats.phase1_refresh_area_us +
                                      dynamicStats.phase2_candidates_us +
                                      dynamicStats.phase2_adjust_us +
                                      dynamicStats.phase2_refresh_area_us;

                leafRunner.applyThreshold(threshold);
                const auto &leafStats = leafRunner.getIterationStats();
                leafIterationUs += leafStats.phase1_candidates_us +
                                   leafStats.phase1_adjust_us +
                                   leafStats.phase1_refresh_area_us +
                                   leafStats.phase2_candidates_us +
                                   leafStats.phase2_adjust_us +
                                   leafStats.phase2_refresh_area_us;

                timingLines.push_back({
                    .threshold = threshold,
                    .dynamicMs = common::totalIterationMillis(dynamicStats),
                    .leafMs = common::totalIterationMillis(leafStats),
                    .naiveMs = common::microsToMillis(naiveIterUs),
                });
            }

            auto dynamicOutput = dynamicRunner.reconstructImage();
            auto leafOutput = leafRunner.reconstructImage();
            const long long dynamicReconstructUs = dynamicRunner.getIterationStats().reconstruct_us;
            const long long leafReconstructUs = leafRunner.getIterationStats().reconstruct_us;
            const long long dynamicTotalUs = dynamicInitUs + dynamicIterationUs + dynamicReconstructUs;
            const long long leafTotalUs = leafInitUs + leafIterationUs + leafReconstructUs;

            if (iterTiming) {
                for (const auto &line : timingLines) {
                    std::cout << "  threshold=" << line.threshold
                              << " | dynamic-subtree=" << line.dynamicMs << " ms"
                              << " | dynamic-leaf=" << line.leafMs << " ms"
                              << " | naive=" << line.naiveMs << " ms\n";
                }
            }

            const bool dynamicMatchesNaive = dynamicOutput->isEqual(naiveOutput);
            const bool leafMatchesNaive = leafOutput->isEqual(naiveOutput);
            const bool dynamicMatchesLeaf = dynamicOutput->isEqual(leafOutput);

            std::cout << "\nCompare summary\n";
            std::cout << "  dynamic-subtree init=" << common::microsToMillis(dynamicInitUs)
                      << " ms iter=" << common::microsToMillis(dynamicIterationUs)
                      << " ms reconstruct=" << common::microsToMillis(dynamicReconstructUs)
                      << " ms total=" << common::microsToMillis(dynamicTotalUs) << " ms\n";
            std::cout << "  dynamic-leaf init=" << common::microsToMillis(leafInitUs)
                      << " ms iter=" << common::microsToMillis(leafIterationUs)
                      << " ms reconstruct=" << common::microsToMillis(leafReconstructUs)
                      << " ms total=" << common::microsToMillis(leafTotalUs) << " ms\n";
            std::cout << "  naive total=" << common::microsToMillis(naiveTotalUs) << " ms\n";
            std::cout << "  dynamic-subtree == naive: " << (dynamicMatchesNaive ? "true" : "false") << "\n";
            std::cout << "  dynamic-leaf == naive: " << (leafMatchesNaive ? "true" : "false") << "\n";
            std::cout << "  dynamic-subtree == dynamic-leaf: " << (dynamicMatchesLeaf ? "true" : "false") << "\n";

            const bool hasRegression = !dynamicMatchesNaive ||
                                       !leafMatchesNaive ||
                                       !dynamicMatchesLeaf;

            if (hasRegression) {
                throw std::runtime_error("compare mode found divergent outputs between implementations");
            }

            output = dynamicOutput;
        }

        if (!noOutput) {
            common::saveGrayImage(outputPath, output);
            std::cout << "\nSaved filtered image to " << outputPath << "\n";
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
