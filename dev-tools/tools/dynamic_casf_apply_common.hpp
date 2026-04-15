#pragma once

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../../morphoTreeAdjust/include/AttributeComputer.hpp"
#include "../../morphoTreeAdjust/include/Common.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTree.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTreeAdjustment.hpp"
#include "../../morphoTreeAdjust/include/DynamicComponentTreeAdjustmentLeaf.hpp"

#include "../external/stb/stb_image.h"
#include "../external/stb/stb_image_write.h"

namespace dynamic_casf_apply_common {

inline const char *attributeName(Attribute attribute) {
    switch (attribute) {
        case AREA: return "area";
        case BOX_WIDTH: return "bbox_width";
        case BOX_HEIGHT: return "bbox_height";
        case DIAGONAL_LENGTH: return "bbox_diagonal";
    }
    return "area";
}

inline Attribute parseAttributeName(const std::string &value) {
    if (value == "area") {
        return AREA;
    }
    if (value == "bbox_width" || value == "box_width") {
        return BOX_WIDTH;
    }
    if (value == "bbox_height" || value == "box_height") {
        return BOX_HEIGHT;
    }
    if (value == "bbox_diagonal" || value == "box_diagonal" || value == "bbox") {
        return DIAGONAL_LENGTH;
    }
    throw std::runtime_error("Invalid attribute: " + value);
}

inline std::unique_ptr<DynamicAttributeComputer> makeAttributeComputer(DynamicComponentTree *tree,
                                                                       Attribute attribute) {
    if (attribute == AREA) {
        return std::make_unique<DynamicAreaComputer>(tree);
    }
    return std::make_unique<DynamicBoundingBoxComputer>(tree, attribute);
}

inline std::vector<float> computeAttributeVector(DynamicComponentTree *tree, Attribute attribute) {
    auto computer = makeAttributeComputer(tree, attribute);
    std::vector<float> buffer(static_cast<std::size_t>(tree->getGlobalIdSpaceSize()), 0.0f);
    if (attribute == AREA) {
        static_cast<DynamicAreaComputer *>(computer.get())->compute(std::span<float>(buffer));
    } else {
        static_cast<DynamicBoundingBoxComputer *>(computer.get())->compute(std::span<float>(buffer));
    }
    return buffer;
}

inline std::vector<int> parseThresholds(int argc, char **argv, int startIndex) {
    std::vector<int> thresholds;
    for (int i = startIndex; i < argc; ++i) {
        char *endPtr = nullptr;
        const long value = std::strtol(argv[i], &endPtr, 10);
        if (endPtr == argv[i] || *endPtr != '\0' || value < 0) {
            throw std::runtime_error(std::string("Invalid threshold: ") + argv[i]);
        }
        thresholds.push_back(static_cast<int>(value));
    }
    return thresholds;
}

inline ImageUInt8Ptr loadGrayImage(const std::string &path) {
    int numCols = 0;
    int numRows = 0;
    int numChannels = 0;
    uint8_t *data = stbi_load(path.c_str(), &numCols, &numRows, &numChannels, 1);
    if (!data) {
        throw std::runtime_error("Failed to load image: " + path);
    }

    ImageUInt8Ptr img = ImageUInt8::create(numRows, numCols);
    std::copy(data, data + (numRows * numCols), img->rawData());
    stbi_image_free(data);
    return img;
}

inline void saveGrayImage(const std::string &path, ImageUInt8Ptr image) {
    if (!image) {
        throw std::runtime_error("Cannot save a null image.");
    }
    const int stride = image->getNumCols();
    if (stbi_write_png(path.c_str(), image->getNumCols(), image->getNumRows(), 1, image->rawData(), stride) == 0) {
        throw std::runtime_error("Failed to save image: " + path);
    }
}

inline long long elapsedMicros(std::chrono::steady_clock::time_point start,
                               std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

inline double microsToMillis(long long micros) {
    return static_cast<double>(micros) / 1000.0;
}

template <typename Factory>
inline auto timedConstruct(long long &elapsedUs, Factory &&factory) {
    const auto start = std::chrono::steady_clock::now();
    auto value = factory();
    elapsedUs = elapsedMicros(start, std::chrono::steady_clock::now());
    return value;
}

template <typename Callback>
inline void timedCall(long long &elapsedUs, Callback &&callback) {
    const auto start = std::chrono::steady_clock::now();
    callback();
    elapsedUs = elapsedMicros(start, std::chrono::steady_clock::now());
}

inline std::vector<NodeId> getNodesThreshold(DynamicComponentTree *tree,
                                             const std::vector<float> &attribute,
                                             int threshold) {
    std::vector<NodeId> out;
    FastQueue<NodeId> queue;
    queue.push(tree->getRoot());

    while (!queue.empty()) {
        const NodeId nodeId = queue.pop();
        if (attribute[nodeId] > threshold) {
            for (NodeId childId : tree->getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }

    return out;
}

inline ImageUInt8Ptr applyNaiveThreshold(ImageUInt8Ptr input,
                                         double radioAdj,
                                         int threshold,
                                         Attribute attribute) {
    AdjacencyRelationPtr adj = std::make_shared<AdjacencyRelation>(
        input->getNumRows(), input->getNumCols(), radioAdj);

    auto maxTreePtr = std::make_shared<DynamicComponentTree>(input, true, adj);
    DynamicComponentTree *maxTree = maxTreePtr.get();
    std::vector<float> maxArea = computeAttributeVector(maxTree, attribute);
    for (NodeId nodeId : getNodesThreshold(maxTree, maxArea, threshold)) {
        maxTree->pruneNode(nodeId);
    }
    ImageUInt8Ptr current = maxTree->reconstructionImage();

    auto minTreePtr = std::make_shared<DynamicComponentTree>(current, false, adj);
    DynamicComponentTree *minTree = minTreePtr.get();
    std::vector<float> minArea = computeAttributeVector(minTree, attribute);
    for (NodeId nodeId : getNodesThreshold(minTree, minArea, threshold)) {
        minTree->pruneNode(nodeId);
    }
    return minTree->reconstructionImage();
}

inline std::vector<int> getDynamicNodesThreshold(const DynamicComponentTree &tree,
                                                 const std::vector<float> &area,
                                                 int threshold) {
    std::vector<int> out;
    FastQueue<int> queue;
    queue.push(tree.getRoot());
    while (!queue.empty()) {
        const int nodeId = queue.pop();
        if (area[nodeId] > threshold) {
            for (int childId : tree.getChildren(nodeId)) {
                queue.push(childId);
            }
        } else {
            out.push_back(nodeId);
        }
    }
    return out;
}

class DynamicComponentTreeCasfRunner {
public:
    struct InitializationStats {
        long long input_clone_us = 0;
        long long adjacency_us = 0;
        long long max_tree_us = 0;
        long long min_tree_us = 0;
        long long max_attribute_computer_us = 0;
        long long min_attribute_computer_us = 0;
        long long max_attribute_values_us = 0;
        long long min_attribute_values_us = 0;
        long long attribute_binding_us = 0;
    };

    struct IterationStats {
        long long phase1_candidates_us = 0;
        long long phase1_adjust_us = 0;
        long long phase1_refresh_area_us = 0;
        long long phase2_candidates_us = 0;
        long long phase2_adjust_us = 0;
        long long phase2_refresh_area_us = 0;
        long long reconstruct_us = 0;

        void reset() {
            phase1_candidates_us = 0;
            phase1_adjust_us = 0;
            phase1_refresh_area_us = 0;
            phase2_candidates_us = 0;
            phase2_adjust_us = 0;
            phase2_refresh_area_us = 0;
            reconstruct_us = 0;
        }
    };

private:
    InitializationStats initStats_{};
    ImageUInt8Ptr inputImage_;
    AdjacencyRelationPtr adj_;
    DynamicComponentTree maxTree_;
    DynamicComponentTree minTree_;
    DynamicComponentTreeAdjustment<AltitudeType> adjust_;
    std::unique_ptr<DynamicAttributeComputer> maxAttributeComputer_;
    std::unique_ptr<DynamicAttributeComputer> minAttributeComputer_;
    Attribute attribute_;
    IterationStats iterStats_{};
    std::vector<float> maxArea_;
    std::vector<float> minArea_;

public:
    explicit DynamicComponentTreeCasfRunner(ImageUInt8Ptr image, double radioAdj, Attribute attribute)
        : inputImage_(timedConstruct(initStats_.input_clone_us, [&]() { return image->clone(); })),
          adj_(timedConstruct(initStats_.adjacency_us, [&]() {
              return std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radioAdj);
          })),
          maxTree_(timedConstruct(initStats_.max_tree_us, [&]() { return DynamicComponentTree(inputImage_, true, adj_); })),
          minTree_(timedConstruct(initStats_.min_tree_us, [&]() { return DynamicComponentTree(inputImage_, false, adj_); })),
          adjust_(&minTree_, &maxTree_, *adj_),
          maxAttributeComputer_(timedConstruct(initStats_.max_attribute_computer_us, [&]() {
              return makeAttributeComputer(&maxTree_, attribute);
          })),
          minAttributeComputer_(timedConstruct(initStats_.min_attribute_computer_us, [&]() {
              return makeAttributeComputer(&minTree_, attribute);
          })),
          attribute_(attribute),
          maxArea_(maxTree_.getGlobalIdSpaceSize(), 0.0f),
          minArea_(minTree_.getGlobalIdSpaceSize(), 0.0f) {
        timedCall(initStats_.max_attribute_values_us, [&]() {
            maxArea_ = computeAttributeVector(&maxTree_, attribute_);
        });
        timedCall(initStats_.min_attribute_values_us, [&]() {
            minArea_ = computeAttributeVector(&minTree_, attribute_);
        });
        timedCall(initStats_.attribute_binding_us, [&]() {
            adjust_.setAttributeComputer(*minAttributeComputer_,
                                         *maxAttributeComputer_,
                                         std::span<float>(minArea_),
                                         std::span<float>(maxArea_));
        });
    }

    void applyThreshold(int threshold) {
        iterStats_.reset();
        auto phaseStart = std::chrono::steady_clock::now();
        auto nodesToPrune = getDynamicNodesThreshold(maxTree_, maxArea_, threshold);
        iterStats_.phase1_candidates_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());

        phaseStart = std::chrono::steady_clock::now();
        adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
        iterStats_.phase1_adjust_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());
        iterStats_.phase1_refresh_area_us = 0;

        phaseStart = std::chrono::steady_clock::now();
        nodesToPrune = getDynamicNodesThreshold(minTree_, minArea_, threshold);
        iterStats_.phase2_candidates_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());

        phaseStart = std::chrono::steady_clock::now();
        adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
        iterStats_.phase2_adjust_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());
        iterStats_.phase2_refresh_area_us = 0;
    }

    ImageUInt8Ptr reconstructImage() const {
        auto phaseStart = std::chrono::steady_clock::now();
        auto out = minTree_.reconstructionImage();
        const_cast<IterationStats &>(iterStats_).reconstruct_us =
            elapsedMicros(phaseStart, std::chrono::steady_clock::now());
        return out;
    }

    const InitializationStats &getInitializationStats() const { return initStats_; }
    const IterationStats &getIterationStats() const { return iterStats_; }
};

class DynamicComponentTreeLeafCasfRunner {
public:
    using InitializationStats = DynamicComponentTreeCasfRunner::InitializationStats;
    using IterationStats = DynamicComponentTreeCasfRunner::IterationStats;

private:
    InitializationStats initStats_{};
    ImageUInt8Ptr inputImage_;
    AdjacencyRelationPtr adj_;
    DynamicComponentTree maxTree_;
    DynamicComponentTree minTree_;
    DynamicComponentTreeAdjustmentLeaf<AltitudeType> adjust_;
    std::unique_ptr<DynamicAttributeComputer> maxAttributeComputer_;
    std::unique_ptr<DynamicAttributeComputer> minAttributeComputer_;
    Attribute attribute_;
    IterationStats iterStats_{};
    std::vector<float> maxArea_;
    std::vector<float> minArea_;

public:
    explicit DynamicComponentTreeLeafCasfRunner(ImageUInt8Ptr image, double radioAdj, Attribute attribute)
        : inputImage_(timedConstruct(initStats_.input_clone_us, [&]() { return image->clone(); })),
          adj_(timedConstruct(initStats_.adjacency_us, [&]() {
              return std::make_shared<AdjacencyRelation>(image->getNumRows(), image->getNumCols(), radioAdj);
          })),
          maxTree_(timedConstruct(initStats_.max_tree_us, [&]() { return DynamicComponentTree(inputImage_, true, adj_); })),
          minTree_(timedConstruct(initStats_.min_tree_us, [&]() { return DynamicComponentTree(inputImage_, false, adj_); })),
          adjust_(&minTree_, &maxTree_, *adj_),
          maxAttributeComputer_(timedConstruct(initStats_.max_attribute_computer_us, [&]() {
              return makeAttributeComputer(&maxTree_, attribute);
          })),
          minAttributeComputer_(timedConstruct(initStats_.min_attribute_computer_us, [&]() {
              return makeAttributeComputer(&minTree_, attribute);
          })),
          attribute_(attribute),
          maxArea_(maxTree_.getGlobalIdSpaceSize(), 0.0f),
          minArea_(minTree_.getGlobalIdSpaceSize(), 0.0f) {
        timedCall(initStats_.max_attribute_values_us, [&]() {
            maxArea_ = computeAttributeVector(&maxTree_, attribute_);
        });
        timedCall(initStats_.min_attribute_values_us, [&]() {
            minArea_ = computeAttributeVector(&minTree_, attribute_);
        });
        timedCall(initStats_.attribute_binding_us, [&]() {
            adjust_.setAttributeComputer(*minAttributeComputer_,
                                         *maxAttributeComputer_,
                                         std::span<float>(minArea_),
                                         std::span<float>(maxArea_));
        });
    }

    void applyThreshold(int threshold) {
        iterStats_.reset();
        auto phaseStart = std::chrono::steady_clock::now();
        auto nodesToPrune = getDynamicNodesThreshold(maxTree_, maxArea_, threshold);
        iterStats_.phase1_candidates_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());

        phaseStart = std::chrono::steady_clock::now();
        adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
        iterStats_.phase1_adjust_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());
        iterStats_.phase1_refresh_area_us = 0;

        phaseStart = std::chrono::steady_clock::now();
        nodesToPrune = getDynamicNodesThreshold(minTree_, minArea_, threshold);
        iterStats_.phase2_candidates_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());

        phaseStart = std::chrono::steady_clock::now();
        adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
        iterStats_.phase2_adjust_us = elapsedMicros(phaseStart, std::chrono::steady_clock::now());
        iterStats_.phase2_refresh_area_us = 0;
    }

    ImageUInt8Ptr reconstructImage() const {
        auto phaseStart = std::chrono::steady_clock::now();
        auto out = minTree_.reconstructionImage();
        const_cast<IterationStats &>(iterStats_).reconstruct_us =
            elapsedMicros(phaseStart, std::chrono::steady_clock::now());
        return out;
    }

    const InitializationStats &getInitializationStats() const { return initStats_; }
    const IterationStats &getIterationStats() const { return iterStats_; }
};

inline long long totalInitializationMicros(const DynamicComponentTreeCasfRunner::InitializationStats &stats) {
    return stats.adjacency_us +
           stats.max_tree_us +
           stats.min_tree_us +
           stats.max_attribute_computer_us +
           stats.min_attribute_computer_us +
           stats.max_attribute_values_us +
           stats.min_attribute_values_us +
           stats.attribute_binding_us;
}

inline double totalIterationMillis(const DynamicComponentTreeCasfRunner::IterationStats &stats) {
    return microsToMillis(stats.phase1_candidates_us +
                          stats.phase1_adjust_us +
                          stats.phase1_refresh_area_us +
                          stats.phase2_candidates_us +
                          stats.phase2_adjust_us +
                          stats.phase2_refresh_area_us +
                          stats.reconstruct_us);
}

inline void printInitializationStats(const DynamicComponentTreeCasfRunner::InitializationStats &stats,
                                     const char *label) {
    std::cout << "  " << label
              << ": adjacency=" << microsToMillis(stats.adjacency_us) << " ms"
              << " max_tree=" << microsToMillis(stats.max_tree_us) << " ms"
              << " min_tree=" << microsToMillis(stats.min_tree_us) << " ms"
              << " max_attribute_computer=" << microsToMillis(stats.max_attribute_computer_us) << " ms"
              << " min_attribute_computer=" << microsToMillis(stats.min_attribute_computer_us) << " ms"
              << " max_attribute_values=" << microsToMillis(stats.max_attribute_values_us) << " ms"
              << " min_attribute_values=" << microsToMillis(stats.min_attribute_values_us) << " ms"
              << " attribute_binding=" << microsToMillis(stats.attribute_binding_us) << " ms"
              << "\n";
}

inline void printDynamicIterationStats(const DynamicComponentTreeCasfRunner::IterationStats &stats,
                                       const char *label) {
    std::cout << "  " << label
              << ": phase1_candidates=" << microsToMillis(stats.phase1_candidates_us) << " ms"
              << " phase1_adjust=" << microsToMillis(stats.phase1_adjust_us) << " ms"
              << " phase1_refresh_area=" << microsToMillis(stats.phase1_refresh_area_us) << " ms"
              << " phase2_candidates=" << microsToMillis(stats.phase2_candidates_us) << " ms"
              << " phase2_adjust=" << microsToMillis(stats.phase2_adjust_us) << " ms"
              << " phase2_refresh_area=" << microsToMillis(stats.phase2_refresh_area_us) << " ms"
              << " reconstruct=" << microsToMillis(stats.reconstruct_us) << " ms"
              << "\n";
}

inline void printDynamicLeafIterationStats(const DynamicComponentTreeLeafCasfRunner::IterationStats &stats,
                                           const char *label) {
    std::cout << "  " << label
              << ": phase1_candidates=" << microsToMillis(stats.phase1_candidates_us) << " ms"
              << " phase1_adjust=" << microsToMillis(stats.phase1_adjust_us) << " ms"
              << " phase1_refresh_area=" << microsToMillis(stats.phase1_refresh_area_us) << " ms"
              << " phase2_candidates=" << microsToMillis(stats.phase2_candidates_us) << " ms"
              << " phase2_adjust=" << microsToMillis(stats.phase2_adjust_us) << " ms"
              << " phase2_refresh_area=" << microsToMillis(stats.phase2_refresh_area_us) << " ms"
              << " reconstruct=" << microsToMillis(stats.reconstruct_us) << " ms"
              << "\n";
}

} // namespace dynamic_casf_apply_common
