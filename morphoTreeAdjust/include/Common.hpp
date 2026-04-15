#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

using NodeId = int;
using PixelId = int;
constexpr NodeId InvalidNode = -1;

/**
 * @brief Utility conversions between linear indices and 2D coordinates.
 */
class ImageUtils {
public:
    /**
     * @brief Converts a linear index `(row * numCols + col)` to `(row, col)`.
     */
    static std::pair<int, int> to2D(int index, int numCols) {
        const int row = index / numCols;
        const int col = index % numCols;
        return std::make_pair(row, col);
    }
};

template <typename PixelType>
class Image {
private:
    int numRows;
    int numCols;
    std::shared_ptr<PixelType[]> data;
    using Ptr = std::shared_ptr<Image<PixelType>>;

public:
    /**
     * @brief Creates an owning image with contiguous storage.
     */
    Image(int rows, int cols)
        : numRows(rows),
          numCols(cols),
          data(new PixelType[rows * cols], std::default_delete<PixelType[]>()) {}

    /**
     * @brief Factory for an uninitialized image.
     */
    static Ptr create(int rows, int cols) {
        return std::make_shared<Image>(rows, cols);
    }

    /**
     * @brief Factory for an image filled with `initValue`.
     */
    static Ptr create(int rows, int cols, PixelType initValue) {
        auto img = create(rows, cols);
        img->fill(initValue);
        return img;
    }

    /**
     * @brief Fills every image pixel with a constant value.
     */
    void fill(PixelType value) {
        std::fill_n(data.get(), numRows * numCols, value);
    }

    /**
     * @brief Compares two images by size and content.
     */
    bool isEqual(const Ptr& other) const {
        if (numRows != other->numRows || numCols != other->numCols) {
            return false;
        }
        const int n = numRows * numCols;
        for (int i = 0; i < n; ++i) {
            if (data[i] != (*other)[i]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Produces a deep copy of the image.
     */
    Ptr clone() const {
        auto newImg = create(numRows, numCols);
        std::copy(data.get(), data.get() + (numRows * numCols), newImg->data.get());
        return newImg;
    }

    /**
     * @brief Shared pointer to the raw contiguous buffer.
     */
    std::shared_ptr<PixelType[]> rawDataPtr() { return data; }
    /**
     * @brief Raw pointer to the first pixel.
     */
    PixelType* rawData() { return data.get(); }
    /**
     * @brief Number of image rows.
     */
    int getNumRows() const { return numRows; }
    /**
     * @brief Number of image columns.
     */
    int getNumCols() const { return numCols; }
    /**
     * @brief Total number of pixels.
     */
    int getSize() const { return numRows * numCols; }
    PixelType& operator[](int index) { return data[index]; }
    const PixelType& operator[](int index) const { return data[index]; }
};

using ImageUInt8 = Image<uint8_t>;
using AltitudeType = uint8_t;
using ImageUInt8Ptr = std::shared_ptr<ImageUInt8>;

template <typename T>
using ImagePtr = std::shared_ptr<Image<T>>;

/**
 * @brief Generation-based mark set for deduplication without O(n) clearing on every use.
 */
struct GenerationStampSet {
    using gen_t = uint32_t;

    std::vector<gen_t> stamp;
    size_t n{0};
    gen_t cur{1};

    GenerationStampSet() = default;

    /**
     * @brief Builds the set with `n` positions.
     */
    explicit GenerationStampSet(size_t n) { resize(n); }

    /**
     * @brief Resizes the set and resets the generation counter.
     */
    void resize(size_t newN) {
        n = newN;
        stamp.assign(n, 0);
        cur = 1;
    }

    /**
     * @brief Marks position `idx` in the current generation.
     */
    inline void mark(size_t idx) noexcept {
        stamp[idx] = cur;
    }

    /**
     * @brief Removes the mark from position `idx` in the current generation.
     */
    inline void unmark(size_t idx) noexcept {
        if (stamp[idx] == cur) {
            stamp[idx] = 0;
        }
    }

    /**
     * @brief Tests whether `idx` is marked in the current generation.
     */
    inline bool isMarked(size_t idx) const noexcept {
        return stamp[idx] == cur;
    }

    /**
     * @brief Advances to a new logical generation.
     *
     * Performs a full physical clear only when the generation counter overflows.
     */
    void resetAll() {
        if (++cur == 0) {
            std::fill(stamp.begin(), stamp.end(), 0);
            cur = 1;
        }
    }
};

/**
 * @brief Lightweight FIFO queue based on `std::vector` and a head index.
 */
template <typename T>
struct FastQueue {
private:
    std::vector<T> data_;
    size_t head_ = 0;

public:
    FastQueue() = default;

    /**
     * @brief Builds the queue with optional initial reserve.
     */
    explicit FastQueue(size_t n) { data_.reserve(n); }

    /**
     * @brief Reserves additional capacity.
     */
    void reserve(size_t n) { data_.reserve(n); }
    /**
     * @brief Removes all elements and resets the logical head.
     */
    void clear() {
        data_.clear();
        head_ = 0;
    }
    /**
     * @brief Tests whether no elements remain.
     */
    bool empty() const { return head_ >= data_.size(); }
    /**
     * @brief Number of elements not yet consumed.
     */
    size_t size() const { return data_.size() - head_; }
    /**
     * @brief Enqueues a copy of `value`.
     */
    void push(const T& value) { data_.push_back(value); }
    /**
     * @brief Enqueues `value` by move.
     */
    void push(T&& value) { data_.push_back(std::move(value)); }
    /**
     * @brief Removes and returns the first logical element in the queue.
     */
    T pop() { return std::move(data_[head_++]); }
    /**
     * @brief Accesses the first logical element without removing it.
     */
    T& front() { return data_[head_]; }
    /**
     * @brief Const access to the first logical element without removing it.
     */
    const T& front() const { return data_[head_]; }
};

/**
 * @brief Simple high-performance stack based on `std::vector`.
 */
template <typename T>
struct FastStack {
private:
    std::vector<T> data_;

public:
    FastStack() = default;

    explicit FastStack(size_t n) {
        data_.reserve(n);
    }

    void reserve(size_t n) { data_.reserve(n); }
    void clear() { data_.clear(); }
    bool empty() const { return data_.empty(); }
    size_t size() const { return data_.size(); }
    void push(const T& value) { data_.push_back(value); }
    void push(T&& value) { data_.push_back(std::move(value)); }

    T pop() {
        T value = std::move(data_.back());
        data_.pop_back();
        return value;
    }

    T& top() { return data_.back(); }
    const T& top() const { return data_.back(); }
};

class Stopwatch {
public:
    using clock = std::chrono::steady_clock;

    /**
     * @brief Starts the stopwatch, clearing any previously accumulated time.
     */
    void start() {
        accumulated_ = clock::duration::zero();
        running_ = true;
        last_ = clock::now();
    }

    /**
     * @brief Pauses timing and accumulates the elapsed interval.
     */
    void pause() {
        if (!running_) return;
        accumulated_ += clock::now() - last_;
        running_ = false;
    }

    /**
     * @brief Resumes timing from the already accumulated value.
     */
    void resume() {
        if (running_) return;
        running_ = true;
        last_ = clock::now();
    }

    clock::duration elapsed() const {
        if (!running_) return accumulated_;
        return accumulated_ + (clock::now() - last_);
    }

    /**
     * @brief Indicates whether the stopwatch is currently running.
     */
    bool running() const { return running_; }

private:
    clock::time_point last_{};
    clock::duration accumulated_{clock::duration::zero()};
    bool running_{false};
};
