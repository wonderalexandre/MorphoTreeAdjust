
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>

class AdjacencyRelation;  // forward declaration
using AdjacencyRelationPtr = std::shared_ptr<AdjacencyRelation>;

/**
 * @brief 2D grid adjacency relation with a symmetric stencil and efficient iteration.
 *
 * Defines neighborhood offsets for a symmetric stencil around the origin. The
 * radius-based constructor preserves the existing API for circular adjacencies
 * (for example, 1.0 -> 4-adj, 1.5 -> 8-adj), but the class also accepts
 * arbitrary symmetric offsets and provides a factory for rectangular
 * adjacency. It also offers a "forward" variant that emits only half of the
 * neighbors, which is useful for asymmetric scans and unique-edge
 * construction.
 */
class AdjacencyRelation {
public:
    using Offset = std::pair<int, int>; // (deltaRow, deltaCol)

private:
    using OffsetKey = std::int64_t;

    int id;
    
    int row;
    int col;    
    int numCols;
    int numRows;
    double radius;
    double radius2;
    int n;
    bool forwardOnly = false; // half stencil only

    std::vector<int> offsetRow;
    std::vector<int> offsetCol;
    std::vector<uint8_t> forwardMask; // "forward" mask per offset i: true if (dy>0) || (dy==0 && dx>0)
    std::unordered_set<OffsetKey> offsetLookup;

    static OffsetKey encodeOffset(int dy, int dx) noexcept;
    static double offsetAngle(int dy, int dx) noexcept;
    static std::vector<Offset> buildCircularOffsets(double radius);
    static std::vector<Offset> buildRectangularOffsets(int halfHeight, int halfWidth);
    static std::vector<Offset> normalizeOffsets(const std::vector<Offset> &offsets);
    void initialize(int numRows, int numCols, const std::vector<Offset> &offsets, double radiusValue);


public:
    /**
     * @brief Builds an adjacency relation for a `numRows`x`numCols` image.
     * @param numRows Number of image rows.
     * @param numCols Number of image columns.
     * @param radius Neighborhood radius (1.0 ~= 4-connectivity, 1.5 ~= 8-connectivity).
     */
    AdjacencyRelation(int numRows, int numCols, double radius);
    /**
     * @brief Builds an adjacency relation from arbitrary symmetric offsets.
     * @param offsets List of offsets (deltaRow, deltaCol), which must contain the origin and be
     * centrally symmetric: if (dy,dx) belongs to the set, (-dy,-dx) must belong as well.
     */
    AdjacencyRelation(int numRows, int numCols, const std::vector<Offset> &offsets);
    /**
     * @brief Convenience helper for building a symmetric rectangular adjacency.
     * @param halfHeight Half-height of the rectangle.
     * @param halfWidth Half-width of the rectangle.
     */
    static AdjacencyRelation rectangular(int numRows, int numCols, int halfHeight, int halfWidth);
    /**
     * @brief Advances to the next valid offset according to bounds and mask.
     * @return Index of the next valid offset, or the size for end-of-range.
     */
    int nextValid();
    /**
     * @brief Returns the number of offsets in the current stencil.
     */
    int getSize();
    /**
     * @brief Configures (row,col) and prepares adjacency iteration without a forward filter. This method includes the origin.
     */
    AdjacencyRelation& getAdjPixels(int row, int col);
    /**
     * @brief Configures by linear index and prepares adjacency iteration without a filter. This method includes the origin.
     */
    AdjacencyRelation& getAdjPixels(int index);
    /**
     * @brief Configures (row,col) and prepares iteration over in-bounds neighbors. This method does NOT include the origin.
     */
    AdjacencyRelation& getNeighborPixels(int row, int col);
    /**
     * @brief Configures by linear index and prepares iteration over in-bounds neighbors. This method does NOT include the origin.
     */
    AdjacencyRelation& getNeighborPixels(int index);
    /**
     * @brief Prepares iteration over only half of the neighbors (forward-only) at (row,col).
     * @note Useful for generating pairs (p,q) without duplication (p->q, never q->p).
     */
    AdjacencyRelation& getNeighborPixelsForward(int row, int col);
    /**
     * @brief Linear-index version of `getNeighborPixelsForward`.
     */
    AdjacencyRelation& getNeighborPixelsForward(int index);

    /**
     * @brief Checks adjacency by linear indices (p,q).
     */
    inline bool isAdjacent(int p, int q)const noexcept;
    /**
     * @brief Checks adjacency by coordinates (px,py) and (qx,qy).
     */
    inline bool isAdjacent(int px, int py, int qx, int qy) const noexcept;
    /**
     * @brief Returns the radius in use. For non-circular adjacencies, returns a negative value.
     */
    double getRadius();
    /**
     * @brief Indicates whether the adjacency was built from a circular radius.
     */
    bool hasRadius() const noexcept;

    int getOffsetRow(int index){
        return offsetRow[index];
    }
    int getOffsetCol(int index){
        return offsetCol[index];
    }
    
    /**
     * @brief Lightweight iterator for traversing neighbors already configured via `get*`.
     *
     * Produces linear indices of valid neighboring pixels while respecting
     * bounds and, when configured, the forward-only mask.
     */
    class IteratorAdjacency { 
    private:
        int index;
        AdjacencyRelation* instance; 

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = int;

        IteratorAdjacency(AdjacencyRelation* obj, int id) :  index(id), instance(obj) { }

        AdjacencyRelation* getInstance() { return instance; } 

        IteratorAdjacency& operator++() { 
            this->index = instance->nextValid();  
            return *this; 
        }

        bool operator==(const IteratorAdjacency& other) const { 
            return index == other.index; 
        }
        bool operator!=(const IteratorAdjacency& other) const { 
            return !(*this == other);
        }

        int operator*() const { 
            return (instance->row + instance->offsetRow[index]) * instance->numCols + (instance->col + instance->offsetCol[index]); 
        }
    };
    /**
     * @brief Start of neighbor iteration according to the current configuration.
     */
    IteratorAdjacency begin();
    /**
     * @brief End marker for neighbor iteration.
     */
   IteratorAdjacency end();	 
};

inline AdjacencyRelation::OffsetKey AdjacencyRelation::encodeOffset(int dy, int dx) noexcept {
    return (static_cast<OffsetKey>(dy) << 32) ^ static_cast<std::uint32_t>(dx);
}

inline double AdjacencyRelation::offsetAngle(int dy, int dx) noexcept {
    constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
    double angle = std::atan2(-static_cast<double>(dy), -static_cast<double>(dx));
    if (angle < 0.0) {
        angle += kTwoPi;
    }
    return angle;
}

inline std::vector<AdjacencyRelation::Offset> AdjacencyRelation::buildCircularOffsets(double radius) {
    if (radius < 0.0) {
        throw std::invalid_argument("The adjacency radius must be non-negative.");
    }

    const int r0 = static_cast<int>(radius);
    const int r2 = static_cast<int>(radius * radius);

    std::vector<Offset> offsets;
    offsets.reserve(static_cast<size_t>((2 * r0 + 1) * (2 * r0 + 1)));

    for (int dy = -r0; dy <= r0; ++dy) {
        for (int dx = -r0; dx <= r0; ++dx) {
            if ((dx * dx) + (dy * dy) <= r2) {
                offsets.emplace_back(dy, dx);
            }
        }
    }

    return offsets;
}

inline std::vector<AdjacencyRelation::Offset> AdjacencyRelation::buildRectangularOffsets(int halfHeight, int halfWidth) {
    if (halfHeight < 0 || halfWidth < 0) {
        throw std::invalid_argument("The rectangle half-dimensions must be non-negative.");
    }

    std::vector<Offset> offsets;
    offsets.reserve(static_cast<size_t>((2 * halfHeight + 1) * (2 * halfWidth + 1)));

    for (int dy = -halfHeight; dy <= halfHeight; ++dy) {
        for (int dx = -halfWidth; dx <= halfWidth; ++dx) {
            offsets.emplace_back(dy, dx);
        }
    }

    return offsets;
}

inline std::vector<AdjacencyRelation::Offset> AdjacencyRelation::normalizeOffsets(const std::vector<Offset> &offsets) {
    if (offsets.empty()) {
        throw std::invalid_argument("The adjacency offset list cannot be empty.");
    }

    std::vector<Offset> uniqueOffsets;
    uniqueOffsets.reserve(offsets.size());

    std::unordered_set<OffsetKey> seen;
    seen.reserve(offsets.size());

    bool hasOrigin = false;
    for (const auto &[dy, dx] : offsets) {
        const OffsetKey key = encodeOffset(dy, dx);
        if (!seen.insert(key).second) {
            continue;
        }

        uniqueOffsets.emplace_back(dy, dx);
        hasOrigin = hasOrigin || (dy == 0 && dx == 0);
    }

    if (!hasOrigin) {
        throw std::invalid_argument("The adjacency must contain the origin (0,0).");
    }

    for (const auto &[dy, dx] : uniqueOffsets) {
        if (seen.find(encodeOffset(-dy, -dx)) == seen.end()) {
            throw std::invalid_argument("The offset list must be centrally symmetric.");
        }
    }

    std::sort(uniqueOffsets.begin(), uniqueOffsets.end(), [](const Offset &lhs, const Offset &rhs) {
        const bool lhsIsOrigin = lhs.first == 0 && lhs.second == 0;
        const bool rhsIsOrigin = rhs.first == 0 && rhs.second == 0;

        if (lhsIsOrigin != rhsIsOrigin) {
            return lhsIsOrigin;
        }

        const double lhsAngle = offsetAngle(lhs.first, lhs.second);
        const double rhsAngle = offsetAngle(rhs.first, rhs.second);
        if (lhsAngle != rhsAngle) {
            return lhsAngle < rhsAngle;
        }

        const int lhsRadius2 = lhs.first * lhs.first + lhs.second * lhs.second;
        const int rhsRadius2 = rhs.first * rhs.first + rhs.second * rhs.second;
        if (lhsRadius2 != rhsRadius2) {
            return lhsRadius2 < rhsRadius2;
        }

        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });

    return uniqueOffsets;
}

inline void AdjacencyRelation::initialize(int numRows, int numCols, const std::vector<Offset> &offsets, double radiusValue) {
    if (numRows <= 0 || numCols <= 0) {
        throw std::invalid_argument("The image dimensions must be positive.");
    }

    const std::vector<Offset> normalizedOffsets = normalizeOffsets(offsets);

    this->numRows = numRows;
    this->numCols = numCols;
    this->radius = radiusValue;
    this->radius2 = radiusValue >= 0.0 ? radiusValue * radiusValue : -1.0;
    this->n = static_cast<int>(normalizedOffsets.size());
    this->id = -1;
    this->row = 0;
    this->col = 0;
    this->forwardOnly = false;

    this->offsetRow.resize(this->n);
    this->offsetCol.resize(this->n);
    this->forwardMask.assign(this->n, 0);
    this->offsetLookup.clear();
    this->offsetLookup.reserve(normalizedOffsets.size());

    for (int i = 0; i < this->n; ++i) {
        this->offsetRow[i] = normalizedOffsets[static_cast<size_t>(i)].first;
        this->offsetCol[i] = normalizedOffsets[static_cast<size_t>(i)].second;
        this->offsetLookup.insert(encodeOffset(this->offsetRow[i], this->offsetCol[i]));

        if (i == 0) {
            continue;
        }

        const int offsetDx = this->offsetCol[i];
        const int offsetDy = this->offsetRow[i];
        this->forwardMask[i] = (offsetDy > 0 || (offsetDy == 0 && offsetDx > 0)) ? 1 : 0;
    }
}

inline AdjacencyRelation::AdjacencyRelation(int numRows, int numCols, double radius) {
    initialize(numRows, numCols, buildCircularOffsets(radius), radius);
}

inline AdjacencyRelation::AdjacencyRelation(int numRows, int numCols, const std::vector<Offset> &offsets) {
    initialize(numRows, numCols, offsets, -1.0);
}

inline AdjacencyRelation AdjacencyRelation::rectangular(int numRows, int numCols, int halfHeight, int halfWidth) {
    return AdjacencyRelation(numRows, numCols, buildRectangularOffsets(halfHeight, halfWidth));
}

inline int AdjacencyRelation::getSize() {
    return this->n;
}

inline double AdjacencyRelation::getRadius() {
    return this->radius;
}

inline bool AdjacencyRelation::hasRadius() const noexcept {
    return this->radius >= 0.0;
}

inline bool AdjacencyRelation::isAdjacent(int px, int py, int qx, int qy) const noexcept {
    const int dx = px - qx;
    const int dy = py - qy;
    return offsetLookup.find(encodeOffset(dy, dx)) != offsetLookup.end();
}

inline bool AdjacencyRelation::isAdjacent(int p, int q) const noexcept {
    const int py = p / numCols;
    const int px = p % numCols;
    const int qy = q / numCols;
    const int qx = q % numCols;

    return isAdjacent(px, py, qx, qy);
}

inline int AdjacencyRelation::nextValid() {
    id += 1;
    while (id < n) {
        if (forwardOnly && !forwardMask[id]) {
            id += 1;
            continue;
        }

        const int newRow = row + offsetRow[id];
        const int newCol = col + offsetCol[id];

        if (newRow >= 0 && newRow < numRows && newCol >= 0 && newCol < numCols) {
            return id;
        }
        id += 1;
    }
    return n;
}

inline AdjacencyRelation::IteratorAdjacency AdjacencyRelation::begin() {
    return IteratorAdjacency(this, nextValid());
}

inline AdjacencyRelation::IteratorAdjacency AdjacencyRelation::end() {
    return IteratorAdjacency(this, this->n);
}

inline AdjacencyRelation& AdjacencyRelation::getAdjPixels(int row, int col) {
    if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Indice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = -1;
    this->forwardOnly = false;

    return *this;
}

inline AdjacencyRelation& AdjacencyRelation::getAdjPixels(int indexVector) {
    return getAdjPixels(indexVector / this->numCols, indexVector % this->numCols);
}

inline AdjacencyRelation& AdjacencyRelation::getNeighborPixels(int row, int col) {
    if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Indice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = 0;
    this->forwardOnly = false;
    return *this;
}

inline AdjacencyRelation& AdjacencyRelation::getNeighborPixels(int indexVector) {
    return getNeighborPixels(indexVector / this->numCols, indexVector % this->numCols);
}

inline AdjacencyRelation& AdjacencyRelation::getNeighborPixelsForward(int row, int col) {
    if (row < 0 || row >= this->numRows || col < 0 || col >= this->numCols) {
        throw std::out_of_range("Indice fora dos limites.");
    }
    this->row = row;
    this->col = col;
    this->id = 0;
    this->forwardOnly = true;
    return *this;
}

inline AdjacencyRelation& AdjacencyRelation::getNeighborPixelsForward(int indexVector) {
    return getNeighborPixelsForward(indexVector / this->numCols, indexVector % this->numCols);
}
