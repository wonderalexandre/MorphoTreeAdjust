#include "include/AdjacencyRelation.hpp"
#include "include/AttributeComputer.hpp"
#include "include/Common.hpp"
#include "include/ComponentTreeCasf.hpp"
#include "include/ComponentTreeCasfNonIncreasing.hpp"
#include "include/DynamicComponentTree.hpp"
#include "include/DualMinMaxTreeIncrementalFilter.hpp"

#include <cctype>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace {

ImageUInt8Ptr image_from_numpy(const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input) {
    const py::buffer_info info = input.request();
    if (info.ndim != 2) {
        throw std::runtime_error("Expected a 2D uint8 NumPy array.");
    }

    const int numRows = static_cast<int>(info.shape[0]);
    const int numCols = static_cast<int>(info.shape[1]);
    auto image = ImageUInt8::create(numRows, numCols);
    std::copy_n(static_cast<const uint8_t *>(info.ptr), static_cast<size_t>(numRows * numCols), image->rawData());
    return image;
}

py::array_t<uint8_t> numpy_from_image(const ImageUInt8Ptr &image) {
    if (image == nullptr) {
        return py::array_t<uint8_t>();
    }

    py::capsule owner(new ImageUInt8Ptr(image), [](void *ptr) {
        delete static_cast<ImageUInt8Ptr *>(ptr);
    });

    return py::array_t<uint8_t>(
        {image->getNumRows(), image->getNumCols()},
        {static_cast<py::ssize_t>(sizeof(uint8_t) * image->getNumCols()), static_cast<py::ssize_t>(sizeof(uint8_t))},
        image->rawData(),
        owner
    );
}

Attribute parse_attribute_string(const std::string &attribute) {
    std::string normalized = attribute;
    for (char &c : normalized) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (normalized == "area") {
        return AREA;
    }
    if (normalized == "bbox_width" || normalized == "bbox-width" ||
        normalized == "box_width" || normalized == "box-width" || normalized == "width") {
        return BOX_WIDTH;
    }
    if (normalized == "bbox_height" || normalized == "bbox-height" ||
        normalized == "box_height" || normalized == "box-height" || normalized == "height") {
        return BOX_HEIGHT;
    }
    if (normalized == "bbox_diagonal" || normalized == "bbox-diagonal" ||
        normalized == "diagonal_length" || normalized == "diagonal-length" || normalized == "diagonal") {
        return DIAGONAL_LENGTH;
    }
    if (normalized == "area_per_box_perimeter" || normalized == "area-per-box-perimeter" ||
        normalized == "box_area_perimeter_ratio" || normalized == "box-area-perimeter-ratio") {
        return AREA_PER_BOX_PERIMETER;
    }
    if (normalized == "rectangularity" ||
        normalized == "bbox_rectangularity" || normalized == "bbox-rectangularity" ||
        normalized == "box_fill_ratio" || normalized == "box-fill-ratio") {
        return RECTANGULARITY;
    }
    if (normalized == "box_aspect_ratio" || normalized == "box-aspect-ratio" ||
        normalized == "bbox_aspect_ratio" || normalized == "bbox-aspect-ratio" ||
        normalized == "aspect_ratio" || normalized == "aspect-ratio") {
        return BOX_ASPECT_RATIO;
    }
    if (normalized == "box_elongation" || normalized == "box-elongation" ||
        normalized == "bbox_elongation" || normalized == "bbox-elongation" ||
        normalized == "elongation") {
        return BOX_ELONGATION;
    }
    if (normalized == "area_per_box_diagonal_squared" || normalized == "area-per-box-diagonal-squared" ||
        normalized == "area_per_bbox_diagonal_squared" || normalized == "area-per-bbox-diagonal-squared" ||
        normalized == "box_diagonal_squared_ratio" || normalized == "box-diagonal-squared-ratio") {
        return AREA_PER_BOX_DIAGONAL_SQUARED;
    }

    throw std::runtime_error(
        "Unknown attribute. Expected one of: area, bbox_width, bbox_height, bbox_diagonal, "
        "area_per_box_perimeter, rectangularity, box_aspect_ratio, box_elongation, "
        "area_per_box_diagonal_squared.");
}

std::vector<NodeId> alive_nodes(const DynamicComponentTree &tree) {
    std::vector<NodeId> nodes;
    nodes.reserve(static_cast<size_t>(tree.getNumNodes()));
    for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
        if (tree.isAlive(nodeId)) {
            nodes.push_back(nodeId);
        }
    }
    return nodes;
}

std::vector<NodeId> children_of(const DynamicComponentTree &tree, NodeId nodeId) {
    std::vector<NodeId> children;
    children.reserve(static_cast<size_t>(std::max(0, tree.getNumChildren(nodeId))));
    for (NodeId childId : tree.getChildren(nodeId)) {
        children.push_back(childId);
    }
    return children;
}

std::vector<NodeId> proper_parts_of(const DynamicComponentTree &tree, NodeId nodeId) {
    std::vector<NodeId> properParts;
    properParts.reserve(static_cast<size_t>(std::max(0, tree.getNumProperParts(nodeId))));
    for (NodeId pixelId : tree.getProperParts(nodeId)) {
        properParts.push_back(pixelId);
    }
    return properParts;
}

std::vector<NodeId> subtree_nodes_of(const DynamicComponentTree &tree, NodeId rootNodeId) {
    std::vector<NodeId> subtree;
    for (NodeId nodeId : tree.getNodeSubtree(rootNodeId)) {
        subtree.push_back(nodeId);
    }
    return subtree;
}

std::vector<NodeId> breadth_first_nodes_of(const DynamicComponentTree &tree) {
    std::vector<NodeId> nodes;
    nodes.reserve(static_cast<size_t>(tree.getNumNodes()));
    for (NodeId nodeId : tree.getIteratorBreadthFirstTraversal()) {
        nodes.push_back(nodeId);
    }
    return nodes;
}

class PyDualMinMaxTreeIncrementalFilter {
private:
    std::shared_ptr<DynamicComponentTree> mintree_;
    std::shared_ptr<DynamicComponentTree> maxtree_;
    DynamicAreaComputer minAreaComputer_;
    DynamicAreaComputer maxAreaComputer_;
    std::vector<float> minArea_;
    std::vector<float> maxArea_;
    DualMinMaxTreeIncrementalFilter<AltitudeType> adjust_;

    static DynamicComponentTree *requireTree(const std::shared_ptr<DynamicComponentTree> &tree, const char *name) {
        if (tree == nullptr) {
            throw std::runtime_error(std::string("DualMinMaxTreeIncrementalFilter requires a valid ") + name + ".");
        }
        return tree.get();
    }

    static AdjacencyRelation &requireAdjacency(const std::shared_ptr<DynamicComponentTree> &mintree,
                                               const std::shared_ptr<DynamicComponentTree> &maxtree) {
        DynamicComponentTree *minTreePtr = requireTree(mintree, "min-tree");
        DynamicComponentTree *maxTreePtr = requireTree(maxtree, "max-tree");
        if (minTreePtr->getAdjacencyRelation() == nullptr || maxTreePtr->getAdjacencyRelation() == nullptr) {
            throw std::runtime_error("DualMinMaxTreeIncrementalFilter requires trees with adjacency information.");
        }
        return *minTreePtr->getAdjacencyRelation();
    }

    static std::vector<float> liveNodeAreaSnapshot(const DynamicComponentTree &tree, const std::vector<float> &area) {
        std::vector<float> snapshot = area;
        snapshot.resize(static_cast<size_t>(tree.getNumInternalNodeSlots()), 0.0f);
        for (NodeId nodeId = 0; nodeId < tree.getNumInternalNodeSlots(); ++nodeId) {
            if (!tree.isAlive(nodeId)) {
                snapshot[static_cast<size_t>(nodeId)] = 0.0f;
            }
        }
        return snapshot;
    }

    void refreshAreaBuffers() {
        minArea_.assign(static_cast<size_t>(mintree_->getNumInternalNodeSlots()), 0.0f);
        maxArea_.assign(static_cast<size_t>(maxtree_->getNumInternalNodeSlots()), 0.0f);
        minAreaComputer_.compute(std::span<float>(minArea_));
        maxAreaComputer_.compute(std::span<float>(maxArea_));
        adjust_.setAttributeComputer(minAreaComputer_, maxAreaComputer_, std::span<float>(minArea_), std::span<float>(maxArea_));
    }

public:
    PyDualMinMaxTreeIncrementalFilter(std::shared_ptr<DynamicComponentTree> mintree, std::shared_ptr<DynamicComponentTree> maxtree)
        : mintree_(std::move(mintree)),
          maxtree_(std::move(maxtree)),
          minAreaComputer_(requireTree(mintree_, "min-tree")),
          maxAreaComputer_(requireTree(maxtree_, "max-tree")),
          adjust_(mintree_.get(), maxtree_.get(), requireAdjacency(mintree_, maxtree_)) {
        refreshAreaBuffers();
    }

    void refreshAttributes() {
        refreshAreaBuffers();
    }

    void updateTree(const std::shared_ptr<DynamicComponentTree> &tree, NodeId subtreeRoot) {
        if (tree.get() != mintree_.get() && tree.get() != maxtree_.get()) {
            throw std::runtime_error("DualMinMaxTreeIncrementalFilter.updateTree requires one of the filter trees.");
        }
        adjust_.updateTree(tree.get(), subtreeRoot);
    }

    void pruneMaxTreeAndUpdateMinTree(std::vector<NodeId> nodesToPrune) {
        adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
    }

    void pruneMinTreeAndUpdateMaxTree(std::vector<NodeId> nodesToPrune) {
        adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
    }

    std::vector<float> getMinArea() const { return liveNodeAreaSnapshot(*mintree_, minArea_); }
    std::vector<float> getMaxArea() const { return liveNodeAreaSnapshot(*maxtree_, maxArea_); }
    std::shared_ptr<DynamicComponentTree> getMinTree() const { return mintree_; }
    std::shared_ptr<DynamicComponentTree> getMaxTree() const { return maxtree_; }
};

class PyComponentTreeCasf : public std::enable_shared_from_this<PyComponentTreeCasf> {
private:
    std::unique_ptr<ComponentTreeCasf<AltitudeType>> casf_;

public:
    PyComponentTreeCasf(const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                        const std::string &attribute,
                        double radiusAdj)
        : casf_(std::make_unique<ComponentTreeCasf<AltitudeType>>(image_from_numpy(input), radiusAdj, parse_attribute_string(attribute))) {}

    PyComponentTreeCasf(const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                        const std::string &attribute,
                        const std::shared_ptr<AdjacencyRelation> &adj)
        : casf_(std::make_unique<ComponentTreeCasf<AltitudeType>>(image_from_numpy(input), adj, parse_attribute_string(attribute))) {}

    py::array_t<uint8_t> filter(const std::vector<int> &thresholds, const std::string &mode = "updating") {
        return numpy_from_image(casf_->filter(thresholds, mode));
    }

    std::shared_ptr<DynamicComponentTree> getMinTree() const {
        return std::shared_ptr<DynamicComponentTree>(this->shared_from_this(), const_cast<DynamicComponentTree *>(&casf_->getMinTree()));
    }

    std::shared_ptr<DynamicComponentTree> getMaxTree() const {
        return std::shared_ptr<DynamicComponentTree>(this->shared_from_this(), const_cast<DynamicComponentTree *>(&casf_->getMaxTree()));
    }
};

class PyComponentTreeCasfNonIncreasing : public std::enable_shared_from_this<PyComponentTreeCasfNonIncreasing> {
private:
    std::unique_ptr<ComponentTreeCasfNonIncreasing<AltitudeType>> casf_;

public:
    PyComponentTreeCasfNonIncreasing(const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                                     const std::string &attribute,
                                     double radiusAdj)
        : casf_(std::make_unique<ComponentTreeCasfNonIncreasing<AltitudeType>>(image_from_numpy(input), radiusAdj, parse_attribute_string(attribute))) {}

    PyComponentTreeCasfNonIncreasing(const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                                     const std::string &attribute,
                                     const std::shared_ptr<AdjacencyRelation> &adj)
        : casf_(std::make_unique<ComponentTreeCasfNonIncreasing<AltitudeType>>(image_from_numpy(input), adj, parse_attribute_string(attribute))) {}

    py::array_t<uint8_t> filter(const std::vector<float> &thresholds) {
        return numpy_from_image(casf_->filter(thresholds));
    }

    std::shared_ptr<DynamicComponentTree> getMinTree() const {
        return std::shared_ptr<DynamicComponentTree>(this->shared_from_this(), const_cast<DynamicComponentTree *>(&casf_->getMinTree()));
    }

    std::shared_ptr<DynamicComponentTree> getMaxTree() const {
        return std::shared_ptr<DynamicComponentTree>(this->shared_from_this(), const_cast<DynamicComponentTree *>(&casf_->getMaxTree()));
    }
};

void init_adjacency_relation(py::module_ &m) {
    py::class_<AdjacencyRelation, std::shared_ptr<AdjacencyRelation>>(
        m,
        "AdjacencyRelation",
        R"doc(Describe a symmetric stencil over a finite 2D image grid.

The stencil always contains the origin.  Configure a pixel with
``getAdjPixels`` and iterate over the returned object to obtain in-bounds,
row-major linear pixel indices.  Configuration is mutable state stored by this
object: a later ``getAdjPixels`` call replaces the active iteration, so nested,
concurrent, or cross-thread iteration over the same instance is unsupported.)doc")
        .def(py::init<int, int, double>(),
             py::arg("numRows"),
             py::arg("numCols"),
             py::arg("radius"),
             R"doc(Create a circular adjacency for an image grid.

``numRows`` and ``numCols`` must be positive and ``radius`` must be
non-negative.  Typical values are 1.0 for 4-connectivity and 1.5 for
8-connectivity.  Raises ``ValueError`` when a dimension or radius is invalid.)doc")
        .def(py::init<int, int, const std::vector<AdjacencyRelation::Offset> &>(),
             py::arg("numRows"),
             py::arg("numCols"),
             py::arg("offsets"),
             R"doc(Create an adjacency from ``(delta_row, delta_col)`` offsets.

Dimensions must be positive.  The non-empty offset collection must contain
``(0, 0)`` and be centrally symmetric: every ``(dy, dx)`` requires
``(-dy, -dx)``.  Duplicate offsets are removed and the stored stencil is put in
a deterministic order.  Invalid input raises ``ValueError``.)doc")
        .def_static("rectangular", [](int numRows, int numCols, int halfHeight, int halfWidth) {
            return std::make_shared<AdjacencyRelation>(
                AdjacencyRelation::rectangular(numRows, numCols, halfHeight, halfWidth)
            );
        },
        py::arg("numRows"),
        py::arg("numCols"),
        py::arg("halfHeight"),
        py::arg("halfWidth"),
        R"doc(Create the full rectangular stencil around the origin.

The result contains every offset in
``[-halfHeight, halfHeight] x [-halfWidth, halfWidth]``.  Image dimensions must
be positive and both half-dimensions non-negative; invalid values raise
``ValueError``.)doc")
        .def_property_readonly(
            "size",
            &AdjacencyRelation::getSize,
            "Number of stored stencil offsets, including the origin.")
        .def("getAdjPixels",
             py::overload_cast<int, int>(&AdjacencyRelation::getAdjPixels),
             py::arg("row"),
             py::arg("col"),
             py::return_value_policy::reference,
             R"doc(Configure iteration around ``(row, col)`` and return ``self``.

Iteration yields every in-bounds stencil position as a row-major linear pixel
index and includes the configured pixel itself because the stencil contains the
origin.  This method mutates the object's active iteration state.)doc")
        .def("getAdjPixels",
             py::overload_cast<int>(&AdjacencyRelation::getAdjPixels),
             py::arg("index"),
             py::return_value_policy::reference,
             R"doc(Configure iteration around a linear pixel index and return ``self``.

The result includes ``index`` itself and all in-bounds stencil positions.  The
index is interpreted in row-major order and must belong to this image grid.
This method mutates the object's active iteration state.)doc")
        .def("getOffsets", [](AdjacencyRelation &self) {
            py::list out;
            for (int i = 0; i < self.getSize(); ++i) {
                out.append(py::make_tuple(self.getOffsetRow(i), self.getOffsetCol(i)));
            }
            return out;
        }, R"doc(Return a new list of stored ``(delta_row, delta_col)`` offsets.

The origin is included.  Changing the returned Python list does not modify the
adjacency.)doc")
        .def("__iter__", [](AdjacencyRelation *adj) {
            return py::make_iterator(adj->begin(), adj->end());
        },
        py::keep_alive<0, 1>(),
        R"doc(Iterate over the pixel configuration selected by ``getAdjPixels``.

The iterator keeps this adjacency alive, but another configuration call on the
same object invalidates its logical iteration state.)doc");
}

void init_dynamic_component_tree(py::module_ &m) {
    py::class_<DynamicComponentTree, std::shared_ptr<DynamicComponentTree>> tree(
        m,
        "DynamicComponentTree",
        R"doc(Mutable pixel-based min-tree or max-tree.

Nodes are represented by integer IDs and pixels by row-major linear indices.
The tree stores direct pixels as a node's *proper parts*; a connected component
also includes the proper parts of every descendant.  Structural operations can
release a node ID and later reuse that same integer for a different node.
Therefore IDs and lists returned by this API are snapshots: check ``isAlive``
and reacquire IDs after mutations.

Most low-level mutation methods intentionally expose the C++ data structure
without defensive validation.  Their documented preconditions are mandatory;
violating them can corrupt the hierarchy rather than raise a Python exception.)doc");

    tree.def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         bool isMaxtree,
                         double radiusAdj) {
            auto img = image_from_numpy(input);
            auto adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radiusAdj);
            return std::make_shared<DynamicComponentTree>(img, isMaxtree, adj);
        }),
        py::arg("image"),
        py::arg("isMaxtree"),
        py::arg("radiusAdj") = 1.5,
        R"doc(Build a component tree using a circular adjacency.

``image`` must be two-dimensional.  It is coerced to C-contiguous ``uint8``
and copied into C++ storage, so later changes to the input do not affect the
tree.  ``isMaxtree=True`` builds a max-tree; ``False`` builds a min-tree.
``radiusAdj`` must be non-negative.)doc")
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         bool isMaxtree,
                         const std::shared_ptr<AdjacencyRelation> &adj) {
            auto image = image_from_numpy(input);
            return std::make_shared<DynamicComponentTree>(image, isMaxtree, adj);
        }),
        py::arg("image"),
        py::arg("isMaxtree"),
        py::arg("adj"),
        R"doc(Build a component tree using an existing adjacency.

``image`` is coerced to C-contiguous ``uint8`` and copied.  ``adj`` is retained
by shared ownership and must have the same grid dimensions as the image.
``isMaxtree`` selects a max-tree when true and a min-tree otherwise.)doc")
        .def("reconstructionImage", [](const DynamicComponentTree &self) {
            return numpy_from_image(self.reconstructionImage());
        }, R"doc(Reconstruct the current hierarchy as a 2D ``uint8`` NumPy array.

Each call returns a newly allocated, writable array whose C++ image storage is
owned by the array.  Editing it does not mutate the tree, and later tree
mutations do not change an array already returned.)doc")
        .def("getNodesThreshold", [](DynamicComponentTree &self, int threshold) {
            return DynamicComponentTree::getNodesThreshold(&self, threshold);
        },
        py::arg("threshold"),
        R"doc(Return the area-threshold frontier as a list of node IDs.

Area is the number of pixels in a node's full connected component.  Traversal
descends through nodes whose area is greater than ``threshold`` and selects the
first nodes at or below it, so selected subtrees do not overlap.  The root may
be returned.  The list is a snapshot and its IDs can become stale or be reused
after mutation.  This query recomputes subtree areas for visited nodes and can
take quadratic time on a deep hierarchy.)doc")
        .def("computeArea", [](DynamicComponentTree &self) {
            DynamicAreaComputer areaComputer(&self);
            return areaComputer.compute();
        }, R"doc(Return a slot-indexed snapshot of connected-component areas.

The returned Python list is indexed by node ID and includes storage slots that
are not currently live; unvisited slots are zero.  It is not updated when the
tree changes.)doc")
        .def("getPixelsOfCC",
             &DynamicComponentTree::getPixelsOfCC,
             py::arg("nodeId"),
             R"doc(Return all row-major pixel indices in a node's component.

Unlike ``getProperParts``, this includes direct proper parts from ``nodeId`` and
from every descendant.  ``nodeId`` must be a currently live node.  A new list
is returned.)doc")
        .def("getChildren", [](const DynamicComponentTree &self, NodeId nodeId) {
            return children_of(self, nodeId);
        },
        py::arg("nodeId"),
        R"doc(Return a new list of the node's direct child IDs.

``nodeId`` must be live.  The list is a snapshot and is not updated after a
topology mutation.)doc")
        .def("getProperParts", [](const DynamicComponentTree &self, NodeId nodeId) {
            return proper_parts_of(self, nodeId);
        },
        py::arg("nodeId"),
        R"doc(Return the node's direct proper-part pixel indices.

Descendant pixels are deliberately excluded; use ``getPixelsOfCC`` for the
full connected component.  ``nodeId`` must be live.  A new list is returned.)doc")
        .def("getNodeSubtree", [](const DynamicComponentTree &self, NodeId nodeId) {
            return subtree_nodes_of(self, nodeId);
        },
        py::arg("nodeId"),
        R"doc(Return a preorder list of IDs in the subtree rooted at ``nodeId``.

The root ID itself is included.  ``nodeId`` must be live.  The returned list is
a snapshot whose IDs may become stale or be reused after mutation.)doc")
        .def("breadthFirstTraversal", [](const DynamicComponentTree &self) {
            return breadth_first_nodes_of(self);
        }, R"doc(Return live node IDs reachable from ``root`` in breadth-first order.

Detached live nodes are not included.  The returned list is a snapshot.)doc")
        .def("moveProperPart",
             &DynamicComponentTree::moveProperPart,
             py::arg("targetNodeId"),
             py::arg("sourceNodeId"),
             py::arg("pixelId"),
             R"doc(Move one direct proper-part pixel from source to target.

Both node IDs must be distinct and live, and ``pixelId`` must be in
``getProperParts(sourceNodeId)``.  The pixel must be in
``[0, numRows * numCols)``.  Preconditions are not checked defensively.)doc")
        .def("moveProperParts",
             &DynamicComponentTree::moveProperParts,
             py::arg("targetNodeId"),
             py::arg("sourceNodeId"),
             R"doc(Move all direct proper parts from source to target.

The distinct node IDs must both be live.  The source remains a live node but
has no direct proper parts afterwards.  Preconditions are not checked
defensively.)doc")
        .def("moveChildren",
             &DynamicComponentTree::moveChildren,
             py::arg("parentId"),
             py::arg("sourceId"),
             R"doc(Reparent every direct child of ``sourceId`` to ``parentId``.

Both distinct IDs must be live, and the move must not create a cycle (in
particular, ``parentId`` cannot be inside a moved child's subtree).  Child order
is preserved and appended after the target's existing children.  Preconditions
are not checked defensively.)doc")
        .def("attachNode",
             &DynamicComponentTree::attachNode,
             py::arg("parentId"),
             py::arg("nodeId"),
             R"doc(Attach a detached live node as the parent's last child.

The IDs must be distinct and live, ``nodeId`` must currently be detached
(``getNodeParent(nodeId) == nodeId``), and ``parentId`` must not lie in its
subtree.  This method does not detach an already attached node and performs no
defensive validation.)doc")
        .def("detachNode",
             &DynamicComponentTree::detachNode,
             py::arg("nodeId"),
             R"doc(Detach an attached non-root node without releasing its ID.

``nodeId`` must be live and have a distinct live parent.  It becomes the root
of an isolated live subtree and is no longer reachable from ``root``.  The
precondition is not checked defensively.)doc")
        .def("removeChild",
             &DynamicComponentTree::removeChild,
             py::arg("parentId"),
             py::arg("childId"),
             py::arg("releaseNodeFlag"),
             R"doc(Remove a direct child, optionally releasing its node ID.

If ``childId`` is not a live direct child of ``parentId``, this method returns
without changing the tree.  With ``releaseNodeFlag=False``, the child becomes a
detached subtree.  With ``True``, the child must already have no children and no
proper parts; its slot is released and may be reused by a later operation.
Releasing a non-empty node discards its links and violates tree invariants.)doc")
        .def_property_readonly("nodes", [](const DynamicComponentTree &self) {
            return alive_nodes(self);
        }, "Snapshot list of all currently live node IDs, including detached nodes.")
        .def_property_readonly(
            "numNodes",
            &DynamicComponentTree::getNumNodes,
            "Number of currently live nodes, including detached live nodes.")
        .def_property_readonly(
            "root",
            &DynamicComponentTree::getRoot,
            "Current hierarchy root ID.")
        .def_property_readonly(
            "isMaxtree",
            &DynamicComponentTree::isMaxtree,
            "Whether this instance is a max-tree (false means min-tree).")
        .def_property_readonly(
            "numRows",
            &DynamicComponentTree::getNumRowsOfImage,
            "Number of rows in the source image grid.")
        .def_property_readonly(
            "numCols",
            &DynamicComponentTree::getNumColsOfImage,
            "Number of columns in the source image grid.")
        .def("getAltitude",
             &DynamicComponentTree::getAltitude,
             py::arg("nodeId"),
             R"doc(Return the level of a currently live node.

The ID is not validated defensively.)doc")
        .def("getNodeParent",
             &DynamicComponentTree::getNodeParent,
             py::arg("nodeId"),
             R"doc(Return the direct parent ID of a live node.

The hierarchy root and a detached subtree root are their own parent.  The ID is
not validated defensively.)doc")
        .def("getNumChildren",
             &DynamicComponentTree::getNumChildren,
             py::arg("nodeId"),
             "Return the number of direct children of a live node.")
        .def("getNumProperParts",
             &DynamicComponentTree::getNumProperParts,
             py::arg("nodeId"),
             R"doc(Return the number of pixels owned directly by a live node.

Descendant pixels are not counted.)doc")
        .def("getSmallestComponent",
             &DynamicComponentTree::getSmallestComponent,
             py::arg("pixelId"),
             R"doc(Return the live node that directly owns a pixel.

``pixelId`` is a row-major index and must be in
``[0, numRows * numCols)``.)doc")
        .def("isAlive",
             &DynamicComponentTree::isAlive,
             py::arg("nodeId"),
             R"doc(Return whether an integer identifies a currently live node.

This check is safe for IDs outside the allocated node-slot domain.  A true
result does not guarantee that the node is reachable from ``root``.)doc")
        .def("isNode",
             &DynamicComponentTree::isNode,
             py::arg("nodeId"),
             R"doc(Return whether an integer is inside the allocated node-slot domain.

A slot can belong to the domain while being released; use ``isAlive`` to test
its current state.)doc")
        .def("isLeaf",
             &DynamicComponentTree::isLeaf,
             py::arg("nodeId"),
             R"doc(Return whether a live node has no direct children.

``nodeId`` must be in the allocated domain and live.)doc")
        .def("hasChild",
             &DynamicComponentTree::hasChild,
             py::arg("parentId"),
             py::arg("childId"),
             R"doc(Return whether ``childId`` has ``parentId`` as its direct parent.

Both IDs must be valid live nodes; they are not checked defensively.)doc");

}

void init_dual_min_max_tree_incremental_filter(py::module_ &m) {
    py::class_<PyDualMinMaxTreeIncrementalFilter>(
        m,
        "DualMinMaxTreeIncrementalFilter",
        R"doc(Keep a mutable min-tree/max-tree pair synchronized during pruning.

The filter retains shared ownership of both Python tree objects and maintains
slot-indexed area buffers for incremental updates.  The trees remain directly
mutable through Python; external structural changes are the caller's
responsibility and may invalidate the paired-tree relationship.  Use
``refreshAttributes`` after a compatible external mutation to recompute area,
but note that recomputing attributes cannot repair structurally inconsistent
min/max trees.

Node IDs are mutable-tree IDs: pruning can release them, and later operations
can reuse the same integers for different nodes.)doc")
        .def(py::init<std::shared_ptr<DynamicComponentTree>, std::shared_ptr<DynamicComponentTree>>(),
             py::arg("minTree"),
             py::arg("maxTree"),
             py::keep_alive<1, 2>(),
             py::keep_alive<1, 3>(),
             R"doc(Create an adjuster for an existing paired tree representation.

``minTree`` must be a min-tree and ``maxTree`` a max-tree for the same image,
dimensions, and adjacency relation.  Both must provide adjacency information.
The filter keeps both objects alive.  Invalid null-like inputs raise
``RuntimeError``; the semantic pairing requirements are preconditions.)doc")
        .def("refreshAttributes",
             &PyDualMinMaxTreeIncrementalFilter::refreshAttributes,
             R"doc(Recompute and re-register both slot-indexed area buffers.

Call this after compatible direct tree mutations that did not pass through this
filter.  Existing ``minArea`` or ``maxArea`` Python lists are snapshots and are
not updated in place.  This operation does not repair a mismatched tree pair.)doc")
        .def("updateTree",
             &PyDualMinMaxTreeIncrementalFilter::updateTree,
             py::arg("tree"),
             py::arg("subtreeRoot"),
             R"doc(Update one dual tree before pruning a subtree from the opposite tree.

``tree`` is the dual hierarchy to mutate and must be exactly ``minTree`` or
``maxTree``.  ``subtreeRoot`` identifies a live non-root subtree in the
*opposite* (primal) tree.  The primal subtree must still be alive because this
method reads its pixels; this method does not prune it.  Passing another tree
raises ``RuntimeError``.  Other preconditions are enforced only by C++ debug
assertions.)doc")
        .def("pruneMaxTreeAndUpdateMinTree",
             &PyDualMinMaxTreeIncrementalFilter::pruneMaxTreeAndUpdateMinTree,
             py::arg("nodesToPrune"),
             R"doc(Update the min-tree and then prune max-tree subtrees in order.

``nodesToPrune`` is copied from the Python iterable.  Each currently live,
non-root max-tree ID is processed sequentially; invalid, root, or already dead
IDs are skipped.  Pruned IDs may subsequently be reused.)doc")
        .def("pruneMinTreeAndUpdateMaxTree",
             &PyDualMinMaxTreeIncrementalFilter::pruneMinTreeAndUpdateMaxTree,
             py::arg("nodesToPrune"),
             R"doc(Update the max-tree and then prune min-tree subtrees in order.

``nodesToPrune`` is copied from the Python iterable.  Each currently live,
non-root min-tree ID is processed sequentially; invalid, root, or already dead
IDs are skipped.  Pruned IDs may subsequently be reused.)doc")
        .def_property_readonly(
            "minTree",
            &PyDualMinMaxTreeIncrementalFilter::getMinTree,
            "The same mutable min-tree object retained by the filter.")
        .def_property_readonly(
            "maxTree",
            &PyDualMinMaxTreeIncrementalFilter::getMaxTree,
            "The same mutable max-tree object retained by the filter.")
        .def_property_readonly(
            "minArea",
            &PyDualMinMaxTreeIncrementalFilter::getMinArea,
            R"doc(A new slot-indexed area snapshot for the min-tree.

The Python list is indexed by node ID; released slots are reported as zero.
The list is not live and must be reacquired after mutations.)doc")
        .def_property_readonly(
            "maxArea",
            &PyDualMinMaxTreeIncrementalFilter::getMaxArea,
            R"doc(A new slot-indexed area snapshot for the max-tree.

The Python list is indexed by node ID; released slots are reported as zero.
The list is not live and must be reacquired after mutations.)doc");
}

void init_component_tree_casf(py::module_ &m) {
    py::class_<PyComponentTreeCasf, std::shared_ptr<PyComponentTreeCasf>>(
        m,
        "ComponentTreeCasf",
        R"doc(Stateful CASF filtering for order-preserving component-tree attributes.

Use this class with ``area``, ``bbox_width``, ``bbox_height``, or
``bbox_diagonal``.  Non-order-preserving descriptors such as rectangularity
belong to ``ComponentTreeCasfNonIncreasing``.  Thresholds mutate one paired
min-tree/max-tree state sequentially; construct a new object for an independent
run.

The ``updating`` mode mutates the existing tree objects.  A non-empty ``naive``
run, or the first step of a non-empty ``hybrid`` run, replaces those objects.
Consequently, previously obtained ``minTree``/``maxTree`` references become
invalid and must never be accessed after such a call.)doc")
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         const std::string &attribute,
                         double radiusAdj) {
                return std::make_shared<PyComponentTreeCasf>(input, attribute, radiusAdj);
            }),
             py::arg("image"),
             py::arg("attribute") = "area",
             py::arg("radiusAdj") = 1.5,
             R"doc(Create a conventional CASF filter with circular adjacency.

``image`` must be two-dimensional.  It is force-cast to C-contiguous ``uint8``
and copied, so the filter never aliases the input array.  ``attribute`` should
be ``area``, ``bbox_width``, ``bbox_height``, or ``bbox_diagonal``; common
hyphenated and short aliases are accepted.  ``radiusAdj`` must be non-negative.
An unknown attribute or invalid image/adjacency raises ``RuntimeError`` or
``ValueError``.)doc")
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         const std::string &attribute,
                         const std::shared_ptr<AdjacencyRelation> &adj) {
                return std::make_shared<PyComponentTreeCasf>(input, attribute, adj);
            }),
             py::arg("image"),
             py::arg("attribute") = "area",
             py::arg("adj"),
             R"doc(Create a conventional CASF filter with an existing adjacency.

``image`` is force-cast to C-contiguous ``uint8`` and copied.  ``adj`` is kept
alive by shared ownership and must match the image dimensions.  Use an
order-preserving attribute: ``area``, ``bbox_width``, ``bbox_height``, or
``bbox_diagonal``.  Unknown attributes and null inputs raise ``RuntimeError``;
the matching-domain requirement is a caller precondition and is not validated.)doc")
        .def("filter",
             &PyComponentTreeCasf::filter,
             py::arg("thresholds"),
             py::arg("mode") = "updating",
             R"doc(Apply integer thresholds sequentially and reconstruct the result.

``mode`` is case-insensitive and accepts ``"updating"``, ``"naive"``, or
``"hybrid"``.  Updating performs every step incrementally.  Naive rebuilds the
pair after every threshold.  Hybrid rebuilds for the first threshold of each
non-empty call and then continues incrementally.  An unknown mode raises
``RuntimeError``.

This method is stateful: each threshold starts from the previous current state.
Within each tree direction, traversal stops at the first node on a path whose
attribute is at or below the threshold.  If the root itself qualifies, that
direction performs no pruning because the root is never removed.

It returns a newly allocated, writable 2D ``uint8`` NumPy array that does not
alias the filter.  Pruning may release and reuse node IDs.  If a naive step
runs (non-empty naive/hybrid input), all tree references obtained before this
call are invalid because their underlying C++ objects were replaced.  The call
is not transactional: an exception does not roll back completed thresholds,
tree edits, or rebuilds.)doc")
        .def("getMinTree",
             &PyComponentTreeCasf::getMinTree,
             R"doc(Return the mutable min-tree owned by this filter.

The returned alias keeps the CASF wrapper alive.  Direct mutation can invalidate
the filter's attribute state and paired-tree invariants.  The reference becomes
invalid after any call that executes a naive rebuild; reacquire it afterwards.)doc")
        .def("getMaxTree",
             &PyComponentTreeCasf::getMaxTree,
             R"doc(Return the mutable max-tree owned by this filter.

The returned alias keeps the CASF wrapper alive.  Direct mutation can invalidate
the filter's attribute state and paired-tree invariants.  The reference becomes
invalid after any call that executes a naive rebuild; reacquire it afterwards.)doc")
        .def_property_readonly(
            "minTree",
            &PyComponentTreeCasf::getMinTree,
            R"doc(The current mutable min-tree.

Reacquire this property after a non-empty naive or hybrid filter call; an older
reference points to an object that was replaced.)doc")
        .def_property_readonly(
            "maxTree",
            &PyComponentTreeCasf::getMaxTree,
            R"doc(The current mutable max-tree.

Reacquire this property after a non-empty naive or hybrid filter call; an older
reference points to an object that was replaced.)doc");
}

void init_component_tree_casf_non_increasing(py::module_ &m) {
    py::class_<PyComponentTreeCasfNonIncreasing, std::shared_ptr<PyComponentTreeCasfNonIncreasing>>(
        m,
        "ComponentTreeCasfNonIncreasing",
        R"doc(Stateful CASF filtering for non-order-preserving flat attributes.

Supported attributes are ``area_per_box_perimeter``, ``rectangularity``,
``box_aspect_ratio``, ``box_elongation``, and
``area_per_box_diagonal_squared``.  At the start of each max-tree or min-tree
phase, nodes whose snapshotted value is at or below the floating-point threshold
become candidates.  Surviving candidates are not threshold-tested again after
earlier merges in that phase.  Local contractions keep the paired tree and
attribute buffers synchronized incrementally.

Filtering mutates the same tree objects, so existing tree references remain
object-valid.  Their node IDs are not stable: merges release slots that later
operations may reuse.)doc")
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         const std::string &attribute,
                         double radiusAdj) {
                return std::make_shared<PyComponentTreeCasfNonIncreasing>(input, attribute, radiusAdj);
            }),
             py::arg("image"),
             py::arg("attribute") = "area_per_box_perimeter",
             py::arg("radiusAdj") = 1.5,
             R"doc(Create a non-increasing CASF filter with circular adjacency.

``image`` must be two-dimensional and is force-cast to C-contiguous ``uint8``
then copied.  The filter does not alias subsequent input-array changes.
``attribute`` must be one of this class's five supported descriptors;
hyphenated and documented aliases are accepted.  ``radiusAdj`` must be
non-negative.  Unsupported attributes and invalid inputs raise ``RuntimeError``
or ``ValueError``.)doc")
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         const std::string &attribute,
                         const std::shared_ptr<AdjacencyRelation> &adj) {
                return std::make_shared<PyComponentTreeCasfNonIncreasing>(input, attribute, adj);
            }),
             py::arg("image"),
             py::arg("attribute") = "area_per_box_perimeter",
             py::arg("adj"),
             R"doc(Create a non-increasing CASF filter with an existing adjacency.

``image`` is force-cast to C-contiguous ``uint8`` and copied.  ``adj`` is kept
alive by shared ownership and must match the image dimensions.  ``attribute``
must be ``area_per_box_perimeter``, ``rectangularity``,
``box_aspect_ratio``, ``box_elongation``, or
``area_per_box_diagonal_squared``.  Null or unsupported inputs raise
``RuntimeError``; the matching-domain requirement is an unchecked precondition.)doc")
        .def("filter",
             &PyComponentTreeCasfNonIncreasing::filter,
             py::arg("thresholds"),
             R"doc(Apply floating-point thresholds sequentially and reconstruct.

For each threshold the max-tree phase runs first, followed by the min-tree
phase.  Every value continues from the state left by the previous value, so use
a new object for an independent run.  A node becomes a merge candidate when
its value in that phase's snapshot is at or below the threshold; ``NaN`` is not
selected.  A surviving candidate is not threshold-tested again after an earlier
merge in the same phase changes its attribute.

The returned value is a newly allocated, writable 2D ``uint8`` NumPy array and
does not alias the filter.  Tree objects remain the same, but merges can release
and later reuse node IDs.  Filtering is not transactional and does not roll
back contractions completed before an exception.)doc")
        .def("getMinTree",
             &PyComponentTreeCasfNonIncreasing::getMinTree,
             R"doc(Return the current mutable min-tree owned by this filter.

The alias keeps the CASF wrapper alive and remains object-valid across filter
calls.  Node IDs can still be released or reused.  Direct mutation can
invalidate incremental attributes and the paired-tree invariants.)doc")
        .def("getMaxTree",
             &PyComponentTreeCasfNonIncreasing::getMaxTree,
             R"doc(Return the current mutable max-tree owned by this filter.

The alias keeps the CASF wrapper alive and remains object-valid across filter
calls.  Node IDs can still be released or reused.  Direct mutation can
invalidate incremental attributes and the paired-tree invariants.)doc")
        .def_property_readonly(
            "minTree",
            &PyComponentTreeCasfNonIncreasing::getMinTree,
            "The current mutable min-tree; node IDs may change after filtering.")
        .def_property_readonly(
            "maxTree",
            &PyComponentTreeCasfNonIncreasing::getMaxTree,
            "The current mutable max-tree; node IDs may change after filtering.");
}

}  // namespace

PYBIND11_MODULE(morphoTreeAdjust, m) {
    m.doc() = R"doc(MorphoTreeAdjust dynamic component-tree bindings.

The module exposes mutable min/max trees, stateful adjacency iteration, paired
incremental adjustment, and CASF filters.  Image inputs are force-cast to
C-contiguous uint8 and copied into C++ storage.  Reconstructed NumPy arrays own
their returned image buffers and do not provide a live view of tree state.

Node IDs belong to mutable slot storage and may be released and reused.  Treat
IDs and returned Python lists as snapshots and revalidate them after every
mutation.)doc";

    init_adjacency_relation(m);
    init_dynamic_component_tree(m);
    init_dual_min_max_tree_incremental_filter(m);
    init_component_tree_casf(m);
    init_component_tree_casf_non_increasing(m);
}
