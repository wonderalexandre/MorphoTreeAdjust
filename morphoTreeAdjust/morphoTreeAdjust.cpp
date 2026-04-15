#include "include/AdjacencyRelation.hpp"
#include "include/AttributeComputer.hpp"
#include "include/Common.hpp"
#include "include/ComponentTreeCasf.hpp"
#include "include/DynamicComponentTree.hpp"
#include "include/DynamicComponentTreeAdjustment.hpp"

#include <memory>
#include <span>
#include <sstream>
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
    if (normalized == "box_width" || normalized == "box-width" || normalized == "width") {
        return BOX_WIDTH;
    }
    if (normalized == "box_height" || normalized == "box-height" || normalized == "height") {
        return BOX_HEIGHT;
    }
    if (normalized == "diagonal_length" || normalized == "diagonal-length" || normalized == "diagonal") {
        return DIAGONAL_LENGTH;
    }

    throw std::runtime_error("Unknown attribute. Expected one of: area, box_width, box_height, diagonal_length.");
}

std::vector<NodeId> alive_nodes(const DynamicComponentTree &tree) {
    std::vector<NodeId> nodes;
    nodes.reserve(static_cast<size_t>(tree.getNumNodes()));
    for (NodeId nodeId = 0; nodeId < tree.getGlobalIdSpaceSize(); ++nodeId) {
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

class PyDynamicComponentTreeAdjustment {
private:
    std::shared_ptr<DynamicComponentTree> mintree_;
    std::shared_ptr<DynamicComponentTree> maxtree_;
    DynamicAreaComputer minAreaComputer_;
    DynamicAreaComputer maxAreaComputer_;
    std::vector<float> minArea_;
    std::vector<float> maxArea_;
    DynamicComponentTreeAdjustment<AltitudeType> adjust_;

    void refreshAreaBuffers() {
        minArea_.assign(static_cast<size_t>(mintree_->getGlobalIdSpaceSize()), 0.0f);
        maxArea_.assign(static_cast<size_t>(maxtree_->getGlobalIdSpaceSize()), 0.0f);
        minAreaComputer_.compute(std::span<float>(minArea_));
        maxAreaComputer_.compute(std::span<float>(maxArea_));
        adjust_.setAttributeComputer(minAreaComputer_, maxAreaComputer_, std::span<float>(minArea_), std::span<float>(maxArea_));
    }

public:
    PyDynamicComponentTreeAdjustment(std::shared_ptr<DynamicComponentTree> mintree,
                                     std::shared_ptr<DynamicComponentTree> maxtree)
        : mintree_(std::move(mintree)),
          maxtree_(std::move(maxtree)),
          minAreaComputer_(mintree_.get()),
          maxAreaComputer_(maxtree_.get()),
          adjust_(mintree_.get(), maxtree_.get(), *mintree_->getAdjacencyRelation()) {
        if (mintree_ == nullptr || maxtree_ == nullptr) {
            throw std::runtime_error("DynamicComponentTreeAdjustment requires valid min-tree and max-tree.");
        }
        if (mintree_->getAdjacencyRelation() == nullptr || maxtree_->getAdjacencyRelation() == nullptr) {
            throw std::runtime_error("DynamicComponentTreeAdjustment requires trees with adjacency information.");
        }
        refreshAreaBuffers();
    }

    void refreshAttributes() {
        refreshAreaBuffers();
    }

    void updateTree(const std::shared_ptr<DynamicComponentTree> &tree, NodeId subtreeRoot) {
        adjust_.updateTree(tree.get(), subtreeRoot);
        refreshAreaBuffers();
    }

    void pruneMaxTreeAndUpdateMinTree(std::vector<NodeId> nodesToPrune) {
        adjust_.pruneMaxTreeAndUpdateMinTree(nodesToPrune);
        refreshAreaBuffers();
    }

    void pruneMinTreeAndUpdateMaxTree(std::vector<NodeId> nodesToPrune) {
        adjust_.pruneMinTreeAndUpdateMaxTree(nodesToPrune);
        refreshAreaBuffers();
    }

    std::vector<float> getMinArea() const { return minArea_; }
    std::vector<float> getMaxArea() const { return maxArea_; }
    std::shared_ptr<DynamicComponentTree> getMinTree() const { return mintree_; }
    std::shared_ptr<DynamicComponentTree> getMaxTree() const { return maxtree_; }
    std::string getOutputLog() const { return adjust_.getOutputLog(); }
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

void init_adjacency_relation(py::module_ &m) {
    py::class_<AdjacencyRelation, std::shared_ptr<AdjacencyRelation>>(m, "AdjacencyRelation")
        .def(py::init<int, int, double>())
        .def(py::init<int, int, const std::vector<AdjacencyRelation::Offset> &>())
        .def_static("rectangular", [](int numRows, int numCols, int halfHeight, int halfWidth) {
            return std::make_shared<AdjacencyRelation>(
                AdjacencyRelation::rectangular(numRows, numCols, halfHeight, halfWidth)
            );
        }, py::arg("numRows"), py::arg("numCols"), py::arg("halfHeight"), py::arg("halfWidth"))
        .def_property_readonly("size", &AdjacencyRelation::getSize)
        .def("getAdjPixels", py::overload_cast<int, int>(&AdjacencyRelation::getAdjPixels), py::return_value_policy::reference)
        .def("getAdjPixels", py::overload_cast<int>(&AdjacencyRelation::getAdjPixels), py::return_value_policy::reference)
        .def("getOffsets", [](AdjacencyRelation &self) {
            py::list out;
            for (int i = 0; i < self.getSize(); ++i) {
                out.append(py::make_tuple(self.getOffsetRow(i), self.getOffsetCol(i)));
            }
            return out;
        })
        .def("__iter__", [](AdjacencyRelation *adj) {
            return py::make_iterator(adj->begin(), adj->end());
        }, py::keep_alive<0, 1>());
}

void init_dynamic_component_tree(py::module_ &m) {
    py::class_<DynamicComponentTree, std::shared_ptr<DynamicComponentTree>> tree(m, "DynamicComponentTree");

    tree.def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         bool isMaxtree,
                         double radiusAdj) {
            auto img = image_from_numpy(input);
            auto adj = std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radiusAdj);
            return std::make_shared<DynamicComponentTree>(img, isMaxtree, adj);
        }),
        py::arg("image"),
        py::arg("isMaxtree"),
        py::arg("radiusAdj") = 1.5)
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         bool isMaxtree,
                         const std::shared_ptr<AdjacencyRelation> &adj) {
            auto image = image_from_numpy(input);
            return std::make_shared<DynamicComponentTree>(image, isMaxtree, adj);
        }),
        py::arg("image"),
        py::arg("isMaxtree"),
        py::arg("adj"))
        .def("reconstructionImage", [](const DynamicComponentTree &self) {
            return numpy_from_image(self.reconstructionImage());
        })
        .def("getNodesThreshold", [](DynamicComponentTree &self, int threshold) {
            return DynamicComponentTree::getNodesThreshold(&self, threshold);
        })
        .def("computeArea", [](DynamicComponentTree &self) {
            DynamicAreaComputer areaComputer(&self);
            return areaComputer.compute();
        })
        .def("getPixelsOfCC", &DynamicComponentTree::getPixelsOfCC)
        .def("getChildren", [](const DynamicComponentTree &self, NodeId nodeId) {
            return children_of(self, nodeId);
        })
        .def("getProperParts", [](const DynamicComponentTree &self, NodeId nodeId) {
            return proper_parts_of(self, nodeId);
        })
        .def("getNodeSubtree", [](const DynamicComponentTree &self, NodeId nodeId) {
            return subtree_nodes_of(self, nodeId);
        })
        .def("breadthFirstTraversal", [](const DynamicComponentTree &self) {
            return breadth_first_nodes_of(self);
        })
        .def("moveProperPart", &DynamicComponentTree::moveProperPart)
        .def("moveProperParts", &DynamicComponentTree::moveProperParts)
        .def("moveChildren", &DynamicComponentTree::moveChildren)
        .def("attachNode", &DynamicComponentTree::attachNode)
        .def("detachNode", &DynamicComponentTree::detachNode)
        .def("removeChild", &DynamicComponentTree::removeChild)
        .def_property_readonly("nodes", [](const DynamicComponentTree &self) {
            return alive_nodes(self);
        })
        .def_property_readonly("numNodes", &DynamicComponentTree::getNumNodes)
        .def_property_readonly("root", &DynamicComponentTree::getRoot)
        .def_property_readonly("isMaxtree", &DynamicComponentTree::isMaxtree)
        .def_property_readonly("numRows", &DynamicComponentTree::getNumRowsOfImage)
        .def_property_readonly("numCols", &DynamicComponentTree::getNumColsOfImage)
        .def("getAltitude", &DynamicComponentTree::getAltitude)
        .def("getNodeParent", &DynamicComponentTree::getNodeParent)
        .def("getNumChildren", &DynamicComponentTree::getNumChildren)
        .def("getNumProperParts", &DynamicComponentTree::getNumProperParts)
        .def("getSmallestComponent", &DynamicComponentTree::getSmallestComponent)
        .def("isAlive", &DynamicComponentTree::isAlive)
        .def("isNode", &DynamicComponentTree::isNode)
        .def("isLeaf", &DynamicComponentTree::isLeaf)
        .def("hasChild", &DynamicComponentTree::hasChild);

}

void init_dynamic_component_tree_adjustment(py::module_ &m) {
    py::class_<PyDynamicComponentTreeAdjustment>(m, "DynamicComponentTreeAdjustment")
        .def(py::init<std::shared_ptr<DynamicComponentTree>, std::shared_ptr<DynamicComponentTree>>(),
             py::keep_alive<1, 2>(),
             py::keep_alive<1, 3>())
        .def("refreshAttributes", &PyDynamicComponentTreeAdjustment::refreshAttributes)
        .def("updateTree", &PyDynamicComponentTreeAdjustment::updateTree)
        .def("pruneMaxTreeAndUpdateMinTree", &PyDynamicComponentTreeAdjustment::pruneMaxTreeAndUpdateMinTree)
        .def("pruneMinTreeAndUpdateMaxTree", &PyDynamicComponentTreeAdjustment::pruneMinTreeAndUpdateMaxTree)
        .def_property_readonly("minTree", &PyDynamicComponentTreeAdjustment::getMinTree)
        .def_property_readonly("maxTree", &PyDynamicComponentTreeAdjustment::getMaxTree)
        .def_property_readonly("minArea", &PyDynamicComponentTreeAdjustment::getMinArea)
        .def_property_readonly("maxArea", &PyDynamicComponentTreeAdjustment::getMaxArea)
        .def("log", &PyDynamicComponentTreeAdjustment::getOutputLog);
}

void init_component_tree_casf(py::module_ &m) {
    py::class_<PyComponentTreeCasf, std::shared_ptr<PyComponentTreeCasf>>(m, "ComponentTreeCasf")
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         const std::string &attribute,
                         double radiusAdj) {
                return std::make_shared<PyComponentTreeCasf>(input, attribute, radiusAdj);
            }),
             py::arg("image"),
             py::arg("attribute") = "area",
             py::arg("radiusAdj") = 1.5)
        .def(py::init([](const py::array_t<uint8_t, py::array::c_style | py::array::forcecast> &input,
                         const std::string &attribute,
                         const std::shared_ptr<AdjacencyRelation> &adj) {
                return std::make_shared<PyComponentTreeCasf>(input, attribute, adj);
            }),
             py::arg("image"),
             py::arg("attribute") = "area",
             py::arg("adj"))
        .def("filter", &PyComponentTreeCasf::filter, py::arg("thresholds"), py::arg("mode") = "updating")
        .def("getMinTree", &PyComponentTreeCasf::getMinTree)
        .def("getMaxTree", &PyComponentTreeCasf::getMaxTree)
        .def_property_readonly("minTree", &PyComponentTreeCasf::getMinTree)
        .def_property_readonly("maxTree", &PyComponentTreeCasf::getMaxTree);
}

}  // namespace

PYBIND11_MODULE(morphoTreeAdjust, m) {
    m.doc() = "MorphoTreeAdjust core dynamic C++/Python bindings.";

    init_adjacency_relation(m);
    init_dynamic_component_tree(m);
    init_dynamic_component_tree_adjustment(m);
    init_component_tree_casf(m);
}
