#include "include/NodeCT.hpp"
#include "pybind/PyBindComponentTree.hpp"
#include "include/AdjacencyRelation.hpp"
#include "pybind/PyBindComponentTreeAdjustmentByLeaf.hpp"
#include "pybind/PyBindComponentTreeAdjustmentBySubtree.hpp"

#include <iterator>
#include <utility>
#include <sstream>

namespace py = pybind11;

void init_NodeCT(py::module &m) {
    py::class_<NodeCT<FlatZones>, std::shared_ptr<NodeCT<FlatZones>>>(m, "NodeCT")
        .def(py::init<>())
        .def_property_readonly("id", &NodeCT<FlatZones>::getIndex)
        .def("__str__", [](NodeCT<FlatZones> &node) {
            std::ostringstream oss;
            oss << "NodeCT(id=" << node.getIndex() 
                << ", level=" << node.getLevel() 
                << ", numCNPs=" << node.getNumCNPs() 
                << ", area=" << node.getArea() 
                << ", threshold1=" << node.getThreshold1() 
                << ", threshold2=" << node.getThreshold2() << ")";
            return oss.str();
        })
        .def("__repr__", [](const NodeCT<FlatZones> &node) { 
            std::ostringstream oss;
            oss << "NodeCT(id=" << node.getIndex() << ", level=" << node.getLevel() << ")";
            return oss.str();
        })
        .def_property_readonly("cnps", [](NodeCT<FlatZones> &node) {
            return node.getCNPs();  // Retorna o iterador
        }, py::return_value_policy::reference_internal)

        .def_property_readonly("level", &NodeCT<FlatZones>::getLevel)
        .def_property_readonly("children", &NodeCT<FlatZones>::getChildren)
        .def_property_readonly("parent", &NodeCT<FlatZones>::getParent)
        .def_property_readonly("numSiblings", &NodeCT<FlatZones>::getNumSiblings)
        .def_property_readonly("numCNPs", &NodeCT<FlatZones>::getNumCNPs)
        .def_property_readonly("area", &NodeCT<FlatZones>::getArea)
        .def_property_readonly("threshold1", &NodeCT<FlatZones>::getThreshold1)
        .def_property_readonly("threshold2", &NodeCT<FlatZones>::getThreshold2)
        .def("pixelsOfCC", &NodeCT<FlatZones>::getPixelsOfCC)
        .def("nodesOfPathToRoot", &NodeCT<FlatZones>::getNodesOfPathToRoot)
        .def("postOrderTraversal", &NodeCT<FlatZones>::getIteratorPostOrderTraversal)
        .def("BFSTraversal", &NodeCT<FlatZones>::getIteratorBreadthFirstTraversal)
        .def("recNode", [](NodeCTPtr<FlatZones> node) {
            return PyBindComponentTree<FlatZones>::recNode(node);
        });
}

void init_NodeCT_Iterators(py::module &m) {
    using NodeType = NodeCT<FlatZones>;


    py::class_<typename NodeType::IteratorPixelsOfCC>(m, "IteratorPixelsOfCC")
        .def(py::init<std::shared_ptr<NodeType>, int>())
        .def("__iter__", [](typename NodeType::IteratorPixelsOfCC &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<typename NodeType::IteratorCNPs>(m, "IteratorCNPs")
        .def(py::init<std::shared_ptr<NodeType>>())  
        .def("__iter__", [](typename NodeType::IteratorCNPs &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>())  // Mantém os dados vivos enquanto o iterador estiver em uso
        .def("__len__", [](typename NodeType::IteratorCNPs &iter) {
            return std::distance(iter.begin(), iter.end()); // Retorna o tamanho do iterador
        })
        .def("__getitem__", [](typename NodeType::IteratorCNPs &iter, int index) {
            auto it = iter.begin();
            std::advance(it, index);
            return *it;
        });

        
    py::class_<typename NodeType::IteratorNodesOfPathToRoot>(m, "IteratorNodesOfPathToRoot")
        .def(py::init<std::shared_ptr<NodeType>>())
        .def("__iter__", [](typename NodeType::IteratorNodesOfPathToRoot &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<typename NodeType::IteratorPostOrderTraversal>(m, "IteratorPostOrderTraversal")
        .def(py::init<std::shared_ptr<NodeType>>())
        .def("__iter__", [](typename NodeType::IteratorPostOrderTraversal &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<typename NodeType::IteratorBreadthFirstTraversal>(m, "IteratorBreadthFirstTraversal")
        .def(py::init<std::shared_ptr<NodeType>>())
        .def("__iter__", [](typename NodeType::IteratorBreadthFirstTraversal &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());
}

void init_ComponentTree(py::module &m) {
    using PyBindTree = PyBindComponentTree<FlatZones>; 
    py::class_<PyBindTree, std::shared_ptr<PyBindTree>>(m, "ComponentTree")
        .def(py::init<py::array_t<uint8_t> &, int, int, bool>())
        //.def(py::init<py::array_t<uint8_t> &, int, int, bool>())
        .def("reconstructionImage", &PyBindTree::reconstructionImage)
        .def("recNode", &PyBindTree::reconstructionNode)
        .def("getSC", [](PyBindTree &self, int p) -> NodeFZPtr {
            NodeFZPtr result = self.getSC(p);
            if (!result) {
                throw py::value_error("getSC: Índice " + std::to_string(p) + " fora dos limites ou nó não inicializado.");
            }
            return result;
        }, "Obtém o nó SC correspondente ao índice p")
        .def("prunning", [](PyBindTree& tree, py::object node) {
            std::shared_ptr<NodeCT<FlatZones>> smart_node = node.cast<std::shared_ptr<NodeCT<FlatZones>>>();
            tree.prunning(smart_node);
            node = py::none();
        })
       /* .def("mergeWithParent", [](PyBindTree& tree, py::object node) {
            NodeCT<FlatZones>* raw_node = node.cast<NodeCT<FlatZones>*>();
            tree.mergeWithParent(raw_node);
            node = py::none();
        })*/
        .def("getNodesThreshold", &PyBindTree::getNodesThreshold)
        .def("leaves", &PyBindTree::getLeaves)
        .def("nodes", &PyBindTree::getNodes)
        .def_property_readonly("numNodes", &PyBindTree::getNumNodes)
        .def_property_readonly("root", &PyBindTree::getRoot)
        .def_property_readonly("isMaxtree", &PyBindTree::isMaxtree)
        .def_property_readonly("numRows", &PyBindTree::getNumRowsOfImage)
        .def_property_readonly("numCols", &PyBindTree::getNumColsOfImage);
}

void init_ComponentTreeAdjustment(py::module &m) {
    using PyBindTree = PyBindComponentTree<FlatZones>; 
    py::class_<PyBindComponentTreeAdjustmentByLeaf>(m, "ComponentTreeAdjustmentByLeaf")
        .def(py::init<std::shared_ptr<PyBindTree>, std::shared_ptr<PyBindTree>>())
        .def("updateTree", [](PyBindComponentTreeAdjustmentByLeaf &self, PyBindComponentTreeFZPtr tree, NodeFZPtr leaf) {
            self.updateTree(tree, leaf);
        })
        .def("buildCollections", &PyBindComponentTreeAdjustmentByLeaf::buildCollections)
        .def("log", &PyBindComponentTreeAdjustmentByLeaf::getOutputLog);

    py::class_<PyBindComponentTreeAdjustmentBySubtree>(m, "ComponentTreeAdjustmentBySubtree")
        .def(py::init<std::shared_ptr<PyBindTree>, std::shared_ptr<PyBindTree>>())
        .def("updateTree", [](PyBindComponentTreeAdjustmentBySubtree &self, PyBindComponentTreeFZPtr tree, NodeFZPtr leaf) {
            self.updateTree(tree, leaf);
        })
        .def("buildCollections", &PyBindComponentTreeAdjustmentBySubtree::buildCollections)
        .def("log", &PyBindComponentTreeAdjustmentBySubtree::getOutputLog);
}

void init_AdjacencyRelation(py::module &m) {
    py::class_<AdjacencyRelation>(m, "AdjacencyRelation")
        .def(py::init<int, int, double>())
        .def_property_readonly("size", &AdjacencyRelation::getSize)
        .def("getAdjPixels", py::overload_cast<int, int>(&AdjacencyRelation::getAdjPixels), py::return_value_policy::reference)
        .def("getAdjPixels", py::overload_cast<int>(&AdjacencyRelation::getAdjPixels), py::return_value_policy::reference)
        .def("__iter__", [](AdjacencyRelation* adj) {
            return py::make_iterator(adj->begin(), adj->end());
        }, py::keep_alive<0, 1>());
}

PYBIND11_MODULE(morphoTreeAdjust, m) {
    m.doc() = "MorphoTreeAdjust is a C++/Python implementation for adjusting the morphological trees.";

    init_AdjacencyRelation(m);
    init_NodeCT(m);
    init_NodeCT_Iterators(m);
    init_ComponentTree(m);
    init_ComponentTreeAdjustment(m);
}