#include "include/NodeCT.hpp"
#include "pybind/PyBindComponentTree.hpp"
#include "include/AdjacencyRelation.hpp"
#include "pybind/PyBindComponentTreeAdjustmentByLeaf.hpp"
#include "pybind/PyBindComponentTreeAdjustmentBySubtree.hpp"
#include "pybind/PyBindFlatZonesGraph.hpp"

#include <sstream>

namespace py = pybind11;

void init_NodeCT(py::module &m) {
    using NodeFZ = NodeCT<FlatZones>;
    py::class_<NodeFZ>(m, "NodeCT")
        .def(py::init<>())
        .def_property_readonly("id", &NodeFZ::getIndex)
        .def("__str__", [](NodeFZ &node) {
            std::ostringstream oss;
            oss << "NodeCT(id=" << node.getIndex() 
                << ", level=" << node.getLevel() 
                << ", parentID=" << node.getParent().getIndex() 
                << ", numCNPs=" << node.getNumCNPs()
                << ", numFlatzones=" << node.getNumFlatzone()
                << ", area=" << node.getArea() << ")";
            return oss.str();
        })
        .def("__repr__", [](NodeFZ &node) { 
            std::ostringstream oss;
            oss << "NodeCT(id=" << node.getIndex() 
                << ", level=" << node.getLevel() 
                << ", parentID=" << node.getParent().getIndex() 
                << ", numCNPs=" << node.getNumCNPs()
                << ", numFlatzones=" << node.getNumFlatzone()
                << ", area=" << node.getArea() << ")";
            return oss.str();
        })
        .def_property_readonly("repCNPs", [](NodeFZ &node) {
            py::list out;
            for (int r : node.getRepCNPs()) 
                out.append(r);
            return out;
        })
        .def_property_readonly("cnps", [](NodeFZ &node) {
            py::list out;
            for (int r : node.getCNPs()) 
                out.append(r);
            return out;
        })
        .def_property_readonly("pixelsOfCC", [](NodeFZ &n){
            py::list out;
            for (int p : n.getPixelsOfCC()) 
                out.append(p);
            return out;
        })
        .def_property_readonly("level", &NodeFZ::getLevel)
        .def_property_readonly("repNode", &NodeFZ::getRepNode)
        .def_property_readonly("parent", &NodeFZ::getParent)
        .def_property_readonly("numCNPs", &NodeFZ::getNumCNPs)
        .def_property_readonly("numFlatzones", &NodeFZ::getNumFlatzone)
        .def_property_readonly("area", &NodeFZ::getArea)
        .def_property_readonly("children", [](NodeFZ &node){
            std::vector<NodeFZ> ch;
            for (auto c : node.getChildren()) ch.push_back(c);
            return ch;
        })
        .def("__bool__", [](const NodeFZ& n){ return static_cast<bool>(n); })

        .def("nodesOfPathToRoot", [](NodeFZ &n){
            py::list out;
            for (auto x : n.getNodesOfPathToRoot()) out.append(x);
            return out;
        })
        .def("postOrderTraversal", [](NodeFZ &n){
            py::list out;
            for (auto x : n.getIteratorPostOrderTraversal()) out.append(x);
            return out;
        })
        .def("BFSTraversal", [](NodeFZ &n){
            py::list out;
            for (auto x : n.getIteratorBreadthFirstTraversal()) out.append(x);
            return out;
        });
}

// iterators de NodeCT agora são expostos diretamente via make_iterator nos métodos do NodeCT

void init_ComponentTree(py::module &m) {
    using PyBindTree = PyBindComponentTree<FlatZones>;  
    py::class_<PyBindTree, std::shared_ptr<PyBindTree>>(m, "ComponentTree")
        .def(py::init<std::shared_ptr<PyBindFlatZonesGraph>&, bool>())
        .def("reconstructionImage", &PyBindTree::reconstructionImage)
        .def("reconstructionNode", &PyBindTree::reconstructionNode)
        .def("getSC", [](PyBindTree &self, int p) -> NodeCT<FlatZones> {
            return self.getSC(p);
        }, "Obtém o nó SC correspondente ao índice p")
        /*
        .def("prunning", [](PyBindTree &self, NodeCT<FlatZones> node) {
            self.prunning(node.getIndex());
        })*/
        .def("getNodesThreshold", &PyBindTree::getNodesThreshold)
        .def("leaves", [](PyBindTree &self){
            py::list out;
            for (auto id : self.getLeaves()) out.append(self.proxy(id));
            return out;
        })
        .def("getPixelsOfFlatzone", [](PyBindTree &self, int repFZ) {
            py::list out;
            for (int p : self.getPixelsOfFlatzone(repFZ)) 
                out.append(p);
            return out;
        })
        .def_property_readonly("nodes", &PyBindTree::getNodes)
        .def_property_readonly("numNodes", &PyBindTree::getNumNodes)
        .def_property_readonly("root", &PyBindTree::getRoot)
        .def_property_readonly("isMaxtree", [](const PyBindTree& self){ return self.isMaxtree(); })
        .def_property_readonly("numRows", [](const PyBindTree& self){ return self.getNumRowsOfImage(); })
        .def_property_readonly("numCols", [](const PyBindTree& self){ return self.getNumColsOfImage(); });
}

void init_ComponentTreeAdjustment(py::module &m) {
    using PyBindTree = PyBindComponentTree<FlatZones>; 
    py::class_<PyBindComponentTreeAdjustmentByLeaf>(m, "ComponentTreeAdjustmentByLeaf")
        .def(py::init<std::shared_ptr<PyBindTree>, std::shared_ptr<PyBindTree>>())
        .def("updateTree", [](PyBindComponentTreeAdjustmentByLeaf &self, PyBindComponentTreeFZPtr tree, NodeCT<FlatZones> leaf) {
            self.updateTree(tree, leaf);
        })
        .def("mergeWithParent", [](PyBindComponentTreeAdjustmentByLeaf &self, PyBindComponentTreeFZPtr tree, NodeCT<FlatZones> node) {
            self.mergeWithParent(tree.get(), node);
        })
        .def("prunning", [](PyBindComponentTreeAdjustmentByLeaf &self, PyBindComponentTreeFZPtr tree, NodeCT<FlatZones> node) {
            self.prunning(tree.get(), node.getIndex());
        })
        .def("buildCollections", &PyBindComponentTreeAdjustmentByLeaf::buildCollections)
        .def("log", &PyBindComponentTreeAdjustmentByLeaf::getOutputLog);

    py::class_<PyBindComponentTreeAdjustmentBySubtree>(m, "ComponentTreeAdjustmentBySubtree")
        .def(py::init<std::shared_ptr<PyBindTree>, std::shared_ptr<PyBindTree>>())
        .def("updateTree", [](PyBindComponentTreeAdjustmentBySubtree &self, PyBindComponentTreeFZPtr tree, NodeCT<FlatZones> leaf) {
            self.updateTree(tree, leaf);
        })
        .def("mergeWithParent", [](PyBindComponentTreeAdjustmentBySubtree &self, PyBindComponentTreeFZPtr tree, NodeCT<FlatZones> node) {
            self.mergeWithParent(tree.get(), node);
        })
        .def("prunning", [](PyBindComponentTreeAdjustmentBySubtree &self, PyBindComponentTreeFZPtr tree, NodeCT<FlatZones> node) {
            self.prunning(tree.get(), node.getIndex());
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
        .def("getOffsets", [](AdjacencyRelation &self){
            py::list out;
            for (int i = 0; i < self.getSize(); ++i) {
                out.append(py::make_tuple(self.getOffsetRow(i), self.getOffsetCol(i)));
            }
            return out;
        })
        .def("__iter__", [](AdjacencyRelation* adj) {
            return py::make_iterator(adj->begin(), adj->end());
        }, py::keep_alive<0, 1>());
}

void init_FlatZonesGraph(py::module &m) {
    
    py::class_<PyBindFlatZonesGraph, std::shared_ptr<PyBindFlatZonesGraph>>(m, "FlatZonesGraph")
        .def(py::init<py::array_t<uint8_t> &, int, int, double>());
}


PYBIND11_MODULE(morphoTreeAdjust, m) {
    m.doc() = "MorphoTreeAdjust is a C++/Python implementation for adjusting the morphological trees.";

    init_AdjacencyRelation(m);
    init_NodeCT(m);
    init_ComponentTree(m);
    init_ComponentTreeAdjustment(m);
    init_FlatZonesGraph(m);
}
