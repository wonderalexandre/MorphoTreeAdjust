
#include "include/NodeCT.hpp"
#include "pybind/PyBindComponentTree.hpp"
#include "include/AdjacencyRelation.hpp"
#include "pybind/PyBindComponentTreeAdjustment.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iterator>
#include <utility>


namespace py = pybind11;


void init_NodeCT(py::module &m){
    py::class_<NodeCT>(m, "NodeCT")
        .def(py::init<>())
        .def_property_readonly("id", &NodeCT::getIndex )
        .def_property_readonly("cnps", &NodeCT::getCNPsToVector )
        .def_property_readonly("level", &NodeCT::getLevel )
        .def_property_readonly("children", &NodeCT::getChildrenToVector )
        .def_property_readonly("parent", &NodeCT::getParent )
        .def_property_readonly("numSiblings", &NodeCT::getNumSiblings )
        .def("nodesOfPathToRoot", &NodeCT::getNodesOfPathToRoot )
        .def("postOrderTraversal", &NodeCT::getIteratorPostOrderTraversal );

    // Configuração para `IteratorNodesOfPathToRoot`
    py::class_<NodeCT::IteratorNodesOfPathToRoot>(m, "IteratorNodesOfPathToRoot")
        .def(py::init<NodeCT *>())
        .def("__iter__", [](NodeCT::IteratorNodesOfPathToRoot &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<NodeCT::InternalIteratorNodesOfPathToRoot>(m, "InternalIteratorNodesOfPathToRoot")
        .def("__iter__", [](NodeCT::InternalIteratorNodesOfPathToRoot &self) -> NodeCT::InternalIteratorNodesOfPathToRoot& { return self; })
        .def("__next__", [](NodeCT::InternalIteratorNodesOfPathToRoot &self) -> NodeCT* {
            if (self == NodeCT::InternalIteratorNodesOfPathToRoot(nullptr))
                throw py::stop_iteration();
            NodeCT* current = *self;
            ++self;
            return current;
        });

    // Configuração para `IteratorPostOrderTraversal`
    py::class_<NodeCT::IteratorPostOrderTraversal>(m, "IteratorPostOrderTraversal")
        .def(py::init<NodeCT *>())
        .def("__iter__", [](NodeCT::IteratorPostOrderTraversal &iter) {
            return py::make_iterator(iter.begin(), iter.end());
        }, py::keep_alive<0, 1>());

    py::class_<NodeCT::InternalIteratorPostOrderTraversal>(m, "InternalIteratorPostOrderTraversal")
        .def("__iter__", [](NodeCT::InternalIteratorPostOrderTraversal &self) -> NodeCT::InternalIteratorPostOrderTraversal& { return self; })
        .def("__next__", [](NodeCT::InternalIteratorPostOrderTraversal &self) -> NodeCT* {
            if (self == NodeCT::InternalIteratorPostOrderTraversal(nullptr))
                throw py::stop_iteration();
            NodeCT* current = *self;
            ++self;
            return current;
        });
}




void init_ComponentTree(py::module &m){
      py::class_<PyBindComponentTree>(m, "ComponentTree")
        .def(py::init<py::array_t<int> &, int, int, bool, double>())
        .def(py::init<py::array_t<int> &, int, int, bool>())
        .def("reconstructionImage", &PyBindComponentTree::reconstructionImage )
        .def("recNode", &PyBindComponentTree::reconstructionNode )
        .def("getSC", &PyBindComponentTree::getSC )
        .def("prunning", [](PyBindComponentTree& tree, NodeCT* node) {
            return tree.prunning(node);
        }, py::arg("node"))
        .def("getNodesThreshold", &PyBindComponentTree::getNodesThreshold)
        .def("leaves", &PyBindComponentTree::getLeaves)
		.def_property_readonly("numNodes", &PyBindComponentTree::getNumNodes )
        .def_property_readonly("root", &PyBindComponentTree::getRoot );

}

void init_ComponentTreeAdjustment(py::module &m){
    py::class_<PyBindComponentTreeAdjustment>(m, "ComponentTreeAdjustment")
    .def(py::init<PyBindComponentTree&, PyBindComponentTree&>())
    .def("updateTree", &PyBindComponentTreeAdjustment::updateTree )
    .def("buildCollections", &PyBindComponentTreeAdjustment::buildCollections);
    



}


void init_AdjacencyRelation(py::module &m){
    	py::class_<AdjacencyRelation>(m, "AdjacencyRelation")
        .def(py::init<int, int, double>())
        .def_property_readonly("size", &AdjacencyRelation::getSize )
        .def("getAdjPixels", py::overload_cast<int, int>( &AdjacencyRelation::getAdjPixels ));

}

PYBIND11_MODULE(morphoTreeAdjust, m) {
    // Optional docstring
    m.doc() = "MorphoTreeAdjust is a C++/Python implementation for adjusting the morpholofical trees.";
    
    init_NodeCT(m);
    init_ComponentTree(m);
    init_AdjacencyRelation(m);
    init_ComponentTreeAdjustment(m);
}
