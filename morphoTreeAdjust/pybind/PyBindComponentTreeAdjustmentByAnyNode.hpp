#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustmentByAnyNode.hpp"
#include "../pybind/PyBindComponentTree.hpp"
#include "../include/Common.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_ANY_NODE_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_ANY_NODE_H

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;

class PyBindComponentTreeAdjustmentByAnyNode: public ComponentTreeAdjustmentByAnyNode {

public:
   
   PyBindComponentTreeAdjustmentByAnyNode(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : ComponentTreeAdjustmentByAnyNode(maxtree, mintree) {} 

   using ComponentTreeAdjustmentByAnyNode::updateTree;
   void updateTree(PyBindComponentTreeFZPtr tree, NodeFZ node);

   //py::tuple buildCollections(PyBindComponentTreeFZPtr tree, NodeFZPtr leaf, int newGrayLevel);

};

#endif