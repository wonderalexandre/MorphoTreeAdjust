#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustmentByLeaf.hpp"
#include "../pybind/PyBindComponentTree.hpp"
#include "../include/Common.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_LEAF_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_LEAF_H

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;

class PyBindComponentTreeAdjustmentByLeaf: public ComponentTreeAdjustmentByLeaf {

public:
   
   PyBindComponentTreeAdjustmentByLeaf(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : ComponentTreeAdjustmentByLeaf(maxtree, mintree) {} 

   using ComponentTreeAdjustmentByLeaf::updateTree;
   void updateTree(PyBindComponentTreeFZPtr tree, NodeFZ L_leaf);

   py::tuple buildCollections(PyBindComponentTreeFZPtr tree, NodeFZ leaf, int newGrayLevel);

};

#endif