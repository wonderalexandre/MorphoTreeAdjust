#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustmentBySubtree.hpp"
#include "../pybind/PyBindComponentTree.hpp"
#include "../include/Common.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_SUBTREE_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_SUBTREE_H

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;

class PyBindComponentTreeAdjustmentBySubtree: public ComponentTreeAdjustmentBySubtree {

public:
   
   PyBindComponentTreeAdjustmentBySubtree(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : ComponentTreeAdjustmentBySubtree(maxtree, mintree) {} 

   using ComponentTreeAdjustmentBySubtree::updateTree;
   void updateTree(PyBindComponentTreeFZPtr tree, NodeFZPtr L_leaf);

   py::tuple buildCollections(PyBindComponentTreeFZPtr tree, std::vector<int> flatZone, int newGrayLevel, bool isMaxtree);

};

#endif