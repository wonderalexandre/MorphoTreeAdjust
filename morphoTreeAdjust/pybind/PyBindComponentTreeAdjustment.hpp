#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustment.hpp"
#include "../pybind/PyBindComponentTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_H

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;

class PyBindComponentTreeAdjustment: public ComponentTreeAdjustment {

public:
   
   PyBindComponentTreeAdjustment(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : ComponentTreeAdjustment(maxtree, mintree) {} 

   void updateTree(PyBindComponentTreeFZPtr tree, NodeFZPtr L_leaf);

   void updateTree2(PyBindComponentTreeFZPtr tree, NodeFZPtr rSubtree);

   void updateTree3(PyBindComponentTreeFZ* tree, NodeFZ* node);


   py::tuple buildCollections(PyBindComponentTreeFZPtr tree, std::vector<int> flatZone, int newGrayLevel, bool isMaxtree);

};

#endif