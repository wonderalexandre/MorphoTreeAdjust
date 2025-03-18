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

class PyBindComponentTreeAdjustment: public ComponentTreeAdjustment {

public:
   
   PyBindComponentTreeAdjustment(PyBindComponentTreeFZ* maxtree, PyBindComponentTreeFZ* mintree)
        : ComponentTreeAdjustment(maxtree, mintree) {} 

   void updateTree(PyBindComponentTreeFZ* tree, NodeFZ* L_leaf);

   void updateTree2(PyBindComponentTreeFZ* tree, NodeFZ* rSubtree);

   py::tuple buildCollections(PyBindComponentTreeFZ* tree, std::vector<int> flatZone, int newGrayLevel, bool isMaxtree);

};

#endif