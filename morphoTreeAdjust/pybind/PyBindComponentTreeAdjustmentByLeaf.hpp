#pragma once

#include "../include/ComponentTreeAdjustmentByLeaf.hpp"
#include "../pybind/PyBindComponentTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;
using PyBindAdjustmentLeafBase = ComponentTreeAdjustmentByLeaf<>;

class PyBindComponentTreeAdjustmentByLeaf: public PyBindAdjustmentLeafBase {

public:
   
   PyBindComponentTreeAdjustmentByLeaf(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : PyBindAdjustmentLeafBase(mintree.get(), maxtree.get()) {} 

   using PyBindAdjustmentLeafBase::updateTree;
   void updateTree(PyBindComponentTreeFZPtr tree, NodeFZ<> L_leaf);

   py::tuple buildCollections(PyBindComponentTreeFZPtr tree, NodeFZ<> leaf, int newGrayLevel);

};
