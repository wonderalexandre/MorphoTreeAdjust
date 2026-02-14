#pragma once

#include "../include/ComponentTreeAdjustmentBySubtree.hpp"
#include "../pybind/PyBindComponentTree.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;
using PyBindAdjustmentSubtreeBase = ComponentTreeAdjustmentBySubtree<>;

class PyBindComponentTreeAdjustmentBySubtree: public PyBindAdjustmentSubtreeBase {

public:
   
   PyBindComponentTreeAdjustmentBySubtree(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : PyBindAdjustmentSubtreeBase(mintree.get(), maxtree.get()) {} 

   using PyBindAdjustmentSubtreeBase::updateTree;
   void updateTree(PyBindComponentTreeFZPtr tree, NodeFZ<> L_leaf);

   py::tuple buildCollections(PyBindComponentTreeFZPtr tree, NodeFZ<> nodeSubtree, int newGrayLevel);

};
