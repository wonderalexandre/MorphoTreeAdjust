#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustmentByFlatzone.hpp"
#include "../pybind/PyBindComponentTree.hpp"
#include "../include/Common.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_FLATZONE_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_FLATZONE_H

namespace py = pybind11;

using PyBindComponentTreeFZ = PyBindComponentTree<FlatZones>;
using PyBindComponentTreeFZPtr = std::shared_ptr<PyBindComponentTreeFZ>;

class PyBindComponentTreeAdjustmentByFlatzone: public ComponentTreeAdjustmentByFlatzone {

public:
   
   PyBindComponentTreeAdjustmentByFlatzone(PyBindComponentTreeFZPtr maxtree, PyBindComponentTreeFZPtr mintree)
        : ComponentTreeAdjustmentByFlatzone(maxtree, mintree) {} 

   using ComponentTreeAdjustmentByFlatzone::updateTree;
   void updateTree(PyBindComponentTreeFZPtr tree, int repFlatzone);

   py::tuple buildCollections(PyBindComponentTreeFZPtr tree, int repFlatzone, int newGrayLevel);

};

#endif
