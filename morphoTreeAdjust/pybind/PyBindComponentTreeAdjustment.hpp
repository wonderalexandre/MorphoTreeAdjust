#include <list>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustment.hpp"
#include "./PyBindComponentTree.hpp"


#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_H

namespace py = pybind11;

class PyBindComponentTreeAdjustment: public ComponentTreeAdjustment {

public:
   
   PyBindComponentTreeAdjustment(PyBindComponentTree* maxtree, PyBindComponentTree* mintree)
        : ComponentTreeAdjustment(maxtree, mintree) {} 

   void updateTree(PyBindComponentTree* tree, NodeCT *L_leaf);

   void updateTree2(PyBindComponentTree* tree, NodeCT *rSubtree);

   py::tuple buildCollections(PyBindComponentTree* tree, std::vector<int> flatZone, int newGrayLevel, bool isMaxtree);

};

#endif