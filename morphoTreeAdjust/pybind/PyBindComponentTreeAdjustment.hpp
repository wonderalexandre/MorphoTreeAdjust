#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTreeAdjustment.hpp"
#include "./PyBindComponentTree.hpp"

#ifndef PYBIND_COMPONENT_TREE_ADJUSTMENT_H
#define PYBIND_COMPONENT_TREE_ADJUSTMENT_H


class PyBindComponentTreeAdjustment: public ComponentTreeAdjustment {

public:
   
   void adjustPyBindMinTree(PyBindComponentTree &mintree, NodeCT *Lmax);

};

#endif