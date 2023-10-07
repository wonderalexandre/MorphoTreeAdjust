#include <list>
#include <vector>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_H
#define PYBIND_COMPONENT_TREE_H


namespace py = pybind11;


class PyBindComponentTree: public ComponentTree {

public:
   
    PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree);

	PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	py::array_t<int> reconstructionImage();

	py::array_t<int> reconstructionNode(NodeCT* node);

};

#endif