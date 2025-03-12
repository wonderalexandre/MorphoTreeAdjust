#include <list>
#include <vector>
#include <map>
#include <stack>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>


#ifndef PYBIND_COMPONENT_TREE_H
#define PYBIND_COMPONENT_TREE_H

namespace py = pybind11;

template<typename CNPsType>
class PyBindComponentTree: public ComponentTree<CNPsType> {

public:

	PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree);
    PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);


	py::array_t<int> reconstructionImage();
    py::array_t<int> reconstructionNode(NodeCT<CNPsType>* node);

    NodeCT<CNPsType>* getNodeByIndex(int index);
    std::map<int, NodeCT<CNPsType>*> getNodes();
    bool isNodesInitialized();

    static py::array_t<int> recNode(NodeCT<CNPsType>* _node);
    
	
};

#include "PyBindComponentTree.tpp"

#endif