#include <list>
#include <vector>
#include <map>
#include <stack>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

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

	
    PyBindComponentTree(py::array_t<uint8_t> &input, int numRows, int numCols, bool isMaxtree);
    PyBindComponentTree(ImageUInt8Ptr img, int numRows, int numCols, bool isMaxtree, double radiusAdj);

	py::array_t<uint8_t> reconstructionImage();

    py::array_t<uint8_t> reconstructionNode(NodeCTPtr<CNPsType> node);

    std::map<int, NodeCTPtr<CNPsType>> getNodes();

    static py::array_t<uint8_t> recNode(NodeCTPtr<CNPsType> _node);
    
	
};

#include "../pybind/PyBindComponentTree.tpp"

#endif