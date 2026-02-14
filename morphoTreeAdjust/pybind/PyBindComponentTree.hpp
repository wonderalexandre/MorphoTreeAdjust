#pragma once

#include <map>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"

#include "../pybind/PyBindFlatZonesGraph.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

template<typename CNPsType>
class PyBindComponentTree: public ComponentTree<CNPsType> {

public:

	
    PyBindComponentTree(py::array_t<uint8_t> &input, int numRows, int numCols, bool isMaxtree);
    PyBindComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radiusAdj);
    PyBindComponentTree(std::shared_ptr<PyBindFlatZonesGraph> graph, bool isMaxtree);

	py::array_t<uint8_t> reconstructionImage();

    py::array_t<uint8_t> reconstructionNode(NodeCT<CNPsType> node);

    std::map<int, NodeCT<CNPsType>> getNodes();

    // Nota: remoção do utilitário estático recNode; use reconstructionNode
    
};

#include "../pybind/PyBindComponentTree.tpp"
