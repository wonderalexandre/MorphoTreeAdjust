#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"
#include "../include/FlatZonesGraph.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>


#ifndef PYBIND_FLATZONES_GRAPH_H
#define PYBIND_FLATZONES_GRAPH_H



namespace py = pybind11;

class PyBindFlatZonesGraph: public FlatZonesGraph {

public:
	
    PyBindFlatZonesGraph(py::array_t<uint8_t> &input, int numRows, int numCols, double radiusAdj);
    PyBindFlatZonesGraph(ImageUInt8Ptr img, double radiusAdj);
	
};

#endif