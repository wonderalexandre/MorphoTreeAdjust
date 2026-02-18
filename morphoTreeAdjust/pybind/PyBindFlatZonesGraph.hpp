#pragma once

#include "../include/FlatZonesGraph.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class PyBindFlatZonesGraph: public DefaultFlatZonesGraph {

public:
	
    PyBindFlatZonesGraph(py::array_t<uint8_t> &input, int numRows, int numCols, double radiusAdj);
    PyBindFlatZonesGraph(ImageUInt8Ptr img, double radiusAdj);
	
};
