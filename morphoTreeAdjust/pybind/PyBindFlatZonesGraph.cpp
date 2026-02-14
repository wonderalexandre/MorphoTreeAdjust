
#include "../pybind/PyBindFlatZonesGraph.hpp"

// Implementação dos construtores
PyBindFlatZonesGraph::PyBindFlatZonesGraph(py::array_t<uint8_t> &input, int numRows, int numCols, double radiusAdj) : 
    PyBindFlatZonesGraph(ImageUInt8::fromExternal(static_cast<uint8_t*>(input.request().ptr), numRows, numCols), radiusAdj) {}

PyBindFlatZonesGraph::PyBindFlatZonesGraph(ImageUInt8Ptr img, double radiusAdj) : DefaultFlatZonesGraph(img, radiusAdj){}
    
