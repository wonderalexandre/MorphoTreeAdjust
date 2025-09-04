
#include "../pybind/PyBindComponentTree.hpp"
#include "../pybind/PybindUtils.hpp"

// Implementação dos construtores
template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(py::array_t<uint8_t> &input, int numRows, int numCols, bool isMaxtree)
    : PyBindComponentTree(ImageUInt8::fromExternal(static_cast<uint8_t*>(input.request().ptr), numRows, numCols), isMaxtree, 1.5){}

template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radiusAdj)
    : PyBindComponentTree(std::make_shared<PyBindFlatZonesGraph>(img, radiusAdj), isMaxtree) {}

template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(std::shared_ptr<PyBindFlatZonesGraph> graph, bool isMaxtree)
    : ComponentTree<CNPsType>(graph, isMaxtree ) {}
    

template <typename CNPsType>
py::array_t<uint8_t> PyBindComponentTree<CNPsType>::reconstructionImage(){
    // usa implementação da árvore e converte para numpy
    ImageUInt8Ptr imgOut = ComponentTree<CNPsType>::reconstructionImage();
    return PybindUtils::toNumpy(imgOut);
}

template <typename CNPsType>
py::array_t<uint8_t> PyBindComponentTree<CNPsType>::reconstructionNode(NodeCT<CNPsType> node) {
    ImageUInt8Ptr imgOut = ImageUInt8::create(this->getNumRowsOfImage(), this->getNumColsOfImage(), 0);
    for(int p: node.getPixelsOfCC()){
        (*imgOut)[p] = 255;
    }
    return PybindUtils::toNumpy(imgOut);
}

template <typename CNPsType>
std::map<int, NodeCT<CNPsType>> PyBindComponentTree<CNPsType>::getNodes() {
    std::map<int, NodeCT<CNPsType>> nodes_map;
    for (auto node : this->getRoot().getIteratorBreadthFirstTraversal()) {
        nodes_map[node.getIndex()] = node;
    }
    return nodes_map;
}


// recNode removido (use reconstructionNode)
