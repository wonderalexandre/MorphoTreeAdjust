
#include "../pybind/PyBindComponentTree.hpp"
#include "../pybind/PybindUtils.hpp"

// Implementação dos construtores
template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(py::array_t<uint8_t> &input, int numRows, int numCols, bool isMaxtree)
    : PyBindComponentTree(
        ImageUInt8::fromExternal(static_cast<uint8_t*>(input.request().ptr), numRows, numCols),
        numRows, numCols, isMaxtree, 1.5)
{}

template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(ImageUInt8Ptr img, int numRows, int numCols, bool isMaxtree, double radiusAdj)
    : ComponentTree<CNPsType>(
        img,
        isMaxtree,
        std::make_shared<AdjacencyRelation>(numRows, numCols, radiusAdj),
        std::make_unique<FlatZonesGraph>(img, std::make_shared<AdjacencyRelation>(numRows, numCols, radiusAdj))
    )
{}
    


template <typename CNPsType>
py::array_t<uint8_t> PyBindComponentTree<CNPsType>::reconstructionImage(){
    ImageUInt8Ptr imgOut = ImageUInt8::create(this->numRows, this->numCols);
    ComponentTreeFZ::reconstruction(this->root, imgOut->rawData());
    return PybindUtils::toNumpy(imgOut);
}

template <typename CNPsType>
py::array_t<uint8_t> PyBindComponentTree<CNPsType>::reconstructionNode(NodeCTPtr<CNPsType> node) {
    

    ImageUInt8Ptr imgOut = ImageUInt8::create(this->getNumRowsOfImage(), this->getNumColsOfImage(), 0);
    for(int p: node->getPixelsOfCC()){
        (*imgOut)[p] = 255;
    }
    return PybindUtils::toNumpy(imgOut);

}

template <typename CNPsType>
std::map<int, NodeCTPtr<CNPsType>> PyBindComponentTree<CNPsType>::getNodes() {
    std::map<int, NodeCTPtr<CNPsType>> nodes_map;
    for (NodeCTPtr<CNPsType> node : this->getRoot()->getIteratorBreadthFirstTraversal()) {
        nodes_map[node->getIndex()] = node;
    }
    return nodes_map;
}


template <typename CNPsType>
py::array_t<uint8_t> PyBindComponentTree<CNPsType>::recNode(NodeCTPtr<CNPsType> _node) {
    int n = _node->getArea();
    auto parent = _node->getParent();
    while (parent != nullptr) {
        n = parent->getArea();
        parent = parent->getParent();
    }

    ImageUInt8Ptr imgOut = ImageUInt8::create(n, 1, 0);
    for(int p: _node->getPixelsOfCC()){
        (*imgOut)[p] = 255;
    }
    return PybindUtils::toNumpy(imgOut);
}