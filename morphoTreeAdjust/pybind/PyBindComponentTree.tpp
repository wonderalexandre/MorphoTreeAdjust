
#include "../pybind/PyBindComponentTree.hpp"


// Implementação dos construtores
template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree)
    : PyBindComponentTree(input, numRows, numCols, isMaxtree, 1.5) {}

template <typename CNPsType>
PyBindComponentTree<CNPsType>::PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
    : ComponentTree<CNPsType>(numRows, numCols, isMaxtree, radiusOfAdjacencyRelation) {
    auto buf_input = input.request();
    int* img = (int *)buf_input.ptr;
    this->build(img);
}

template <typename CNPsType>
py::array_t<int> PyBindComponentTree<CNPsType>::reconstructionImage() {
    int n = this->getNumRowsOfImage() * this->getNumColsOfImage();
    auto img_numpy = py::array(py::buffer_info(
        nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {n}, {sizeof(int)}
    ));
    auto buf_img = img_numpy.request();
    int *imgOut = (int *)buf_img.ptr;
    this->reconstruction(this->getRoot(), imgOut);
    return img_numpy;
}

template <typename CNPsType>
py::array_t<int> PyBindComponentTree<CNPsType>::reconstructionNode(NodeCT<CNPsType>* node) {
    int n = this->getNumRowsOfImage() * this->getNumColsOfImage();
    auto img_numpy = py::array(py::buffer_info(
        nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {n}, {sizeof(int)}
    ));
    auto buf_img = img_numpy.request();
    int* imgOut = (int*)buf_img.ptr;

    for (int i = 0; i < n; i++) imgOut[i] = 0;

    std::stack<NodeCT<CNPsType>*> s;
    s.push(node);

    while (!s.empty()) {
        NodeCT<CNPsType>* n = s.top();
        s.pop();

        for (int p : n->getCNPs()) {
            imgOut[p] = 255;
        }

        for (NodeCT<CNPsType>* child : n->getChildren()) {
            s.push(child);
        }
    }

    return img_numpy;
}

template <typename CNPsType>
std::map<int, NodeCT<CNPsType>*> PyBindComponentTree<CNPsType>::getNodes() {
    std::map<int, NodeCT<CNPsType>*> nodes_map;
    for (NodeCT<CNPsType>* node : this->getRoot()->getIteratorBreadthFirstTraversal()) {
        nodes_map[node->getIndex()] = node;
    }
    return nodes_map;
}


template <typename CNPsType>
py::array_t<int> PyBindComponentTree<CNPsType>::recNode(NodeCT<CNPsType>* _node) {
    int n = _node->getArea();
    NodeCT<CNPsType>* parent = _node->getParent();
    while (parent != nullptr) {
        n = parent->getArea();
        parent = parent->getParent();
    }

    auto img_numpy = py::array(py::buffer_info(
        nullptr, sizeof(int), py::format_descriptor<int>::value,
        1, {n}, {sizeof(int)}
    ));
    auto buf_img = img_numpy.request();
    int* imgOut = (int*) buf_img.ptr;

    for (int p = 0; p < n; p++)
        imgOut[p] = 0;

    std::stack<NodeCT<CNPsType>*> s;
    s.push(_node);

    while (!s.empty()) {
        NodeCT<CNPsType>* node = s.top();
        s.pop();
        for (int p : node->getCNPs()) {
            imgOut[p] = 255;
        }
        for (NodeCT<CNPsType>* child : node->getChildren()) {
            s.push(child);
        }
    }
    return img_numpy;
}
