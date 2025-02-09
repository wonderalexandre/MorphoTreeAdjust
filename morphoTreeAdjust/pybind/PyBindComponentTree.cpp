#include <list>
#include <vector>
#include <stack>
#include <queue>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "./PyBindComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <iostream>


PyBindComponentTree::PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree) 
	: PyBindComponentTree(input, numRows, numCols, isMaxtree, 1.5){ }
 
PyBindComponentTree::PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
	: ComponentTree(numRows, numCols, isMaxtree, radiusOfAdjacencyRelation){
	
	auto buf_input = input.request();
		
	int* img = (int *) buf_input.ptr;
	build(img);
} 

py::array_t<int> PyBindComponentTree::reconstructionImage(){
	int n = this->numRows * this->numCols;
	auto img_numpy = py::array(py::buffer_info(
			nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
			sizeof(int),     /* Size of one item */
			py::format_descriptor<int>::value, /* Buffer format */
			1,          /* How many dimensions? */
			{ ( n ) },  /* Number of elements for each dimension */
			{ sizeof(int) }  /* Strides for each dimension */
	));
	auto buf_img = img_numpy.request();
	int *imgOut = (int *) buf_img.ptr;
	this->reconstruction(this->root, imgOut);
	return img_numpy;

}

bool PyBindComponentTree::isNodesInitialized() {
	for (int i = 0; i < this->numRows * this->numCols; i++) {
		if (this->nodes[i] == nullptr) {
			std::cerr << "Erro: nodes[" << i << "] não foi inicializado!" << std::endl;
			return false;
		}
	}
	return true;
}
 

py::array_t<int> PyBindComponentTree::reconstructionNode(NodeCT* _node){
	int n = this->numRows * this->numCols;
	auto img_numpy = py::array(py::buffer_info(
			nullptr,            /* Pointer to data (nullptr -> ask NumPy to allocate!) */
			sizeof(int),     /* Size of one item */
			py::format_descriptor<int>::value, /* Buffer format */
			1,          /* How many dimensions? */
			{ ( n ) },  /* Number of elements for each dimension */
			{ sizeof(int) }  /* Strides for each dimension */
	));
	auto buf_img = img_numpy.request();
	int *imgOut = (int *) buf_img.ptr;
	
	for(int p=0; p < n; p++)
		imgOut[p] = 0;
	
	std::stack<NodeCT*> s;
	s.push(_node);
	while(!s.empty()){
    	NodeCT* node = s.top(); s.pop();
		for(int p: node->getCNPs())
			imgOut[p] = 255;
		for (NodeCT *child: node->getChildren()){
			s.push(child);
		}
	}
	return img_numpy;
}

std::map<int, NodeCT*> PyBindComponentTree::getNodes() {
	std::map<int, NodeCT*> nodes;
	for (NodeCT* node : this->root->getIteratorBreadthFirstTraversal()) {
		nodes[node->getIndex()] = node;
	}
	return nodes;
}

NodeCT* PyBindComponentTree::getNodeByIndex(int index){
	for (NodeCT* node : this->root->getIteratorBreadthFirstTraversal()) {
		if(node->getIndex() == index){
			return node;
		}
	}
	return nullptr;
}