#include <list>
#include <vector>
#include <map>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <pybind11/numpy.h>

#ifndef PYBIND_COMPONENT_TREE_H
#define PYBIND_COMPONENT_TREE_H


namespace py = pybind11;


class PyBindComponentTree: public ComponentTree {

public:
   
    PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree);

	PyBindComponentTree(py::array_t<int> &input, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	py::array_t<int> reconstructionImage();

	py::array_t<int> reconstructionNode(NodeCT* node);


	static py::array_t<int> recNode(NodeCT* _node){
		int n = _node->getArea();
		NodeCT* parent = _node->getParent();
		while(parent != nullptr){
			n = parent->getArea();
			parent = parent->getParent();
		}
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
	std::map<int, NodeCT*> getNodes() ;


};

#endif