#include <list>
#include <vector>
#include <array>

#include "../include/NodeCT.hpp"
#include "../include/AdjacencyRelation.hpp"

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H

class ComponentTree {

protected:
		
	NodeCT* root;
	int numNodes;
	NodeCT** nodes;
	
	int numCols;
	int numRows;
	bool maxtreeTreeType;
	AdjacencyRelation* adj;
	
	int* countingSort(int* img);
	int* createTreeByUnionFind(int* orderedPixels, int* img);
	int findRoot(int *zPar, int x);
	void reconstruction(NodeCT* node, int* imgOut);

public:
   
	ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

	~ComponentTree();

	void build(int* img);

	NodeCT* getSC(int p);
	
	NodeCT* getRoot();

	void setRoot(NodeCT* n);

	bool isMaxtree();

	int getNumNodes();

	int getNumRowsOfImage();

	int getNumColsOfImage();

	void computerArea(NodeCT* node);

	int* reconstructionImage();

	AdjacencyRelation* getAdjacencyRelation();

	void setSC(int p, NodeCT* n);

	bool prunning(NodeCT* node);

	std::list<NodeCT*> getLeaves();

	std::list<NodeCT*> getDescendantsInPostOrder(NodeCT* rootSubtree);

    std::list<NodeCT*> getNodesThreshold(int threshold);
};

#endif