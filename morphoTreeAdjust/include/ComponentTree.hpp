#include <list>
#include <vector>
#include <array>

#include <optional>
#include <functional>

#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H



class NodeCT;  //Forward declaration


class ComponentTree {

protected:
		
	NodeCT* root;
	int numNodes;
	int maxIndex;
	NodeCT** nodes;
	
	std::vector<std::list<int>*> pixelToFlatzone;  

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

	std::list<int>& getFlatzoneRef(int p);
	std::list<int>* getFlatzonePointer(int p);
	
	void updatePixelToFlatzone(int p, std::list<int>* newFlatzone);

	void build(int* img);

	NodeCT* getSC(int p);
	
	NodeCT* getRoot();

	void setRoot(NodeCT* n);

	bool isMaxtree();

	int getNumNodes();

	void setNumNodes(int numNodes){this->numNodes = numNodes;}

	int getNumRowsOfImage();

	int getNumColsOfImage();

	void computerArea(NodeCT* node);

	int* reconstructionImage();

	AdjacencyRelation* getAdjacencyRelation();

	void setSC(int p, NodeCT* n);

	void prunning(NodeCT* node);

	std::vector<NodeCT*> getLeaves();

    std::vector<NodeCT*> getNodesThreshold(int threshold);

	void assignCNPs();


	
};


#endif