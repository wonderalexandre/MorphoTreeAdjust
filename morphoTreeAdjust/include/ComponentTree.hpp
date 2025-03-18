#include <list>
#include <vector>
#include <array>
#include <unordered_set>
#include <utility>
#include <optional>
#include <functional>
#include <variant>

#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"

#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H

template <typename CNPsType>
class NodeCT;  // Forward declaration

template <typename CNPsType>
class ComponentTree {
protected:
    NodeCT<CNPsType>* root;
    NodeCT<CNPsType>** pixelToNode; //mapping from pixel to node
    
    int numNodes;
    int maxIndex;
    bool maxtreeTreeType; //maxtree is true; mintree is false

    int numCols;
    int numRows;
    int numPixels;
    
    AdjacencyRelation* adj; //disk of a given ratio: ratio(1) for 4-connect and ratio(1.5) for 8-connect 
    
    int* countingSort(int* img);
    int* createTreeByUnionFind(int* orderedPixels, int* img);
    int findRoot(int* zPar, int x);
    void reconstruction(NodeCT<CNPsType>* node, int* imgOut);

public:

    // Define `flatzoneGraph` apenas para `FlatZones`
    using FlatzoneGraphType = std::conditional_t<std::is_same_v<CNPsType, FlatZones>,
        FlatzoneGraph, 
        std::monostate>;
    FlatzoneGraphType flatzoneGraph;
    
    ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

    ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

    ~ComponentTree();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
	std::list<int>& getFlatzoneByID(int p);
	
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void updateGraphAfterPruning(std::list<FlatZone>& flatZoneList, FlatZone& unifiedFlatzone, NodeFZ* nodeStar);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void updateGraph(std::list<FlatZoneNode>& flatZoneNodeList, FlatZone& unifiedFlatzone, NodeFZ* nodeStar);

    void assignCNPs();

    void build(int* img);

    NodeCT<CNPsType>* getSC(int p);

    NodeCT<CNPsType>* getRoot();

    void setRoot(NodeCT<CNPsType>* n);

    bool isMaxtree();

    int getNumNodes();

    void setNumNodes(int numNodes) { this->numNodes = numNodes; }

    int getNumRowsOfImage();

    int getNumColsOfImage();

    void computerArea(NodeCT<CNPsType>* node);

    int* reconstructionImage();

    AdjacencyRelation* getAdjacencyRelation();

    void setSC(int p, NodeCT<CNPsType>* n);

    void prunning(NodeCT<CNPsType>* node);
    
    std::vector<NodeCT<CNPsType>*> getLeaves();

    std::vector<NodeCT<CNPsType>*> getNodesThreshold(int threshold);
    

};


#include "../include/ComponentTree.tpp"

#endif