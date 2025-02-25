#include <list>
#include <vector>
#include <array>
#include <unordered_set>
#include <utility>
#include <optional>
#include <functional>

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
    int numNodes;
    int maxIndex;
    NodeCT<CNPsType>** nodes;
    
    //std::vector<FlatZone*> pixelToFlatzone;
    // Define `pixelToFlatzone` apenas para `FlatZones`
    using PixelToFlatzoneType = std::conditional_t<std::is_same_v<CNPsType, FlatZones>, 
        std::vector<FlatZone*>, 
        std::monostate>;
    PixelToFlatzoneType pixelToFlatzone;
    
    int numCols;
    int numRows;
    bool maxtreeTreeType;
    AdjacencyRelation* adj;
    
    int* countingSort(int* img);
    int* createTreeByUnionFind(int* orderedPixels, int* img);
    int findRoot(int* zPar, int x);
    void reconstruction(NodeCT<CNPsType>* node, int* imgOut);

public:

    //std::unordered_map<std::reference_wrapper<std::list<int>>, std::unordered_set<std::reference_wrapper<std::list<int>>, ListRefHash, ListRefEqual>, ListRefHash, ListRefEqual> flatzoneGraph;
  // Define `flatzoneGraph` apenas para `FlatZones`
    using FlatzoneGraphType = 
        std::conditional_t<std::is_same_v<CNPsType, FlatZones>,
        std::unordered_map<FlatZoneRef, std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual>, ListRefHash, ListRefEqual>, 
        std::monostate>;
    FlatzoneGraphType flatzoneGraph;
    
    ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

    ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

    ~ComponentTree();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
	std::list<int>& getFlatzoneRef(int p);
	
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    std::list<int>* getFlatzonePointer(int p);
	
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void updatePixelToFlatzone(int p, std::list<int>* newFlatzone);
	
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void buildFlatzoneGraph();

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


#include "ComponentTree.tpp"


#endif