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

class FlatZonesGraph;

template <typename CNPsType>
class NodeCT;  // Forward declaration

template <typename CNPsType>
class ComponentTree : public std::enable_shared_from_this<ComponentTree<CNPsType>> {
protected:
    NodeCTPtr<CNPsType> root;
    std::vector<NodeCTPtr<CNPsType>> pixelToNode; //mapping from pixel to node
    
    int numNodes;
    int maxIndex;
    bool maxtreeTreeType; //maxtree is true; mintree is false

    int numCols;
    int numRows;
    int numPixels;
    
    AdjacencyRelation adj; //disk of a given ratio: ratio(1) for 4-connect and ratio(1.5) for 8-connect 
    
    int* countingSort(int* img);
    int* createTreeByUnionFind(int* orderedPixels, int* img);
    int findRoot(int* zPar, int x);
    void reconstruction(NodeCTPtr<CNPsType> node, int* imgOut);


    // Define `flatzoneGraph` apenas para `FlatZones`
    using FlatzoneGraphType = std::conditional_t<std::is_same_v<CNPsType, FlatZones>,
        std::unique_ptr<FlatZonesGraph>, 
        std::monostate>;
    FlatzoneGraphType flatzoneGraph;

public:

    
    ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

    ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation);

    virtual ~ComponentTree() = default;

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
	std::list<int>& getFlatzoneByID(int p);
	
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    ListOfAdjacentFlatzones& getListOfAdjacentFlatzones();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    std::unique_ptr<FlatZonesGraph>& getFlatZonesGraph();

    void assignCNPs(int* img);

    void build(int* img);

    NodeCTPtr<CNPsType> getSC(int p);

    NodeCTPtr<CNPsType> getRoot();

    void setRoot(NodeCTPtr<CNPsType> n);

    bool isMaxtree();

    int getNumNodes();

    void setNumNodes(int numNodes) { this->numNodes = numNodes; }

    int getNumRowsOfImage();

    int getNumColsOfImage();

    void computerArea(NodeCTPtr<CNPsType> node);

    int* reconstructionImage();

    AdjacencyRelation& getAdjacencyRelation();

    void setSC(int p, NodeCTPtr<CNPsType> n);

    void prunning(NodeCTPtr<CNPsType> node);
    
    std::vector<NodeCTPtr<CNPsType>> getLeaves();

    std::vector<NodeCTPtr<CNPsType>> getNodesThreshold(int threshold);
    

};


#include "../include/ComponentTree.tpp"

#endif