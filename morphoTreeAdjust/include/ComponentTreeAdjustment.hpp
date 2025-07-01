#include <iterator>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_H
#define COMPONENT_TREE_ADJUSTMENT_H

#include <array>
#include <vector>
#include <unordered_set>

class MergedNodesCollection {
protected:
    std::array<std::vector<NodeFZPtr>, 256> collectionF;
    std::vector<bool> visited;
    std::vector<bool> visitedAdj;

    int maxIndex;
    std::vector<int> lambdaList; // Lista ordenada de lambdas (sempre crescente)
    int currentIndex = 0; // Índice atual
    bool isMaxtree; // Se true, percorre de forma decrescente
    std::vector<NodeFZPtr> nodesNL;

public:
    
    // Construtor permite definir a ordem de iteração
    MergedNodesCollection(int maxIndex): maxIndex(maxIndex) {
        this->visited.resize(this->maxIndex, false);  // Inicializa com false
        this->visitedAdj.resize(this->maxIndex, false);
    }

    std::vector<NodeFZPtr>& getMergedNodes(int level) {
        return collectionF[level]; 
    }

    void computerAdjacentNodes(ComponentTreeFZPtr tree, const std::vector<int>& flatZonesID) {
        bool isMaxtree = tree->isMaxtree();
        FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
        for (int flatZoneID_P : flatZonesID) {   
            int grayFlatZoneP = tree->getSC(flatZoneID_P)->getLevel(); //is same that: f(p)
    
            for (int flatZoneID_Q : graph->getAdjacentFlatzones(flatZoneID_P)) {
                NodeFZPtr node = tree->getSC(flatZoneID_Q);
                if ( (isMaxtree && node->getLevel() > grayFlatZoneP) || (!isMaxtree && node->getLevel() < grayFlatZoneP) ) {
                    if(!visitedAdj[node->getIndex()]){
                        nodesNL.push_back(node);  
                        visitedAdj[node->getIndex()] = true;
                    }
                }
            }
        }
        resetAdjacentNode();
    }

    void computerAdjacentNodes(ComponentTreeFZPtr tree, int flatZoneID_P) {
        bool isMaxtree = tree->isMaxtree();
        FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
        
        int grayFlatZoneP = tree->getSC(flatZoneID_P)->getLevel(); //is same that: f(p)
        for (int flatZoneID_Q : graph->getAdjacentFlatzones(flatZoneID_P)) {
            NodeFZPtr node = tree->getSC(flatZoneID_Q);
            if ( (isMaxtree && node->getLevel() > grayFlatZoneP) || (!isMaxtree && node->getLevel() < grayFlatZoneP) ) {
                if(!visitedAdj[node->getIndex()]){
                    nodesNL.push_back(node);  
                    visitedAdj[node->getIndex()] = true;
                }
            }
        }
        resetAdjacentNode();
    }

    std::vector<NodeFZPtr>& getAdjacentNodes(){
        return nodesNL;
    }


    std::array<std::vector<NodeFZPtr>, 256>& getCollectionF(){
        return collectionF;
    }

    void resetAdjacentNode(){
        for(NodeFZPtr node: nodesNL){
            visitedAdj[node->getIndex()] = false;
        }
    }

    void resetCollection(bool descendingOrder) {
        this->isMaxtree = descendingOrder;
        for (auto& vec : collectionF) {
            vec.clear();
        }
        nodesNL.clear();
        lambdaList.clear();
        std::fill(visited.begin(), visited.end(), false);
        currentIndex = 0;
    }

    void addNodesOfPath(NodeFZPtr nodeNL, NodeFZPtr nodeTauL) {
        if(!visited[nodeNL->getIndex()])
            for (NodeFZPtr n : nodeNL->getNodesOfPathToRoot()) {
                int index = n->getIndex();
                if(!visited[index]) {
                    collectionF[n->getLevel()].push_back(n);
                    visited[index] = true;
                }else{
                    break;
                }
                
                if (n == nodeTauL) {
                    break;
                }
            }
    }

    int firstLambda() {
        lambdaList.clear();
        for (int i = 0; i < 256; ++i) {
            if (!collectionF[i].empty()) {
                lambdaList.push_back(i);
            }
        }
        currentIndex = isMaxtree ? lambdaList.size() - 1 : 0;
        return lambdaList[currentIndex];
    }


    int nextLambda() {
        if (isMaxtree) {
            return lambdaList[--currentIndex];
        } else {
            return lambdaList[++currentIndex];
        }
    }
};


class ComponentTreeAdjustment {

protected:
    ComponentTreeFZPtr mintree;
    ComponentTreeFZPtr maxtree; 
    int maxIndex; 
    MergedNodesCollection F;
    std::unordered_set<NodeFZPtr> Fb;
    std::ostringstream outputLog;

    

    void disconnect(NodeFZPtr node, bool isFreeMemory=false) {
        if(node->getParent() != nullptr){
	        node->getParent()->getChildren().remove(node);
		    node->setParent(nullptr);
            if(isFreeMemory){
               // delete node; 
                node = nullptr;
            }
        }
    }

    void mergedParentAndChildren(NodeFZPtr nodeUnion, NodeFZPtr n){
        for (NodeFZPtr son : n->getChildren()) {
            son->setParent(nodeUnion);
        }
        nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());
    }
    
    
public:

    std::string getOutputLog() {
        return outputLog.str();
    }

    ComponentTreeAdjustment(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree): mintree(mintree), maxtree(maxtree), maxIndex(std::max(maxtree->getNumNodes(), mintree->getNumNodes())), F(maxIndex)  { }

    ComponentTreeAdjustment() = delete; 

    virtual ~ComponentTreeAdjustment() = default;
 
    virtual void buildMergedAndNestedCollections(ComponentTreeFZPtr, std::vector<int>&, int, int, bool) = 0;
    
    virtual void buildMergedAndNestedCollections(ComponentTreeFZPtr, int, int, int, bool) = 0;
    
    virtual void updateTree(ComponentTreeFZPtr tree, NodeFZPtr node) = 0;

    virtual void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning) = 0;
    
    virtual void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning) = 0;
  
};


#endif