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
    ComponentTreeFZPtr maxtree;     
    ComponentTreeFZPtr mintree;
    
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
    
    ComponentTreeFZPtr getOtherTree(bool isMaxtree){
        return isMaxtree ? mintree : maxtree;
    }

    
public:

    std::string getOutputLog() {
        return outputLog.str();
    }

    ComponentTreeAdjustment(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree): maxtree(maxtree), mintree(mintree), maxIndex(std::max(maxtree->getNumNodes(), mintree->getNumNodes())), F(maxIndex)  { }

    ComponentTreeAdjustment() = delete; 

    virtual ~ComponentTreeAdjustment() = default;
 
    //virtual void updateTree(ComponentTreeFZPtr tree, NodeFZPtr node) = 0;

    //virtual void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning) = 0;
    
    //virtual void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning) = 0;

    
    void buildMergedAndNestedCollections(ComponentTreeFZPtr tree, std::vector<int>& flatZonesID,  int pixelUpperBound, int newGrayLevel, bool isMaxtree){
        Fb.clear();
        F.resetCollection(isMaxtree);
        F.computerAdjacentNodes(tree, flatZonesID);
        NodeFZPtr nodeTauStar = tree->getSC(pixelUpperBound); //pixel de tauStar ou tauL, para termos o node (limite) mais proximo de root

        for (NodeFZPtr nodeNL: F.getAdjacentNodes()) {
            if( (isMaxtree && nodeNL->getLevel() <= newGrayLevel) || (!isMaxtree &&  nodeNL->getLevel() >= newGrayLevel)) { //o nodeNL está entre g(p) e f(p)
                F.addNodesOfPath(nodeNL, nodeTauStar); 
            } 
            else { 
                //o nodeNL está abaixo de g(p). É armazenado somente a raiz da subtree antes de atingir o nivel g(p)
                NodeFZPtr nodeSubtree = nodeNL;
                for (NodeFZPtr n : nodeNL->getNodesOfPathToRoot()) {
                    if ( (isMaxtree && newGrayLevel > n->getLevel()) || (!isMaxtree && newGrayLevel < n->getLevel())) {
                        break;
                    }
                    nodeSubtree = n; 
                }
                // se a subtree tiver level = g(p), então ela entra em F[\lambda]
                if (nodeSubtree->getLevel() == newGrayLevel) {
                    F.addNodesOfPath(nodeSubtree, nodeTauStar); //F_lambda
                } 
                else {
                    Fb.insert(nodeSubtree); //F_{lambda} > b
                }
            }
            
        }
    }
            
    void buildMergedAndNestedCollections(ComponentTreeFZPtr tree, int flatZoneID, int pixelUpperBound, int newGrayLevel, bool isMaxtree){
        Fb.clear();
        F.resetCollection(isMaxtree);
        F.computerAdjacentNodes(tree, flatZoneID);
        NodeFZPtr nodeTauL = tree->getSC(pixelUpperBound); //pixel de tauStar ou tauL, para termos o node (limite) mais proximo de root

        for (NodeFZPtr nodeNL: F.getAdjacentNodes()) {
            if( (isMaxtree && nodeNL->getLevel() <= newGrayLevel) || (!isMaxtree &&  nodeNL->getLevel() >= newGrayLevel)) { //o nodeNL está entre g(p) e f(p)
                F.addNodesOfPath(nodeNL, nodeTauL); 
            } 
            else { 
                //o nodeNL está abaixo de g(p). É armazenado somente a raiz da subtree antes de atingir o nivel g(p)
                NodeFZPtr nodeSubtree = nodeNL;
                for (NodeFZPtr n : nodeNL->getNodesOfPathToRoot()) {
                    if ( (isMaxtree && newGrayLevel > n->getLevel()) || (!isMaxtree && newGrayLevel < n->getLevel())) {
                        break;
                    }
                    nodeSubtree = n; 
                }
                // se a subtree tiver level = g(p), então ela entra em F[\lambda]
                if (nodeSubtree->getLevel() == newGrayLevel) {
                    F.addNodesOfPath(nodeSubtree, nodeTauL); //F_lambda
                } 
                else {
                    if(nodeSubtree->getParent() && nodeSubtree->getParent()->getIndex() != nodeTauL->getIndex() ){ // caso raro: se o pai de nodeSubtree for diferente de tauL
                        F.addNodesOfPath(nodeSubtree->getParent(), nodeTauL); // Adiciona os nodes do caminho até tauL
                    }else{
                        Fb.insert(nodeSubtree); //F_{lambda} > b
                    }
                }
            }
            
        }
    }
  
};


#endif