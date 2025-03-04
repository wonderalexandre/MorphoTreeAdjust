#include <iterator>
#include <functional>
#include <iostream>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_H
#define COMPONENT_TREE_ADJUSTMENT_H

#include <array>
#include <vector>

class MergedNodesCollection {
protected:
    std::array<std::vector<NodeFZ*>, 256> collectionF;
    bool* visited;
    
    int maxIndex;
    std::vector<int> lambdaList; // Lista ordenada de lambdas (sempre crescente)
    int currentIndex = 0; // Índice atual
    bool isMaxtree; // Se true, percorre de forma decrescente
    std::vector<NodeFZ*> nodesNL;

public:
    bool* visitedAdj;
    // Construtor permite definir a ordem de iteração
    MergedNodesCollection(int maxIndex): maxIndex(maxIndex) {
        this->visited = new bool[this->maxIndex]();  // Inicializa com false
        this->visitedAdj = new bool[this->maxIndex]();
    }

    ~MergedNodesCollection() {
        delete[] this->visited;  // Libera memória
        delete[] this->visitedAdj;
    }

    std::vector<NodeFZ*>& getMergedNodes(int level) {
        return collectionF[level]; 
    }

    void computerAdjacentNodes(ComponentTreeFZ* tree, std::vector<FlatZoneRef>& flatZones) {
        bool isMaxtree = tree->isMaxtree();
        
        for (FlatZoneRef flatZonePRef : flatZones) {   
            FlatZone& flatZoneP = flatZonePRef.get();
            int flatZoneID_P = tree->getIdFlatZone(flatZoneP);
            int grayFlatZoneP = tree->getSC(flatZoneID_P)->getLevel(); //is same that: f(p)
    
            for (int flatZoneID_Q : *tree->flatzoneGraph[flatZoneID_P]) {
                NodeFZ* node = tree->getSC(flatZoneID_Q);
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

    std::vector<NodeFZ*>& getAdjacentNodes(){
        return nodesNL;
    }


    std::array<std::vector<NodeFZ*>, 256>& getCollectionF(){
        return collectionF;
    }

    void resetAdjacentNode(){
        for(NodeFZ* node: nodesNL){
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
        std::fill(visited, visited + maxIndex, false);
        currentIndex = 0;
    }

    void addNodesOfPath(NodeFZ* nodeNL, NodeFZ* nodeTauL) {
        if(!visited[nodeNL->getIndex()])
            for (NodeFZ* n : nodeNL->getNodesOfPathToRoot()) {
                int index = n->getIndex();
                if(!visited[index]) {
                    collectionF[n->getLevel()].push_back(n);
                    visited[index] = true;
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

class UnionNodes {
protected:
    bool isMaxtree; // Se true, percorre de forma decrescente
    std::list<FlatZoneNode> flatZoneNodeList;
    std::list<int> listFlatzoneIDs;
    NodeFZ* nodeStar;
public:

    // Construtor permite definir a ordem de iteração
    UnionNodes(bool isMaxtree):  isMaxtree(isMaxtree) {}
    
    
    std::list<FlatZoneNode>& getFlatzoneNodeList() {
        return flatZoneNodeList;
    }

    void setNodeStar(NodeFZ* nodeStar){
       this->nodeStar = nodeStar; 
    }

    std::vector<FlatZoneRef> getFlatzones() {
        std::vector<FlatZoneRef> flatzones;
        for(FlatZoneNode fzNode: flatZoneNodeList){
            flatzones.push_back(*fzNode.flatzone);
        }
        return flatzones;
    }


    void addCNPsToConnectedFlatzone(NodeFZ* nodeUnion, ComponentTreeFZ* tree) {
        if (flatZoneNodeList.size() > 1) {
            FlatZone unifiedFlatzone;
            tree->updateGraph(flatZoneNodeList, unifiedFlatzone, nodeStar);
            nodeUnion->addCNPsToConnectedFlatzone(std::move(unifiedFlatzone), tree);
        } else {
            nodeUnion->addCNPsToConnectedFlatzone(std::move(*flatZoneNodeList.front().flatzone), tree);
        }
    }

    void removeFlatzones(ComponentTreeFZ* tree) {
        for(int idFlatZone: listFlatzoneIDs){
            NodeFZ* node = tree->getSC(idFlatZone);
            node->removeFlatzone(idFlatZone);  
        }
    }
    
    void resetCollection(bool isMaxtree) {
        this->isMaxtree = isMaxtree;
        flatZoneNodeList.clear();
        listFlatzoneIDs.clear();
    }
    
    void addNode(NodeFZ* node, std::list<int>& fzTau) {
        flatZoneNodeList.emplace_back(node, fzTau);
        listFlatzoneIDs.push_back(fzTau.front());
    }

};


class ComponentTreeAdjustment {

protected:
    ComponentTreeFZ* mintree;
    ComponentTreeFZ* maxtree; 
    int maxIndex; 
    int pixelUpperBound=-1;    
    UnionNodes unionNodeTauSubtree;
    MergedNodesCollection F;
    std::vector<NodeFZ*> B_L;

    void disconnect(NodeFZ* node, bool isFreeMemory=false) {
        if(node->getParent() != nullptr){
	        node->getParent()->getChildren().remove(node);
		    node->setParent(nullptr);
            if(isFreeMemory){
                delete node; 
                node = nullptr;
            }
        }
    }

    void mergedParentAndChildren(NodeFZ* nodeUnion, NodeFZ* n){
        for (NodeFZ* son : n->getChildren()) {
            son->setParent(nodeUnion);
        }
        nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());
    }
    
    
public:

    ComponentTreeAdjustment(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree); 

    ~ComponentTreeAdjustment(); 
 
    void buildMergedAndNestedCollections(ComponentTreeFZ* tree, std::vector<FlatZoneRef>& flatZone, int newGrayLevel, bool isMaxtree);
    
    void updateTree(ComponentTreeFZ* tree, NodeFZ *L_leaf);

    void updateTree2(ComponentTreeFZ* tree, NodeFZ *rSubtree);


    void adjustMinTree(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning);
    
    void adjustMaxTree(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning);

    void adjustMinTree2(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning);
    
    void adjustMaxTree2(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning);



    
};


#endif