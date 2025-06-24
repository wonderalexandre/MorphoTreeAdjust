#include <iterator>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"
#include "../include/ComponentTreeAdjustment.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_SUBTREE_H
#define COMPONENT_TREE_ADJUSTMENT_SUBTREE_H

class UnionNodes {
protected:
    bool isMaxtree; // Se true, percorre de forma decrescente
    std::list<FlatZoneNode> flatZoneNodeList;
    std::list<int> listFlatzoneIDs;
    FlatZoneNode nodeTauStar; //nodeTauStar é o nó correspondente da folha da sub-arvore a ser podada com maior (ou menor, para min-tree) intensidade
    std::unordered_set<NodeFZPtr> nodesToBeRemoved; //nodes que foram removido devido ao fusão das zonas planas quando newGrayLevel = \lambda

public:

    // Construtor permite definir a ordem de iteração
    UnionNodes(bool isMaxtree):  isMaxtree(isMaxtree) {}
    
    
    std::list<FlatZoneNode>& getFlatzoneNodeList() {
        return flatZoneNodeList;
    }

    std::vector<FlatZonePtr> getFlatzones() {
        std::vector<FlatZonePtr> flatzones;
        for(FlatZoneNode& fzNode: flatZoneNodeList){
            flatzones.push_back(fzNode.flatzone);
        }
        return flatzones;
    }

    FlatZoneNode& getNodeTauStar() {
        return nodeTauStar;
    }
   
    std::unordered_set<NodeFZPtr>& getNodesToBeRemoved() {
        return nodesToBeRemoved;
    }

    void addCNPsToConnectedFlatzone(NodeFZPtr nodeUnion, ComponentTreeFZPtr tree) {
        if (flatZoneNodeList.size() > 1) {
            FlatZone unifiedFlatzone;
            std::unique_ptr<FlatZonesGraph>& flatzoneGraph = tree->getFlatZonesGraph();
            flatzoneGraph->updateGraph(flatZoneNodeList, unifiedFlatzone, nodeTauStar.node, tree);
            nodeUnion->addCNPsToConnectedFlatzone(std::move(unifiedFlatzone), tree);
        }else{
            nodeUnion->addCNPsToConnectedFlatzone(std::move(*nodeTauStar.flatzone), tree);
        }
        
    }

    void removeFlatzones() {
        for(FlatZoneNode& fzNode: flatZoneNodeList){
            NodeFZPtr node = fzNode.node;
            node->removeFlatzone(fzNode.idFlatZone);  
            if(fzNode.node->getNumCNPs() == 0){
                nodesToBeRemoved.insert(fzNode.node);
            }
        }
    }

    bool isRemoved(NodeFZPtr node){
        return nodesToBeRemoved.find(node) == nodesToBeRemoved.end();
    }
    
    void resetCollection(bool isMaxtree) {
        this->isMaxtree = isMaxtree;
        this->flatZoneNodeList.clear();
        this->listFlatzoneIDs.clear();
        this->nodesToBeRemoved.clear();

        this->nodeTauStar.node = nullptr;
        this->nodeTauStar.flatzone = nullptr;
        this->nodeTauStar.idFlatZone = -1;
    }
    
    void addNode(NodeFZPtr nodeTau, std::list<int>& fzTau) {
        flatZoneNodeList.emplace_back(nodeTau, fzTau);
        listFlatzoneIDs.push_back(fzTau.front());
        
        if (!this->nodeTauStar.node || ( (!isMaxtree && nodeTau->getLevel() > nodeTauStar.node->getLevel()) || (isMaxtree && nodeTau->getLevel() < nodeTauStar.node->getLevel()))) {
            this->nodeTauStar.node = nodeTau;
            this->nodeTauStar.flatzone = &fzTau;
        }
    }

};


class ComponentTreeAdjustmentBySubtree: public ComponentTreeAdjustment {

protected:
    UnionNodes unionNodeTauSubtree;

    
public:

    ComponentTreeAdjustmentBySubtree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree) : ComponentTreeAdjustment(maxtree, mintree), unionNodeTauSubtree(maxtree->isMaxtree()) { }
    
    void buildMergedAndNestedCollections(ComponentTreeFZPtr tree, std::vector<FlatZonePtr>& flatZone, int pixelUpperBound, int newGrayLevel, bool isMaxtree);
      
    void updateTree(ComponentTreeFZPtr tree, NodeFZPtr node);
    
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning);
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning);


    
};


#endif