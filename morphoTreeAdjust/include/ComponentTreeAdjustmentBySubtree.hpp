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
    
    NodeFZPtr nodeTauStar; //nodeTauStar é o nó correspondente da folha da sub-arvore a ser podada com maior (ou menor, para min-tree) intensidade
    int fzTauStarID;
    
    std::vector<int> flatZonesID; // Lista de flatzones unificadas
    
    std::vector<NodeFZPtr> nodesWithFZToBeRemoved; // nodes que possuem flatzones a serem removidas
    std::vector<bool> isNodesToBeRemoved; // vetor de booleanos que indica se o nó deve ser removido ou não
    
    FlatZone* unifiedFlatzone;
public:

    // Construtor permite definir a ordem de iteração
    UnionNodes(bool isMaxtree, int numNodes) : isMaxtree(isMaxtree), isNodesToBeRemoved(numNodes, false) { }
    

    std::vector<int>& getFlatzonesID() {
        return flatZonesID;
    }

    NodeFZPtr getNodeTauStar() {
        return nodeTauStar;
    }
   
    int getFlatzoneIDTauStar() {
        return fzTauStarID;
    }

    //esse método é bem lento, usar somente para debug
    std::vector<NodeFZPtr> getNodesToBeRemoved() {
        std::vector<NodeFZPtr> nodesToBeRemovedTmp;
        for(NodeFZPtr node: nodesWithFZToBeRemoved) {
            if (isNodesToBeRemoved[node->getIndex()]) {
                nodesToBeRemovedTmp.push_back(node);
            }
        }
        return nodesToBeRemovedTmp;
    }

    void addCNPsToConnectedFlatzone(NodeFZPtr nodeUnion, ComponentTreeFZPtr tree) {
        if (flatZonesID.size() > 1) {
            std::shared_ptr<FlatZonesGraph>& flatzoneGraph = tree->getFlatZonesGraph();
            flatzoneGraph->updateGraph(flatZonesID, unifiedFlatzone->front(), nodeTauStar, tree);
        }

        nodeUnion->addCNPsToConnectedFlatzone(std::move(*unifiedFlatzone), tree);
        this->removeFlatzones();
    }

    void removeFlatzones() {
        for(size_t i = 0; i < flatZonesID.size(); ++i) {
            NodeFZPtr node = nodesWithFZToBeRemoved[i];
            int idFlatZone = flatZonesID[i];
            node->removeFlatzone(idFlatZone);  
            if(node->getNumCNPs() == 0){
                isNodesToBeRemoved[node->getIndex()] = true; // Marca o nó para remoção
            }
        }
    }

    bool isRemoved(NodeFZPtr node){
        return isNodesToBeRemoved[node->getIndex()];
    }
    
    void resetCollections(bool isMaxtree) {
        this->isMaxtree = isMaxtree;
        for (auto& node : nodesWithFZToBeRemoved) {
            this->isNodesToBeRemoved[node->getIndex()] = false; 
        }
        this->fzTauStarID = -1; 
        this->nodesWithFZToBeRemoved.clear();
        this->flatZonesID.clear();    

        this->nodeTauStar = nullptr;
        this->unifiedFlatzone = nullptr;
    }
    
    void addNode(NodeFZPtr nodeTau, std::list<int>& fzTau) {
        int idFzTau = fzTau.front();
        flatZonesID.push_back(idFzTau);
        nodesWithFZToBeRemoved.push_back(nodeTau);
        
        if (!this->nodeTauStar || ( (!isMaxtree && nodeTau->getLevel() > nodeTauStar->getLevel()) || (isMaxtree && nodeTau->getLevel() < nodeTauStar->getLevel()))) {
            this->nodeTauStar = nodeTau;
            this->fzTauStarID = idFzTau; 
        }

        // unifiedFlatzone é a flatzone que será unificada, ou seja, ela terá os pixels de todas as flatzones que estão sendo unificadas formando assim um grande componente conectado
        if(!this->unifiedFlatzone) {
            unifiedFlatzone = &fzTau;
        }else{
            if(unifiedFlatzone->front() < idFzTau) {
                unifiedFlatzone->splice(unifiedFlatzone->end(), fzTau); 
            }else{
                unifiedFlatzone->splice(unifiedFlatzone->begin(), fzTau); 
            }
        }
    }

};


class ComponentTreeAdjustmentBySubtree: public ComponentTreeAdjustment {

protected:
    UnionNodes unionNodeTauSubtree;

public:

    ComponentTreeAdjustmentBySubtree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree) : ComponentTreeAdjustment(maxtree, mintree), unionNodeTauSubtree(maxtree->isMaxtree(), std::max(mintree->getNumNodes(), maxtree->getNumNodes())) { }
      
    void updateTree(ComponentTreeFZPtr tree, NodeFZPtr node);
    
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning) ;
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning) ;


};


#endif