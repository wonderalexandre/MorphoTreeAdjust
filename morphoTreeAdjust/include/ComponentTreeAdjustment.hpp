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

public:
    // Construtor permite definir a ordem de iteração
    MergedNodesCollection(int maxIndex): maxIndex(maxIndex) {
        this->visited = new bool[this->maxIndex]();  // Inicializa com false
    }

    ~MergedNodesCollection() {
        delete[] this->visited;  // Libera memória
    }

    std::vector<NodeFZ*>& getMergedNodes(int level) {
        return collectionF[level]; 
    }

    std::array<std::vector<NodeFZ*>, 256> getCollectionF(){
        return collectionF;
    }


    void resetCollection(bool descendingOrder) {
        this->isMaxtree = descendingOrder;
        for (auto& vec : collectionF) {
            vec.clear();
        }
        lambdaList.clear();
        std::fill(visited, visited + maxIndex, false);
        currentIndex = 0;
    }

    void addNodesOfPath(NodeFZ* nodeNL, NodeFZ* nodeTauL) {
        for (NodeFZ* n : nodeNL->getNodesOfPathToRoot()) {
            int index = n->getIndex();
            if (!visited[index]) {
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
    //std::array<std::vector<NodeFZ*>, 256> unionNodes;
    //std::array<std::vector<bool>, 256> cnpsIsEqualsMap;
    bool isMaxtree; // Se true, percorre de forma decrescente



public:
std::vector<FlatZoneRef> flatZonesList;
std::vector<NodeFZ*> nodesList;


    // Construtor permite definir a ordem de iteração
    UnionNodes(bool isMaxtree):  isMaxtree(isMaxtree) {}
    
    
    std::vector<FlatZoneRef>& getFlatzones() {
        return flatZonesList;
    }

    std::vector<NodeFZ*> getNodes() {
        return nodesList;
    }

    void addCNPsToConnectedFlatzone(NodeFZ* nodeUnion, ComponentTreeFZ* tree){
        if (flatZonesList.size() > 1) {
            std::list<int> unifiedFZ;
            std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual> flatzonesToMergeSet;

            for (std::list<int>& fz : flatZonesList) {
                flatzonesToMergeSet.insert(fz);
            }
            
            // Remover do grafo as flatzones que serão fundidas
            std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual> neighborsToTransfer;
            for (std::list<int>& flatzoneRef : flatZonesList) {
                
                // Coletar vizinhos que não fazem parte da fusão
                for (const auto& neighborRef : tree->flatzoneGraph[flatzoneRef]) {
                    if (flatzonesToMergeSet.find(neighborRef) == flatzonesToMergeSet.end()) {
                        tree->flatzoneGraph[neighborRef].erase(flatzoneRef);
                        neighborsToTransfer.insert(neighborRef);
                    }
                }

                // Remover do grafo
                tree->flatzoneGraph.erase(flatzoneRef);

                // Adicionar pixels ao cnpsCC
                unifiedFZ.splice(unifiedFZ.end(), flatzoneRef);
            }

            assert(!unifiedFZ.empty() && "ERRO: unifiedFZ está vazio após a fusão!");

            // Adicionar `cnpsCC` ao grafo com suas conexões
            FlatZoneRef unifiedFZRef = unifiedFZ;
            tree->flatzoneGraph[unifiedFZRef] = neighborsToTransfer;

            for (const auto& neighborRef : neighborsToTransfer) {
                tree->flatzoneGraph[neighborRef].insert(unifiedFZRef);
            }
            assert(tree->flatzoneGraph.find(unifiedFZRef) != tree->flatzoneGraph.end() && "Erro: unifiedFZRef não está registrada no grafo!");
            nodeUnion->addCNPsToConnectedFlatzone(std::move(unifiedFZ), tree);
        }else{
            nodeUnion->addCNPsToConnectedFlatzone(std::move(flatZonesList[0].get()), tree);
        }

    }

    void removeFlatzones() {
        for(int i=0; i < flatZonesList.size(); i++){
            std::list<int>& flatzone = flatZonesList[i].get();
            NodeFZ* node = nodesList[i];
            node->removeFlatzone(flatzone);  
        }
    }
    
    /*
    std::vector<NodeFZ*>& getNodes(int level) {
        return unionNodes[level]; 
    }

    std::vector<bool>& getCnpsIsEquals(int level) {
        return cnpsIsEqualsMap[level]; 
    }
*/
    void resetCollection(bool isMaxtree) {
        this->isMaxtree = isMaxtree;
        flatZonesList.clear();
        nodesList.clear();
    }
    
    void addNode(NodeFZ* node, std::list<int>& fzTau) {
        
        flatZonesList.push_back(fzTau);
        nodesList.push_back(node);  
        //cnpsIsEquals.push_back(removeNode);
        //cnpsIsEquals.push_back(nodeSubtree->getNumFlatzone() == node->getNumFlatzone());
    }
    /*
    class Iterator {
        private:
            UnionNodes& container;
            int level;
            size_t index;
    
            void advanceToNextValid() {
                while (level >= 0 && level < 256) {
                    if (index < container.unionNodes[level].size()) {
                        return;
                    }
                    index = 0;
                    level += (container.isMaxtree ? 1 : -1);
                }
            }
        
        public:
            Iterator(UnionNodes& container, int startLevel) : container(container), level(startLevel), index(0) {
                advanceToNextValid();
            }
    
            std::pair<NodeFZ*, bool> operator*() const {
                return {container.unionNodes[level][index], container.cnpsIsEqualsMap[level][index]};
            }
    
            Iterator& operator++() {
                ++index;
                advanceToNextValid();
                return *this;
            }
    
            bool operator!=(const Iterator& other) const {
                return level != other.level || index != other.index;
            }
        };
    
        Iterator begin() {
            return Iterator(*this, isMaxtree ? 0 : 255);
        }
    
        Iterator end() {
            return Iterator(*this, isMaxtree ? 256 : -1);
        }
    
        class IterableWrapper {
        private:
            UnionNodes& container;
        public:
            IterableWrapper(UnionNodes& container) : container(container) {}
            Iterator begin() { return container.begin(); }
            Iterator end() { return container.end(); }
        };
    
        IterableWrapper getIterator() {
            return IterableWrapper(*this);
        }
     */

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

    void disconnect(NodeFZ* node) {
        if(node->getParent() != nullptr){
	        node->getParent()->getChildren().remove(node);
		    node->setParent(nullptr);
        }
    }
    
    
public:

    ComponentTreeAdjustment(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree); 

    ~ComponentTreeAdjustment(); 
 
    void buildMergedAndNestedCollections(ComponentTreeFZ* tree, std::vector<FlatZoneRef>& flatZone, int newGrayLevel, bool isMaxtree);
    
    std::vector<NodeFZ*> getAdjacentNodes(ComponentTreeFZ* tree, std::vector<FlatZoneRef>& flatZone);

    void updateTree(ComponentTreeFZ* tree, NodeFZ *L_leaf);

    void updateTree2(ComponentTreeFZ* tree, NodeFZ *rSubtree);

    void adjustMinTree(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning);
    
    void adjustMaxTree(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning);

    void adjustMinTree2(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning);
    
    void adjustMaxTree2(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning);

};


#endif