#include <iterator>
#include <functional>
#include <iostream>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_H
#define COMPONENT_TREE_ADJUSTMENT_H

#include <array>
#include <vector>

class MergedNodesCollection {
protected:
    std::array<std::vector<NodeCT*>, 256> collectionF;
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

    std::vector<NodeCT*>& getMergedNodes(int level) {
        return collectionF[level]; 
    }

    std::array<std::vector<NodeCT*>, 256> getCollectionF(){
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

    void addNodesOfPath(NodeCT* nodeNL, NodeCT* nodeTauL) {
        for (NodeCT* n : nodeNL->getNodesOfPathToRoot()) {
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
    std::array<std::vector<NodeCT*>, 256> unionNodes;
    std::array<std::vector<bool>, 256> cnpsIsEqualsMap;
    bool isMaxtree; // Se true, percorre de forma decrescente
    int level;

public:
    // Construtor permite definir a ordem de iteração
    UnionNodes(bool isMaxtree):  isMaxtree(isMaxtree) {}
    
    std::vector<NodeCT*>& getNodes(int level) {
        return unionNodes[level]; 
    }

    std::vector<bool>& getCnpsIsEquals(int level) {
        return cnpsIsEqualsMap[level]; 
    }

    void resetCollection(bool isMaxtree) {
        this->isMaxtree = isMaxtree;
        for (auto& vec : unionNodes) {
            vec.clear();
        }
    }

    void addNode(NodeCT* node, NodeCT* nodeSubtree) {
        auto& nodeList = unionNodes[node->getLevel()];

        for (NodeCT* existingNode : nodeList) {
            if (existingNode->getIndex() == node->getIndex()) {
                return;  
            }
        }
        
        auto& cnpsIsEquals = cnpsIsEqualsMap[node->getLevel()];
        nodeList.push_back(node);  
        cnpsIsEquals.push_back(nodeSubtree->getCNPs().size() == node->getCNPs().size());
    }

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
    
            std::pair<NodeCT*, bool> operator*() const {
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
     

};




class ComponentTreeAdjustment {

protected:
    ComponentTree* mintree;
    ComponentTree* maxtree; 
    int maxIndex; 
    bool* visited = nullptr; 
    
    UnionNodes unionNodes;
    MergedNodesCollection F;
    std::vector<NodeCT*> B_L;

    void disconnect(NodeCT* node) {
        if(node->getParent() != nullptr){
	        node->getParent()->getChildren().remove(node);
		    node->setParent(nullptr);
        }
    }
    
    
public:

    ComponentTreeAdjustment(ComponentTree* maxtree, ComponentTree* mintree); 

    ~ComponentTreeAdjustment(); 
 
    void buildMergedAndNestedCollections(ComponentTree* tree, std::list<int> flatZone, int newGrayLevel, bool isMaxtree);
    
    std::vector<NodeCT*> getAdjacentNodes(ComponentTree* tree, std::list<int> flatZone);

    void updateTree(ComponentTree* tree, NodeCT *L_leaf);

    void updateTree2(ComponentTree* tree, NodeCT *rSubtree);

    void adjustMinTree(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning);
    
    void adjustMaxTree(ComponentTree* maxtree, ComponentTree* mintree, std::vector<NodeCT*> nodesToPruning);

    void adjustMinTree2(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning);
    
    void adjustMaxTree2(ComponentTree* maxtree, ComponentTree* mintree, std::vector<NodeCT*> nodesToPruning);

};

#endif