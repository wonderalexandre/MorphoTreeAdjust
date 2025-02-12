#include <iterator>
#include <functional>

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

        // Verifica se o nó já existe no vetor antes de adicioná-lo
        for (NodeCT* existingNode : nodeList) {
            if (existingNode == node) {
                return;  
            }
        }
        auto& cnpsIsEquals = cnpsIsEqualsMap[node->getLevel()];
        nodeList.push_back(node);  
        cnpsIsEquals.push_back(nodeSubtree->getCNPs().size() == node->getCNPs().size());
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