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
    bool descendingOrder; // Se true, percorre de forma decrescente

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
        this->descendingOrder = descendingOrder;
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

    void prepareLambdaList() {
        lambdaList.clear();
        for (int i = 0; i < 256; ++i) {
            if (!collectionF[i].empty()) {
                lambdaList.push_back(i);
            }
        }
        currentIndex = descendingOrder ? lambdaList.size() - 1 : 0;
    }

    int firstLambda() const {
        return lambdaList.empty() ? -1 : lambdaList[currentIndex];
    }

    int nextLambda() {
        if (lambdaList.empty()) return -1;

        if (descendingOrder) {
            return (currentIndex > 0) ? lambdaList[--currentIndex] : -1;
        } else {
            return (currentIndex + 1 < lambdaList.size()) ? lambdaList[++currentIndex] : -1;
        }
    }
};


class ComponentTreeAdjustment {

protected:
    ComponentTree* mintree;
    ComponentTree* maxtree; 
    int maxIndex; 
    bool* visited = nullptr; 
    
    
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
    
    std::vector<NodeCT*> getAdjacentNodes(ComponentTree* tree, std::list<int> flatZone, int grayFlatZone);

    void updateTree(ComponentTree* tree, NodeCT *L_leaf);

    void updateTree2(ComponentTree* tree, NodeCT *L_leaf);

    void adjustMinTree(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning);
    
    void adjustMaxTree(ComponentTree* maxtree, ComponentTree* mintree, std::vector<NodeCT*> nodesToPruning);

};

#endif