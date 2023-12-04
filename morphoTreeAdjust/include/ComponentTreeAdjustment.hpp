#include <unordered_set>
#include <map>
#include<iterator>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_H
#define COMPONENT_TREE_ADJUSTMENT_H

class ComponentTreeAdjustment {

protected:
    

    std::map <int, std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* > collectionF;

    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* getF(int level){
        if(existsInF(level))
            return collectionF[level];
        else{
            std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F = new std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>();
            collectionF.insert({ level, F }); 
            return collectionF[level];
        }
    }
    bool existsInF(int level){
        return collectionF.count(level) == 1;
    }

    void clearCollectionF(){
        std::map <int, std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* >::iterator it;
        for (it = collectionF.begin(); it != collectionF.end();  ++it){
            std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F = it->second;
            delete F;
        }
        collectionF.clear();
    }

    void disconnect(NodeCT* node) {
        if(node->getParent() != nullptr){
	        node->getParent()->getChildren().remove(node);
		    node->setParent(nullptr);
        }
    }
    
    void addNodesOfPath(NodeCT* nodeNL, NodeCT* nodeTauL) {
        for (NodeCT* n : nodeNL->getNodesOfPathToRoot()) {
            std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F = getF(n->getLevel());
            F->insert(n);
            if (n == nodeTauL){ 
                break; 
            }
        }
    }

public:

    
    void adjustMinTree(ComponentTree &mintree, NodeCT *Lmax);

    void adjustMaxTree(ComponentTree &maxtree, NodeCT *Lmin);

    void adjustMinTree(ComponentTree &mintree, ComponentTree &maxtree, std::list<NodeCT*> nodesToPruning);
    
    void adjustMaxTree(ComponentTree &maxtree, ComponentTree &mintree, std::list<NodeCT*> nodesToPruning);

};

#endif