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
#include "../include/ComponentTreeAdjustmentBySubtree.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_ANY_NODE_H
#define COMPONENT_TREE_ADJUSTMENT_ANY_NODE_H


class ComponentTreeAdjustmentByAnyNode: public  ComponentTreeAdjustment {

protected:
    UnionNodes unionNodeTauSubtree;
    
    
public:

    ComponentTreeAdjustmentByAnyNode(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree) : ComponentTreeAdjustment(maxtree, mintree), unionNodeTauSubtree(maxtree->isMaxtree(), std::max(mintree->getNumNodes(), maxtree->getNumNodes())) { }
      
    void updateTree(ComponentTreeFZPtr tree, NodeFZPtr node);
    
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning);
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning);
    
    


    
};


#endif