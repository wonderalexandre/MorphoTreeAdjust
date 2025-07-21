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

#ifndef COMPONENT_TREE_ADJUSTMENT_LEAF_H
#define COMPONENT_TREE_ADJUSTMENT_LEAF_H


class ComponentTreeAdjustmentByLeaf: public ComponentTreeAdjustment {
    
public:

    ComponentTreeAdjustmentByLeaf(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree) : ComponentTreeAdjustment(mintree, maxtree) { }
    
    void updateTree(ComponentTreeFZPtr tree, NodeFZPtr L_leaf);

    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning) ;
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning) ;
  
};


#endif