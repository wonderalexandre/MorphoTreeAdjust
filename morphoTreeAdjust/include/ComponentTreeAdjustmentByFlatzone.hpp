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

#ifndef COMPONENT_TREE_ADJUSTMENT_FLATZONE_H
#define COMPONENT_TREE_ADJUSTMENT_FLATZONE_H


class ComponentTreeAdjustmentByFlatzone: public  ComponentTreeAdjustment {
private:    
    UnionNodes unionNodeTauSubtree;
    
public:

    ComponentTreeAdjustmentByFlatzone(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree) : ComponentTreeAdjustment(mintree, maxtree), unionNodeTauSubtree(maxtree->isMaxtree(), std::max(mintree->getNumNodes(), maxtree->getNumNodes())) { }
          
    void updateTree(ComponentTreeFZPtr tree, FlatZone* flatzone);
    
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, FlatZone* flatzone);
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, FlatZone* flatzone);
    
    


    
};


#endif