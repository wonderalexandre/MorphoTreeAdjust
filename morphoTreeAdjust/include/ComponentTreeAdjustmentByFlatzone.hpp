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

#ifndef COMPONENT_TREE_ADJUSTMENT_FLATZONE_H
#define COMPONENT_TREE_ADJUSTMENT_FLATZONE_H


class ComponentTreeAdjustmentByFlatzone: public  ComponentTreeAdjustment {

    
public:

    ComponentTreeAdjustmentByFlatzone(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree) : ComponentTreeAdjustment(maxtree, mintree) { }
      
    
    void updateTree(ComponentTreeFZPtr tree, FlatZone* flatzone);
    
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, FlatZone* flatzone);
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, FlatZone* flatzone);
    
    


    
};


#endif