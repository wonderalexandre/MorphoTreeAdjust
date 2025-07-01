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

    ComponentTreeAdjustmentByLeaf(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree) : ComponentTreeAdjustment(maxtree, mintree) { }

    void buildMergedAndNestedCollections(ComponentTreeFZPtr, std::vector<int>&, int, int, bool) override {
        throw std::runtime_error("Método com vetor de FlatZones não suportado nesta subclasse.");
    }

    void buildMergedAndNestedCollections(ComponentTreeFZPtr tree, int flatZoneID, int pixelUpperBound, int newGrayLevel, bool isMaxtree) override;
    
    void updateTree(ComponentTreeFZPtr tree, NodeFZPtr L_leaf) override;

    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning) override;
    
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning) override;
  
};


#endif