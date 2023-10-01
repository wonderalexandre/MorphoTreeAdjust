#include <unordered_set>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_H
#define COMPONENT_TREE_ADJUSTMENT_H




class ComponentTreeAdjustment {


public:

    void addNodesOfPath(NodeCT* nodeNL, NodeCT* nodeTauL, std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F);

    void adjustMinTree(ComponentTree &mintree, NodeCT *Lmax);

};

#endif