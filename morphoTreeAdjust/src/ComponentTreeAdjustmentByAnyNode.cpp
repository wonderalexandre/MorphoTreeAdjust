#include "../include/ComponentTreeAdjustmentByAnyNode.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>


 

void ComponentTreeAdjustmentByAnyNode::updateTree(ComponentTreeFZPtr tree, NodeFZ node) {
    long int areaFZsRemovedTmp = 0;
    for(int rep: node.getRepCNPs()) {
        ComponentTreeAdjustmentByFlatzone::updateTree(tree, rep);
        areaFZsRemovedTmp += areaFZsRemoved;
    }
    areaFZsRemoved = areaFZsRemovedTmp;
}


void ComponentTreeAdjustmentByAnyNode::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToRemoved) {
    for (NodeId nodeId : nodesToRemoved) {  	
        NodeFZ node = maxtree->proxy(nodeId);	
        assert(node != maxtree->getRoot() && "node is root");
        updateTree(mintree, node); 
        mergeWithParent(maxtree, node); 
    }
}

void ComponentTreeAdjustmentByAnyNode::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToRemoved) {
    for (NodeId nodeId : nodesToRemoved) {  
        NodeFZ node = mintree->proxy(nodeId);	
        assert(node != mintree->getRoot() && "node is root");

        updateTree(maxtree, node);   
        mergeWithParent(mintree, node); 
    }
}
