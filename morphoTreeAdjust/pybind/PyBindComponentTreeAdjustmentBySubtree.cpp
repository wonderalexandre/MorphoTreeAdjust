#include "PyBindComponentTreeAdjustmentBySubtree.hpp"



void PyBindComponentTreeAdjustmentBySubtree::updateTree(PyBindComponentTreeFZPtr tree, NodeFZ nodeSubtree) {
   ComponentTreeAdjustmentBySubtree::updateTree(tree, nodeSubtree);
}

py::tuple PyBindComponentTreeAdjustmentBySubtree::buildCollections(PyBindComponentTreeFZPtr tree, NodeFZ rootSubtree, int newGrayLevel) { 
    ComponentTreeFZPtr otherTree = tree->isMaxtree()? this->mintree : this->maxtree;
    bool isMaxtree = tree->isMaxtree();
    
    properPartsCollector.resetCollections(isMaxtree);
    for (auto nSubtree : rootSubtree.getIteratorBreadthFirstTraversal()) {
        for(int repFZ : nSubtree.getRepCNPs()) {   
            properPartsCollector.addNode(tree, tree->getSC(repFZ), repFZ);
        }
    }
    int pixelUpperBound = properPartsCollector.getRepFZTauStar();
    
    ComponentTreeAdjustmentBySubtree::buildMergedAndNestedCollections(tree,  properPartsCollector.getRepsFZ(), pixelUpperBound, newGrayLevel, isMaxtree);
    
    auto &collectionF = this->F.getCollectionF();
    std::map<int, std::vector<NodeFZ>> mapCollectionF;
    for (int i = 0; i < 256; ++i) {
        if (!collectionF[i].empty()) {
            std::vector<NodeFZ> nodes;
            nodes.reserve(collectionF[i].size());
            for (NodeId nid : collectionF[i]) nodes.push_back(tree->proxy(nid));
            mapCollectionF[i] = std::move(nodes);
        }
    }

    return py::make_tuple(mapCollectionF, this->F.getFb());
}
