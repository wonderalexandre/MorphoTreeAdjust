#include "PyBindComponentTreeAdjustmentBySubtree.hpp"



void PyBindComponentTreeAdjustmentBySubtree::updateTree(PyBindComponentTreeFZPtr tree, NodeFZPtr nodeSubtree) {
   ComponentTreeAdjustmentBySubtree::updateTree(tree, nodeSubtree);
}

py::tuple PyBindComponentTreeAdjustmentBySubtree::buildCollections(PyBindComponentTreeFZPtr tree, NodeFZPtr rootSubtree, int newGrayLevel) { 
    
    ComponentTreeFZPtr otherTree = tree->isMaxtree()? this->mintree : this->maxtree;

    bool isMaxtree = tree->isMaxtree();
    
    unionNodeTauSubtree.resetCollections(isMaxtree);
    for (NodeFZPtr nSubtree : rootSubtree->getIteratorBreadthFirstTraversal()) {
        for(auto& [idFlatZoneNSubtree, fzSubtree]: nSubtree->getCNPsByFlatZone()){    
            NodeFZPtr nodeTau = tree->getSC(idFlatZoneNSubtree);
            FlatZone& fzTau = nodeTau->getFlatZone(idFlatZoneNSubtree); //tree->getFlatzoneByID(idFlatZoneNSubtree);
            unionNodeTauSubtree.addNode(nodeTau, fzTau); //, fzSubtree.size() == fzTau.size()  
        }
    }
    int pixelUpperBound = unionNodeTauSubtree.getFlatzoneIDTauStar();
    
    ComponentTreeAdjustmentBySubtree::buildMergedAndNestedCollections(tree,  unionNodeTauSubtree.getFlatzonesID(), pixelUpperBound, newGrayLevel, isMaxtree);
    
    std::array<std::vector<NodeFZPtr>, 256>& collectionF = this->F.getCollectionF();
    std::map<int, std::vector<NodeFZPtr>> mapCollectionF;
    for (int i = 0; i < 256; ++i) {
        if (!collectionF[i].empty()) {
            std::vector<NodeFZPtr> nodes = collectionF[i];
            mapCollectionF[i] = nodes;
        }
    }

    return py::make_tuple(mapCollectionF, this->Fb);
}