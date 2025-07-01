#include "PyBindComponentTreeAdjustmentByLeaf.hpp"



void PyBindComponentTreeAdjustmentByLeaf::updateTree(PyBindComponentTreeFZPtr tree, NodeFZPtr L_leaf) {
   ComponentTreeAdjustmentByLeaf::updateTree(tree, L_leaf);
}

py::tuple PyBindComponentTreeAdjustmentByLeaf::buildCollections(PyBindComponentTreeFZPtr tree, NodeFZPtr leaf, int newGrayLevel) { 

    int idLeaf = leaf->getCNPsByFlatZone().begin()->second.front(); //pixel (id) of flatzone 
    int pixelUpperBound = idLeaf; 
    NodeFZPtr nodeTauL = tree->getSC(idLeaf); //node of correspondence flatzone in other treee
    FlatZonePtr flatzoneTauL = &tree->getFlatzoneByID(idLeaf); 
    
    ComponentTreeAdjustmentByLeaf::buildMergedAndNestedCollections(tree, flatzoneTauL->front(), pixelUpperBound, newGrayLevel, tree->isMaxtree());

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