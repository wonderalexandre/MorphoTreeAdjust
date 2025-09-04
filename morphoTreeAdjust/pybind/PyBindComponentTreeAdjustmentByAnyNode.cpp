#include "PyBindComponentTreeAdjustmentByAnyNode.hpp"



void PyBindComponentTreeAdjustmentByAnyNode::updateTree(PyBindComponentTreeFZPtr tree, NodeFZ node) {
   ComponentTreeAdjustmentByAnyNode::updateTree(tree, node);
}

/*
py::tuple PyBindComponentTreeAdjustmentByAnyNode::buildCollections(PyBindComponentTreeFZPtr tree, NodeFZPtr node, int newGrayLevel) { 

    int idLeaf = leaf->getCNPsByFlatZone().begin()->second.front(); //pixel (id) of flatzone 
    int pixelUpperBound = idLeaf; 
    NodeFZPtr nodeTauL = tree->getSC(idLeaf); //node of correspondence flatzone in other treee
    FlatZone* flatzoneTauL = &tree->getFlatzoneByID(idLeaf); 
    
    ComponentTreeAdjustmentByAnyNode::buildMergedAndNestedCollections(tree, flatzoneTauL->front(), pixelUpperBound, newGrayLevel, tree->isMaxtree());

    std::array<std::vector<NodeFZPtr>, 256>& collectionF = this->F.getCollectionF();
    std::map<int, std::vector<NodeFZPtr>> mapCollectionF;
    for (int i = 0; i < 256; ++i) {
        if (!collectionF[i].empty()) {
            std::vector<NodeFZPtr> nodes = collectionF[i];
            mapCollectionF[i] = nodes;
        }
    }

    return py::make_tuple(mapCollectionF, this->Fb);
}*/