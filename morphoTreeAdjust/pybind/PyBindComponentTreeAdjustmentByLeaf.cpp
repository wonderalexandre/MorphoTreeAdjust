#include "PyBindComponentTreeAdjustmentByLeaf.hpp"



void PyBindComponentTreeAdjustmentByLeaf::updateTree(PyBindComponentTreeFZPtr tree, NodeFZPtr L_leaf) {
   ComponentTreeAdjustmentByLeaf::updateTree(tree, L_leaf);
}

py::tuple PyBindComponentTreeAdjustmentByLeaf::buildCollections(PyBindComponentTreeFZPtr tree, std::vector<int> vflatZone, int newGrayLevel, bool isMaxtree) { 
    std::shared_ptr<FlatZone> flatzone = std::make_shared<FlatZone>(vflatZone.begin(), vflatZone.end());
    std::vector<FlatZonePtr> flatZonesList;
    flatZonesList.push_back(flatzone.get());

    ComponentTreeAdjustmentByLeaf::buildMergedAndNestedCollections(tree, flatZonesList, vflatZone.front(), newGrayLevel, isMaxtree);
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