#include "PyBindComponentTreeAdjustment.hpp"



void PyBindComponentTreeAdjustment::updateTree(PyBindComponentTreeFZPtr tree, NodeFZPtr L_leaf) {
   ComponentTreeAdjustment::updateTree(tree, L_leaf);
}


void PyBindComponentTreeAdjustment::updateTree2(PyBindComponentTreeFZPtr tree, NodeFZPtr rSubtree) {
   ComponentTreeAdjustment::updateTree2(tree, rSubtree);
}


py::tuple PyBindComponentTreeAdjustment::buildCollections(PyBindComponentTreeFZPtr tree, std::vector<int> vflatZone, int newGrayLevel, bool isMaxtree) {
    FlatZone flatZone(vflatZone.begin(), vflatZone.end()); 
    std::vector<FlatZoneRef> flatZonesList;
    flatZonesList.push_back(flatZone);

    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, flatZonesList, newGrayLevel, isMaxtree);
    std::array<std::vector<NodeFZPtr>, 256> collectionF = this->F.getCollectionF();

    std::map<int, std::vector<NodeFZPtr>> mapCollectionF;
    for (int i = 0; i < 256; ++i) {
        if (!collectionF[i].empty()) {
            std::vector<NodeFZPtr> nodes = collectionF[i];
            mapCollectionF[i] = nodes;
        }
    }

    return py::make_tuple(mapCollectionF, this->Fb);
}