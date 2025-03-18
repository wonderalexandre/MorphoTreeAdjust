#include "PyBindComponentTreeAdjustment.hpp"



void PyBindComponentTreeAdjustment::updateTree(PyBindComponentTreeFZ* tree, NodeFZ* L_leaf) {
   ComponentTreeAdjustment::updateTree(tree, L_leaf);
}


void PyBindComponentTreeAdjustment::updateTree2(PyBindComponentTreeFZ* tree, NodeFZ* rSubtree) {
   ComponentTreeAdjustment::updateTree2(tree, rSubtree);
}


py::tuple PyBindComponentTreeAdjustment::buildCollections(PyBindComponentTreeFZ* tree, std::vector<int> vflatZone, int newGrayLevel, bool isMaxtree) {
    FlatZone flatZone(vflatZone.begin(), vflatZone.end()); 
    std::vector<FlatZoneRef> flatZonesList;
    flatZonesList.push_back(flatZone);

    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, flatZonesList, newGrayLevel, isMaxtree);
    std::array<std::vector<NodeFZ*>, 256> collectionF = this->F.getCollectionF();

    std::map<int, std::vector<NodeFZ*>> mapCollectionF;
    for (int i = 0; i < 256; ++i) {
        if (!collectionF[i].empty()) {
            std::vector<NodeCT<FlatZones>*> nodes = collectionF[i];
            mapCollectionF[i] = nodes;
        }
    }

    return py::make_tuple(mapCollectionF, this->Fb);
}