#include "PyBindComponentTreeAdjustmentByFlatzone.hpp"



void PyBindComponentTreeAdjustmentByFlatzone::updateTree(PyBindComponentTreeFZPtr tree, int repFlatzone) {
   ComponentTreeAdjustmentByFlatzone::updateTree(tree, repFlatzone);
}

py::tuple PyBindComponentTreeAdjustmentByFlatzone::buildCollections(PyBindComponentTreeFZPtr tree, int repFlatzone, int newGrayLevel) { 
    int idFZ = repFlatzone;
    ComponentTreeAdjustmentByFlatzone::buildMergedAndNestedCollections(tree, idFZ, newGrayLevel, tree->isMaxtree());

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
