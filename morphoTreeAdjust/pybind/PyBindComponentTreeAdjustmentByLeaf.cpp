#include "PyBindComponentTreeAdjustmentByLeaf.hpp"



void PyBindComponentTreeAdjustmentByLeaf::updateTree(PyBindComponentTreeFZPtr tree, NodeFZ<> L_leaf) {
   PyBindAdjustmentLeafBase::updateTree(tree.get(), L_leaf.getIndex());
}

py::tuple PyBindComponentTreeAdjustmentByLeaf::buildCollections(PyBindComponentTreeFZPtr tree, NodeFZ<> leaf, int newGrayLevel) { 
    // representante (pixel) da flat-zone associada à folha
    int idLeaf = leaf.getRepCNPs().empty() ? -1 : leaf.getRepCNPs().front();
    if (idLeaf < 0) {
        return py::make_tuple(std::map<int, std::vector<NodeFZ<>>> {}, std::vector<NodeId>{});
    }
    (void)tree->getSC(idLeaf); // força validação do SC (não utilizado aqui diretamente)
    PyBindAdjustmentLeafBase::buildMergedAndNestedCollections(tree.get(), idLeaf, newGrayLevel, tree->isMaxtree());

    auto &collectionF = this->F.getCollectionF();
    std::map<int, std::vector<NodeFZ<>>> mapCollectionF;
    for (int i = 0; i < 256; ++i) {
        if (!collectionF[i].empty()) {
            std::vector<NodeFZ<>> nodes;
            nodes.reserve(collectionF[i].size());
            for (NodeId nid : collectionF[i]) nodes.push_back(tree->proxy(nid));
            mapCollectionF[i] = std::move(nodes);
        }
    }
    
    std::vector<NodeFZ<>> Fb;
    for(NodeId nid : this->F.getFb()){
        Fb.push_back(tree->proxy(nid));
    }

    return py::make_tuple(mapCollectionF, Fb);
}
