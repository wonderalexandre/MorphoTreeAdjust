
#include "./PyBindComponentTreeAdjustment.hpp"

void PyBindComponentTreeAdjustment::updateTree(PyBindComponentTree* tree, NodeCT *L_leaf){
   ComponentTreeAdjustment::updateTree(tree, L_leaf);
}

py::tuple PyBindComponentTreeAdjustment::buildCollections(PyBindComponentTree* tree, std::vector<int> flatZone, int newGrayLevel, bool isMaxtree){
      std::list<int> flatZoneList(flatZone.begin(), flatZone.end()); 
      ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, flatZoneList, newGrayLevel, isMaxtree);

      std::array<std::vector<NodeCT*>, 256> collectionF = this->F.getCollectionF();

      std::map<int, std::vector<NodeCT*>> mapCollectionF;
      for (int i = 0; i < 256; ++i) {
         if (!collectionF[i].empty()) {
            std::vector<NodeCT*> nodes = collectionF[i];
            mapCollectionF[i] = nodes;
         }
      }

      return py::make_tuple(mapCollectionF, this->B_L);

}