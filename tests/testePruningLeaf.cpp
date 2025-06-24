
/* Debugar
1. cmake .. -DCMAKE_BUILD_TYPE=Debug
2. make
*/

#include <iostream>
#include <list>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"

#include "../tests/Tests.hpp"
#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"

int main()
{

    std::cout << "Iniciando...\n";
    std::string filename = "/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/dat/icdar_resized/val_091.png";
    
    
    int numCols, numRows, nchannels;
    uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    
    if (!data) {
        std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
        return 1;
    }
    std::cout << "Resolution: " << numCols << "x" << numRows << std::endl;
    ImageUInt8Ptr img = ImageUInt8::fromRaw(data, numRows, numCols);
    //img = getCharImage();

    double radioAdj = 1.5;
    
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::unique_ptr<FlatZonesGraph> maxTreeFZGraph = std::make_unique<FlatZonesGraph>(img, adj);
    std::unique_ptr<FlatZonesGraph> minTreeFZGraph = std::make_unique<FlatZonesGraph>(*maxTreeFZGraph);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, std::move(maxTreeFZGraph));
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, std::move(minTreeFZGraph));

    testComponentTreeFZ(maxtree, "max-tree", img);
    testComponentTreeFZ(mintree, "min-tree", img);

    //printMappingSC(maxtree);
    //std::cout <<"\n\n" << std::endl;
    //printMappingSC(mintree);
    
    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);


    //for(int i=0; i < numNodes; i++){
    //    NodeCT* L_leaf = mintree->getLeaves().front();
       // if(L_leaf->getIndex() != 293) continue;
    auto nodesToPruning = maxtree->getNodesThreshold(60);
    adjust.adjustMinTree(mintree, maxtree, nodesToPruning);

    for(NodeFZPtr node : mintree->getNodesThreshold(60)){
        std::cout << "\n\n---Pruning node:" << node->getIndex() << ", level:" << node->getLevel() << ", |cnps|:" << node->getNumCNPs() << std::endl;
        testComponentTreeFZ(maxtree, "max-tree", maxtree->reconstructionImage());
        testComponentTreeFZ(mintree, "min-tree", mintree->reconstructionImage());
        for (NodeFZPtr L_leaf : node->getIteratorPostOrderTraversal()) { 
            if(L_leaf == mintree->getRoot()) break;
            
            adjust.updateTree(maxtree, L_leaf);
            mintree->prunning(L_leaf);            
        }

        std::cout << "\n\n+++Pruning node:" << node->getIndex() << ", level:" << node->getLevel() << ", |cnps|:" << node->getNumCNPs() << std::endl;
        testComponentTreeFZ(maxtree, "max-tree", maxtree->reconstructionImage());
        testComponentTreeFZ(mintree, "min-tree", mintree->reconstructionImage());

    }
	std::cout << "\n\nFim do teste...\n\n";

    return 0;
}
