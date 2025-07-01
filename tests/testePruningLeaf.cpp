
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
    
    
    
    ImageUInt8Ptr img = getPeppersImage();

    double radioAdj = 1.5;
    
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);
    

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, graph);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, graph);

    testComponentTreeFZ(maxtree, "max-tree", img);
    testComponentTreeFZ(mintree, "min-tree", img);

    
    
    std::cout <<"\nMintree:" << std::endl;
    printTree(mintree->getRoot());
    printMappingSC(mintree);

    std::cout <<"\n\nMaxtree:" << std::endl;
    printTree(maxtree->getRoot());
    printMappingSC(maxtree);

    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);

    int numNodes = maxtree->getNumNodes();
    for(int i=0; i < numNodes; i++){
        NodeFZPtr L_leaf = mintree->getLeaves().front();   
        if(L_leaf == mintree->getRoot()) break;
        std::cout << "\n\n---Pruning (min-tree) node:" << L_leaf->getIndex() << ", level:" << L_leaf->getLevel() << ", |cnps|:" << L_leaf->getNumCNPs() << std::endl;
        adjust.updateTree(maxtree, L_leaf);
        mintree->prunning(L_leaf);            

        testComponentTreeFZ(maxtree, "max-tree", maxtree->reconstructionImage());
        testComponentTreeFZ(mintree, "min-tree", mintree->reconstructionImage());

        
    
      
    }
	std::cout << "\n\nFim do teste...\n\n";

    return 0;
}
