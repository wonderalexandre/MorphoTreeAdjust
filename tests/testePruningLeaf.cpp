
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
#include "../morphoTreeAdjust/include/ComponentTreeAdjustment.hpp"

#include "../tests/Tests.hpp"


int main()
{

    std::cout << "Iniciando...\n";
    
    auto img = getLenaCropImage();
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
    

    ComponentTreeAdjustment adjust(mintree, maxtree);


    //for(int i=0; i < numNodes; i++){
    //    NodeCT* L_leaf = mintree->getLeaves().front();
       // if(L_leaf->getIndex() != 293) continue;
    for(NodeFZPtr L_leaf : mintree->getRoot()->getIteratorPostOrderTraversal()){
        if(L_leaf == mintree->getRoot()) break;
        
        int id = L_leaf->getIndex();
        std::cout <<"\nPruning L_leaf:" << id << ", level:" << L_leaf->getLevel() << ", |cnps|:" << L_leaf->getNumCNPs() <<  std::endl;
        adjust.updateTree(maxtree, L_leaf);
        mintree->prunning(L_leaf);

        auto imgOutMaxtree = maxtree->reconstructionImage();
        auto imgOutMintree = mintree->reconstructionImage();
        
        testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
        testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
        if(imgOutMaxtree->isEqual(imgOutMintree))
            std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        else
            std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;

    }
	std::cout << "\n\nFim do teste...\n\n";

    return 0;
}
