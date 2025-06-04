
/* Debugar
1. cmake .. -DCMAKE_BUILD_TYPE=Debug
2. make
*/

#include <iostream>
#include <list>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>


#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustment.hpp"

#include "../tests/Tests.hpp"


int main()
{

    std::cout << "Iniciando...\n";
    

    //int* img = getPassatImage(numRows, numCols);
    ImageUInt8Ptr img = getPassatImage();
    double radioAdj = 1.5;

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, radioAdj);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, radioAdj);
    
    testComponentTreeFZ(maxtree, "max-tree", img);
    testComponentTreeFZ(mintree, "min-tree", img);

    //printMappingSC(maxtree);
    //std::cout <<"\n\n" << std::endl;
    //printMappingSC(mintree);
    
    NodeFZPtr N = nullptr;
    int index = 5;
    for (NodeFZPtr node : mintree->getRoot()->getIteratorBreadthFirstTraversal()) {
        if(node->getIndex() == index){
            N = node;
            break;
        }
    }    

    ComponentTreeAdjustment adjust(mintree, maxtree);
     

    //for(int i=0; i < numNodes; i++){
    //    NodeCT* L_leaf = mintree->getLeaves().front();
       // if(L_leaf->getIndex() != 293) continue;
    //for(NodeFZPtr L_leaf : mintree->getRoot()->getIteratorPostOrderTraversal()){
        
        std::cout <<"\nN:" << N->getIndex() << ", level:" << N->getLevel() << ", |cnps|:" << N->getNumCNPs() <<  std::endl;
        adjust.updateTree3(maxtree, N);
        mintree->mergeWithParent(N);

        auto imgOutMaxtree = maxtree->reconstructionImage();
        auto imgOutMintree = mintree->reconstructionImage();
        
        testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
        testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
        if(imgOutMaxtree->isEqual(imgOutMintree))
            std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        else
            std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;

    
    //}
    std::cout <<"\n\nApos mudancas\n" << std::endl;
    printImage(imgOutMaxtree);
    std::cout <<"\n\n" << std::endl;
    printImage(imgOutMintree);
    
	std::cout << "\n\nFim do teste...\n\n";
    
    return 0;
}
