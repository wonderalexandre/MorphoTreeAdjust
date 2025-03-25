
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
    int numRows, numCols;

    //int* img = getPassatImage(numRows, numCols);
    int* img = getLenaCropImage(numRows, numCols);
    int n = numRows * numCols;
    double radioAdj = 1.5;

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, numRows, numCols, true, radioAdj);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, numRows, numCols, false, radioAdj);

    testComponentTreeFZ(maxtree, "max-tree", img, numRows, numCols);
    testComponentTreeFZ(mintree, "min-tree", img, numRows, numCols);

    //printMappingSC(maxtree);
    //std::cout <<"\n\n" << std::endl;
    //printMappingSC(mintree);
    

    ComponentTreeAdjustment adjust(mintree, maxtree);
    int numNodes = mintree->getNumNodes()-1;


    //for(int i=0; i < numNodes; i++){
    //    NodeCT* L_leaf = mintree->getLeaves().front();
       // if(L_leaf->getIndex() != 293) continue;
    for(NodeFZPtr L_leaf : mintree->getRoot()->getIteratorPostOrderTraversal()){
        if(L_leaf == mintree->getRoot()) break;
        
        int id = L_leaf->getIndex();
        std::cout <<"\nPruning L_leaf:" << id << ", level:" << L_leaf->getLevel() << ", |cnps|:" << L_leaf->getNumCNPs() <<  std::endl;
        adjust.updateTree2(maxtree, L_leaf);
        mintree->prunning(L_leaf);

        int* imgOutMaxtree = maxtree->reconstructionImage();
        int* imgOutMintree = mintree->reconstructionImage();
        
        testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree, numRows, numCols);
        testComponentTreeFZ(mintree, "min-tree", imgOutMintree, numRows, numCols);
        if(isEquals(imgOutMaxtree, imgOutMintree, n))                
            std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        else
            std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;

        delete[] imgOutMaxtree;
        delete[] imgOutMintree;
    
    }
   // std::cout <<"\n\nApos mudancas\n" << std::endl;
   // printMappingSC(maxtree);
   // std::cout <<"\n\n" << std::endl;
   // printMappingSC(mintree);
   // delete maxtree;
   // delete mintree;    
	std::cout << "\n\nFim do teste...\n\n";

    delete[] img;
    
    return 0;
}
