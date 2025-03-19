
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
    int* img = getPassatImage(numRows, numCols);
    int n = numRows * numCols;
    double radioAdj = 1.5;

    ComponentTreeFZ* maxtree = new ComponentTreeFZ(img, numRows, numCols, true, radioAdj);
    ComponentTreeFZ* mintree = new ComponentTreeFZ(img, numRows, numCols, false, radioAdj);
    
    testComponentTreeFZ(maxtree, "max-tree", img, numRows, numCols);
    testComponentTreeFZ(mintree, "min-tree", img, numRows, numCols);

    //printMappingSC(maxtree);
    //std::cout <<"\n\n" << std::endl;
    //printMappingSC(mintree);
    

    ComponentTreeAdjustment adjust(mintree, maxtree);
    int numNodes = mintree->getNumNodes()-1;

    NodeFZ* N = mintree->getSC(2327);

    //for(int i=0; i < numNodes; i++){
    //    NodeCT* L_leaf = mintree->getLeaves().front();
       // if(L_leaf->getIndex() != 293) continue;
    //for(NodeFZ* L_leaf : mintree->getRoot()->getIteratorPostOrderTraversal()){
        
        std::cout <<"\nN:" << N->getIndex() << ", level:" << N->getLevel() << ", |cnps|:" << N->getNumCNPs() <<  std::endl;
        adjust.updateTree2(maxtree, N);
        mintree->mergeParent(N);

        int* imgOutMaxtree = maxtree->reconstructionImage();
        int* imgOutMintree = mintree->reconstructionImage();
        
        testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree, numRows, numCols);
        testComponentTreeFZ(mintree, "min-tree", imgOutMintree, numRows, numCols);
        if(isEquals(imgOutMaxtree, imgOutMintree, n))                
            std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        else
            std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;

    
    //}
    std::cout <<"\n\nApos mudancas\n" << std::endl;
    printImage(imgOutMaxtree, numRows, numCols);
    std::cout <<"\n\n" << std::endl;
    printImage(imgOutMintree, numRows, numCols);
    
    delete[] imgOutMaxtree;
    delete[] imgOutMintree;

    
    delete maxtree;
    delete mintree;    
	std::cout << "\n\nFim do teste...\n\n";

    delete[] img;
    
    return 0;
}
