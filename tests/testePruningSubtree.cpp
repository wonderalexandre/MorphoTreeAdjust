
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

    int numRows,numCols;
    int* img = getPassatImage(numRows, numCols);

    int n = numRows * numCols;
    double radioAdj = 1.5;

    ComponentTreeFZ* maxtree = new ComponentTreeFZ(img, numRows, numCols, true, radioAdj);
    ComponentTreeFZ* mintree = new ComponentTreeFZ(img, numRows, numCols, false, radioAdj);
    std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    printMappingSC(mintree);
    
    std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    printMappingSC(maxtree);
    

    
    NodeFZ* TauS_r = getNodeByIndex(maxtree, 7);
    int id_TauS_r = TauS_r->getIndex();
    ComponentTreeAdjustment adjust(mintree, maxtree);
    

    adjust.updateTree2(mintree, TauS_r);
    maxtree->prunning(TauS_r);

    int* imgOutMaxtree = maxtree->reconstructionImage();
    int* imgOutMintree = mintree->reconstructionImage();
        
        
    std::cout << "NunNodes (mintree):" << mintree->getNumNodes() << "\tNumNodes (maxtree):" << maxtree->getNumNodes();
    std::cout <<"TauS_r(id): " << id_TauS_r << "\t Rec(maxtree) = Rec(mintree):" << isEquals(imgOutMaxtree, imgOutMintree, n) << std::endl;
    

    std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    printMappingSC(mintree);
    std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    printMappingSC(maxtree);
    
    
    std::cout <<"\n===========  Rec(max-tree) ===========\n" << std::endl;
    printImage(imgOutMaxtree, numRows, numCols);
    std::cout <<"\n===========  Rec(min-tree) ===========\n" << std::endl;
    printImage(imgOutMintree, numRows, numCols);
    
    delete[] imgOutMaxtree;
    delete[] imgOutMintree;
    delete maxtree;
    delete mintree;    
	std::cout << "\n\nFim do teste...\n\n";

    delete[] img;
    
    return 0;
}
