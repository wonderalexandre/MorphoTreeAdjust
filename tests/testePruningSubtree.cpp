
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
    //int* img = getPassatImage(numRows, numCols);
    //int* img = getCharImage(numRows, numCols);
    int* img = getPeppersImage(numRows, numCols);
    //int* img = getLenaCropImage(numRows, numCols);

    int n = numRows * numCols;
    double radioAdj = 1.5;

    ComponentTreeFZ* maxtree = new ComponentTreeFZ(img, numRows, numCols, true, radioAdj);
    ComponentTreeFZ* mintree = new ComponentTreeFZ(img, numRows, numCols, false, radioAdj);
    testComponentTreeFZ(maxtree, "max-tree", img, numRows, numCols);
    testComponentTreeFZ(mintree, "min-tree", img, numRows, numCols);
    //std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    //printMappingSC(mintree);
    //std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    //printMappingSC(maxtree);
    ComponentTreeAdjustment adjust(mintree, maxtree);
    for(int threshold = 10; threshold <= 1000; threshold += 10){
        int cont = 1;
        for(NodeFZ* rootSubtree: mintree->getNodesThreshold(threshold)){
            
            std::cout<< "\n" << cont++ << " - Processing the subtree rooted (mintree) at node id: " << rootSubtree->getIndex()<< ", numDescendants: "<< rootSubtree->computerNumDescendants()  << ", area: "<< rootSubtree->getArea() << std::endl;
            
            adjust.updateTree2(maxtree, rootSubtree);
            mintree->prunning(rootSubtree);

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
   
    
     std::cout << "\n\n\n====================================================\n\n\n";
        cont = 1;
        for(NodeFZ* rootSubtree: maxtree->getNodesThreshold(threshold)){
            std::cout << "\n" << cont++ << " - Processing the subtree rooted (maxtree) at node id: " << rootSubtree->getIndex()<< ", numDescendants: "<< rootSubtree->computerNumDescendants()  << ", area: "<< rootSubtree->getArea() << std::endl;
   
            adjust.updateTree2(mintree, rootSubtree);
            maxtree->prunning(rootSubtree);

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
    }
        
   /*
    std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    printMappingSC(mintree);
    std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    printMappingSC(maxtree);
    
    std::cout <<"\n===========  Rec(max-tree) ===========\n" << std::endl;
    printImage(imgOutMaxtree, numRows, numCols);
    std::cout <<"\n===========  Rec(min-tree) ===========\n" << std::endl;
    printImage(imgOutMintree, numRows, numCols);
    */

    delete maxtree;
    delete mintree;    
	std::cout << "\n\nFim do teste...\n\n";

    delete[] img;
    
    return 0;
}
