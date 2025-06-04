
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

    
    //auto img = getPassatImage();
    //auto img = getCharImage();
    auto img = getPeppersImage();
    //auto img = getLenaCropImage();

    int n = img->getSize();
    double radioAdj = 1.5;

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, radioAdj);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, radioAdj);
    testComponentTreeFZ(maxtree, "max-tree", img);
    testComponentTreeFZ(mintree, "min-tree", img);
    //std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    //printMappingSC(mintree);
    //std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    //printMappingSC(maxtree);
    ComponentTreeAdjustment adjust(mintree, maxtree);
    for(int threshold = 10; threshold <= 1000; threshold += 10){
        int cont = 1;
        for(NodeFZPtr rootSubtree: mintree->getNodesThreshold(threshold)){
            
            std::cout<< "\n" << cont++ << " - Processing the subtree rooted (mintree) at node id: " << rootSubtree->getIndex()<< ", numDescendants: "<< rootSubtree->computerNumDescendants()  << ", area: "<< rootSubtree->getArea() << std::endl;
            
            adjust.updateTree2(maxtree, rootSubtree);
            mintree->prunning(rootSubtree);

            auto imgOutMaxtree = maxtree->reconstructionImage();
            auto imgOutMintree = mintree->reconstructionImage();
            
            testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
            testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
            if(imgOutMaxtree->isEqual(imgOutMintree)) 
              std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
            else
              std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;
        
        }
   
    
     std::cout << "\n\n\n====================================================\n\n\n";
        cont = 1;
        for(NodeFZPtr rootSubtree: maxtree->getNodesThreshold(threshold)){
            std::cout << "\n" << cont++ << " - Processing the subtree rooted (maxtree) at node id: " << rootSubtree->getIndex()<< ", numDescendants: "<< rootSubtree->computerNumDescendants()  << ", area: "<< rootSubtree->getArea() << std::endl;
   
            adjust.updateTree2(mintree, rootSubtree);
            maxtree->prunning(rootSubtree);

            auto imgOutMaxtree = maxtree->reconstructionImage();
            auto imgOutMintree = mintree->reconstructionImage();
                    
            testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
            testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
            if(imgOutMaxtree->isEqual(imgOutMintree))
                std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
            else
                std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;

        }
    }
        
   /*
    std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    printMappingSC(mintree);
    std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    printMappingSC(maxtree);
    
    std::cout <<"\n===========  Rec(max-tree) ===========\n" << std::endl;
    printImage(imgOutMaxtree);
    std::cout <<"\n===========  Rec(min-tree) ===========\n" << std::endl;
    printImage(imgOutMintree);
    */

    //delete maxtree;
    //delete mintree;    
	std::cout << "\n\nFim do teste...\n\n";

    
    
    return 0;
}
