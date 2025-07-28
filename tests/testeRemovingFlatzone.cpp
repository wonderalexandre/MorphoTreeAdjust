
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
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByFlatzone.hpp"

#include "../tests/Tests.hpp"


int main()
{

    std::cout << "Iniciando...\n";
    

    //int* img = getPassatImage(numRows, numCols);
    ImageUInt8Ptr img = getPeppersImage();
    double radioAdj = 1.5;

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, graph);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, graph);
    
    testComponentTreeFZ(maxtree, "max-tree", img);
    testComponentTreeFZ(mintree, "min-tree", img);

    //printMappingSC(maxtree);
    //std::cout <<"\n\n" << std::endl;
    //printMappingSC(mintree);
    
    ComponentTreeAdjustmentByFlatzone adjust(mintree, maxtree);

    for (NodeFZPtr node : mintree->getRoot()->getIteratorBreadthFirstTraversal()) {
        if(node == mintree->getRoot()){//} || node->getIndex() != 6){
            continue;
        }
            
        auto& flatzones = node->getCNPsByFlatZone(); 
        while (!flatzones.empty()) {
            FlatZone& flatzone = flatzones.begin()->second;
            std::cout << "\nN:" << node->getIndex()
                    << ", level:" << node->getLevel()
                    << ", numFlatzone:" << node->getNumFlatzone()
                    << ", |cnps|:" << node->getNumCNPs()
                    << std::endl;

            adjust.updateTree(maxtree, &flatzone);
            std::cout << adjust.getOutputLog() << std::endl;
            adjust.clearOutputLog();

            // Atualização da mintree
            NodeFZPtr nodeTmp = mintree->getSC(flatzone.front());
            nodeTmp->setArea(nodeTmp->getArea() - flatzone.size());
            mintree->mergeWithParent(&flatzone);        
        }    

        auto imgOutMaxtree = maxtree->reconstructionImage();
        auto imgOutMintree = mintree->reconstructionImage();
        
        testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
        testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
        if(imgOutMaxtree->isEqual(imgOutMintree))
            std::cout <<"\n✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        else
            std::cout <<"\n❌ Rec(maxtree) != Rec(mintree)" << std::endl;

        /*
        std::cout <<"\n\nApos mudancas\n" << std::endl;
        printImage(imgOutMaxtree);
        std::cout <<"\n\n" << std::endl;
        printImage(imgOutMintree);
        */    
    }
    
	std::cout << "\n\nFim do teste...\n\n";
    
    return 0;
}
