
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
    ImageUInt8Ptr img = getSmallImage();
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
    
    NodeFZPtr N = nullptr;
    int index = 2;
    for (NodeFZPtr node : mintree->getRoot()->getIteratorBreadthFirstTraversal()) {
        if(node->getIndex() == index){
            N = node;
            break;
        }
    }    
    

    ComponentTreeAdjustmentByFlatzone adjust(mintree, maxtree);
     

    FlatZone flatzone = N->getCNPsByFlatZone().begin()->second;
    std::cout <<"\nN:" << N->getIndex() << ", level:" << N->getLevel() << ", |cnps|:" << N->getNumCNPs() <<  std::endl;
    adjust.updateTree(maxtree, &flatzone);
    std::cout << adjust.getOutputLog() << std::endl;   
    mintree->mergeWithParent(&flatzone);

    auto imgOutMaxtree = maxtree->reconstructionImage();
    auto imgOutMintree = mintree->reconstructionImage();
    
    testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
    testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
    if(imgOutMaxtree->isEqual(imgOutMintree))
        std::cout <<"\n✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    else
        std::cout <<"\n❌ Rec(maxtree) != Rec(mintree)" << std::endl;


    std::cout <<"\n\nApos mudancas\n" << std::endl;
    printImage(imgOutMaxtree);
    std::cout <<"\n\n" << std::endl;
    printImage(imgOutMintree);
    
	std::cout << "\n\nFim do teste...\n\n";
    
    return 0;
}
