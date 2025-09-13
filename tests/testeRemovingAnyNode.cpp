
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
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByAnyNode.hpp"

#include "../tests/Tests.hpp"


int main()
{

    std::cout << "Iniciando...\n";
    bool PRINTS = false;

    ImageUInt8Ptr img = getCameramanImage();
    double radioAdj = 1.5;

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(graph, true);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(graph, false);
    if(PRINTS){
        printImage(img);
        testComponentTreeFZ(maxtree, "max-tree", img);
        testComponentTreeFZ(mintree, "min-tree", img);
    }
    //printMappingSC(maxtree);
    //std::cout <<"\n\n" << std::endl;
    //printMappingSC(mintree);


    ComponentTreeAdjustmentByAnyNode adjust(mintree, maxtree);
    for (int rep: graph->getFlatzoneRepresentatives()){
        NodeFZ node = mintree->getSC(rep);
        if(node == mintree->getRoot()){//} || node.getIndex() != 6){
            continue;
        }
        if(PRINTS){
            std::cout << "\n\nImprimindo as árvores antes do ajuste:\n\n";
            std::cout << "Max-tree:\n";
            printTree(maxtree->getRoot());
            std::cout << "Min-tree:\n";
            printTree(mintree->getRoot());
            
            printConnectedComponent(maxtree->getSC(rep), maxtree);
            std::cout << "\nN:" << node.getIndex()
                        << ", level:" << node.getLevel()
                        << ", numFlatzone:" << node.getNumFlatzone()
                        << ", |cnps|:" << node.getNumCNPs()
                        << std::endl;
            printConnectedComponent(node, mintree);

            std::cout << "Max-tree:\n";
            printTree(maxtree->getRoot());
            std::cout << "Min-tree:\n";
            printTree(mintree->getRoot());
        }

        adjust.updateTree(maxtree, node);
        adjust.mergeWithParent(mintree, node);
        
        auto imgOutMaxtree = maxtree->reconstructionImage();
        auto imgOutMintree = mintree->reconstructionImage();
        if(PRINTS){
            testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
            testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
        }
        if(imgOutMaxtree->isEqual(imgOutMintree)){
            std::cout <<"\n✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        }
        else{
            std::cout <<"\n❌ Rec(maxtree) != Rec(mintree)" << std::endl;
            //checkRepresentatives(maxtree, graph);
            //checkRepresentatives(mintree, graph);
            printImage(imgOutMaxtree);
            printImage(imgOutMintree);
            std::cout << "Max-tree:\n";
            printTree(maxtree->getRoot());
            std::cout << "Min-tree:\n";
            printTree(mintree->getRoot());

            break;
        }

    
    }
    
	std::cout << "\n\nFim do teste...\n\n";
    
    return 0;
}
