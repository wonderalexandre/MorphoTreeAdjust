
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
    bool PRINTS = false;

    ImageUInt8Ptr img = getCameramanImage();
    double radioAdj = 1.5;

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(graph, true);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(graph, false);
    if(PRINTS){
        testComponentTreeFZ(maxtree, "max-tree", img);
        testComponentTreeFZ(mintree, "min-tree", img);
    }
    
    ComponentTreeAdjustmentByFlatzone adjust(mintree, maxtree);
    for (int repFlatzone: graph->getFlatzoneRepresentatives()){
        NodeFZ node = mintree->getSC(repFlatzone);
        if(node == mintree->getRoot()){//} || node.getIndex() != 6){
            continue;
        }
        //if(node.getIndex() == 10)
        //    std::cout << "AQUI";
        if(PRINTS){
            std::cout << "\n\nImprimindo as árvores antes do ajuste:\n\n";
            std::cout << "Max-tree:\n";
            printTree(maxtree->getRoot());
            std::cout << "Min-tree:\n";
            printTree(mintree->getRoot());
            
            int areaFZ = graph->getNumPixelInFlatzone(repFlatzone);
            std::cout << "\nN:" << node.getIndex()
                        << ", level:" << node.getLevel()
                        << ", numFlatzone:" << node.getNumFlatzone()
                        << ", |cnps|:" << node.getNumCNPs()
                        << ", rep:" << repFlatzone 
                        << ", areaFZ:" << areaFZ
                        << std::endl;
            printConnectedComponent(node, mintree);
            
            std::cout << "\nTauL:\n";
            printConnectedComponent(maxtree->getSC(repFlatzone), maxtree);
        }

        adjust.updateTree(maxtree, repFlatzone);
        //std::cout << adjust.getOutputLog() << std::endl;
        //adjust.clearOutputLog();

        // Atualização da mintree
        adjust.mergeWithParent(mintree, repFlatzone);

        if(PRINTS){
            std::cout << "max-tree (update):\n";
            printTree(maxtree->getRoot());
            
            std::cout << "min-tree (mergeWithParent):\n";
            printTree(mintree->getRoot());
        }

        auto imgOutMaxtree = maxtree->reconstructionImage();
        auto imgOutMintree = mintree->reconstructionImage();
        
        if(PRINTS){
            testComponentTreeFZ(maxtree, "max-tree (update)", imgOutMaxtree);
            testComponentTreeFZ(mintree, "min-tree (mergeWithParent)", imgOutMintree);
        }
        if(imgOutMaxtree->isEqual(imgOutMintree))
            std::cout <<"\n✅ Rec(maxtree) = Rec(mintree)" << std::endl;
        else{
            std::cout <<"\n❌ Rec(maxtree) != Rec(mintree)" << std::endl;
            //checkRepresentatives(maxtree, graph);
            //checkRepresentatives(mintree, graph);
            std::cout << "Max-tree:\n";
            printImage(imgOutMaxtree);
            std::cout << "Min-tree:\n";
            printImage(imgOutMintree);
            std::cout << "Max-tree:\n";
            printTree(maxtree->getRoot());
            std::cout << "Min-tree:\n";
            printTree(mintree->getRoot());

            break;
        }

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
