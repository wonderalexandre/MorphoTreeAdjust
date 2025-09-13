
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
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentBySubtree.hpp"


#include "../tests/Tests.hpp"

int main()
{

    
    //auto img = getPassatImage();
    //auto img = getCharImage();
    auto img = getWonderImage();
    printImage(img);

    double radioAdj = 1.5;

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(graph, true);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(graph, false);
    //testComponentTreeFZ(maxtree, "max-tree", img);
    //testComponentTreeFZ(mintree, "min-tree", img);
    //std::cout <<"\n=========== mapIDs min-tree ===========\n" << std::endl;
    //printMappingSC(mintree);
    std::cout <<"\n=========== mapIDs max-tree ===========\n" << std::endl;
    printMappingSC(maxtree);

    std::cout <<"\nMintree:" << std::endl;
    printTree(mintree->getRoot());

    std::cout <<"\n\nMaxtree:" << std::endl;
    printTree(maxtree->getRoot());


    ComponentTreeAdjustmentBySubtree adjust(mintree, maxtree);


/*
    for (NodeFZPtr node : mintree->getRoot()->getIteratorBreadthFirstTraversal()) {
        std::cout<< "Processing the subtree rooted (mintree) at node id: " << node->getIndex()<< ", numDescendants: "<< node->computerNumDescendants()  << ", area: "<< node->getArea() << std::endl;
        printConnectedComponent(node, mintree);

    }

    NodeFZPtr N = nullptr;
    int index = 11;
    for (NodeFZPtr node : mintree->getRoot()->getIteratorBreadthFirstTraversal()) {
        if(node->getIndex() == index){
            N = node;
            break;
        }
    }    
    
    adjust.updateTree(maxtree, N);
    std::cout << adjust.getOutputLog() << std::endl;   
    mintree->prunning(N);

    auto imgOutMaxtree = maxtree->reconstructionImage();
    auto imgOutMintree = mintree->reconstructionImage();
        
    testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
    testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
    if(imgOutMaxtree->isEqual(imgOutMintree))
        std::cout <<"\n✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    else
        std::cout <<"\n❌ Rec(maxtree) != Rec(mintree)" << std::endl;
*/
    
    for(int threshold = 50; threshold <= 100; threshold += 30){
        int cont = 1;
        for(NodeId rootSubtreeId: ComponentTreeFZ::getNodesThreshold(mintree, threshold)){
            NodeFZ rootSubtree = mintree->proxy(rootSubtreeId);
            std::cout<< "\n" << cont++ << " - Processing the subtree rooted (mintree) at node id: " << rootSubtree.getIndex()<< ", numDescendants: "<< rootSubtree.computerNumDescendants()  << ", area: "<< rootSubtree.getArea() << std::endl;
            printConnectedComponent(rootSubtree, mintree);

            
            adjust.updateTree(maxtree, rootSubtree);
            std::cout << adjust.getOutputLog() << std::endl;
            adjust.prunning(mintree, rootSubtree);

            checkRepresentatives(mintree, graph);
            checkRepresentatives(maxtree, graph);

            auto imgOutMaxtree = maxtree->reconstructionImage();
            auto imgOutMintree = mintree->reconstructionImage();
            
            //testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
            //testComponentTreeFZ(mintree, "min-tree", imgOutMintree);
            if(imgOutMaxtree->isEqual(imgOutMintree)) 
              std::cout <<"✅ Rec(maxtree) = Rec(mintree)" << std::endl;
            else{
              std::cout <<"❌ Rec(maxtree) != Rec(mintree)" << std::endl;
                ComponentTreeFZPtr maxtreeCorrect = std::make_shared<ComponentTreeFZ>(std::make_shared<FlatZonesGraph>(imgOutMintree, adj), true);
                std::cout << "Correct tree\n";
                printTree(maxtreeCorrect->getRoot());
            }
        


            std::cout <<"\n Min-tree" << std::endl;
            printTree(mintree->getRoot());
            
            std::cout <<"\n\nMax-tree" << std::endl;
            printTree(maxtree->getRoot());

            std::cout << "Maxtree reconstruction image:\n";
            printImage(imgOutMaxtree);
            std::cout << "\n\nMintree reconstruction image:\n";
            printImage(imgOutMintree);

            break;
        }
   
    
     std::cout << "\n\n\n====================================================\n\n\n";
        cont = 1;
        for(NodeId rootSubtreeId: ComponentTreeFZ::getNodesThreshold(maxtree, threshold)){
            NodeFZ rootSubtree = maxtree->proxy(rootSubtreeId);
            std::cout << "\n" << cont++ << " - Processing the subtree rooted (maxtree) at node id: " << rootSubtree.getIndex()<< ", numDescendants: "<< rootSubtree.computerNumDescendants()  << ", area: "<< rootSubtree.getArea() << std::endl;
   
            adjust.updateTree(mintree, rootSubtree);
            adjust.prunning(maxtree, rootSubtree);

            auto imgOutMaxtree = maxtree->reconstructionImage();
            auto imgOutMintree = mintree->reconstructionImage();
            
            printTree(maxtree->getRoot());                    
            testComponentTreeFZ(maxtree, "max-tree", imgOutMaxtree);
            
            printTree(mintree->getRoot());                    
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
