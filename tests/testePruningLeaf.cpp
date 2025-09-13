
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
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"

#include "../tests/Tests.hpp"
#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"



int main()
{

    std::cout << "Iniciando...\n";
    
    
    
    ImageUInt8Ptr img = getWonderImage();
    printImage(img);
    
    double radioAdj = 1.5;
    
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);
    

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(graph, true);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(graph, false);

    
    

    std::cout <<"\nMintree:" << std::endl;
    printTree(mintree->getRoot());
    testComponentTreeFZ(mintree, "min-tree", img);
    //printMappingSC(mintree);

    std::cout <<"\n\nMaxtree:" << std::endl;
    printTree(maxtree->getRoot());
    testComponentTreeFZ(maxtree, "max-tree", img);
   // printMappingSC(maxtree);

    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);
    ComponentTreeFZPtr tree1, tree2;
    tree1 = tree2 = nullptr;
    int numNodes = mintree->getNumNodes();
    for(int i=0; i < numNodes; i++){
        if(i % 2 == 0){
            tree1 = mintree;
            tree2 = maxtree;
        }else{
            tree1 = maxtree;
            tree2 = mintree;
        }
        
        NodeId L_leafId = tree1->getLeaves().front();   
        NodeCT L_leaf = tree1->proxy(L_leafId);
        if(L_leaf == tree1->getRoot()) break;
        
        if(L_leaf.getIndex() == 11){
            std::cout << "===============================================\n";
            testComponentTreeFZ(tree2, "tree (pruning)", tree2->reconstructionImage());
            printImage(tree2->reconstructionImage());
            printTree(tree2->getRoot());    
            std::cout << "===============================================\n";

        }

        std::cout << "\n\n---Pruning (" << (tree1->isMaxtree()? "max-tree":"min-tree") << ") node:" << L_leaf.getIndex() << ", level:" << L_leaf.getLevel() << ", |cnps|:" << L_leaf.getNumCNPs() << ", idFlatzone:"<< L_leaf.getRepCNPs().front()  << std::endl;
        //graph->printUnionFind();
        adjust.updateTree(tree2, L_leaf);
        //graph->printUnionFind();
        adjust.prunning(tree1, L_leaf);            
        //checkRepresentatives(tree1, graph);
        //checkRepresentatives(tree2, graph);
        

        ImageUInt8Ptr imgTree2 = tree2->reconstructionImage();
        ImageUInt8Ptr imgTree1 = tree1->reconstructionImage();
        
        printImage(imgTree2);
        printTree(tree2->getRoot());                    
        testComponentTreeFZ(tree2, "tree (update)", imgTree2);
        
        printImage(imgTree1);
        printTree(tree1->getRoot());                    
        testComponentTreeFZ(tree1, "tree (pruning)", imgTree1);
        

        bool isEqual =  imgTree2->isEqual(imgTree1);
        std::cout << "The images are equals: " << ( isEqual? "True":"False") << "\n";
        printImage(imgTree1);
        if(!isEqual){
            
            AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
            ComponentTreePPtr maxtreeCorrect = std::make_shared<ComponentTreeP>(imgTree1, true, adj);
            std::cout << "Correct tree\n";
            printTree(maxtreeCorrect->getRoot());

            std::cout << "Maxtree reconstruction image:\n";
            printImage(imgTree2);
            std::cout << "\n\n========Falhou==========\n\n";
            break;
        }
    
        std::cout <<"\n Tree1 (type: " << (tree1->isMaxtree()? "max-tree":"min-tree") << ")" << std::endl;
        //printTree(tree1->getRoot());
        
        std::cout <<"\n\nTree2 (" << (tree2->isMaxtree()? "max-tree":"min-tree") << ")" << std::endl;
        //printTree(tree2->getRoot());
      
    }
	std::cout << "\n\nFim do teste...\n\n";

    return 0;
}
