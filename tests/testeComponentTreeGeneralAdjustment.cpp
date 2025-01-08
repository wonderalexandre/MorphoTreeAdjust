
/* Debugar
1. cmake .. -DCMAKE_BUILD_TYPE=Debug
2. make
*/

#include <iostream>
#include <list>


#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeGeneralAdjustment.hpp"


void printTree(NodeCT* root, int indent = 0) {
    
    // Imprime o nó atual com indentação
    for (int i = 0; i < indent; ++i) {
        std::cout << "|-";
    }
    std::cout << "Node: " << root->getIndex() <<  ", Level: " << root->getLevel()<< std::endl;

    // Chama recursivamente a função para cada filho
    for (NodeCT* child : root->getChildren()) {
        printTree(child, indent + 1);
    }
}


int main()
{

    std::cout << "Iniciando...\n";
    
    int numRows=9;
    int numCols=9;
    int* img=new int[numRows * numCols]{
        7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 1, 1, 1, 1, 1, 7, 7, 7,
        7, 1, 1, 1, 3, 1, 7, 7, 7,
        7, 1, 1, 1, 1, 1, 7, 7, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7,
        7, 7, 7, 0, 0, 0, 0, 0, 7,
        7, 7, 7, 0, 2, 0, 2, 0, 7,
        7, 7, 7, 0, 0, 0, 0, 0, 7,
        7, 7, 7, 7, 7, 7, 7, 7, 7};
    
    int n = numRows * numCols;
    double radioAdj = 1.5;

    ComponentTree maxtree(img, numRows, numCols, true, radioAdj);
    std::cout << " --- NunNodes:" << maxtree.getNumNodes() << " ---\n";
    printTree(maxtree.getRoot());
    

    int indexPixelFZ = 48;
    std::list<int> flatZone = maxtree.getSC(indexPixelFZ)->getCNPs(); //pixels de level 0
    std::cout << "\n\nflatZone de level:" << maxtree.getSC(indexPixelFZ)->getLevel() << "\n";
    
    
    ComponentTreeGeneralAdjustment alg1;
    int newGrayLevel = maxtree.getSC(indexPixelFZ)->getLevel() + 1;
    alg1.adjustMaxTree(maxtree, flatZone, newGrayLevel);
    
    


	std::cout << "\n\nFim do teste...\n\n";

    delete[] img;
    
    return 0;
}
