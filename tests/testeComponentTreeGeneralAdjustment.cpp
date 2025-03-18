
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


int main(){

    std::cout << "Iniciando...\n";
    return 0;
}
    /*
    int numRows=10;
    int numCols=17;
    int* img=new int[numRows * numCols]{
        6, 6, 6, 1, 1, 0, 1, 6, 6, 5, 2, 1, 1, 1, 0, 7, 7,
        6, 1, 6, 1, 0, 0, 1, 6, 1, 5, 2, 2, 1, 1, 0, 3, 3,
        6, 4, 0, 0, 6, 0, 1, 5, 5, 5, 1, 2, 1, 1, 0, 3, 3,
        6, 1, 0, 6, 6, 6, 6, 0, 2, 1, 5, 3, 1, 1, 0, 5, 5,
        6, 1, 0, 6, 0, 0, 6, 2, 2, 2, 8, 8, 1, 0, 1, 5, 5,
        6, 1, 0, 6, 6, 6, 1, 2, 2, 2, 3, 8, 0, 0, 4, 1, 1,
        1, 1, 3, 0, 0, 0, 7, 7, 3, 8, 8, 0, 4, 4, 1, 4, 4,
        6, 3, 4, 3, 7, 0, 7, 6, 6, 6, 0, 2, 4, 1, 6, 1, 4,
        6, 3, 4, 3, 7, 0, 2, 3, 3, 0, 2, 2, 4, 1, 5, 1, 4,
        1, 3, 3, 3, 1, 0, 2, 3, 3, 0, 2, 2, 4, 1, 1, 1, 4
    };
    
    int n = numRows * numCols;
    double radioAdj = 1.5;

    ComponentTree maxtree(img, numRows, numCols, true, radioAdj);
    std::cout << " --- NunNodes:" << maxtree.getNumNodes() << " ---\n";
    printTree(maxtree.getRoot());
    

    int indexPixelFZ = 76; //particao do paper
    std::list<int> flatZone = maxtree.getSC(indexPixelFZ)->getCNPs(); //pixels de level 0
    std::cout << "\n\nflatZone de level:" << maxtree.getSC(indexPixelFZ)->getLevel() << "\n";
    
    
    ComponentTreeGeneralAdjustment alg1;
    int newGrayLevel = 7;
    alg1.adjustMaxTree(maxtree, flatZone, newGrayLevel);
    
    


	std::cout << "\n\nFim do teste...\n\n";

    delete[] img;
    */
    
