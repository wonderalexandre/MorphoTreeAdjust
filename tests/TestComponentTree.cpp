#include "../tests/Tests.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustment.hpp"

int main() {
    // Definição da imagem e parâmetros
    int numRows, numCols;
    int* img = getPassatImage(numRows, numCols);
    int n = numRows * numCols;
    double radioAdj = 1.5;

    // Criação das Component Trees
    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, numRows, numCols, true, radioAdj);
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, numRows, numCols, false, radioAdj);

    // Executar testes
    testComponentTreeFZ(mintree, "Min-Tree", mintree->reconstructionImage(), numRows, numCols);
    testComponentTreeFZ(maxtree, "Max-Tree", maxtree->reconstructionImage(), numRows, numCols);

    NodeFZPtr Lmin_leaf = mintree->getLeaves().front();
    std::cout << "---- Pruning leaf id: " << Lmin_leaf->getIndex() << " parent id: " << Lmin_leaf->getParent()->getIndex() << " ----" << std::endl;
    mintree->prunning(Lmin_leaf);
    int* imgMintree = mintree->reconstructionImage();
    testComponentTreeFZ(mintree, "Min-Tree after pruning", imgMintree, numRows, numCols);

    NodeFZPtr Lmax_leaf = maxtree->getLeaves().front();
    std::cout << "---- Pruning leaf id: " << Lmax_leaf->getIndex() << " parent id: " << Lmax_leaf->getParent()->getIndex() << " ----" << std::endl;
    maxtree->prunning(Lmax_leaf);
    int* imgMaxtree = maxtree->reconstructionImage();
    testComponentTreeFZ(maxtree, "Max-Tree after pruning", imgMaxtree, numRows, numCols);

    delete[] imgMaxtree;
    delete[] imgMintree;
    //delete maxtree;
    //delete mintree;

    maxtree = std::make_shared<ComponentTreeFZ>(img, numRows, numCols, true, radioAdj);
    mintree = std::make_shared<ComponentTreeFZ>(img, numRows, numCols, false, radioAdj);
    ComponentTreeAdjustment adjust(mintree, maxtree);

    NodeFZPtr Lmin_leaf1 = mintree->getLeaves().front();
    std::cout << "Pruning id: " << Lmin_leaf1->getIndex() << std::endl;
    adjust.updateTree(maxtree, Lmin_leaf1);
    mintree->prunning(Lmin_leaf1);
    imgMaxtree = maxtree->reconstructionImage();
    imgMintree = mintree->reconstructionImage();

    if (isEquals(imgMaxtree, imgMintree, n)) {
        std::cout << "✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    } else {
        std::cout << "❌ Rec(maxtree) ≠ Rec(mintree)" << std::endl;
    }

    testComponentTreeFZ(maxtree, "(1) Max-Tree after adjustment", imgMaxtree, numRows, numCols);
    testComponentTreeFZ(mintree, "(1) Min-Tree after pruning", imgMintree, numRows, numCols);

    NodeFZPtr Lmax_leaf1 = maxtree->getLeaves().front();
    std::cout << "Pruning id: " << Lmax_leaf1->getIndex() << std::endl;
    adjust.updateTree(mintree, Lmax_leaf1);
    maxtree->prunning(Lmax_leaf1);
    imgMaxtree = maxtree->reconstructionImage();
    imgMintree = mintree->reconstructionImage();

    if (isEquals(imgMaxtree, imgMintree, n)) {
        std::cout << "✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    } else {
        std::cout << "❌ Rec(maxtree) ≠ Rec(mintree)" << std::endl;
    }

    testComponentTreeFZ(mintree, "(1) Min-Tree after adjustment", imgMintree, numRows, numCols);
    testComponentTreeFZ(maxtree, "(1) Max-Tree after pruning", imgMaxtree, numRows, numCols);

    // Executar mais operações de pruning e ajuste
    NodeFZPtr node = maxtree->getRoot()->getChildren().front();
    maxtree->prunning(node);
    testComponentTreeFZ(maxtree, "(2) Max-Tree after big pruning", maxtree->reconstructionImage(), numRows, numCols);

    // Liberação de memória
    delete[] imgMaxtree;
    delete[] imgMintree;
    //delete maxtree;
    //delete mintree;
    delete[] img;

    return 0;
}