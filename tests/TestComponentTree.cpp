#include "../tests/Tests.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"

#include <cassert>
#include <vector>
#include <iostream>


int main(){
    auto img = getWonderImage();
    double radioAdj = 1.5;
    printImage(img);

    // Criação das Component Trees
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, false, adj);

    auto imgMaxtree = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());

    NodeId nodeId = maxtree->getLeaves().front();// maxtree->getSC(28);
    NodeFZ node = maxtree->proxy(nodeId);
    std::cout << "\nNode - ID: " << node.getIndex() << ", Level: " << node.getLevel() << ", Area: " << node.getArea() << "\n" << std::endl;
    maxtree->prunning(node);


    auto imgPrunned = maxtree->reconstructionImage();
    printTree(maxtree->getRoot());
    printImage(imgPrunned);
    testComponentTreeFZ(maxtree, "maxtreeFZ sem grafo", imgPrunned);
    return 0;
}

int test1() {
    // Definição da imagem e parâmetros
    auto img = getPassatImage();
    double radioAdj = 1.5;

    // Criação das Component Trees
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::unique_ptr<FlatZonesGraph> graph = std::make_unique<FlatZonesGraph>(img, adj);
    std::unique_ptr<FlatZonesGraph> maxTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);
    std::unique_ptr<FlatZonesGraph> minTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);

    ComponentTreeFZPtr maxtreePtr = std::make_shared<ComponentTreeFZ>(std::move(maxTreeFZGraph), true);
    ComponentTreeFZPtr mintreePtr = std::make_shared<ComponentTreeFZ>(std::move(minTreeFZGraph), false);
    ComponentTreeFZ* maxtree = maxtreePtr.get();
    ComponentTreeFZ* mintree = mintreePtr.get();

    // Executar testes
    testComponentTreeFZ(mintreePtr, "Min-Tree", mintree->reconstructionImage());
    testComponentTreeFZ(maxtreePtr, "Max-Tree", maxtree->reconstructionImage());


    maxTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);
    minTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);

    maxtreePtr = std::make_shared<ComponentTreeFZ>(std::move(maxTreeFZGraph), true);
    mintreePtr = std::make_shared<ComponentTreeFZ>(std::move(minTreeFZGraph), false);
    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);

    NodeId Lmin_leaf1Id = mintree->getLeaves().front();
    NodeFZ Lmin_leaf1 = mintree->proxy(Lmin_leaf1Id);
    std::cout << "Pruning id: " << Lmin_leaf1.getIndex() << std::endl;
    adjust.updateTree(maxtree, Lmin_leaf1);
    adjust.prunning(mintree, Lmin_leaf1);

    
    auto imgMaxtree = maxtree->reconstructionImage();
    auto imgMintree = mintree->reconstructionImage();

    if (imgMaxtree->isEqual(imgMintree)) {
        std::cout << "✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    } else {
        std::cout << "❌ Rec(maxtree) ≠ Rec(mintree)" << std::endl;
    }

    testComponentTreeFZ(maxtreePtr, "(1) Max-Tree after adjustment", imgMaxtree);
    testComponentTreeFZ(mintreePtr, "(1) Min-Tree after pruning", imgMintree);

    NodeId Lmax_leaf1Id = maxtree->getLeaves().front();
    NodeFZ Lmax_leaf1 = maxtree->proxy(Lmax_leaf1Id);
    std::cout << "Pruning id: " << Lmax_leaf1.getIndex() << std::endl;
    adjust.updateTree(mintree, Lmax_leaf1);
    adjust.prunning(maxtree, Lmax_leaf1);
    imgMaxtree = maxtree->reconstructionImage();
    imgMintree = mintree->reconstructionImage();
    
    if (imgMaxtree->isEqual(imgMintree)) {
        std::cout << "✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    } else {
        std::cout << "❌ Rec(maxtree) ≠ Rec(mintree)" << std::endl;
    }

    testComponentTreeFZ(mintreePtr, "(1) Min-Tree after adjustment", imgMintree);
    testComponentTreeFZ(maxtreePtr, "(1) Max-Tree after pruning", imgMaxtree);

    return 0;
}