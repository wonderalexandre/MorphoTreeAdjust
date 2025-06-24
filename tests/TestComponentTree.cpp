#include "../tests/Tests.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"

int main() {
    // Definição da imagem e parâmetros
    auto img = getPassatImage();
    double radioAdj = 1.5;

    // Criação das Component Trees
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::unique_ptr<FlatZonesGraph> graph = std::make_unique<FlatZonesGraph>(img, adj);
    std::unique_ptr<FlatZonesGraph> maxTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);
    std::unique_ptr<FlatZonesGraph> minTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, std::move(maxTreeFZGraph));
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, std::move(minTreeFZGraph));

    // Executar testes
    testComponentTreeFZ(mintree, "Min-Tree", mintree->reconstructionImage());
    testComponentTreeFZ(maxtree, "Max-Tree", maxtree->reconstructionImage());

    NodeFZPtr Lmin_leaf = mintree->getLeaves().front();
    std::cout << "---- Pruning leaf id: " << Lmin_leaf->getIndex() << " parent id: " << Lmin_leaf->getParent()->getIndex() << " ----" << std::endl;
    mintree->prunning(Lmin_leaf);
    auto imgMintree = mintree->reconstructionImage();
    testComponentTreeFZ(mintree, "Min-Tree after pruning", imgMintree);

    NodeFZPtr Lmax_leaf = maxtree->getLeaves().front();
    std::cout << "---- Pruning leaf id: " << Lmax_leaf->getIndex() << " parent id: " << Lmax_leaf->getParent()->getIndex() << " ----" << std::endl;
    maxtree->prunning(Lmax_leaf);
    auto imgMaxtree = maxtree->reconstructionImage();
    testComponentTreeFZ(maxtree, "Max-Tree after pruning", imgMaxtree);

    maxTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);
    minTreeFZGraph = std::make_unique<FlatZonesGraph>(*graph);

    maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, std::move(maxTreeFZGraph));
    mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, std::move(minTreeFZGraph));
    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);

    NodeFZPtr Lmin_leaf1 = mintree->getLeaves().front();
    std::cout << "Pruning id: " << Lmin_leaf1->getIndex() << std::endl;
    adjust.updateTree(maxtree, Lmin_leaf1);
    mintree->prunning(Lmin_leaf1);
    imgMaxtree = maxtree->reconstructionImage();
    imgMintree = mintree->reconstructionImage();

    if (imgMaxtree->isEqual(imgMintree)) {
        std::cout << "✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    } else {
        std::cout << "❌ Rec(maxtree) ≠ Rec(mintree)" << std::endl;
    }

    testComponentTreeFZ(maxtree, "(1) Max-Tree after adjustment", imgMaxtree);
    testComponentTreeFZ(mintree, "(1) Min-Tree after pruning", imgMintree);

    NodeFZPtr Lmax_leaf1 = maxtree->getLeaves().front();
    std::cout << "Pruning id: " << Lmax_leaf1->getIndex() << std::endl;
    adjust.updateTree(mintree, Lmax_leaf1);
    maxtree->prunning(Lmax_leaf1);
    imgMaxtree = maxtree->reconstructionImage();
    imgMintree = mintree->reconstructionImage();
    
    if (imgMaxtree->isEqual(imgMintree)) {
        std::cout << "✅ Rec(maxtree) = Rec(mintree)" << std::endl;
    } else {
        std::cout << "❌ Rec(maxtree) ≠ Rec(mintree)" << std::endl;
    }

    testComponentTreeFZ(mintree, "(1) Min-Tree after adjustment", imgMintree);
    testComponentTreeFZ(maxtree, "(1) Max-Tree after pruning", imgMaxtree);

    // Executar mais operações de pruning e ajuste
    NodeFZPtr node = maxtree->getRoot()->getChildren().front();
    maxtree->prunning(node);
    testComponentTreeFZ(maxtree, "(2) Max-Tree after big pruning", maxtree->reconstructionImage());

    return 0;
}