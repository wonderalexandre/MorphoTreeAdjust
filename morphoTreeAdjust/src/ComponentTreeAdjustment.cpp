#include "../include/ComponentTreeAdjustment.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <utility>
#include "../../tests/Tests.hpp"


ComponentTreeAdjustment::ComponentTreeAdjustment(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree) 
    : maxtree(maxtree), mintree(mintree), 
      maxIndex(std::max(maxtree->getNumNodes(), mintree->getNumNodes())),
      F(maxIndex),
      unionNodeTauSubtree(maxtree->isMaxtree()) 
{    }

ComponentTreeAdjustment::~ComponentTreeAdjustment() { }



void ComponentTreeAdjustment::buildMergedAndNestedCollections(ComponentTreeFZ* tree, std::vector<FlatZoneRef>& flatZone, int newGrayLevel, bool isMaxtree){
	B_L.clear();
	F.resetCollection(isMaxtree);
    if(this->pixelUpperBound == -1){
        this->pixelUpperBound = flatZone.front().get().front();
    }
    F.computerAdjacentNodes(tree, flatZone);
    NodeFZ* nodeTauL = tree->getSC(this->pixelUpperBound); //pixel de tauStar ou tauL, para termos o node (limite) mais proximo de root

    for (NodeFZ* nodeNL: F.getAdjacentNodes()) {
	    if( (isMaxtree && nodeNL->getLevel() <= newGrayLevel) || (!isMaxtree &&  nodeNL->getLevel() >= newGrayLevel)) { //o nodeNL está entre g(p) e f(p)
		    F.addNodesOfPath(nodeNL, nodeTauL); 
	    } 
	    else { 
            //o nodeNL está abaixo de g(p). É armazenado somente a raiz da subtree antes de atingir o nivel g(p)
	        NodeFZ* nodeSubtree = nodeNL;
            for (NodeFZ *n : nodeNL->getNodesOfPathToRoot()) {
                if ( (isMaxtree && newGrayLevel > n->getLevel()) || (!isMaxtree && newGrayLevel < n->getLevel())) {
                    break;
                }
                nodeSubtree = n; 
            }
	        // se a subtree tiver level = g(p), então ela entra em F[\lambda]
            if (nodeSubtree->getLevel() == newGrayLevel) {
                F.addNodesOfPath(nodeSubtree, nodeTauL);
            } 
            else {
                B_L.push_back(nodeSubtree);
            }
	    }
        
	}
}


void ComponentTreeAdjustment::updateTree2(ComponentTreeFZ* tree, NodeFZ* rSubtree) {
    assert(rSubtree != nullptr && "rSubtree is nullptr"); 
    ComponentTreeFZ* otherTree = tree->isMaxtree()? this->mintree : this->maxtree;

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = rSubtree->getParent()->getLevel();  // g(p)
    
    NodeFZ* nodeTauStar = tree->getSC(rSubtree->getRepresentativeCNPs());
    unionNodeTauSubtree.resetCollection(isMaxtree);
    for (NodeFZ* n : rSubtree->getIteratorBreadthFirstTraversal()) {
        for(auto& [idFlatZoneSubtree, fzSubtree]: n->getCNPsByFlatZone()){
            //int idFlatZoneSubtree = fzSubtree.front();
            NodeFZ* nodeTau = tree->getSC(idFlatZoneSubtree);
            FlatZone& fzTau = tree->getFlatzoneByID(idFlatZoneSubtree);
            unionNodeTauSubtree.addNode(nodeTau, fzTau); //, fzSubtree.size() == fzTau.size()           
            if ((!isMaxtree && n->getLevel() > nodeTauStar->getLevel()) || (isMaxtree && n->getLevel() < nodeTauStar->getLevel())) {
                nodeTauStar = nodeTau;
            }
        }
    }
    unionNodeTauSubtree.setNodeStar(nodeTauStar);
    this->pixelUpperBound = nodeTauStar->getRepresentativeCNPs(); 
    int grayTauStar = nodeTauStar->getLevel();  // f(pixelUpperBound)
    std::cout << "newGrayLevel: " << newGrayLevel  << std::endl; //g(p)
    std::cout << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
    std::cout << "Area(rSubtree)= " << rSubtree->getArea() <<  std::endl;
    std::cout << "nodeTauStar: Id:" << nodeTauStar->getIndex() << "; level:" << nodeTauStar->getLevel() <<"; |cnps|:" << nodeTauStar->getNumCNPs() << std::endl;
    std::cout << "unionNodes: [";
    for(FlatZoneNode fzNode: unionNodeTauSubtree.getFlatzoneNodeList()){
        NodeFZ* nodeTau = fzNode.node;
        FlatZone* fzTau = fzNode.flatzone;
        std::cout << "(id:" << nodeTau->getIndex() << ", level:" << nodeTau->getLevel() << ", |cnps|:"<< nodeTau->getNumCNPs() << ", idFZ:" << fzTau->front() <<", |fz|:" << fzTau->size() << "), ";
    }
    std::cout << "]"<< std::endl;

    bool nodeTauCNPsIsEqualL = (nodeTauStar->getNumFlatzone() ==  1);
    std::vector<FlatZoneRef>  flatzonesTauSubtree = unionNodeTauSubtree.getFlatzones();
    this->buildMergedAndNestedCollections(tree,  flatzonesTauSubtree, newGrayLevel, isMaxtree);
    
    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    NodeFZ* nodeUnion = nullptr;
    NodeFZ* nodeUnionPrevious = nullptr;
    NodeFZ* nodeTauParentSubtree = nullptr;


    
    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        std::vector<NodeFZ*>& F_lambda = F.getMergedNodes(lambda);
        std::cout << "Lambda:" << lambda << "; size(F_lambda):" << F_lambda.size() << std::endl;

        nodeUnion = F_lambda.front();
        disconnect(nodeUnion);

        for (NodeFZ* n : F_lambda) {
            if (n != nodeUnion) {
                nodeUnion->addCNPsOfDisjointFlatzones(n->moveCNPsByFlatZone(), tree);
                mergedParentAndChildren(nodeUnion, n);
                disconnect(n, true);
                tree->setNumNodes(tree->getNumNodes() - 1);
            }
        }
        if (lambda == newGrayLevel) {
            unionNodeTauSubtree.addCNPsToConnectedFlatzone(nodeUnion, tree); // Os mapeamentos são atualizados
            unionNodeTauSubtree.removeFlatzones(tree);
            for (NodeFZ* n : B_L) {
                disconnect(n);
                nodeUnion->addChild(n);
                n->setParent(nodeUnion);
            }
            nodeTauParentSubtree = nodeUnion;

        }
        if (nodeUnionPrevious != nullptr) {
            nodeUnionPrevious->setParent(nodeUnion);
            nodeUnion->addChild(nodeUnionPrevious);
        }

        // Atualiza atributo de área
        nodeUnion->setArea(nodeUnion->getNumCNPs());
        for (NodeFZ* n : nodeUnion->getChildren()) {
            nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }
    //std::cout << "==> Fez os merges" << std::endl;

    if (nodeTauCNPsIsEqualL) {
       // std::cout << "==> nodeTauCNPsIsEqualL";

        NodeFZ* parentNodeTauL = nodeTauStar->getParent();
        nodeUnion->setParent(parentNodeTauL);

        if (parentNodeTauL != nullptr) {
            parentNodeTauL->addChild(nodeUnion);
            for (NodeFZ* n : nodeTauStar->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } else {  // Novo root
            NodeFZ* newRoot = nodeUnion;
            if (!nodeTauStar->getChildren().empty()) {
                for (NodeFZ* n : nodeTauStar->getChildren()) {
                    if ((isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeFZ* n : nodeTauStar->getChildren()) {
                    if (n != newRoot && !nodeUnion->isChild(n)) {
                        newRoot->addChild(n);
                        n->setParent(newRoot);
                    }
                }
            }
            
            newRoot->setArea(nodeTauStar->getArea());
            newRoot->setParent(nullptr);
            tree->setRoot(nodeUnion);
            
        }
        tree->setNumNodes( tree->getNumNodes() -1 );  
        disconnect(nodeTauStar);
        nodeTauStar->getChildren().clear();
        //nodeTauStar->getCNPs().clear();
        delete nodeTauStar;
        nodeTauStar = nullptr;
       // std::cout << " OK"<< std::endl;
    } else {
       // std::cout << "==> NOT nodeTauCNPsIsEqualL";
        if (nodeUnion != nodeTauStar) {
            nodeUnion->setParent(nodeTauStar);
            nodeTauStar->addChild(nodeUnion);
        }
        //nodeTauStar->removeCNPs(cnpsSubtree);
        //nodeTauStar->getCNPs().remove_if([&tree, &nodeTauStar](int p) {
        //   return tree->getSC(p) != nodeTauStar;
        //});

    }
       
    
}


void ComponentTreeAdjustment::updateTree(ComponentTreeFZ* tree, NodeFZ* L_leaf) {
    assert(L_leaf != nullptr && "L_leaf is nullptr"); 

    bool isMaxtree = tree->isMaxtree();
    int newGrayL = L_leaf->getParent()->getLevel();  // g(p)
    int grayL = L_leaf->getLevel();  // f(p)
    int idLeaf = L_leaf->getRepresentativeCNPs(); //sao CNPs  pois L_leaf é uma folha
    
    NodeFZ* nodeTauL = tree->getSC(idLeaf);
    this->pixelUpperBound = idLeaf;

    bool nodeTauCNPsIsEqualL = nodeTauL->getNumFlatzone() == 1;

    FlatZone& flatzoneTauL = tree->getFlatzoneByID(idLeaf); 
    std::vector<FlatZoneRef> flatZonesTauL = {flatzoneTauL};
    
    assert(L_leaf->getNumCNPs() == flatzoneTauL.size() && "O número de CNPs de L_leaf é diferente do número de pixels da flatzone de tauL");

    this->buildMergedAndNestedCollections(tree, flatZonesTauL, newGrayL, isMaxtree);

    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    NodeFZ* nodeUnion = nullptr;
    NodeFZ* nodeUnionPrevious = nullptr;

    // Definição da direção do loop
    while ( (isMaxtree && lambda > grayL) || (!isMaxtree && lambda < grayL)) {
        std::vector<NodeFZ*>& F_lambda = F.getMergedNodes(lambda);

        nodeUnion = F_lambda.front();
        disconnect(nodeUnion);

        for (NodeFZ* n : F_lambda) {
            if (n != nodeUnion) {
                nodeUnion->addCNPsOfDisjointFlatzones(n->moveCNPsByFlatZone(), tree);
                mergedParentAndChildren(nodeUnion, n);
                disconnect(n, true);
                tree->setNumNodes(tree->getNumNodes() - 1);
            }
        }
        if (lambda == newGrayL) {
            int idFlatzoneTauL = flatzoneTauL.front();
            nodeUnion->addCNPsToConnectedFlatzone(std::move(flatzoneTauL), tree); // Os mapeamentos são atualizados
            nodeTauL->removeFlatzone(idFlatzoneTauL);
            for (NodeFZ* n : B_L) {
                disconnect(n);
                nodeUnion->addChild(n);
                n->setParent(nodeUnion);
            }
        }
        if (nodeUnionPrevious != nullptr) {
            nodeUnionPrevious->setParent(nodeUnion);
            nodeUnion->addChild(nodeUnionPrevious);
        }

        // Atualiza atributo de área
        nodeUnion->setArea(nodeUnion->getNumCNPs());
        for (NodeFZ* n : nodeUnion->getChildren()) {
            nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }

    if (nodeTauCNPsIsEqualL) {
        NodeFZ* parentNodeTauL = nodeTauL->getParent();
        nodeUnion->setParent(parentNodeTauL);

        if (parentNodeTauL != nullptr) {
            parentNodeTauL->addChild(nodeUnion);
            for (NodeFZ* n : nodeTauL->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } 
        else {  // Novo root
            NodeFZ* newRoot = nodeUnion;
            if (!nodeTauL->getChildren().empty()) {
                for (NodeFZ* n : nodeTauL->getChildren()) {
                    if ((isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeFZ* n : nodeTauL->getChildren()) {
                    if (n != newRoot && !nodeUnion->isChild(n)) {
                        newRoot->addChild(n);
                        n->setParent(newRoot);
                    }
                }
            }
            
            newRoot->setArea(nodeTauL->getArea());
            newRoot->setParent(nullptr);
            tree->setRoot(nodeUnion);
        }
        tree->setNumNodes(tree->getNumNodes() - 1);
        disconnect(nodeTauL, true);
    } else {
        if (nodeUnion != nodeTauL) {
            nodeUnion->setParent(nodeTauL);
            nodeTauL->addChild(nodeUnion);
        }
    }
}


void ComponentTreeAdjustment::adjustMinTree(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* node : nodesToPruning) {
        for (NodeFZ* Lmax : node->getIteratorPostOrderTraversal()) { 
            assert(Lmax != mintree->getRoot() && "Lmax is root");
            assert(Lmax->isLeaf() && "Lmax não é uma folha");
            updateTree(mintree, Lmax); 
            maxtree->prunning(Lmax);
        }
    }
}

void ComponentTreeAdjustment::adjustMaxTree(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* node : nodesToPruning) {  	
        for (NodeFZ* Lmin : node->getIteratorPostOrderTraversal()) {
            assert(Lmin != mintree->getRoot() && "Lmin is root");
            assert(Lmin->isLeaf() && "Lmin não é uma folha");

            updateTree(maxtree, Lmin);             
            mintree->prunning(Lmin);
        }
    }
}

void ComponentTreeAdjustment::adjustMinTree2(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* rSubtree : nodesToPruning) {  	
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");
        updateTree2(mintree, rSubtree); 
        maxtree->prunning(rSubtree);
    }
}

void ComponentTreeAdjustment::adjustMaxTree2(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* rSubtree : nodesToPruning) {  
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");

        updateTree(maxtree, rSubtree);   
        mintree->prunning(rSubtree);
    }
}