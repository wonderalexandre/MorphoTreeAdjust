#include "../include/ComponentTreeAdjustment.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <utility>
#include "../../tests/Tests.hpp"

ComponentTreeAdjustment::ComponentTreeAdjustment(ComponentTree* maxtree, ComponentTree* mintree) 
    : maxtree(maxtree), mintree(mintree), 
      maxIndex(std::max(maxtree->getNumNodes(), mintree->getNumNodes())),
      visited(new bool[maxIndex]()), 
      F(maxIndex),
      unionNodes(maxtree->isMaxtree()) {  
    
    //std::cout << "maxIndex: " << maxIndex << std::endl; 
}


ComponentTreeAdjustment::~ComponentTreeAdjustment() {
    delete[] visited;
}


std::vector<NodeCT*> ComponentTreeAdjustment::getAdjacentNodes(ComponentTree* tree, NodeCT::IteratorCNPs flatZone) {
    bool isMaxtree = tree->isMaxtree();
    AdjacencyRelation* adj = tree->getAdjacencyRelation();

    std::vector<NodeCT*> nodesNL;
    for (int p : flatZone) {
        int grayFlatZone = tree->getSC(p)->getLevel();
        for (int q : adj->getAdjPixels(p)) {
            NodeCT* node = tree->getSC(q);
            int index = node->getIndex();
            if (!this->visited[index]) { 
                if ((!isMaxtree && node->getLevel() < grayFlatZone) || (isMaxtree && node->getLevel() > grayFlatZone)) {
                    nodesNL.push_back(node);  
                    this->visited[index] = true;
                }
            }
        }
    }
     
    return nodesNL;
}







void ComponentTreeAdjustment::buildMergedAndNestedCollections(ComponentTree* tree, NodeCT::IteratorCNPs flatZone, int newGrayLevel, bool isMaxtree){
	B_L.clear();
	F.resetCollection(isMaxtree);

    std::vector<NodeCT*> nodesNL = this->getAdjacentNodes(tree, flatZone);

    NodeCT* nodeTauL = tree->getSC(flatZone.front());

    for (NodeCT *nodeNL : nodesNL) {
        int index = nodeNL->getIndex();
        this->visited[index] = false; //mandendo o vetor visited inicializado com false para ser usado novamente no metodo getAdjacentNodes

	    if ((!isMaxtree && newGrayLevel <= nodeNL->getLevel()) || (isMaxtree && newGrayLevel >= nodeNL->getLevel())) { //o nodeNL está entre g(p) e f(p)
		    F.addNodesOfPath(nodeNL, nodeTauL);
	    } 
	    else { //o nodeNL está abaixo de g(p)
            // é armazenado somente a raiz da subtree antes de atingir o nivel g(p)
	        NodeCT* nodeSubtree = nodeNL;
            for (NodeCT *n : nodeNL->getNodesOfPathToRoot()) {
                if ((!isMaxtree && n->getLevel() > newGrayLevel) || (isMaxtree && n->getLevel() < newGrayLevel) ) {
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



void ComponentTreeAdjustment::updateTree2(ComponentTree* tree, NodeCT *rSubtree) {
   // std::cout << "\n==> start updateTree" << std::endl;
    if (rSubtree == nullptr) {
        std::cout << "tauRoot is nullptr" << std::endl;
        return;
    }

/*

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = rSubtree->getParent()->getLevel();  // g(p)

    
    std::list<int> flatzoneSubtree; 
    NodeCT* nodeTauStar = tree->getSC(rSubtree->getCNP(0));
    unionNodes.resetCollection(isMaxtree);
    
    for (NodeCT* n : rSubtree->getIteratorBreadthFirstTraversal()) {
        int n_cnp = n->getCNP(0);
        NodeCT* nodeTau = tree->getSC(n_cnp);
        
        unionNodes.addNode(nodeTau, n);
        
        std::list<int> flatzone = tree->getFlatzone(n_cnp); //TODO:
    
        if ((!isMaxtree && n->getLevel() > nodeTauStar->getLevel()) || (isMaxtree && n->getLevel() < nodeTauStar->getLevel())) {
            nodeTauStar = nodeTau;
            flatzoneSubtree.splice(flatzoneSubtree.begin(), flatzone);  
        } else {
            flatzoneSubtree.splice(flatzoneSubtree.end(), flatzone);
        }
    }
    

    int grayTauStar = nodeTauStar->getLevel();  // f(p)
    std::cout << "newGrayLevel: " << newGrayLevel  << std::endl;
    std::cout << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
    std::cout << "|cnpsSubtree| = " << flatzoneSubtree.size() << " = Area(rSubtree)= " << rSubtree->getArea() <<  std::endl;
    std::cout << "nodeTauStar: Id:" << nodeTauStar->getIndex() << "; level:" << nodeTauStar->getLevel() <<"; |cnps|:" << nodeTauStar->getNumCNPs() << std::endl;
    std::cout << "unionNodes: [";
    for (auto [nodeTau, isEqual] : unionNodes.getIterator()) {
        std::cout << nodeTau->getIndex() << ":" << nodeTau->getLevel() << ":" << isEqual << ", ";
    }
    std::cout << "]"<< std::endl;


    bool nodeTauCNPsIsEqualL = (nodeTauStar->getNumFlatzone() ==  1);
    
    
    this->buildMergedAndNestedCollections(tree,  NodeCT::IteratorCNPs{std::vector<std::list<int>>{std::move(flatzoneSubtree)}}, newGrayLevel, isMaxtree);
    
    
    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    NodeCT* nodeUnion = nullptr;
    NodeCT* nodeUnionPrevious = nullptr;
    NodeCT* nodeTauParentSubtree = nullptr;
    //std::cout << "==> Construiu as estruturas" << std::endl;

    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        
        std::vector<NodeCT*>& F_lambda = F.getMergedNodes(lambda);
        std::cout << "Lambda:" << lambda << "; size(F_lambda):" << F_lambda.size() << std::endl;

        nodeUnion = F_lambda.front();
        disconnect(nodeUnion);

        for (NodeCT* n : F_lambda) {
            if (n != nodeUnion) {
                nodeUnion->addCNPsOfAllFlatZone(std::move(n->getCNPsByFlatZone()), tree); //os mapeamentos sao atualizados
                
                for (NodeCT* son : n->getChildren()) {
                    son->setParent(nodeUnion);
                }
                nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());

                disconnect(n);
                n->getChildren().clear();
                //n->getCNPs().clear();
                delete n;
                n = nullptr;
                tree->setNumNodes( tree->getNumNodes() -1 );
            }
        }
        if (lambda == newGrayLevel) {
            nodeUnion->addCNPsToConnectedFlatzone(std::move(flatzoneSubtree), tree); //os mapeamentos sao atualizados
            
            for (NodeCT* n : B_L) {
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
        for (NodeCT* n : nodeUnion->getChildren()) {
            nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }
    //std::cout << "==> Fez os merges" << std::endl;

    if (nodeTauCNPsIsEqualL) {
       // std::cout << "==> nodeTauCNPsIsEqualL";

        NodeCT* parentNodeTauL = nodeTauStar->getParent();
        nodeUnion->setParent(parentNodeTauL);

        if (parentNodeTauL != nullptr) {
            parentNodeTauL->addChild(nodeUnion);
            for (NodeCT* n : nodeTauStar->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } else {  // Novo root
            NodeCT* newRoot = nodeUnion;
            if (!nodeTauStar->getChildren().empty()) {
                for (NodeCT* n : nodeTauStar->getChildren()) {
                    if ((isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeCT* n : nodeTauStar->getChildren()) {
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
    */
}



void ComponentTreeAdjustment::updateTree(ComponentTree* tree, NodeCT *L_leaf) {
   // std::cout << "\n==> start updateTree" << std::endl;
   assert(L_leaf != nullptr && "L_leaf is nullptr"); 
   


    bool isMaxtree = tree->isMaxtree();
    int newGrayL = L_leaf->getParent()->getLevel();  // g(p)
    int grayL = L_leaf->getLevel();  // f(p)
    NodeCT::IteratorCNPs cnpsL = L_leaf->getCNPs();
    
    NodeCT* nodeTauL = tree->getSC(cnpsL.front());
    bool nodeTauCNPsIsEqualL = nodeTauL->getNumFlatzone() == 1;
    std::list<int>& flatzoneTauL = tree->getFlatzoneRef(cnpsL.front()); 

    assert(L_leaf->getNumCNPs() == flatzoneTauL.size() && "O número de CNPs de L_leaf é diferente do numero de pixels da faltzone de tauL");
    
    
    this->buildMergedAndNestedCollections(tree, cnpsL, newGrayL, isMaxtree);
    
    
    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    NodeCT* nodeUnion = nullptr;
    NodeCT* nodeUnionPrevious = nullptr;
    
    //std::cout << "==> Construiu as estruturas" << std::endl;

    // Definição da direção do loop
    while ((isMaxtree && lambda > grayL) || (!isMaxtree && lambda < grayL)) {
        std::vector<NodeCT*>& F_lambda = F.getMergedNodes(lambda);

        nodeUnion = F_lambda.front();
        disconnect(nodeUnion);

        for (NodeCT* n : F_lambda) {
            if (n != nodeUnion) {
                nodeUnion->addCNPsOfDisjointFlatzones(n->moveCNPsByFlatZone(), tree);
                for (NodeCT* son : n->getChildren()) {
                    son->setParent(nodeUnion);
                }
                nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());
                disconnect(n);
                delete n;
                n = nullptr;
                tree->setNumNodes( tree->getNumNodes() -1 );
            }
        }
        if (lambda == newGrayL) {
            nodeUnion->addCNPsToConnectedFlatzone(std::move(flatzoneTauL), tree); //os mapeamentos sao atualizados
            nodeTauL->removeFlatzone(flatzoneTauL);
            for (NodeCT* n : B_L) {
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
        for (NodeCT* n : nodeUnion->getChildren()) {
            nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }
    //std::cout << "==> Fez os merges" << std::endl;

    if (nodeTauCNPsIsEqualL) {
       // std::cout << "==> nodeTauCNPsIsEqualL";

        NodeCT* parentNodeTauL = nodeTauL->getParent();
        nodeUnion->setParent(parentNodeTauL);

        if (parentNodeTauL != nullptr) {
            parentNodeTauL->addChild(nodeUnion);
            for (NodeCT* n : nodeTauL->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } else {  // Novo root
            NodeCT* newRoot = nodeUnion;
            if (!nodeTauL->getChildren().empty()) {
                for (NodeCT* n : nodeTauL->getChildren()) {
                    if ((isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeCT* n : nodeTauL->getChildren()) {
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
        tree->setNumNodes( tree->getNumNodes() -1 );  
        disconnect(nodeTauL);
        nodeTauL->getChildren().clear();
        delete nodeTauL;
        nodeTauL = nullptr;
       // std::cout << " OK"<< std::endl;
    } else {
       // std::cout << "==> NOT nodeTauCNPsIsEqualL";
        if (nodeUnion != nodeTauL) {
            nodeUnion->setParent(nodeTauL);
            nodeTauL->addChild(nodeUnion);
        }
        
    }
}


void ComponentTreeAdjustment::adjustMinTree(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning) {
	for(NodeCT* node: nodesToPruning) {
	    for (NodeCT* Lmax: node->getIteratorPostOrderTraversal()) { 
            assert(Lmax != mintree->getRoot() && "Lmax is root");
            assert(Lmax->isLeaf() && "Lmax não é uma folha");
            updateTree(mintree, Lmax); 
			maxtree->prunning(Lmax);

	    }
	}
	
}
    
void ComponentTreeAdjustment::adjustMaxTree(ComponentTree* maxtree, ComponentTree* mintree, std::vector<NodeCT*> nodesToPruning) {
	for(NodeCT* node: nodesToPruning) {  	
	    for (NodeCT* Lmin: node->getIteratorPostOrderTraversal()) {
			assert(Lmin != mintree->getRoot() && "Lmin is root");
            assert(Lmin->isLeaf() && "Lmin não é uma folha");

            updateTree(maxtree, Lmin);             
			mintree->prunning(Lmin);
            
	    }
	}
}


void ComponentTreeAdjustment::adjustMinTree2(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning) {
    
    for(NodeCT* rSubtree: nodesToPruning) {  	
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");
        updateTree2(mintree, rSubtree); 
        maxtree->prunning(rSubtree);
    }
	    
}

    
void ComponentTreeAdjustment::adjustMaxTree2(ComponentTree* maxtree, ComponentTree* mintree,std::vector<NodeCT*> nodesToPruning) {
    for(NodeCT* rSubtree: nodesToPruning) {  
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");

        updateTree(maxtree, rSubtree);   
        mintree->prunning(rSubtree);


    }
}