#include "../include/ComponentTreeAdjustment.hpp"
#include <unordered_set>
#include <list>
#include <iostream>
#include <functional>

// compare `<` std::less<int>() e  `>` std::greater<int>()
/*std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> ComponentTreeAdjustment::getAdjacentNodes(ComponentTree &tree, std::list<int> flatZone){
	bool isMaxtree = tree.isMaxtree();
	
	AdjacencyRelation* adj = tree.getAdjacencyRelation();
	std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> nodesNL;
	int grayFlatZone = tree.getSC(flatZone.front())->getLevel();
	for (int p : flatZone) {
		for (int q : adj->getAdjPixels(p)) {
			if ( (!isMaxtree && tree.getSC(q)->getLevel() < grayFlatZone) || ((isMaxtree && tree.getSC(q)->getLevel() > grayFlatZone)) ){ // p in Lmax and q notin Lmax
				nodesNL.insert( tree.getSC(q) );
			}
		}
	}
	return nodesNL;
}*/

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


std::vector<NodeCT*> ComponentTreeAdjustment::getAdjacentNodes(ComponentTree* tree, std::list<int> flatZone) {
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







void ComponentTreeAdjustment::buildMergedAndNestedCollections(ComponentTree* tree, std::list<int> flatZone, int newGrayLevel, bool isMaxtree){
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

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = rSubtree->getParent()->getLevel();  // g(p)
    
    
    std::list<int> cnpsSubtree;
    NodeCT* nodeTauStar = tree->getSC(rSubtree->getCNPs().front());
    unionNodes.resetCollection(isMaxtree);
    
    for(NodeCT* n: rSubtree->getIteratorBreadthFirstTraversal()){
        NodeCT* nodeTau = tree->getSC(n->getCNPs().front());
        unionNodes.addNode( nodeTau, n ); 
        std::list<int> n_cnps = n->getCNPsCopy();
        if( (!isMaxtree && n->getLevel() > nodeTauStar->getLevel()) || (isMaxtree && n->getLevel() < nodeTauStar->getLevel() ) ){
            nodeTauStar = tree->getSC(n->getCNPs().front());
            cnpsSubtree.splice(cnpsSubtree.begin(), n_cnps);
        }else{
            cnpsSubtree.splice(cnpsSubtree.end(), n_cnps);
        } 
        nodeTau->setArea(nodeTau->getArea() - n_cnps.size());   
        nodeTau->removeCNPs(n_cnps);
    }
    
    int grayTauStar = nodeTauStar->getLevel();  // f(p)
    std::cout << "newGrayLevel: " << newGrayLevel  << std::endl;
    std::cout << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
    std::cout << "|cnpsSubtree| = " << cnpsSubtree.size() << " = Area(rSubtree)= " << rSubtree->getArea() <<  std::endl;
    std::cout << "nodeTauStar: Id:" << nodeTauStar->getIndex() << "; level:" << nodeTauStar->getLevel() << std::endl;
    std::cout << "unionNodes: [";
    for (auto [nodeTau, isEqual] : unionNodes.getIterator()) {
        std::cout << nodeTau->getIndex() << ":" << nodeTau->getLevel() << ":" << isEqual << ", ";
    }
    std::cout << "]"<< std::endl;


    bool nodeTauCNPsIsEqualL = (cnpsSubtree.size() == nodeTauStar->getCNPs().size());
    
    
    this->buildMergedAndNestedCollections(tree, cnpsSubtree, newGrayLevel, isMaxtree);
    
    
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
                for (int p : n->getCNPs()) {  
                    if(tree->getSC(p)->getIndex() != nodeTauParentSubtree->getIndex())
                        tree->setSC(p, nodeUnion);
                }
                nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), n->getCNPs());
                for (NodeCT* son : n->getChildren()) {
                    son->setParent(nodeUnion);
                }
                  nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());

                disconnect(n);
                n->getChildren().clear();
                n->getCNPs().clear();
                delete n;
                n = nullptr;
                tree->setNumNodes( tree->getNumNodes() -1 );
            }
        }
        if (lambda == newGrayLevel) {
            for (int p : cnpsSubtree) {
                tree->setSC(p, nodeUnion);
            }
            nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), cnpsSubtree);

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
        nodeUnion->setArea(nodeUnion->getCNPs().size());
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
        nodeTauStar->getCNPs().clear();
        delete nodeTauStar;
        nodeTauStar = nullptr;
       // std::cout << " OK"<< std::endl;
    } else {
       // std::cout << "==> NOT nodeTauCNPsIsEqualL";
        if (nodeUnion != nodeTauStar) {
            nodeUnion->setParent(nodeTauStar);
            nodeTauStar->addChild(nodeUnion);
        }
        nodeTauStar->getCNPs().remove_if([&tree, &nodeTauStar](int p) {
            return tree->getSC(p) != nodeTauStar;
        });

    }
}



void ComponentTreeAdjustment::updateTree(ComponentTree* tree, NodeCT *L_leaf) {
   // std::cout << "\n==> start updateTree" << std::endl;
    if (L_leaf == nullptr) {
        std::cout << "L_leaf is nullptr" << std::endl;
        return;
    }

    bool isMaxtree = tree->isMaxtree();
    int newGrayL = L_leaf->getParent()->getLevel();  // g(p)
    int grayL = L_leaf->getLevel();  // f(p)
    std::list<int> cnpsL = L_leaf->getCNPsCopy();
    NodeCT* nodeTauL = tree->getSC(cnpsL.front());
    bool nodeTauCNPsIsEqualL = (cnpsL.size() == nodeTauL->getCNPs().size());
    
    
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
                for (int p : n->getCNPs()) {  // Atualiza mapeamento SC
                    tree->setSC(p, nodeUnion);
                }
                nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), n->getCNPs());
                   for (NodeCT* son : n->getChildren()) {
                        son->setParent(nodeUnion);
                }
                  nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());

                disconnect(n);
                n->getChildren().clear();
                n->getCNPs().clear();
                delete n;
                n = nullptr;
                tree->setNumNodes( tree->getNumNodes() -1 );
            }
        }
        if (lambda == newGrayL) {
            for (int p : cnpsL) {
                tree->setSC(p, nodeUnion);
            }
            nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), cnpsL);

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
        nodeUnion->setArea(nodeUnion->getCNPs().size());
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
        nodeTauL->getCNPs().clear();
        delete nodeTauL;
        nodeTauL = nullptr;
       // std::cout << " OK"<< std::endl;
    } else {
       // std::cout << "==> NOT nodeTauCNPsIsEqualL";
        if (nodeUnion != nodeTauL) {
            nodeUnion->setParent(nodeTauL);
            nodeTauL->addChild(nodeUnion);
        }
        nodeTauL->getCNPs().remove_if([&tree, &nodeTauL](int p) {
            return tree->getSC(p) != nodeTauL;
        });

    }
}


void ComponentTreeAdjustment::adjustMinTree(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning) {
	for(NodeCT* node: nodesToPruning) {
	    for (NodeCT* Lmax: node->getIteratorPostOrderTraversal()) { 
		 	if(Lmax == maxtree->getRoot()){
                std::cout << " Lmax is root"<< std::endl;
                exit(0);
            }
            updateTree(mintree, Lmax); 
			maxtree->prunning(Lmax);
	    }
	}
	
}
    
void ComponentTreeAdjustment::adjustMaxTree(ComponentTree* maxtree, ComponentTree* mintree, std::vector<NodeCT*> nodesToPruning) {
	for(NodeCT* node: nodesToPruning) {  	
	    for (NodeCT* Lmin: node->getIteratorPostOrderTraversal()) {
			if(Lmin == mintree->getRoot()){
                std::cout << " Lmin is root"<< std::endl;
                exit(0);
            }
            updateTree(maxtree, Lmin);   
			mintree->prunning(Lmin);
	    }
	}
}


void ComponentTreeAdjustment::adjustMinTree2(ComponentTree* mintree, ComponentTree* maxtree, std::vector<NodeCT*> nodesToPruning) {
    for(NodeCT* rSubtree: nodesToPruning) {  	
        updateTree2(mintree, rSubtree); 
        maxtree->prunning(rSubtree);
    }
	    
}
    
void ComponentTreeAdjustment::adjustMaxTree2(ComponentTree* maxtree, ComponentTree* mintree,std::vector<NodeCT*> nodesToPruning) {
    for(NodeCT* rSubtree: nodesToPruning) {  	
        updateTree(maxtree, rSubtree);   
        mintree->prunning(rSubtree);
    }
}