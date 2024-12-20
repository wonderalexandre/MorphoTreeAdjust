#include "../include/ComponentTreeAdjustment.hpp"
#include <unordered_set>
#include <list>
#include <iostream>


void ComponentTreeAdjustment::adjustMinTree(ComponentTree &mintree, NodeCT *Lmax) {
	if(Lmax == nullptr){
		std::cout << "Lmax is nullptr";
		return;
	}

	AdjacencyRelation* adj = mintree.getAdjacencyRelation();

    int newGrayLmax = Lmax->getParent()->getLevel(); //g(p)
    int grayLmax = Lmax->getLevel(); //f(p)
	std::list<int> cnpsL = Lmax->getCNPsCopy();
	NodeCT* nodeTauL = mintree.getSC( cnpsL.front() );
	bool nodeTauCNPsIsIgualsL = (cnpsL.size() == nodeTauL->getCNPs().size());

    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> nodesNL;
	for (int p : cnpsL) {
	    for (int q : adj->getAdjPixels(p)) {
		    if (mintree.getSC(q)->getLevel() < grayLmax) { // p in Lmax and q notin Lmax
				nodesNL.insert( mintree.getSC(q) );
		    }
	    }
	}
	
	std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> B_L;
    for (NodeCT *nodeNL : nodesNL) {
	    if (newGrayLmax <= nodeNL->getLevel()) { //o nodeNL está entre g(p) e f(p)
		    this->addNodesOfPath(nodeNL, nodeTauL);
	    } 
	    else { //o nodeNL está abaixo de g(p)
            // é armazenado somente a raiz da subtree antes de atingir o nivel g(p)
	        NodeCT* nodeSubtree = nodeNL;
            for (NodeCT *n : nodeNL->getNodesOfPathToRoot()) {
                if (n->getLevel() > newGrayLmax) {
                    break;
                }
                nodeSubtree = n;
            }
	        // se a subtree tiver level = g(p), então ela entra em F[\lambda]
            if (nodeSubtree->getLevel() == newGrayLmax) {
                this->addNodesOfPath(nodeSubtree, nodeTauL);
            } 
            else {
                B_L.insert(nodeSubtree);
            }
	    }
	}

    // merge
	int lambda = newGrayLmax; //g(p)
	NodeCT* nodeUnion = nullptr;
	NodeCT* nodeUnionPrevious = nullptr;
	while (lambda < grayLmax) {
	    if (existsInF(lambda)) {
			std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F_lambda = getF(lambda);
		    nodeUnion = *F_lambda->begin();
			disconnect(nodeUnion);

		    for (NodeCT* n : *F_lambda) {
                if (n != nodeUnion) {
                    for (int p : n->getCNPs()) { //para manter atualizado o mapeamento SC: pixels -> nodes
                        mintree.setSC(p, nodeUnion);
                    }
					nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), n->getCNPs());

                    for (NodeCT* son : n->getChildren()) {
                        son->setParent(nodeUnion);
                    }
					nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());

					disconnect(n);
					delete n; n = nullptr;
                }
		    }
		    if (lambda == newGrayLmax) {
				for (int p : cnpsL) { //para manter atualizado o mapeamento SC: pixels -> nodes
                    mintree.setSC(p, nodeUnion);
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
		    
			//update: area attribute
			nodeUnion->setArea( nodeUnion->getCNPs().size() );
			for(NodeCT* n: nodeUnion->getChildren()){
				nodeUnion->setArea( nodeUnion->getArea() + n->getArea());
			}

			nodeUnionPrevious = nodeUnion;
	    }
	    lambda++;
	}
    if (nodeTauCNPsIsIgualsL) {
		NodeCT* parentNodeTauL = nodeTauL->getParent();
		nodeUnion->setParent(parentNodeTauL);

		if (parentNodeTauL != nullptr) {
			parentNodeTauL->addChild(nodeUnion);

			for (NodeCT* n : nodeTauL->getChildren()) {
		    	if (n != nodeUnion && !nodeUnion->isChild(n)){
		        	nodeUnion->getChildren().push_back(n);
					nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
		       		n->setParent(nodeUnion);
		    	}
	    	}
		}else{ 
			NodeCT* newRoot = nodeUnion;
			if(nodeTauL->getChildren().size() != 0){
				for (NodeCT* n : nodeTauL->getChildren()) {
					if(n->getLevel() > newRoot->getLevel())
						newRoot = n;
				}
				if(newRoot != nodeUnion){
					newRoot->addChild(nodeUnion);
					nodeUnion->setParent(newRoot);
				}
				for (NodeCT* n : nodeTauL->getChildren()) {
					if (n != newRoot && !nodeUnion->isChild(n)){
						newRoot->addChild(n);
						n->setParent(newRoot);
					}
				}
			}
			// novo root
			newRoot->setArea(nodeTauL->getArea());
			newRoot->setParent(nullptr);
			mintree.setRoot(nodeUnion);
		}
		
		disconnect(nodeTauL);
		delete nodeTauL; nodeTauL = nullptr;

	} 
	else {
		if(nodeUnion != nodeTauL){
	    	nodeUnion->setParent(nodeTauL);
	    	nodeTauL->addChild(nodeUnion);
		}

		std::list<int> newCNPs;
  		for (int p : nodeTauL->getCNPs()) {
			if(mintree.getSC(p) == nodeTauL)
				newCNPs.push_back(p);
			}
	    nodeTauL->setCNPs(newCNPs);
	}
	clearCollectionF();
}




void ComponentTreeAdjustment::adjustMaxTree(ComponentTree &maxtree, NodeCT *Lmin) {
	if(Lmin == nullptr){
		std::cout << "Lmin is nullptr";
		return;
	}

	AdjacencyRelation* adj = maxtree.getAdjacencyRelation();

    int newGrayLmax = Lmin->getParent()->getLevel(); //g(p)
    int grayLmax = Lmin->getLevel(); //f(p)

	std::list<int> cnpsL = Lmin->getCNPsCopy();

	NodeCT* nodeTauL = maxtree.getSC( cnpsL.front() );
	
	bool nodeTauCNPsIsIgualsL = (cnpsL.size() == nodeTauL->getCNPs().size());

    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> nodesNL;
	for (int p : cnpsL) {
	    for (int q : adj->getAdjPixels(p)) {
		    if (maxtree.getSC(q)->getLevel() > grayLmax) { // p in Lmin and q notin Lmin
				nodesNL.insert( maxtree.getSC(q) );
		    }
	    }
	}
	
	std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> B_L;
    for (NodeCT *nodeNL : nodesNL) {
	    if (newGrayLmax >= nodeNL->getLevel()) { //o nodeNL está entre g(p) e f(p)
		    this->addNodesOfPath(nodeNL, nodeTauL);
	    } 
	    else { //o nodeNL está acima de g(p)
            // é armazenado somente a raiz da subtree antes de atingir o nivel g(p)
	        NodeCT* nodeSubtree = nodeNL;
            for (NodeCT *n : nodeNL->getNodesOfPathToRoot()) {
                if (n->getLevel() < newGrayLmax) {
                    break;
                }
                nodeSubtree = n;
            }
	        // se a subtree tiver level = g(p), então ela entra em F[\lambda]
            if (nodeSubtree->getLevel() == newGrayLmax) {
                this->addNodesOfPath(nodeSubtree, nodeTauL);
            } 
            else {
                B_L.insert(nodeSubtree);
            }
	    }
	}

    // merge
	int lambda = newGrayLmax; //g(p)
	NodeCT* nodeUnion = nullptr;
	NodeCT* nodeUnionPrevious = nullptr;
	while (lambda > grayLmax) {
	    if (existsInF(lambda)) {
			std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F_lambda = getF(lambda);
		    nodeUnion = *F_lambda->begin();
			disconnect(nodeUnion);

		    for (NodeCT* n : *F_lambda) {
                if (n != nodeUnion) {
                    for (int p : n->getCNPs()) { //para manter atualizado o mapeamento SC: pixels -> nodes
                        maxtree.setSC(p, nodeUnion);
                    }
					nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), n->getCNPs());

                    for (NodeCT* son : n->getChildren()) {
                        son->setParent(nodeUnion);
                    }
					nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());

					disconnect(n);
					delete n; n = nullptr;
                }
		    }
		    if (lambda == newGrayLmax) {
				for (int p : cnpsL) { //para manter atualizado o mapeamento SC: pixels -> nodes
                    maxtree.setSC(p, nodeUnion);
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
		    
			//update: area attribute
			nodeUnion->setArea( nodeUnion->getCNPs().size() );
			for(NodeCT* n: nodeUnion->getChildren()){
				nodeUnion->setArea( nodeUnion->getArea() + n->getArea());
			}

			nodeUnionPrevious = nodeUnion;
	    }
	    lambda--;
	}
    if (nodeTauCNPsIsIgualsL) {
		NodeCT* parentNodeTauL = nodeTauL->getParent();
		nodeUnion->setParent(parentNodeTauL);

		if (parentNodeTauL != nullptr) {
			parentNodeTauL->addChild(nodeUnion);

			for (NodeCT* n : nodeTauL->getChildren()) {
		    	if (n != nodeUnion && !nodeUnion->isChild(n)){
		        	nodeUnion->getChildren().push_back(n);
					nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
		       		n->setParent(nodeUnion);
		    	}
	    	}
		}else{ 
			NodeCT* newRoot = nodeUnion;
			if(nodeTauL->getChildren().size() != 0){
				for (NodeCT* n : nodeTauL->getChildren()) {
					if(n->getLevel() < newRoot->getLevel())
						newRoot = n;
				}
				if(newRoot != nodeUnion){
					newRoot->addChild(nodeUnion);
					nodeUnion->setParent(newRoot);
				}
				for (NodeCT* n : nodeTauL->getChildren()) {
					if (n != newRoot && !nodeUnion->isChild(n)){
						newRoot->addChild(n);
						n->setParent(newRoot);
					}
				}
			}
			// novo root
			newRoot->setArea(nodeTauL->getArea());
			newRoot->setParent(nullptr);
			maxtree.setRoot(nodeUnion);
		}
		
		disconnect(nodeTauL);
		delete nodeTauL; nodeTauL = nullptr;

	} 
	else {
		if(nodeUnion != nodeTauL){
	    	nodeUnion->setParent(nodeTauL);
	    	nodeTauL->addChild(nodeUnion);
		}
		
		std::list<int> newCNPs;
  		for (int p : nodeTauL->getCNPs()) {
			if(maxtree.getSC(p) == nodeTauL)
				newCNPs.push_back(p);
			}
	    nodeTauL->setCNPs(newCNPs);
	}
	clearCollectionF();
}



void ComponentTreeAdjustment::adjustMinTree(ComponentTree &mintree, ComponentTree &maxtree, std::list<NodeCT*> nodesToPruning) {
	for(NodeCT* node: nodesToPruning) {
	    std::list<NodeCT*> descendants = maxtree.getDescendantsInPostOrder(node);
	    for (NodeCT* Lmax: descendants) { 
			
		 	adjustMinTree(mintree, Lmax);   
			maxtree.prunning(Lmax);
	    }
	}
	
}
    
void ComponentTreeAdjustment::adjustMaxTree(ComponentTree &maxtree, ComponentTree &mintree, std::list<NodeCT*> nodesToPruning) {
	for(NodeCT* node: nodesToPruning) {
	    std::list<NodeCT*> descendants = mintree.getDescendantsInPostOrder(node);
	    	
	    for (NodeCT* Lmin: descendants) {
			adjustMaxTree(maxtree, Lmin);   
			mintree.prunning(Lmin);
	    }
	}
}