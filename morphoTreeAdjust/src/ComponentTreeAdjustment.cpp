#include "../include/ComponentTreeAdjustment.hpp"
#include <unordered_set>
#include <list>
#include <iostream>


void ComponentTreeAdjustment::addNodesOfPath(NodeCT* nodeNL, NodeCT* nodeTauL, std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F) {
    std::cout << "=> adjustMinTree [addNodesOfPath - 1] \t ";
	for (NodeCT* n : nodeNL->getNodesOfPathToRoot()) {
		/*if(mapF[n->getLevel()]){
			mapF[n->getLevel()] = true;
			std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> set;
			F[n->getLevel()] = set;
		}*/
		F[n->getLevel()].insert(n);
        if (n == nodeTauL){ 
            break;
        }
    }
	std::cout << "=> adjustMinTree [addNodesOfPath - 2]\n";

}



void ComponentTreeAdjustment::adjustMinTree(ComponentTree &mintree, NodeCT *Lmax) {
	std::cout << "=> adjustMinTree\n";
	AdjacencyRelation* adj = mintree.getAdjacencyRelation();

    int newGrayLmax = Lmax->getParent()->getLevel(); //g(p)
    int grayLmax = Lmax->getLevel(); //f(p)

	std::list<int> cnpsL = Lmax->getCNPs();


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
	std::cout << "=> adjustMinTree [nodesNL - OK]:" << nodesNL.size() << "\n";
	
	
	std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>* F = new std::unordered_set<NodeCT*, NodeCT::NodeHashFunction>[255];
	
	std::cout << "=> adjustMinTree [alocacao F]\n";

	
    std::list<NodeCT*> subtrees;
    for (NodeCT *nodeNL : nodesNL) {
		std::cout << "=> adjustMinTree [for (NodeCT *nodeNL : nodesNL)] 1 \n";
	    if (newGrayLmax <= nodeNL->getLevel()) { //o nodeNL está entre g(p) e f(p)
			std::cout << "=> adjustMinTree [for (NodeCT *nodeNL : nodesNL)] 2 \n";
		    this->addNodesOfPath(nodeNL, nodeTauL, F);
	    } 
	    else { //o nodeNL está abaixo de g(p)
            // é armazenado somente a raiz da subtree antes de atingir o nivel g(p)
			std::cout << "=> adjustMinTree [for (NodeCT *nodeNL : nodesNL)] 3 \n";
            NodeCT* nodeSubtree = nodeNL;
            for (NodeCT *n : nodeNL->getNodesOfPathToRoot()) {
                if (n->getLevel() > newGrayLmax) {
                    break;
                }
                nodeSubtree = n;
            }
			std::cout << "=> adjustMinTree [for (NodeCT *nodeNL : nodesNL)] 3.1 \n";
            // se a subtree tiver level = g(p), então ela entra em F[\lambda]
            if (nodeSubtree->getLevel() == newGrayLmax) {
				std::cout << "=> adjustMinTree [for (NodeCT *nodeNL : nodesNL)] 3.2 \n";
				std::cout <<"nodeSubtree:" <<nodeSubtree->getIndex() << "\tnodeTauL:" << nodeTauL->getIndex()<< "\n";
				

                this->addNodesOfPath(nodeSubtree, nodeTauL, F);
            } 
            else {
				std::cout << "=> adjustMinTree [for (NodeCT *nodeNL : nodesNL)] 3.3 \n";
                subtrees.push_back(nodeSubtree);
            }
	    }

	}
    std::cout << "=> adjustMinTree [F[i] e subtrees - OK]\n";
	
    // merge
	int lambda = newGrayLmax; //g(p)
	NodeCT* nodeUnion = nullptr;
	NodeCT* nodeUnionPrevious = nullptr;
	std::cout << "lambda:" <<lambda << "\tF[lambda]:" << F[lambda].size() << "\n";
	while (lambda < grayLmax) {
	    if (!F[lambda].empty()) {
		    nodeUnion = *F[lambda].begin();
		    nodeUnion->getParent()->getChildren().remove(nodeUnion);
		    nodeUnion->setParent(nullptr);
		    for (NodeCT* n : F[lambda]) {
                if (n != nodeUnion) {
                    nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), n->getCNPs());
                    for (int p : n->getCNPs()) { //para manter atualizado o mapeamento SC: pixels -> nodes
                        mintree.setSC(p, nodeUnion);
                    }
                    nodeUnion->getChildren().splice(nodeUnion->getChildren().end(), n->getChildren());
                    for (NodeCT* son : n->getChildren()) {
                        son->setParent(nodeUnion);
                    }
                    n->getParent()->getChildren().remove(n);
                }
		    }
		    if (lambda == newGrayLmax) {
		        nodeUnion->getCNPs().splice(nodeUnion->getCNPs().end(), cnpsL);
		        for (NodeCT* n : subtrees) {
			        nodeUnion->addChild(n);
			        n->getParent()->getChildren().remove(n);
			        n->setParent(nodeUnion);
		        }
		    }
		    if (nodeUnionPrevious != nullptr) {
		        nodeUnionPrevious->setParent(nodeUnion);
		        nodeUnion->addChild(nodeUnionPrevious);
		    }
		    nodeUnionPrevious = nodeUnion;
	    }
	    lambda++;
	}
	std::cout << "=> adjustMinTree [merges - OK]\n";
    if (nodeTauCNPsIsIgualsL) {
	    nodeUnion->setParent(nodeTauL->getParent());
	    if (nodeTauL->getParent() != nullptr) {
		    nodeTauL->getParent()->addChild(nodeUnion);
		    nodeTauL->getParent()->getChildren().remove(nodeTauL);
	    } 
	    else { // novo root
		    mintree.setRoot(nodeUnion);
	    }
	    for (NodeCT* n : nodeTauL->getChildren()) {
		    if (n != nodeUnion && !nodeUnion->isChild(n)){
		        nodeUnion->getChildren().push_back(n);
		        n->setParent(nodeUnion);
		    }
	    }
		delete nodeTauL;

	} else {
	    nodeUnion->setParent(nodeTauL);
	    nodeTauL->addChild(nodeUnion);
	}


	delete[] F;
	F = nullptr;
}