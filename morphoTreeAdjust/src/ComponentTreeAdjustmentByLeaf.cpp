#include "../include/ComponentTreeAdjustmentByLeaf.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>



void ComponentTreeAdjustmentByLeaf::updateTree(ComponentTreeFZPtr tree, NodeFZPtr leaf) {
    assert(leaf != nullptr && "L_leaf is nullptr"); 
    assert(leaf->isLeaf() && "L_leaf is not leaf"); 

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = leaf->getParent()->getLevel();  // b = g(p)
    int oldGrayLevel = leaf->getLevel();  // a = f(p)
    int idLeaf = leaf->getCNPsByFlatZone().begin()->second.front(); //pixel (id) of flatzone 
    
    NodeFZPtr nodeTauL = tree->getSC(idLeaf); //node of correspondence flatzone in other treee
    int pixelUpperBound = idLeaf; 

    bool nodeTauCNPsIsEqualL = nodeTauL->getNumFlatzone() == 1;
    FlatZone* flatzoneTauL = &tree->getFlatzoneByID(idLeaf); 
    //std::vector<FlatZonePtr> flatZonesTauL = {flatzoneTauL};
    
    assert(leaf->getNumCNPs() == flatzoneTauL->size() && "O número de CNPs de L_leaf é diferente do número de pixels da flatzone de tauL");

    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, flatzoneTauL->front(), pixelUpperBound, newGrayLevel, isMaxtree);

    int lambda = F.firstLambda(); //start with b = newGrayLevel
    NodeFZPtr nodeUnion = nullptr; // tau_{lambda}
    NodeFZPtr nodeUnionPrevious = nullptr; //maxtree: tau_{\lambda+1}, mintree:  tau_{\lambda-1}

    // Definição da direção do loop
    while ( (isMaxtree && lambda > oldGrayLevel) || (!isMaxtree && lambda < oldGrayLevel)) {
        std::vector<NodeFZPtr>& F_lambda = F.getMergedNodes(lambda);

        nodeUnion = F_lambda.front();
        disconnect(nodeUnion);

        for (NodeFZPtr n : F_lambda) {
            if (n != nodeUnion) {
                nodeUnion->addCNPsOfDisjointFlatzones(n->moveCNPsByFlatZone(), tree);
                mergedParentAndChildren(nodeUnion, n);
                disconnect(n, true);
                tree->setNumNodes(tree->getNumNodes() - 1);
            }
        }
        if (lambda == newGrayLevel) {
            int idFlatzoneTauL = flatzoneTauL->front();
            nodeUnion->addCNPsToConnectedFlatzone(std::move(*flatzoneTauL), tree); // Os mapeamentos são atualizados
            nodeTauL->removeFlatzone(idFlatzoneTauL);
            for (NodeFZPtr n : this->Fb) {
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
        for (NodeFZPtr n : nodeUnion->getChildren()) {
            nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
        }
        
        if(PRINT_LOG){
            outputLog << "\tnodeUnion = union(F_{" << lambda << "}) = " << 
                " id:" << nodeUnion->getIndex() << 
                ", level: " << nodeUnion->getLevel() << 
                ", |cnps|: " << nodeUnion->getNumCNPs() <<  
                ", |children|: " << nodeUnion->getChildren().size() << std::endl;
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }

    if (nodeTauCNPsIsEqualL) {
        NodeFZPtr parentNodeTauL = nodeTauL->getParent();
        nodeUnion->setParent(parentNodeTauL);

        if (parentNodeTauL != nullptr) {
            parentNodeTauL->addChild(nodeUnion);
            for (NodeFZPtr n : nodeTauL->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } 
        else {  // Novo root
            NodeFZPtr newRoot = nodeUnion;
            if (!nodeTauL->getChildren().empty()) {
                for (NodeFZPtr n : nodeTauL->getChildren()) {
                    if ( (isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeFZPtr n : nodeTauL->getChildren()) {
                    if (n != newRoot && !nodeUnion->isChild(n)) {
                        newRoot->addChild(n);
                        n->setParent(newRoot);
                    }
                }
            }
            
            newRoot->setArea(nodeTauL->getArea());
            newRoot->setParent(nullptr);
            tree->setRoot(newRoot);
            
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


void ComponentTreeAdjustmentByLeaf::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToPruning) {
    for (NodeFZPtr node : nodesToPruning) {        
        
        for (NodeFZPtr Lmax : node->getIteratorPostOrderTraversal()) { 
            assert(Lmax != maxtree->getRoot() && "Lmax is root");
            assert(Lmax->isLeaf() && "Lmax não é uma folha");
            
            updateTree(mintree, Lmax);     
            maxtree->prunning(Lmax);
        }
    }
}

void ComponentTreeAdjustmentByLeaf::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToPruning) {
    for (NodeFZPtr node : nodesToPruning) {  	
        for (NodeFZPtr Lmin : node->getIteratorPostOrderTraversal()) {
            assert(Lmin != mintree->getRoot() && "Lmin is root");
            assert(Lmin->isLeaf() && "Lmin não é uma folha");
            
            updateTree(maxtree, Lmin);             
            mintree->prunning(Lmin);
        }
    }
}
