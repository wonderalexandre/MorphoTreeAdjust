#include "../include/ComponentTreeAdjustmentByFlatzone.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>




void ComponentTreeAdjustmentByFlatzone::updateTree(ComponentTreeFZPtr tree, FlatZone* flatzone) {

    bool isMaxtree = tree->isMaxtree();
    ComponentTreeFZPtr otherTree = getOtherTree(isMaxtree);

    int idNode = flatzone->front(); //pixel (id) of flatzone 
    NodeFZPtr node = tree->getSC(idNode);
    int newGrayLevel =  otherTree->getSC(idNode)->getParent()->getLevel(); //b = g(p)
    int oldGrayLevel = node->getLevel();  // a = f(p)
    
    NodeFZPtr nodeTauL = tree->getSC(idNode); //node of correspondence flatzone in other treee
    int pixelUpperBound = idNode; 

    bool nodeTauCNPsIsEqualL = nodeTauL->getNumFlatzone() == 1;
    FlatZone* flatzoneTauL = &tree->getFlatzoneByID(idNode); 
    if (PRINT_LOG) {
        outputLog << "Updating tree: " << (isMaxtree ? "Maxtree" : "Mintree") << std::endl;
        outputLog << "\tNode to update: id:" << node->getIndex() << ", level: " << node->getLevel() 
                  << ", |cnps|: " << node->getNumCNPs() << ", |children|: " << node->getChildren().size() 
                  << ", pixelUpperBound: " << pixelUpperBound 
                  << ", newGrayLevel: " << newGrayLevel 
                  << ", oldGrayLevel: " << oldGrayLevel 
                  << ", idFlatzone: " << idNode 
                  << ", nodeTauL id: " << nodeTauL->getIndex() 
                  << ", nodeTauL level: " << nodeTauL->getLevel() 
                  << ", nodeTauL numCNPs: " << nodeTauL->getNumCNPs() 
                  << std::endl;
        std::cout << outputLog.str() << std::endl;
    }
    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, flatzoneTauL->front(), pixelUpperBound, newGrayLevel, isMaxtree);

    int lambda = F.firstLambda(); //inicia com b = newGrayLevel
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



void ComponentTreeAdjustmentByFlatzone::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, FlatZone* flatzone) {
    updateTree(mintree, flatzone); 
    maxtree->mergeWithParent(flatzone); 
}

void ComponentTreeAdjustmentByFlatzone::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, FlatZone* flatzone) {
    updateTree(maxtree, flatzone);   
    mintree->mergeWithParent(flatzone); 
}
