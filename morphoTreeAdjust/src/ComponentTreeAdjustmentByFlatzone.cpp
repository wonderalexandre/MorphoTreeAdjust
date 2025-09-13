#include "../include/ComponentTreeAdjustmentByFlatzone.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>



void ComponentTreeAdjustmentByFlatzone::updateTree(ComponentTreeFZPtr tree, int repFZ) {
    bool isMaxtree = tree->isMaxtree();
    ComponentTreeFZPtr otherTree = getOtherTree(isMaxtree);
    NodeFZ nodeTauL = tree->getSC(repFZ); //node of correspondence flatzone in other treee
    this->areaFZsRemoved = tree->getFlatZonesGraph()->getNumPixelInFlatzone(repFZ);

    int newGrayLevel =  otherTree->getSC(repFZ).getParent().getLevel(); //b = g(p)
    int oldGrayLevel = nodeTauL.getLevel();  // a = f(p)
    
    bool nodeTauCNPsIsEqualL = nodeTauL.getNumFlatzone() == 1;
    if (PRINT_LOG) {
        outputLog << "Updating tree: " << (isMaxtree ? "Maxtree" : "Mintree") << std::endl;
        outputLog << "\tNode to update: id:" << nodeTauL.getIndex() << ", level: " << nodeTauL.getLevel() 
                  << ", |cnps|: " << nodeTauL.getNumCNPs() << ", |children|: " << nodeTauL.getNumChildren() 
                  << ", newGrayLevel: " << newGrayLevel 
                  << ", oldGrayLevel: " << oldGrayLevel 
                  << ", repFZ: " << repFZ 
                  << std::endl;
        //std::cout << outputLog.str() << std::endl;
    }
    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, repFZ, newGrayLevel, isMaxtree);
    int lambda = F.firstLambda(); //inicia com b = newGrayLevel

    
    NodeFZ newNode;
    if(lambda != newGrayLevel){ //Não tem nenhum nó em F_{\lambda} no nível b, então precisamos criar um novo nó
        
        int threshold1 = (isMaxtree) ? nodeTauL.getLevel()+1 : nodeTauL.getLevel()-1; 
        newNode = tree->createNode(repFZ, nodeTauL, threshold1, newGrayLevel);
        newNode.addRepCNPs(repFZ);
        newNode.setArea(areaFZsRemoved);

        nodeTauL.addChild(newNode);
        nodeTauL.removeFlatzone(repFZ); 

        //anexa Fb no newNode
        for (NodeId nodeId : F.getFb()) {
            newNode.setArea(newNode.getArea() + tree->getAreaById(nodeId)); //atualiza a area do newNode            
            disconnect(tree, nodeId, false);
            tree->addChildById(newNode, nodeId);
        }
    }
    

    NodeFZ nodeUnion = newNode; // tau_{lambda}
    NodeFZ nodeUnionPrevious; //maxtree: tau_{\lambda+1}, mintree:  tau_{\lambda-1}

    // Definição da direção do loop
    while ((lambda!=-1) && ((isMaxtree && lambda > oldGrayLevel) || (!isMaxtree && lambda < oldGrayLevel))) {
        std::vector<NodeId>& F_lambda = F.getMergedNodes(lambda);

        nodeUnion = tree->proxy(F_lambda.front());
        disconnect(tree, nodeUnion, false);
    
        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion && nodeId != newNode) {
                nodeUnion.addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId), tree);
                mergedParentAndChildren(tree, nodeUnion, nodeId);
                disconnect(tree, nodeId, true);
                //tree->releaseNode(nodeId);
            }
        }
        
        if (lambda == newGrayLevel) {
            if(tree->getFlatZonesGraph()->isAdjacent(repFZ, nodeUnion.getRepCNPs())){
                nodeUnion.addCNPsToConnectedFlatzone(repFZ, tree); 
            } 
            else {
                nodeUnion.addRepCNPs(repFZ);
            }
            //disconnect(tree, newNode, true);
            //tree->releaseNode(newNode);
            nodeTauL.removeFlatzone(repFZ); 
            for (NodeId nodeId : F.getFb()) {
                disconnect(tree, nodeId, false);
                tree->addChildById(nodeUnion, nodeId);
            }
        }
        if (nodeUnionPrevious) {
            nodeUnionPrevious.setParent(nodeUnion);
            nodeUnion.addChild(nodeUnionPrevious);
        }

        // Atualiza atributo de área
        nodeUnion.setArea(nodeUnion.getNumCNPs());
        for(NodeId nodeId: tree->getChildrenById(nodeUnion)){
            nodeUnion.setArea(nodeUnion.getArea() + tree->getAreaById(nodeId));
        }
    

        
        if(PRINT_LOG){
            outputLog << "\tnodeUnion = union(F_{" << lambda << "}) = " << 
                " id:" << nodeUnion.getIndex() << 
                ", level: " << nodeUnion.getLevel() << 
                ", |cnps|: " << nodeUnion.getNumCNPs() <<  
                ", |children|: " << nodeUnion.getNumChildren() << std::endl;
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }

    if (nodeTauCNPsIsEqualL) {
        NodeFZ parentNodeTauL = nodeTauL.getParent();
        nodeUnion.setParent(parentNodeTauL);

        if (parentNodeTauL) {
            parentNodeTauL.addChild(nodeUnion);
            for(NodeId n: tree->getChildrenById(nodeTauL)){
                if (n != nodeUnion && !tree->hasChildById(nodeUnion, n)) {
                    tree->addChildById(nodeUnion, n);
                    nodeUnion.setArea(nodeUnion.getArea() + tree->getAreaById(n));
                }
            }
            disconnect(tree, nodeTauL, true);
        } 
        else {  // Novo root
            NodeId newRoot = nodeUnion;
            if (nodeTauL.getNumChildren() > 0) {
                for(NodeId n: tree->getChildrenById(nodeTauL)){
                    if ( (isMaxtree && tree->getLevelById(n) < tree->getLevelById(newRoot)) || (!isMaxtree && tree->getLevelById(n) > tree->getLevelById(newRoot))) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    tree->addChildById(newRoot, nodeUnion);
                }
                for(NodeId n: tree->getChildrenById(nodeTauL)){
                    if (n != newRoot && !tree->hasChildById(nodeUnion, n)) {
                        tree->addChildById(newRoot, n);
                    }
                }

            }
            
            tree->setAreaById(newRoot, nodeTauL.getArea());
            tree->setRootById(newRoot);
            tree->releaseNode(nodeTauL);
        }
        
    } else {
        if (nodeUnion != nodeTauL) {
            nodeUnion.setParent(nodeTauL);
            nodeTauL.addChild(nodeUnion);
        }

    }
}



void ComponentTreeAdjustmentByFlatzone::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, int repFlatzone) {
    updateTree(mintree, repFlatzone); 
    mergeWithParent(maxtree, repFlatzone); 
}

void ComponentTreeAdjustmentByFlatzone::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, int repFlatzone) {
    updateTree(maxtree, repFlatzone);   
    mergeWithParent(mintree, repFlatzone); 
}
