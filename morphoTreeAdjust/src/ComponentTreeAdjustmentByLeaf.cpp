#include "../include/ComponentTreeAdjustmentByLeaf.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>



void ComponentTreeAdjustmentByLeaf::updateTree(ComponentTreeFZPtr tree, NodeFZ leaf) {
    assert(leaf != nullptr && "L_leaf is nullptr"); 
    assert(leaf.isLeaf() && "L_leaf is not leaf"); 
    assert(leaf.getParent() && "L_leaf is root"); 
    
    this->areaFZsRemoved = leaf.getArea();

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = leaf.getParent().getLevel();  // b = g(p)
    int oldGrayLevel = leaf.getLevel();  // a = f(p)
    int idLeaf = leaf.getRepCNPs().front(); //pixel (id) of flatzone 
    NodeFZ nodeTauL = tree->getSC(idLeaf); //node of correspondence flatzone in other treee
    
    bool nodeTauCNPsIsEqualL = nodeTauL.getNumFlatzone() == 1;
    int idFlatzoneTauL = idLeaf;//flatzoneTauL.front();
    
    assert(leaf.getNumCNPs(tree->getFlatZonesGraph()) == flatzoneTauL && "O número de CNPs de L_leaf é diferente do número de pixels da flatzone de tauL");

    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree, idFlatzoneTauL, newGrayLevel, isMaxtree);

    int lambda = F.firstLambda(); //start with b = newGrayLevel
    NodeFZ nodeUnion; // tau_{lambda}
    NodeFZ nodeUnionPrevious; //maxtree: tau_{\lambda+1}, mintree:  tau_{\lambda-1}

    // Definição da direção do loop
    while ( (isMaxtree && lambda > oldGrayLevel) || (!isMaxtree && lambda < oldGrayLevel)) {
        std::vector<NodeId>& F_lambda = F.getMergedNodes(lambda);

        nodeUnion = tree->proxy(F_lambda.front());    
        disconnect(tree, nodeUnion);

        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion) {
                nodeUnion.addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId), tree);
                mergedParentAndChildren(tree, nodeUnion, nodeId);
                disconnect(tree, nodeId);
                //tree->setNumNodes(tree->getNumNodes() - 1);
                tree->releaseNode(nodeId);
            }
        }
        if (lambda == newGrayLevel) {
            nodeUnion.addCNPsToConnectedFlatzone(idFlatzoneTauL, tree); 
            nodeTauL.removeFlatzone(idFlatzoneTauL);
            for (NodeId nodeId : F.getFb()) {
                disconnect(tree, nodeId);
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
            
        }
        
        disconnect(tree, nodeTauL);
        //tree->setNumNodes(tree->getNumNodes() - 1);
        tree->releaseNode(nodeTauL);
    } else {
        if (nodeUnion != nodeTauL) {
            nodeUnion.setParent(nodeTauL);
            nodeTauL.addChild(nodeUnion);
        }

    }
}


void ComponentTreeAdjustmentByLeaf::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId nodeId : nodesToPruning) {        
        NodeFZ node = maxtree->proxy(nodeId);
        for (NodeFZ Lmax : node.getIteratorPostOrderTraversal()) { 
            assert(Lmax != maxtree->getRoot() && "Lmax is root");
            assert(Lmax.isLeaf() && "Lmax não é uma folha");
            
            if(!Lmax.getParent()) {
                continue; //Lmax é root, não pode ser ajustado
            }
            updateTree(mintree, Lmax);     
            prunning(maxtree, Lmax);
        }
    }
}

void ComponentTreeAdjustmentByLeaf::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId nodeId : nodesToPruning) {  
        NodeFZ node = mintree->proxy(nodeId);	
        for (NodeFZ Lmin : node.getIteratorPostOrderTraversal()) { 
            assert(Lmin != mintree->getRoot() && "Lmin is root");
            assert(Lmin.isLeaf() && "Lmin não é uma folha");
            
            if(! Lmin.getParent()) {
                continue; //Lmax é root, não pode ser ajustado
            }

            updateTree(maxtree, Lmin);             
            prunning(mintree, Lmin);
        }
    }
}
