#include "../include/ComponentTreeAdjustmentByLeaf.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>



template<typename Computer>
void ComponentTreeAdjustmentByLeaf<Computer>::updateTree(ComponentTreeFZ* tree, NodeFZ leaf) {
    assert(leaf && "L_leaf is nullptr"); 
    assert(leaf.isLeaf() && "L_leaf is not leaf"); 
    assert(leaf.getParent() && "L_leaf is root"); 
    
    this->areaFZsRemoved = leaf.getArea();

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = leaf.getParent().getLevel();  // b = g(p)
    int oldGrayLevel = leaf.getLevel();  // a = f(p)
    int repLeaf = leaf.getRepCNPs().front(); //pixel (id) of flatzone 
    NodeFZ nodeTauL = tree->getSC(repLeaf); //node of correspondence flatzone in other treee
    
    bool nodeTauCNPsIsEqualL = nodeTauL.getNumFlatzone() == 1;
    int idFlatzoneTauL = repLeaf;
    
    ComponentTreeAdjustment<Computer>::buildMergedAndNestedCollections(tree, idFlatzoneTauL, newGrayLevel, isMaxtree);

    int lambda = this->F.firstLambda(oldGrayLevel, newGrayLevel); //start with b = newGrayLevel
    NodeFZ nodeUnion; // tau_{lambda}
    NodeFZ nodeUnionPrevious; //maxtree: tau_{\lambda+1}, mintree:  tau_{\lambda-1}

    // Definição da direção do loop
    while ( (isMaxtree && lambda > oldGrayLevel) || (!isMaxtree && lambda < oldGrayLevel)) {
        std::vector<NodeId>& F_lambda = this->F.getMergedNodes(lambda);

        nodeUnion = tree->proxy(F_lambda.front());    
        this->disconnect(tree, nodeUnion, false);

        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion) {
                nodeUnion.addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId));
                this->mergedParentAndChildren(tree, nodeUnion, nodeId);
                this->disconnect(tree, nodeId, true);
            }
        }
        if (lambda == newGrayLevel) {
            nodeUnion.addCNPsToConnectedFlatzone(idFlatzoneTauL); 
            nodeTauL.removeFlatzone(idFlatzoneTauL);
            for (NodeId nodeId : this->F.getFb()) {
                this->disconnect(tree, nodeId, false);
                tree->addChildById(nodeUnion, nodeId);
            }
        }
        if (nodeUnionPrevious) {                
            nodeUnionPrevious.setParent(nodeUnion);
            nodeUnion.addChild(nodeUnionPrevious);
        }
        

        //O atributo area SEMPRE é atualizado
        nodeUnion.setArea(nodeUnion.getNumCNPs()); //1. pre-processing
        for(NodeId nodeId: tree->getChildrenById(nodeUnion)){
            nodeUnion.setArea(nodeUnion.getArea() + tree->getAreaById(nodeId)); //2. merge-processing
        }
        
        /********************************************************
         **  Atualização incremental de atributos do nodeUnion **
         ********************************************************/
        if(isMaxtree){
            this->attrComputerMax->preProcessing(nodeUnion, this->bufferMax);
            for(NodeId nodeId: tree->getChildrenById(nodeUnion)){
                this->attrComputerMax->mergeProcessing(nodeUnion, nodeId, this->bufferMax);
            }
            this->attrComputerMax->postProcessing(nodeUnion, this->bufferMax);
        } else {
            this->attrComputerMin->preProcessing(nodeUnion, this->bufferMin);
            for(NodeId nodeId: tree->getChildrenById(nodeUnion)){
                this->attrComputerMin->mergeProcessing(nodeUnion, nodeId, this->bufferMin);
            }
            this->attrComputerMin->postProcessing(nodeUnion, this->bufferMin);
        }
        /********************************************************/

        if(PRINT_LOG){
            this->outputLog << "\tnodeUnion = union(F_{" << lambda << "}) = " << 
                " id:" << nodeUnion.getIndex() << 
                ", level: " << nodeUnion.getLevel() << 
                ", |cnps|: " << nodeUnion.getNumCNPs() <<  
                ", |children|: " << nodeUnion.getNumChildren() << std::endl;
        }

        nodeUnionPrevious = nodeUnion;
        lambda = this->F.nextLambda();  // Avança para o próximo lambda
    }

    if (nodeTauCNPsIsEqualL) {
        NodeFZ parentNodeTauL = nodeTauL.getParent();
        nodeUnion.setParent(parentNodeTauL);

        if (parentNodeTauL) {
            parentNodeTauL.addChild(nodeUnion);
            for(NodeId n: tree->getChildrenById(nodeTauL)){
                if (n != nodeUnion && !tree->hasChildById(nodeUnion, n)) {
                    tree->addChildById(nodeUnion, n);
                    
                    nodeUnion.setArea(nodeUnion.getArea() + tree->getAreaById(n)); //atualização de atributo: merge-processing
                    
                    /*****************************************************************
                     **  Atualização incremental de atributos do (ultimo) nodeUnion **
                    ******************************************************************/
                    if(isMaxtree){
                        this->attrComputerMax->mergeProcessing(nodeUnion, n, this->bufferMax);
                    } else {
                        this->attrComputerMin->mergeProcessing(nodeUnion, n, this->bufferMin);
                    }
                    /******************************************************************/
                }
            }
            this->disconnect(tree, nodeTauL, true);
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
            

            /*******************************************************
             **  Atualização incremental de atributos do newRoot  **
            ********************************************************/
            if(isMaxtree){
                this->attrComputerMax->preProcessing(newRoot, this->bufferMax);                
                for(NodeId nodeId: tree->getChildrenById(newRoot)){
                    this->attrComputerMax->mergeProcessing(newRoot, nodeId, this->bufferMax);
                }
                this->attrComputerMax->postProcessing(newRoot, this->bufferMax);
            }else{
                this->attrComputerMin->preProcessing(newRoot, this->bufferMin);                
                for(NodeId nodeId: tree->getChildrenById(newRoot)){
                    this->attrComputerMin->mergeProcessing(newRoot, nodeId, this->bufferMin);
                }
                this->attrComputerMin->postProcessing(newRoot, this->bufferMin);
            }
            /********************************************************/



            tree->setAreaById(newRoot, nodeTauL.getArea()); //nesse caso: O novo root pode herda de nodeTauL  
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


template<typename Computer>
void ComponentTreeAdjustmentByLeaf<Computer>::adjustMinTree(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId nodeId : nodesToPruning) {        
        NodeFZ node = maxtree->proxy(nodeId);
        for (NodeFZ Lmax : node.getIteratorPostOrderTraversal()) { 
            assert(Lmax != maxtree->getRoot() && "Lmax is root");
            assert(Lmax.isLeaf() && "Lmax não é uma folha");
            
            if(!Lmax.getParent()) {
                continue; //Lmax é root, não pode ser ajustado
            }
            updateTree(mintree, Lmax);     
            this->prunning(maxtree, Lmax);
        }
    }
}

template<typename Computer>
void ComponentTreeAdjustmentByLeaf<Computer>::adjustMaxTree(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId nodeId : nodesToPruning) {  
        NodeFZ node = mintree->proxy(nodeId);	
        for (NodeFZ Lmin : node.getIteratorPostOrderTraversal()) { 
            assert(Lmin != mintree->getRoot() && "Lmin is root");
            assert(Lmin.isLeaf() && "Lmin não é uma folha");
            
            if(! Lmin.getParent()) {
                continue; //Lmax é root, não pode ser ajustado
            }

            updateTree(maxtree, Lmin);             
            this->prunning(mintree, Lmin);
        }
    }
}

template class ComponentTreeAdjustmentByLeaf<DefaultAttributeComputer>;
