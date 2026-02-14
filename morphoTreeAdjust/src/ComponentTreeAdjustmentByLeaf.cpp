#include "../include/ComponentTreeAdjustmentByLeaf.hpp"



template<typename Computer, typename GraphT>
void ComponentTreeAdjustmentByLeaf<Computer, GraphT>::updateTree(ComponentTreeFZ<GraphT>* tree, NodeId leafId) {
    assert(leafId!=InvalidNode && "leafId is invalid");
    assert(tree->isLeafById(leafId) && "leafId não é uma folha");
    bool isMaxtree = tree->isMaxtree();

    ComponentTreeFZ<GraphT>* otherTree = isMaxtree ? this->mintree : this->maxtree;
    this->areaFZsRemoved = otherTree->getAreaById(leafId);

    int newGrayLevel = otherTree->getLevelById( otherTree->getParentById(leafId) );  // b = g(p)    
    
    
    int repLeaf = otherTree->getRepCNPsById(leafId).front(); //pixel (id) of flatzone
    NodeId nodeTauL = tree->getSCById(repLeaf); //node of correspondence flatzone in other treee
    int oldGrayLevel = tree->getLevelById(nodeTauL);  // a = f(p)
    
    bool nodeTauCNPsIsEqualL = tree->getNumFlatzoneById(nodeTauL) == 1;
    int idFlatzoneTauL = repLeaf;
    
    this->buildMergedAndNestedCollections(tree, idFlatzoneTauL, newGrayLevel, isMaxtree);

    int lambda = this->F.firstLambda(); //start with b = newGrayLevel
    NodeId nodeUnion; // tau_{lambda}
    NodeId nodeUnionPrevious = InvalidNode; //maxtree: tau_{\lambda+1}, mintree:  tau_{\lambda-1}

    // Definição da direção do loop
    while ( (isMaxtree && lambda > oldGrayLevel) || (!isMaxtree && lambda < oldGrayLevel)) {
        std::vector<NodeId>& F_lambda = this->F.getMergedNodes(lambda);

        nodeUnion = F_lambda.front(); // Inicializa com o primeiro nó de F_λ
        this->disconnect(tree, nodeUnion, false);

        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion) {
                NodeFZ<GraphT>::addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId), nodeUnion, tree);    
                this->mergedParentAndChildren(tree, nodeUnion, nodeId);
                this->disconnect(tree, nodeId, true);
            }
        }
        if (lambda == newGrayLevel) {
            NodeFZ<GraphT>::addCNPsToConnectedFlatzone(idFlatzoneTauL, nodeUnion, tree);
            NodeFZ<GraphT>::removeFlatzone(idFlatzoneTauL, nodeTauL, tree);

            for (NodeId nodeId : this->F.getFb()) {
                this->disconnect(tree, nodeId, false);
                tree->addChildById(nodeUnion, nodeId);
            }
        }
        if (nodeUnionPrevious != InvalidNode) {
            tree->setParentById(nodeUnionPrevious, nodeUnion);
            tree->addChildById(nodeUnion, nodeUnionPrevious);
        }
        

        // O atributo área SEMPRE é atualizado
        tree->setAreaById(nodeUnion, tree->getNumCNPsById(nodeUnion));
        for(NodeId nodeId: tree->getChildrenById(nodeUnion)){
            tree->setAreaById(nodeUnion, tree->getAreaById(nodeUnion) + tree->getAreaById(nodeId));
        }
        
        /********************************************************
         **   Atualização incremental de atributos do nodeUnion **
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
                " id:" << nodeUnion << 
                ", level: " << tree->getLevelById(nodeUnion) << 
                ", |cnps|: " << tree->getNumCNPsById(nodeUnion) <<  
                ", |children|: " << tree->getNumChildrenById(nodeUnion) << std::endl;
        }

        nodeUnionPrevious = nodeUnion;
        lambda = this->F.nextLambda();  // Avança para o próximo lambda
    }

    if (nodeTauCNPsIsEqualL) {
        NodeId parentNodeTauL = tree->getParentById(nodeTauL); 
        tree->setParentById(nodeUnion, parentNodeTauL);
        
        if (parentNodeTauL != InvalidNode) {
            tree->addChildById(parentNodeTauL, nodeUnion);
            for(NodeId n: tree->getChildrenById(nodeTauL)){
                if (n != nodeUnion && !tree->hasChildById(nodeUnion, n)) {
                    tree->addChildById(nodeUnion, n);
                    tree->setAreaById(nodeUnion, tree->getAreaById(nodeUnion) + tree->getAreaById(n));
                    
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
            if (tree->getNumChildrenById(nodeTauL) > 0) {
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



            tree->setAreaById(newRoot, tree->getAreaById(nodeTauL)); //nesse caso: O novo root pode herda de nodeTauL  
            tree->setRootById(newRoot);
            tree->releaseNode(nodeTauL);
            
        }
        
    } else {
        if (nodeUnion != nodeTauL) {
            tree->setParentById(nodeUnion, nodeTauL);
            tree->addChildById(nodeTauL, nodeUnion);
        }

    }
}


template<typename Computer, typename GraphT>
void ComponentTreeAdjustmentByLeaf<Computer, GraphT>::adjustMinTree(ComponentTreeFZ<GraphT>* mintree,
                                                                    ComponentTreeFZ<GraphT>* maxtree,
                                                                    std::vector<NodeId>& nodesToPruning) {
    for (NodeId nodeId : nodesToPruning) {        
        for (NodeId LmaxId : maxtree->getIteratorPostOrderTraversalById(nodeId)) { 
            assert(LmaxId != maxtree->getRootById(LmaxId) && "Lmax is root");
            assert(maxtree->isLeafById(LmaxId) && "Lmax não é uma folha");
            
            if (LmaxId == InvalidNode) {
                continue; //Lmax é root, não pode ser ajustado
            }
            updateTree(mintree, LmaxId);     
            this->prunning(maxtree, LmaxId);
        }
    }
}

template<typename Computer, typename GraphT>
void ComponentTreeAdjustmentByLeaf<Computer, GraphT>::adjustMaxTree(ComponentTreeFZ<GraphT>* maxtree,
                                                                    ComponentTreeFZ<GraphT>* mintree,
                                                                    std::vector<NodeId>& nodesToPruning) {
    for (NodeId nodeId : nodesToPruning) {  
        for (NodeId LminId : mintree->getIteratorPostOrderTraversalById(nodeId)) { 
            assert(LminId != mintree->getRootById(LminId) && "Lmin is root");
            assert(mintree->isLeafById(LminId) && "Lmin não é uma folha");
            if (LminId == InvalidNode) {
                continue; //Lmin é root, não pode ser ajustado
            }
            updateTree(maxtree, LminId);             
            this->prunning(mintree, LminId);
        }
    }
}

template class ComponentTreeAdjustmentByLeaf<DefaultAttributeComputer, DefaultFlatZonesGraph>;
template class ComponentTreeAdjustmentByLeaf<BoundingBoxComputerFZ, DefaultFlatZonesGraph>;
template class ComponentTreeAdjustmentByLeaf<DefaultAttributeComputerT<FlatZonesGraphFullEdges>, FlatZonesGraphFullEdges>;
template class ComponentTreeAdjustmentByLeaf<BoundingBoxComputerFZT<FlatZonesGraphFullEdges>, FlatZonesGraphFullEdges>;
template class ComponentTreeAdjustmentByLeaf<DefaultAttributeComputerT<FlatZonesGraphOnDemandEdgesByPixel>, FlatZonesGraphOnDemandEdgesByPixel>;
template class ComponentTreeAdjustmentByLeaf<BoundingBoxComputerFZT<FlatZonesGraphOnDemandEdgesByPixel>, FlatZonesGraphOnDemandEdgesByPixel>;
