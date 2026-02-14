#include "../include/ComponentTreeAdjustmentBySubtree.hpp"


template<typename Computer, typename GraphT>
void ComponentTreeAdjustmentBySubtree<Computer, GraphT>::updateTree(ComponentTreeFZ<GraphT>* tree, NodeId rootIdSubtree) {
    assert(rootIdSubtree!=InvalidNode && "rootSubtree is invalid"); 
    bool isMaxtree = tree->isMaxtree();

    ComponentTreeFZ<GraphT>* otherTree = isMaxtree ? this->mintree : this->maxtree;
    this->areaFZsRemoved = otherTree->getAreaById(rootIdSubtree);
    int newGrayLevel = otherTree->getLevelById( otherTree->getParentById(rootIdSubtree) );  // b = g(p)
    
    properPartsCollector.resetCollections(isMaxtree);
    for (NodeId nSubtree : otherTree->getIteratorBreadthFirstTraversalById(rootIdSubtree)) {
        for(int repFZ : otherTree->getRepCNPsById(nSubtree)) {
            properPartsCollector.addNode(tree, tree->getSCById(repFZ), repFZ); 
        }
    }
    
    NodeId nodeIdTauStar = properPartsCollector.getNodeTauStar();
    int pixelUpperBound = properPartsCollector.getRepFZTauStar();
    int grayTauStar = tree->getLevelById(nodeIdTauStar);
    if(PRINT_LOG){
        this->outputLog.clear();
        this->outputLog << "Area(rSubtree)= " << otherTree->getAreaById(rootIdSubtree) << ", level(rSubtree)= " << otherTree->getLevelById(rootIdSubtree) << ", level(parent(rSubtree))= " << otherTree->getLevelById(otherTree->getParentById(rootIdSubtree)) << std::endl;
        this->outputLog << "newGrayLevel: " << newGrayLevel  << std::endl; //g(p)
        this->outputLog << "Proper parts: (Tau_S): [";
        bool flagPrint = false;
        for(int repFZTau: properPartsCollector.getRepsFZ()){
            NodeFZ<GraphT> nodeTau = tree->getSC(repFZTau);
            if(flagPrint){
                this->outputLog << "\t";
            }
            flagPrint = true;
            this->outputLog << "\t(id:" << nodeTau.getIndex() << ", level:" << nodeTau.getLevel() << ", |cnps|:"<< nodeTau.getNumCNPs() << ", repFZ:" << repFZTau << "), \n";
        }
        this->outputLog << "]"<< std::endl;
        if(tree->isMaxtree())
            this->outputLog << "Intervalo: [" << grayTauStar << ", " << newGrayLevel << "]" << std::endl;
        else
            this->outputLog << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
        this->outputLog << "nodeTauStar: Id:" << nodeIdTauStar << "; level:" << grayTauStar <<"; |cnps|:" << tree->getNumCNPsById(nodeIdTauStar) << std::endl;
    }    


    
    this->buildMergedAndNestedCollections(tree,  properPartsCollector.getRepsFZ(), pixelUpperBound, newGrayLevel, isMaxtree);

    if(PRINT_LOG){
        this->outputLog << "F_λ = { ";
        for(int lambda=newGrayLevel; lambda != grayTauStar; ){
            std::vector<NodeId>& F_lambda = this->F.getMergedNodes(lambda);
            if(!F_lambda.empty()){
                this->outputLog << lambda << ":[ ";
                for(NodeId node: F_lambda){
                    this->outputLog << "Id:" << node << " ";
                }
                this->outputLog << "] ";
            }
            if(tree->isMaxtree()) lambda--; else lambda++;       
        }
        this->outputLog << "}\nF_{λ>b} = {";
        for (NodeId node : this->F.getFb()) {
            this->outputLog << " Id:" << node << " ";
        }
        this->outputLog << "}\n" << std::endl;
        
    }    
    

    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = this->F.firstLambda();
    NodeId nodeUnion;
    NodeId nodeUnionPrevious = InvalidNode;
    size_t mergedCalls = 0;
    size_t totalMerged = 0;
    long long areaSum = 0;
    long long areaMax = 0;
    size_t areaCount = 0;
    size_t adjacentCount = 0;
    long long areaTauStar = 0;
    size_t loopIterations = 0;

    auto mergedParentAndChildrenTracked = [&](ComponentTreeFZ<GraphT>* targetTree, NodeId parentId, NodeId childId) {
        this->mergedParentAndChildren(targetTree, parentId, childId);
        if (this->metrics_) {
            ++mergedCalls;
        }
    };

    if (this->metrics_) {
        this->pauseMetrics();
        areaTauStar = (nodeIdTauStar != InvalidNode) ? tree->getAreaById(nodeIdTauStar) : 0;
        for (int l = 0; l <= 255; ++l) {
            const auto& nodes = this->F.getMergedNodes(l);
            totalMerged += nodes.size();
            for (NodeId nodeId : nodes) {
                long long area = tree->getAreaById(nodeId);
                areaSum += area;
                if (area > areaMax) {
                    areaMax = area;
                }
                ++areaCount;
            }
        }
        adjacentCount = this->F.getAdjacentNodes().size();
        this->resumeMetrics();
    }


    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        if (this->metrics_) {
            ++loopIterations;
        }
        std::vector<NodeId>& F_lambda = this->F.getMergedNodes(lambda);
        
        // Encontrar um nodeUnion que NÃO esteja em nodesToBeRemoved
        nodeUnion = InvalidNode;
        for (NodeId nodeId : F_lambda) {
            if (!properPartsCollector.isRemoved(nodeId)) { 
                nodeUnion = nodeId;
                break; 
            }
        }

        
        // Se não encontrou nenhum node válido, continua para a próxima iteração
        if (nodeUnion == InvalidNode) {
            for (NodeId nodeId : F_lambda) {
                mergedParentAndChildrenTracked(tree, tree->getParentById(nodeId), nodeId);
                this->disconnect(tree, nodeId, true);
            }
            lambda = this->F.nextLambda();  
            nodeUnion = nodeUnionPrevious;
            continue;
        }
        if(PRINT_LOG){
            this->outputLog << "F_{" << lambda << "} = \n";
            this->outputLog << "\t(Id:" << nodeUnion << "; level:" << tree->getLevelById(nodeUnion) <<"; |cnps|:" << tree->getNumCNPsById(nodeUnion) << ") " << std::endl;
        }

        this->disconnect(tree, nodeUnion, false);

        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion) {
                if(PRINT_LOG){
                    this->outputLog << "\t(Id:" << nodeId << "; level:" << tree->getLevelById(nodeId) <<"; |cnps|:" << tree->getNumCNPsById(nodeId) << ") " << std::endl;
                }

                if(!properPartsCollector.isRemoved(nodeId)){ //node não foi removido
                    NodeFZ<GraphT>::addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId), nodeUnion, tree);    
                }else{
                    if(PRINT_LOG)
                        this->outputLog << "\t\twas removed" << std::endl;
                }
                    
                mergedParentAndChildrenTracked(tree, nodeUnion, nodeId);
                this->disconnect(tree, nodeId, true);
            }
        }
        if (lambda == newGrayLevel) {
            properPartsCollector.addCNPsToConnectedFlatzone(nodeUnion, tree); // Os mapeamentos são atualizados
            
            if(PRINT_LOG){
                this->outputLog << "\t\tAfter add CNPs of S: (Id:" << nodeUnion << "; level:" << tree->getLevelById(nodeUnion) <<"; |cnps|:" << tree->getNumCNPsById(nodeUnion) << ") " << std::endl;
            }
            for (NodeId nodeId : this->F.getFb()) {
                this->disconnect(tree, nodeId, false);
                tree->addChildById(nodeUnion, nodeId);
            }
            
            if(PRINT_LOG){
                auto nodesToBeRemoved = properPartsCollector.getNodesToBeRemoved(); 
                if(!nodesToBeRemoved.empty()){
                    this->outputLog << "\tNodes to be removed from tree: ";
                    for(NodeId nodeId: nodesToBeRemoved){
                        this->outputLog << 
                            "(id:" << nodeId << 
                            ", level: " << tree->getLevelById(nodeId) << 
                            ", |cnps|: " << tree->getNumCNPsById(nodeId) <<  
                            ", |children|: " << tree->getNumChildrenById(nodeId) << "), ";
                    }
                    this->outputLog << std::endl;
                }
            }
        }
        
        if (nodeUnionPrevious != InvalidNode) {
            tree->setParentById(nodeUnionPrevious, nodeUnion);
            tree->addChildById(nodeUnion, nodeUnionPrevious);
        }
        
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
        
        if(PRINT_LOG)
            this->outputLog << "\tnodeUnion = union(F_{" << lambda << "}) = " << 
                " id:" << nodeUnion << 
                ", level: " << tree->getLevelById(nodeUnion) << 
                ", |cnps|: " << tree->getNumCNPsById(nodeUnion) <<  
                ", |children|: " << tree->getNumChildrenById(nodeUnion) << std::endl;
            
        nodeUnionPrevious = nodeUnion;
        lambda = this->F.nextLambda();  
    }
    

    if (properPartsCollector.isRemoved(nodeIdTauStar)) { //nodeTauStar foi removido
        NodeId parentIdNodeTauStar = tree->getParentById(nodeIdTauStar);
        tree->setParentById(nodeUnion, parentIdNodeTauStar);
        
        if (parentIdNodeTauStar != InvalidNode) {
            tree->addChildById(parentIdNodeTauStar, nodeUnion);
            for(NodeId n: tree->getChildrenById(nodeIdTauStar)){
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
            this->disconnect(tree, nodeIdTauStar, true);
        } 
        else {  // Novo root
            NodeId newRoot = nodeUnion;
            if (tree->getNumChildrenById(nodeIdTauStar) > 0) {
                for(NodeId n: tree->getChildrenById(nodeIdTauStar)){
                    if ( (isMaxtree && tree->getLevelById(n) < tree->getLevelById(newRoot)) || (!isMaxtree && tree->getLevelById(n) > tree->getLevelById(newRoot))) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    tree->addChildById(newRoot, nodeUnion);
                }
                for(NodeId n: tree->getChildrenById(nodeIdTauStar)){
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

            tree->setAreaById(newRoot, tree->getAreaById(nodeIdTauStar));
            tree->setRootById(newRoot);
            tree->releaseNode(nodeIdTauStar);
        }
        
    } else {
        if (nodeUnion != nodeIdTauStar) {
            tree->setParentById(nodeUnion, nodeIdTauStar);
            tree->addChildById(nodeIdTauStar, nodeUnion);
        }   
    }

    if (this->metrics_) {
        this->pauseMetrics();
        this->metrics_->reset();
        this->metrics_->a = newGrayLevel;
        this->metrics_->b = grayTauStar;
        this->metrics_->area_substree = this->areaFZsRemoved;
        this->metrics_->count_proper_parts = properPartsCollector.getRepsFZ().size();
        this->metrics_->count_nodes_Fb = this->F.getFb().size();
        this->metrics_->count_total_nodes_merged = totalMerged;
        this->metrics_->count_total_nodes_removed = properPartsCollector.getRemovedCount();
        this->metrics_->count_total_nodes_and_children_merged = mergedCalls;
        this->metrics_->count_adjacent_nodes = adjacentCount;
        this->metrics_->area_tau_star = areaTauStar;
        this->metrics_->loop_iterations = loopIterations;
        this->metrics_->avg_area_nodes_merged = areaCount ? (static_cast<double>(areaSum) / static_cast<double>(areaCount)) : 0.0;
        this->metrics_->max_area_nodes_merged = areaMax;
        this->metrics_->valid = true;
        this->resumeMetrics();
    }
}


template<typename Computer, typename GraphT>
void ComponentTreeAdjustmentBySubtree<Computer, GraphT>::adjustMinTree(ComponentTreeFZ<GraphT>* mintree,
                                                                       ComponentTreeFZ<GraphT>* maxtree,
                                                                       std::vector<NodeId>& nodesToPruning) {
    for (NodeId rSubtreeId : nodesToPruning) {  	
        assert(rSubtreeId != maxtree->getRootById() && "rSubtree is root");
        if (rSubtreeId == InvalidNode) {
            continue; //rSubtree é root, não pode ser ajustado
        }
        updateTree(mintree, rSubtreeId); 
        this->prunning(maxtree, rSubtreeId);
    }
}

template<typename Computer, typename GraphT>
void ComponentTreeAdjustmentBySubtree<Computer, GraphT>::adjustMaxTree(ComponentTreeFZ<GraphT>* maxtree,
                                                                       ComponentTreeFZ<GraphT>* mintree,
                                                                       std::vector<NodeId>& nodesToPruning) {
    for (NodeId rSubtreeId : nodesToPruning) {  
        assert(rSubtreeId != mintree->getRootById() && "rSubtree is root");
        if (rSubtreeId == InvalidNode) {
            continue; //rSubtree é root, não pode ser ajustado
        }
        updateTree(maxtree, rSubtreeId);   
        this->prunning(mintree, rSubtreeId);
    }
}

template class ComponentTreeAdjustmentBySubtree<DefaultAttributeComputer, DefaultFlatZonesGraph>;
template class ComponentTreeAdjustmentBySubtree<BoundingBoxComputerFZ, DefaultFlatZonesGraph>;
template class ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<FlatZonesGraphFullEdges>, FlatZonesGraphFullEdges>;
template class ComponentTreeAdjustmentBySubtree<BoundingBoxComputerFZT<FlatZonesGraphFullEdges>, FlatZonesGraphFullEdges>;
template class ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<FlatZonesGraphOnDemandEdgesByPixel>, FlatZonesGraphOnDemandEdgesByPixel>;
template class ComponentTreeAdjustmentBySubtree<BoundingBoxComputerFZT<FlatZonesGraphOnDemandEdgesByPixel>, FlatZonesGraphOnDemandEdgesByPixel>;
