#include "../include/ComponentTreeAdjustmentBySubtree.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>


template<typename Computer>
void ComponentTreeAdjustmentBySubtree<Computer>::updateTree(ComponentTreeFZPtr tree, NodeFZ rootSubtree) {
    assert(rootSubtree && "rootSubtree is invalid"); 
    ComponentTreeFZPtr otherTree = tree->isMaxtree()? this->mintree : this->maxtree;
    
    this->areaFZsRemoved = rootSubtree.getArea();
    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = rootSubtree.getParent().getLevel();  // g(p)
    
    properPartsCollector.resetCollections(isMaxtree);
    
    for (NodeId nSubtree : otherTree->getIteratorBreadthFirstTraversalById(rootSubtree)) {
        for(int repFZ : otherTree->getRepCNPsById(nSubtree)) {
            properPartsCollector.addNode(tree, tree->getSCById(repFZ), repFZ); 
        }
    }
    
    NodeFZ nodeTauStar = tree->proxy(properPartsCollector.getNodeTauStar());
    int pixelUpperBound = properPartsCollector.getRepFZTauStar();
    int grayTauStar = nodeTauStar.getLevel();  // f(pixelUpperBound)
    if(PRINT_LOG){
        this->outputLog.clear();
        this->outputLog << "Area(rSubtree)= " << rootSubtree.getArea() << ", level(rSubtree)= " << rootSubtree.getLevel() << ", level(parent(rSubtree))= " << rootSubtree.getParent().getLevel() << std::endl;
        this->outputLog << "newGrayLevel: " << newGrayLevel  << std::endl; //g(p)
        this->outputLog << "Proper parts: (Tau_S): [";
        bool flagPrint = false;
        for(int repFZTau: properPartsCollector.getRepsFZ()){
            NodeFZ nodeTau = tree->getSC(repFZTau);
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
        this->outputLog << "nodeTauStar: Id:" << nodeTauStar.getIndex() << "; level:" << nodeTauStar.getLevel() <<"; |cnps|:" << nodeTauStar.getNumCNPs() << std::endl;
    }    


    
    ComponentTreeAdjustment<Computer>::buildMergedAndNestedCollections(tree,  properPartsCollector.getRepsFZ(), pixelUpperBound, newGrayLevel, isMaxtree);

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
    int lambda = this->F.firstLambda(grayTauStar, newGrayLevel);
    NodeFZ nodeUnion;
    NodeFZ nodeUnionPrevious;
    NodeFZ nodeTauParentSubtree;
    
    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        std::vector<NodeId>& F_lambda = this->F.getMergedNodes(lambda);
        
        // Encontrar um nodeUnion que NÃO esteja em nodesToBeRemoved
        nodeUnion = NodeFZ();
        for (NodeId nodeId : F_lambda) {
            if (!properPartsCollector.isRemoved(nodeId)) { 
                nodeUnion = tree->proxy(nodeId);;
                break; 
            }
        }
        // Se não encontrou nenhum node válido, continua para a próxima iteração
        if (!nodeUnion) {
            for (NodeId nodeId : F_lambda) {
                this->mergedParentAndChildren(tree, tree->getParentById(nodeId), nodeId);
                this->disconnect(tree, nodeId, true);
                //tree->releaseNode(nodeId);
            }
            lambda = this->F.nextLambda();  
            nodeUnion = nodeUnionPrevious;
            continue;
        }
        if(PRINT_LOG){
            this->outputLog << "F_{" << lambda << "} = \n";
            this->outputLog << "\t(Id:" << nodeUnion.getIndex() << "; level:" << nodeUnion.getLevel() <<"; |cnps|:" << nodeUnion.getNumCNPs() << ") " << std::endl;
        }

        this->disconnect(tree, nodeUnion, false);

        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion) {
                if(PRINT_LOG){
                    this->outputLog << "\t(Id:" << nodeId << "; level:" << tree->getLevelById(nodeId) <<"; |cnps|:" << tree->getNumCNPsById(nodeId) << ") " << std::endl;
                }

                if(!properPartsCollector.isRemoved(nodeId)){ //node não foi removido
                    nodeUnion.addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId), tree);    
                }else{
                    if(PRINT_LOG)
                        this->outputLog << "\t\twas removed" << std::endl;
                }
                    
                this->mergedParentAndChildren(tree, nodeUnion, nodeId);
                this->disconnect(tree, nodeId, true);
                //tree->releaseNode(nodeId);
            }
        }
        if (lambda == newGrayLevel) {
            properPartsCollector.addCNPsToConnectedFlatzone(nodeUnion, tree); // Os mapeamentos são atualizados
            
            if(PRINT_LOG){
                this->outputLog << "\t\tAfter add CNPs of S: (Id:" << nodeUnion.getIndex() << "; level:" << nodeUnion.getLevel() <<"; |cnps|:" << nodeUnion.getNumCNPs() << ") " << std::endl;
            }
            for (NodeId nodeId : this->F.getFb()) {
                this->disconnect(tree, nodeId, false);
                tree->addChildById(nodeUnion, nodeId);
            }
            nodeTauParentSubtree = nodeUnion;
            
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

        
        if (nodeUnionPrevious) {
            nodeUnionPrevious.setParent(nodeUnion);
            nodeUnion.addChild(nodeUnionPrevious);
        }
        
        
        nodeUnion.setArea(nodeUnion.getNumCNPs());
        for(NodeId nodeId: tree->getChildrenById(nodeUnion)){
            nodeUnion.setArea(nodeUnion.getArea() + tree->getAreaById(nodeId));
        }
        /********************************************************
         **   Atualização inremental de atributos do nodeUnion **
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
                " id:" << nodeUnion.getIndex() << 
                ", level: " << nodeUnion.getLevel() << 
                ", |cnps|: " << nodeUnion.getNumCNPs() <<  
                ", |children|: " << nodeUnion.getNumChildren() << std::endl;
            
        nodeUnionPrevious = nodeUnion;
        lambda = this->F.nextLambda();  
    }
    

    if (properPartsCollector.isRemoved(nodeTauStar)) { //nodeTauStar foi removido
        NodeFZ parentNodeTauStar = nodeTauStar.getParent();
        nodeUnion.setParent(parentNodeTauStar);

        if (parentNodeTauStar) {
            parentNodeTauStar.addChild(nodeUnion);
            for(NodeId n: tree->getChildrenById(nodeTauStar)){
                if (n != nodeUnion && !tree->hasChildById(nodeUnion, n)) {
                    tree->addChildById(nodeUnion, n);
                    nodeUnion.setArea(nodeUnion.getArea() + tree->getAreaById(n));
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
            this->disconnect(tree, nodeTauStar, true);
        } 
        else {  // Novo root
            NodeId newRoot = nodeUnion;
            if (nodeTauStar.getNumChildren() > 0) {
                for(NodeId n: tree->getChildrenById(nodeTauStar)){
                    if ( (isMaxtree && tree->getLevelById(n) < tree->getLevelById(newRoot)) || (!isMaxtree && tree->getLevelById(n) > tree->getLevelById(newRoot))) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    tree->addChildById(newRoot, nodeUnion);
                }
                for(NodeId n: tree->getChildrenById(nodeTauStar)){
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
            
            tree->setAreaById(newRoot, nodeTauStar.getArea());
            tree->setRootById(newRoot);
            tree->releaseNode(nodeTauStar);
        }
        
    } else {
        if (nodeUnion != nodeTauStar) {
            nodeUnion.setParent(nodeTauStar);
            nodeTauStar.addChild(nodeUnion);
        }   
    }
}


template<typename Computer>
void ComponentTreeAdjustmentBySubtree<Computer>::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId rSubtreeId : nodesToPruning) {  	
        NodeFZ rSubtree = maxtree->proxy(rSubtreeId);	
        assert(rSubtree != maxtree->getRoot() && "rSubtree is root");
        if (!rSubtree.getParent()) {
            continue; //rSubtree é root, não pode ser ajustado
        }
        updateTree(mintree, rSubtree); 
        this->prunning(maxtree, rSubtree);
    }
}

template<typename Computer>
void ComponentTreeAdjustmentBySubtree<Computer>::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId rSubtreeId : nodesToPruning) {  
        NodeFZ rSubtree = mintree->proxy(rSubtreeId);	
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");
        if (!rSubtree.getParent()) {
            continue; //rSubtree é root, não pode ser ajustado
        }
        updateTree(maxtree, rSubtree);   
        this->prunning(mintree, rSubtree);
    }
}

template class ComponentTreeAdjustmentBySubtree<DefaultAttributeComputer>;
