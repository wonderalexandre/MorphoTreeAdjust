#include "../include/ComponentTreeAdjustmentBySubtree.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>


void ComponentTreeAdjustmentBySubtree::updateTree(ComponentTreeFZPtr tree, NodeFZ rootSubtree) {
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
        outputLog.clear();
        outputLog << "Area(rSubtree)= " << rootSubtree.getArea() << ", level(rSubtree)= " << rootSubtree.getLevel() << ", level(parent(rSubtree))= " << rootSubtree.getParent().getLevel() << std::endl;
        outputLog << "newGrayLevel: " << newGrayLevel  << std::endl; //g(p)
        outputLog << "Proper parts: (Tau_S): [";
        bool flagPrint = false;
        for(int repFZTau: properPartsCollector.getRepsFZ()){
            NodeFZ nodeTau = tree->getSC(repFZTau);
            if(flagPrint){
                outputLog << "\t";
            }
            flagPrint = true;
            outputLog << "\t(id:" << nodeTau.getIndex() << ", level:" << nodeTau.getLevel() << ", |cnps|:"<< nodeTau.getNumCNPs() << ", repFZ:" << repFZTau << "), \n";
        }
        outputLog << "]"<< std::endl;
        if(tree->isMaxtree())
            outputLog << "Intervalo: [" << grayTauStar << ", " << newGrayLevel << "]" << std::endl;
        else
            outputLog << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
        outputLog << "nodeTauStar: Id:" << nodeTauStar.getIndex() << "; level:" << nodeTauStar.getLevel() <<"; |cnps|:" << nodeTauStar.getNumCNPs() << std::endl;
    }    


    
    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree,  properPartsCollector.getRepsFZ(), pixelUpperBound, newGrayLevel, isMaxtree);

    if(PRINT_LOG){
        outputLog << "F_λ = { ";
        for(int lambda=newGrayLevel; lambda != grayTauStar; ){
            std::vector<NodeId>& F_lambda = F.getMergedNodes(lambda);
            if(!F_lambda.empty()){
                outputLog << lambda << ":[ ";
                for(NodeId node: F_lambda){
                    outputLog << "Id:" << node << " ";
                }
                outputLog << "] ";
            }
            if(tree->isMaxtree()) lambda--; else lambda++;       
        }
        outputLog << "}\nF_{λ>b} = {";
        for (NodeId node : F.getFb()) {
            outputLog << " Id:" << node << " ";
        }
        outputLog << "}\n" << std::endl;
        
    }    


    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    NodeFZ nodeUnion;
    NodeFZ nodeUnionPrevious;
    NodeFZ nodeTauParentSubtree;
    
    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        std::vector<NodeId>& F_lambda = F.getMergedNodes(lambda);
        
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
                mergedParentAndChildren(tree, tree->getParentById(nodeId), nodeId);
                disconnect(tree, nodeId, true);
                //tree->releaseNode(nodeId);
            }
            lambda = F.nextLambda();  
            nodeUnion = nodeUnionPrevious;
            continue;
        }
        if(PRINT_LOG){
            outputLog << "F_{" << lambda << "} = \n";
            outputLog << "\t(Id:" << nodeUnion.getIndex() << "; level:" << nodeUnion.getLevel() <<"; |cnps|:" << nodeUnion.getNumCNPs() << ") " << std::endl;
        }

        disconnect(tree, nodeUnion, false);

        for (NodeId nodeId : F_lambda) {
            if (nodeId != nodeUnion) {
                if(PRINT_LOG){
                    outputLog << "\t(Id:" << nodeId << "; level:" << tree->getLevelById(nodeId) <<"; |cnps|:" << tree->getNumCNPsById(nodeId) << ") " << std::endl;
                }

                if(!properPartsCollector.isRemoved(nodeId)){ //node não foi removido
                    nodeUnion.addCNPsOfDisjointFlatzones(tree->getRepCNPsById(nodeId), tree);    
                }else{
                    if(PRINT_LOG)
                        outputLog << "\t\twas removed" << std::endl;
                }
                    
                mergedParentAndChildren(tree, nodeUnion, nodeId);
                disconnect(tree, nodeId, true);
                //tree->releaseNode(nodeId);
            }
        }
        if (lambda == newGrayLevel) {
            properPartsCollector.addCNPsToConnectedFlatzone(nodeUnion, tree); // Os mapeamentos são atualizados
            
            if(PRINT_LOG){
                outputLog << "\t\tAfter add CNPs of S: (Id:" << nodeUnion.getIndex() << "; level:" << nodeUnion.getLevel() <<"; |cnps|:" << nodeUnion.getNumCNPs() << ") " << std::endl;
            }
            for (NodeId nodeId : F.getFb()) {
                disconnect(tree, nodeId, false);
                tree->addChildById(nodeUnion, nodeId);
            }
            nodeTauParentSubtree = nodeUnion;
            
            if(PRINT_LOG){
                auto nodesToBeRemoved = properPartsCollector.getNodesToBeRemoved(); 
                if(!nodesToBeRemoved.empty()){
                    outputLog << "\tNodes to be removed from tree: ";
                    for(NodeId nodeId: nodesToBeRemoved){
                        outputLog << 
                            "(id:" << nodeId << 
                            ", level: " << tree->getLevelById(nodeId) << 
                            ", |cnps|: " << tree->getNumCNPsById(nodeId) <<  
                            ", |children|: " << tree->getNumChildrenById(nodeId) << "), ";
                    }
                    outputLog << std::endl;
                }
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
        if(PRINT_LOG)
            outputLog << "\tnodeUnion = union(F_{" << lambda << "}) = " << 
                " id:" << nodeUnion.getIndex() << 
                ", level: " << nodeUnion.getLevel() << 
                ", |cnps|: " << nodeUnion.getNumCNPs() <<  
                ", |children|: " << nodeUnion.getNumChildren() << std::endl;
            
        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  
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
                }
            }
            disconnect(tree, nodeTauStar, true);
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


void ComponentTreeAdjustmentBySubtree::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId rSubtreeId : nodesToPruning) {  	
        NodeFZ rSubtree = maxtree->proxy(rSubtreeId);	
        assert(rSubtree != maxtree->getRoot() && "rSubtree is root");
        if (!rSubtree) {
            continue; //rSubtree é root, não pode ser ajustado
        }
        updateTree(mintree, rSubtree); 
        prunning(maxtree, rSubtree);
    }
}

void ComponentTreeAdjustmentBySubtree::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToPruning) {
    for (NodeId rSubtreeId : nodesToPruning) {  
        NodeFZ rSubtree = mintree->proxy(rSubtreeId);	
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");
        if (!rSubtree) {
            continue; //rSubtree é root, não pode ser ajustado
        }
        updateTree(maxtree, rSubtree);   
        prunning(mintree, rSubtree);
    }
}