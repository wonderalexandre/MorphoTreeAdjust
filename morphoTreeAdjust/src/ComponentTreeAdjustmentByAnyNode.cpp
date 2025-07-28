#include "../include/ComponentTreeAdjustmentByAnyNode.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <algorithm> 
#include <utility>




void ComponentTreeAdjustmentByAnyNode::updateTree(ComponentTreeFZPtr tree, NodeFZPtr node) {
    assert(node != nullptr && "rootSubtree is nullptr"); 
    ComponentTreeFZPtr otherTree = tree->isMaxtree()? this->mintree : this->maxtree;

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = node->getParent()->getLevel();  // g(p)

    unionNodeTauSubtree.resetCollections(isMaxtree);
    for (NodeFZPtr nSubtree : node->getIteratorBreadthFirstTraversal()) {
        for(auto& [idFlatZoneNSubtree, fzSubtree]: nSubtree->getCNPsByFlatZone()){    
            NodeFZPtr nodeTau = tree->getSC(idFlatZoneNSubtree);
            FlatZone& fzTau = nodeTau->getFlatZone(idFlatZoneNSubtree); 
            unionNodeTauSubtree.addNode(nodeTau, fzTau); 
        }
    }
    NodeFZPtr nodeTauStar = unionNodeTauSubtree.getNodeTauStar();
    int pixelUpperBound = unionNodeTauSubtree.getFlatzoneIDTauStar();
    int grayTauStar = nodeTauStar->getLevel();  // f(pixelUpperBound)
    if(PRINT_LOG){
        outputLog.clear();
        outputLog << "Area(node)= " << node->getArea() << ", level(node)= " << node->getLevel() << ", level(parent(node))= " << node->getParent()->getLevel() << std::endl;
        outputLog << "newGrayLevel: " << newGrayLevel  << std::endl; //g(p)
        outputLog << "unionNodes (Tau_S): [";
        bool flagPrint = false;
        for(int fzTauID: unionNodeTauSubtree.getFlatzonesID()){
            NodeFZPtr nodeTau = tree->getSC(fzTauID);
            if(flagPrint){
                outputLog << "\t";
            }
            flagPrint = true;
            outputLog << "\t(id:" << nodeTau->getIndex() << ", level:" << nodeTau->getLevel() << ", |cnps|:"<< nodeTau->getNumCNPs() << ", idFZ:" << fzTauID << "), \n";
        }
        outputLog << "]"<< std::endl;
        if(tree->isMaxtree())
            outputLog << "Intervalo: [" << grayTauStar << ", " << newGrayLevel << "]" << std::endl;
        else
            outputLog << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
        outputLog << "nodeTauStar: Id:" << nodeTauStar->getIndex() << "; level:" << nodeTauStar->getLevel() <<"; |cnps|:" << nodeTauStar->getNumCNPs() << std::endl;
    }    


    
    ComponentTreeAdjustment::buildMergedAndNestedCollections(tree,  unionNodeTauSubtree.getFlatzonesID(), pixelUpperBound, newGrayLevel, isMaxtree);

    if(PRINT_LOG){
        outputLog << "F_λ = { ";
        for(int lambda=newGrayLevel; lambda != grayTauStar; ){
            std::vector<NodeFZPtr>& F_lambda = F.getMergedNodes(lambda);
            if(!F_lambda.empty()){
                outputLog << lambda << ":[ ";
                for(NodeFZPtr node: F_lambda){
                    outputLog << "Id:" << node->getIndex() << " ";
                }
                outputLog << "] ";
            }
            if(tree->isMaxtree()) lambda--; else lambda++;       
        }
        outputLog << "}\nF_{λ>b} = {";
        for (NodeFZPtr node : Fb) {
            outputLog << " Id:" << node->getIndex() << " ";
        }
        outputLog << "}\n" << std::endl;
    }    


    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    if(lambda == -1) { 
        //se não houver lambdas, não há o que atualizar
        std::cout << "=========Não há lambdas =======" << std::endl;
        return;
    }

    NodeFZPtr nodeUnion = nullptr;
    NodeFZPtr nodeUnionPrevious = nullptr;
    NodeFZPtr nodeTauParentSubtree = nullptr;
    
    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        std::vector<NodeFZPtr>& F_lambda = F.getMergedNodes(lambda);
        
        // Encontrar um nodeUnion que NÃO esteja em nodesToBeRemoved
        nodeUnion = nullptr;
        for (NodeFZPtr node : F_lambda) {
            if (!unionNodeTauSubtree.isRemoved(node)) { 
                nodeUnion = node;
                break; 
            }
        }
        // Se não encontrou nenhum node válido, continua para a próxima iteração
        if (!nodeUnion) {
            for (NodeFZPtr node : F_lambda) {
                mergedParentAndChildren(node->getParent(), node);
                disconnect(node, true);
                tree->setNumNodes(tree->getNumNodes() - 1);
            }
            lambda = F.nextLambda();  
            nodeUnion = nodeUnionPrevious;
            continue;
        }
        if(PRINT_LOG){
            outputLog << "F_{" << lambda << "} = \n";
            outputLog << "\t(Id:" << nodeUnion->getIndex() << "; level:" << nodeUnion->getLevel() <<"; |cnps|:" << nodeUnion->getNumCNPs() << ") " << std::endl;
        }

        disconnect(nodeUnion);

        for (NodeFZPtr n : F_lambda) {
            
            if (n != nodeUnion) {
                if(PRINT_LOG){
                    outputLog << "\t(Id:" << n->getIndex() << "; level:" << n->getLevel() <<"; |cnps|:" << n->getNumCNPs() << ") " << std::endl;
                }

                if(!unionNodeTauSubtree.isRemoved(n)){ //node não foi removido
                    nodeUnion->addCNPsOfDisjointFlatzones(n->moveCNPsByFlatZone(), tree);    
                }else{
                    if(PRINT_LOG)
                        outputLog << "\t\twas removed" << std::endl;
                }
                    
                mergedParentAndChildren(nodeUnion, n);
                disconnect(n, true);
                tree->setNumNodes(tree->getNumNodes() - 1);
            }
        }
        if (lambda == newGrayLevel) {
            unionNodeTauSubtree.addCNPsToConnectedFlatzone(nodeUnion, tree); // Os mapeamentos são atualizados
            
            if(PRINT_LOG){
                outputLog << "\t\tAfter add CNPs of S: (Id:" << nodeUnion->getIndex() << "; level:" << nodeUnion->getLevel() <<"; |cnps|:" << nodeUnion->getNumCNPs() << ") " << std::endl;
            }
            for (NodeFZPtr n : Fb) {
                disconnect(n);
                nodeUnion->addChild(n);
                n->setParent(nodeUnion);
            }
            nodeTauParentSubtree = nodeUnion;
            
            if(PRINT_LOG){
                auto nodesToBeRemoved = unionNodeTauSubtree.getNodesToBeRemoved(); 
                if(!nodesToBeRemoved.empty()){
                    outputLog << "\tNodes to be removed from tree: ";
                    for(NodeFZPtr node: nodesToBeRemoved){
                        outputLog << 
                            "(id:" << node->getIndex() << 
                            ", level: " << node->getLevel() << 
                            ", |cnps|: " << node->getNumCNPs() <<  
                            ", |children|: " << node->getChildren().size() << "), ";
                    }
                    outputLog << std::endl;
                }
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
        if(PRINT_LOG)
            outputLog << "\tnodeUnion = union(F_{" << lambda << "}) = " << 
                " id:" << nodeUnion->getIndex() << 
                ", level: " << nodeUnion->getLevel() << 
                ", |cnps|: " << nodeUnion->getNumCNPs() <<  
                ", |children|: " << nodeUnion->getChildren().size() << std::endl;
            
        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  
    }
    

    if (unionNodeTauSubtree.isRemoved(nodeTauStar)) { //nodeTauStar foi removido
        NodeFZPtr parentNodeTauStar = nodeTauStar->getParent();
        nodeUnion->setParent(parentNodeTauStar);

        if (parentNodeTauStar != nullptr) {
            parentNodeTauStar->addChild(nodeUnion);
            for (NodeFZPtr n : nodeTauStar->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } else {  // Novo root
            NodeFZPtr newRoot = nodeUnion;
            if (!nodeTauStar->getChildren().empty()) {
                for (NodeFZPtr n : nodeTauStar->getChildren()) {
                    if ((isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeFZPtr n : nodeTauStar->getChildren()) {
                    if (n != newRoot && !nodeUnion->isChild(n)) {
                        newRoot->addChild(n);
                        n->setParent(newRoot);
                    }
                }
            }
            
            newRoot->setArea(nodeTauStar->getArea());
            newRoot->setParent(nullptr);
            tree->setRoot(newRoot);
            
        }
        tree->setNumNodes( tree->getNumNodes() -1 );  
        disconnect(nodeTauStar, true);
    } else {
        if (nodeUnion != nodeTauStar) {
            nodeUnion->setParent(nodeTauStar);
            nodeTauStar->addChild(nodeUnion);
        }   
    }
}




void ComponentTreeAdjustmentByAnyNode::adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeFZPtr>& nodesToRemoved) {
    for (NodeFZPtr node : nodesToRemoved) {  	
        assert(node != maxtree->getRoot() && "node is root");
        updateTree(mintree, node); 
        maxtree->mergeWithParent(node); 
    }
}

void ComponentTreeAdjustmentByAnyNode::adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeFZPtr>& nodesToRemoved) {
    for (NodeFZPtr node : nodesToRemoved) {  
        assert(node != mintree->getRoot() && "node is root");

        updateTree(maxtree, node);   
        mintree->mergeWithParent(node); 
    }
}
