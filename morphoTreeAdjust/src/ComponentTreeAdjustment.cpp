#include "../include/ComponentTreeAdjustment.hpp"
#include <unordered_set>
#include <list>
#include <vector>
#include <iostream>
#include <functional>
#include <utility>


ComponentTreeAdjustment::ComponentTreeAdjustment(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree) 
    : maxtree(maxtree), mintree(mintree), 
      maxIndex(std::max(maxtree->getNumNodes(), mintree->getNumNodes())),
      F(maxIndex),
      unionNodeTauSubtree(maxtree->isMaxtree()) 
{    }

ComponentTreeAdjustment::~ComponentTreeAdjustment() { }



void ComponentTreeAdjustment::buildMergedAndNestedCollections(ComponentTreeFZ* tree, std::vector<FlatZoneRef>& flatZone, int newGrayLevel, bool isMaxtree){
	Fb.clear();
	F.resetCollection(isMaxtree);
    if(this->pixelUpperBound == -1){
        this->pixelUpperBound = flatZone.front().get().front();
    }
    F.computerAdjacentNodes(tree, flatZone);
    NodeFZ* nodeTauL = tree->getSC(this->pixelUpperBound); //pixel de tauStar ou tauL, para termos o node (limite) mais proximo de root

    for (NodeFZ* nodeNL: F.getAdjacentNodes()) {
	    if( (isMaxtree && nodeNL->getLevel() <= newGrayLevel) || (!isMaxtree &&  nodeNL->getLevel() >= newGrayLevel)) { //o nodeNL está entre g(p) e f(p)
		    F.addNodesOfPath(nodeNL, nodeTauL); 
	    } 
	    else { 
            //o nodeNL está abaixo de g(p). É armazenado somente a raiz da subtree antes de atingir o nivel g(p)
	        NodeFZ* nodeSubtree = nodeNL;
            for (NodeFZ *n : nodeNL->getNodesOfPathToRoot()) {
                if ( (isMaxtree && newGrayLevel > n->getLevel()) || (!isMaxtree && newGrayLevel < n->getLevel())) {
                    break;
                }
                nodeSubtree = n; 
            }
	        // se a subtree tiver level = g(p), então ela entra em F[\lambda]
            if (nodeSubtree->getLevel() == newGrayLevel) {
                F.addNodesOfPath(nodeSubtree, nodeTauL); //F_lambda
            } 
            else {
                Fb.insert(nodeSubtree); //F_{lambda} > b
            }
	    }
        
	}
}


void ComponentTreeAdjustment::updateTree2(ComponentTreeFZ* tree, NodeFZ* rootSubtree) {
    assert(rootSubtree != nullptr && "rootSubtree is nullptr"); 
    ComponentTreeFZ* otherTree = tree->isMaxtree()? this->mintree : this->maxtree;

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = rootSubtree->getParent()->getLevel();  // g(p)
    
    unionNodeTauSubtree.resetCollection(isMaxtree);
    for (NodeFZ* nSubtree : rootSubtree->getIteratorBreadthFirstTraversal()) {
        for(auto& [idFlatZoneNSubtree, fzSubtree]: nSubtree->getCNPsByFlatZone()){    
            NodeFZ* nodeTau = tree->getSC(idFlatZoneNSubtree);
            FlatZone& fzTau = nodeTau->getFlatZone(idFlatZoneNSubtree); //tree->getFlatzoneByID(idFlatZoneNSubtree);
            unionNodeTauSubtree.addNode(nodeTau, fzTau); //, fzSubtree.size() == fzTau.size()  
        }
    }
    NodeFZ* nodeTauStar = unionNodeTauSubtree.getNodeTauStar().node;
    FlatZone* fzTauStar = unionNodeTauSubtree.getNodeTauStar().flatzone;
    this->pixelUpperBound = fzTauStar->front();
    int grayTauStar = nodeTauStar->getLevel();  // f(pixelUpperBound)
    if(PRINT_LOG){
        outputLog.clear();
        outputLog << "Area(rSubtree)= " << rootSubtree->getArea() << ", level(rSubtree)= " << rootSubtree->getLevel() << ", level(parent(rSubtree))= " << rootSubtree->getParent()->getLevel() << std::endl;
        outputLog << "newGrayLevel: " << newGrayLevel  << std::endl; //g(p)
        outputLog << "unionNodes (Tau_S): [";
        bool flagPrint = false;
        for(FlatZoneNode& fzNode: unionNodeTauSubtree.getFlatzoneNodeList()){
            NodeFZ* nodeTau = fzNode.node;
            FlatZone* fzTau = fzNode.flatzone;
            if(flagPrint){
                outputLog << "\t";
            }
            flagPrint = true;
            outputLog << "\t(id:" << nodeTau->getIndex() << ", level:" << nodeTau->getLevel() << ", |cnps|:"<< nodeTau->getNumCNPs() << ", idFZ:" << fzTau->front() <<", |fz|:" << fzTau->size() << "), \n";
        }
        outputLog << "]"<< std::endl;
        if(tree->isMaxtree())
            outputLog << "Intervalo: [" << grayTauStar << ", " << newGrayLevel << "]" << std::endl;
        else
            outputLog << "Intervalo: [" << newGrayLevel << ", " << grayTauStar << "]" << std::endl;
        outputLog << "nodeTauStar: Id:" << nodeTauStar->getIndex() << "; level:" << nodeTauStar->getLevel() <<"; |cnps|:" << nodeTauStar->getNumCNPs() <<"; |flatZoneTauStar|:" << fzTauStar->size() << std::endl;
    }    


    //bool nodeTauCNPsIsEqualL = (nodeTauStar->getNumFlatzone() ==  1);
    std::vector<FlatZoneRef>  flatzonesTauSubtree = unionNodeTauSubtree.getFlatzones();
    this->buildMergedAndNestedCollections(tree,  flatzonesTauSubtree, newGrayLevel, isMaxtree);

    if(PRINT_LOG){
        outputLog << "F_λ = { ";
        for(int lambda=newGrayLevel; lambda != grayTauStar; ){
            std::vector<NodeFZ*>& F_lambda = F.getMergedNodes(lambda);
            if(!F_lambda.empty()){
                outputLog << lambda << ":[ ";
                for(NodeFZ* node: F_lambda){
                    outputLog << "Id:" << node->getIndex() << " ";
                }
                outputLog << "] ";
            }
            if(tree->isMaxtree()) lambda--; else lambda++;       
        }
        outputLog << "}\nF_{λ>b} = {";
        for (NodeFZ* node : Fb) {
            outputLog << " Id:" << node->getIndex() << " ";
        }
        outputLog << "}\n" << std::endl;
    }    


    // Ordenação dos lambdas (crescente para Min-Tree, decrescente para Max-Tree)
    int lambda = F.firstLambda();
    NodeFZ* nodeUnion = nullptr;
    NodeFZ* nodeUnionPrevious = nullptr;
    NodeFZ* nodeTauParentSubtree = nullptr;
    
    // Definição da direção do loop
    while ((isMaxtree && lambda > grayTauStar) || (!isMaxtree && lambda < grayTauStar)) {
        std::vector<NodeFZ*>& F_lambda = F.getMergedNodes(lambda);
        
        // Encontrar um nodeUnion que NÃO esteja em nodesToBeRemoved
        nodeUnion = nullptr;
        for (NodeFZ* node : F_lambda) {
            if (unionNodeTauSubtree.isRemoved(node)) { 
                nodeUnion = node;
                break; 
            }
        }
        // Se não encontrou nenhum node válido, continua para a próxima iteração
        if (!nodeUnion) {
            for (NodeFZ* node : F_lambda) {
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

        for (NodeFZ* n : F_lambda) {
            
            if (n != nodeUnion) {
                if(PRINT_LOG){
                    outputLog << "\t(Id:" << n->getIndex() << "; level:" << n->getLevel() <<"; |cnps|:" << n->getNumCNPs() << ") " << std::endl;
                }

                if(unionNodeTauSubtree.isRemoved(n)){ //node não foi removido
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
            unionNodeTauSubtree.removeFlatzones(tree);
            if(PRINT_LOG){
                outputLog << "\t\tAfter add CNPs of S: (Id:" << nodeUnion->getIndex() << "; level:" << nodeUnion->getLevel() <<"; |cnps|:" << nodeUnion->getNumCNPs() << ") " << std::endl;
            }
            for (NodeFZ* n : Fb) {
                disconnect(n);
                nodeUnion->addChild(n);
                n->setParent(nodeUnion);
            }
            nodeTauParentSubtree = nodeUnion;
            //nodesToBeRemoved = unionNodeTauSubtree.getNodesToBeRemoved();
            
            if(PRINT_LOG){
                if(!unionNodeTauSubtree.getNodesToBeRemoved().empty()){
                    outputLog << "\tNodes to be removed from tree: ";
                    for(NodeFZ* node: unionNodeTauSubtree.getNodesToBeRemoved()){
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
        for (NodeFZ* n : nodeUnion->getChildren()) {
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
    

    if (!unionNodeTauSubtree.isRemoved(nodeTauStar)) { //nodeTauStar não foi removido
        NodeFZ* parentNodeTauStar = nodeTauStar->getParent();
        nodeUnion->setParent(parentNodeTauStar);

        if (parentNodeTauStar != nullptr) {
            parentNodeTauStar->addChild(nodeUnion);
            for (NodeFZ* n : nodeTauStar->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } else {  // Novo root
            NodeFZ* newRoot = nodeUnion;
            if (!nodeTauStar->getChildren().empty()) {
                for (NodeFZ* n : nodeTauStar->getChildren()) {
                    if ((isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeFZ* n : nodeTauStar->getChildren()) {
                    if (n != newRoot && !nodeUnion->isChild(n)) {
                        newRoot->addChild(n);
                        n->setParent(newRoot);
                    }
                }
            }
            
            newRoot->setArea(nodeTauStar->getArea());
            newRoot->setParent(nullptr);
            tree->setRoot(nodeUnion);
            
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


void ComponentTreeAdjustment::updateTree(ComponentTreeFZ* tree, NodeFZ* leaf) {
    assert(leaf != nullptr && "L_leaf is nullptr"); 

    bool isMaxtree = tree->isMaxtree();
    int newGrayLevel = leaf->getParent()->getLevel();  // b = g(p)
    int oldGrayLevel = leaf->getLevel();  // a = f(p)
    int idLeaf = leaf->getRepresentativeCNPs(); //pixel (id) of flatzone 
    
    NodeFZ* nodeTauL = tree->getSC(idLeaf); //node of correspondence flatzone in other treee
    this->pixelUpperBound = idLeaf; 

    bool nodeTauCNPsIsEqualL = nodeTauL->getNumFlatzone() == 1;
    FlatZone& flatzoneTauL = tree->getFlatzoneByID(idLeaf); 
    std::vector<FlatZoneRef> flatZonesTauL = {flatzoneTauL};
    
    assert(leaf->getNumCNPs() == flatzoneTauL.size() && "O número de CNPs de L_leaf é diferente do número de pixels da flatzone de tauL");

    this->buildMergedAndNestedCollections(tree, flatZonesTauL, newGrayLevel, isMaxtree);

    int lambda = F.firstLambda(); //star with b = newGrayLevel
    NodeFZ* nodeUnion = nullptr; // tau_{lambda}
    NodeFZ* nodeUnionPrevious = nullptr; //maxtree: tau_{\lambda+1}, mintree:  tau_{\lambda-1}

    // Definição da direção do loop
    while ( (isMaxtree && lambda > oldGrayLevel) || (!isMaxtree && lambda < oldGrayLevel)) {
        std::vector<NodeFZ*>& F_lambda = F.getMergedNodes(lambda);

        nodeUnion = F_lambda.front();
        disconnect(nodeUnion);

        for (NodeFZ* n : F_lambda) {
            if (n != nodeUnion) {
                nodeUnion->addCNPsOfDisjointFlatzones(n->moveCNPsByFlatZone(), tree);
                mergedParentAndChildren(nodeUnion, n);
                disconnect(n, true);
                tree->setNumNodes(tree->getNumNodes() - 1);
            }
        }
        if (lambda == newGrayLevel) {
            int idFlatzoneTauL = flatzoneTauL.front();
            nodeUnion->addCNPsToConnectedFlatzone(std::move(flatzoneTauL), tree); // Os mapeamentos são atualizados
            nodeTauL->removeFlatzone(idFlatzoneTauL);
            for (NodeFZ* n : Fb) {
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
        for (NodeFZ* n : nodeUnion->getChildren()) {
            nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
        }

        nodeUnionPrevious = nodeUnion;
        lambda = F.nextLambda();  // Avança para o próximo lambda
    }


    if (nodeTauCNPsIsEqualL) {
        NodeFZ* parentNodeTauL = nodeTauL->getParent();
        nodeUnion->setParent(parentNodeTauL);

        if (parentNodeTauL != nullptr) {
            parentNodeTauL->addChild(nodeUnion);
            for (NodeFZ* n : nodeTauL->getChildren()) {
                if (n != nodeUnion && !nodeUnion->isChild(n)) {
                    nodeUnion->getChildren().push_back(n);
                    nodeUnion->setArea(nodeUnion->getArea() + n->getArea());
                    n->setParent(nodeUnion);
                }
            }
        } 
        else {  // Novo root
            NodeFZ* newRoot = nodeUnion;
            if (!nodeTauL->getChildren().empty()) {
                for (NodeFZ* n : nodeTauL->getChildren()) {
                    if ( (isMaxtree && n->getLevel() < newRoot->getLevel()) || (!isMaxtree && n->getLevel() > newRoot->getLevel())) {
                        newRoot = n;
                    }
                }
                if (newRoot != nodeUnion) {
                    newRoot->addChild(nodeUnion);
                    nodeUnion->setParent(newRoot);
                }
                for (NodeFZ* n : nodeTauL->getChildren()) {
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


void ComponentTreeAdjustment::adjustMinTree(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* node : nodesToPruning) {
        for (NodeFZ* Lmax : node->getIteratorPostOrderTraversal()) { 
            assert(Lmax != maxtree->getRoot() && "Lmax is root");
            assert(Lmax->isLeaf() && "Lmax não é uma folha");
   
            updateTree(mintree, Lmax); 
            maxtree->prunning(Lmax);
        }
    }
}

void ComponentTreeAdjustment::adjustMaxTree(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* node : nodesToPruning) {  	
        for (NodeFZ* Lmin : node->getIteratorPostOrderTraversal()) {
            assert(Lmin != mintree->getRoot() && "Lmin is root");
            assert(Lmin->isLeaf() && "Lmin não é uma folha");

            updateTree(maxtree, Lmin);             
            mintree->prunning(Lmin);
        }
    }
}

void ComponentTreeAdjustment::adjustMinTree2(ComponentTreeFZ* mintree, ComponentTreeFZ* maxtree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* rSubtree : nodesToPruning) {  	
        assert(rSubtree != maxtree->getRoot() && "rSubtree is root");
        updateTree2(mintree, rSubtree); 
        maxtree->prunning(rSubtree);
    }
}

void ComponentTreeAdjustment::adjustMaxTree2(ComponentTreeFZ* maxtree, ComponentTreeFZ* mintree, std::vector<NodeFZ*> nodesToPruning) {
    for (NodeFZ* rSubtree : nodesToPruning) {  
        assert(rSubtree != mintree->getRoot() && "rSubtree is root");

        updateTree2(maxtree, rSubtree);   
        mintree->prunning(rSubtree);
    }
}