#include <iterator>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_H
#define COMPONENT_TREE_ADJUSTMENT_H

#include <array>
#include <vector>
#include <unordered_set>

/**
 * @brief Estrutura auxiliar para gerenciar coleções de nós mesclados e adjacentes em árvores de componentes.
 *
 * `MergedNodesCollection` organiza nós por nível (0–255), mantendo listas de nós mesclados (collectionF),
 * adjacentes (nodesNL), e conjuntos auxiliares (Fb). Também controla a ordem de iteração (maxtree/min-tree)
 * e evita duplicações através de marcações (`visited`, `visitedAdj`).
 */
class MergedNodesCollection {
protected:
    std::array<std::vector<NodeId>, 256> collectionF;
    
    int maxIndex;
    GenerationStampSet visited;
    GenerationStampSet visitedAdj;
    
    std::vector<int> lambdaList; // Lista ordenada de lambdas (sempre crescente)
    int currentIndex = 0; // Índice atual
    bool isMaxtree; // Se true, percorre de forma decrescente
    std::vector<NodeId> nodesNL;
    std::vector<NodeId> Fb;

public:
    /// @brief Construtor permite definir a ordem de iteração
    MergedNodesCollection(int maxIndex): maxIndex(maxIndex), visited(maxIndex), visitedAdj(maxIndex) { }

    /// @brief Retorna a lista de nós mesclados armazenados no nível dado.
    std::vector<NodeId>& getMergedNodes(int level) {
        return collectionF[level]; 
    }

    /// @brief Adiciona um nó mesclado (versão NodeFZ).
    void addMergedNode(NodeFZ nodeNL) {
        //if(!visited.isMarked(nodeNL->getIndex())) {
            collectionF[nodeNL.getLevel()].push_back(nodeNL);
            visited.mark(nodeNL.getIndex());
        //}
    }

    /// @brief Adiciona um nó mesclado especificando árvore e id.
    void addMergedNode(ComponentTreeFZPtr tree, NodeId nodeId) {
        //if(!visited.isMarked(nodeNL->getIndex())) {
            collectionF[tree->getLevelById(nodeId)].push_back(nodeId);
            visited.mark(nodeId);
        //}
    }

    /// @brief Computa e armazena nós adjacentes a uma lista de flat-zones.
    void computerAdjacentNodes(ComponentTreeFZPtr tree, const std::vector<int>& flatZonesID) {
        bool isMaxtree = tree->isMaxtree();
        FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
        for (int flatZoneID_P : flatZonesID) {   
            int grayFlatZoneP =  tree->getLevelById(tree->getSCById(flatZoneID_P)); //is same that: level de "a"
    
            for (int flatZoneID_Q : graph->getAdjacentFlatzonesFromPixel(flatZoneID_P)) {
                NodeId nodeId = tree->getSCById(flatZoneID_Q);
                if ( (isMaxtree && tree->getLevelById(nodeId) > grayFlatZoneP) || (!isMaxtree && tree->getLevelById(nodeId) < grayFlatZoneP) ) {
                    if(!visitedAdj.isMarked(nodeId)) {
                        nodesNL.push_back(nodeId);  
                        visitedAdj.mark(nodeId);
                    }
                }
            }
        }
    }

    /// @brief Computa e armazena nós adjacentes a uma flat-zone específica.
    void computerAdjacentNodes(ComponentTreeFZPtr tree, int flatZoneID_P) {
        bool isMaxtree = tree->isMaxtree();
        FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
        
        int grayFlatZoneP =  tree->getLevelById(tree->getSCById(flatZoneID_P)); //is same that: level de "a"

        for (int flatZoneID_Q : graph->getAdjacentFlatzonesFromPixel(flatZoneID_P)) {
            NodeId nodeId = tree->getSCById(flatZoneID_Q);
            if ( (isMaxtree && tree->getLevelById(nodeId) > grayFlatZoneP) || (!isMaxtree && tree->getLevelById(nodeId) < grayFlatZoneP) ) {
                if(!visitedAdj.isMarked(nodeId)) {
                    nodesNL.push_back(nodeId);  
                    visitedAdj.mark(nodeId);
                }
            }
        }
    }

    /// @brief Retorna a lista de nós adjacentes coletados.
    std::vector<NodeId>& getAdjacentNodes(){
        return nodesNL;
    }

    /// @brief Retorna a coleção completa de nós por nível.
    std::array<std::vector<NodeId>, 256>& getCollectionF(){
        return collectionF;
    }

    /// @brief Reinicializa a coleção, limpando listas, visitas e Fb.
    void resetCollection(bool descendingOrder) {
        this->isMaxtree = descendingOrder;

        for(int index: lambdaList) {
            collectionF[index].clear();
        }

        lambdaList.clear();
        visited.resetAll();
        visitedAdj.resetAll();

        nodesNL.clear();
        Fb.clear();

        currentIndex = 0;
    }

    /// @brief Adiciona um nó à lista Fb, se ainda não marcado.
    void addNodesInFb(NodeId nodeNL) {
        if(!visited.isMarked(nodeNL)) {
            Fb.push_back(nodeNL);
            visited.mark(nodeNL);
        }
    }

    /// @brief Retorna a lista Fb de nós armazenados.
    std::vector<NodeId>& getFb() {
        return Fb;
    }

    /// @brief Adiciona todos os nós do caminho até a raiz, parando em nodeTauL.
    void addNodesOfPath(ComponentTreeFZPtr tree, NodeId nodeNL, NodeId nodeTauL) {
        if(!visited.isMarked(nodeNL)){
            for (NodeId nodeId : tree->getNodesOfPathToRootById(nodeNL)) {
                if(!visited.isMarked(nodeId)) {
                    collectionF[tree->getLevelById(nodeId)].push_back(nodeId);
                    visited.mark(nodeId);
                }else{
                    break;
                }
                
                if (nodeId == nodeTauL) {
                    break;
                }
            }
        }
    }
    
    /// @brief Retorna o primeiro valor lambda disponível na coleção.
    int firstLambda() {
        lambdaList.clear();
        for (int i = 0; i < 256; ++i) {
            if (!collectionF[i].empty()) {
                lambdaList.push_back(i);
            }
        }
        if(lambdaList.empty()) {
            return -1; // Retorna -1 se não houver lambdas
        }
        currentIndex = isMaxtree ? lambdaList.size() - 1 : 0;
        return lambdaList[currentIndex];
    }

    /// @brief Avança para o próximo valor lambda.
    int nextLambda() {
        currentIndex = (isMaxtree)? currentIndex-1 : currentIndex+1; 
        
        if(currentIndex < 0 || currentIndex >= static_cast<int>(lambdaList.size())) {
            return -1; // Retorna -1 se não houver mais lambdas
        }else {
            return lambdaList[currentIndex];
        }
    }
};

/**
 * @brief Algoritmos de ajuste em árvores de componentes (min-tree e max-tree).
 *
 * `ComponentTreeAdjustment` gerencia operações de fusão, poda (prunning),
 * e construção de coleções aninhadas/mescladas em árvores de componentes
 * baseadas em flat-zones. É utilizada em algoritmos de ajuste e atualização
 * entre max-tree e min-tree.
 */
class ComponentTreeAdjustment {

protected: 
    ComponentTreeFZPtr mintree;
    ComponentTreeFZPtr maxtree;    

    int maxIndex; 
    MergedNodesCollection F;
    std::ostringstream outputLog;
    long int areaFZsRemoved;
    
    void disconnect(ComponentTreeFZPtr tree, NodeId nodeId, bool releaseNode) {
        NodeId parentId = tree->getParentById(nodeId);
        tree->removeChildById(parentId, nodeId, releaseNode);
        
    }

    void mergedParentAndChildren(ComponentTreeFZPtr tree, NodeId nodeUnion, NodeId n){
        for(NodeId son: tree->getChildrenById(n)){
            tree->setParentById(son, nodeUnion);
        }
        tree->spliceChildrenById(nodeUnion, n);
    }
    
    ComponentTreeFZPtr getOtherTree(bool isMaxtree){
        return isMaxtree ? mintree : maxtree;
    }

    
public:
    /// @brief Retorna o log de saída acumulado.
    std::string getOutputLog() {
        return outputLog.str();
    }

    /// @brief Limpa o log de saída.
    void clearOutputLog() {
        outputLog.str(""); // Limpa o conteúdo do stream
        outputLog.clear();
    }

    /// @brief Construtor principal a partir de min-tree e max-tree.
    ComponentTreeAdjustment(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree): mintree(mintree), maxtree(maxtree), maxIndex(maxtree->getFlatZonesGraph()->getNumFlatZones()), F(maxIndex)  { }

    ComponentTreeAdjustment() = delete; 

    virtual ~ComponentTreeAdjustment() = default;
    
    /// @brief Constrói coleções mescladas e aninhadas a partir de uma lista de flat-zones.
    void buildMergedAndNestedCollections(ComponentTreeFZPtr tree, std::vector<int>& flatZonesID,  int pixelUpperBound, int newGrayLevel, bool isMaxtree){
        
        F.resetCollection(isMaxtree);
        F.computerAdjacentNodes(tree, flatZonesID);
        NodeId nodeTauStar = tree->getSCById(pixelUpperBound); //node tauL: nó que contem L na maxtree

        for (NodeId nodeNL: F.getAdjacentNodes()) { //todos os vizinho de flatzone (folha) L com level > "a" 
            if( (isMaxtree && tree->getLevelById(nodeNL) <= newGrayLevel) || (!isMaxtree &&  tree->getLevelById(nodeNL) >= newGrayLevel)) { 
                //o nodeNL está no intervalo [a+1, b] =>Note que: newGrayLevel é a variavel "b" do paper
                F.addNodesOfPath(tree, nodeNL, nodeTauStar);  
            } 
            else { 
                //O level(nodeNL) > b 
                NodeId nodeSubtree = nodeNL;
                for (NodeId n : tree->getNodesOfPathToRootById(nodeNL)) {
                    if ( (isMaxtree && newGrayLevel > tree->getLevelById(n)) || (!isMaxtree && newGrayLevel < tree->getLevelById(n))) {
                        break;
                    }
                    nodeSubtree = n; //nodeSubtree é a raiz da subarvore antes de atingir o nivel b
                }
                // se a subtree tiver level = b, então ela entra em F[\lambda]
                if ( (tree->getLevelById(nodeSubtree) == newGrayLevel) || (tree->getNumFlatzoneById(nodeTauStar)==1 && tree->getParentById(nodeSubtree) && tree->getParentById(nodeSubtree) != nodeTauStar)) {
                //if(nodeSubtree->getParent() && nodeSubtree->getParent()->getIndex() != nodeTauStar->getIndex() ){ // caso raro: se o pai de nodeSubtree for diferente de tauL
                    F.addNodesOfPath(tree, nodeSubtree, nodeTauStar); //F_lambda
                } 
                else {
                    F.addNodesInFb(nodeSubtree); //F_{lambda} > b
                }
            }
            
        }
    }
    
    /// @brief Constrói coleções mescladas e aninhadas para uma flat-zone específica.
    void buildMergedAndNestedCollections(ComponentTreeFZPtr tree, int flatZoneID, int newGrayLevel, bool isMaxtree){
        F.resetCollection(isMaxtree);
        F.computerAdjacentNodes(tree, flatZoneID);
        NodeId nodeTauL = tree->getSCById(flatZoneID);

        for (NodeId nodeNL: F.getAdjacentNodes()) { //todos os vizinho de flatzone (folha) L com level > "a" 
            if( (isMaxtree && tree->getLevelById(nodeNL) <= newGrayLevel) || (!isMaxtree &&  tree->getLevelById(nodeNL) >= newGrayLevel)) {
                //o nodeNL está no intervalo [a+1, b] =>Note que: newGrayLevel é a variavel "b" do paper
                F.addNodesOfPath(tree, nodeNL, nodeTauL); 
            } 
            else { 
                //O level(nodeNL) > b
                NodeId nodeSubtree = nodeNL;
                for (NodeId n : tree->getNodesOfPathToRootById(nodeNL)) {
                    if ( (isMaxtree && newGrayLevel > tree->getLevelById(n)) || (!isMaxtree && newGrayLevel < tree->getLevelById(n))) {
                        break;
                    }
                    nodeSubtree = n; //nodeSubtree é a raiz da subarvore antes de atingir o nivel b
                }
                
                if (tree->getLevelById(nodeSubtree) == newGrayLevel) {
                    // se nodeSubtree tiver level = b, então ela entra em F[\lambda]
                    F.addNodesOfPath(tree, nodeSubtree, nodeTauL); //F_lambda
                } 
                else {
                    // caso raro: se o pai de nodeSubtree for diferente de tauL
                    if(tree->getParentById(nodeSubtree) && tree->getParentById(nodeSubtree) != nodeTauL ){ 
                        F.addNodesOfPath(tree, tree->getParentById(nodeSubtree), nodeTauL); // Adiciona os nodes do caminho até tauL
                    }else{
                        F.addNodesInFb(nodeSubtree); //F_{lambda} > b
                    }
                }
            }
            
        }
    }

    /// @brief Realiza poda (prunning) de um nó e seus descendentes, atualizando SCs.
    void prunning(ComponentTreeFZPtr tree, NodeFZ node) {
        assert(node && "Erro: node is nullptr");
        assert(node.getParent() && "Erro: node é a raiz");

        if (node != tree->getRoot()) {
            NodeFZ parent = node.getParent();
            
            FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
            //tratamento dos representantes de flatzones que já foram fundidos 
            std::vector<int>& repsParent = parent.getRepCNPs();        
            //1 normaliza os representantes
            for (int& r : repsParent) r = graph->findRepresentative(r); 
            //2 remove os repetidos
            std::sort(repsParent.begin(), repsParent.end()); 
            repsParent.erase(std::unique(repsParent.begin(), repsParent.end()), repsParent.end()); 
            
            FastQueue<NodeFZ> s;
            s.push(node);
            while (!s.empty()) {
                NodeFZ child = s.pop();
                for (NodeFZ n : child.getChildren()) {
                    s.push(n);
                }
                for(int repBase: child.getRepCNPs()) {
                    tree->setSC(repBase, parent); // atualiza o SC para apontar para o pai
                    
                    //Normaliza o representante base
                    int winner = graph->findRepresentative(repBase); // normaliza o representante
                        
                    //Remove todos que colapsaram no mesmo winner 
                    std::erase_if(repsParent, [&](int r){ return r == winner; });

                    //Adiciona o vencedor 1x no final
                    repsParent.push_back(winner);
                }
                tree->removeChildById(child.getParent(), child, true);
                //tree->releaseNode(child);
            }

        }

    }

    /// @brief Funde um nó com seu pai, movendo filhos e ajustando representantes.
    void mergeWithParent(ComponentTreeFZPtr tree, NodeFZ node) {
        if (!node || node.getParent().getIndex() == -1) return;
        FlatZonesGraphPtr& flatzoneGraph = tree->getFlatZonesGraph();
        NodeFZ parent = node.getParent();
        
        // Retira 'node' da cadeia do pai e move todos os seus filhos para o pai
        parent.removeChild(node, false);
        parent.spliceChildren(node);   // reparenta todos os filhos de 'node' para 'parent'
        
        //tratamento dos representantes de flatzones que já foram fundidos 
        std::vector<int>& repsParent = parent.getRepCNPs();        

        //1 normaliza os representantes
        for (int& r : repsParent) r = flatzoneGraph->findRepresentative(r); 

        //2 remove os repetidos
        std::sort(repsParent.begin(), repsParent.end()); 
        repsParent.erase(std::unique(repsParent.begin(), repsParent.end()), repsParent.end()); 

        for(int repBase: node.getRepCNPs()) {
            tree->setSC(repBase, parent); // atualiza o SC para apontar para o pai
            
            //Normaliza o representante base
            int winner = flatzoneGraph->findRepresentative(repBase); // normaliza o representante
                
            //Remove todos que colapsaram no mesmo winner 
            std::erase_if(repsParent, [&](int r){ return r == winner; });

            //Adiciona o vencedor 1x no final
            repsParent.push_back(winner);
        }
        tree->releaseNode(node);
    }

    /// @brief Funde o nó SC associado a um representante com seu pai.
    void mergeWithParent(ComponentTreeFZPtr tree, int repBase) {
        NodeFZ node = tree->getSC(repBase);
        if (!node || node.getParent().getIndex() == -1) return;
        if(node.getNumFlatzone() == 1) {
            mergeWithParent(tree, node);
            return; 
        }
        
        node.removeFlatzone(repBase); 
        node.setArea(node.getArea() - areaFZsRemoved); //a area é o unico atributo que mantido pela arena

        FlatZonesGraphPtr& flatzoneGraph = tree->getFlatZonesGraph();
        NodeFZ parent = node.getParent();

        //tratamento dos representantes de flatzones que já foram fundidos 
        std::vector<int>& repsParent = parent.getRepCNPs();        

        //1 normaliza os representantes
        for (int& r : repsParent) {
            r = flatzoneGraph->findRepresentative(r); 
     //       tree->setSC(r, parent); // atualiza o SC para apontar para o pai
        }

        //2 remove os repetidos
        std::sort(repsParent.begin(), repsParent.end()); 
        repsParent.erase(std::unique(repsParent.begin(), repsParent.end()), repsParent.end()); 

        tree->setSC(repBase, parent); // atualiza o SC para apontar para o pai
            
        //Normaliza o representante base
        int winner = flatzoneGraph->findRepresentative(repBase); // normaliza o representante
                
        //Remove todos que colapsaram no mesmo winner 
        std::erase_if(repsParent, [&](int r){ return r == winner; });

        //Adiciona o vencedor 1x no final
        repsParent.push_back(winner);
        
        
    }
  
};


#endif