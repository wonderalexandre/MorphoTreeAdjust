#include <iterator>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"
#include "../include/ComponentTreeAdjustment.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_SUBTREE_H
#define COMPONENT_TREE_ADJUSTMENT_SUBTREE_H

/**
 * @brief Coletor auxiliar para união de partes próprias (proper parts) em ajustes por subárvore.
 *
 * Dado o nó raiz da subárvore na árvore complementar, este coletor mantém:
 *  - o conjunto de representantes de flat-zones envolvidos (repFlatzones);
 *  - os nós da árvore corrente que possuem flat-zones a remover;
 *  - marcações rápidas de remoção (isNodesToBeRemoved) por NodeId;
 *  - o nó τ⋆ (nodeTauStar) que delimita o intervalo [g(p), f(τ⋆)] do ajuste;
 *  - o representante da FZ de τ⋆ e o vencedor final para união das FZs (fzWinner).
 */
class ProperPartsMergeCollector {
protected:
    bool isMaxtree; // Se true, percorre de forma decrescente
    
    NodeId nodeTauStar; //nodeTauStar é o nó correspondente da sub-arvore a ser podada com maior (ou menor, para min-tree) intensidade
    int repFZTauStar;
    
    std::vector<int> repFlatzones; // Lista de flatzones que serão unidas
    
    std::vector<NodeId> nodesWithFZToBeRemoved; // nodes que possuem flatzones a serem removidas
    std::vector<uint8_t> isNodesToBeRemoved; // vetor de booleanos que indica se o nó deve ser removido ou não
    
    int fzWinner;
public:
    /**
     * @brief Constrói o coletor definindo direção (max/min) e dimensão da marcação.
     * @param isMaxtree Verdadeiro se for para max-tree (ordem decrescente de níveis).
     * @param numNodes Número máximo de nós para dimensionar a estrutura de marcação.
     */
    ProperPartsMergeCollector(bool isMaxtree, int numNodes) : isMaxtree(isMaxtree), isNodesToBeRemoved(numNodes, 0) { }
    
    /**
     * @brief Acessa os representantes de flat-zones coletados.
     * @return Referência mutável ao vetor de representantes.
     */
    std::vector<int>& getRepsFZ() {
        return repFlatzones;
    }

    /**
     * @brief Retorna o NodeId de τ⋆ (raiz da subárvore-limite no intervalo do ajuste).
     */
    NodeId getNodeTauStar() {
        return nodeTauStar;
    }
   
    /**
     * @brief Retorna o representante da flat-zone de τ⋆.
     */
    int getRepFZTauStar() {
        return repFZTauStar;
    }

    /**
     * @brief Lista nós marcados para remoção. O(nós marcados).
     * @note Método para debug/inspeção (varre a lista). Em produção, prefira `isRemoved`.
     */
    std::vector<NodeId> getNodesToBeRemoved() {
        std::vector<NodeId> nodesToBeRemovedTmp;
        for(NodeId node: nodesWithFZToBeRemoved) {
            if (isNodesToBeRemoved[node]) {
                nodesToBeRemovedTmp.push_back(node);
            }
        }
        return nodesToBeRemovedTmp;
    }

    /**
     * @brief Une todas as FZs coletadas a uma base conectada dentro de `nodeUnion`.
     * @param nodeUnion Nó que receberá a FZ conectada (vencedor).
     * @param tree Árvore na qual a união ocorre.
     */
    void addCNPsToConnectedFlatzone(NodeFZ nodeUnion, ComponentTreeFZPtr tree) {
       nodeUnion.addCNPsToConnectedFlatzone(repFlatzones, fzWinner, tree);
       this->removeFlatzones(tree);
    }

    /**
     * @brief Remove das FZs dos nós coletados aquelas que foram agregadas ao vencedor.
     * @param tree Árvore-alvo.
     */
    void removeFlatzones(ComponentTreeFZPtr tree){
        for(size_t i = 0; i < repFlatzones.size(); ++i) {
            NodeFZ node = tree->proxy( nodesWithFZToBeRemoved[i] );
            int repFZ = repFlatzones[i];
            node.removeFlatzone(repFZ);  
            if(node.getNumFlatzone() == 0){ 
                isNodesToBeRemoved[node] = 1; // Marca o nó para remoção
            }
        }
    }

    /**
     * @brief Testa se um nó foi completamente removido após a união.
     */
    bool isRemoved(NodeId nodeId){
        return isNodesToBeRemoved[nodeId];
    }
    
    /**
     * @brief Limpa coleções e reconfigura direção (max/min) para um novo ajuste.
     * @param isMaxtree Verdadeiro para max-tree.
     */
    void resetCollections(bool isMaxtree) {
        this->isMaxtree = isMaxtree;
        for (NodeId nodeId : nodesWithFZToBeRemoved) {
            this->isNodesToBeRemoved[nodeId] = 0; 
        }
        this->repFZTauStar = -1; 
        this->nodesWithFZToBeRemoved.clear();
        this->repFlatzones.clear();    

        this->nodeTauStar = -1;
        this->fzWinner = -1;
    }
    
    /**
     * @brief Adiciona um nó e a FZ correspondente ao conjunto a ser unido.
     * @param tree Árvore referência para níveis.
     * @param nodeTau NodeId do nó a ser considerado.
     * @param repFZ Representante da FZ ligada a `nodeTau`.
     */
    void addNode(ComponentTreeFZPtr tree, NodeId nodeTau, int repFZ) {
        repFlatzones.push_back(repFZ);
        nodesWithFZToBeRemoved.push_back(nodeTau);
        
        if (this->nodeTauStar == -1 || ( (!isMaxtree && tree->getLevelById(nodeTau) > tree->getLevelById(nodeTauStar)) || (isMaxtree && tree->getLevelById(nodeTau) < tree->getLevelById(nodeTauStar)))) {
            this->nodeTauStar = nodeTau;
            this->repFZTauStar = repFZ; 
        }

        if(this->fzWinner == -1) {
            fzWinner = repFZ;
        }else if(fzWinner > repFZ) {
            fzWinner = repFZ;
        }
    }

};

/**
 * @brief Ajuste de árvores guiado por uma subárvore (raiz rSubtree) na árvore complementar.
 *
 * Coleta as partes próprias de rSubtree na outra árvore (min ↔ max), constrói
 * as coleções F_λ/F_{λ>b} e executa as fusões/reatamentos conforme o intervalo
 * [g(p), f(τ⋆)]. Opera simetricamente para min-tree e max-tree.
 */
class ComponentTreeAdjustmentBySubtree: public ComponentTreeAdjustment {

protected:
    ProperPartsMergeCollector properPartsCollector;

public:
    /**
     * @brief Constrói o ajustador configurando o coletor de proper parts.
     * @param mintree Ponteiro para a min-tree.
     * @param maxtree Ponteiro para a max-tree.
     */
    ComponentTreeAdjustmentBySubtree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree) : ComponentTreeAdjustment(mintree, maxtree), properPartsCollector(maxtree->isMaxtree(), std::max(mintree->getNumNodes(), maxtree->getNumNodes())) { }
      
    /**
     * @brief Atualiza a árvore `tree` a partir da subárvore enraizada em `rootSubtree` na outra árvore.
     * @param tree Árvore (min ou max) a ser atualizada.
     * @param rootSubtree Raiz da subárvore removida/mesclada na árvore complementar.
     */
    void updateTree(ComponentTreeFZPtr tree, NodeFZ node);
    
    /**
     * @brief Ajusta a min-tree após podas por subárvore na max-tree.
     * @param mintree Árvore mínima a ser ajustada.
     * @param maxtree Árvore máxima onde subárvores foram removidas/mescladas.
     * @param nodesToPruning Lista de NodeId raízes de subárvore na max-tree.
     */
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToPruning) ;
    
    /**
     * @brief Ajusta a max-tree após podas por subárvore na min-tree.
     * @param maxtree Árvore máxima a ser ajustada.
     * @param mintree Árvore mínima onde subárvores foram removidas/mescladas.
     * @param nodesToPruning Lista de NodeId raízes de subárvore na min-tree.
     */
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToPruning) ;


};


#endif
