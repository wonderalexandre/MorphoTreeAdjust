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

#ifndef COMPONENT_TREE_ADJUSTMENT_LEAF_H
#define COMPONENT_TREE_ADJUSTMENT_LEAF_H

/**
 * @brief Ajuste de árvores de componentes guiado por uma folha.
 *
 * Implementa o caso clássico de ajuste quando uma folha L é podada/mesclada
 * em uma árvore (min ↔ max). O algoritmo constrói coleções de nós a serem
 * fundidos por níveis (F_λ) e reata conexões de subárvores conforme a
 * mudança de nível do elemento folha na árvore complementar.
 */
template<typename Computer = DefaultAttributeComputer>
class ComponentTreeAdjustmentByLeaf: public ComponentTreeAdjustment<Computer> {
    
public:
    /**
     * @brief Constrói o ajustador para árvores mínimo/máximo emparelhadas.
     * @param mintree Ponteiro para a min-tree.
     * @param maxtree Ponteiro para a max-tree.
     */
    ComponentTreeAdjustmentByLeaf(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree)
        : ComponentTreeAdjustment<Computer>(mintree, maxtree) { }
    
    /**
     * @brief Atualiza a árvore `tree` após a remoção/mescla da folha `L_leaf` na outra árvore.
     *
     * Constrói as coleções F_λ e F_{λ>b}, mescla nós por nível e integra a
     * flat-zone da folha no nível correspondente (λ = g(p)).
     *
     * @param tree Árvore (min ou max) a ser atualizada.
     * @param L_leaf Folha removida/mesclada na árvore complementar.
     */
    void updateTree(ComponentTreeFZPtr tree, NodeFZ L_leaf);

    /**
     * @brief Ajusta a min-tree após poda/mescla de folhas na max-tree.
     * @param mintree Árvore mínima a ser ajustada.
     * @param maxtree Árvore máxima onde folhas foram removidas/mescladas.
     * @param nodesToPruning Lista de NodeId na max-tree a percorrer (pós-ordem).
     */
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToPruning) ;
    
    /**
     * @brief Ajusta a max-tree após poda/mescla de folhas na min-tree.
     * @param maxtree Árvore máxima a ser ajustada.
     * @param mintree Árvore mínima onde folhas foram removidas/mescladas.
     * @param nodesToPruning Lista de NodeId na min-tree a percorrer (pós-ordem).
     */
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToPruning) ;
  
};

#endif
