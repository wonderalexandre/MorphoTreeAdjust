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

#ifndef COMPONENT_TREE_ADJUSTMENT_FLATZONE_H
#define COMPONENT_TREE_ADJUSTMENT_FLATZONE_H

/**
 * @brief Ajuste de árvores guiado por um representante de flat-zone.
 *
 * Atualiza a árvore complementar (min ↔ max) quando a flat-zone identificada por
 * `repFlatzone` sofre remoção/elevação de nível. Constrói coleções F_λ e F_{λ>b},
 * mescla nós por níveis, integra a flat-zone ao nível alvo e reconecta subárvores.
 * 
 */
class ComponentTreeAdjustmentByFlatzone: public  ComponentTreeAdjustment {
    
public:
    /**
     * @brief Constrói o ajustador para árvores mínimo/máximo emparelhadas.
     * @param mintree Ponteiro para a min-tree.
     * @param maxtree Ponteiro para a max-tree.
     */
    ComponentTreeAdjustmentByFlatzone(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree) : ComponentTreeAdjustment(mintree, maxtree) { }
          
    /**
     * @brief Atualiza a árvore `tree` com base na flat-zone representada por `repFlatzone`.
     * @param tree Árvore (min ou max) a ser atualizada.
     * @param repFlatzone Representante da flat-zone cujo efeito é propagado.
     */
    void updateTree(ComponentTreeFZPtr tree, int repFlatzone);
    
    /**
     * @brief Ajusta a min-tree após alteração de uma flat-zone na max-tree.
     * @param mintree Árvore mínima a ser ajustada.
     * @param maxtree Árvore máxima (origem da alteração).
     * @param repFlatzone Representante da flat-zone que orienta o ajuste.
     */
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, int repFlatzone);
    
    /**
     * @brief Ajusta a max-tree após alteração de uma flat-zone na min-tree.
     * @param maxtree Árvore máxima a ser ajustada.
     * @param mintree Árvore mínima (origem da alteração).
     * @param repFlatzone Representante da flat-zone que orienta o ajuste.
     */
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, int repFlatzone);
    
    


    
};


#endif
