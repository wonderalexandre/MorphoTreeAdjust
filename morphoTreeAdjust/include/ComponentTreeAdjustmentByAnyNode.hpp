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
#include "../include/ComponentTreeAdjustmentByFlatzone.hpp"

#ifndef COMPONENT_TREE_ADJUSTMENT_ANY_NODE_H
#define COMPONENT_TREE_ADJUSTMENT_ANY_NODE_H

/**
 * @brief Ajuste de árvores de componentes a partir de qualquer nó (não apenas folhas).
 *
 * Esta classe especializa a atualização quando um nó arbitrário (com uma ou mais
 * flat-zones) é removido/podado na árvore complementar (min ↔ max). O ajuste é
 * executado propagando o efeito de cada flat-zone do nó nos níveis correspondentes
 * da outra árvore, por meio de ComponentTreeAdjustmentByFlatzone.
 *
 */
class ComponentTreeAdjustmentByAnyNode: public ComponentTreeAdjustmentByFlatzone {

    
public:
    /**
     * @brief Constrói o ajustador para árvores mínimo/máximo emparelhadas.
     * @param mintree Ponteiro para a min-tree.
     * @param maxtree Ponteiro para a max-tree.
     */
    ComponentTreeAdjustmentByAnyNode(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree) : ComponentTreeAdjustmentByFlatzone(mintree, maxtree) { }
      
    /**
     * @brief Atualiza a árvore `tree` considerando as flat-zones do nó fornecido.
     *
     * Para cada representante de flat-zone armazenado no nó, delega a
     * ComponentTreeAdjustmentByFlatzone::updateTree a atualização pontual.
     *
     * @param tree Árvore (min ou max) a ser atualizada.
     * @param node Nó arbitrário da árvore complementar que orienta o ajuste.
     */      
    void updateTree(ComponentTreeFZPtr tree, NodeFZ node);
    
    /**
     * @brief Ajusta a min-tree após remover/mesclar nós arbitrários na max-tree.
     * @param mintree Árvore mínima a ser ajustada.
     * @param maxtree Árvore máxima onde os nós foram removidos/mesclados.
     * @param nodesToPruning Lista de NodeId na max-tree a serem processados.
     */
    void adjustMinTree(ComponentTreeFZPtr mintree, ComponentTreeFZPtr maxtree, std::vector<NodeId>& nodesToPruning);
    
    /**
     * @brief Ajusta a max-tree após remover/mesclar nós arbitrários na min-tree.
     * @param maxtree Árvore máxima a ser ajustada.
     * @param mintree Árvore mínima onde os nós foram removidos/mesclados.
     * @param nodesToPruning Lista de NodeId na min-tree a serem processados.
     */
    void adjustMaxTree(ComponentTreeFZPtr maxtree, ComponentTreeFZPtr mintree, std::vector<NodeId>& nodesToPruning);
    
    


    
};


#endif
