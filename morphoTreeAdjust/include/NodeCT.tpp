#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <functional> 
#include <stdexcept>  



// --- FlatZones: adiciona reps de FZ disjuntas
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsOfDisjointFlatzones(std::vector<int>& repsFlatZones, ComponentTreeFZPtr tree) {
    for (int rep: repsFlatZones) {
        tree->setSCById(rep, this->getIndex());
    }
    auto& reps = this->getRepCNPs();
    reps.insert(reps.end(), repsFlatZones.begin(), repsFlatZones.end());
}

// --- FlatZones: remove uma FZ do nó
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::removeFlatzone(int repFlatZone) {
    auto& reps = this->getRepCNPs();
    std::erase(reps, repFlatZone);
}

// --- FlatZones: conecta reps adicionais à FZ base (uma base) 
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(int repFlatZone, ComponentTreeFZPtr tree) { //usado na versão ByLeaf
    auto& reps = this->getRepCNPs();
    assert(!reps.empty() && "Erro: conjunto de FZs do nó está vazio!");
    
    FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
    int repWinner = graph->mergeAdjacentCandidatesInPlace(repFlatZone, reps);
    tree->setSCById(repWinner, this->getIndex()); // atualiza pixelToNode
}

// --- FlatZones: conecta reps adicionais à FZ base (várias bases)
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(const std::vector<int>& repBases, int repBaseWinner, ComponentTreeFZPtr tree) { //usado na versão BySubtree
    auto& reps = this->getRepCNPs();
    assert(!reps.empty() && "Erro: conjunto de FZs do nó está vazio!");
    FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
    int repWinner = graph->mergeBasesWithAdjacentCandidatesInPlace(repBases, reps, repBaseWinner);

    tree->setSCById(repWinner, this->getIndex()); // atualiza pixelToNode
}

