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
void NodeFZ::addCNPsOfDisjointFlatzones(std::vector<int>& repsFlatZones) {
    NodeFZ::addCNPsOfDisjointFlatzones(repsFlatZones, this->getIndex(), this->tree);
}

// --- FlatZones: remove uma FZ do nó
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::removeFlatzone(int repFlatZone) {
    NodeFZ::removeFlatzone(repFlatZone, this->getIndex(), this->tree);
}

// --- FlatZones: conecta reps adicionais à FZ base (uma base) 
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(int repFlatZone) { //usado na versão ByLeaf
    NodeFZ::addCNPsToConnectedFlatzone(repFlatZone, this->getIndex(), this->tree);
}

// --- FlatZones: conecta reps adicionais à FZ base (várias bases)
template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(const std::vector<int>& repBases, int repBaseWinner) { //usado na versão BySubtree
    NodeFZ::addCNPsToConnectedFlatzone(repBases, repBaseWinner, this->getIndex(), this->tree);
}
