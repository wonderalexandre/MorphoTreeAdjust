
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

// Construtor padrão
template <typename CNPsType>
NodeCT<CNPsType>::NodeCT() : children(), cnps() {}

// Construtor parametrizado
template <typename CNPsType>
NodeCT<CNPsType>::NodeCT(int index, NodeCT* parent, int threshold1, int threshold2) {
    this->index = index;
    this->parent = parent;
    this->threshold2 = threshold2;
    this->threshold1 = threshold1;
    this->areaCC = 0;
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, Pixels>::value, int>>
void NodeP::addCNPs(int p){
    this->cnps.push_back(p);
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
int NodeFZ::getNumFlatzone() {
    return this->cnps.size();  // Número de flatzones armazenadas
}


template <typename CNPsType>
int NodeCT<CNPsType>::getNumCNPs() const {
    if constexpr (std::is_same<CNPsType, Pixels>::value) {
        return this->cnps.size();  // Retorna diretamente o número de pixels
    } else {
        int count = 0;
        for (const auto& flatzone : this->cnps) {
            count += flatzone.size();  // Soma os pixels de todas as flatzones
        }
        return count;
    }
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZones NodeFZ::moveCNPsByFlatZone() {
    return std::exchange(this->cnps, {});  // Move os CNPs e deixa `cnps` vazio
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZones& NodeFZ::getCNPsByFlatZone() {
    return this->cnps;
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsOfDisjointFlatzone(FlatZone&& flatZone, ComponentTreeFZ* tree) {
    // Adiciona a nova flatzone à lista de flatzones do nó
    this->cnps.push_back(std::move(flatZone));
    FlatZone& newFlatzone = this->cnps.back();  // Obtém referência à nova flatzone

    // Atualiza o mapeamento de pixels para flatzone
    for (int p : newFlatzone) { 
        tree->updatePixelToFlatzone(p, &newFlatzone);
    }
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsOfDisjointFlatzones(FlatZones&& flatZones, ComponentTreeFZ* tree) {
    size_t oldSize = this->cnps.size();

    // Movemos todas as flatzones para `cnps` sem copiar elementos
    this->cnps.splice(this->cnps.end(), std::move(flatZones));

    // Atualiza `pixelToFlatzone` e `tree->setSC(p, this)`
    auto it = this->cnps.begin();
    std::advance(it, oldSize);  // Move iterador para o início das novas flatzones

    for (; it != this->cnps.end(); ++it) {
        FlatZone& flatzone = *it;  // `flatzone` é `std::list<int>`

        for (int p : flatzone) {
            tree->updatePixelToFlatzone(p, &flatzone);
            tree->setSC(p, this);
        }
    }
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::removeFlatzone(FlatZone& flatzone) {
    // Encontrar e remover a flatzone da lista cnps
    for (auto it = cnps.begin(); it != cnps.end(); ++it) {
        if (&(*it) == &flatzone) {  // Verifica se a referência coincide
            cnps.erase(it);
            return;
        }
    }

    assert(false && "Erro: A flatzone a ser removida não está presente na lista de flatzones do nó!");
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(FlatZone&& flatZone, ComponentTreeFZ* tree) {
    assert(!flatZone.empty() && "Erro: flatZone passada está vazia!");
    
    std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual> flatzonesToMergeSet;
    std::list<int>* targetFlatzone = nullptr;

    assert(tree->flatzoneGraph.find(flatZone) != tree->flatzoneGraph.end() && "Erro: flatZone não está registrada no grafo!");
    assert(!tree->flatzoneGraph[flatZone].empty() && "Erro: flatZone não tem vizinhos registrados no grafo!");

    // Encontrar as flatzones vizinhas usando o grafo
    for (const auto& neighborRef : tree->flatzoneGraph[flatZone]) {
        FlatZone& neighborFlatzone = neighborRef.get();

        if (tree->getSC(neighborFlatzone.front()) == this) {
            if (flatzonesToMergeSet.insert(neighborRef).second) {
                if (!targetFlatzone || neighborFlatzone.size() > targetFlatzone->size()) {
                    targetFlatzone = &neighborFlatzone;
                }
            }
        }
    }
    
    assert(targetFlatzone && "Erro: nenhuma flatzone válida para fusão!");

    // Atualizar SC para os novos pixels antes de modificar `flatZone`
    for (int p : flatZone) {
        tree->setSC(p, this);
    }

    // === Atualizar o grafo ANTES da fusão ===
    std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual> neighborsToReinsert;
    for (const auto& flatzoneRef : flatzonesToMergeSet) {
        for (const auto& neighborRef : tree->flatzoneGraph[flatzoneRef]) {
            if (flatzonesToMergeSet.find(neighborRef) == flatzonesToMergeSet.end()) {
                neighborsToReinsert.insert(neighborRef);
            }
            tree->flatzoneGraph[neighborRef].erase(flatzoneRef);  // Remover a referência do vizinho
        }
        // Remover a flatzone do grafo
        tree->flatzoneGraph.erase(flatzoneRef);
    }
    // Remover `flatZone` do grafo
    for (const auto& neighborRef : tree->flatzoneGraph[flatZone]) {
        neighborsToReinsert.insert(neighborRef);
        tree->flatzoneGraph[neighborRef].erase(flatZone);
    }
    tree->flatzoneGraph.erase(flatZone);
    // ==========================================
    
    // Guardamos o tamanho antes da fusão, pois `targetFlatzone` pode crescer após os splices
    size_t targetFlatzoneOriginalSize = targetFlatzone->size();

    // Fundir todas as flatzones na `targetFlatzone`
    for (auto it = this->cnps.begin(); it != this->cnps.end(); ) {
        if (flatzonesToMergeSet.count(*it) > 0 && &(*it) != targetFlatzone) {
            targetFlatzone->splice(targetFlatzone->end(), *it);
            it = this->cnps.erase(it);  
        } else {
            ++it;
        }
    }

    // Fundir `flatZone` na `targetFlatzone`
    targetFlatzone->splice(targetFlatzone->end(), flatZone);

    // Atualizar `pixelToFlatzone` apenas para os novos pixels adicionados
    size_t newPixelsToMap = targetFlatzone->size() - targetFlatzoneOriginalSize;
    for (auto it = targetFlatzone->rbegin(); newPixelsToMap > 0; ++it, --newPixelsToMap) {
        tree->updatePixelToFlatzone(*it, targetFlatzone);
    }

    assert([&]() { 
        for (int p : *targetFlatzone) {
            if (tree->getFlatzoneRef(p) != *targetFlatzone) {
                return false;
            }
        }
        return true;
    }() && "Erro: mapeamento pixelToFlatzone está errado para os pixels de targetFlatzone");
    

    // === Atualizar o grafo APÓS a fusão ===
    //Adicionar ao grafo a flatzone fundida e suas arestas
    FlatZoneRef targetFlatzoneRef = *targetFlatzone;
    tree->flatzoneGraph[targetFlatzoneRef] = {};  // Criar a entrada para a nova flatzone fundida

    for (const auto& neighborRef : neighborsToReinsert) {
        if(!neighborRef.get().empty()){
            tree->flatzoneGraph[targetFlatzoneRef].insert(neighborRef);
            tree->flatzoneGraph[neighborRef].insert(targetFlatzoneRef);  
        }
    }
    // ==========================================

    assert([&]() {
        if (targetFlatzone->empty()) {
            std::cerr << "ERRO: targetFlatzone está vazia após a fusão!" << std::endl;
            return false;
        }
        
        if (tree->flatzoneGraph.find(targetFlatzoneRef) == tree->flatzoneGraph.end()) {
            std::cerr << "ERRO: targetFlatzone não está registrada no grafo!" << std::endl;
            return false;
        }

        for (const auto& neighborRef : tree->flatzoneGraph[targetFlatzoneRef]) {
            const FlatZone& neighborFlatzone = neighborRef.get();

            if (neighborFlatzone.empty()) {
                std::cerr << "ERRO: Flatzone vizinha de targetFlatzone está vazia!" << std::endl;
                return false;
            }

            if (tree->flatzoneGraph[neighborRef].find(targetFlatzoneRef) == tree->flatzoneGraph[neighborRef].end()) {
                std::cerr << "ERRO: Conexão assimétrica entre targetFlatzone e seu vizinho!" << std::endl;
                return false;
            }
        }

        return true;
    }() && "Erro: Grafo de flatzones inconsistente após a fusão!");
}

template <typename CNPsType>
void NodeCT<CNPsType>::addChild(NodeCT<CNPsType>* child) {
    this->children.push_back(child);
} 

template <typename CNPsType>
bool NodeCT<CNPsType>::isChild(NodeCT<CNPsType>* child) const {
    return std::find(this->children.begin(), this->children.end(), child) != children.end();
}

template <typename CNPsType>
int NodeCT<CNPsType>::getIndex() const { 
    return this->index; 
}

template <typename CNPsType>
int NodeCT<CNPsType>::getThreshold1() const { 
    return this->threshold1; 
}

template <typename CNPsType>
int NodeCT<CNPsType>::getThreshold2() const { 
    return this->threshold2; 
}

template <typename CNPsType>
int NodeCT<CNPsType>::getLevel() const { 
    return this->threshold2; 
}

template <typename CNPsType>
void NodeCT<CNPsType>::setLevel(int level) { 
    this->threshold2 = level; 
}

template <typename CNPsType>
void NodeCT<CNPsType>::setArea(long int area) { 
    this->areaCC = area; 
}

template <typename CNPsType>
long int NodeCT<CNPsType>::getArea() const { 
    return this->areaCC; 
}

template <typename CNPsType>
NodeCT<CNPsType>* NodeCT<CNPsType>::getParent() {  
    return this->parent; 
}

template <typename CNPsType>
void NodeCT<CNPsType>::setParent(NodeCT<CNPsType>* parent) { 
    this->parent = parent; 
}

template <typename CNPsType>
bool NodeCT<CNPsType>::isLeaf() const { 
    return this->children.empty(); 
}

template <typename CNPsType>
std::list<NodeCT<CNPsType>*>& NodeCT<CNPsType>::getChildren() {  
    return this->children; 
}

template <typename CNPsType>
int NodeCT<CNPsType>::getNumSiblings() const {
    return (this->parent != nullptr) ? this->parent->getChildren().size() : 0;
}

template <typename CNPsType>
int NodeCT<CNPsType>::getRepresentativeCNPs() const{
    if constexpr (std::is_same<CNPsType, Pixels>::value) {
        return this->cnps.front();
    }else{
        return this->cnps.front().front();
    }
}