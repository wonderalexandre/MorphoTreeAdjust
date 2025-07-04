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
NodeCT<CNPsType>::NodeCT(int index, NodeCTPtr<CNPsType> parent, int threshold1, int threshold2) {
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


//template <>
//template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
template <typename CNPsType>
int NodeCT<CNPsType>::getNumFlatzone() {
    return this->cnps.size();  //TODO: Esse método só está correto para FlatZones type
}


template <typename CNPsType>
int NodeCT<CNPsType>::getNumCNPs() {
    if constexpr (std::is_same<CNPsType, Pixels>::value) {
        return this->cnps.size();  // Retorna diretamente o número de pixels
    } else {
       // if(this->numCNPs != -1) return this->numCNPs;
        int numCNPs = 0; //recomputando o cache contendo o numero de cnps
        for (const auto& [id, flatzone] : this->cnps) {
            numCNPs += flatzone.size();  // Soma os pixels de todas as flatzones
        }
        return numCNPs;
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
    assert([&]() {
        for (const auto& [idFlatZone, flatzone] : this->cnps) {
            if(flatzone.empty()){
                std::cerr << "\n\nIndex: " << index << std::endl;
                std::cerr << "idFlatZone: " << idFlatZone << std::endl;
                std::cerr << "|flatzone|: " << flatzone.size() << std::endl;
                std::cerr << "cnps.size: " << cnps.size() << std::endl;
                return false;
            }
        }
        return true;
    }() && "Erro: Existem flatzones vazias");

    return this->cnps;
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZone& NodeFZ::getFlatZone(int idFlatZone) {
    auto it = this->cnps.find(idFlatZone);
    if (it != this->cnps.end()) {
        return it->second;
    }
    std::string msgErro = "Erro: Nenhuma FlatZone encontrada com idFlatZone = " + std::to_string(idFlatZone) + "\n\n";
    msgErro += "Index: " + std::to_string(index) + "\n";
    msgErro += "NumCNPs: " + std::to_string(this->getNumCNPs()) + "\n";
    msgErro += "cnps.size: " + std::to_string(cnps.size()) + "\n";
    for (const auto& [id, flatzone] : this->cnps) {
        if(!flatzone.empty()){
            msgErro += "\tidFlatZone: " + std::to_string(id) + "\n";
            msgErro += "\t|FlatZone|: " + std::to_string(flatzone.size()) + "\n";
            for(int p: flatzone){
                if(p == idFlatZone){
                    msgErro += "\tO pixel " + std::to_string(p) + " está na flatzone mas o id da flatzone é: " + std::to_string(id) + "\n";
                    break;
                }
            }
        }
    }
    throw std::runtime_error(msgErro);
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsOfDisjointFlatzone(FlatZone&& flatZone, ComponentTreeFZPtr tree, int capacity) {
    if(capacity != -1){
        this->cnps.reserve(capacity);
    }
    int id = flatZone.front();
    this->cnps[id] = std::move(flatZone);
    if (tree) {
        auto ptr = this->shared_from_this();
        for (int p : this->cnps[id]) {
            tree->setSC(p, ptr);
        }
    }
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsOfDisjointFlatzones(FlatZones&& flatZones, ComponentTreeFZPtr tree) {
    for (auto& [id, flatzone] : flatZones) {  
        this->cnps[id] = std::move(flatzone); 
        auto ptr = this->shared_from_this();
        for (int p : this->cnps[id]) {
            tree->setSC(p, ptr);
        }
    }
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::removeFlatzone(int idFlatZone) {
    auto it = cnps.find(idFlatZone);

    assert([&]() {
        if (it == cnps.end()) {
            std::cerr << "\n\nIndex: " << index << std::endl;
            std::cerr << "idFlatZone: " << idFlatZone << std::endl;
            std::cerr << "cnps.size: " << cnps.size() << std::endl;
            return false;
        }
        return true;
    }() && "Erro: A flatzone a ser removida não está presente no nó!");

    if (it != cnps.end()) {
        cnps.erase(it);  
       // this->numCNPs = -1; //o cache contendo o numero de cnps será recomputado
    } 
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(FlatZone&& flatZone, ComponentTreeFZPtr tree) {
   assert(!flatZone.empty() && "Erro: flatZone passada está vazia!");
   //this->numCNPs = -1; //o cache contendo o numero de cnps será recomputado
   int flatZoneID = flatZone.front();    

   FlatZonesGraphPtr& graph = tree->getFlatZonesGraph();
   auto [unifiedFlatzoneID, flatzonesToMergeList] = graph->mergeConnectedFlatzone(flatZoneID, this->shared_from_this(), tree) ;
   FlatZone& unifiedFlatzone = this->getFlatZone(unifiedFlatzoneID); //tree->getSC(unifiedFlatzoneID)->getFlatZone(unifiedFlatzoneID); 

    // Atualizar SC para os novos pixels antes de modificar `flatZone`
    auto ptr = this->shared_from_this();
    for (int p : flatZone) {
        tree->setSC(p, ptr);
    }

    // Iterar diretamente sobre `flatzonesToMergeList` para fundir seus CNPs na `unifiedFlatzone`
    for (const int fzID : flatzonesToMergeList){
        auto it = this->cnps.find(fzID);
        unifiedFlatzone.splice(unifiedFlatzone.end(), it->second);  // Fundir na unifiedFlatzone
        this->cnps.erase(it);  // Remove do unordered_map
    }

    // Fundir `flatZone` na `unifiedFlatzone`
    if(unifiedFlatzoneID < flatZoneID){
        unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone);
    }else{        
        graph->remapFlatzoneIDInGraph(unifiedFlatzoneID, flatZoneID);

        flatZone.splice(flatZone.end(), unifiedFlatzone);
        this->cnps[flatZoneID] = std::move(flatZone);
        this->cnps.erase(unifiedFlatzoneID);

    }   

}



template <typename CNPsType>
void NodeCT<CNPsType>::addChild(NodeCTPtr<CNPsType> child) {
    this->children.push_back(child);
} 

template <typename CNPsType>
bool NodeCT<CNPsType>::isChild(NodeCTPtr<CNPsType> child) const {
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
NodeCTPtr<CNPsType> NodeCT<CNPsType>::getParent() {  
    return this->parent; 
}

template <typename CNPsType>
void NodeCT<CNPsType>::setParent(NodeCTPtr<CNPsType> parent) { 
    this->parent = parent; 
}

template <typename CNPsType>
bool NodeCT<CNPsType>::isLeaf() const { 
    return this->children.empty(); 
}

template <typename CNPsType>
std::list<NodeCTPtr<CNPsType>>& NodeCT<CNPsType>::getChildren() {  
    return this->children; 
}

template <typename CNPsType>
int NodeCT<CNPsType>::getNumSiblings() const {
    return (this->parent != nullptr) ? this->parent->getChildren().size() : 0;
}


template <typename CNPsType>
int NodeCT<CNPsType>::computerNumDescendants() {
    int numDescendants = 0;
    for(NodeCTPtr<CNPsType> desc: this->getIteratorBreadthFirstTraversal()){
        if(desc != this->shared_from_this()){
            numDescendants++;
        }
    }
    return numDescendants;
}
