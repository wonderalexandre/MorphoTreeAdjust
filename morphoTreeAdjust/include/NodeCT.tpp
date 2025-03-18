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


//template <>
//template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
template <typename CNPsType>
int NodeCT<CNPsType>::getNumFlatzone() {
    return this->cnps.size();  //TODO: Esse método só está correto para FlatZones type
}


template <typename CNPsType>
int NodeCT<CNPsType>::getNumCNPs() const {
    if constexpr (std::is_same<CNPsType, Pixels>::value) {
        return this->cnps.size();  // Retorna diretamente o número de pixels
    } else {
        int count = 0;
        for (const auto& [id, flatzone] : this->cnps) {
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
void NodeFZ::addCNPsOfDisjointFlatzone(FlatZone&& flatZone, ComponentTreeFZ* tree) {
    this->cnps[flatZone.front()] = std::move(flatZone);
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsOfDisjointFlatzones(FlatZones&& flatZones, ComponentTreeFZ* tree) {
    for (auto& [id, flatzone] : flatZones) {  
        this->cnps[id] = std::move(flatzone); 

        for (int p : this->cnps[id]) {
            tree->setSC(p, this);
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
    } 
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void NodeFZ::addCNPsToConnectedFlatzone(FlatZone&& flatZone, ComponentTreeFZ* tree) {
    assert(!flatZone.empty() && "Erro: flatZone passada está vazia!");
    int flatZoneID = flatZone.front();    
    assert(tree->flatzoneGraph[flatZoneID] != nullptr && "Erro: flatZone não está registrada no grafo!");
    assert(!tree->flatzoneGraph[flatZoneID]->empty() && "Erro: flatZone não tem vizinhos registrados no grafo!");

    AdjacentFlatzones* flatZoneNeighbors = tree->flatzoneGraph[flatZoneID];
    std::list<int> flatzonesToMergeList; //flatzones que serão fundidas e removidas
    int unifiedFlatzoneID = std::numeric_limits<int>::max();
    for (int neighborID : *flatZoneNeighbors) {
        if (tree->getSC(neighborID) == this) {
            unifiedFlatzoneID = std::min(unifiedFlatzoneID, neighborID);
        }
    }
    std::list<int>* unifiedFlatzone = &tree->getFlatzoneByID(unifiedFlatzoneID);
    

    flatZoneNeighbors->erase(unifiedFlatzoneID);
    tree->flatzoneGraph[flatZoneID]->erase(unifiedFlatzoneID);
    tree->flatzoneGraph[unifiedFlatzoneID]->erase(flatZoneID);
    
    for (int flatzonMergedID : *flatZoneNeighbors) {
        if (tree->getSC(flatzonMergedID) == this) {
            flatzonesToMergeList.push_back(flatzonMergedID); 
        }
    }
    
    for (int flatzonMergedID : flatzonesToMergeList) { //vizinhos de flatzoneID
        for (int neighborID : *tree->flatzoneGraph[flatzonMergedID]) {
            if ( flatZoneID != neighborID) {
                tree->flatzoneGraph[unifiedFlatzoneID]->insert(neighborID);
                tree->flatzoneGraph[neighborID]->insert(unifiedFlatzoneID); 
            }
            tree->flatzoneGraph[neighborID]->erase(flatzonMergedID);
        }
        delete tree->flatzoneGraph[flatzonMergedID];
        tree->flatzoneGraph[flatzonMergedID] = nullptr;
    }

    // Remover `flatZone` do grafo
    for (int neighborID : *flatZoneNeighbors) {
        tree->flatzoneGraph[neighborID]->erase(flatZoneID);
        tree->flatzoneGraph[unifiedFlatzoneID]->insert(neighborID);
        tree->flatzoneGraph[neighborID]->insert(unifiedFlatzoneID);  
    }
    delete tree->flatzoneGraph[flatZoneID];
    tree->flatzoneGraph[flatZoneID] = nullptr; 
    

    // Atualizar SC para os novos pixels antes de modificar `flatZone`
    for (int p : flatZone) {
        tree->setSC(p, this);
    }

    // Iterar diretamente sobre `flatzonesToMergeList` para fundir seus CNPs na `unifiedFlatzone`
    for (int fzID : flatzonesToMergeList) {
        auto it = this->cnps.find(fzID);
        unifiedFlatzone->splice(unifiedFlatzone->end(), it->second);  // Fundir na unifiedFlatzone
        this->cnps.erase(it);  // Remove do unordered_map
    }

    // Fundir `flatZone` na `unifiedFlatzone`
    if(unifiedFlatzoneID < flatZoneID){
        unifiedFlatzone->splice(unifiedFlatzone->end(), flatZone);
    }else{
        //para manter a propriedade do id da flatzone ser o menor pixel
        for (int neighborID : *tree->flatzoneGraph[unifiedFlatzoneID]) {
            tree->flatzoneGraph[neighborID]->erase(unifiedFlatzoneID);
            tree->flatzoneGraph[neighborID]->insert(flatZoneID);
        }
        tree->flatzoneGraph[flatZoneID] = tree->flatzoneGraph[unifiedFlatzoneID];
        tree->flatzoneGraph[unifiedFlatzoneID] = nullptr;        
        
        flatZone.splice(flatZone.end(), *unifiedFlatzone);
        
        this->cnps[flatZoneID] = std::move(flatZone);
        this->cnps.erase(unifiedFlatzoneID);

        unifiedFlatzoneID = flatZoneID;
        unifiedFlatzone = &tree->getFlatzoneByID(unifiedFlatzoneID);

    }


    assert([&]() {
        int minPixel = *std::min_element(unifiedFlatzone->begin(), unifiedFlatzone->end());
        return minPixel == unifiedFlatzoneID && unifiedFlatzoneID == unifiedFlatzone->front();
    }() && "ERRO: O menor pixel da flatzone unificada não é o seu ID!");
    
    assert([&]() {
        
        
        if (unifiedFlatzone->empty()) {
            std::cerr << "ERRO: unifiedFlatzone está vazia após a fusão!" << std::endl;
            return false;
        }
        
        if (tree->flatzoneGraph[unifiedFlatzoneID] == nullptr) {
            std::cerr << "ERRO: unifiedFlatzone não está registrada no grafo!" << std::endl;
            return false;
        }

        for (int neighborID : *tree->flatzoneGraph[unifiedFlatzoneID]) {

            if (tree->flatzoneGraph[neighborID] == nullptr) {
                std::cerr << "ERRO: Conexão assimétrica entre unifiedFlatzone e seu vizinho!" << std::endl;
                return false;
            }
            
            const FlatZone& neighborFlatzone = tree->getFlatzoneByID(neighborID);

            if (neighborFlatzone.empty()) {
                std::cerr << "neighborID: " << neighborID << std::endl;
                std::cerr << "ERRO: Flatzone vizinha de unifiedFlatzone está vazia APOS fusão!" << std::endl;
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
        int minPixel = std::numeric_limits<int>::max();
        for (const auto& [idFlatZone, flatZone] : this->cnps) {
            minPixel = std::min(minPixel, idFlatZone);
        }
        return minPixel;
    }
}

template <typename CNPsType>
int NodeCT<CNPsType>::computerNumDescendants() {
    int numDescendants = 0;
    for(NodeCT<CNPsType>* desc: this->getIteratorBreadthFirstTraversal()){
        if(desc != this){
            numDescendants++;
        }
    }
    return numDescendants;
}
