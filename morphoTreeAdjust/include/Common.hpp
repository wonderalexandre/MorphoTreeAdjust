#ifndef COMMONS_HPP  
#define COMMONS_HPP  


#define NDEBUG  // Remove os asserts do código
#include <cassert>
#include <list>

// Forward declaration dos templates
template <typename T>
class ComponentTree;

template <typename T>
class NodeCT;

// Definição de tipos para CNPs
using Pixels = std::list<int>;                 // Representa uma lista de pixels (CNPs)
using FlatZone = std::list<int>;                 // Representa uma flatzone
using FlatZones = std::list<FlatZone>;   // Representa uma lista de flatzones (CNPs ficarão separados em flatzones)

//Alias em função do tipo dos CNPs
using FlatZoneRef = std::reference_wrapper<std::list<int>>;   // Representa uma flatzone
using FlatZonesRef = std::list<FlatZoneRef>;   // Representa uma lista de flatzones
using ComponentTreeFZ = ComponentTree<FlatZones>; //representa um component tree por flatzones
using ComponentTreeP = ComponentTree<Pixels>; //representa um component tree sem tratamento de flatzones
using NodeFZ = NodeCT<FlatZones>; //representa um node com separação dos cnps em flatzones
using NodeP = NodeCT<Pixels>; //representa um node sem tratamento de flatzones


struct ListRefHash {
    size_t operator()(const std::reference_wrapper<std::list<int>>& ref) const {
        return reinterpret_cast<size_t>(&ref.get());  // Usa o endereço como hash
    }
};

struct ListRefEqual {
    bool operator()(const std::reference_wrapper<std::list<int>>& lhs, 
                    const std::reference_wrapper<std::list<int>>& rhs) const {
        return &lhs.get() == &rhs.get();  // Compara pelos endereços
    }
};



#endif 
