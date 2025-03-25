#ifndef COMMONS_HPP  
#define COMMONS_HPP  


#define NDEBUG  // Remove os asserts do código
#include <cassert>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <memory>
#include <limits>

#define PRINT_LOG 1 

#define PRINT_DEBUG 0 

// Forward declaration dos templates
template <typename T>
class ComponentTree;

template <typename T>
class NodeCT;




// Definição de tipos para CNPs
using Pixels = std::list<int>;                 // Representa uma lista de pixels (CNPs)
using FlatZone = std::list<int>;                 // Representa uma flatzone
using FlatZones = std::unordered_map<int, FlatZone>;   // Representa uma lista de flatzones (CNPs ficarão separados em flatzones)

//Alias em função do tipo dos CNPs
using FlatZoneRef = std::reference_wrapper<FlatZone>;   // Representa uma flatzone
using FlatZonesRef = std::list<FlatZoneRef>;   // Representa uma lista de flatzones
using ComponentTreeFZ = ComponentTree<FlatZones>; //representa um component tree por flatzones
using ComponentTreeP = ComponentTree<Pixels>; //representa um component tree sem tratamento de flatzones
using NodeFZ = NodeCT<FlatZones>; //representa um node com separação dos cnps em flatzones
using NodeP = NodeCT<Pixels>; //representa um node sem tratamento de flatzones


template <typename CNPsType>
using NodeCTPtr = std::shared_ptr<NodeCT<CNPsType>>;
using NodeFZPtr = std::shared_ptr<NodeCT<FlatZones>>;
using NodePPtr  = std::shared_ptr<NodeCT<Pixels>>;

template <typename CNPsType>
using ComponentTreePtr = std::shared_ptr<ComponentTree<CNPsType>>; 
using ComponentTreeFZPtr = std::shared_ptr<ComponentTreeFZ>; //representa um component tree por flatzones
using ComponentTreePPtr = std::shared_ptr<ComponentTreeP>; //representa um component tree sem tratamento de flatzones


using AdjacentFlatzones = std::unordered_set<int>;
using FlatzoneGraph = std::vector<std::unique_ptr<AdjacentFlatzones>>;


struct FlatZoneNode {
    NodeFZPtr node = nullptr;
    FlatZone* flatzone = nullptr;
    int idFlatZone;

    FlatZoneNode(){}

    // Construtor para mover a FlatZone
    FlatZoneNode(NodeFZPtr n, FlatZone& fz) : node(n), flatzone(&fz), idFlatZone(fz.front()) {} 
    

};



#endif 
