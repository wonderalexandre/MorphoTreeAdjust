

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


NodeCT::NodeCT(): children(), cnpsByFlatzone(){}

NodeCT::NodeCT(int index,  NodeCT* parent, int threshold1, int threshold2) {
		this->index = index;
        this->parent = parent;
        this->threshold2 = threshold2;
        this->threshold1 = threshold1;
        this->areaCC = 0;
}

int NodeCT::getNumFlatzone(){
    return this->cnpsByFlatzone.size();
}

int NodeCT::getCNP(size_t index) {
    size_t count = 0;

    for (auto flatZoneIt = cnpsByFlatzone.begin(); flatZoneIt != cnpsByFlatzone.end(); ++flatZoneIt) {
        size_t flatZoneSize = flatZoneIt->size();

        // Se o índice está dentro desta flatzone
        if (index < count + flatZoneSize) {
            size_t localIndex = index - count;

            //  Se estiver na primeira metade, percorremos do início
            if (localIndex < flatZoneSize / 2) {
                auto it = flatZoneIt->begin();
                std::advance(it, localIndex);
                return *it;
            } else {
                //  Se estiver na segunda metade, percorremos a partir do final
                auto it = flatZoneIt->end();
                std::advance(it, localIndex - flatZoneSize);
                return *it;
            }
        }

        count += flatZoneSize;  // Soma os pixels já percorridos
    }

    throw std::out_of_range("Índice fora do limite dos CNPs.");
}


int NodeCT::getNumCNPs() const{
    int count = 0;
    for (auto& flatzone : this->cnpsByFlatzone) 
        count += flatzone.size();
    return count;
}

std::list<std::list<int>> NodeCT::moveCNPsByFlatZone() {
    return std::exchange(this->cnpsByFlatzone, {});  
}

std::list<std::list<int>> NodeCT::getCopyCNPsByFlatZone(){
    return this->cnpsByFlatzone;
}


void NodeCT::addCNPsOfDisjointFlatzone(std::list<int>&& flatZone, ComponentTree* tree) {
    
    this->cnpsByFlatzone.push_back(std::move(flatZone));
    std::list<int>& newFlatzone = this->cnpsByFlatzone.back();  // ✅ CORRETO

    for (int p : newFlatzone) { 
        tree->updatePixelToFlatzone(p, &newFlatzone);  // ✅ CORRETO (passa referência)
    }
}




void NodeCT::addCNPsOfDisjointFlatzones(std::list<std::list<int>>&& flatZones, ComponentTree* tree) {
    if (flatZones.empty()) {
        return;  // Nada para adicionar, evita operações desnecessárias
    }

    size_t oldSize = this->cnpsByFlatzone.size();

    // Movemos todas as flatzones para `cnpsByFlatzone` sem copiar elementos
    this->cnpsByFlatzone.splice(
        this->cnpsByFlatzone.end(),
        std::move(flatZones)
    );

    // Atualiza `pixelToFlatzone` e `tree->setSC(p, this)`
    auto it = this->cnpsByFlatzone.begin();
    std::advance(it, oldSize);  // Move iterador para o início das novas flatzones

    for (; it != this->cnpsByFlatzone.end(); ++it) {
        std::list<int>& flatzone = *it;  

        for (int p : flatzone) {
            tree->updatePixelToFlatzone(p, &flatzone);
            tree->setSC(p, this);
        }
    }
}


void NodeCT::removeFlatzone(std::list<int>& flatzone) {
    // Encontrar e remover a flatzone da lista cnpsByFlatzone
    for (auto it = cnpsByFlatzone.begin(); it != cnpsByFlatzone.end(); ++it) {
        if (&(*it) == &flatzone) {  // Verifica se a referência coincide
            cnpsByFlatzone.erase(it);
            return;
        }
    }

    assert(false && "Erro: A flatzone a ser removida não está presente na lista de flatzones do nó!");
}


void NodeCT::addCNPsToConnectedFlatzone(std::list<int>&& flatZone, ComponentTree* tree) {
    assert(!flatZone.empty() && "Erro: flatZone passada está vazia!");

    std::unordered_set<std::reference_wrapper<std::list<int>>, ListRefHash, ListRefEqual> flatzonesToMergeSet;
    AdjacencyRelation* adj = tree->getAdjacencyRelation();
    std::list<int>* targetFlatzone = nullptr;
    
    // Encontrar as flatzones adjacentes e definir a maior como `targetFlatzone`
    for (int p : flatZone) {
        for (int np : adj->getAdjPixels(p)) {
            if (p != np && tree->getSC(np) == this) {
                std::list<int>* adjFlatzone = tree->getFlatzonePointer(np);
                if (flatzonesToMergeSet.insert(*adjFlatzone).second) {
                    if (!targetFlatzone || adjFlatzone->size() > targetFlatzone->size()) {
                        targetFlatzone = adjFlatzone;
                    }
                }
            }
        }
    }

    assert(targetFlatzone && "Erro: nenhuma flatzone válida para fusão!");

    // Atualizar SC para os novos pixels antes de modificar `flatZone`
    for (int p : flatZone) {
        tree->setSC(p, this);
    }

    // Guardamos o tamanho antes da fusão, pois `targetFlatzone` pode crescer após os splices
    size_t targetFlatzoneOriginalSize = targetFlatzone->size();

    
    // Fundir todas as flatzones na `targetFlatzone`
    for (auto it = cnpsByFlatzone.begin(); it != cnpsByFlatzone.end(); ) {
        if (flatzonesToMergeSet.count(*it) > 0 && &(*it) != targetFlatzone) {
            targetFlatzone->splice(targetFlatzone->end(), *it);
            it = cnpsByFlatzone.erase(it);  
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
    


}




void NodeCT::addChild(NodeCT* child) {
	this->children.push_back(child);
} 

bool NodeCT::isChild(NodeCT* child){
    auto it = std::find(this->children.begin(), this->children.end(), child);
    return it != children.end();
}

int NodeCT::getIndex() const{ return this->index; }

int NodeCT::getThreshold1() const{ return this->threshold1; }

int NodeCT::getThreshold2() const{ return this->threshold2; }

int NodeCT::getLevel() const{ return this->threshold2; }

void NodeCT::setLevel(int level){ this->threshold2 = level; }

void NodeCT::setArea(long int area){this->areaCC = area;}

long int NodeCT::getArea() const{return this->areaCC;} 

NodeCT* NodeCT::getParent(){  return this->parent; }

void NodeCT::setParent(NodeCT* parent){ this->parent = parent; }

bool NodeCT::isLeaf(){ return this->children.empty(); }

std::list<NodeCT*>& NodeCT::getChildren(){  return this->children; }

int NodeCT::getNumSiblings() {
    if(this->parent != nullptr)
		return this->parent->getChildren().size();
	else
		return 0;
}
