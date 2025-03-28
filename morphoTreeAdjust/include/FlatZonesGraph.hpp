#ifndef FLATZONEGRAPH_HPP
#define FLATZONEGRAPH_HPP

#include <memory>
#include <map>
#include <list>
#include <vector>
#include <queue>
#include <cassert>
#include "../include/Common.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"

class FlatZonesGraph {
private:
    ListOfAdjacentFlatzones listOfAdjacentFlatzones;
    std::list<FlatZone> listFlatZones;

    inline static std::unique_ptr<FlatZonesGraph> instance = nullptr;
    inline static bool forceRebuild = false;
    
    FlatZonesGraph() = default;

    FlatZonesGraph(int* img, int numRows, int numCols, AdjacencyRelation& adj) {
        int numPixels = numRows * numCols;
        int* pixelToFlatzone = new int[numPixels];
        bool* visited = new bool[numPixels]();
        bool* isCountor = new bool[numPixels]();
        listOfAdjacentFlatzones.resize(numPixels);

        for (int p = 0; p < numPixels; p++) {
            if (visited[p]) continue;

            int grayLevel = img[p];
            FlatZone flatZone;
            std::queue<int> queue;
            queue.push(p);
            visited[p] = true;
            listOfAdjacentFlatzones[p] = std::make_unique<AdjacentFlatzones>();

            while (!queue.empty()) {
                int q = queue.front(); queue.pop();
                flatZone.push_back(q);

                for (int nq : adj.getAdjPixels(q)) {
                    if (!visited[nq] && img[nq] == grayLevel) {
                        visited[nq] = true;
                        queue.push(nq);
                    } else if (img[nq] != grayLevel) {
                        isCountor[nq] = true;
                        isCountor[q] = true;
                        pixelToFlatzone[q] = p;
                    }
                }
            }

            assert(!flatZone.empty() && "ERRO: Existem flatzones vazias!");
            listFlatZones.push_back (std::move(flatZone) );
        }

        for (int p = 0; p < numPixels; p++) {
            if (!isCountor[p]) continue;
            int flatZoneID_P = pixelToFlatzone[p];
            AdjacentFlatzones& setP = *listOfAdjacentFlatzones[flatZoneID_P];
            for (int np : adj.getAdjPixels(p)) {
                if (!isCountor[np]) continue;
                int flatZoneID_NP = pixelToFlatzone[np];
                if (flatZoneID_NP != flatZoneID_P) {
                    AdjacentFlatzones& setNP = *listOfAdjacentFlatzones[flatZoneID_NP];
                    setP.insert(flatZoneID_NP);
                    setNP.insert(flatZoneID_P);
                }
            }
        }

        delete[] pixelToFlatzone;
        delete[] isCountor;
        delete[] visited;
    }

    std::unique_ptr<FlatZonesGraph> clone() const {
        std::unique_ptr<FlatZonesGraph> copy = std::unique_ptr<FlatZonesGraph>(new FlatZonesGraph());

        copy->listFlatZones = this->listFlatZones;
        copy->listOfAdjacentFlatzones.resize(this->listOfAdjacentFlatzones.size());
        for (size_t i = 0; i < this->listOfAdjacentFlatzones.size(); ++i) {
            if (this->listOfAdjacentFlatzones[i]) {
                copy->listOfAdjacentFlatzones[i] = std::make_unique<AdjacentFlatzones>(*this->listOfAdjacentFlatzones[i]);
            }
        }
        return copy;
    }

public:
    FlatZonesGraph(const FlatZonesGraph&) = delete;
    FlatZonesGraph& operator=(const FlatZonesGraph&) = delete;
    
    static std::unique_ptr<FlatZonesGraph> createInstance(int* img, int numRows, int numCols, AdjacencyRelation& adj) {
        if (!instance) {
            instance = std::unique_ptr<FlatZonesGraph> ( new FlatZonesGraph(img, numRows, numCols, adj) );
            forceRebuild = false;
        }

        if (forceRebuild) {
            return instance->clone();
        } else {
            forceRebuild = true;
            return std::move(instance);
        }
    }

    std::list<std::list<int>>& getFlatzones() {
        return listFlatZones;
    }
    
    ListOfAdjacentFlatzones& getGraph() {
        return listOfAdjacentFlatzones;
    }
    
    void remapFlatzoneIDInGraph(int oldID, int newID){
        //para manter a propriedade do id da flatzone ser o menor pixel
        for (int neighborID : *listOfAdjacentFlatzones[oldID]) {
            listOfAdjacentFlatzones[neighborID]->erase(oldID);
            listOfAdjacentFlatzones[neighborID]->insert(newID);
        }
        listOfAdjacentFlatzones[newID] = std::move(listOfAdjacentFlatzones[oldID]);
    }

    bool isAdjacent(int flatZoneID, NodeFZPtr node) {
        // Obtém o conjunto de flatzones adjacentes ao nodeFlatZoneID
        const std::unordered_set<int>& adjacentFlatzones = *(listOfAdjacentFlatzones[flatZoneID]);
    
        // Percorre os IDs das flatzones do parent para verificar se há alguma adjacente
        for (auto& [id, flatzone]: node->getCNPsByFlatZone()) {
            if (adjacentFlatzones.find(id) != adjacentFlatzones.end()) {
                return true;  
            }
        }
    
        return false; 
    }

    /**
     * Esse método unifica a lista de flatzones flatZoneNodeList em uma unica flatzone e atualiza o grafo.
     * A propriedade do menor pixel ser o id da flatzone é mantido
     */
    void updateGraphAfterPruning(std::list<FlatZone>& flatZoneNodeList, FlatZone& unifiedFlatzone, NodeFZPtr node, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");
        int unifiedFlatzoneID = node->getRepresentativeCNPs(); //retorna o menor pixel entre as flatzones de node
        FlatZone* unifiedFlatzoneInitial = &node->getFlatZone(unifiedFlatzoneID);
        for(FlatZone& flatZone: flatZoneNodeList)   {
            if(unifiedFlatzoneID > flatZone.front()){
                unifiedFlatzoneID = flatZone.front();
                unifiedFlatzoneInitial = &flatZone;
            }
            
        }
        unifiedFlatzone.splice(unifiedFlatzone.end(), *unifiedFlatzoneInitial); 

        for(FlatZone& flatZone: flatZoneNodeList)   {
            if(flatZone.empty()) continue; //pode acontecer devido ao splice
            //atualiza o grafo e coleta os cnps em uma unica flatzone
            int flatZoneID = flatZone.front();
            std::unordered_set<int> neighborsCopy = *(listOfAdjacentFlatzones[flatZoneID]);
            if(flatZoneID == unifiedFlatzoneID){
                for (int neighborID : neighborsCopy) {
                    //vizinhos que serão unificados
                    if( (tree->isMaxtree() && node->getLevel() <= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && node->getLevel() >= tree->getSC(neighborID)->getLevel())){
                        this->listOfAdjacentFlatzones[neighborID]->erase(unifiedFlatzoneID);
                        this->listOfAdjacentFlatzones[unifiedFlatzoneID]->erase(neighborID);
                    }
                }
            }else{
                for (int neighborID : neighborsCopy) {
                    if(neighborID != unifiedFlatzoneID && this->listOfAdjacentFlatzones[neighborID] ){
                        //vizinhos que NÃO serão unificados
                        if( (tree->isMaxtree() && node->getLevel() >= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && node->getLevel() <= tree->getSC(neighborID)->getLevel())){
                            this->listOfAdjacentFlatzones[neighborID]->erase(flatZoneID);
                            this->listOfAdjacentFlatzones[flatZoneID]->erase(neighborID);
                                
                            this->listOfAdjacentFlatzones[neighborID]->insert(unifiedFlatzoneID);
                            this->listOfAdjacentFlatzones[unifiedFlatzoneID]->insert(neighborID);
                        }
                    }
                }
                unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); // Adicionar pixels ao cnpsCC

                this->listOfAdjacentFlatzones[unifiedFlatzoneID]->erase(flatZoneID);
                //delete this->flatzoneGraph[flatZoneID];
                this->listOfAdjacentFlatzones[flatZoneID] = nullptr;
            }
        }
        
        assert(!unifiedFlatzone.empty() && "ERRO: unifiedFlatzone está vazio após a fusão!");
        assert([&]() {
            int minPixel = *std::min_element(unifiedFlatzone.begin(), unifiedFlatzone.end());
            return minPixel == unifiedFlatzoneID && unifiedFlatzoneID == unifiedFlatzone.front();
        }() && "ERRO: O menor pixel da flatzone unificada não é o seu ID!");
        
    }




    /*
    Esses metodo fundira as flatzones de flatZoneNodeList em uma unica flatzone.
    O primeiro pixel de unifiedFlatzone será o id da nova flatzone.
    A união dos pixels das flatzones flatZoneNodeList formam um componente conexo
    */
    void updateGraph(std::list<FlatZoneNode>& flatZoneNodeList,  FlatZone& unifiedFlatzone, NodeFZPtr tauStar, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");

        int levelTauStar = tauStar->getLevel();
        int unifiedFlatzoneID = std::numeric_limits<int>::max();
        FlatZone* unifiedFlatzoneInitial = nullptr;
        for(FlatZoneNode& flatZoneNode: flatZoneNodeList)   {
            FlatZone& flatZone = *flatZoneNode.flatzone;
            if(!unifiedFlatzoneInitial || unifiedFlatzoneID > flatZone.front()){
                unifiedFlatzoneID = flatZone.front();
                unifiedFlatzoneInitial = &flatZone;
            }
            
        }
        unifiedFlatzone.splice(unifiedFlatzone.end(), *unifiedFlatzoneInitial); 


        for(FlatZoneNode& flatZoneNode: flatZoneNodeList)   {
            FlatZone& flatZone = *flatZoneNode.flatzone;
            if(flatZone.empty())
                continue; //pode acontecer devido ao splice em unifiedFlatzone com unifiedFlatzoneInitial
            

            int flatZoneID = flatZone.front();
            std::unordered_set<int> neighborsCopy = *(listOfAdjacentFlatzones[flatZoneID]);    
            for (int neighborID : neighborsCopy) {
                if(neighborID != unifiedFlatzoneID && this->listOfAdjacentFlatzones[neighborID] ){

                    //vizinhos que NÃO serão unificados
                    if( (tree->isMaxtree() && levelTauStar <= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && levelTauStar >= tree->getSC(neighborID)->getLevel())){
                        this->listOfAdjacentFlatzones[neighborID]->erase(flatZoneID);
                        this->listOfAdjacentFlatzones[flatZoneID]->erase(neighborID);
                            
                        this->listOfAdjacentFlatzones[unifiedFlatzoneID]->insert(neighborID);
                        this->listOfAdjacentFlatzones[neighborID]->insert(unifiedFlatzoneID);
                            
                                
                    }
                }
            }
            
            this->listOfAdjacentFlatzones[unifiedFlatzoneID]->erase(flatZoneID);
            //delete this->flatzoneGraph[flatZoneID];
            this->listOfAdjacentFlatzones[flatZoneID] = nullptr;
            unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); // Adicionar pixels ao cnpsCC
        }
        
        assert([&]() {
            int minPixel = *std::min_element(unifiedFlatzone.begin(), unifiedFlatzone.end());
            return minPixel == unifiedFlatzoneID && unifiedFlatzoneID == unifiedFlatzone.front();
        }() && "ERRO: O menor pixel da flatzone unificada não é o seu ID!");
        
        assert([&]() {
            if (unifiedFlatzone.empty()) {
                std::cerr << "ERRO: unifiedFlatzone está vazia após a fusão!" << std::endl;
                return false;
            }
            
            if (this->listOfAdjacentFlatzones[unifiedFlatzoneID] == nullptr) {
                std::cerr << "ERRO: unifiedFlatzone não está registrada no grafo!" << std::endl;
                return false;
            }

            for (int neighborID : *this->listOfAdjacentFlatzones[unifiedFlatzoneID]) {
                if (this->listOfAdjacentFlatzones[neighborID] == nullptr) {
                    std::cerr << "ERRO: Conexão assimétrica entre unifiedFlatzone e seu vizinho!" << std::endl;
                    std::cerr << "neighborID: " << neighborID << std::endl;
                    return false;
                }
                
                const FlatZone& neighborFlatzone = this->getFlatzoneByID(neighborID);

                if (neighborFlatzone.empty()) {
                    std::cerr << "ERRO: Flatzone vizinha de unifiedFlatzone está vazia APOS fusão!" << std::endl;
                    std::cerr << "neighborID: " << neighborID << std::endl;
                    return false;
                }


            }

            return true;
        }() && "Erro: Grafo de flatzones inconsistente após a fusão!");
    }

    /**
     * Esse método fundirá a flatzoneID com uma outra flatzone adjacentes presente na lista de flatzones do node 
     */
    std::tuple<int, std::unique_ptr<std::list<int>>> mergeConnectedFlatzone(int flatZoneID, NodeFZPtr node, std::shared_ptr<ComponentTreeFZ> tree) 
    {
        assert(listOfAdjacentFlatzones[flatZoneID] && "Erro: flatZone não está registrada no grafo!");
        assert(!listOfAdjacentFlatzones[flatZoneID]->empty() && "Erro: flatZone não tem vizinhos registrados no grafo!");
        
        // aloca a lista dinamicamente
        auto flatzonesToMergeList = std::make_unique<std::list<int>>();
        int unifiedFlatzoneID = std::numeric_limits<int>::max();
    
        for (int neighborID : *listOfAdjacentFlatzones[flatZoneID]) {
            if (tree->getSC(neighborID) == node) {
                unifiedFlatzoneID = std::min(unifiedFlatzoneID, neighborID);
            }
        }
    
        std::list<int>* unifiedFlatzone = &tree->getFlatzoneByID(unifiedFlatzoneID);
    
        listOfAdjacentFlatzones[flatZoneID]->erase(unifiedFlatzoneID);
        listOfAdjacentFlatzones[unifiedFlatzoneID]->erase(flatZoneID);
    
        for (int flatzonMergedID : *listOfAdjacentFlatzones[flatZoneID]) {
            if (tree->getSC(flatzonMergedID) == node) {
                flatzonesToMergeList->push_back(flatzonMergedID); 
            }
        }
    
        for (int flatzonMergedID : *flatzonesToMergeList) {
            for (int neighborID : *listOfAdjacentFlatzones[flatzonMergedID]) {
                if (flatZoneID != neighborID) {
                    listOfAdjacentFlatzones[unifiedFlatzoneID]->insert(neighborID);
                    listOfAdjacentFlatzones[neighborID]->insert(unifiedFlatzoneID); 
                }
                listOfAdjacentFlatzones[neighborID]->erase(flatzonMergedID);
            }
            listOfAdjacentFlatzones[flatzonMergedID] = nullptr;
        }
    
        for (int neighborID : *listOfAdjacentFlatzones[flatZoneID]) {
            listOfAdjacentFlatzones[neighborID]->erase(flatZoneID);
            listOfAdjacentFlatzones[unifiedFlatzoneID]->insert(neighborID);
            listOfAdjacentFlatzones[neighborID]->insert(unifiedFlatzoneID);  
        }
        listOfAdjacentFlatzones[flatZoneID] = nullptr;
    
        // retorno usando std::move no ponteiro
        return std::make_tuple(unifiedFlatzoneID, std::move(flatzonesToMergeList));
    }
    
};

#endif // FLATZONEGRAPH_HPP
