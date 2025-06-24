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
    ListOfAdjacentFlatZones listOfAdjacentFlatZones;
    std::list<FlatZone> listFlatZones;
    std::unordered_map<int, size_t> flatZoneToIndex;
public:

    FlatZonesGraph(const FlatZonesGraph&) = default;
    FlatZonesGraph& operator=(const FlatZonesGraph&) = default;


    FlatZonesGraph(ImageUInt8Ptr imgPtr, AdjacencyRelationPtr adj) {
        int numPixels = imgPtr->getNumRows() * imgPtr->getNumCols();
        auto img = imgPtr->rawData();
        std::unique_ptr<int[]> pixelToFlatzone(new int[numPixels]);
        std::unique_ptr<bool[]> visited(new bool[numPixels]());
        std::unique_ptr<bool[]> isContour(new bool[numPixels]());
        
        for (int p = 0; p < numPixels; p++) {
            if (visited[p]) continue;

            int grayLevel = img[p];
            FlatZone flatZone;
            std::queue<int> queue;
            queue.push(p);
            visited[p] = true;
            flatZoneToIndex[p] = listFlatZones.size();

            while (!queue.empty()) {
                int q = queue.front(); queue.pop();
                flatZone.push_back(q);

                for (int nq : adj->getAdjPixels(q)) {
                    if (!visited[nq] && img[nq] == grayLevel) {
                        visited[nq] = true;
                        queue.push(nq);
                    } else if (img[nq] != grayLevel) {
                        isContour[nq] = true;
                        isContour[q] = true;
                        pixelToFlatzone[q] = p;
                    }
                }
            }

            assert(!flatZone.empty() && "ERRO: Existem flatzones vazias!");
            listFlatZones.push_back (std::move(flatZone) );
        }

        // Inicialização
        listOfAdjacentFlatZones.resize(listFlatZones.size());
        
        for (int p = 0; p < numPixels; p++) {
            if (!isContour[p]) continue;
            int flatZoneID_P = pixelToFlatzone[p];
            int keyFlatZoneID_P = flatZoneToIndex[flatZoneID_P];
            for (int np : adj->getAdjPixels(p)) {
                if (!isContour[np]) continue;
                int flatZoneID_NP = pixelToFlatzone[np];
                if (flatZoneID_NP != flatZoneID_P) {
                    listOfAdjacentFlatZones[ keyFlatZoneID_P ].insert(flatZoneID_NP);
                    listOfAdjacentFlatZones[ flatZoneToIndex[flatZoneID_NP] ].insert(flatZoneID_P);
                }
            }
        }

    }

    AdjacentFlatZones& getAdjacentFlatzones(int flatZoneID) {
        return listOfAdjacentFlatZones[flatZoneToIndex[flatZoneID]];
    }

    std::list<FlatZone>& getFlatzones() {
        return listFlatZones;
    }
    
    void remapFlatzoneIDInGraph(int oldID, int newID){
        int keyOldID = flatZoneToIndex[oldID];
        int keyNewID = flatZoneToIndex[newID];

        //para manter a propriedade do id da flatzone ser o menor pixel
        for (int neighborID : listOfAdjacentFlatZones[keyOldID]) {
            int keyNeighborID = flatZoneToIndex[neighborID];
            auto& neighborSet = listOfAdjacentFlatZones[keyNeighborID];

            neighborSet.erase(oldID); 
            neighborSet.insert(newID);
        }

        listOfAdjacentFlatZones[keyNewID] = std::move(listOfAdjacentFlatZones[keyOldID]);
    }

    bool isAdjacent(int flatZoneID, NodeFZPtr node) {
        // Obtém o conjunto de flatzones adjacentes ao nodeFlatZoneID
        const std::unordered_set<int>& adjacentFlatzones = listOfAdjacentFlatZones[flatZoneToIndex[flatZoneID]];
    
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
    void updateGraphAfterPruning(std::list<FlatZone>&& flatZoneNodeList, FlatZone& unifiedFlatzone, NodeFZPtr node, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");
        FlatZone* unifiedFlatzoneInitial = &flatZoneNodeList.front();
        int unifiedFlatzoneID = unifiedFlatzoneInitial->front();
        for(FlatZone& flatZone: flatZoneNodeList)   {
            if(unifiedFlatzoneID > flatZone.front()){
                unifiedFlatzoneID = flatZone.front();
                unifiedFlatzoneInitial = &flatZone;
            }
            
        }
        unifiedFlatzone.splice(unifiedFlatzone.end(), *unifiedFlatzoneInitial); 
        int keyUnifiedFlatzoneID = flatZoneToIndex[unifiedFlatzoneID];

        for(FlatZone& flatZone: flatZoneNodeList)   {
            if(flatZone.empty()) continue; //pode acontecer devido ao splice

            //atualiza o grafo e coleta os cnps em uma unica flatzone
            int flatZoneID = flatZone.front();
            int keyFlatZoneID = flatZoneToIndex[flatZoneID];
            std::vector<int> neighborsCopy(listOfAdjacentFlatZones[keyFlatZoneID].begin(), listOfAdjacentFlatZones[keyFlatZoneID].end());
            if(flatZoneID == unifiedFlatzoneID){
                for (int neighborID : neighborsCopy) {
                    //vizinhos que serão unificados
                    if( (tree->isMaxtree() && node->getLevel() <= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && node->getLevel() >= tree->getSC(neighborID)->getLevel())){
                        this->listOfAdjacentFlatZones[flatZoneToIndex[neighborID]].erase(unifiedFlatzoneID);
                        this->listOfAdjacentFlatZones[keyUnifiedFlatzoneID].erase(neighborID);
                    }
                }
            }else{
                for (int neighborID : neighborsCopy) {
                    int keyNeighborID = flatZoneToIndex[neighborID];
                    if(neighborID != unifiedFlatzoneID && !this->listOfAdjacentFlatZones[keyNeighborID].empty() ){
                        //vizinhos que NÃO serão unificados
                        if( (tree->isMaxtree() && node->getLevel() >= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && node->getLevel() <= tree->getSC(neighborID)->getLevel())){
                            this->listOfAdjacentFlatZones[keyNeighborID].erase(flatZoneID);
                            this->listOfAdjacentFlatZones[keyFlatZoneID].erase(neighborID);
                                
                            this->listOfAdjacentFlatZones[keyNeighborID].insert(unifiedFlatzoneID);
                            this->listOfAdjacentFlatZones[keyUnifiedFlatzoneID].insert(neighborID);
                        }
                    }
                }
                unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); // Adicionar pixels ao cnpsCC

                this->listOfAdjacentFlatZones[keyUnifiedFlatzoneID].erase(flatZoneID);
                //delete this->flatzoneGraph[flatZoneID];
                this->listOfAdjacentFlatZones[keyFlatZoneID].clear();

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
            FlatZone* flatZone = flatZoneNode.flatzone;
            if(!unifiedFlatzoneInitial || unifiedFlatzoneID > flatZone->front()){
                unifiedFlatzoneID = flatZone->front();
                unifiedFlatzoneInitial = flatZone;
            }
            
        }
        unifiedFlatzone.splice(unifiedFlatzone.end(), *unifiedFlatzoneInitial); 
        int keyUnifiedFlatzoneID = flatZoneToIndex[unifiedFlatzoneID];

        for(FlatZoneNode& flatZoneNode: flatZoneNodeList)   {
            FlatZone* flatZone = flatZoneNode.flatzone;
            if(flatZone->empty())
                continue; //pode acontecer devido ao splice em unifiedFlatzone com unifiedFlatzoneInitial
            
            int flatZoneID = flatZone->front();
            int keyFlatZoneID = flatZoneToIndex[flatZoneID];
            std::vector<int> neighborsCopy(listOfAdjacentFlatZones[keyFlatZoneID].begin(), listOfAdjacentFlatZones[keyFlatZoneID].end());
            for (int neighborID : neighborsCopy) {
                int keyNeighborID = flatZoneToIndex[neighborID];
                if(neighborID != unifiedFlatzoneID && !this->listOfAdjacentFlatZones[keyNeighborID].empty() ){

                    //vizinhos que NÃO serão unificados
                    if( (tree->isMaxtree() && levelTauStar <= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && levelTauStar >= tree->getSC(neighborID)->getLevel())){
                        this->listOfAdjacentFlatZones[keyNeighborID].erase(flatZoneID);
                        this->listOfAdjacentFlatZones[keyFlatZoneID].erase(neighborID);
                            
                        this->listOfAdjacentFlatZones[keyUnifiedFlatzoneID].insert(neighborID);
                        this->listOfAdjacentFlatZones[keyNeighborID].insert(unifiedFlatzoneID);
                                 
                    }
                }
            }

            this->listOfAdjacentFlatZones[flatZoneToIndex[unifiedFlatzoneID]].erase(flatZoneID);
            this->listOfAdjacentFlatZones[flatZoneToIndex[flatZoneID]].clear();
            unifiedFlatzone.splice(unifiedFlatzone.end(), *flatZone); // Adicionar pixels ao cnpsCC

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
            
            if (this->listOfAdjacentFlatZones[flatZoneToIndex[unifiedFlatzoneID]].empty()) {
                std::cerr << "ERRO: unifiedFlatzone não está registrada no grafo!" << std::endl;
                return false;
            }

            for (int neighborID : this->listOfAdjacentFlatZones[flatZoneToIndex[unifiedFlatzoneID]]) {
                if (this->listOfAdjacentFlatZones[flatZoneToIndex[neighborID]].empty()) {
                    std::cerr << "ERRO: Conexão assimétrica entre unifiedFlatzone e seu vizinho!" << std::endl;
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
    std::tuple<int, std::list<int>> mergeConnectedFlatzone(int flatZoneID, NodeFZPtr node, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(flatZoneToIndex.count(flatZoneID) && "flatZoneID inválido");
        auto& neighborsFlatZoneID = listOfAdjacentFlatZones[flatZoneToIndex[flatZoneID]];
        assert(!neighborsFlatZoneID.empty() && "Erro: flatZone não tem vizinhos registrados no grafo!");
        
        // aloca a lista dinamicamente
        std::list<int> flatzonesToMergeList;
        int unifiedFlatzoneID = std::numeric_limits<int>::max();
        for (int neighborID : neighborsFlatZoneID) {
            if (tree->getSC(neighborID) == node) {
                unifiedFlatzoneID = std::min(unifiedFlatzoneID, neighborID);
            }
        }
        
        for (int flatzonMergedID : neighborsFlatZoneID) {
            if (unifiedFlatzoneID != flatzonMergedID && tree->getSC(flatzonMergedID) == node) {
                flatzonesToMergeList.push_back(flatzonMergedID); 
            }
        }

        listOfAdjacentFlatZones[flatZoneToIndex[flatZoneID]].erase(unifiedFlatzoneID);
        listOfAdjacentFlatZones[flatZoneToIndex[unifiedFlatzoneID]].erase(flatZoneID);
        

        for (int flatzonMergedID : flatzonesToMergeList) {
            auto& adjFZ = listOfAdjacentFlatZones[flatZoneToIndex[flatzonMergedID]];
            for (int neighborID : adjFZ) {
                if (flatZoneID != neighborID) {
                    listOfAdjacentFlatZones[flatZoneToIndex[unifiedFlatzoneID]].insert(neighborID);
                    listOfAdjacentFlatZones[flatZoneToIndex[neighborID]].insert(unifiedFlatzoneID); 
                }
                listOfAdjacentFlatZones[flatZoneToIndex[neighborID]].erase(flatzonMergedID);
            }
            adjFZ.clear(); //listOfAdjacentFlatZones[flatZoneToIndex[flatzonMergedID]].clear();
        }

        
        for (int neighborID : std::vector<int>(neighborsFlatZoneID.begin(), neighborsFlatZoneID.end())) {
            listOfAdjacentFlatZones[flatZoneToIndex[neighborID]].erase(flatZoneID);
            listOfAdjacentFlatZones[flatZoneToIndex[unifiedFlatzoneID]].insert(neighborID);
            listOfAdjacentFlatZones[flatZoneToIndex[neighborID]].insert(unifiedFlatzoneID);  
        }
        listOfAdjacentFlatZones[flatZoneToIndex[flatZoneID]].clear();
    
        // retorno usando std::move no ponteiro
        return std::make_tuple(unifiedFlatzoneID, std::move(flatzonesToMergeList));
    }

};

using FlatZonesGraphPtr = std::unique_ptr<FlatZonesGraph>;

#endif // FLATZONEGRAPH_HPP
