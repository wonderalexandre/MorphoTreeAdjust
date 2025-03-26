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

class GraphOfFlatZones {
private:
    FlatzoneGraph graph;
    std::list<FlatZone> listFlatZones;

    inline static std::unique_ptr<GraphOfFlatZones> instance = nullptr;
    inline static bool forceRebuild = false;
    
    GraphOfFlatZones() = default;

    GraphOfFlatZones(int* img, int numRows, int numCols, AdjacencyRelation& adj) {
        int numPixels = numRows * numCols;
        int* pixelToFlatzone = new int[numPixels];
        bool* visited = new bool[numPixels]();
        bool* isCountor = new bool[numPixels]();
        graph.resize(numPixels);

        for (int p = 0; p < numPixels; p++) {
            if (visited[p]) continue;

            int grayLevel = img[p];
            FlatZone flatZone;
            std::queue<int> queue;
            queue.push(p);
            visited[p] = true;
            graph[p] = std::make_unique<AdjacentFlatzones>();

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
            AdjacentFlatzones& setP = *graph[flatZoneID_P];
            for (int np : adj.getAdjPixels(p)) {
                if (!isCountor[np]) continue;
                int flatZoneID_NP = pixelToFlatzone[np];
                if (flatZoneID_NP != flatZoneID_P) {
                    AdjacentFlatzones& setNP = *graph[flatZoneID_NP];
                    setP.insert(flatZoneID_NP);
                    setNP.insert(flatZoneID_P);
                }
            }
        }

        delete[] pixelToFlatzone;
        delete[] isCountor;
        delete[] visited;
    }

    std::unique_ptr<GraphOfFlatZones> clone() const {
        std::unique_ptr<GraphOfFlatZones> copy = std::unique_ptr<GraphOfFlatZones>(new GraphOfFlatZones());

        copy->listFlatZones = this->listFlatZones;
        copy->graph.resize(this->graph.size());
        for (size_t i = 0; i < this->graph.size(); ++i) {
            if (this->graph[i]) {
                copy->graph[i] = std::make_unique<AdjacentFlatzones>(*this->graph[i]);
            }
        }
        return copy;
    }

public:
    GraphOfFlatZones(const GraphOfFlatZones&) = delete;
    GraphOfFlatZones& operator=(const GraphOfFlatZones&) = delete;
    



    static std::unique_ptr<GraphOfFlatZones> createInstance(int* img, int numRows, int numCols, AdjacencyRelation& adj) {
        if (!instance) {
            instance = std::unique_ptr<GraphOfFlatZones> ( new GraphOfFlatZones(img, numRows, numCols, adj) );
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
    
    /*FlatzoneGraph&& releaseflatzoneGraph() {
        return std::move(flatzoneGraph);
    }*/

    FlatzoneGraph& getGraph() {
        return graph;
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
            std::unordered_set<int> neighborsCopy = *(graph[flatZoneID]);
            if(flatZoneID == unifiedFlatzoneID){
                for (int neighborID : neighborsCopy) {
                    //vizinhos que serão unificados
                    if( (tree->isMaxtree() && node->getLevel() <= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && node->getLevel() >= tree->getSC(neighborID)->getLevel())){
                        this->graph[neighborID]->erase(unifiedFlatzoneID);
                        this->graph[unifiedFlatzoneID]->erase(neighborID);
                    }
                }
            }else{
                for (int neighborID : neighborsCopy) {
                    if(neighborID != unifiedFlatzoneID && this->graph[neighborID] ){
                        //vizinhos que NÃO serão unificados
                        if( (tree->isMaxtree() && node->getLevel() >= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && node->getLevel() <= tree->getSC(neighborID)->getLevel())){
                            this->graph[neighborID]->erase(flatZoneID);
                            this->graph[flatZoneID]->erase(neighborID);
                                
                            this->graph[neighborID]->insert(unifiedFlatzoneID);
                            this->graph[unifiedFlatzoneID]->insert(neighborID);
                        }
                    }
                }
                unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); // Adicionar pixels ao cnpsCC

                this->graph[unifiedFlatzoneID]->erase(flatZoneID);
                //delete this->flatzoneGraph[flatZoneID];
                this->graph[flatZoneID] = nullptr;
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
            std::unordered_set<int> neighborsCopy = *(graph[flatZoneID]);    
            for (int neighborID : neighborsCopy) {
                if(neighborID != unifiedFlatzoneID && this->graph[neighborID] ){

                    //vizinhos que NÃO serão unificados
                    if( (tree->isMaxtree() && levelTauStar <= tree->getSC(neighborID)->getLevel()) || (!tree->isMaxtree() && levelTauStar >= tree->getSC(neighborID)->getLevel())){
                        this->graph[neighborID]->erase(flatZoneID);
                        this->graph[flatZoneID]->erase(neighborID);
                            
                        this->graph[unifiedFlatzoneID]->insert(neighborID);
                        this->graph[neighborID]->insert(unifiedFlatzoneID);
                            
                                
                    }
                }
            }
            
            this->graph[unifiedFlatzoneID]->erase(flatZoneID);
            //delete this->flatzoneGraph[flatZoneID];
            this->graph[flatZoneID] = nullptr;
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
            
            if (this->graph[unifiedFlatzoneID] == nullptr) {
                std::cerr << "ERRO: unifiedFlatzone não está registrada no grafo!" << std::endl;
                return false;
            }

            for (int neighborID : *this->graph[unifiedFlatzoneID]) {
                if (this->graph[neighborID] == nullptr) {
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



    
};

#endif // FLATZONEGRAPH_HPP
