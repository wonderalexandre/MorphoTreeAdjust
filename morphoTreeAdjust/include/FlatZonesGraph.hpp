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
    std::vector<int> indexToFlatZone;
    std::vector<int> parent; // DSU (representantes das flatzones)
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


        int numFZ = listFlatZones.size();
        
        // Inicialização do union-find
        parent.resize(numFZ);
        for (int i = 0; i < numFZ; ++i) 
            parent[i] = i;

        // Mapear índice para flatzoneID
        indexToFlatZone.resize(listFlatZones.size());

        // Inicialização da lista de adjacência
        listOfAdjacentFlatZones.resize(numFZ);

        //criação da lista de adjacência
        for (int p = 0; p < numPixels; p++) {
            if (!isContour[p]) continue;
            int flatZoneID_P = pixelToFlatzone[p];
            int keyFlatZoneID_P = flatZoneToIndex[flatZoneID_P];
            indexToFlatZone[keyFlatZoneID_P] = flatZoneID_P;
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

    int find(int i) {
        if (parent[i] != i) {
            parent[i] = find(parent[i]); // compressão de caminho
        }
        return parent[i];
    }

    void unite(int idxUnified, int idxMerged) {
        int rootUnified = find(idxUnified);
        int rootMerged = find(idxMerged);

        if (rootUnified == rootMerged) return;

        parent[rootMerged] = rootUnified;  // rootUnified será o representante
        
    }

    int findRepresentative(int flatZoneID) {
        int idx = flatZoneToIndex[flatZoneID];
        int rep = find(idx);
        return indexToFlatZone[rep];
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

        // Atualiza union-find
        parent[keyOldID] = find(keyNewID);
        int repIdx = find(keyNewID);
        indexToFlatZone[keyOldID] = indexToFlatZone[repIdx];
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
     * Esse método fundirá a flatzoneID com outras flatzones adjacentes presente na lista de flatzones do node 
     */
    std::tuple<int, std::list<int>> mergeConnectedFlatzone(int flatZoneID, NodeFZPtr node, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(flatZoneToIndex.count(flatZoneID) && "flatZoneID inválido");

        int idxFlatZone = flatZoneToIndex[flatZoneID];
        auto& adjFZ = listOfAdjacentFlatZones[idxFlatZone];
        assert(!adjFZ.empty() && "Erro: flatZone não tem vizinhos registrados!");

        std::list<int> flatzonesToMergeList;
        int unifiedFlatzoneID = std::numeric_limits<int>::max();

        // Encontrar o menor ID entre os vizinhos pertencentes ao mesmo node
        for (int neighborID : adjFZ) {
            if (tree->getSC(neighborID) == node) {
                unifiedFlatzoneID = std::min(unifiedFlatzoneID, neighborID);
            }
        }

        int idxUnified = flatZoneToIndex[unifiedFlatzoneID];
        if (unifiedFlatzoneID < flatZoneID) {
            unite(idxUnified, idxFlatZone);
        }

        // União com os demais
        for (int mergedID : adjFZ) {
            if (mergedID != unifiedFlatzoneID && tree->getSC(mergedID) == node) {
                flatzonesToMergeList.push_back(mergedID);
                unite(idxUnified, flatZoneToIndex[mergedID]);
            }
        }

        // Atualizar adjacências
        for (int mergedID : flatzonesToMergeList) {
            int idxMerged = flatZoneToIndex[mergedID];
            auto& neighborsMerged = listOfAdjacentFlatZones[idxMerged];

            for (int neighborID : neighborsMerged) {
                if (neighborID != flatZoneID) {
                    int idxNeighbor = flatZoneToIndex[neighborID];
                    listOfAdjacentFlatZones[idxUnified].insert(neighborID);
                    listOfAdjacentFlatZones[idxNeighbor].insert(unifiedFlatzoneID);
                    listOfAdjacentFlatZones[idxNeighbor].erase(mergedID);
                }
            }

            // Libera memória interna do conjunto
            std::unordered_set<int>().swap(neighborsMerged); 
        }

        // Remoção da relação flatZoneID ↔ unified
        listOfAdjacentFlatZones[idxFlatZone].erase(unifiedFlatzoneID);
        listOfAdjacentFlatZones[idxUnified].erase(flatZoneID);

        // Atualizar adjacência dos antigos vizinhos do flatZoneID
        for (int neighborID : adjFZ) {
        //for (int neighborID : std::vector<int>(adjFZ.begin(), adjFZ.end())) {
            int idxNeighbor = flatZoneToIndex[neighborID];
            listOfAdjacentFlatZones[idxNeighbor].erase(flatZoneID);
            listOfAdjacentFlatZones[idxNeighbor].insert(unifiedFlatzoneID);
            listOfAdjacentFlatZones[idxUnified].insert(neighborID);
        }

        // Limpa adjacência da flatZoneID
        std::unordered_set<int>().swap(adjFZ); 

        return { unifiedFlatzoneID, std::move(flatzonesToMergeList) };
    }

    /**
     * Esse método fundirá a flatzoneID com uma outra flatzone adjacentes presente na lista de flatzones do node 
     *
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
        
        int idxUnified = flatZoneToIndex[unifiedFlatzoneID];
        if(unifiedFlatzoneID < flatZoneID){
            int idxflatZoneID = flatZoneToIndex[flatZoneID];
            unite(idxUnified, idxflatZoneID);
        }

        for (int flatzonMergedID : neighborsFlatZoneID) {
            if (unifiedFlatzoneID != flatzonMergedID && tree->getSC(flatzonMergedID) == node) {
                flatzonesToMergeList.push_back(flatzonMergedID); //aplicar union

                //: aplica union no índice
                int idxMerged = flatZoneToIndex[flatzonMergedID];
                unite(idxUnified, idxMerged);
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
    }*/


    /*
    Esses metodo fundira as flatzones de flatZoneNodeList em uma unica flatzone.
    A união dos pixels das flatzones flatZoneNodeList formam um componente conexo
    */
    void updateGraph(std::vector<int>& flatZonesID, int unifiedFlatzoneID, NodeFZPtr tauStar, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");
        int levelTauStar = tauStar->getLevel();
        
        
        //int unifiedFlatzoneID = unifiedFlatzone->front();
        int keyUnified = flatZoneToIndex[unifiedFlatzoneID];
        int idxUnified = find(keyUnified);

        std::unordered_set<int>& adjUnified = listOfAdjacentFlatZones[idxUnified];
        for(size_t i=0; i<flatZonesID.size(); i++){    
            int flatZoneID = flatZonesID[i];
            int keyFused = flatZoneToIndex[flatZoneID];
            int idxFused = find(keyFused);
            if (idxFused == idxUnified) continue;
            
            std::unordered_set<int>& adjFused = listOfAdjacentFlatZones[idxFused];

            // Remover self-loop
            adjFused.erase(unifiedFlatzoneID);
            adjUnified.erase(flatZoneID);

            for (int neighborID : adjFused) {
                int neighborLevel = tree->getSC(neighborID)->getLevel();
                if ((tree->isMaxtree() && levelTauStar <= neighborLevel) || (!tree->isMaxtree() && levelTauStar >= neighborLevel)) {
                    int idxNeighbor = flatZoneToIndex[neighborID];
                    listOfAdjacentFlatZones[idxNeighbor].erase(flatZoneID);
                    listOfAdjacentFlatZones[idxNeighbor].insert(unifiedFlatzoneID);

                    adjUnified.insert(neighborID);
                }
            }
            
            std::unordered_set<int>().swap(adjFused); // desaloca
            //unifiedFlatzone->splice(unifiedFlatzone->end(), *flatZones[i]); 

            if (idxFused != idxUnified) { //atualiza union-find, mesmo que: unite(idxUnified, idxFused);
                parent[idxFused] = idxUnified;
            }
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


        

        //return unifiedFlatzone;
    }



    /*
    Esses metodo fundira as flatzones de flatZoneNodeList em uma unica flatzone.
    A união dos pixels das flatzones flatZoneNodeList formam um componente conexo
    *
    void updateGraph(std::vector<FlatZonePtr>& flatZoneNodeList,  FlatZone& unifiedFlatzone, NodeFZPtr tauStar, std::shared_ptr<ComponentTreeFZ> tree) {
        assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");

        int levelTauStar = tauStar->getLevel();
        int unifiedFlatzoneID = std::numeric_limits<int>::max();
        FlatZone* unifiedFlatzoneInitial = nullptr;
        for(FlatZonePtr& flatZone: flatZoneNodeList)   {
            if(!unifiedFlatzoneInitial || unifiedFlatzoneID > flatZone->front()){
                unifiedFlatzoneID = flatZone->front();
                unifiedFlatzoneInitial = flatZone;
            }
            
        }
        unifiedFlatzone.splice(unifiedFlatzone.end(), *unifiedFlatzoneInitial); 
        int keyUnifiedFlatzoneID = flatZoneToIndex[unifiedFlatzoneID];
        int idxUnified = find(flatZoneToIndex[unifiedFlatzoneID]);

        *
        for(FlatZonePtr& flatZone: flatZoneNodeList)   {
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

            // atualiza union-find
            int idxFused = find(flatZoneToIndex[flatZoneID]);
            if (idxFused != idxUnified) {
                parent[idxFused] = idxUnified;
            }
        }*

        std::unordered_set<int>& adjUnified = listOfAdjacentFlatZones[idxUnified];

        for (FlatZonePtr& flatZone : flatZoneNodeList) {
            if (flatZone->empty()) continue;

            int flatZoneID = flatZone->front();
            int idxFused = find(flatZoneToIndex[flatZoneID]);
            if (idxFused == idxUnified) continue;

            std::unordered_set<int>& adjFused = listOfAdjacentFlatZones[idxFused];

            // Remover self-loop
            adjFused.erase(unifiedFlatzoneID);
            adjUnified.erase(flatZoneID);

            for (int neighborID : adjFused) {
                int neighborLevel = tree->getSC(neighborID)->getLevel();
                bool keepAsNeighbor =
                    (tree->isMaxtree() && levelTauStar <= neighborLevel) ||
                    (!tree->isMaxtree() && levelTauStar >= neighborLevel);

                if (keepAsNeighbor) {
                    int idxNeighbor = flatZoneToIndex[neighborID];
                    listOfAdjacentFlatZones[idxNeighbor].erase(flatZoneID);
                    listOfAdjacentFlatZones[idxNeighbor].insert(unifiedFlatzoneID);

                    adjUnified.insert(neighborID);
                }
            }
            unifiedFlatzone.splice(unifiedFlatzone.end(), *flatZone); // Adicionar pixels ao cnpsCC

            adjFused.clear(); // Não mais necessário
            if (idxFused != idxUnified) {
                parent[idxFused] = idxUnified;
            }
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
*/



    /**
     * Esse método unifica a lista de flatzones flatZoneNodeList em uma unica flatzone e atualiza o grafo.
     * A propriedade do menor pixel ser o id da flatzone é mantido
     *
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
        
    }*/



};

using FlatZonesGraphPtr = std::shared_ptr<FlatZonesGraph>;

#endif // FLATZONEGRAPH_HPP
