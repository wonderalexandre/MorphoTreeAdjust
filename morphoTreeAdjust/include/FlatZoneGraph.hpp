#include "../include/Common.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include<map>

class FlatZoneGraph{

    private:
        std::vector<std::unique_ptr<AdjacentFlatzones>> adjacencyRelations;
        std::map<int, std::list<int>> flatzones;
        
    public:

        FlatZoneGraph() = default;

        FlatZoneGraph(int* img, int numRows, int numCols, AdjacencyRelation& adj) {
            int numPixels = numRows * numCols;
            int* pixelToFlatzone = new int[numPixels];
            bool* visited = new bool[numPixels]();
            bool* isCountor = new bool[numPixels]();
            adjacencyRelations.resize(numPixels);  
        
            for (int p = 0; p < numPixels; p++) {
                if (visited[p]) continue; 
                assert(pixelToNode[p] != nullptr && "Falha no mapeamento SC");
                    
                int grayLevel = img[p];
                std::list<int> flatZone;
                std::queue<int> queue;
                queue.push(p);
                visited[p] = true;
                adjacencyRelations[p] = std::make_unique<AdjacentFlatzones>();
        
        
                while (!queue.empty()) {
                    int q = queue.front(); queue.pop();
                    flatZone.push_back(q);
                    
                    for (int nq : adj.getAdjPixels(q)) {
                        if (!visited[nq] && img[nq] == grayLevel) {
                            visited[nq] = true;
                            queue.push(nq);
                        }
                        else if (img[nq] != grayLevel){
                            isCountor[nq] = true;
                            isCountor[q] = true;
                            pixelToFlatzone[q] = p; //id da flatzone
                        }
                    }
        
                }
                assert(flatZone.size() > 0 && "ERRO: Existem flatzones vazias!");
                flatzones[p] = std::move(flatZone);
            }
        
            //build graph
            for (int p = 0; p < numPixels; p++) {
                if(!isCountor[p]) continue;
                int flatZoneID_P = pixelToFlatzone[p];
                AdjacentFlatzones& setP = *adjacencyRelations[flatZoneID_P];
                for (int np : adj.getAdjPixels(p)) {
                    if(!isCountor[np]) continue;
                    int flatZoneID_NP = pixelToFlatzone[np];
                    if (flatZoneID_NP != flatZoneID_P) {
                        AdjacentFlatzones& setNP = *adjacencyRelations[flatZoneID_NP];
                        setP.insert(flatZoneID_NP);
                        setNP.insert(flatZoneID_P);
                    }
                }
            }
        
            delete[] pixelToFlatzone;
            delete[] isCountor;
            delete[] visited;
        
        }


        FlatZoneGraph clone() const {
            FlatZoneGraph copy;
            copy.flatzones = this->flatzones;
            copy.adjacencyRelations.resize(adjacencyRelations.size());
            for (size_t i = 0; i < adjacencyRelations.size(); ++i) {
                if (adjacencyRelations[i]) {
                    copy.adjacencyRelations[i] = std::make_unique<AdjacentFlatzones>(*adjacencyRelations[i]);
                }
            }
            return copy;
        }
        

};