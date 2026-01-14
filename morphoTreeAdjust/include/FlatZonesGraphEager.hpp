

#ifndef FLATZONEGRAPH_HPP
#define FLATZONEGRAPH_HPP

#include <memory>
#include <vector>
#include <utility> 
#include <span>
#include <cassert>
#include <climits> 
#include <cmath>

#include <numeric>
#include <algorithm> 
#include <iostream>
#include <iomanip> // para formatar saída

#include "../include/Common.hpp"
#include "../include/AdjacencyRelation.hpp"


/**
 * @brief Grafo de flat-zones com união eficiente e adjacências deduplicadas.
 *
 * `FlatZonesGraphEager` constrói e mantém um grafo cujos vértices são flat-zones
 * (componentes conexos de nível de cinza constante) extraídas de uma imagem.
 * Ele provê:
 *  - construção das flat-zones por BFS sobre `AdjacencyRelation`;
 *  - representação compacta das zonas via `PixelSetManager` (listas circulares);
 *  - união/fusão de zonas com **splice O(1)** das listas de pixels;
 *  - manutenção de adjacências sem duplicatas através de `AdjacentFlatZonesSet`;
 *  - estrutura Union-Find (`parent`) com compressão de caminho para consultas rápidas.
 *
 * ## Estruturas internas
 * - `pixelBuffer`/`pixelView` (`PixelSetManager`): mapeia `pixel ↔ índice de set`,
 *   tamanho do set e lista circular `pixelsNext` para percorrer os pixels de cada flat-zone.
 * - `parent` (`std::vector<int>`): Union-Find sobre **índices de flat-zones** (slots)
 *   para rastrear fusões (find com path compression).
 * - `listOfAdjacentFlatZones` (`std::vector<AdjacentFlatZonesSet>`): lista de vizinhos
 *   por slot (armazenados por **representantes de pixel**), com espelhamento e
 *   deduplicação via `mirrorAndFinalize`.
 * 
 * ## Complexidade (típica)
 * - Construção: O(N) nos pixels + O(E) nas arestas entre zonas (com prefiltro local).
 * - `find` (UF): α(N) (quase O(1)).
 * - `unite` (UF + splice de listas): O(1) amortizado.
 * - Operações de adjacência: proporcionais ao grau da zona envolvida.
 *
 * ## Invariantes
 * - Slots invalidados após fusões têm `indexToPixel[i] == -1` e `sizeSets[i] == 0`.
 * - `pixelToIndex[rep]` aponta para o slot (raiz após `find`) da sua flat-zone.
 * - Não existem auto-laços nas adjacências; as listas são simétricas.
 *
 * ## Exemplo
 * @code
 * auto img = ImageUInt8::create(rows, cols); // preenchida previamente
 * auto adj = std::make_shared<AdjacencyRelation>(rows, cols, 1.5);
 * FlatZonesGraphEager fzg(img, adj);
 *
 * // Percorre pixels de uma flat-zone por seu representante de pixel
 * int rep =  menor pixel da FZ ;
 * for (int p : fzg.getPixelsByFlatzone(rep)) {
 *     // processa p
 * }
 *
 * // Funde base com candidatos realmente adjacentes
 * std::vector<int> base   = {repA, repB};
 * std::vector<int> cands  = {repC, repD, repE};
 * int winner = fzg.mergeBasesWithAdjacentCandidatesInPlace(base, cands);
 * @endcode
 */
class FlatZonesGraphEager {
private:

    ImageUInt8Ptr imgPtr;
    AdjacencyRelationPtr adj;
    std::vector<AdjacentFlatZones> listOfAdjacentFlatZones;

    std::shared_ptr<PixelSetManager> pixelBuffer;
    PixelSetManager::View pixelView;

    std::vector<int> parent; // unin-find: registrará a fusão das flatzones. Tamanho: numFlatZones
    std::vector<int> tmpCandidatesBuffer;
    int tmpCandidatesSize{0};
    

    /**
     * Método auxiliar do union-find (DSU) para encontrar o representante de um conjunto.
     * @param i Índice do elemento.
     * @return Índice do representante do conjunto ao qual o elemento pertence.
     */
    int find(int index) {
        if (parent[index] != index) {
            parent[index] = find(parent[index]); // compressão de caminho
        }
        return parent[index];
    }

    // --- helpers defensivos ---
    inline bool isRootIndex_(int idx) const {
        return idx >= 0 && idx < (int)parent.size() && parent[idx] == idx;
    }

    inline int repIndexFromPixel_(int repPixel) const {
        // pixelView.pixelToIndex[repPixel] existe no seu código atual
        if (repPixel < 0 || repPixel >= (int)pixelView.pixelToIndex.size()) return -1;
        return pixelView.pixelToIndex[repPixel];
    }

    /**
     * @brief Une duas flat zones no grafo de flat zones.
     *
     * Este método realiza a fusão completa entre duas flat zones, escolhendo uma como
     * vencedora (pré-condição: o chamador já fornece o índice do vencedor, que deve
     * corresponder ao menor pixel representante).
     *
     * A fusão envolve:
     *   - Atualização da estrutura union-find (`parent` e `size`);
     *   - Splice em O(1) das listas circulares de pixels (`pixelsNext`);
     *   - Consolidação do tamanho da flat zone vencedora;
     *   - Invalidação do slot perdedor em `indexToPixel` (não mais raiz ativa);
     *   - Redirecionamento de consultas antigas (`pixelToIndex`) para a raiz vencedora.
     *
     * Pré-condições:
     *   - `winnerIndex` e `loserIndex` são índices válidos de flat zones;
     *   - O chamador já garantiu que o representante do vencedor é o menor pixel.
     *
     * @param idxWinner Índice (slot) da flat zone vencedora.
     * @param idxLoser  Índice (slot) da flat zone perdedora.
     */
    void unite(int idxWinner, int idxLoser) {
        // 1. Normaliza para as raízes
        int idxRootWinner = find(idxWinner);
        int idxRootLoser  = find(idxLoser);
        if (idxRootWinner == idxRootLoser) return;

        // 2. Recupera representantes (pixels)
        int repWinner = pixelView.indexToPixel[idxRootWinner];
        int repLoser  = pixelView.indexToPixel[idxRootLoser];

        // 3. Atualiza union-find: loser -> winner
        parent[idxRootLoser] = idxRootWinner;
        pixelView.sizeSets[idxRootWinner] += pixelView.sizeSets[idxRootLoser];

        // 4. Splice O(1) das listas circulares (pixels)
        int nextWinner = pixelView.pixelsNext[repWinner];
        int nextLoser  = pixelView.pixelsNext[repLoser];
        pixelView.pixelsNext[repWinner] = nextLoser;
        pixelView.pixelsNext[repLoser]  = nextWinner;

        // 5. Invalida slot perdedor
        pixelView.sizeSets[idxRootLoser]  = 0;
        pixelView.indexToPixel[idxRootLoser] = -1;

        // 6. Redireciona lookups pelo antigo rep pixel
        pixelView.pixelToIndex[repLoser] = idxRootWinner;
    }

    void ensureTmpCandidatesCapacity(int required) {
        if (required > (int)tmpCandidatesBuffer.size()) {
            tmpCandidatesBuffer.resize(required);
        }
    }

public:
    

    FlatZonesGraphEager(const FlatZonesGraphEager&) = default;
    FlatZonesGraphEager& operator=(const FlatZonesGraphEager&) = default;    
    FlatZonesGraphEager(ImageUInt8Ptr imgPtr, int radiusAdj) : FlatZonesGraphEager(imgPtr, std::make_shared<AdjacencyRelation>(imgPtr->getNumRows(), imgPtr->getNumCols(), radiusAdj)) {}
    
    /**
     * Construtor que cria o grafo de flatzones a partir de uma imagem e uma relação de adjacência.
     * @param imgPtr Ponteiro para a imagem de entrada.
     * @param adj Ponteiro para a relação de adjacência que define os vizinhos dos pixels.
     */
    FlatZonesGraphEager(ImageUInt8Ptr imgPtr, AdjacencyRelationPtr adj)
        : imgPtr(imgPtr), adj(adj), pixelBuffer(std::make_shared<PixelSetManager>(imgPtr->getSize())){

        int numPixels = imgPtr->getSize();
        auto img = imgPtr->rawData();
        std::vector<uint8_t> visited(numPixels, 0);
        std::vector<uint8_t> isBoundary(numPixels, 0); 
        
        auto& [pixelToIndex, indexToPixel, sizeSets, pixelsNext] = *pixelBuffer;
        int numFZ = 0;
        FastQueue<int> queue(numPixels/4);
        for (int p = 0; p < numPixels; ++p) {
            if (visited[p]) continue;

            int tail = p; 
            int sizeFlatzone = 0;
            int idxFZ = numFZ++;  // índice desta flatzone
            queue.push(p);
            visited[p] = 1;
            pixelToIndex[p] = idxFZ;   // já marca o seed
            indexToPixel[idxFZ] = p;
            
            while (!queue.empty()) {
                int q = queue.pop();
                sizeFlatzone++;
                bool hasDiff = false;
                for (int nq : adj->getNeighborPixels(q)) {
                    if (!visited[nq] && img[nq] == img[p]) {
                        visited[nq] = 1;
                        queue.push(nq);
                        pixelToIndex[nq] = idxFZ;           
                        pixelsNext[tail] = nq;
                        tail = nq;
                    }
                    else if (img[nq] != img[p]) {
                        hasDiff = true; // q é de borda
                    }
                }
                if (hasDiff) {
                    isBoundary[q] = 1;
                }
                
            }
            pixelsNext[tail] = p; // fecha ciclo da lista circular
            sizeSets[idxFZ] = sizeFlatzone;
        }
        pixelBuffer->shrinkToNumSets(numFZ);
        pixelView = pixelBuffer->view();

        parent.resize(numFZ);
        std::iota(parent.begin(), parent.end(), 0);

        listOfAdjacentFlatZones.resize(numFZ);
        
        auto guess_degree = [](int area) -> int { // Heurística: grau ~ O(sqrt(area)). Bons defaults para 4/8-adj.
            int g = 10 + int(2.2 * std::sqrt(double(area)));
            if (g > 64) g = 64;
            return g;
        };
        for (int idx = 0; idx < numFZ; ++idx) listOfAdjacentFlatZones[idx].reserve( guess_degree(sizeSets[idx]) );
        

        
        //pre-filtro: evita repetir o mesmo (idxP,idxQ) dentro da mesma FZ
        LocalPrefilter64 prefilter;
        
        // Gera adjacências: para cada aresta (p,n) com níveis diferentes
        for (int p = 0; p < numPixels; ++p) {
            if (!isBoundary[p]) continue; // p não é de borda: ignora

            int idxP = pixelToIndex[p];
            
            if (idxP != prefilter.currentFZ) prefilter.reset(idxP); // reset do cache só quando muda a FZ base

            for (int q : adj->getNeighborPixelsForward(p)) {
                if (img[q] == img[p]) continue; // mesma flatone: ignora
                int idxQ = pixelToIndex[q];
                
                
                if (prefilter.contains(idxQ)) continue;  // já emitido para esta FZ
                prefilter.insert(idxQ);
                
                // adicionar as aresta em uma direção
                if (idxP < idxQ) {
                    listOfAdjacentFlatZones[idxP].appendUnchecked(indexToPixel[idxQ]);
                } else {
                    listOfAdjacentFlatZones[idxQ].appendUnchecked(indexToPixel[idxP]);
                }
            }
        }
        // espelha as areasta para outra direção
        AdjacentFlatZonesSet::mirrorAndFinalize(listOfAdjacentFlatZones, pixelView);        
        
    }

    ImageUInt8Ptr getImage() { 
        return imgPtr; 
    }

    AdjacencyRelationPtr getAdjacencyRelation() { 
        return adj; 
    }
    
    std::shared_ptr<PixelSetManager> getPixelSetManager() { 
        return pixelBuffer; 
    }
    
    float averageDegree() const {
        size_t totalDegree = 0;
        int numFZ = getNumActiveFlatZones();
        for (int i = 0; i < getNumFlatZones(); ++i) {
            totalDegree += listOfAdjacentFlatZones[i].size();
        }
        return static_cast<float>(totalDegree) / static_cast<float>(numFZ);
    }

    int getNumEdges() const {
        size_t totalDegree = 0;
        for (int i = 0; i < getNumFlatZones(); ++i) {
            totalDegree += listOfAdjacentFlatZones[i].size();
        }
        return static_cast<int>(totalDegree / 2); // cada aresta contada duas vezes
    }
    int getNumActiveFlatZones() const {
        int cnt = 0;
        for (int idx = 0; idx < (int)parent.size(); ++idx) {
            if (parent[idx] != idx) continue; // raiz DSU
            cnt++;
        }
        return cnt;
    }

    int getNumFlatZones() const { 
        return pixelBuffer->numSets();
    }

    int getNumPixelInFlatzone(int rep){ 
        return pixelBuffer->numPixelsInSet(rep);
    }
    
    int findRepresentative(int rep) {
        int idx = pixelView.pixelToIndex[rep];
        return pixelView.indexToPixel[find(idx)];
    }

    /**
     * Retorna a lista de adjacências de uma dada flatzone.
     * @param rep da flatzone para a qual se deseja obter as adjacências.
     * @return Referência para o conjunto de flatzones adjacentes.
     * @throws std::out_of_range se o flatZoneID não existir no grafo
     */
    AdjacentFlatZones& getAdjacentFlatzonesFromPixel(int rep) {
        return listOfAdjacentFlatZones[pixelView.pixelToIndex[rep]];
    }

    template <class Func>
    void forEachAdjacentFlatzoneFromPixel(int rep, Func&& func) {
        int idx = pixelView.pixelToIndex[rep];
        assert(idx >= 0);
        int rootIdx = find(idx);
        for (int neighRep : listOfAdjacentFlatZones[rootIdx]) {
            func(neighRep);
        }
    }

    template <class Func>
    inline void forEachAdjacentFlatzoneFromPixelStatic(int rep, Func&& func) {
        int idx = pixelView.pixelToIndex[rep];
        assert(idx >= 0);
        for (int neighRep : listOfAdjacentFlatZones[idx]) {
            func(neighRep);
        }
    }
    
    /**
     * Verifica se uma flatzone é adjacente a outras.
     * @param repBase representante da flatzone a ser verificada.
     * @param repsFlatzones lista de flatzones a serem verificadas.
     * @return true se a flatzone for adjacente a alguma flazone de repsFlatzones, false caso contrário.
     */
    bool isAdjacent(int repBase, std::vector<int>& repsFlatzones) {
        // Obtém o conjunto de flatzones adjacentes ao nodeFlatZoneID
        AdjacentFlatZones& adjacentFlatzones = listOfAdjacentFlatZones[pixelView.pixelToIndex[repBase]];
    
        // Percorre os IDs das flatzones do parent para verificar se há alguma adjacente
        for (int rep: repsFlatzones) {
            if (adjacentFlatzones.find(rep) != adjacentFlatzones.end()) {
                return true;  
            }
        }
    
        return false; 
    }

    

    /**
     * @brief Funde um conjunto de flatzones adjacentes em torno de um representante base.
     *
     * Dado um representante de flatzone (repFlatzone) e um conjunto de representantes (repsFlatzones),
     * este método identifica quais elementos de repsFlatzones são realmente adjacentes a repFlatzone,
     * determina o representante vencedor (menor pixel entre repFlatzone e os candidatos adjacentes)
     * e realiza a fusão dessas flatzones no grafo:
     *   - Atualiza as relações de adjacência (removendo os arcos dos perdedores e ligando vizinhos ao vencedor);
     *   - Atualiza a estrutura union-find para refletir a unificação;
     *   - Integra o conteúdo (predessores/pixels) das flatzones perdedoras no vencedor.
     *
     * Após a fusão, repsFlatzones é modificada para remover os elementos realmente fundidos,
     * garantindo a presença do vencedor.
     *
     * @param repFlatzone Representante base da flatzone a partir do qual a fusão é feita.
     * @param repsFlatzones Lista de representantes candidatos à fusão (pode ser modificada).
     * @return O representante vencedor (ID) da flatzone unificada.
     */
    int mergeAdjacentCandidatesInPlace(int repFlatzone, std::vector<int>& repsFlatzones) {
        // 1) Índice e adjacências do repFlatzone
        int idxFlatzone = pixelView.pixelToIndex[repFlatzone];
        const AdjacentFlatZones& adjFlatzone = listOfAdjacentFlatZones[idxFlatzone];

        // 2) Filtrar apenas candidatos realmente adjacentes e escolher o vencedor (menor pixel)
        tmpCandidatesSize = 0;
        ensureTmpCandidatesCapacity((int)repsFlatzones.size() + 1);
        int winnerRep = repFlatzone;
        for (int rep : repsFlatzones) {
            if (adjFlatzone.find(rep) != adjFlatzone.end()) {
                tmpCandidatesBuffer[tmpCandidatesSize++] = rep;
                if (rep < winnerRep) winnerRep = rep;
            }
        }

        // Inclui a própria base se ela perder para o winner
        if (repFlatzone != winnerRep) {
            tmpCandidatesBuffer[tmpCandidatesSize++] = repFlatzone;
        }else{
            repsFlatzones.push_back(winnerRep); // garante que o winnerRep esteja na lista
        }
        int idxWinner = pixelView.pixelToIndex[winnerRep];

        // 3) Rewire das adjacências (cada loser -> winner), simétrico e sem self-loop
        for (int i = 0; i < tmpCandidatesSize; ++i) {
            int loserRep = tmpCandidatesBuffer[i];
            if (loserRep == winnerRep) continue;
            int idxL = pixelView.pixelToIndex[loserRep];

            AdjacentFlatZones& adjLoser = listOfAdjacentFlatZones[idxL];

            for (int neighRep : adjLoser) {
                if (neighRep == winnerRep) continue; // evita self-loop
                int idxN = pixelView.pixelToIndex[neighRep];

                // winner <-> neighbor
                listOfAdjacentFlatZones[idxWinner].insert(neighRep);
                listOfAdjacentFlatZones[idxN].insert(winnerRep);

                // remove loser <-> neighbor
                listOfAdjacentFlatZones[idxN].erase(loserRep);
            }

            // remover arco winner <-> loser, se existir
            listOfAdjacentFlatZones[idxWinner].erase(loserRep);

            // limpa adjacências do perdedor
            AdjacentFlatZones().swap(adjLoser);

            // atualiza union-find
            unite(idxWinner, idxL);
        }

        // 4) Atualiza repsFlatzones:
        for (int i = 0; i < tmpCandidatesSize; ++i) {
            int x = tmpCandidatesBuffer[i];
            if (x != winnerRep) std::erase(repsFlatzones, x);
        }
        
        return winnerRep;
    }
    
    /**
     * @brief Funde um conjunto de flatzones base e seus candidatos adjacentes em uma única flat zone.
     *
     * Dado um conjunto de representantes de flatzones base (baseReps), que se sabe formar
     * um componente conexo, este método:
     *   - Determina o representante vencedor (menor pixel entre todas as bases e candidatos);
     *   - Realiza a fusão de todas as flatzones base em torno do vencedor;
     *   - Identifica e integra os candidatos realmente adjacentes a qualquer uma das bases;
     *   - Atualiza as relações de adjacência no grafo (removendo arcos dos perdedores e ligando vizinhos ao vencedor);
     *   - Atualiza a estrutura union-find para refletir a unificação;
     *   - Integra os ciclos de predecessores/pixels das flatzones fundidas.
     *
     * Ao final, todas as flatzones base e suas adjacentes são unificadas em torno do representante vencedor.
     *
     * @param baseReps Lista de representantes que formam um único componente conexo (não modificada).
     * @param repsFlatzones Lista de representantes candidatos à fusão (pode ser modificada).
     * @param winnerRep (opcional) Representante vencedor (menor pixel) já foi computado; se -1, o menor pixel será computado.
     * @return O representante vencedor (ID) da flat zone unificada.
     */
    int mergeBasesWithAdjacentCandidatesInPlace(const std::vector<int>& baseReps, std::vector<int>& repsFlatzones, int winnerRep = -1){

        // 1) vencedor = menor representante
        if (winnerRep == -1)
            winnerRep = *std::min_element(baseReps.begin(), baseReps.end());
        int idxWinner = pixelView.pixelToIndex[winnerRep];

        // 2) funde todas as bases diretamente no winner
        for (int rep : baseReps) {
            if (rep == winnerRep) continue;

            int idxL = pixelView.pixelToIndex[rep];
            AdjacentFlatZones& adjL = listOfAdjacentFlatZones[idxL];

            // rewire: neighbors(loser) -> winner (simétrico, sem self-loop)
            for (int neighRep : adjL) {
                if (neighRep == winnerRep) continue;
                int idxN = pixelView.pixelToIndex[neighRep];

                listOfAdjacentFlatZones[idxWinner].insert(neighRep);
                listOfAdjacentFlatZones[idxN].insert(winnerRep);
                listOfAdjacentFlatZones[idxN].erase(rep);
            }

            // remove arco winner <-> loser (se existir) e limpa adj do loser
            listOfAdjacentFlatZones[idxWinner].erase(rep);
            AdjacentFlatZones().swap(adjL);

            // UF + predecessores
            unite(idxWinner, idxL);
        }

        // 3) agora anexa candidatos adjacentes ao winner
        return mergeAdjacentCandidatesInPlace(winnerRep, repsFlatzones);
    }

    
    //Um único representante
    auto getPixelsByFlatzone(int rep) {
        return PixelSetManager::PixelsBySetRange<std::array<int,1>>(pixelView, std::array<int,1>{rep});
    }

    // Qualquer range (vector<int>, span<const int>, RepsOfCCRange, …)
    template<class Range>
    auto getPixelsByFlatzone(Range reps) {
        return PixelSetManager::PixelsBySetRange<Range>(pixelView, std::move(reps));
    }

    auto getFlatzoneRepresentatives() const {
        return pixelBuffer->getFlatzoneRepresentatives();
    }


    void printUnionFind() {
        std::cout << "===== Union-Find State =====" << std::endl;

        for (size_t id = 0; id < parent.size(); ++id) {
            int rep = pixelView.indexToPixel[find(id)];
            std::cout << "Flatzone " << std::setw(4) << pixelView.indexToPixel[id] 
                    << " -> Representative " << std::setw(4) << rep
                    << std::endl;
        }
    }

    void printUnionFindGroups()  {
        std::unordered_map<int, std::vector<int>> groups;

        for (size_t id = 0; id < parent.size(); ++id) {
            int rep = pixelView.indexToPixel[find(id)];
            groups[rep].push_back( pixelView.indexToPixel[id] );
        }

        std::cout << "===== Union-Find Groups =====" << std::endl;
        for (auto& [rep, members] : groups) {
            std::cout << "Rep " << rep << ": ";
            for (int id : members) {
                std::cout << id << " ";
            }
            std::cout << std::endl;
        }
    }

    /*
    struct FZGraphCheckReport {
        long long activeRoots = 0;       // |R|
        long long halfEdges   = 0;       // sum_{r in R} |Adj(r)|
        long long selfLoops   = 0;       // count of (r in Adj(r))
        long long nonRootNeigh = 0;      // count of neighbors that are not DSU roots
        long long duplicateNeigh = 0;    // duplicates inside Adj(r)
        long long asymEdges   = 0;       // count of directed adjacencies lacking reverse
        long long invalidNeigh = 0;      // neighbor rep not mapped / out of range (defensive)
        
        bool ok() const {
            return selfLoops == 0 && nonRootNeigh == 0 && duplicateNeigh == 0 && asymEdges == 0 && invalidNeigh == 0;
        }
    };


    **
     * @brief Valida invariantes do grafo após merges.
     *
     * @param checkSymmetry   se true, testa u in Adj(r) => r in Adj(u)
     * @param checkDuplicates se true, testa duplicatas em Adj(r)
     * @param maxRootsToCheck se >0, amostra no máximo esse número de raízes (debug rápido)
     *
    FZGraphCheckReport checkConsistency(bool checkSymmetry = true,
                                        bool checkDuplicates = true,
                                        int maxRootsToCheck = 0) const
    {
        FZGraphCheckReport rep;

        const int n = (int)parent.size();
        if (n == 0) return rep;

        // stamp para detectar duplicatas por raiz: marca por índice do DSU (vizinho)
        static thread_local std::vector<uint32_t> stamp;
        static thread_local uint32_t token = 1;

        if ((int)stamp.size() != n) { stamp.assign(n, 0); token = 1; }
        if (++token == 0) { std::fill(stamp.begin(), stamp.end(), 0); token = 1; }

        int checkedRoots = 0;

        for (int idxR = 0; idxR < n; ++idxR) {
            if (parent[idxR] != idxR) continue; // não é root (loser)
            rep.activeRoots++;
            checkedRoots++;
            if (maxRootsToCheck > 0 && checkedRoots > maxRootsToCheck) break;

            const int repR = pixelView.indexToPixel[idxR];
            // se o índice->pixel for inválido, pule (ou conte, dependendo da sua modelagem)
            if (repR < 0) continue;

            const AdjacentFlatZones& adjR = listOfAdjacentFlatZones[idxR];
            rep.halfEdges += (long long)adjR.size();

            // token exclusivo por raiz (para duplicatas)
            uint32_t tokR = token;
            if (checkDuplicates) {
                // usa um novo token por raiz (se token overflow, já tratamos acima)
                if (++token == 0) { std::fill(stamp.begin(), stamp.end(), 0); token = 1; }
                tokR = token;
            }

            for (int neighPixel : adjR) {
                if (neighPixel == repR) {
                    rep.selfLoops++;
                    continue;
                }

                const int idxN = repIndexFromPixel_(neighPixel);
                if (idxN < 0 || idxN >= n) {
                    rep.invalidNeigh++;
                    continue;
                }

                // vizinho deve ser raiz DSU (canônico)
                if (parent[idxN] != idxN) {
                    rep.nonRootNeigh++;
                }

                if (checkDuplicates) {
                    // marca por índice do vizinho (DSU index)
                    if (stamp[idxN] == tokR) rep.duplicateNeigh++;
                    else stamp[idxN] = tokR;
                }

                if (checkSymmetry) {
                    // simetria só faz sentido se o vizinho for um root "válido"
                    // (se não for root, a assimetria já está explicada por nonRootNeigh)
                    if (parent[idxN] == idxN) {
                        const AdjacentFlatZones& adjN = listOfAdjacentFlatZones[idxN];

                        // AdjacentFlatZonesSet tem find() no seu Common.hpp.
                        // Se o seu tipo não tiver, troque por busca linear.
                        if (adjN.find(repR) == adjN.end()) {
                            rep.asymEdges++;
                        }
                    }
                }
            }
        }

        return rep;
    }*/

};


#endif // FLATZONEGRAPH_HPP
