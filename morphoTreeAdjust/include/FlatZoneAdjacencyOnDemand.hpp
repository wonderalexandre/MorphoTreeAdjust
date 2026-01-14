#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include "Common.hpp"
#include "AdjacencyRelation.hpp"

/**
 * @brief Estrutura de adjacência de flatzones sob demanda (sem grafo explícito).
 *
 * Esta classe mantém apenas os vértices (flatzones) e calcula adjacências
 * quando solicitado. A base do funcionamento é:
 * - DSU (union-find) para rastrear merges e encontrar representantes canônicos;
 * - PixelSetManager para percorrer os pixels de cada flatzone;
 * - lista ligada de pixels de borda por flatzone (borderHead/borderNext/borderTail),
 *   evitando varrer pixels internos nas consultas.
 *
 * Fluxo de construção:
 * - BFS por componente conexo (mesmo nível de cinza) cria cada flatzone;
 * - pixels de borda são detectados (vizinho com nível diferente) e inseridos na
 *   lista ligada da flatzone.
 *
 * Consultas de adjacência:
 * - percorre apenas pixels de borda da flatzone base;
 * - para cada vizinho de pixel, obtém o root via DSU;
 * - deduplica roots com a estrutura de "stamp".
 *
 * Merges:
 * - une conjuntos no DSU;
 * - faz splice O(1) das listas circulares de pixels;
 * - concatena listas de borda em O(1) e marca a borda como "suja".
 *
 * Borda suja e refiltro parcial:
 * - após merges, a lista de borda pode conter pixels que viraram internos;
 * - a limpeza é incremental: a cada consulta, percorre-se apenas um número limitado
 *   de pixels (budget) removendo os que não são mais borda;
 * - quando o cursor alcança o fim da lista, a borda volta a ser considerada limpa.
 *
 * Complexidade típica:
 * - merge: O(1) amortizado (DSU + concatenação de borda);
 * - consulta: O(|borda| * grau), com deduplicação por stamp;
 * - refiltro parcial: O(budget * grau) por consulta.
 *
 * Invariantes:
 * - borderNext forma uma lista ligada sem ciclos por flatzone;
 * - borderCount reflete o tamanho atual da lista de borda;
 * - indexToPixel == -1 indica slot inativo após merge.
 */
class FlatZoneAdjacencyOnDemand {
public:
    /**
     * @brief Construtor de cópia padrão.
     */
    FlatZoneAdjacencyOnDemand(const FlatZoneAdjacencyOnDemand&) = default;
    /**
     * @brief Atribuição de cópia padrão.
     */
    FlatZoneAdjacencyOnDemand& operator=(const FlatZoneAdjacencyOnDemand&) = default;

    /**
     * @brief Constrói a estrutura a partir de uma imagem e raio de adjacência.
     */
    FlatZoneAdjacencyOnDemand(ImageUInt8Ptr imagePtr, int radiusAdj)
        : FlatZoneAdjacencyOnDemand(imagePtr,
              std::make_shared<AdjacencyRelation>(imagePtr->getNumRows(),
                                                  imagePtr->getNumCols(),
                                                  radiusAdj)) {}

    /**
     * @brief Constrói a estrutura a partir de uma imagem e relação de adjacência.
     */
    FlatZoneAdjacencyOnDemand(ImageUInt8Ptr imagePtr, AdjacencyRelationPtr adjRelation)
        : imgPtr(std::move(imagePtr)),
          adj(std::move(adjRelation)),
          pixelBuffer(std::make_shared<PixelSetManager>(imgPtr->getSize())),
          pixelView(pixelBuffer->view()) {

        // Cria vértices (flatzones) a partir da imagem.
        buildFlatZones();
    }

    /**
     * @brief Retorna a imagem de entrada.
     */
    ImageUInt8Ptr getImage() { return imgPtr; }
    /**
     * @brief Retorna a relação de adjacência.
     */
    AdjacencyRelationPtr getAdjacencyRelation() { return adj; }
    /**
     * @brief Retorna o PixelSetManager usado nas flatzones.
     */
    std::shared_ptr<PixelSetManager> getPixelSetManager() { return pixelBuffer; }

    /**
     * @brief Número total de flatzones ativas (componentes DSU ativos).
     * Aqui usamos o número de "slots" criados (pixelBuffer->numSets()) e
     * contamos slots válidos (indexToPixel != -1) para refletir merges.
     */
    int getNumFlatZones() const {
        int acc = 0;
        const auto& idxToPix = pixelBuffer->indexToPixel;
        for (int p : idxToPix) if (p != -1) ++acc;
        return acc;
    }


    /**
     * @brief Retorna grau médio (indisponível no modo on-demand).
     */
    float averageDegree() const {
        return -1;
    }

    /**
     * @brief Retorna número de arestas (indisponível no modo on-demand).
     */
    int getNumEdges() const {
        return -1;
    }
    /**
     * @brief Retorna número de raízes DSU ativas (flatzones atuais).
     */
    int getNumActiveFlatZones() const {
        int cnt = 0;
        for (int idx = 0; idx < (int)parent.size(); ++idx) {
            // Apenas raízes DSU contam como flatzones ativas.
            if (parent[idx] != idx) continue;
            cnt++;
        }
        return cnt;
    }

    /**
     * @brief DSU: encontra o representante canônico (pixel-id) da flatzone
     * que contém o pixel/rep informado.
     *
     * Observação: "rep" pode ser qualquer pixel dentro da flatzone.
     */
    int findRepresentative(int rep) {
        int idx = pixelView.pixelToIndex[rep];
        if (idx < 0) return rep;
        int root = findIndex(idx);
        return pixelView.indexToPixel[root];
    }

    /**
     * @brief Merge de duas flatzones (por representantes em pixel-id).
     * Não rewire adjacências (não existem arestas armazenadas).
     *
     * @return pixel-id do vencedor (representante canônico).
     */
    int mergeFlatZones(int repA, int repB) {
        int idxA = pixelView.pixelToIndex[repA];
        int idxB = pixelView.pixelToIndex[repB];
        if (idxA < 0 || idxB < 0) return findRepresentative(repA);

        int rootA = findIndex(idxA);
        int rootB = findIndex(idxB);
        if (rootA == rootB) return pixelView.indexToPixel[rootA];

        // União por tamanho (usa sizeSets do PixelSetManager).
        int repRootA = pixelView.indexToPixel[rootA];
        int repRootB = pixelView.indexToPixel[rootB];

        int sizeA = pixelView.sizeSets[rootA];
        int sizeB = pixelView.sizeSets[rootB];

        int winnerRoot = rootA;
        int loserRoot  = rootB;
        int winnerRep  = repRootA;
        int loserRep   = repRootB;

        if (sizeB > sizeA) {
            winnerRoot = rootB; loserRoot = rootA;
            winnerRep  = repRootB; loserRep  = repRootA;
        }

        // Aponta o perdedor para o vencedor no DSU.
        parent[loserRoot] = winnerRoot;

        // Splice O(1) das listas circulares e invalida o slot perdedor.
        pixelBuffer->mergeSetsByRep(winnerRep, loserRep);

        // Concatena listas de borda e marca para refiltro incremental.
        concatBorderLists(winnerRoot, loserRoot);
        borderDirty[winnerRoot] = 1u;
        resetBorderRefilterCursor(winnerRoot);

        return winnerRep;
    }

    /**
     * @brief Itera vizinhos de rep e chama func(repVizinho).
     */
    template <class Func>
    void forEachAdjacentFlatzoneFromPixel(int rep, Func&& func) {
        int baseIdxRaw = pixelView.pixelToIndex[rep];
        int baseIdx = -1;
        int baseRep = rep;
        if (baseIdxRaw >= 0 && parent[baseIdxRaw] == baseIdxRaw) {
            baseIdx = baseIdxRaw;
        } else {
            // Normaliza para o representante canônico.
            baseRep = findRepresentative(rep);
            baseIdxRaw = pixelView.pixelToIndex[baseRep];
            assert(baseIdxRaw >= 0);
            baseIdx = findIndex(baseIdxRaw);
        }

        if (borderDirty[baseIdx]) {
            // Limpa parcialmente a borda antes de consultar.
            partialRefilterStep(baseIdx, computePartialRefilterBudget(baseIdx));
            if (borderCursor[baseIdx] == -1) {
                borderDirty[baseIdx] = 0u;
            }
        }

        beginToken();
        forEachBorderPixel(baseIdx, [&](int p) {
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0) continue;
                int rootQ = (parent[idxQ] == idxQ) ? idxQ : findIndex(idxQ);
                if (rootQ == baseIdx) continue;

                // Dedup: evita repetir o mesmo root na mesma consulta.
                if (stamp[rootQ] != token) {
                    stamp[rootQ] = token;
                    int neighRep = pixelView.indexToPixel[rootQ];
                    if (neighRep != -1) func(neighRep);
                }
            }
        });
    }

    /**
     * @brief Itera vizinhos de rep sem DSU find (assume ausência de merges).
     */
    template <class Func>
    inline void forEachAdjacentFlatzoneFromPixelStatic(int rep, Func&& func) {
        const int baseIdx = pixelView.pixelToIndex[rep];
        assert(baseIdx >= 0);

        beginToken();
        forEachBorderPixel(baseIdx, [&](int p) {
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0 || idxQ == baseIdx) continue;

                // Dedup: índice já é raiz (sem find).
                if (stamp[idxQ] != token) {
                    stamp[idxQ] = token;
                    int neighRep = pixelView.indexToPixel[idxQ];
                    if (neighRep != -1) func(neighRep);
                }
            }
        });
    }

    /**
     * @brief Teste direto de adjacência entre duas flatzones usando pixels.
     * Implementado varrendo a menor flatzone (por sizeSets).
     */
    bool areAdjacentByPixels(int repA, int repB) {
        int a = findRepresentative(repA);
        int b = findRepresentative(repB);
        if (a == b) return false;

        int idxA = findIndex(pixelView.pixelToIndex[a]);
        int idxB = findIndex(pixelView.pixelToIndex[b]);

        int sizeA = pixelView.sizeSets[idxA];
        int sizeB = pixelView.sizeSets[idxB];

        int smallRep = a;
        int otherIdx = idxB;

        if (sizeB < sizeA) {
            smallRep = b; 
            otherIdx = idxA;
        }

        int p = smallRep;
        bool first = true;
        while (first || p != smallRep) {
            first = false;
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0) continue;
                int rootQ = (parent[idxQ] == idxQ) ? idxQ : findIndex(idxQ);
                if (rootQ == otherIdx) return true;
            }
            p = pixelView.pixelsNext[p];
            if (p < 0) break;
        }

        return false;
    }

    /**
     * @brief Versão "naive" do merge de candidatos adjacentes:
     * filtra candidatos adjacentes (via pixels) e faz unions no DSU.
     *
     * Mantém a mesma assinatura do seu FlatZonesGraph eager.
     */
    int mergeAdjacentCandidatesInPlace(int repFlatzone, std::vector<int>& repsFlatzones) {
        const int baseCanon = findRepresentative(repFlatzone);
        const int baseIdx = markAdjacentRootsFromRep(baseCanon);
        if (baseIdx < 0) return baseCanon;

        tmpCandidatesSize = 0;
        ensureTmpCandidatesCapacity((int)repsFlatzones.size() + 1);

        int winnerRep = baseCanon;

        // Canoniza e filtra quem realmente é vizinho.
        for (int r : repsFlatzones) {
            int c = findRepresentative(r);
            if (c == baseCanon) continue;
            int idxC = pixelView.pixelToIndex[c];
            if (idxC < 0) continue;
            int rootC = findIndex(idxC);
            if (stamp[rootC] != token) continue;
            tmpCandidatesBuffer[tmpCandidatesSize++] = c;
            if (c < winnerRep) winnerRep = c;
        }

        // Garante que a base participe se não for o vencedor.
        if (winnerRep != baseCanon) tmpCandidatesBuffer[tmpCandidatesSize++] = baseCanon;

        // Faz merges: todos -> vencedor.
        for (int i = 0; i < tmpCandidatesSize; ++i) {
            int loser = tmpCandidatesBuffer[i];
            if (loser == winnerRep) continue;
            mergeFlatZones(winnerRep, loser);
        }

        // Atualiza repsFlatzones: remove perdedores colapsados no vencedor.
        const int wCanon = findRepresentative(winnerRep);
        repsFlatzones.erase(std::remove_if(repsFlatzones.begin(), repsFlatzones.end(),
                                           [&](int x){ return findRepresentative(x) == wCanon; }),
                            repsFlatzones.end());
        repsFlatzones.push_back(wCanon);
        return wCanon;
    }

    /**
     * @brief Funde um conjunto de bases no winner e depois tenta anexar candidatos adjacentes.
     * Assinatura compatível com a versão eager.
     */
    int mergeBasesWithAdjacentCandidatesInPlace(const std::vector<int>& baseReps,
                                               std::vector<int>& repsFlatzones,
                                               int winnerRep = -1) {
        if (baseReps.empty()) return -1;

        if (winnerRep == -1) {
            winnerRep = *std::min_element(baseReps.begin(), baseReps.end());
        }
        winnerRep = findRepresentative(winnerRep);

        for (int r : baseReps) {
            int c = findRepresentative(r);
            if (c == winnerRep) continue;
            mergeFlatZones(winnerRep, c);
        }

        return mergeAdjacentCandidatesInPlace(winnerRep, repsFlatzones);
    }

private:
    ImageUInt8Ptr imgPtr;
    AdjacencyRelationPtr adj;

    std::shared_ptr<PixelSetManager> pixelBuffer;
    PixelSetManager::View pixelView;

    // DSU em índices de flatzones.
    std::vector<int> parent;

    // Borda por flatzone (lista encadeada por pixel).
    std::vector<int> borderHead;
    std::vector<int> borderTail;
    // Por pixel, aponta o próximo na lista de borda.
    std::vector<int> borderNext;
    std::vector<int> borderCount;
    std::vector<int> borderCursor;
    std::vector<int> borderCursorPrev;
    std::vector<uint8_t> borderDirty;

    std::vector<int> tmpCandidatesBuffer;
    int tmpCandidatesSize{0};

    // Deduplicação por índice-root (stamp).
    std::vector<uint32_t> stamp;
    uint32_t token{1};

private:
    /**
     * @brief Encontra a raiz do DSU com compressão de caminho.
     */
    int findIndex(int i) {
        int root = i;
        while (parent[root] != root) root = parent[root];

        // Compressão de caminho para acelerar próximas buscas.
        while (parent[i] != i) {
            int next = parent[i];
            parent[i] = root;
            i = next;
        }
        return root;
    }

    /**
     * @brief Avança a geração de stamp e limpa ao estourar.
     */
    void beginToken() {
        if (++token == 0u) {
            std::fill(stamp.begin(), stamp.end(), 0u);
            token = 1u;
            std::cout << "  [FlatZoneAdjacencyOnDemand] Stamp reset.\n";
        }
    }

    /**
     * @brief Marca roots adjacentes a partir do representante base.
     */
    int markAdjacentRootsFromRep(int baseRep) {
        int baseIdxRaw = pixelView.pixelToIndex[baseRep];
        if (baseIdxRaw < 0) return -1;
        int baseIdx = -1;
        if (parent[baseIdxRaw] == baseIdxRaw) {
            baseIdx = baseIdxRaw;
        } else {
            baseIdx = findIndex(baseIdxRaw);
        }

        beginToken();

        // Percorre apenas pixels de borda e coleta vizinhos.
        forEachBorderPixel(baseIdx, [&](int p) {
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0) continue;
                int rootQ = (parent[idxQ] == idxQ) ? idxQ : findIndex(idxQ);
                if (rootQ == baseIdx) continue;

                // Dedup por stamp para não repetir roots.
                if (stamp[rootQ] != token) {
                    stamp[rootQ] = token;
                }
            }
        });

        return baseIdx;
    }

    /**
     * @brief Constrói flatzones e listas de borda a partir da imagem.
     */
    void buildFlatZones() {
        int numPixels = imgPtr->getSize();
        const auto* imageData = imgPtr->rawData();

        std::vector<uint8_t> visited(numPixels, 0);
        FastQueue<int> queue(numPixels / 4);

        auto& [pixelToIndex, indexToPixel, sizeSets, pixelsNext] = *pixelBuffer;
        int numFZ = 0;

        // Estruturas de borda (tamanho máximo: numPixels).
        borderHead.assign(numPixels, -1);
        borderTail.assign(numPixels, -1);
        borderNext.assign(numPixels, -1);
        borderCount.assign(numPixels, 0);

        // Durante o build, preenche pixelToIndex e constrói pixelsNext.
        for (int p = 0; p < numPixels; ++p) {
            if (visited[p]) continue;

            int tail = p;
            int sizeFlatzone = 0;
            // Índice desta flatzone.
            int idxFZ = numFZ++;
            queue.push(p);
            visited[p] = 1;
            pixelToIndex[p] = idxFZ;
            indexToPixel[idxFZ] = p;

            // BFS na componente conexa de mesmo nível.
            while (!queue.empty()) {
                int q = queue.pop();
                sizeFlatzone++;
                bool hasDiff = false;

                for (int nq : adj->getNeighborPixels(q)) {
                    if (!visited[nq] && imageData[nq] == imageData[p]) {
                        visited[nq] = 1;
                        queue.push(nq);
                        pixelToIndex[nq] = idxFZ;
                        pixelsNext[tail] = nq;
                        tail = nq;
                    } else if (imageData[nq] != imageData[p]) {
                        hasDiff = true;
                    }
                }

                // Adiciona pixel à borda se houver vizinho com nível diferente.
                if (hasDiff) appendBorderPixel(idxFZ, q);
            }

            // Fecha ciclo da lista circular de pixels.
            pixelsNext[tail] = p;
            sizeSets[idxFZ] = sizeFlatzone;
        }

        // Reduz vetores para numFZ (economia de memória).
        pixelBuffer->shrinkToNumSets(numFZ);
        pixelView = pixelBuffer->view();

        borderHead.resize(numFZ);
        borderTail.resize(numFZ);
        borderCount.resize(numFZ);
        borderCursor.resize(numFZ);
        borderCursorPrev.resize(numFZ);
        borderDirty.resize(numFZ, 0u);
        for (int i = 0; i < numFZ; ++i) {
            borderCursor[i] = borderHead[i];
            borderCursorPrev[i] = -1;
        }

        // DSU inicial: cada set é raiz.
        parent.resize(numFZ);
        std::iota(parent.begin(), parent.end(), 0);

        // Stamps por índice de set (dedup em consultas).
        stamp.assign(numFZ, 0u);
        token = 1u;

        tmpCandidatesBuffer.resize(numFZ);
        tmpCandidatesSize = 0;
    }

    /**
     * @brief Garante capacidade do buffer temporário de candidatos.
     */
    void ensureTmpCandidatesCapacity(int required) {
        if (required > (int)tmpCandidatesBuffer.size()) {
            tmpCandidatesBuffer.resize(required);
        }
    }

    /**
     * @brief Anexa um pixel à lista de borda de uma raiz.
     */
    void appendBorderPixel(int rootIdx, int p) {
        if (borderHead[rootIdx] == -1) {
            borderHead[rootIdx] = borderTail[rootIdx] = p;
        } else {
            borderNext[borderTail[rootIdx]] = p;
            borderTail[rootIdx] = p;
        }
        // Garante terminação da lista ligada.
        borderNext[p] = -1;
        borderCount[rootIdx] += 1;
    }

    /**
     * @brief Concatena listas de borda após um merge.
     */
    void concatBorderLists(int winnerRoot, int loserRoot) {
        const int loserHead = borderHead[loserRoot];
        if (loserHead == -1) {
            borderCount[loserRoot] = 0;
            borderHead[loserRoot] = -1;
            borderTail[loserRoot] = -1;
            borderCursor[loserRoot] = -1;
            borderCursorPrev[loserRoot] = -1;
            borderDirty[loserRoot] = 0u;
            return;
        }

        // Concatenação O(1): liga a cauda do vencedor ao início do perdedor.
        if (borderHead[winnerRoot] == -1) {
            borderHead[winnerRoot] = loserHead;
            borderTail[winnerRoot] = borderTail[loserRoot];
        } else {
            borderNext[borderTail[winnerRoot]] = loserHead;
            borderTail[winnerRoot] = borderTail[loserRoot];
        }

        borderHead[loserRoot] = -1;
        borderTail[loserRoot] = -1;
        borderCount[winnerRoot] += borderCount[loserRoot];
        borderCount[loserRoot] = 0;
        borderCursor[loserRoot] = -1;
        borderCursorPrev[loserRoot] = -1;
        borderDirty[loserRoot] = 0u;
    }

    /**
     * @brief Reseta o cursor de refiltro parcial de uma raiz.
     */
    void resetBorderRefilterCursor(int rootIdx) {
        borderCursor[rootIdx] = borderHead[rootIdx];
        borderCursorPrev[rootIdx] = -1;
    }

    /**
     * @brief Itera todos os pixels de borda de uma raiz.
     */
    template <class Func>
    void forEachBorderPixel(int rootIdx, Func&& func) {
        int p = borderHead[rootIdx];
        while (p != -1) {
            int next = borderNext[p];
            func(p);
            p = next;
        }
    }

    /**
     * @brief Refiltra parcialmente a borda até o limite de budget.
     */
    void partialRefilterStep(int rootIdx, int budget) {
        if (budget <= 0) return;
        int current = borderCursor[rootIdx];
        int prev = borderCursorPrev[rootIdx];
        if (current == -1) {
            current = borderHead[rootIdx];
            prev = -1;
        }

        int processed = 0;
        while (current != -1 && processed < budget) {
            int next = borderNext[current];
            // Remove pixels que deixaram de ser borda.
            if (!isBorderPixelForRoot(current, rootIdx)) {
                if (prev == -1) {
                    borderHead[rootIdx] = next;
                } else {
                    borderNext[prev] = next;
                }
                if (borderTail[rootIdx] == current) {
                    borderTail[rootIdx] = prev;
                }
                borderNext[current] = -1;
                borderCount[rootIdx] -= 1;
            } else {
                prev = current;
            }
            current = next;
            processed += 1;
        }

        borderCursor[rootIdx] = current;
        borderCursorPrev[rootIdx] = prev;
        if (borderHead[rootIdx] == -1) {
            borderTail[rootIdx] = -1;
            borderCursor[rootIdx] = -1;
            borderCursorPrev[rootIdx] = -1;
        } else if (current == -1) {
            borderCursorPrev[rootIdx] = -1;
        }
    }

    /**
     * @brief Retorna true se o pixel ainda é borda para rootIdx.
     */
    bool isBorderPixelForRoot(int p, int rootIdx) {
        for (int q : adj->getNeighborPixels(p)) {
            int idxQ = pixelView.pixelToIndex[q];
            if (idxQ < 0) continue;
            int rootQ = (parent[idxQ] == idxQ) ? idxQ : findIndex(idxQ);
            if (rootQ != rootIdx) return true;
        }
        return false;
    }

    /**
     * @brief Computa quantos pixels de borda refiltrar por consulta.
     */
    int computePartialRefilterBudget(int rootIdx) const {
        int count = borderCount[rootIdx];
        if (count <= 0) return 0;
        // Heurística simples: controla custo por consulta vs. limpeza da borda.
        int budget = 0;
        if (count < 64) {
            budget = 8;
        } else if (count < 256) {
            budget = 16;
        } else if (count < 1024) {
            budget = count / 8;
        } else {
            budget = count / 4;
        }
        if (budget < 8) budget = 8;
        if (budget > 4096) budget = 4096;
        return budget;
    }


};
