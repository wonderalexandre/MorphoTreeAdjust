#pragma once


#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>
#include <utility>
#include <cmath>

#include "Common.hpp"
#include "AdjacencyRelation.hpp"

template <typename GraphT>
concept FlatZonesGraphCommonInterface = requires(GraphT& graph,
                                                 const GraphT& constGraph,
                                                 int rep,
                                                 std::vector<int>& repsFlatzones,
                                                 const std::vector<int>& baseReps) {
    { graph.getImage() } -> std::same_as<ImageUInt8Ptr>;
    { graph.getAdjacencyRelation() } -> std::same_as<AdjacencyRelationPtr>;
    { graph.getPixelSetManager() } -> std::same_as<std::shared_ptr<PixelSetManager>>;
    { constGraph.getNumFlatZones() } -> std::convertible_to<int>;
    { graph.findRepresentative(rep) } -> std::convertible_to<int>;
    { graph.forEachAdjacentFlatzoneFromPixel(rep, [](int){}) } -> std::same_as<void>;
    { graph.forEachAdjacentFlatzoneFromPixelStatic(rep, [](int){}) } -> std::same_as<void>;
    { graph.mergeAdjacentCandidatesInPlace(rep, repsFlatzones) } -> std::convertible_to<int>;
    { graph.mergeBasesWithAdjacentCandidatesInPlace(baseReps, repsFlatzones, rep) } -> std::convertible_to<int>;
};

class FlatZonesGraphBase {
protected:
    ImageUInt8Ptr imgPtr;
    AdjacencyRelationPtr adj;
    std::shared_ptr<PixelSetManager> pixelBuffer;
    PixelSetManager::View pixelView;
    std::vector<int> parent;

    FlatZonesGraphBase(ImageUInt8Ptr imagePtr, AdjacencyRelationPtr adjRelation)
        : imgPtr(std::move(imagePtr)),
          adj(std::move(adjRelation)),
          pixelBuffer(std::make_shared<PixelSetManager>(imgPtr->getSize())),
          pixelView(pixelBuffer->view()) {}

    int findIndex(int i) {
        int root = i;
        while (parent[root] != root) root = parent[root];
        while (parent[i] != i) {
            int next = parent[i];
            parent[i] = root;
            i = next;
        }
        return root;
    }

public:
    ImageUInt8Ptr getImage() { return imgPtr; }
    AdjacencyRelationPtr getAdjacencyRelation() { return adj; }
    std::shared_ptr<PixelSetManager> getPixelSetManager() { return pixelBuffer; }
    int getNumFlatZones() const { return pixelBuffer->numSets(); }

    int findRepresentative(int rep) {
        int idx = pixelView.pixelToIndex[rep];
        if (idx < 0) return rep;
        int root = findIndex(idx);
        return pixelView.indexToPixel[root];
    }
};

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
class FlatZonesGraphOnDemandEdgesByBoundary : public FlatZonesGraphBase {
public:
    /**
     * @brief Construtor de cópia padrão.
     */
    FlatZonesGraphOnDemandEdgesByBoundary(const FlatZonesGraphOnDemandEdgesByBoundary&) = default;
    /**
     * @brief Atribuição de cópia padrão.
     */
    FlatZonesGraphOnDemandEdgesByBoundary& operator=(const FlatZonesGraphOnDemandEdgesByBoundary&) = default;

    /**
     * @brief Constrói a estrutura a partir de uma imagem e raio de adjacência.
     */
    FlatZonesGraphOnDemandEdgesByBoundary(ImageUInt8Ptr imagePtr, int radiusAdj)
        : FlatZonesGraphOnDemandEdgesByBoundary(imagePtr,
              std::make_shared<AdjacencyRelation>(imagePtr->getNumRows(),
                                                  imagePtr->getNumCols(),
                                                  radiusAdj)) {}

    /**
     * @brief Constrói a estrutura a partir de uma imagem e relação de adjacência.
     */
    FlatZonesGraphOnDemandEdgesByBoundary(ImageUInt8Ptr imagePtr, AdjacencyRelationPtr adjRelation)
        : FlatZonesGraphBase(std::move(imagePtr), std::move(adjRelation)) {

        // Cria vértices (flatzones) a partir da imagem.
        buildFlatZones();
    }

private:

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

        int winnerRoot = rootA;
        int loserRoot  = rootB;
        int winnerRep  = repRootA;
        int loserRep   = repRootB;

        // Mantém o menor pixel como representante (propriedade de rep mínima).
        if (repRootB < repRootA) {
            winnerRoot = rootB;
            loserRoot  = rootA;
            winnerRep  = repRootB;
            loserRep   = repRootA;
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

public:

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
     * @brief Versão "naive" do merge de candidatos adjacentes:
     * filtra candidatos adjacentes (via pixels) e faz unions no DSU.
     *
     * Mantém a mesma assinatura do seu FlatZonesGraphFullEdges.
     */
    int mergeAdjacentCandidatesInPlace(int repFlatzone, std::vector<int>& repsFlatzones) {
        // Complexidade (on-demand): O(|Borda_base| * d + R * α(N)).
        // - Borda_base = numero de pixels de borda da flatzone base.
        // - d = grau da adjacencia (4/8 etc.).
        // - R = repsFlatzones.size(); α(N) ~ constante (DSU).
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
        // Complexidade (on-demand): O(B * α(N) + |Borda_winner| * d + R * α(N)).
        // - B = baseReps.size().
        if (baseReps.empty()) return -1;

        if (winnerRep == -1) {
            winnerRep = std::numeric_limits<int>::max();
            for (int r : baseReps) {
                int c = findRepresentative(r);
                if (c < winnerRep) {
                    winnerRep = c;
                }
            }
        } else {
            winnerRep = findRepresentative(winnerRep);
        }
        if (winnerRep == std::numeric_limits<int>::max()) {
            return -1;
        }

        for (int r : baseReps) {
            int c = findRepresentative(r);
            if (c == winnerRep) continue;
            mergeFlatZones(winnerRep, c);
        }

        return mergeAdjacentCandidatesInPlace(winnerRep, repsFlatzones);
    }

private:
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
     * @brief Avança a geração de stamp e limpa ao estourar.
     */
    void beginToken() {
        if (++token == 0u) {
            std::fill(stamp.begin(), stamp.end(), 0u);
            token = 1u;
            std::cout << "  [FlatZonesGraphOnDemandEdgesByBoundary] Stamp reset.\n";
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

        auto& pixelToIndex = pixelBuffer->pixelToIndex;
        auto& indexToPixel = pixelBuffer->indexToPixel;
        auto& sizeSets = pixelBuffer->sizeSets;
        auto& pixelsNext = pixelBuffer->pixelsNext;
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



/**
 * @brief Grafo "naive" de flatzones sem lista de adjacências.
 *
 * Esta versão mantém apenas:
 * - PixelSetManager (listas circulares de pixels);
 * - DSU (parent) para merges e representantes canônicos;
 * - deduplicação por "stamp" nas consultas.
 *
 * As adjacências são computadas sob demanda, varrendo TODOS os pixels
 * da flatzone base e testando seus vizinhos.
 *
 * Complexidade:
 * - merge: O(1) amortizado (DSU + splice O(1) de listas circulares);
 * - consulta: O(|FZ| * grau).
 */
class FlatZonesGraphOnDemandEdgesByPixel : public FlatZonesGraphBase {
public:
    FlatZonesGraphOnDemandEdgesByPixel(const FlatZonesGraphOnDemandEdgesByPixel&) = default;
    FlatZonesGraphOnDemandEdgesByPixel& operator=(const FlatZonesGraphOnDemandEdgesByPixel&) = default;

    FlatZonesGraphOnDemandEdgesByPixel(ImageUInt8Ptr imagePtr, int radiusAdj)
        : FlatZonesGraphOnDemandEdgesByPixel(imagePtr,
              std::make_shared<AdjacencyRelation>(imagePtr->getNumRows(),
                                                  imagePtr->getNumCols(),
                                                  radiusAdj)) {}

    FlatZonesGraphOnDemandEdgesByPixel(ImageUInt8Ptr imagePtr, AdjacencyRelationPtr adjRelation)
        : FlatZonesGraphBase(std::move(imagePtr), std::move(adjRelation)) {
        buildFlatZones();
    }

private:

    int mergeFlatZones(int repA, int repB) {
        int idxA = pixelView.pixelToIndex[repA];
        int idxB = pixelView.pixelToIndex[repB];
        if (idxA < 0 || idxB < 0) return findRepresentative(repA);

        int rootA = findIndex(idxA);
        int rootB = findIndex(idxB);
        if (rootA == rootB) return pixelView.indexToPixel[rootA];

        int repRootA = pixelView.indexToPixel[rootA];
        int repRootB = pixelView.indexToPixel[rootB];

        int sizeA = pixelView.sizeSets[rootA];
        int sizeB = pixelView.sizeSets[rootB];

        int winnerRoot = rootA;
        int loserRoot = rootB;
        int winnerRep = repRootA;
        int loserRep = repRootB;

        if (sizeB > sizeA) {
            winnerRoot = rootB;
            loserRoot = rootA;
            winnerRep = repRootB;
            loserRep = repRootA;
        }

        parent[loserRoot] = winnerRoot;
        pixelBuffer->mergeSetsByRep(winnerRep, loserRep);
        return winnerRep;
    }

public:

    template <class Func>
    void forEachAdjacentFlatzoneFromPixel(int rep, Func&& func) {
        int baseRep = findRepresentative(rep);
        int baseIdxRaw = pixelView.pixelToIndex[baseRep];
        if (baseIdxRaw < 0) return;
        int baseIdx = findIndex(baseIdxRaw);

        beginToken();
        for (int p : pixelBuffer->getPixelsBySet(baseRep)) {
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0) continue;
                int rootQ = (parent[idxQ] == idxQ) ? idxQ : findIndex(idxQ);
                if (rootQ == baseIdx) continue;
                if (stamp[rootQ] != token) {
                    stamp[rootQ] = token;
                    int neighRep = pixelView.indexToPixel[rootQ];
                    if (neighRep != -1) func(neighRep);
                }
            }
        }
    }

    template <class Func>
    inline void forEachAdjacentFlatzoneFromPixelStatic(int rep, Func&& func) {
        const int baseIdx = pixelView.pixelToIndex[rep];
        if (baseIdx < 0) return;

        beginToken();
        for (int p : pixelBuffer->getPixelsBySet(rep)) {
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0 || idxQ == baseIdx) continue;
                if (stamp[idxQ] != token) {
                    stamp[idxQ] = token;
                    int neighRep = pixelView.indexToPixel[idxQ];
                    if (neighRep != -1) func(neighRep);
                }
            }
        }
    }

    int mergeAdjacentCandidatesInPlace(int repFlatzone, std::vector<int>& repsFlatzones) {
        // Complexidade (naive): O(|FZ_base| * d + R * α(N)).
        // - |FZ_base| = numero de pixels da flatzone base.
        // - d = grau da adjacencia (4/8 etc.).
        // - R = repsFlatzones.size(); α(N) ~ constante (DSU).
        const int baseCanon = findRepresentative(repFlatzone);
        const int baseIdx = markAdjacentRootsFromRep(baseCanon);
        if (baseIdx < 0) return baseCanon;

        tmpCandidatesSize = 0;
        ensureTmpCandidatesCapacity((int)repsFlatzones.size() + 1);

        int winnerRep = baseCanon;
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

        if (winnerRep != baseCanon) {
            tmpCandidatesBuffer[tmpCandidatesSize++] = baseCanon;
        }

        for (int i = 0; i < tmpCandidatesSize; ++i) {
            int loser = tmpCandidatesBuffer[i];
            if (loser == winnerRep) continue;
            mergeFlatZones(winnerRep, loser);
        }

        const int wCanon = findRepresentative(winnerRep);
        repsFlatzones.erase(
            std::remove_if(repsFlatzones.begin(), repsFlatzones.end(),
                           [&](int x) { return findRepresentative(x) == wCanon; }),
            repsFlatzones.end());
        repsFlatzones.push_back(wCanon);
        return wCanon;
    }

    int mergeBasesWithAdjacentCandidatesInPlace(const std::vector<int>& baseReps,
                                               std::vector<int>& repsFlatzones,
                                               int winnerRep = -1) {
        // Complexidade (naive): O(B * α(N) + |FZ_winner| * d + R * α(N)).
        // - B = baseReps.size().
        // - |FZ_winner| pode crescer apos unir bases.
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
    std::vector<int> tmpCandidatesBuffer;
    int tmpCandidatesSize{0};

    std::vector<uint32_t> stamp;
    uint32_t token{1};
    GenerationStampSet repCNPsStamp;

private:
    void beginToken() {
        if (++token == 0u) {
            std::fill(stamp.begin(), stamp.end(), 0u);
            token = 1u;
        }
    }

    int markAdjacentRootsFromRep(int baseRep) {
        int baseIdxRaw = pixelView.pixelToIndex[baseRep];
        if (baseIdxRaw < 0) return -1;
        int baseIdx = findIndex(baseIdxRaw);

        beginToken();
        for (int p : pixelBuffer->getPixelsBySet(baseRep)) {
            for (int q : adj->getNeighborPixels(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (idxQ < 0) continue;
                int rootQ = (parent[idxQ] == idxQ) ? idxQ : findIndex(idxQ);
                if (rootQ == baseIdx) continue;
                if (stamp[rootQ] != token) {
                    stamp[rootQ] = token;
                }
            }
        }

        return baseIdx;
    }

    void buildFlatZones() {
        int numPixels = imgPtr->getSize();
        const auto* imageData = imgPtr->rawData();

        std::vector<uint8_t> visited(numPixels, 0);
        FastQueue<int> queue(numPixels / 4);

        auto& pixelToIndex = pixelBuffer->pixelToIndex;
        auto& indexToPixel = pixelBuffer->indexToPixel;
        auto& sizeSets = pixelBuffer->sizeSets;
        auto& pixelsNext = pixelBuffer->pixelsNext;
        int numFZ = 0;

        for (int p = 0; p < numPixels; ++p) {
            if (visited[p]) continue;

            int tail = p;
            int sizeFlatzone = 0;
            int idxFZ = numFZ++;
            queue.push(p);
            visited[p] = 1;
            pixelToIndex[p] = idxFZ;
            indexToPixel[idxFZ] = p;

            const uint8_t level = imageData[p];
            while (!queue.empty()) {
                int q = queue.pop();
                sizeFlatzone++;

                for (int nq : adj->getNeighborPixels(q)) {
                    if (!visited[nq] && imageData[nq] == level) {
                        visited[nq] = 1;
                        queue.push(nq);
                        pixelToIndex[nq] = idxFZ;
                        pixelsNext[tail] = nq;
                        tail = nq;
                    }
                }
            }

            pixelsNext[tail] = p;
            sizeSets[idxFZ] = sizeFlatzone;
        }

        pixelBuffer->shrinkToNumSets(numFZ);
        pixelView = pixelBuffer->view();

        parent.resize(numFZ);
        std::iota(parent.begin(), parent.end(), 0);

        stamp.assign(numFZ, 0u);
        token = 1u;

        tmpCandidatesBuffer.resize(numFZ);
        tmpCandidatesSize = 0;
    }

    void ensureTmpCandidatesCapacity(int required) {
        if (required > (int)tmpCandidatesBuffer.size()) {
            tmpCandidatesBuffer.resize(required);
        }
    }
};



/**
 * @brief Grafo de flat-zones com união eficiente e adjacências deduplicadas.
 *
 * `FlatZonesGraphFullEdges` constrói e mantém um grafo cujos vértices são flat-zones
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
 * FlatZonesGraphFullEdges fzg(img, adj);
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
class FlatZonesGraphFullEdges : public FlatZonesGraphBase {
private:
    std::vector<AdjacentFlatZones> listOfAdjacentFlatZones;
    std::vector<int> tmpCandidatesBuffer;
    int tmpCandidatesSize{0};
    

    /**
     * Método auxiliar do union-find (DSU) para encontrar o representante de um conjunto.
     * @param i Índice do elemento.
     * @return Índice do representante do conjunto ao qual o elemento pertence.
     */
    int find(int index) {
        return findIndex(index);
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

    int canonicalRep(int repPixel) {
        int idx = repIndexFromPixel_(repPixel);
        if (idx < 0) return repPixel;
        int root = find(idx);
        return pixelView.indexToPixel[root];
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
        pixelBuffer->mergeSetsByRep(repWinner, repLoser);
    }

    void ensureTmpCandidatesCapacity(int required) {
        if (required > (int)tmpCandidatesBuffer.size()) {
            tmpCandidatesBuffer.resize(required);
        }
    }

public:
    

    FlatZonesGraphFullEdges(const FlatZonesGraphFullEdges&) = default;
    FlatZonesGraphFullEdges& operator=(const FlatZonesGraphFullEdges&) = default;    
    FlatZonesGraphFullEdges(ImageUInt8Ptr imgPtr, int radiusAdj) : FlatZonesGraphFullEdges(imgPtr, std::make_shared<AdjacencyRelation>(imgPtr->getNumRows(), imgPtr->getNumCols(), radiusAdj)) {}
    
    /**
     * Construtor que cria o grafo de flatzones a partir de uma imagem e uma relação de adjacência.
     * @param imgPtr Ponteiro para a imagem de entrada.
     * @param adj Ponteiro para a relação de adjacência que define os vizinhos dos pixels.
     */
    FlatZonesGraphFullEdges(ImageUInt8Ptr imagePtr, AdjacencyRelationPtr adjRelation)
        : FlatZonesGraphBase(std::move(imagePtr), std::move(adjRelation)) {

        int numPixels = imgPtr->getSize();
        auto img = imgPtr->rawData();
        std::vector<uint8_t> visited(numPixels, 0);
        std::vector<uint8_t> isBoundary(numPixels, 0); 
        
        auto& pixelToIndex = pixelBuffer->pixelToIndex;
        auto& indexToPixel = pixelBuffer->indexToPixel;
        auto& sizeSets = pixelBuffer->sizeSets;
        auto& pixelsNext = pixelBuffer->pixelsNext;
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
                if (img[q] == img[p]) continue; // mesma flatzone: ignora
                int idxQ = pixelToIndex[q];
                
                
                if (prefilter.contains(idxQ)) continue;  // já emitido para esta FZ
                prefilter.insert(idxQ);
                
                // adiciona as arestas em uma direção
                if (idxP < idxQ) {
                    listOfAdjacentFlatZones[idxP].appendUnchecked(indexToPixel[idxQ]);
                } else {
                    listOfAdjacentFlatZones[idxQ].appendUnchecked(indexToPixel[idxP]);
                }
            }
        }
        // espelha as arestas para a outra direção
        AdjacentFlatZonesSet::mirrorAndFinalize(listOfAdjacentFlatZones, pixelView);        
        
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
     * @brief Funde um conjunto de flatzones adjacentes em torno de um representante base.
     *
     * Dado um representante de flatzone (repFlatzone) e um conjunto de representantes (repsFlatzones),
     * este método identifica quais elementos de repsFlatzones são realmente adjacentes a repFlatzone,
     * determina o representante vencedor (menor pixel entre repFlatzone e os candidatos adjacentes)
     * e realiza a fusão dessas flatzones no grafo:
     *   - Atualiza as relações de adjacência (removendo os arcos dos perdedores e ligando vizinhos ao vencedor);
     *   - Atualiza a estrutura union-find para refletir a unificação;
     *   - Integra o conteúdo (predecessores/pixels) das flatzones perdedoras no vencedor.
     *
     * Após a fusão, repsFlatzones é modificada para remover os elementos realmente fundidos,
     * garantindo a presença do vencedor.
     *
     * @param repFlatzone Representante base da flatzone a partir do qual a fusão é feita.
     * @param repsFlatzones Lista de representantes candidatos à fusão (pode ser modificada).
     * @return O representante vencedor (ID) da flatzone unificada.
     */
    int mergeAdjacentCandidatesInPlace(int repFlatzone, std::vector<int>& repsFlatzones) {
        // Complexidade (eager): O(R * lookup + sum(deg(loser))).
        // - R = repsFlatzones.size().
        // - lookup em AdjacentFlatZones e normalmente O(1) médio.
        // - sum(deg(loser)) vem do rewire das adjacências dos perdedores.
        // 1) Índice e adjacências do repFlatzone (canonical)
        int baseRep = canonicalRep(repFlatzone);
        int idxFlatzone = pixelView.pixelToIndex[baseRep];
        const AdjacentFlatZones& adjFlatzone = listOfAdjacentFlatZones[idxFlatzone];

        // 2) Filtrar apenas candidatos realmente adjacentes e escolher o vencedor (menor pixel)
        tmpCandidatesSize = 0;
        ensureTmpCandidatesCapacity((int)repsFlatzones.size() + 1);
        int winnerRep = baseRep;
        for (int rep : repsFlatzones) {
            int canon = canonicalRep(rep);
            if (adjFlatzone.find(canon) != adjFlatzone.end()) {
                tmpCandidatesBuffer[tmpCandidatesSize++] = canon;
                if (canon < winnerRep) winnerRep = canon;
            }
        }

        // Inclui a própria base se ela perder para o winner
        if (baseRep != winnerRep) {
            tmpCandidatesBuffer[tmpCandidatesSize++] = baseRep;
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

        // 4) Atualiza repsFlatzones: remove perdedores (inclui reps colapsados)
        repsFlatzones.erase(std::remove_if(repsFlatzones.begin(), repsFlatzones.end(),
                                           [&](int x){ return canonicalRep(x) == winnerRep; }),
                            repsFlatzones.end());
        repsFlatzones.push_back(winnerRep);
        
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
        // Complexidade (eager): O(B * α(N) + sum(deg(base_loser)) + custo de mergeAdjacentCandidatesInPlace).
        // - B = baseReps.size().
        // - sum(deg(base_loser)) vem do rewire ao fundir as bases.

        // 1) vencedor = menor representante (canonical)
        std::vector<int> baseCanon;
        baseCanon.reserve(baseReps.size());
        for (int rep : baseReps) {
            baseCanon.push_back(canonicalRep(rep));
        }
        std::sort(baseCanon.begin(), baseCanon.end());
        baseCanon.erase(std::unique(baseCanon.begin(), baseCanon.end()), baseCanon.end());
        if (baseCanon.empty()) return -1;

        if (winnerRep == -1) {
            winnerRep = baseCanon.front();
        } else {
            winnerRep = canonicalRep(winnerRep);
            if (winnerRep != baseCanon.front()) {
                winnerRep = baseCanon.front();
            }
        }
        int idxWinner = pixelView.pixelToIndex[winnerRep];

        // 2) funde todas as bases diretamente no winner
        for (int rep : baseCanon) {
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

};

static_assert(FlatZonesGraphCommonInterface<FlatZonesGraphOnDemandEdgesByBoundary>,
              "FlatZonesGraphOnDemandEdgesByBoundary must satisfy FlatZonesGraphCommonInterface");
static_assert(FlatZonesGraphCommonInterface<FlatZonesGraphOnDemandEdgesByPixel>,
              "FlatZonesGraphOnDemandEdgesByPixel must satisfy FlatZonesGraphCommonInterface");
static_assert(FlatZonesGraphCommonInterface<FlatZonesGraphFullEdges>,
              "FlatZonesGraphFullEdges must satisfy FlatZonesGraphCommonInterface");
