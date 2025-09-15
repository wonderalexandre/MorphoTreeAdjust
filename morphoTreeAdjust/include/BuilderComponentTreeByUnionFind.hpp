#ifndef BUILDER_COMPONENT_TREE_BY_UNION_FIND_HPP
#define BUILDER_COMPONENT_TREE_BY_UNION_FIND_HPP

#include "../include/ComponentTree.hpp"

// Builder externo para construir ComponentTree via Union-Find.
// Mantém a mesma lógica do antigo tipo aninhado, mas agora como classe template independente.
template <typename CNPsType>
class BuilderComponentTreeByUnionFind {
private:
    ComponentTree<CNPsType>* tree;

public:
    explicit BuilderComponentTreeByUnionFind(ComponentTree<CNPsType>* t) : tree(t) {}

    // Ordenação estável dos pixels por nível de cinza
    std::vector<int> countingSort(ImageUInt8Ptr imgPtr) {
        int n = tree->getNumRowsOfImage() * tree->getNumColsOfImage();
        auto img = imgPtr->rawData();
        int maxvalue = img[0];
        for (int i = 1; i < n; i++) if (maxvalue < img[i]) maxvalue = img[i];

        std::vector<uint32_t> counter(maxvalue + 1, 0);
        std::vector<int> orderedPixels(n);

        if (tree->isMaxtree()) {
            for (int i = 0; i < n; i++) counter[img[i]]++;
            for (int i = 1; i < maxvalue; i++) counter[i] += counter[i - 1];
            counter[maxvalue] += counter[maxvalue - 1];
            for (int i = n - 1; i >= 0; --i) orderedPixels[--counter[img[i]]] = i;
        } else {
            for (int i = 0; i < n; i++) counter[maxvalue - img[i]]++;
            for (int i = 1; i < maxvalue; i++) counter[i] += counter[i - 1];
            counter[maxvalue] += counter[maxvalue - 1];
            for (int i = n - 1; i >= 0; --i) orderedPixels[--counter[maxvalue - img[i]]] = i;
        }
        return orderedPixels;
    }

    // Versão FlatZones: ordena representantes de FZ por nível de cinza
    template<typename U = CNPsType, std::enable_if_t<std::is_same_v<U, FlatZones>, int> = 0>
    std::vector<int> countingSort() {
        auto& pixelView = tree->pixelView; // acesso concedido via friend

        int numFZ = tree->flatzoneGraph->getNumFlatZones();
        auto img = tree->flatzoneGraph->getImage()->rawData();

        uint8_t maxvalue = img[pixelView.indexToPixel[0]];
        for (int i = 1; i < numFZ; ++i) {
            int pixelID = pixelView.indexToPixel[i];
            uint8_t val = img[pixelID];
            if (val > maxvalue) maxvalue = val;
        }
        std::vector<uint32_t> counter(maxvalue + 1);
        std::vector<int> orderedFlatzones(numFZ);

        if (tree->isMaxtree()) {
            for (int i = 0; i < numFZ; ++i) {
                int pixelID = pixelView.indexToPixel[i];
                uint8_t gray = img[pixelID];
                counter[gray]++;
            }
            for (int i = 1; i < maxvalue; ++i) counter[i] += counter[i - 1];
            counter[maxvalue] += counter[maxvalue - 1];
            for (int i = numFZ - 1; i >= 0; --i) {
                int pixelID = pixelView.indexToPixel[i];
                uint8_t gray = img[pixelID];
                orderedFlatzones[--counter[gray]] = pixelID;
            }
        } else {
            for (int i = 0; i < numFZ; ++i) {
                int pixelID = pixelView.indexToPixel[i];
                uint8_t gray = img[pixelID];
                counter[maxvalue - gray]++;
            }
            for (int i = 1; i < maxvalue; ++i) counter[i] += counter[i - 1];
            counter[maxvalue] += counter[maxvalue - 1];
            for (int i = numFZ - 1; i >= 0; --i) {
                int pixelID = pixelView.indexToPixel[i];
                uint8_t gray = img[pixelID];
                orderedFlatzones[--counter[maxvalue - gray]] = pixelID;
            }
        }
        return orderedFlatzones;
    }

    // FlatZones: constrói por união de flat-zones pré-computadas
    template<typename U = CNPsType, std::enable_if_t<std::is_same_v<U, FlatZones>, int> = 0>
    void createTreeByUnionFind() {
        std::vector<int> orderedPixelFlatzones = countingSort();
        auto& pixelView = tree->pixelView;
        auto* flatzoneGraph = tree->flatzoneGraph.get();
        auto& pixelToNodeId = tree->pixelToNodeId;

        int numFZ = flatzoneGraph->getNumFlatZones();
        std::vector<int> zPar(numFZ, -1);
        std::vector<int> parent(numFZ, -1);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };

        for (int i = numFZ - 1; i >= 0; i--) {
            int p = orderedPixelFlatzones[i];
            int idxP = pixelView.pixelToIndex[p];
            zPar[idxP] = idxP;
            parent[idxP] = idxP;
            for (int q : flatzoneGraph->getAdjacentFlatzonesFromPixel(p)) {
                int idxQ = pixelView.pixelToIndex[q];
                if (zPar[idxQ] != -1) {
                    int idxR = findRoot(idxQ);
                    if (idxP != idxR) {
                        parent[idxR] = idxP;
                        zPar[idxR] = idxP;
                    }
                }
            }
        }

        auto img = flatzoneGraph->getImage()->rawData();
        for (int i = 0; i < numFZ; i++) {
            int p = orderedPixelFlatzones[i];
            int idxP = pixelView.pixelToIndex[p];
            int idxPParent = parent[idxP];
            int pParent = pixelView.indexToPixel[idxPParent];

            if (idxP == idxPParent) {
                int threshold1 = tree->maxtreeTreeType ? 0 : 255;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->root = tree->makeNode(p, -1, threshold1, threshold2);
            } else if (img[p] != img[pParent]) {
                int threshold1 = tree->maxtreeTreeType ? img[pParent] + 1 : img[pParent] - 1;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->makeNode(p, pixelToNodeId[pParent], threshold1, threshold2);
            } else {
                pixelToNodeId[p] = pixelToNodeId[pParent];
            }
        }
    }

    // FlatZones (pixel-driven): fluxo alternativo com BFS por platôs
    template<typename U = CNPsType, std::enable_if_t<std::is_same_v<U, FlatZones>, int> = 0>
    void createTreeByUnionFind(ImageUInt8Ptr imgPtr) {
        std::vector<int> orderedPixels = countingSort(imgPtr);
        int numPixels = tree->getNumRowsOfImage() * tree->getNumColsOfImage();
        auto& pixelToNodeId = tree->pixelToNodeId;
        auto* adj = tree->adj.get();
        auto& arena = tree->arena;

        std::vector<int> zPar(numPixels, -1);
        std::vector<int> parent(numPixels, -1);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };
        auto img = imgPtr->rawData();

        for (int i = numPixels - 1; i >= 0; i--) {
            int p = orderedPixels[i];
            parent[p] = p;
            zPar[p] = p;
            for (int q : adj->getNeighborPixels(p)) {
                if (zPar[q] != -1) {
                    int r = findRoot(q);
                    if (p != r) { parent[r] = p; zPar[r] = p; }
                }
            }
        }

        int numNodes = 0;
        for (int i = 0; i < numPixels; i++) {
            int p = orderedPixels[i];
            int q = parent[p];
            if (img[parent[q]] == img[q]) parent[p] = parent[q];
            if (p == parent[p] || img[p] != img[parent[p]]) ++numNodes;
        }

        tree->reserveNodes(numNodes);

        for (int i = 0; i < numPixels; i++) {
            int p = orderedPixels[i];
            if (p == parent[p]) {
                int threshold1 = tree->maxtreeTreeType ? 0 : 255;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->root = tree->makeNode(p, -1, threshold1, threshold2);
            } else if (img[p] != img[parent[p]]) {
                int threshold1 = tree->maxtreeTreeType ? img[parent[p]] + 1 : img[parent[p]] - 1;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->makeNode(p, pixelToNodeId[parent[p]], threshold1, threshold2);
            } else {
                pixelToNodeId[p] = pixelToNodeId[parent[p]];
            }
        }

        tree->pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numPixels);
        tree->pixelView = tree->pixelBuffer->view();
        auto& pixelView = tree->pixelView;

        auto isVisited = [&](int p) { return zPar[p] == -2; };
        auto setVisited = [&](int p) { zPar[p] = -2; };
        FastQueue<int> queue(numPixels / 3);
        int nextIdxFZ = 0;

        for (int p0 = 0; p0 < numPixels; ++p0) {
            if (isVisited(p0)) continue;
            int repNode = (parent[p0] == p0 || img[parent[p0]] != img[p0]) ? p0 : parent[p0];
            int L = img[p0];
            int headFZ = p0;

            int idx = nextIdxFZ++;
            pixelView.indexToPixel[idx] = headFZ;
            pixelView.sizeSets[idx] = 0;
            pixelView.pixelsNext[headFZ] = headFZ;

            setVisited(p0);
            queue.push(p0);

            while (!queue.empty()) {
                int p = queue.pop();
                pixelView.pixelsNext[p] = pixelView.pixelsNext[headFZ];
                pixelView.pixelsNext[headFZ] = p;
                pixelView.pixelToIndex[p] = idx;
                pixelView.sizeSets[idx]++;
                for (int q : adj->getNeighborPixels(p)) {
                    if (!isVisited(q) && img[q] == L) {
                        int rq = (parent[q] == q || img[parent[q]] != img[q]) ? q : parent[q];
                        if (rq == repNode) { setVisited(q); queue.push(q); }
                    }
                }
            }

            arena.repCNPs[pixelToNodeId[repNode]].push_back(headFZ);
        }

        tree->pixelBuffer->shrinkToNumSets(nextIdxFZ);
        tree->pixelView = tree->pixelBuffer->view();
    }

    // Pixels: fluxo padrão por UF
    template<typename U = CNPsType, std::enable_if_t<std::is_same_v<U, Pixels>, int> = 0>
    void createTreeByUnionFind(ImageUInt8Ptr imgPtr) {
        std::vector<int> orderedPixels = countingSort(imgPtr);
        auto& pixelToNodeId = tree->pixelToNodeId;
        auto* adj = tree->adj.get();

        int numPixels = tree->getNumRowsOfImage() * tree->getNumColsOfImage();
        std::vector<int> zPar(numPixels, -1);
        std::vector<int> parent(numPixels, -1);
        auto findRoot = [&](int p) {
            while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
            return p;
        };
        auto img = imgPtr->rawData();

        for (int i = numPixels - 1; i >= 0; i--) {
            int p = orderedPixels[i];
            parent[p] = p;
            zPar[p] = p;
            for (int q : adj->getNeighborPixels(p)) {
                if (zPar[q] != -1) {
                    int r = findRoot(q);
                    if (p != r) { parent[r] = p; zPar[r] = p; }
                }
            }
        }

        int numNodes = 0;
        for (int i = 0; i < numPixels; i++) {
            int p = orderedPixels[i];
            int q = parent[p];
            if (img[parent[q]] == img[q]) parent[p] = parent[q];
            if (parent[p] == p || img[parent[p]] != img[p]) ++numNodes;
        }

        tree->reserveNodes(numNodes);
        tree->pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numNodes);
        tree->pixelView = tree->pixelBuffer->view();
        auto& pixelView = tree->pixelView;
        int indice = 0;
        for (int i = 0; i < numPixels; i++) {
            int p = orderedPixels[i];

            //contrução da árvore e arena
            if (p == parent[p]) {
                int threshold1 = tree->maxtreeTreeType ? 0 : 255;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->root = tree->makeNode(p, -1, threshold1, threshold2);
            } else if (img[p] != img[parent[p]]) {
                int threshold1 = tree->maxtreeTreeType ? img[parent[p]] + 1 : img[parent[p]] - 1;
                int threshold2 = img[p];
                pixelToNodeId[p] = tree->makeNode(p, pixelToNodeId[parent[p]], threshold1, threshold2);
            } else {
                pixelToNodeId[p] = pixelToNodeId[parent[p]];
            }

            //construção de PixelSetManager
            if (p == parent[p] || img[p] != img[parent[p]]) {
                pixelView.indexToPixel[indice] = p;
                pixelView.pixelToIndex[p] = indice;
                pixelView.sizeSets[indice] = 1;
                pixelView.pixelsNext[p] = p;
                indice++;
            } else {
                pixelView.pixelsNext[p] = pixelView.pixelsNext[parent[p]];
                pixelView.pixelsNext[parent[p]] = p;
                int idx = pixelView.pixelToIndex[parent[p]];
                pixelView.sizeSets[idx]++;
            }
        }
        
        assert((indice == numNodes) && "Erro na contagem de sets");
    }
};

#endif // BUILDER_COMPONENT_TREE_BY_UNION_FIND_HPP
