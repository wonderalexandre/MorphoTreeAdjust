#include <list>
#include <vector>
#include <array>
#include <unordered_set>
#include <utility>
#include <optional>
#include <functional>
#include <variant>
#include <span>

#include "../include/AdjacencyRelation.hpp"
#include "../include/Common.hpp"
#include "../include/FlatZonesGraph.hpp"


#ifndef COMPONENT_TREE_H
#define COMPONENT_TREE_H


/**
 * @brief Arena de nós para Component Trees com armazenamento contíguo e acesso rápido.
 *
 * A `NodeArena<CNPsType>` gerencia os atributos estruturais e dados de cada nó da árvore
 * de componentes (component tree) em vetores paralelos, proporcionando eficiência de cache
 * e operações O(1) para acesso a campos.
 *
 * ## Estrutura de dados
 * Cada vetor armazena um atributo específico dos nós:
 * - `repNode`: representante do conjunto de pixels (Union-Find).
 * - `threshold2`: nível (gray-level ou threshold máximo).
 * - `areaCC`: área acumulada do componente conexo.
 * - `repCNPs`: representantes de CNPs (pixels ou flat zones, dependendo do tipo template).
 * - Ponteiros estruturais: `parentId`, `firstChildId`, `nextSiblingId`, `prevSiblingId`, `lastChildId`.
 * - `childCount`: número de filhos diretos (cache).
 *
 * ## Operações principais
 * - `allocate(rep, thr1, thr2)`: cria um novo nó e inicializa seus campos com valores padrão.
 * - `reserve(n)`: reserva espaço para `n` nós, evitando realocações.
 * - `size()`: retorna o número de nós já alocados.
 * - `children(id)`: retorna um range leve (`ChildRange`) para iterar sobre os filhos de um nó.
 * - `getRepsOfCC(id)`: retorna um range BFS (`RepsOfCCRangeById`) que percorre todos os representantes
 *   da subárvore enraizada em `id`.
 * - **Reutilização de IDs**: `releaseNode(id)` libera o slot e `allocate(...)` reaproveita IDs disponíveis sem crescer os vetores.
 *
 * ## Iteradores
 * - `ChildRange`: itera os filhos diretos de um nó via `range-for`.
 * - `RepsOfCCIteratorById`: percorre em BFS os representantes armazenados nos nós de uma CC.
 *
 * ## Exemplo
 * @code
 * NodeArena<std::vector<int>> arena;
 * arena.reserve(1000);
 * NodeId root = arena.allocate(rep, thr1, thr2);
 * for (NodeId c : arena.children(root)) {
 *     // processa cada filho
 * }
 * for (int rep : arena.getRepsOfCC(root)) {
 *     // processa reps da subárvore
 * }
 * @endcode
 *
 * @tparam CNPsType Tipo de dados armazenado em `repCNPs` (ex.: Pixels=int, FlatZones=std::vector<int>).
 */
template <typename CNPsType>
class NodeArena {
    private:
    template<typename T> friend class ComponentTree; 
    template<typename T> friend class NodeCT; 
    
    std::vector<int>      repNode;       // representante do UF
    std::vector<int>      threshold2;    // level (max threshold)
    //std::vector<int>    threshold1;    // min threshold
    std::vector<int32_t>  areaCC;
    std::vector<CNPsType> repCNPs;       // Pixels:int; FlatZones: std::vector<int>
    
    std::vector<NodeId>   parentId;      // -1 = raiz
    std::vector<NodeId>   firstChildId;  // -1 = sem filhos
    std::vector<NodeId>   nextSiblingId; // -1 = sem próximo
    std::vector<NodeId>   prevSiblingId; // -1 = sem anterior
    std::vector<NodeId>   lastChildId;   // -1 = sem filhos
    std::vector<int>      childCount;    // cache para ter acesso a quantidade de filhos diretos
    // Lista de IDs livres para reutilização (LIFO)
    std::vector<NodeId>   freeIds;

    public:
   
    // --- gerenciamento ---
    inline NodeId allocate(int rep, int thr1, int thr2) {
        // 1) Reutiliza ID livre, se houver
        if (!freeIds.empty()) {
            NodeId id = freeIds.back();
            freeIds.pop_back();

            // Reinicializa o slot existente
            repNode[id]      = rep;              // representante UF
            threshold2[id]   = thr2;             // level
            (void)thr1; //threshold1[id] = thr1;
            areaCC[id]       = 0;
            if constexpr (std::is_same_v<CNPsType, Pixels>) {
                repCNPs[id] = 0;
            } else {
                repCNPs[id].clear();
            }

            parentId[id]     = -1;
            firstChildId[id] = -1;
            nextSiblingId[id]= -1;
            prevSiblingId[id]= -1;
            lastChildId[id]  = -1;
            childCount[id]   = 0;
            return id;
        }

        // 2) Caso contrário, aloca um novo slot no final (comportamento antigo)
        NodeId id = static_cast<NodeId>(repNode.size());
        repNode.push_back(rep);
        threshold2.push_back(thr2);
        (void)thr1; //threshold1.push_back(thr1);
        areaCC.push_back(0);
        repCNPs.emplace_back();   // default de CNPsType
        parentId.push_back(-1);
        firstChildId.push_back(-1);
        nextSiblingId.push_back(-1);
        prevSiblingId.push_back(-1);
        lastChildId.push_back(-1);
        childCount.push_back(0);
        return id;
    }

    inline void reserve(size_t n) {
        repNode.reserve(n); threshold2.reserve(n); /*threshold1.reserve(n);*/ areaCC.reserve(n);
        repCNPs.reserve(n); parentId.reserve(n); firstChildId.reserve(n); nextSiblingId.reserve(n); 
        prevSiblingId.reserve(n); lastChildId.reserve(n); childCount.reserve(n);
    }

    // Marca o ID como livre para reutilização. Pré-condição: nó já está desconectado (sem pai e sem filhos).
    inline void releaseNode(NodeId id) noexcept {
        // Zera/normaliza campos observáveis
        repNode[id]      = -1;      // marcador de slot livre (não conflita com IDs válidos >= 0)
        threshold2[id]   = 0;
        areaCC[id]       = 0;
        if constexpr (std::is_same_v<CNPsType, Pixels>) {
            repCNPs[id] = 0;
        } else {
            repCNPs[id].clear();
        }
        parentId[id]     = -1;
        firstChildId[id] = -1;
        nextSiblingId[id]= -1;
        prevSiblingId[id]= -1;
        lastChildId[id]  = -1;
        childCount[id]   = 0;

        freeIds.push_back(id);
    }

    // Consulta simples: retorna true se o slot parece livre
    inline bool isFree(NodeId id) const noexcept {
        return id >= 0 && id < static_cast<NodeId>(repNode.size()) && repNode[id] == -1;
    }

    size_t size() const { return repNode.size(); }


    // ============================================
    // ADIÇÃO: Range leve para filhos (range-for)
    // ============================================
    class ChildRange {
    public:
        // iterador 'input' mínimo para range-for
        class iterator {
        public:
            iterator(NodeId cur, const NodeArena* arena) noexcept
            : cur_(cur), arena_(arena) {}

            NodeId operator*() const noexcept { return cur_; }

            iterator& operator++() noexcept {
                cur_ = (cur_ == -1) ? -1 : arena_->nextSiblingId[cur_];
                return *this;
            }

            bool operator!=(const iterator& other) const noexcept {
                return cur_ != other.cur_;
            }

        private:
            NodeId cur_;
            const NodeArena* arena_;
        };

        ChildRange(NodeId first, const NodeArena* arena) noexcept
        : first_(first), arena_(arena) {}

        iterator begin() const noexcept { return iterator(first_, arena_); }
        iterator end()   const noexcept { return iterator(-1,    arena_); }

        // açúcares úteis:
        bool empty() const noexcept { return first_ == -1; }
        NodeId front() const noexcept { return first_; }

    private:
        NodeId first_;
        const NodeArena* arena_;
    };

    // uso: for (NodeId c : arena.children(parentId)) { ... }
    inline ChildRange children(NodeId id) const noexcept {
        return ChildRange(firstChildId[id], this);
    }



    // Itera os representantes (ints) em BFS sobre a CC (subárvore) do nó.
    class RepsOfCCIteratorById {
    private:
        const NodeArena* arena_ = nullptr;
        FastQueue<int> q_;               // BFS por IDs
        const int* curPtr_  = nullptr;    // ponteiro p/ bloco atual
        const int* curEnd_  = nullptr;
        int singleBuf_      = -1;         // buffer p/ CNPsType==Pixels

        void enqueueChildrenOf(int nid) {
            for(int c: arena_->children(nid)) 
                q_.push(c);
        }
        void loadBlockFromId(int nid) {
            if constexpr (std::is_same_v<CNPsType, Pixels>) {
                singleBuf_ = arena_->repNode[nid];
                curPtr_ = &singleBuf_;
                curEnd_ = &singleBuf_ + 1;
            } else {
                curPtr_ = arena_->repCNPs[nid].data();
                curEnd_ = arena_->repCNPs[nid].data() + arena_->repCNPs[nid].size();
            }
        }
        void advanceToNextNodeWithReps() {
            curPtr_ = curEnd_ = nullptr;
            while (!q_.empty()) {
                int nid = q_.pop();
                // enfileira filhos primeiro (BFS)
                enqueueChildrenOf(nid);
                // carrega reps do nó corrente
                loadBlockFromId(nid);
                if (curPtr_ != curEnd_) return; // achou bloco não-vazio
            }
            // fim
            curPtr_ = curEnd_ = nullptr;
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = int;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const int*;
        using reference         = const int&;

        RepsOfCCIteratorById(const NodeArena* arena, NodeId root, bool isEnd) {
            if (!isEnd && root >= 0) {
                arena_ = arena;
                // inicia fila com a RAIZ da CC (inclui reps do próprio nó)
                q_.push(root);
                advanceToNextNodeWithReps();
            }
        }

        reference operator*()  const { return *curPtr_; }
        pointer   operator->() const { return  curPtr_; }

        RepsOfCCIteratorById& operator++() {
            if (!curPtr_) return *this;         // já no fim
            ++curPtr_;
            if (curPtr_ == curEnd_) {           // acabou bloco do nó → próximo nó com reps
                advanceToNextNodeWithReps();
            }
            return *this;
        }

        RepsOfCCIteratorById operator++(int) { auto tmp = *this; ++(*this); return tmp; }

        bool operator==(const RepsOfCCIteratorById& other) const {
            const bool endA = (curPtr_ == nullptr && curEnd_ == nullptr && q_.empty());
            const bool endB = (other.curPtr_ == nullptr && other.curEnd_ == nullptr && other.q_.empty());
            if (endA || endB) return endA == endB;
            // estado interno diferente => não igual (suficiente p/ uso como input iterator)
            return false;
        }
        bool operator!=(const RepsOfCCIteratorById& other) const { return !(*this == other); }
    };

    class RepsOfCCRangeById {
    private:
        const NodeArena* arena_;
        NodeId root_;
    public:
        explicit RepsOfCCRangeById(const NodeArena* arena, NodeId root): arena_(arena), root_(root) {}

        RepsOfCCIteratorById begin() const { return RepsOfCCIteratorById(arena_, root_, false); }
        RepsOfCCIteratorById end()   const { return RepsOfCCIteratorById(arena_, root_, true ); }
    };

    // Exponha um método público para uso direto em range-for:
    inline RepsOfCCRangeById getRepsOfCC(NodeId id) const {
        return RepsOfCCRangeById(this, id);
    }

};


/**
 * @brief Estrutura de árvore de componentes (Component Tree) para imagens, com suporte a Pixels ou FlatZones.
 *
 * `ComponentTree<CNPsType>` organiza hierarquicamente regiões conexas de uma imagem
 * (componentes) em uma estrutura de árvore, permitindo análise multiescala e operações
 * de filtragem baseadas em atributos. O parâmetro de template `CNPsType` define se a
 * árvore é construída diretamente por pixels (`Pixels`) ou por flat-zones (`FlatZones`).
 *
 * ## Características principais
 * - **Construção**: via Union-Find (otimizado) usando pixels ou flat-zones.
 * - **Estrutura interna**: dados de nós armazenados em `NodeArena`, com acesso rápido O(1).
 * - **Proxy**: interface de nó exposta via `NodeCT<CNPsType>`, que encapsula acesso e
 *   relações pai/filho.
 *
 * ## Exemplo mínimo
 * @code
 * ImageUInt8Ptr img = ...;
 * auto adj = std::make_shared<AdjacencyRelation>(rows, cols, 1.5);
 * ComponentTree<Pixels> T(img, true, adj); // max-tree por pixels
 *
 * NodeCT<Pixels> root = T.getRoot();
 * for (auto c : root.getChildren()) {
 *     int area = c.getArea();
 * }
 * auto recon = T.reconstructionImage();
 * @endcode
 *
 * @tparam CNPsType Define o tipo de construção da árvore: `Pixels` ou `FlatZones`.
 */
template <typename CNPsType>
class ComponentTree : public std::enable_shared_from_this<ComponentTree<CNPsType>> {
protected:
    template<typename T> friend class NodeCT; 
    NodeId root;
    int numRows;
    int numCols;
    bool maxtreeTreeType; //maxtree is true; mintree is false
    AdjacencyRelationPtr adj; //disk of a given ratio: ratio(1) for 4-connect and ratio(1.5) for 8-connect 
    int numNodes;
    
    std::vector<NodeId> pixelToNodeId; //Mapeamento dos pixels representantes para NodeID. Para adquirir todos os representantes valido utilize o método getRepCNPs
    std::shared_ptr<PixelSetManager> pixelBuffer; PixelSetManager::View pixelView; //gerenciamento de pixels da arvore
    NodeArena<CNPsType> arena; // armazenamento indexado dos dados de todos os nós 
    
    NodeId makeNode(int repNode, NodeId parentId, int threshold1, int threshold2);
    void reserveNodes(int expected) { arena.reserve(expected);}
    
    // Define `flatzoneGraph` apenas para `FlatZones`
    using FlatzoneGraphType = std::conditional_t<std::is_same_v<CNPsType, FlatZones>,
        std::shared_ptr<FlatZonesGraph>, 
        std::monostate>;
    FlatzoneGraphType flatzoneGraph;

    std::vector<int> countingSort(ImageUInt8Ptr img);
    
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    std::vector<int> countingSort();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void createTreeByUnionFind(std::vector<int>& orderedPixels, ImageUInt8Ptr img);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void createTreeByUnionFind(std::vector<int>& orderedFlatzones);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, Pixels>::value, int> = 0>
    void createTreeByUnionFind(std::vector<int>& orderedPixels, ImageUInt8Ptr img);

    void reconstruction(NodeId node, uint8_t* data);
    void build(ImageUInt8Ptr img);
    void build();

public:


    ComponentTree(ImageUInt8Ptr img, bool isMaxtree, AdjacencyRelationPtr adj);
    
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    ComponentTree(std::shared_ptr<FlatZonesGraph> graph, bool isMaxtree);

    virtual ~ComponentTree() = default;

    NodeCT<CNPsType> proxy(NodeId id) const;

    NodeCT<CNPsType> createNode(int repNode,NodeCT<CNPsType> parent,int threshold1,int threshold2);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    std::shared_ptr<FlatZonesGraph>& getFlatZonesGraph();

    bool isMaxtree()const noexcept{ return maxtreeTreeType;}
    int getNumNodes()const noexcept{ return numNodes; }
    int getNumRowsOfImage()const noexcept{ return numRows;}
    int getNumColsOfImage()const noexcept{ return numCols;}
    AdjacencyRelationPtr getAdjacencyRelation() noexcept {return adj;}
    inline NodeArena<CNPsType>::RepsOfCCRangeById getRepCNPs() const { return arena.getRepsOfCC(root); } //iterador

    inline NodeArena<CNPsType>::RepsOfCCRangeById getRepCNPsOfCCById(NodeId id) const { return arena.getRepsOfCC(id); } //iterador
    inline auto getPixelsOfCCById(NodeId id) const{ return pixelBuffer->getPixelsBySet(arena.getRepsOfCC(id));  } //iterador
    inline auto getPixelsOfFlatzone(int repFZ) const{ return pixelBuffer->getPixelsBySet(repFZ); } //iterador
    inline int getLevelById(NodeId id) const noexcept{return arena.threshold2[id];}
    inline void setLevelById(NodeId id, int level){ arena.threshold2[id] = level;}
    inline CNPsType& getRepCNPsById(NodeId id) noexcept {return arena.repCNPs[id];}
    inline const CNPsType& getRepCNPsById(NodeId id) const{ return arena.repCNPs[id]; }
    inline void addRepCNPsById(NodeId id, int rep) { 
        pixelToNodeId[rep] = id; 
        if constexpr (std::is_same_v<CNPsType, int>) 
            arena.repCNPs[id] = rep;  
        else  
            arena.repCNPs[id].push_back(rep); 
    }

    void releaseNode(NodeId id) noexcept { arena.releaseNode(id); numNodes--; }
    inline int32_t getAreaById(NodeId id) const noexcept{ return arena.areaCC[id];}
    inline void setAreaById(NodeId id, int32_t area) noexcept { arena.areaCC[id] = area;}

    inline NodeId getParentById(NodeId id) const noexcept {return arena.parentId[id];}
    inline int getNumChildrenById(NodeId id) const noexcept{ return arena.childCount[id]; }
    inline bool hasChildById(NodeId nodeId, NodeId childId) const { return arena.parentId[childId] == nodeId;}
    inline bool isLeaf(NodeId id) const { return arena.childCount[id] == 0; }
    
    inline auto getChildrenById(NodeId id) const {return arena.children(id);}    
    inline NodeId getRootById()const noexcept { return this->root;}
    inline void setRootById(NodeId n){ setParentById(n, -1); this->root = n; }
    inline NodeId getSCById(int p) const noexcept{ return this->pixelToNodeId[p];}
    inline void setSCById(int p, NodeId id) noexcept { this->pixelToNodeId[p] = id;}
    
    NodeCT<CNPsType> getSC(int p) const noexcept { return proxy(this->pixelToNodeId[p]);}
    void setSC(int p, NodeCT<CNPsType> node){ setSCById(p, node); }
    NodeCT<CNPsType> getRoot() { return proxy(this->root); }
    void setRoot(NodeCT<CNPsType> n){ setRootById(n); }

    inline int getNumFlatzoneById(NodeId id) const;
    inline auto getCNPsById(NodeId id) const;
    inline int getNumCNPsById(NodeId id) const;

    inline void removeChildById(int parentId, int childId, bool release);
    inline void addChildById(int parentId, int childId);
    inline void spliceChildrenById(int toId, int fromId);
    inline void setParentById(NodeId nodeId, NodeId parentId);

    void computerArea(NodeId node);
    ImageUInt8Ptr reconstructionImage();
    std::vector<NodeId> getLeaves();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, Pixels>::value, int> = 0>
    void prunning(NodeId nodeId);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void prunning(NodeId nodeId);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, Pixels>::value, int> = 0>
    void mergeWithParent(NodeCT<CNPsType> node);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, Pixels>::value, int> = 0>
    void mergeWithParent(std::vector<int>& flatzone);
    
    // Soma de "flat-zones" nos descendentes (exclui o próprio nó)
    // (usa n->getNumFlatzone(), que você já especializa por CNPsType)
    int computerNumFlatzoneDescendants(NodeId id) {
        int acc = 0;
        FastQueue<int> st;
        for (NodeId c : this->arena.children(id)) 
            st.push(c);
        while (!st.empty()) {
            NodeId u = st.pop();
            if constexpr (std::is_same_v<CNPsType, Pixels>) {
                // Para Pixels (CNPsType == int), cada nó contribui com 1 flat-zone
                acc += 1;
            } else {
                // Para FlatZones (CNPsType == std::vector<int>), soma a quantidade de reps no nó
                acc += static_cast<int>(this->arena.repCNPs[u].size());
            }
            for (NodeId v : this->arena.children(u)) 
                st.push(v);
        }
        return acc;
    }

    // Conta descendentes (exclui o próprio nó)
    int computerNumDescendants(NodeId id) {
        int cont = 0;
        // Empilha apenas os filhos diretos; não conta o próprio nó
        FastQueue<int> st;
        for(int c: this->arena.children(id)) st.push(c);
        while (!st.empty()) {
            int u = st.pop();
            ++cont; // conta nó u
            // empilha filhos de u
            for(int v: this->arena.children(u)) st.push(v);
        }
        return cont;
    }


    static bool validateStructure(ComponentTreePtr<CNPsType> tree)  {
        return validateStructure(tree.get());
    }
    static bool validateStructure(const ComponentTree* tree){
        const auto erroMsg = [&](std::string msg){ 
            std::cerr << "❌ Erro " << msg << "\n";
            return false;
        };
        const auto infoMsg = [&](std::string msg){ 
            std::cout << "✅ " << msg << "\n";
        };

        if (tree->root < 0 || tree->root >= (int)tree->arena.size())
            return erroMsg("root inválido");

        // 1) Exigir exatamente 1 raiz (DESCONSIDERANDO slots liberados pela free-list)
        int roots = 0;
        for (int id = 0; id < (int)tree->arena.size(); ++id) {
            if (!tree->arena.isFree(id) && tree->arena.parentId[id] == -1) {
                ++roots;
            }
        }
        if (roots != 1)
            return erroMsg("1: A árvore NÃO possui exatamente uma raiz; a soma de parentId == -1 (desconsiderando slots livres) é "+ std::to_string(roots) +" mas deveria ser 1");
        else
            infoMsg("A árvore contém exatamente 1 raiz (excluindo slots livres)");
    

        // 2) Pai consistente
        for (int id = 0; id < (int)tree->arena.size(); ++id) {
            if (tree->arena.parentId[id] != -1) {
                if (tree->arena.parentId[id] < 0 || tree->arena.parentId[id] >= (int)tree->arena.size()) 
                    return erroMsg("2: O parentId="+ std::to_string(id)+" está fora do range [0, "+ std::to_string(tree->arena.size()) + "]");
                if (id == tree->arena.parentId[id]) 
                    return erroMsg("3: O parentId="+std::to_string(id)+" está apontando para si mesmo");
            }
        }
        infoMsg("A estrutura de parentesco arena.parentId está consistente");

        // 3) Encadeamento filhos/irmãos + lastChildId + childCount
        for (int u = 0; u < (int)tree->arena.size(); ++u) {
            int cnt = 0, last = -1;
            if (tree->arena.firstChildId[u] == -1) {
                if (tree->arena.lastChildId[u] != -1 || tree->arena.childCount[u] != 0) 
                    return erroMsg("4: Nó sem filhos mas last/childCount incoerentes");
            } else {
                // Caminha pela lista de filhos de u
                for (int c = tree->arena.firstChildId[u]; c != -1; c = tree->arena.nextSiblingId[c]) {
                    if (c < 0 || c >= (int)tree->arena.size()) 
                        return erroMsg("5: Estrutura de filhos e irmãos (firstChildId/nextSiblingId) fora do range [0, "+ std::to_string(tree->arena.size()) + "]");
                    if (tree->arena.parentId[c] != u) 
                        return erroMsg("6: Filho com parentId diferente do pai");
                    
                        // Checa simetria prev/next
                    if (tree->arena.prevSiblingId[c] != last) 
                        return erroMsg("7: Estrutura de irmãos prevSiblingId está inconsistente");
                    if (last != -1 && tree->arena.nextSiblingId[last] != c) 
                        return erroMsg("8: Estrutura de irmãos nextSiblingId está inconsistente");
                    last = c; 
                    ++cnt;
                }
                if (last != tree->arena.lastChildId[u]) 
                    return erroMsg("8: Estrutura de filhos lastChildId não bate com último encadeado");
                if (cnt != tree->arena.childCount[u]) 
                    return erroMsg("9: Estrutura que armazena a quantidade de filhos childCount não bate com encadeamento");
            }
        }
        infoMsg("A estrutura de filhos/irmãos está consistente");
        return true;
    }
    
    static std::vector<NodeId> getNodesThreshold(ComponentTreePtr<CNPsType> tree, int areaThreshold, bool ENABLE_LOG = false){
        std::vector<NodeId> lista;
        FastQueue<NodeId> queue;
        queue.push(tree->root);

        int sumArea = 0; // estatística
        int numFlatZones = 0;
        int numNodes = 0;
        while (!queue.empty()) {
            NodeId id = queue.pop();
            if (tree->arena.areaCC[id] > areaThreshold) {
                for(NodeId c: tree->arena.children(id)){
                    queue.push(c);
                }
            } else {
                if (ENABLE_LOG) {
                    sumArea += tree->arena.areaCC[id];
                    // Para estatísticas abaixo, usamos os helpers já existentes no NodeCT (cria handle só aqui):
                    auto h = tree->proxy(id);
                    numFlatZones += (h.computerNumFlatzoneDescendants() + h.getNumFlatzone());
                    numNodes     += h.computerNumDescendants() + 1;
                }
                lista.push_back(id);
            }
        }
        if (ENABLE_LOG) {
            int areaImage = tree->getNumColsOfImage() * tree->getNumRowsOfImage();
            std::cout << "\tArea threshold: " << areaThreshold
                    << ", #Nodes: " << numNodes
                    << ", #FlatZones: " << numFlatZones
                    << ", #InputTreeNodes: " << tree->getNumNodes()
                    << ", |Pruning Area|: " << sumArea
                    << " (" << std::fixed << std::setprecision(2)
                    << (static_cast<double>(sumArea) / areaImage) * 100.0 << "% of the image area)"
                    << std::endl;
        }
        return lista;
    }
    
    
    // ====================== Iteradores por ID (sem proxy) ====================== //


    /**
     * @brief Iterador interno que salta slots vazios e retorna IDs válidos.
     */
    class InternalIteratorValidNodeIds {
    private:
        const int* rep_;        // ponteiro p/ arena.repNode[0]
        NodeId cur_;            // posição atual
        NodeId end_;            // N = arena.repNode.size()

        // avança cur_ até um id válido ou coloca cur_ = end_
        inline void settle_() noexcept {
            while (cur_ < end_ && rep_[cur_] == InvalidNode) ++cur_;
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        inline InternalIteratorValidNodeIds(ComponentTree<CNPsType>* T, NodeId start) noexcept
        : rep_(T ? T->arena.repNode.data() : nullptr),
        cur_(T ? start : 0),
        end_(T ? static_cast<NodeId>(T->arena.repNode.size()) : 0) {
            settle_();
        }

        inline InternalIteratorValidNodeIds& operator++() noexcept {
            ++cur_;
            settle_();
            return *this;
        }

        inline NodeId operator*() const noexcept { return cur_; }

        // iterador input: comparar só posição é suficiente
        inline bool operator==(const InternalIteratorValidNodeIds& other) const noexcept { return cur_ == other.cur_; }
        inline bool operator!=(const InternalIteratorValidNodeIds& other) const noexcept { return cur_ != other.cur_; }
    };

    /**
     * @brief Range wrapper para percorrer todos os NodeId ativos da arena.
     */
    class IteratorValidNodeIds {
    private:
        ComponentTree<CNPsType>* T_ = nullptr;
    public:
        inline explicit IteratorValidNodeIds(ComponentTree<CNPsType>* T) noexcept : T_(T) {}

        inline InternalIteratorValidNodeIds begin() const noexcept { return InternalIteratorValidNodeIds(T_, 0); }
        inline InternalIteratorValidNodeIds end() const noexcept {
            // end iterator shares the same end_ (size) as begin();
            // if T_ is null, both begin/end will be empty
            return InternalIteratorValidNodeIds(T_, T_ ? static_cast<NodeId>(T_->arena.repNode.size()) : 0);
        }
    };

    /** Range para iterar NodeId válidos (exclui slots livres) */
    inline IteratorValidNodeIds getNodeIds() noexcept { return IteratorValidNodeIds(this); }
    inline IteratorValidNodeIds getNodeIds() const noexcept { return IteratorValidNodeIds(const_cast<ComponentTree<CNPsType>*>(this)); }


    // Pós-ordem (retorna NodeId)
    class InternalIteratorPostOrderTraversalId {
    private:
        struct Item { int id; bool expanded; };

        ComponentTree<CNPsType>* T_ = nullptr;
        FastStack<Item> st_;
        NodeId current_ = -1;

        // Avança até o próximo nó a ser emitido (ou deixa current_ = -1 se acabou)
        void settle_() noexcept {
            while (!st_.empty()) {
                Item &top = st_.top();
                if (!top.expanded) {
                    top.expanded = true;

                    // A ordem resultante não é garantida (costuma ser direita->esquerda).
                    for (int c = T_->arena.firstChildId[top.id]; c != -1; c = T_->arena.nextSiblingId[c]) {
                        st_.push(Item{c, false});
                    }
                    // volta ao loop: agora o topo será algum filho
                } else {
                    current_ = top.id;      // todos os filhos já emitidos
                    return;
                }
            }
            current_ = -1; // fim
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        InternalIteratorPostOrderTraversalId(ComponentTree<CNPsType>* T, NodeId rootId) noexcept : T_(T) {
            if (T_ && rootId >= 0) {
                st_.push(Item{rootId, false});
                settle_(); // posiciona no primeiro elemento
            } else {
                current_ = -1;
            }
        }

        // pré-incremento
        InternalIteratorPostOrderTraversalId& operator++() noexcept {
            if (!st_.empty()) st_.pop();  // consome o atual
            settle_();                    // acha o próximo
            return *this;
        }

        // desreferência → NodeId atual
        NodeId operator*() const noexcept { return current_; }

        bool operator==(const InternalIteratorPostOrderTraversalId& other) const noexcept {
            return current_ == other.current_;
        }
        bool operator!=(const InternalIteratorPostOrderTraversalId& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorPostOrderTraversalId {
    private:
        ComponentTree<CNPsType>* T_ = nullptr;
        int rootId_ = -1;
    public:
        explicit IteratorPostOrderTraversalId(ComponentTree<CNPsType>* T, int rootId) noexcept : T_(T), rootId_(rootId) {}

        InternalIteratorPostOrderTraversalId begin() const noexcept {
            return InternalIteratorPostOrderTraversalId(T_, rootId_);
        }
        InternalIteratorPostOrderTraversalId end() const noexcept {
            return InternalIteratorPostOrderTraversalId(nullptr, -1);
        }
    };

    auto getIteratorPostOrderTraversalById(int id) {        
        return IteratorPostOrderTraversalId(this, id);
    }
    auto getIteratorPostOrderTraversalById() {        
        return IteratorPostOrderTraversalId(this, root);
    }

    // Largura (BFS) — retorna NodeId
    class InternalIteratorBreadthFirstTraversalId {
    private:
        ComponentTree<CNPsType>* T_ = nullptr;
        FastQueue<int> q_;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = int;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const int*;
        using reference         = const int&;

        InternalIteratorBreadthFirstTraversalId(ComponentTree<CNPsType>* T, int rootId) noexcept : T_(T) {
            if (T_ && rootId >= 0) q_.push(rootId);
        }

        InternalIteratorBreadthFirstTraversalId& operator++() noexcept {
            if (!q_.empty()) {
                int u = q_.pop();
                for (int c: T_->arena.children(u)) {
                    q_.push(c);
                }
            }
            return *this;
        }

        int operator*() const noexcept { return q_.front(); }

        bool operator==(const InternalIteratorBreadthFirstTraversalId& other) const noexcept {
            return q_.empty() == other.q_.empty();
        }
        bool operator!=(const InternalIteratorBreadthFirstTraversalId& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorBreadthFirstTraversalId {
    private:
        ComponentTree<CNPsType>* T_ = nullptr;
        int rootId_ = -1;
    public:
        explicit IteratorBreadthFirstTraversalId(ComponentTree<CNPsType>* T, int rootId) noexcept : T_(T), rootId_(rootId) {}

        InternalIteratorBreadthFirstTraversalId begin() const noexcept {
            return InternalIteratorBreadthFirstTraversalId(T_, rootId_);
        }
        InternalIteratorBreadthFirstTraversalId end() const noexcept {
            return InternalIteratorBreadthFirstTraversalId(nullptr, -1);
        }
    };

    IteratorBreadthFirstTraversalId getIteratorBreadthFirstTraversalById(NodeId id) noexcept {
        return IteratorBreadthFirstTraversalId(this, id);
    }
    IteratorBreadthFirstTraversalId getIteratorBreadthFirstTraversalById() noexcept {
        return IteratorBreadthFirstTraversalId(this, root);
    }
    // ================== Fim dos iteradores por ID (sem proxy) ================== //

    // ================== Iterador para caminho até a raiz por NodeId (sem proxy) ================== //
    class InternalIteratorNodesOfPathToRootId {
    private:
        ComponentTree<CNPsType>* T_ = nullptr;
        NodeId currentId_ = -1;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeId;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const NodeId*;
        using reference         = const NodeId&;

        InternalIteratorNodesOfPathToRootId(ComponentTree<CNPsType>* T, NodeId startId) noexcept : T_(T), currentId_(startId) {}

        InternalIteratorNodesOfPathToRootId& operator++() noexcept {
            if (T_ && currentId_ != -1) {
                currentId_ = T_->arena.parentId[currentId_];
            }
            return *this;
        }

        NodeId operator*() const noexcept { return currentId_; }

        bool operator==(const InternalIteratorNodesOfPathToRootId& other) const noexcept {
            return currentId_ == other.currentId_;
        }
        bool operator!=(const InternalIteratorNodesOfPathToRootId& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorNodesOfPathToRootId {
    private:
        ComponentTree<CNPsType>* T_ = nullptr;
        NodeId startId_ = -1;

    public:
        IteratorNodesOfPathToRootId(ComponentTree<CNPsType>* T, NodeId startId) noexcept : T_(T), startId_(startId) {}

        InternalIteratorNodesOfPathToRootId begin() const noexcept {
            return InternalIteratorNodesOfPathToRootId(T_, startId_);
        }
        InternalIteratorNodesOfPathToRootId end() const noexcept {
            return InternalIteratorNodesOfPathToRootId(nullptr, -1);
        }
    };

    IteratorNodesOfPathToRootId getNodesOfPathToRootById(NodeId id) noexcept {
        return IteratorNodesOfPathToRootId(this, id);
    }
};


#include "../include/ComponentTree.tpp"

#endif
