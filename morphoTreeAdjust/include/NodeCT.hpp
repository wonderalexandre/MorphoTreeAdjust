#include <list>
#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <utility>
#include <functional> 
#include <optional>
#include <stdexcept>
#include <type_traits>
#include "../include/Common.hpp"
#include "../include/ComponentTree.hpp"

#ifndef NODECT_H
#define NODECT_H


/**
 * @brief Handle (proxy leve) para nós da Component Tree com acesso O(1) via NodeArena.
 *
 * `NodeCT<CNPsType>` é um *wrapper* leve que referencia um nó da árvore de componentes
 * (identificado por `id`) e delega leituras/escritas para a `NodeArena` da
 * árvore dona (`ComponentTree<CNPsType>* tree`). Não possui propriedade do nó nem
 * alocação própria — serve como *view* handle seguro e barato de copiar.
 *
 * ## Principais capacidades
 * - **Acesso O(1)** aos atributos do nó (nível, área, representante, CNPs).
 * - **Relações estruturais**: `getParent()`, `addChild()`, `removeChild()`,
 *   `spliceChildren()` e as variantes por `Id` otimizadas (operam só na arena).
 * - **Ranges/Iteradores** para percursos:
 *   - `getChildren()` → range de filhos (por *proxy* `NodeCT`).
 *   - `getNodesOfPathToRoot()` → caminho até a raiz.
 *   - `getIteratorPostOrderTraversal()` → pós-ordem.
 *   - `getIteratorBreadthFirstTraversal()` → largura (BFS).
 *   - `getRepsOfCC()` → itera representantes (CNPs) em BFS na subárvore.
 * - **Métricas**: `computerNumDescendants()` e `computerNumFlatzoneDescendants()`.
 *
 * ## Semântica de tipo (CNPsType)
 * - `Pixels` (alias para `int`): cada nó guarda um único representante; métodos
 *   como `addRepCNPs(int)` substituem o valor.
 * - `FlatZones` (ex.: `std::vector<int>`): cada nó pode armazenar vários reps;
 *   há métodos especializados (habilitados por SFINAE) para adicionar/remover
 *   e consultar contagens de flat-zones.
 *
 * ## Invariantes
 * - Este handle é válido se e somente se `tree != nullptr` e `id >= 0`.
 * - As operações que alteram a estrutura (ex.: `addChild`, `setParent`) são
 *   encaminhadas à `ComponentTree`, que mantém a coerência dos vetores da arena.
 *
 * ## Complexidade (típica)
 * - Acesso/atribuição de campo: O(1).
 * - Inserção/remoção de relação pai↔filho (arena-based): O(1).
 * - Iteração por ranges/iteradores: O(k) no número de nós/itens visitados.
 *
 * ## Exemplo mínimo
 * @code
 * ComponentTree<Pixels> T = ... ;
 * NodeCT<Pixels> n = T.proxy(rootId);
 * if (n) {
 *     int lvl = n.getLevel();
 *     for (auto c : n.getChildren()) {
 *         // processa filhos por proxy
 *     }
 *     for (int rep : n.getRepsOfCC()) {
 *         // reps em BFS da subárvore de n
 *     }
 * }
 * @endcode
 *
 * @tparam CNPsType Tipo das CNPs armazenadas por nó (ex.: `Pixels=int`, `FlatZones=std::vector<int>`).
 */
template <typename CNPsType>
class NodeCT{
private:
    template <typename> friend class ComponentTree;

    NodeId id = -1;
    ComponentTree<CNPsType>* tree = nullptr;  // árvore dona
    
public:
    
    // Construtor padrão
    NodeCT() = default;
    NodeCT(ComponentTree<CNPsType>* owner, NodeId i) : id(i), tree(owner) {}

    //id >= 0 é o operador booleno
    explicit operator bool() const noexcept { return tree != nullptr && id >= 0; }

    //id é o operador == e !=
    bool operator==(const NodeCT& other) const noexcept { return tree == other.tree && id == other.id; }
    bool operator!=(const NodeCT& other) const noexcept { return !(*this == other); }
    
    //id é o operador int
    operator int() const noexcept { return id; }

    //acesso aos dados da arena
    inline CNPsType& getRepCNPs() noexcept { return tree->arena.repCNPs[id]; }
    inline const CNPsType& getRepCNPs() const noexcept{ return tree->arena.repCNPs[id]; }
    inline NodeArena<CNPsType>::RepsOfCCRangeById getRepCNPsOfCC() const { return tree->arena.getRepsOfCC(id); } //iterador
    inline int getRepNode() const noexcept{ return tree->arena.repNode[id]; }
    inline int getIndex() const noexcept{ return id; }
    //inline int getThreshold1() const noexcept{ return tree->arena.threshold1[id]; }
    //inline int getThreshold2() const noexcept{ return tree->arena.threshold2[id]; }
    inline int getLevel() const noexcept{ return tree->arena.threshold2[id]; }
    inline int32_t getArea() const noexcept{ return tree->arena.areaCC[id]; }
    inline int getNumChildren() const noexcept{ return tree->arena.childCount[id]; }
    inline void setArea(int32_t a) noexcept { tree->arena.areaCC[id] = a; }
    inline void setLevel(int lv) noexcept { tree->arena.threshold2[id] = lv; }
    inline void setRepCNPs(const CNPsType& reps) noexcept { tree->arena.repCNPs[id] = reps; }
    
    //propaga para versão por id 
    inline auto getPixelsOfCC() const { return tree->getPixelsOfCCById(id); } //iterador
    inline auto getCNPs() const { return tree->getCNPsById(id); } //iterador 
    inline int getNumFlatzone() const { return tree->getNumFlatzoneById(id); }
    inline int getNumCNPs() const { return tree->getNumCNPsById(id);} 
    inline bool isChild(NodeCT<CNPsType> node) const noexcept { return tree && node && hasChildById(id, node.getIndex()); }
    inline  NodeCT<CNPsType> getParent() { if (!tree || tree->arena.parentId[id] < 0) return {}; return NodeCT<CNPsType>(tree, tree->arena.parentId[id]);}
    inline  void setParent(NodeCT<CNPsType> node) { if (!tree) return;  tree->setParentById(id, node.getIndex()); }
    inline void addChild(NodeCT<CNPsType> child) { if (!tree || !child) return; tree->addChildById(id, child.getIndex()); }
    inline void removeChild(NodeCT<CNPsType> child, bool releaseNode) { if (!tree || !child) return; tree->removeChildById(id, child.getIndex(), releaseNode); }
    inline void spliceChildren(NodeCT<CNPsType> from) { if (!tree || !from || from.getIndex() == id) return; tree->spliceChildrenById(id, from.getIndex()); }
    inline void addRepCNPs(int rep) { if(!tree) return; tree->addRepCNPsById(id, rep);}
    inline bool isLeaf() const { return tree->isLeaf(id); }


    ///Métodos disponíveis SOMENTE para `FlatZones`
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzones(std::vector<int>& repsFlatZones, ComponentTreeFZPtr tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsToConnectedFlatzone(int repFlatZone, ComponentTreeFZPtr tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsToConnectedFlatzone(const std::vector<int>& repsBase, int repWinner, ComponentTreeFZPtr tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void removeFlatzone(int idFlatZone);
          
    // Conta descendentes (exclui o próprio nó)
    int computerNumDescendants() {
        if (!tree) return 0;
        int cont = 0;
        // Empilha apenas os filhos diretos; não conta o próprio nó
        FastQueue<int> st;
        for(int c: tree->arena.children(id)) st.push(c);
        while (!st.empty()) {
            int u = st.pop();
            ++cont; // conta nó u
            // empilha filhos de u
            for(int v: tree->arena.children(u)) st.push(v);
        }
        return cont;
    }

    // Soma de "flat-zones" nos descendentes (exclui o próprio nó)
    // (usa n->getNumFlatzone(), que você já especializa por CNPsType)
    int computerNumFlatzoneDescendants() {
        if (!tree) return 0;
        int acc = 0;
        FastQueue<int> st;
        for (NodeId c : tree->arena.children(id)) 
            st.push(c);
        while (!st.empty()) {
            NodeId u = st.pop();
            if constexpr (std::is_same_v<CNPsType, Pixels>) {
                // Para Pixels (CNPsType == int), cada nó contribui com 1 flat-zone
                acc += 1;
            } else {
                // Para FlatZones (CNPsType == std::vector<int>), soma a quantidade de reps no nó
                acc += static_cast<int>(tree->arena.repCNPs[u].size());
            }
            for (NodeId v : tree->arena.children(u)) 
                st.push(v);
        }
        return acc;
    }


    // Ranges existentes que devolvem filhos (por ponteiro lógico) continuam,
    // mas internamente usam `tree->proxy(id)` (que agora devolve handle)
    class ChildIdRange {
        int cur; ComponentTree<CNPsType>* T;
    public:
        struct It {
            int id; ComponentTree<CNPsType>* T;
            bool operator!=(const It& o) const { return id != o.id; }
            void operator++() { id = (id == -1) ? -1 : T->arena.nextSiblingId[id]; }
            NodeCT<CNPsType> operator*() const { return T->proxy(id); }
        };
        It begin() const { return {cur, T}; }
        It end()   const { return {-1,  T}; }
        ChildIdRange(int first, ComponentTree<CNPsType>* t) : cur(first), T(t) {}
    };
    auto getChildren() const {
        return ChildIdRange(tree ? tree->arena.firstChildId[id] : -1, tree);
    }

    //============= Iterator para iterar os nodes do caminho até o root==============//
    class InternalIteratorNodesOfPathToRoot {
    private:
        NodeCT<CNPsType> currentNode;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCT<CNPsType>;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCT<CNPsType>*;
        using reference = NodeCT<CNPsType>&;
    
        InternalIteratorNodesOfPathToRoot(NodeCT<CNPsType> obj) : currentNode(obj) {}
    
        InternalIteratorNodesOfPathToRoot& operator++() {
            if (currentNode) {
                currentNode = currentNode.getParent();  // Retorna outro shared_ptr
            }
            return *this;
        }
            
        bool operator==(const InternalIteratorNodesOfPathToRoot& other) const {
            const bool aEnd = !currentNode;
            const bool bEnd = !other.currentNode;
            if (aEnd || bEnd) return aEnd == bEnd;
            return currentNode.getIndex() == other.currentNode.getIndex();
        }
        bool operator!=(const InternalIteratorNodesOfPathToRoot& other) const { return !(*this == other); }
    
        reference operator*() {
            return currentNode;
        }
    };
    
    class IteratorNodesOfPathToRoot {
    private:
        NodeCT<CNPsType> instance;
    
    public:
        explicit IteratorNodesOfPathToRoot(NodeCT<CNPsType> obj) : instance(obj) {}
    
        InternalIteratorNodesOfPathToRoot begin() const { return InternalIteratorNodesOfPathToRoot(instance); }
        InternalIteratorNodesOfPathToRoot end() const {
            return InternalIteratorNodesOfPathToRoot(NodeCT<CNPsType>()); 
        }
    };
    
    // Chamador usa this como shared_ptr:
    IteratorNodesOfPathToRoot getNodesOfPathToRoot() {
        return IteratorNodesOfPathToRoot(NodeCT<CNPsType>(tree, id));
    }
    

/////////
// **Iterador para coletar ramos em pós-ordem**
// Classe do Iterador Pós-Ordem por Ramos
class InternalIteratorBranchPostOrderTraversal {
    private:
        std::stack<NodeCT<CNPsType>> processingStack;
        std::stack<NodeCT<CNPsType>> postOrderStack;
        std::list<std::list<NodeCT<CNPsType>>> branches;
        typename std::list<std::list<NodeCT<CNPsType>>>::iterator branchIterator;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::list<NodeCT<CNPsType>>;
        using pointer = std::list<NodeCT<CNPsType>>*;
        using reference = std::list<NodeCT<CNPsType>>&;
    
        InternalIteratorBranchPostOrderTraversal(NodeCT<CNPsType> root) {
            if (!root) return;
    
            std::stack<NodeCT<CNPsType>> tempStack;
            tempStack.push(root);
    
            while (!tempStack.empty()) {
                auto current = tempStack.top();
                tempStack.pop();
                postOrderStack.push(current);
    
                for (auto& child : current.getChildren()) {
                    tempStack.push(child);
                }
            }
    
            std::list<NodeCT<CNPsType>> currentBranch;
            while (!postOrderStack.empty()) {
                auto node = postOrderStack.top();
                postOrderStack.pop();
    
                if (!currentBranch.empty()) {
                    auto lastNode = currentBranch.back();
                    auto parent   = lastNode.getParent();
                    if (parent) {
                        int last = -1;
                        for (int c = parent.getFirstChildId(); c != -1; c = parent.tree->arena.nextSiblingId[c]) {
                            last = c;
                        }
                        const bool lastIsLastChild = (lastNode.getIndex() == last);
                        if (!lastIsLastChild) {
                            branches.push_back(currentBranch);
                            currentBranch.clear();
                        }
                    }
                }
    
                currentBranch.push_back(node);
            }
    
            if (!currentBranch.empty()) {
                branches.push_back(currentBranch);
            }
    
            branchIterator = branches.begin();
        }
    
        InternalIteratorBranchPostOrderTraversal& operator++() {
            if (branchIterator != branches.end()) {
                ++branchIterator;
            }
            return *this;
        }
    
        reference operator*() {
            return *branchIterator;
        }
    
        bool operator==(const InternalIteratorBranchPostOrderTraversal& other) const {
            return branchIterator == other.branchIterator;
        }
    
        bool operator!=(const InternalIteratorBranchPostOrderTraversal& other) const {
            return !(*this == other);
        }
    };
    
    // Classe externa
    class IteratorBranchPostOrderTraversal {
    private:
        NodeCT<CNPsType> root;
    public:
        explicit IteratorBranchPostOrderTraversal(NodeCT<CNPsType> root) : root(root) {}
    
        InternalIteratorBranchPostOrderTraversal begin() { return InternalIteratorBranchPostOrderTraversal(root); }
        InternalIteratorBranchPostOrderTraversal end() {
            return InternalIteratorBranchPostOrderTraversal(NodeCT<CNPsType>()); // default-vazio
        }
    };
    
    // Método da classe
    IteratorBranchPostOrderTraversal getIteratorBranchPostOrderTraversal() {
        return IteratorBranchPostOrderTraversal(NodeCT<CNPsType>(tree, id));
    }
    


    //============= Iterator para iterar os nodes de um percuso em pos-ordem ==============//
    class InternalIteratorPostOrderTraversal {
    private:
        struct Item { int id; bool expanded; };

        ComponentTree<CNPsType>* T_ = nullptr;
        FastStack<Item> st_;
        NodeId current_ = -1;

        inline void settle_() noexcept {
            while (!st_.empty()) {
                Item &top = st_.top();
                if (!top.expanded) {
                    top.expanded = true;
                    // empilha filhos (direita→esquerda na prática; inverta se quiser L→R)
                    for (int c = T_->arena.firstChildId[top.id]; c != -1; c = T_->arena.nextSiblingId[c]) {
                        st_.push(Item{c, false});
                    }
                } else {
                    current_ = st_.top().id;  // todos os filhos já emitidos
                    return;
                }
            }
            current_ = -1; // fim
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeCT<CNPsType>;
        using difference_type   = std::ptrdiff_t;
        using pointer           = NodeCT<CNPsType>*;     // apenas para conformidade
        using reference         = NodeCT<CNPsType>;      // retornaremos por valor (handle leve)

        InternalIteratorPostOrderTraversal(NodeCT<CNPsType> root) noexcept {
            if (root) {
                T_ = root.tree;
                st_.push({root.getIndex(), false});
                settle_();
            } else {
                current_ = -1;
            }
        }

        inline InternalIteratorPostOrderTraversal& operator++() noexcept {
            if (!st_.empty()) st_.pop();  // consome o atual
            settle_();                    // posiciona no próximo
            return *this;
        }

        // devolve o PROXY (handle leve) para o nó atual
        inline value_type operator*() const noexcept {
            return (current_ >= 0) ? NodeCT<CNPsType>(T_, current_) : NodeCT<CNPsType>();
        }

        inline bool operator==(const InternalIteratorPostOrderTraversal& other) const noexcept {
            const bool aEnd = (current_ == -1);
            const bool bEnd = (other.current_ == -1);
            if (aEnd || bEnd) return aEnd == bEnd;
            // fora do fim, basta comparar o id atual
            return current_ == other.current_;
        }
        inline bool operator!=(const InternalIteratorPostOrderTraversal& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorPostOrderTraversal {
        NodeCT<CNPsType> root_;
    public:
        explicit IteratorPostOrderTraversal(NodeCT<CNPsType> root) : root_(root) {}
        InternalIteratorPostOrderTraversal begin() const { return InternalIteratorPostOrderTraversal(root_); }
        InternalIteratorPostOrderTraversal end()   const { return InternalIteratorPostOrderTraversal(NodeCT<CNPsType>()); }
    };

    IteratorPostOrderTraversal getIteratorPostOrderTraversal() {
        return IteratorPostOrderTraversal(NodeCT<CNPsType>(tree, id));
    }



    //============= Iterator para iterar os nodes de um percuso em largura ==============//
    class InternalIteratorBreadthFirstTraversal {
    private:
        ComponentTree<CNPsType>* T_ = nullptr; // apenas leitura/encaminhamento
        FastQueue<int> q_;                     // guarda somente NodeIds

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type        = NodeCT<CNPsType>;    // devolvemos PROXY
        using difference_type   = std::ptrdiff_t;
        using pointer           = void;               // não expomos ponteiro real
        using reference         = NodeCT<CNPsType>;   // retornamos por valor (handle leve)

        // Constrói a partir de um proxy raiz. Enfileira apenas o id da raiz.
        explicit InternalIteratorBreadthFirstTraversal(NodeCT<CNPsType> root) noexcept {
            if (root) {
                T_ = root.tree;
                q_.push(root.getIndex());
            }
        }

        // Pré-incremento: consome o nó da frente e enfileira seus filhos por ID
        InternalIteratorBreadthFirstTraversal& operator++() noexcept {
            if (!q_.empty()) {
                int u = q_.front();
                q_.pop();
                for(int c: T_->arena.children(u))
                    q_.push(c);
            }
            return *this;
        }

        // Desreferencia: devolve um PROXY (NodeCT) para o nó atual (frente da fila)
        value_type operator*() const noexcept {
            return q_.empty() ? NodeCT<CNPsType>() : NodeCT<CNPsType>(T_, q_.front());
        }

        bool operator==(const InternalIteratorBreadthFirstTraversal& other) const noexcept {
            return q_.empty() == other.q_.empty();
        }
        bool operator!=(const InternalIteratorBreadthFirstTraversal& other) const noexcept {
            return !(*this == other);
        }
    };

    class IteratorBreadthFirstTraversal {
    private:
        NodeCT<CNPsType> root_;

    public:
        explicit IteratorBreadthFirstTraversal(NodeCT<CNPsType> root) : root_(root) {}

        InternalIteratorBreadthFirstTraversal begin() const { return InternalIteratorBreadthFirstTraversal(root_); }
        InternalIteratorBreadthFirstTraversal end()   const { return InternalIteratorBreadthFirstTraversal(NodeCT<CNPsType>()); }
    };

    // Método para expor o iterador na classe NodeCT
    IteratorBreadthFirstTraversal getIteratorBreadthFirstTraversal() {
        return IteratorBreadthFirstTraversal(NodeCT<CNPsType>(tree, id));
    }
    
};


#include "../include/NodeCT.tpp"

#endif


