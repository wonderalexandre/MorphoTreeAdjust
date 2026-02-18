#pragma once


#define NDEBUG  // Remove os asserts do código
#include <array>
#include <chrono>
#include <cstdint>
#include <algorithm>
#include <iterator>
#include <memory>
#include <span>
#include <utility>
#include <vector>


#define PRINT_LOG 0
#define PRINT_DEBUG 0

// Forward declaration dos templates
class FlatZonesGraphOnDemandEdgesByBoundary;
class FlatZonesGraphOnDemandEdgesByPixel;
class FlatZonesGraphFullEdges;

using DefaultFlatZonesGraph = FlatZonesGraphOnDemandEdgesByBoundary;

template <typename GraphT = DefaultFlatZonesGraph>
using FlatZonesGraphPtr = std::shared_ptr<GraphT>;

template <typename T, typename GraphT = DefaultFlatZonesGraph>
class ComponentTree;

template <typename T, typename GraphT = DefaultFlatZonesGraph>
class NodeCT;


using NodeId = int; 
constexpr NodeId InvalidNode = -1; //-1 indica nó inválido
inline bool isValidNode(NodeId id) noexcept { return id != InvalidNode;}
inline bool isInvalid(NodeId id) noexcept { return id == InvalidNode; }

// Definição de tipos para CNPs
using Pixels = int;                 // representante dos CNPs do nó. Esse representante serve para acessar os CNPs no PixelSetManager. É o mesmo que repNode.
using FlatZones = std::vector<int>; // lista de representantes das flatzones. Esses representantes servem para acessar os pixels da FZ no PixelSetManager.


template <typename GraphT = DefaultFlatZonesGraph>
using ComponentTreeFZ = ComponentTree<FlatZones, GraphT>; // representa uma component tree por flatzones
using ComponentTreeP = ComponentTree<Pixels, DefaultFlatZonesGraph>; // representa uma component tree sem tratamento de flatzones

template <typename GraphT = DefaultFlatZonesGraph>
using NodeFZ = NodeCT<FlatZones, GraphT>;
using NodeP  = NodeCT<Pixels, DefaultFlatZonesGraph>;

template <typename CNPsType, typename GraphT = DefaultFlatZonesGraph>
using ComponentTreePtr = std::shared_ptr<ComponentTree<CNPsType, GraphT>>; 
template <typename GraphT = DefaultFlatZonesGraph>
using ComponentTreeFZPtr = std::shared_ptr<ComponentTreeFZ<GraphT>>; // representa uma component tree por flatzones
using ComponentTreePPtr = std::shared_ptr<ComponentTreeP>; // representa uma component tree sem tratamento de flatzones


class ImageUtils{
public:
    // Converte (row, col) para índice 1D (row-major)
    static int to1D(int row, int col, int numCols) {
        return row * numCols + col;
    }

    // Converte índice 1D para (row, col) (row-major)
    static std::pair<int, int> to2D(int index, int numCols) {
        int row = index / numCols;
        int col = index % numCols;
        return std::make_pair(row, col);
    }
};

/**
 * @brief Classe de imagem genérica 2D com armazenamento contíguo e controle de vida via std::shared_ptr.
 *
 * A `Image<PixelType>` representa uma imagem 2D em ordem row-major, encapsulando
 * largura, altura e um buffer contíguo gerenciado por `std::shared_ptr<PixelType[]>`.
 * Fornece utilitários para criação/copiar/preencher e acesso indexado em 1D.
 *
 * ## Semântica de propriedade do buffer
 * - `create(rows, cols)`: aloca um novo buffer e o gerencia (deleter padrão).
 * - `create(rows, cols, initValue)`: idem, preenchendo com o valor inicial.
 * - `fromExternal(rawPtr, rows, cols)`: **não** assume a propriedade
 *   (o deleter é vazio). Útil quando o ciclo de vida do ponteiro bruto é externo.
 * - `fromRaw(rawPtr, rows, cols)`: **assume** a propriedade do array
 *   (deleter padrão de array). Use quando a instância deve gerenciar a memória.
 *
 * ## Layout de memória
 * - Acesso linear (row-major): índice `i = row * numCols + col`.
 * - Operador `operator[](int)` fornece acesso por índice linear.
 *
 *
 * ## Exemplo de uso
 * @code
 * using ImgU8 = Image<uint8_t>;
 * auto img = ImgU8::create(480, 640, 0);     // aloca e zera
 * img->fill(255);                            // preenche com 255
 * int idx = ImageUtils::to1D(10, 20, img->getNumCols());
 * (*img)[idx] = 128;                         // acesso linear
 * auto clone = img->clone();                  // deep copy
 * bool eq = img->isEqual(clone);              // true
 * @endcode
 *
 * @tparam PixelType tipo do pixel armazenado (ex.: uint8_t, int32_t, float).
 */
template <typename PixelType>
class Image {
    private:
        int numRows;
        int numCols;
        std::shared_ptr<PixelType[]> data;
        using Ptr = std::shared_ptr<Image<PixelType>>;

    public:
    
    Image(int rows, int cols): numRows(rows), numCols(cols), data(new PixelType[rows * cols], std::default_delete<PixelType[]>()) {}

    static Ptr create(int rows, int cols) {
        return std::make_shared<Image>(rows, cols);
    }

    static Ptr create(int rows, int cols, PixelType initValue) {
        auto img = create(rows, cols);
        img->fill(initValue);
        return img;
    }

    static Ptr fromExternal(PixelType* rawPtr, int rows, int cols) {
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, [](PixelType*) {
            // deleter vazio: não libera o ponteiro
        });
        return img;
    }

    static Ptr fromRaw(PixelType* rawPtr, int rows, int cols) {
        auto img = create(rows, cols);
        img->data = std::shared_ptr<PixelType[]>(rawPtr, std::default_delete<PixelType[]>());
        return img;
    }

    
    void fill(PixelType value) {
        std::fill_n(data.get(), numRows * numCols, value);
    }

    bool isEqual(const Ptr& other) const {
        if (numRows != other->numRows || numCols != other->numCols)
            return false;
        int n = numRows * numCols;
        for (int i = 0; i < n; ++i) {
            if (data[i] != (*other)[i])
                return false;
        }
        return true;
    }
    
    Ptr clone() const {
        auto newImg = create(numRows, numCols);
        std::copy(data.get(), data.get() + (numRows * numCols), newImg->data.get());
        return newImg;
    }

    std::shared_ptr<PixelType[]> rawDataPtr(){ return data; }
    PixelType* rawData() { return data.get(); }
    int getNumRows() const { return numRows; }
    int getNumCols() const { return numCols; }
    int getSize() const { return numRows * numCols; }
    PixelType& operator[](int index) { return data[index]; }
    const PixelType& operator[](int index) const { return data[index]; }


};

// Aliases
using ImageUInt8 = Image<uint8_t>;
using ImageInt32 = Image<int32_t>;
using ImageFloat = Image<float>;

using ImageUInt8Ptr = std::shared_ptr<ImageUInt8>;
using ImageInt32Ptr = std::shared_ptr<ImageInt32>;
using ImageFloatPtr = std::shared_ptr<ImageFloat>;

template <typename T>
using ImagePtr = std::shared_ptr<Image<T>>;





/**
 * @brief Estrutura de dados para marcação eficiente (visited flags) usando carimbos de geração.
 *
 * A `GenerationStampSet` mantém um array de inteiros (stamps), cada posição
 * associada a um índice de elemento (ex.: nó de grafo). Em vez de limpar o
 * array inteiro a cada iteração, um contador de geração (`cur`) é incrementado
 * e usado como "marca lógica". 
 *
 *
 * @code
 * GenerationStampSet visited(numNodes);
 *
 * visited.mark(nodeIdx);
 *
 * if (!visited.isMarked(otherIdx)) {
 *     // processa nó não visitado
 * }
 *
 * visited.resetAll();  // O(1) para preparar nova iteração
 * @endcode
 */
struct GenerationStampSet {
    using gen_t = uint32_t;

    std::vector<gen_t> stamp; // array de carimbos
    size_t n{0};              // tamanho
    gen_t cur{1};             // geração atual (0 = “limpo”)

    GenerationStampSet() = default;
    explicit GenerationStampSet(size_t n) { resize(n); }

    void resize(size_t newN) {
        n = newN;
        stamp.assign(n, 0);
        cur = 1;
    }

    inline void mark(size_t idx) noexcept {
        stamp[idx] = cur;
    }

    inline bool isMarked(size_t idx) const noexcept {
        return stamp[idx] == cur;
    }

    // reset lógico em O(1)
    void resetAll() {
        if (++cur == 0) {
            std::fill(stamp.begin(), stamp.end(), 0);
            cur = 1;
        }
    }

    // limpeza forçada em O(N)
    void clearAll() {
        std::fill(stamp.begin(), stamp.end(), 0);
        cur = 1;
    }

    gen_t generation() const noexcept { return cur; }
};


/**
 * @brief Gerenciador de conjuntos disjuntos de pixels (flat zones ou CNPs) com listas circulares e mapeamentos O(1).
 *
 * O `PixelSetManager` mantém a relação entre pixels e seus conjuntos (flat zones ou CNPs)
 * usando quatro vetores paralelos: `pixelToIndex`, `indexToPixel`, `sizeSets` e
 * `pixelsNext`. O desenho provê operações O(1) para consultas e splices
 * (concatenação de listas circulares) durante fusões de conjuntos, além de
 * *views* baseadas em `std::span` para iteração sem cópias.
 *
 * ## Estruturas internas
 * - `pixelToIndex[p]` → índice do conjunto ao qual o pixel representante `p` pertence.
 * - `indexToPixel[i]` → representante (pixel head) do conjunto de índice `i`.
 * - `sizeSets[i]` → tamanho (número de pixels) do conjunto `i`.
 * - `pixelsNext[p]` → próximo pixel na lista circular do conjunto ao qual `p` pertence.
 *
 * ## Operações principais
 * - `numSets()` (conjuntos ativos), `numPixelsInSet(rep)`, `numPixelsInSets(reps)` — consultas O(1)/O(k).
 * - `mergeSetsByRep(repWinner, repLoser)` — fusão O(1) com splice de listas circulares.
 * - `shrinkToNumSets(n)` — reduz vetores de conjuntos ao número real de FZs ou |numNodes|
 * - *Views*: `view()`, `viewOf*()` expõem `std::span` para zero-cópia.
 * - Iteração de pixels por set: `getPixelsBySet(...)` retorna um range lazy.
 * - Iteração de representantes válidos: `getFlatzoneRepresentatives()` retorna um range.
 *
 * ## Complexidade
 * - Acesso por mapeamento: O(1).
 * - Fusão (`mergeSetsByRep`): O(1) (atualiza contadores e faz splice das listas).
 * - Iteração: O(|S|) proporcional ao número de pixels/sets percorridos.
 *
 * ## Invariantes
 * - Se `indexToPixel[i] == -1` então o slot do conjunto `i` é inválido (após fusões).
 * - Para qualquer representante `rep`: `pixelToIndex[rep]` aponta para um índice `i`
 *   tal que `indexToPixel[i] == headRep` e `sizeSets[i] > 0`.
 * - As listas de pixels de um mesmo conjunto formam um ciclo por `pixelsNext`.
 *
 * ## Exemplo mínimo
 * @code
 * PixelSetManager psm(numPixels);
 * // inicialização dos reps/estruturas omitida
 * int a = repA, b = repB;
 * psm.mergeSetsByRep(a, b);       // b fundido em a
 * for (int px : psm.getPixelsBySet(a)) {
 *   // processa pixels do conjunto resultante
 * }
 * for (int rep : psm.getFlatzoneRepresentatives()) {
 *   // percorre reps válidos
 * }
 * @endcode
 */
struct PixelSetManager{
    
    std::vector<int> pixelToIndex; //mapeamento do pixel representante para índice na lista de conjuntos disjuntos. Tamanho: numPixels
    std::vector<int> indexToPixel; // mapeamento de índice para pixel representante. Tamanho: numSets
    std::vector<int> sizeSets; // usada para armazenar o tamanho dos conjuntos disjuntos. Tamanho: numSets
    std::vector<int> pixelsNext; // mapa de pixels dos conjuntos disjuntos: Tamanho: numPixels
    int activeSetsCount;

    PixelSetManager(int numPixels, int numSets)
        : pixelToIndex(numPixels, -1),
          indexToPixel(numSets, -1),
          sizeSets(numSets, 0),
          pixelsNext(numPixels, -1),
          activeSetsCount(numSets) { }
        
    PixelSetManager(int numPixels)
        : pixelToIndex(numPixels, -1),
          indexToPixel(numPixels, -1),
          sizeSets(numPixels, 0),
          pixelsNext(numPixels, -1),
          activeSetsCount(numPixels) { }
    
    int numSets() const { return activeSetsCount; }

    int numPixelsInSet(int rep){ return sizeSets[pixelToIndex[rep]]; }

    int numPixelsInSets(const std::vector<int>& reps){
        int sum = 0;
        for (int rep : reps) {
            sum += sizeSets[pixelToIndex[rep]];
        }
        return sum;
    }
    
    int indexOfPixel(int pixel) const {
        return pixelToIndex[pixel];
    }
    
    int pixelOfIndex(int idx) const {
        return indexToPixel[idx];
    }

    /**
     * @brief Redimensiona os vetores relacionados a conjuntos (flat zones)
     * para refletir o número real de FZs criadas.
     *
     * @param newNumSets Número real de conjuntos encontrados.
     */
    void shrinkToNumSets(int newNumSets) {
        indexToPixel.resize(newNumSets);
        sizeSets.resize(newNumSets);
        activeSetsCount = newNumSets;
    }
    
    void mergeSetsByRep(int repWinner, int repLoser) {
        
        // 1. Recupera índices dos representantes
        int idxRootWinner = pixelToIndex[repWinner];
        int idxRootLoser  = pixelToIndex[repLoser];
        if (idxRootWinner < 0 || idxRootLoser < 0 || idxRootWinner == idxRootLoser) return;

        sizeSets[idxRootWinner] += sizeSets[idxRootLoser];

        // 2. Splice O(1) das listas circulares (pixels)
        int nextWinner = pixelsNext[repWinner];
        int nextLoser  = pixelsNext[repLoser];
        pixelsNext[repWinner] = nextLoser;
        pixelsNext[repLoser]  = nextWinner;

        // 3. Invalida slot perdedor
        sizeSets[idxRootLoser]  = 0;
        indexToPixel[idxRootLoser] = -1;

        // 4. Redireciona lookups pelo antigo rep pixel
        pixelToIndex[repLoser] = idxRootWinner;
        if (activeSetsCount > 0) --activeSetsCount;
    }


    struct View {
        std::span<int> pixelToIndex;
        std::span<int> indexToPixel;
        std::span<int> sizeSets;
        std::span<int> pixelsNext;
    };

    View view() noexcept {
        return {std::span<int>(pixelToIndex), std::span<int>(indexToPixel), std::span<int>(sizeSets), std::span<int>(pixelsNext) };
    }

    
    std::span<int> viewOfPixelToIndex(){ return std::span<int>(pixelToIndex); }
    std::span<int> viewOfIndexToPixel(){ return std::span<int>(indexToPixel); }
    std::span<int> viewOfSizeSets(){ return std::span<int>(sizeSets); }
    std::span<int> viewOfPixelsNext(){ return std::span<int>(pixelsNext); }
    
    
    /**
     * @brief Faixa iterável sobre os pixels de uma ou mais sets.
     *
     * Esta classe encapsula um range de representantes de sets (ex.: `int`, 
     * `std::vector<int>`, `std::span<const int>`, `RepsOfCCRange` etc.) e fornece 
     * iteradores (`PixelsBySetIterator`) que percorrem todos os pixels desses sets. 
     *
     * O funcionamento é baseado nos seguintes princípios:
     *  - Cada representante identifica um set.
     *  - O iterador avança sobre todos os pixels de cada set usando a lista circular 
     *    interna (`pixelsNext`) e o tamanho registrado em `sizeSets`.
     *  - Ao término de um set, o iterador passa automaticamente para o próximo 
     *    set indicada no range de representantes.
     *  - O range é apenas uma "view" sobre os reps. Ele não copia os pixels, apenas 
     *    percorre dinamicamente as listas já armazenadas.
     *
     * Exemplos de uso:
     * @code
     *
     * // 1) Único representante
     * for (int px : v.getPixelsBySet(rep)) {
     *     // processa cada pixel do set de 'rep'
     * }
     *
     * // 2) Vários reps (std::vector<int>)
     * std::vector<int> reps = {rep1, rep2};
     * for (int px : v.getPixelsBySet(reps)) {
     *     // processa pixels dos sets de rep1 e rep2
     * }
     *
     * // 3) Usando span (sem cópia do vetor)
     * std::span<const int> s(reps.data(), reps.size());
     * for (int px : v.getPixelsBySet(s)) {
     *     // processa pixels, sem overhead de cópia
     * }
     *
     * // 4) Usando um range custom (ex.: RepsOfCCRange de um NodeCT)
     * for (int px : g.getPixelsBySet(node->getRepsOfCC())) {
     *     // processa pixels de todos os sets alcançadas na BFS
     * }
     * @endcode
     *
     * @tparam Range Tipo do container ou view que contém os reps.
     *               Pode ser `std::array<int,N>`, `std::vector<int>`, 
     *               `std::span<const int>`, ou um range custom compatível
     *               com `std::begin`/`std::end` que produza `int`.
     */
    template<class Range>
    class PixelsBySetRange {
    private:
        PixelSetManager::View   v_;  // guarda o View (spans) por valor
        Range reps_;   // range de representantes

        using RepIt = decltype(std::begin(std::declval<const Range&>()));

    public:
        // -------- Iterator ----------
        class PixelsBySetIterator {
            PixelSetManager::View   v_;
            RepIt it_, last_;
            int   cur_{-1};      // pixel atual
            int   remaining_{0}; // pixels restantes na FZ atual

            void startNextSegment() {
                cur_ = -1;
                remaining_ = 0;
                while (it_ != last_) {
                    const int rep = *it_;
                    const int idx = v_.pixelToIndex[rep];
                    if (idx >= 0) {
                        const int head = v_.indexToPixel[idx];
                        const int sz   = v_.sizeSets[idx];
                        if (head != -1 && sz > 0) {
                            cur_ = head;
                            remaining_ = sz;
                            return;
                        }
                    }
                    ++it_; // tenta próximo representante
                }
            }

        public:
            using iterator_category = std::forward_iterator_tag;
            using value_type        = int;
            using difference_type   = std::ptrdiff_t;
            using pointer           = const int*;
            using reference         = const int&;

            PixelsBySetIterator(PixelSetManager::View v, RepIt first, RepIt last) : v_(v), it_(first), last_(last) {
                startNextSegment();
            }

            reference operator*()  const { return cur_; }
            pointer   operator->() const { return &cur_; }

            PixelsBySetIterator& operator++() {
                if (remaining_ > 0) {
                    --remaining_;
                    if (remaining_ == 0) {
                        ++it_;
                        startNextSegment();        // próxima FZ
                    } else {
                        cur_ = v_.pixelsNext[cur_]; // próximo pixel na FZ
                    }
                }
                return *this;
            }

            PixelsBySetIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

            friend bool operator==(const PixelsBySetIterator& a, const PixelsBySetIterator& b) {
                return a.cur_ == b.cur_ && a.it_ == b.it_ && a.last_ == b.last_;
            }
            friend bool operator!=(const PixelsBySetIterator& a, const PixelsBySetIterator& b) {
                return !(a == b);
            }

        };

        // -------- Range ----------
        PixelsBySetRange(PixelSetManager::View v, Range reps)
            : v_(v), reps_(std::move(reps)) {}

        PixelsBySetIterator begin() const { return PixelsBySetIterator(v_, std::begin(reps_), std::end(reps_)); }
        PixelsBySetIterator end()   const { return PixelsBySetIterator(v_, std::end(reps_),  std::end(reps_));  }
    };

    //Um único representante
    auto getPixelsBySet(int rep) {
        return PixelsBySetRange<std::array<int,1>>(this->view(), std::array<int,1>{rep});
    }

    // Qualquer range (vector<int>, span<const int>, RepsOfCCRange, …)
    template<class Range>
    auto getPixelsBySet(Range reps) {
        return PixelsBySetRange<Range>(this->view(), std::move(reps));
    }




    /**
     * @brief Iterador para percorrer todos os representantes de sets válidos.
     *
     * Este iterador percorre o array `indexToPixel`, saltando entradas inválidas
     * (marcadas com -1). Cada elemento retornado é o pixel representante da set.
     *
     * Uso típico:
     * @code
     * for (int rep : getFlatzoneRepresentatives()) {
     *     // rep é um representante válido
     * }
     * @endcode
     */
    class RepresentativeIterator {
    private:
        std::span<int> indexToPixel_;
        size_t size_;
        size_t idx_;

        void skipInvalid() {
            while (idx_ < size_ && indexToPixel_[idx_] == -1) {
                ++idx_;
            }
        }

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type        = int;
        using difference_type   = std::ptrdiff_t;
        using pointer           = const int*;
        using reference         = const int&;

        RepresentativeIterator(std::span<int> data, size_t size, size_t startIdx)
            : indexToPixel_(data), size_(size), idx_(startIdx) {
            skipInvalid();
        }

        reference operator*()  const { return indexToPixel_[idx_]; }
        pointer   operator->() const { return &indexToPixel_[idx_]; }

        RepresentativeIterator& operator++() {
            ++idx_;
            skipInvalid();
            return *this;
        }

        RepresentativeIterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const RepresentativeIterator& a, const RepresentativeIterator& b) {
            return a.idx_ == b.idx_;
        }
        friend bool operator!=(const RepresentativeIterator& a, const RepresentativeIterator& b) {
            return !(a == b);
        }
    };


    /**
     * @brief Faixa iterável de representantes de flat zones válidos.
     *
     * Retorna um objeto que pode ser usado em range-based for loops
     * para iterar apenas sobre os representantes ativos.
     */
    class RepresentativeRange {
    private:
        std::span<int> indexToPixel_;
        size_t size_;

    public:
        explicit RepresentativeRange(std::span<int> data, size_t size)
            : indexToPixel_(data), size_(size) {}

        RepresentativeIterator begin() const { return RepresentativeIterator(indexToPixel_, size_, 0); }
        RepresentativeIterator end()   const { return RepresentativeIterator(indexToPixel_, size_, size_); }
    };

    /**
     * @brief Obtém um range iterável de representantes de flat zones válidos.
     *
     * Uso:
     * @code
     * for (int rep : getFlatzoneRepresentatives()) {
     *     // processar rep
     * }
     * @endcode
     */
    RepresentativeRange getFlatzoneRepresentatives()  {
        return RepresentativeRange(viewOfIndexToPixel(), indexToPixel.size());
    }


};


/**
 * @brief Fila linear baseada em std::vector para alto desempenho em BFS.
 *
 * Essa estrutura encapsula um std::vector<T> e um índice de leitura (`head_`),
 * funcionando como uma fila FIFO.  
 * Ao contrário de std::queue, não tem overhead de alocações/push/pop, pois:
 *   - `push` adiciona ao fim do vetor.
 *   - `pop` apenas avança o índice `head_`.
 *   - `clear` reseta o índice e o tamanho para reutilização sem desalocar memória.
 *
 * Uso típico:
 * @code
 * FastQueue<int> q;
 * q.reserve(2048);
 * q.push(10);
 * q.push(20);
 * while (!q.empty()) {
 *     int x = q.pop();
 *     // processa x
 * }
 * @endcode
 */
template <typename T>
struct FastQueue {
private:
    std::vector<T> data_;
    size_t head_ = 0; // índice do próximo elemento a ser lido

public:
    FastQueue() = default;

    FastQueue(size_t n){
        data_.reserve(n); 
    } 

    /// Reserva espaço inicial (opcional, para evitar realocações)
    void reserve(size_t n) { data_.reserve(n); }

    /// Remove todos os elementos, reseta o índice de leitura
    void clear() { data_.clear(); head_ = 0; }

    /// Retorna se a fila está vazia
    bool empty() const { return head_ >= data_.size(); }

    /// Retorna o tamanho atual da fila
    size_t size() const { return data_.size() - head_; }

    /// Adiciona um elemento ao fim
    void push(const T& value) { data_.push_back(value); }

    void push(T&& value) { data_.push_back(std::move(value)); }

    /// Remove e retorna o próximo elemento
    T pop() { return std::move(data_[head_++]); }

    /// Acesso ao próximo elemento sem remover
    T& front() { return data_[head_]; }
    const T& front() const { return data_[head_]; }
};

/**
 * @brief Pilha (stack) simples e performática baseada em `std::vector`.
 *
 * `FastStack<T>` provê a interface essencial de uma pilha LIFO com
 * operações de custo amortizado O(1) e controle de capacidade via
 * `reserve`. É útil em rotinas de DFS, processamento de componentes,
 * e estruturas auxiliares onde o overhead de `std::stack` e alocações
 * frequentes deve ser evitado.
 *
 * ## Operações
 * - `push(const T&)`, `push(T&&)` — insere no topo (amortizado O(1)).
 * - `pop()` — remove e retorna o topo (amortizado O(1)).
 * - `top()` — acesso ao elemento do topo (O(1)).
 * - `empty()`, `size()` — consultas O(1).
 * - `reserve(n)`, `clear()` — gestão de capacidade e limpeza.
 *
 * ## Exemplo de uso
 * @code
 * FastStack<int> st;
 * st.reserve(1024);
 * st.push(3);
 * st.push(7);
 * int x = st.top();   // 7
 * x = st.pop();       // 7; agora o topo é 3
 * @endcode
 */
template <typename T>
struct FastStack {
private:
    std::vector<T> data_;

public:
    FastStack() = default;

    explicit FastStack(size_t n) {
        data_.reserve(n);
    }

    /// Reserva espaço inicial (opcional)
    void reserve(size_t n) { data_.reserve(n); }

    /// Remove todos os elementos
    void clear() { data_.clear(); }

    /// Retorna se a pilha está vazia
    bool empty() const { return data_.empty(); }

    /// Retorna o tamanho atual da pilha
    size_t size() const { return data_.size(); }

    /// Adiciona um elemento ao topo
    void push(const T& value) { data_.push_back(value); }

    void push(T&& value) { data_.push_back(std::move(value)); }

    /// Remove e retorna o elemento do topo
    T pop() {
        T value = std::move(data_.back());
        data_.pop_back();
        return value;
    }

    /// Acesso ao topo sem remover
    T& top() { return data_.back(); }
    const T& top() const { return data_.back(); }
};

/**
 * @brief Prefiltro local de adjacências com 64 buckets (1 máscara de 64 bits).
 *
 * Ajuda a evitar inserir repetidamente o mesmo vizinho `b` enquanto
 * você processa os pixels de borda de uma única flat zone `a`.
 *
 * Estratégia:
 *  - Uma máscara de 64 bits (`mask`) indica quais buckets já foram vistos.
 *  - Uma lista curta `small[64]` guarda os valores exatos para resolver
 *    colisões de hash (quando dois valores caem no mesmo bucket).
 *
 * Uso típico:
 *  @code
 *  LocalPrefilter64 pf;
 *  ...
 *  if (idxP != pf.currentFZ) pf.reset(idxP); // ao mudar de FZ base
 *  if (!pf.contains(idxQ)) {
 *      pf.insert(idxQ);
 *      // emitir aresta (apenas um lado); o espelhamento virá depois
 *  }
 *  @endcode
 *
 * Complexidade:
 *  - reset: O(1)
 *  - contains: O(1) sem colisão; O(k) na lista curta quando há colisão (k pequeno)
 *  - insert: O(1)
 *
 * Observações:
 *  - Nunca dá falso negativo (se `contains(x)` é true, x foi inserido).
 *  - Pode haver falso positivo no bucket, mas a lista curta confirma exatidão.
 *  - Não é thread-safe; use uma instância por thread.
 */
struct LocalPrefilter64 {
    uint64_t mask = 0ull;   // 64 buckets (bitset)
    int      small[64];     // valores para confirmar colisões
    int      sz = 0;        // tamanho efetivo de 'small'
    int      currentFZ = -1;

    /// Reinicia o filtro para uma nova FZ base.
    inline void reset(int fz) noexcept {
        mask = 0ull;
        sz = 0;
        currentFZ = fz;
    }

    /// Retorna true se 'x' já foi visto nesta FZ.
    inline bool contains(int x) const noexcept {
        const uint64_t bit = 1ull << (unsigned(x) & 63u);
        if ((mask & bit) == 0ull) return false;   // bucket vazio ⇒ certeza de não visto
        // bucket marcado: confirmar na lista curta
        for (int i = 0; i < sz; ++i) {
            if (small[i] == x) return true;
        }
        return false;


    }

    /// Marca 'x' como visto (atualiza bucket e lista curta).
    inline void insert(int x) noexcept {
        const uint64_t bit = 1ull << (unsigned(x) & 63u);
        mask |= bit;
        if (sz < 64) small[sz++] = x; // truncamento conservador acima de 64 entradas
    }
};

/**
 * Esta estrutura implementa a interface essencial de um “set de inteiros” usando
 * um `std::vector<int>` como armazenamento subjacente.
 *
 * **Dois aceleradores internos**:
 *  1. *Modo ordenado* (`sorted_`): após `finalize()`, o vetor é ordenado e deduplicado;
 *     buscas usam `lower_bound` (O(log d)) e inserções/remoções preservam a ordenação.
 *  2. *Tiny Bloom* (64 bits, 2 hashes): pré-teste O(1) para descartar buscas negativas
 *     sem tocar no vetor. Ativado somente quando `size() >= kBloomMin` (padrão 12).
 *     O Bloom **não** gera falsos negativos (apenas alguns falsos positivos) e é
 *     reconstruído em `finalize()`. Remoções não “desmarcam” bits (conservador).
 *
 * **API exposta (drop-in para os seus métodos críticos)**:
 *  - Iteração: `begin()/end()` (const e não-const).
 *  - Consulta: `find(int)` (const e não-const) → iterador ou `end()`.
 *  - Atualização: `insert(int)`, `erase(int)`, `swap(...)`, `reserve(size_t)`.
 *  - Manutenção: `finalize()` (ordena + deduplica + reconstrói Bloom).
 *  - Estado: `size()`, `empty()`.
 *
 * **Complexidade amortizada** (\f$d =\f$ grau do nó):
 *  - `find(x)`: 
 *      - sem ordenação: O(d) (linear desenrolado) — Bloom pode abortar em O(1) para misses;
 *      - com ordenação: O(log d) via `lower_bound` — Bloom idem.
 *  - `insert(x)`:
 *      - sem ordenação: O(d) (checa duplicata) + O(1) push_back;
 *      - com ordenação: O(d) no pior caso (inserção estável) com custo pequeno em graus baixos.
 *  - `erase(x)`:
 *      - sem ordenação: O(d) (swap-pop, não preserva ordem);
 *      - com ordenação: O(d) (erase preservando ordem).
 *  - `finalize()`: O(d log d) para ordenar + O(d) para deduplicar e reconstruir Bloom.
 *
 * **Corretude & Invariantes**
 *  - Representa um conjunto de inteiros (sem duplicatas) após `finalize()`.
 *  - `find(x)!=end()` ⇒ `x` está presente; Bloom **nunca** introduz falso negativo.
 *  - Em modo não-ordenado, `insert` evita duplicatas por busca linear; `finalize()`
 *    reforça unicidade e define `sorted_=true`.
 *
 * **Quando usar `finalize()`**:
 *  - Ao término da construção inicial das adjacências do grafo.
 *  - Após **lotes** de fusões; evite chamar por operação individual.
 *
 * **Notas práticas**
 *  - Chame `reserve(k)` após `resize(numFZ)` da lista de sets; use k≈12–24 conforme o p95 do grau.
 *  - O Bloom é aplicado apenas quando `size() >= kBloomMin` para evitar overhead em listas muito pequenas.
 *  - `swap` permite o “swap trick” para liberar memória em O(1) (opcional no seu fluxo).
 *
 * @warning Não é *thread-safe*; proteja externamente se houver acesso concorrente.
 * @see FlatZonesGraphFullEdges::getAdjacentFlatzonesFromPixel
 * @see FlatZonesGraphFullEdges::mergeAdjacentCandidatesInPlace
 * @see FlatZonesGraphFullEdges::mergeBasesWithAdjacentCandidatesInPlace
 *
 * @par Exemplo
 * @code
 * AdjacentFlatZonesSet adj;
 * adj.reserve(16);
 * adj.insert(10);
 * adj.insert(7);
 * adj.finalize();                 // ordena: {7,10}
 * if (adj.find(8) == adj.end()) { // Bloom pode abortar em O(1)
 *   adj.insert(8);                // insere preservando ordenação: {7,8,10}
 * }
 * adj.erase(8);                   // remove preservando ordenação
 * @endcode
 */
class AdjacentFlatZonesSet {
    std::vector<int> v_;
    bool     sorted_{false};   // true => v_ ordenado e sem duplicatas
    

    /**
     * @brief Bloom filter mínimo (64 bits, 2 hashes) para acelerar consultas negativas.
     *
     * Sem falsos negativos (se diz que não tem, não tem); com raros falsos positivos.
     * Ideal para compor estruturas como `AdjacentFlatZonesSet` após `finalize()`.
     *
     * Uso:
     *  - `add(x)` para marcar,
     *  - `maybeHas(x, sizeHint)` para peneirar consultas (O(1)),
     *  - `clear()` para zerar,
     *  - `rebuild(span)` para reconstruir a partir de um vetor ordenado.
     */
    struct TinyBloom64 {
        uint64_t bits{0};

        // splitmix-like hashes (rápidos e com boa dispersão)
        static inline uint32_t h1(uint32_t x) noexcept {
            x ^= x >> 16; x *= 0x7feb352dU; x ^= x >> 15; x *= 0x846ca68bU; x ^= x >> 16; return x;
        }
        static inline uint32_t h2(uint32_t x) noexcept {
            x ^= x >> 15; x *= 0x2c1b3c6dU; x ^= x >> 12; x *= 0x297a2d39U; x ^= x >> 15; return x;
        }
        static inline uint64_t bit(uint32_t h) noexcept { return 1ull << (h & 63u); }

        inline void clear() noexcept { bits = 0; }

        inline void add(int v) noexcept {
            uint32_t a = h1((uint32_t)v), b = h2((uint32_t)v);
            bits |= bit(a) | bit(b);
        }

        /// Retorna false ⇒ certeza de NÃO conter; true ⇒ "talvez".
        inline bool maybeHas(int v, int sizeHint) const noexcept {
            if (sizeHint < 8) return true; // desliga em listas pequenas (ajustável)
            uint32_t a = h1((uint32_t)v), b = h2((uint32_t)v);
            uint64_t m = bit(a) | bit(b);
            return (bits & m) == m;
        }

        template<class RangeLike>
        inline void rebuild(const RangeLike& sortedUnique) noexcept {
            clear();
            for (int x : sortedUnique) add(x);
        }
    };
    TinyBloom64 bloom_;

    inline bool bloom_maybe_has(int x) const noexcept {
        return bloom_.maybeHas(x, (int)v_.size());
    }

public:
    using value_type     = int;
    using iterator       = std::vector<int>::iterator;
    using const_iterator = std::vector<int>::const_iterator;

    AdjacentFlatZonesSet() = default;

    // capacidade
    void reserve(size_t n) { v_.reserve(n); }

    // tamanho/estado
    int  size()  const { return static_cast<int>(v_.size()); }
    bool empty() const { return v_.empty(); }
    
    bool isSorted() const noexcept { return sorted_; }

    // Força operações em modo não-ordenado (evita shifts em vector::insert/erase).
    // Mantém o Bloom (pode ter falsos positivos, mas nunca falsos negativos).
    void markUnsorted() noexcept { sorted_ = false; }

    // Reconstrói o Bloom a partir do conteúdo atual (ordem não importa).
    void rebuildBloom() { bloom_.rebuild(v_); }
    void clear() {
        v_.clear();
        sorted_ = false;
        bloom_.clear();
    }
    // iteração
    iterator begin() { return v_.begin(); }
    iterator end()   { return v_.end();   }
    const_iterator begin() const { return v_.begin(); }
    const_iterator end()   const { return v_.end();   }

    // ---- find (const) ----
    const_iterator find(int x) const {
        if (!bloom_.maybeHas(x, (int)v_.size())) return v_.cend();
        if (!sorted_) {
            for (auto it = v_.cbegin(), e = v_.cend(); it != e; ++it) if (*it == x) return it;
            return v_.cend();
        } else {
            auto it = std::lower_bound(v_.cbegin(), v_.cend(), x);
            return (it != v_.cend() && *it == x) ? it : v_.cend();
        }
    }

    // ---- find (não-const) ----
    iterator find(int x) {
        if (!bloom_.maybeHas(x, (int)v_.size())) return v_.end();
        if (!sorted_) {
            for (auto it = v_.begin(), e = v_.end(); it != e; ++it) if (*it == x) return it;
            return v_.end();
        } else {
            auto it = std::lower_bound(v_.begin(), v_.end(), x);
            return (it != v_.end() && *it == x) ? it : v_.end();
        }
    }

    // ---- insert ----
    void insert(int x) {
        if (sorted_) {
            auto it = std::lower_bound(v_.begin(), v_.end(), x);
            if (it == v_.end() || *it != x) {
                v_.insert(it, x);
                bloom_.add(x);           // <<<< marca no bloom
            }
        } else {
            // construção não-ordenada
            for (int y : v_) if (y == x) return;
            v_.push_back(x);
            // opcional: não marcar no bloom aqui (como antes), ou:
            bloom_.add(x);               // pode marcar, sem custo relevante
        }
    }

    void appendUnchecked(int x) { v_.push_back(x); }

    // ---- finalize: ordena, deduplica, reconstrói Bloom ----
    void finalize(bool rebuildBloom = true) {
        std::sort(v_.begin(), v_.end());
        v_.erase(std::unique(v_.begin(), v_.end()), v_.end());
        sorted_ = true;
        if (rebuildBloom) bloom_.rebuild(v_); 
    }

    // Ordena apenas quando o tamanho justifica; caso contrário só reconstrói Bloom.
    void finalizeMaybeSorted(int sortThreshold, bool rebuildBloom = true) {
        if ((int)v_.size() >= sortThreshold) {
            finalize(rebuildBloom);
            return;
        }
        sorted_ = false;
        if (rebuildBloom) bloom_.rebuild(v_);
    }

    
    // ---- erase ----
    void erase(int x) {
        if (sorted_) {
            // remove preservando ordenação
            auto it = std::lower_bound(v_.begin(), v_.end(), x);
            if (it != v_.end() && *it == x) v_.erase(it);
            return;
        }
        // modo não-ordenado: swap-pop (O(1))
        auto it = find(x);
        if (it != end()) {
            *it = v_.back();
            v_.pop_back();
        }
    }

    // ---- swap ----
    void swap(AdjacentFlatZonesSet& other) noexcept {
        v_.swap(other.v_);
        std::swap(sorted_, other.sorted_);
        std::swap(bloom_,  other.bloom_);
    }


    /**
     * @brief Espelha arestas não-direcionadas e finaliza todas as listas.
     *
     * Uso típico:
     *  - No loop de construção, adicione **apenas um lado** (regra idxP < idxQ) via appendUnchecked().
     *  - Ao final, chame mirrorAndFinalize(listOfAdjacentFlatZones, pixelView).
     *
     * Passos:
     *  1) Para cada FZ a, faz uma deduplicação local dos reps (sem tocar no Bloom),
     *  2) Espelha a aresta para B (appendUnchecked do repA em B),
     *  3) Finaliza tudo: sort+unique e reconstrução do Bloom.
     *
     * @param lists vetor de listas de adjacência por FZ (representantes inteiros).
     * @param v     PixelSetManager::View com mapeamentos (indexToPixel, pixelToIndex, ...).
     */
    static inline void mirrorAndFinalize(std::vector<AdjacentFlatZonesSet>& listOfAdjacentFlatZones, PixelSetManager::View pixelView){
        
        for (auto& adjFZ : listOfAdjacentFlatZones) adjFZ.finalize(false);
        
        // espelha as arestas para o outro lado
        for (int a = 0; a < (int)listOfAdjacentFlatZones.size(); ++a) {
            const int repA = pixelView.indexToPixel[a];
            for (int repB : listOfAdjacentFlatZones[a]) {
                const int b = pixelView.pixelToIndex[repB];
                // só faltará o lado b->a (garantido pela regra idxP<idxQ)
                listOfAdjacentFlatZones[b].appendUnchecked(repA);
            }
        }

        // finalize final (ordena/dedup + reconstrói Bloom)
        for (auto& adjFZ : listOfAdjacentFlatZones) adjFZ.finalize(true);
        
    }
};

using AdjacentFlatZones = AdjacentFlatZonesSet; //Tipo de dado que armazena as aresta do grafo




class Stopwatch {
public:
    using clock = std::chrono::steady_clock;

    void start() {
        accumulated_ = clock::duration::zero();
        running_ = true;
        last_ = clock::now();
    }

    void pause() {
        if (!running_) return;
        accumulated_ += clock::now() - last_;
        running_ = false;
    }

    void resume() {
        if (running_) return;
        running_ = true;
        last_ = clock::now();
    }

    // tempo decorrido "ativo" até agora
    clock::duration elapsed() const {
        if (!running_) return accumulated_;
        return accumulated_ + (clock::now() - last_);
    }

    bool running() const { return running_; }

private:
    clock::time_point last_{};
    clock::duration accumulated_{clock::duration::zero()};
    bool running_{false};
};
