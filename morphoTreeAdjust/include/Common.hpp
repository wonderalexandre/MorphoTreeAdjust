#ifndef COMMONS_HPP  
#define COMMONS_HPP  


#define NDEBUG  // Remove os asserts do código
#include <cassert>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <memory>
#include <limits>

#define PRINT_LOG 0 
#define PRINT_DEBUG 0 

// Forward declaration dos templates
template <typename T>
class ComponentTree;

template <typename T>
class NodeCT;




// Definição de tipos para CNPs
using Pixels = std::list<int>;                 // Representa uma lista de pixels (CNPs)
using FlatZone = std::list<int>;                 // Representa uma flatzone
using FlatZones = std::unordered_map<int, FlatZone>;   // Representa uma lista de flatzones (CNPs ficarão separados em flatzones)

//Alias em função do tipo dos CNPs
using FlatZoneRef = std::reference_wrapper<FlatZone>;   // Representa uma flatzone
using FlatZonesRef = std::list<FlatZoneRef>;   // Representa uma lista de flatzones

using ComponentTreeFZ = ComponentTree<FlatZones>; //representa uma component tree por flatzones
using ComponentTreeP = ComponentTree<Pixels>; //representa uma component tree sem tratamento de flatzones
using NodeFZ = NodeCT<FlatZones>; //representa um node com separação dos cnps em flatzones
using NodeP = NodeCT<Pixels>; //representa um node sem tratamento de flatzones


template <typename CNPsType>
using NodeCTPtr = std::shared_ptr<NodeCT<CNPsType>>;
using NodeFZPtr = std::shared_ptr<NodeCT<FlatZones>>;
using NodePPtr  = std::shared_ptr<NodeCT<Pixels>>;

template <typename CNPsType>
using ComponentTreePtr = std::shared_ptr<ComponentTree<CNPsType>>; 
using ComponentTreeFZPtr = std::shared_ptr<ComponentTreeFZ>; //representa um component tree por flatzones
using ComponentTreePPtr = std::shared_ptr<ComponentTreeP>; //representa um component tree sem tratamento de flatzones


using AdjacentFlatzones = std::unordered_set<int>;
using ListOfAdjacentFlatzones = std::vector<std::unique_ptr<AdjacentFlatzones>>;


struct FlatZoneNode {
    NodeFZPtr node = nullptr;
    FlatZone* flatzone = nullptr;
    int idFlatZone;

    FlatZoneNode(){}

    // Construtor para mover a FlatZone
    FlatZoneNode(NodeFZPtr n, FlatZone& fz) : node(n), flatzone(&fz), idFlatZone(fz.front()) {} 
    

};




/*
Uma classe de imagem genérica que pode ser usada para armazenar imagens de diferentes tipos de pixel.
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


#endif 
