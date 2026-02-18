#pragma once

#include <iterator>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <memory>
#include <unordered_map>
#include <vector>
#include <span>

#include "../include/AdjacencyRelation.hpp"
#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/Common.hpp"

enum class Attribute {
    AREA,
    VOLUME,
	RELATIVE_VOLUME,
    LEVEL,
    GRAY_HEIGHT,
    
    BOX_WIDTH,
    BOX_HEIGHT,
	DIAGONAL_LENGTH,
    
    INERTIA
};
using enum Attribute;


/**
 * @brief Funções utilitárias para computar atributos em árvores de forma incremental.
 *
 * A classe expõe algoritmos genéricos de travessia pós-ordem permitindo
 * compor etapas de pré-processamento, mesclagem e pós-processamento sem criar
 * estruturas auxiliares temporárias. 
 */
template <typename CNPsType, typename GraphT = DefaultFlatZonesGraph>
class AttributeComputedIncrementally {
public:
    
    template<class PreProcessing, class MergeProcessing, class PostProcessing>
    static void computerAttribute(ComponentTree<CNPsType, GraphT>* tree, NodeId root, 
                                    PreProcessing&& preProcessing, 
                                    MergeProcessing&& mergeProcessing, 
                                    PostProcessing&& postProcessing) {
        preProcessing(root);
        for (NodeId child : tree->getChildrenById(root)) {
            AttributeComputedIncrementally::computerAttribute(tree, child, preProcessing, mergeProcessing, postProcessing); // passar por ref (sem cópia)
            mergeProcessing(root, child);
        }
        postProcessing(root);
    }


    
};



/**
 * @brief Resolve nomes de atributos em índices lineares usados pelos buffers.
 */
class AttributeNames {
public:
    std::unordered_map<Attribute, int> indexMap;
    const int NUM_ATTRIBUTES;

    AttributeNames(std::unordered_map<Attribute, int>&& map)
        : indexMap(std::move(map)), NUM_ATTRIBUTES(static_cast<int>(indexMap.size())) {}

    static AttributeNames fromList(int n, const std::vector<Attribute>& attributes) {
        std::unordered_map<Attribute, int> map;
        int i = 0;
        for (auto attr : attributes) {
            map[attr] = i++ * n;
        }
        return AttributeNames(std::move(map));
    }

    int getIndex(Attribute attr) const {
        return indexMap.at(attr);
    }

    int linearIndex(int nodeIndex, Attribute attr) const {
        return nodeIndex * NUM_ATTRIBUTES + getIndex(attr);
    }
};

template<typename PreProcessing, typename MergeProcessing, typename PostProcessing>
struct IncrementalFunctions {
    PreProcessing preProcessing;
    MergeProcessing mergeProcessing;
    PostProcessing postProcessing;
};

template<typename CNPsType, typename GraphT, typename IncrementalFunctionsT, typename Derived>
class AttributeComputer {
protected:
    static_assert(!std::is_same_v<CNPsType, FlatZones> || FlatZonesGraphCommonInterface<GraphT>,
                  "GraphT must satisfy FlatZonesGraphCommonInterface when CNPsType is FlatZones");
    ComponentTree<CNPsType, GraphT>* tree;
    Attribute attribute;
    IncrementalFunctionsT incrementalFunctions;
    mutable Derived* selfPtr = nullptr;

public:
    AttributeComputer(ComponentTree<CNPsType, GraphT>* tree,
                      Attribute attribute,
                      IncrementalFunctionsT incrementalFunctions)
        : tree(tree), attribute(attribute), incrementalFunctions(std::move(incrementalFunctions)) {}

    virtual ~AttributeComputer() = default;

    ComponentTree<CNPsType, GraphT>* getTree() { return tree; }
    const ComponentTree<CNPsType, GraphT>* getTree() const { return tree; }

    // Permite registrar explicitamente o ponteiro para o objeto derivado ao final do ctor do derivado.
    void registerDerived(Derived* self) const {
        selfPtr = self;
    }

    std::vector<float> compute() {
        std::vector<float> buffer(tree->getNumNodes(), 0.0f);
        compute(buffer);
        return buffer;
    }

    void compute(std::span<float> buffer) {
        
        AttributeComputedIncrementally<CNPsType, GraphT>::computerAttribute(
            tree, tree->getRootById(),
            [this, buffer](NodeId idx) mutable { preProcessing(idx, buffer); },
            [this, buffer](NodeId parentId, NodeId childId) mutable { mergeProcessing(parentId, childId, buffer); },
            [this, buffer](NodeId idx) mutable { postProcessing(idx, buffer); }
        );
        
    }

    Attribute getAttribute() const {
        return attribute;
    }

    // Acesso direto aos functors configurados, com as assinaturas desejadas.
    void preProcessing(NodeId idx, std::span<float> buf) {
        auto self = selfPtr ? selfPtr : static_cast<Derived*>(this);
        (self->*incrementalFunctions.preProcessing)(idx, buf);
    }

    void mergeProcessing(NodeId parentId, NodeId childId, std::span<float> buf) {
        auto self = selfPtr ? selfPtr : static_cast<Derived*>(this);
        (self->*incrementalFunctions.mergeProcessing)(parentId, childId, buf);
    }

    void postProcessing(NodeId idx, std::span<float> buf) {
        auto self = selfPtr ? selfPtr : static_cast<Derived*>(this);
        (self->*incrementalFunctions.postProcessing)(idx, buf);
    }
};

/**
 * @brief Computa a área (número de pixels) de cada nó da árvore.
 */
template<typename CNPsType, typename GraphT = DefaultFlatZonesGraph>
class AreaComputer : public AttributeComputer<
    CNPsType,
    GraphT,
    IncrementalFunctions<
        void (AreaComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>),
        void (AreaComputer<CNPsType, GraphT>::*)(NodeId, NodeId, std::span<float>),
        void (AreaComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>)>,
    AreaComputer<CNPsType, GraphT>> {

public:
    using Functions = IncrementalFunctions<
        void (AreaComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>),
        void (AreaComputer<CNPsType, GraphT>::*)(NodeId, NodeId, std::span<float>),
        void (AreaComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>)>;
    using Base = AttributeComputer<CNPsType, GraphT, Functions, AreaComputer<CNPsType, GraphT>>;

    Functions makeIncrementalFunctions() const {
        return Functions{
            &AreaComputer::preProcessingImpl,
            &AreaComputer::mergeProcessingImpl,
            &AreaComputer::postProcessingImpl
        };
    }

    explicit AreaComputer(ComponentTree<CNPsType, GraphT>* tree)
        : Base(tree, AREA, makeIncrementalFunctions() )  {
        this->registerDerived(this);
    }

    void preProcessingImpl(NodeId idx, std::span<float> buf) {
        buf[idx] = static_cast<float>(this->tree->getNumCNPsById(idx));
    }

    void mergeProcessingImpl(NodeId parentId, NodeId childId, std::span<float> buf) {
        buf[parentId] += buf[childId];
    }

    void postProcessingImpl(NodeId, std::span<float>) {
        /* no-op */
    }

};
using AreaComputerP = AreaComputer<Pixels>;
using AreaComputerFZ = AreaComputer<FlatZones>;
template <typename GraphT = DefaultFlatZonesGraph>
using AreaComputerFZT = AreaComputer<FlatZones, GraphT>;





/**
 * @brief Computa atributos derivadas do boundind box
 */
template<typename CNPsType, typename GraphT = DefaultFlatZonesGraph>
class BoundingBoxComputer : public AttributeComputer<
    CNPsType,
    GraphT,
    IncrementalFunctions<
        void (BoundingBoxComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>),
        void (BoundingBoxComputer<CNPsType, GraphT>::*)(NodeId, NodeId, std::span<float>),
        void (BoundingBoxComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>)>,
    BoundingBoxComputer<CNPsType, GraphT>> {

private:
    std::vector<int> xmin;
    std::vector<int> xmax;
    std::vector<int> ymin;
    std::vector<int> ymax;
    int numCols;
    int numRows;


public:
    using Functions = IncrementalFunctions<
        void (BoundingBoxComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>),
        void (BoundingBoxComputer<CNPsType, GraphT>::*)(NodeId, NodeId, std::span<float>),
        void (BoundingBoxComputer<CNPsType, GraphT>::*)(NodeId, std::span<float>)>;
    using Base = AttributeComputer<CNPsType, GraphT, Functions, BoundingBoxComputer<CNPsType, GraphT>>;

    Functions makeIncrementalFunctions() const {
        return Functions{
            &BoundingBoxComputer::preProcessingImpl,
            &BoundingBoxComputer::mergeProcessingImpl,
            &BoundingBoxComputer::postProcessingImpl
        };
    }



    explicit BoundingBoxComputer(ComponentTree<CNPsType, GraphT>* tree, Attribute attribute) 
        : Base(tree, attribute, makeIncrementalFunctions()), 
        xmin(tree->getNumNodes(), tree->getNumColsOfImage()), 
        xmax(tree->getNumNodes(), 0), 
        ymin(tree->getNumNodes(), tree->getNumRowsOfImage()), 
        ymax(tree->getNumNodes(), 0),
        numCols(tree->getNumColsOfImage()),
        numRows(tree->getNumRowsOfImage()) {
        
            this->registerDerived(this);
    }

    void preProcessingImpl(NodeId idx, std::span<float>) {
        xmin[idx] = numCols;
        xmax[idx] = 0;
        ymin[idx] = numRows;
        ymax[idx] = 0;

        for (int p : this->tree->getCNPsById(idx)) {
            auto [y, x] = ImageUtils::to2D(p, numCols);
            xmin[idx] = std::min(xmin[idx], x);
            xmax[idx] = std::max(xmax[idx], x);
            ymin[idx] = std::min(ymin[idx], y);
            ymax[idx] = std::max(ymax[idx], y);
        }
    }

    void mergeProcessingImpl(NodeId pid, NodeId cid, std::span<float>) {
        xmin[pid] = std::min(xmin[pid], xmin[cid]);
        xmax[pid] = std::max(xmax[pid], xmax[cid]);
        ymin[pid] = std::min(ymin[pid], ymin[cid]);
        ymax[pid] = std::max(ymax[pid], ymax[cid]);
    }

    void postProcessingImpl(NodeId idx, std::span<float> buffer) {
        if(this->attribute == BOX_WIDTH)
            buffer[idx]  = xmax[idx] - xmin[idx] + 1;
        else if(this->attribute == BOX_HEIGHT)
            buffer[idx] = ymax[idx] - ymin[idx] + 1;
        else if(this->attribute == DIAGONAL_LENGTH) {
            float width  = xmax[idx] - xmin[idx] + 1;
            float height = ymax[idx] - ymin[idx] + 1;
            buffer[idx] = std::sqrt(width*width + height*height);
        }else{
            std::cout << "Attribute not supported!";
            exit(1);
        }
    }

};
using BoundingBoxComputerP = BoundingBoxComputer<Pixels>;
using BoundingBoxComputerFZ = BoundingBoxComputer<FlatZones>;
template <typename GraphT = DefaultFlatZonesGraph>
using BoundingBoxComputerFZT = BoundingBoxComputer<FlatZones, GraphT>;
