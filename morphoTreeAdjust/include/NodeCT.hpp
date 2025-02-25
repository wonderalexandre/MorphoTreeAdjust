#include <list>
#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <utility>
#include <functional> 
#include "../include/Common.hpp"

#ifndef NODECT_H
#define NODECT_H

template <typename CNPsType>
class ComponentTree;  //Forward declaration


template <typename CNPsType>
class NodeCT {
private:
	int index=-1; 
    int threshold2; //for maxtree: maximal threshold, same that "level"
	int threshold1;  //for maxtree: minimal threshold
	long int areaCC;
	
	NodeCT* parent;
	CNPsType cnps; //pixels of the proper part 
    std::list<NodeCT*> children;


public:
	
    NodeCT();
    NodeCT(int index, NodeCT* parent, int threshold1, int threshold2);
	~NodeCT() {
        parent = nullptr;
    }
    
    
    ///Métodos disponíveis SOMENTE para `FlatZones`
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    CNPsType moveCNPsByFlatZone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    CNPsType& getCNPsByFlatZone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    int getNumFlatzone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzone(typename CNPsType::value_type&& flatZone, ComponentTree<CNPsType>* tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzones(CNPsType&& flatZones, ComponentTree<CNPsType>* tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsToConnectedFlatzone(typename CNPsType::value_type&& flatZone, ComponentTree<CNPsType>* tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void removeFlatzone(typename CNPsType::value_type& flatzone);

    ///Métodos disponíveis SOMENTE para `Pixels`
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, Pixels>::value, int> = 0>
    void addCNPs(int p);

    void setArea(long int area);
	long int getArea() const;
    void addChild(NodeCT<CNPsType>* child);
    int getIndex() const;
	int getThreshold1() const;
	int getThreshold2() const;
    int getNumCNPs() const;
	int getLevel() const;
	void setLevel(int level);
	bool isChild(NodeCT<CNPsType>* node) const;
    bool isLeaf() const;
	void setNumDescendants(int num);
	NodeCT<CNPsType>* getParent();
	void setParent(NodeCT<CNPsType>* parent);
	std::list<NodeCT<CNPsType>*>& getChildren();
	int getNumSiblings() const;
    int getRepresentativeCNPs() const;



//============= Iterator para iterar os nodes do caminho até o root==============//
	class InternalIteratorNodesOfPathToRoot {
    private:
        NodeCT* currentNode;
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCT*;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCT*;
        using reference = NodeCT*; 

        InternalIteratorNodesOfPathToRoot(NodeCT* obj) : currentNode(obj) {}

        InternalIteratorNodesOfPathToRoot& operator++() {
            if (currentNode) {
                currentNode = currentNode->getParent();
            }
            return *this;
        }

        bool operator==(const InternalIteratorNodesOfPathToRoot& other) const {
            return currentNode == other.currentNode;
        }

        bool operator!=(const InternalIteratorNodesOfPathToRoot& other) const {
            return !(*this == other);
        }

        reference operator*() {  
            return currentNode;  // Retorna ponteiro para o nó atual
        }
    };

	class IteratorNodesOfPathToRoot {
    private:
        NodeCT* instance;
    public:
        explicit IteratorNodesOfPathToRoot(NodeCT* obj) : instance(obj) {}

        InternalIteratorNodesOfPathToRoot begin() const { return InternalIteratorNodesOfPathToRoot(instance); }
        InternalIteratorNodesOfPathToRoot end() const { return InternalIteratorNodesOfPathToRoot(nullptr); }
    };

    IteratorNodesOfPathToRoot getNodesOfPathToRoot() { return IteratorNodesOfPathToRoot(this); }



//============= Iterator para iterar os nodes de um percuso em pos-ordem ==============//
	class InternalIteratorPostOrderTraversal {
    private:
        std::stack<NodeCT*> nodeStack;
        std::stack<NodeCT*> outputStack;
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCT*;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCT*;
        using reference = NodeCT*; // Retorna ponteiro!

        InternalIteratorPostOrderTraversal(NodeCT* root) {
            if (root) {
                nodeStack.push(root);
                while (!nodeStack.empty()) {
                    NodeCT* current = nodeStack.top();nodeStack.pop();
                    outputStack.push(current);
                    for (NodeCT* child : current->getChildren()) {
                        nodeStack.push(child);
                    }
                }
            }
        }

        InternalIteratorPostOrderTraversal& operator++() {
            if (!outputStack.empty()) {
                outputStack.pop();
            }
            return *this;
        }

        reference operator*() {  
            return outputStack.top();  // Retorna ponteiro!
        }

        bool operator==(const InternalIteratorPostOrderTraversal& other) const {
            return (outputStack.empty() == other.outputStack.empty());
        }

        bool operator!=(const InternalIteratorPostOrderTraversal& other) const {
            return !(*this == other);
        }
    };

	class IteratorPostOrderTraversal {
    private:
        NodeCT* root;
    public:
        explicit IteratorPostOrderTraversal(NodeCT* root) : root(root) {}

        InternalIteratorPostOrderTraversal begin() { return InternalIteratorPostOrderTraversal(root); }
        InternalIteratorPostOrderTraversal end() { return InternalIteratorPostOrderTraversal(nullptr); }
    };

    IteratorPostOrderTraversal getIteratorPostOrderTraversal() { return IteratorPostOrderTraversal(this); }



//============= Iterator para iterar os nodes de um percuso em largura ==============//
    class InternalIteratorBreadthFirstTraversal {
    private:
        std::queue<NodeCT*> nodeQueue;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCT*;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCT*;
        using reference = NodeCT*; // Retorna ponteiro!

        InternalIteratorBreadthFirstTraversal(NodeCT* root) {
            if (root) {
                nodeQueue.push(root);
            }
        }

        InternalIteratorBreadthFirstTraversal& operator++() {
            if (!nodeQueue.empty()) {
                NodeCT* current = nodeQueue.front();
                nodeQueue.pop();
                for (NodeCT* child : current->getChildren()) {
                    nodeQueue.push(child);
                }
            }
            return *this;
        }

        reference operator*() {
            return nodeQueue.front();
        }

        bool operator==(const InternalIteratorBreadthFirstTraversal& other) const {
            return nodeQueue.empty() == other.nodeQueue.empty();
        }

        bool operator!=(const InternalIteratorBreadthFirstTraversal& other) const {
            return !(*this == other);
        }
    };

    class IteratorBreadthFirstTraversal {
    private:
        NodeCT* root;

    public:
        explicit IteratorBreadthFirstTraversal(NodeCT* root) : root(root) {}

        InternalIteratorBreadthFirstTraversal begin() { return InternalIteratorBreadthFirstTraversal(root); }
        InternalIteratorBreadthFirstTraversal end() { return InternalIteratorBreadthFirstTraversal(nullptr); }
    };

    // Método para expor o iterador na classe NodeCT
    IteratorBreadthFirstTraversal getIteratorBreadthFirstTraversal() { 
        return IteratorBreadthFirstTraversal(this); 
    }


//////////////////    

//============= Iterator para iterar os pixels de um CC==============//
class InternalIteratorPixelsOfCC {
    private:
        NodeCT<CNPsType>* currentNode;
        std::stack<NodeCT<CNPsType>*> s;
        typename CNPsType::iterator iter;  // Iterador principal (pode ser Pixels ou FlatZones)
        typename CNPsType::iterator endIter;  // Final da lista principal
        typename CNPsType::value_type::iterator subIter;  //  Iterador secundário (caso seja FlatZones)
        bool isFlatZones;  // Indica se estamos lidando com `FlatZones`
        int countArea;

        using iterator_category = std::input_iterator_tag;
        using value_type = int;

    public:
        InternalIteratorPixelsOfCC(NodeCT<CNPsType>* obj, int area) : currentNode(obj), countArea(area) {
            isFlatZones = std::is_same<CNPsType, std::list<std::list<int>>>::value;

            if (!this->currentNode->cnps.empty()) {
                iter = this->currentNode->cnps.begin();
                endIter = this->currentNode->cnps.end();

                if (isFlatZones && iter != endIter) {
                    subIter = iter->begin();
                }
            }

            for (NodeCT<CNPsType>* child : this->currentNode->getChildren()) {
                s.push(child);
            }
        }

        InternalIteratorPixelsOfCC& operator++() {
            if (isFlatZones) {
                // Iterando dentro de uma flatzone
                ++subIter;
                while (iter != endIter && subIter == iter->end()) {
                    ++iter;
                    if (iter != endIter) {
                        subIter = iter->begin();
                    }
                }
            } else {
                // Iterando em uma lista simples de pixels
                ++iter;
            }

            // Se acabou, ir para próximo nó na árvore
            if (iter == endIter && !s.empty()) {
                this->currentNode = s.top();
                s.pop();

                if (!this->currentNode->cnps.empty()) {
                    iter = this->currentNode->cnps.begin();
                    endIter = this->currentNode->cnps.end();

                    if (isFlatZones && iter != endIter) {
                        subIter = iter->begin();
                    }
                }

                for (NodeCT<CNPsType>* child : this->currentNode->getChildren()) {
                    s.push(child);
                }
            }

            this->countArea++;
            return *this;
        }

        bool operator==(const InternalIteratorPixelsOfCC& other) const {
            return this->countArea == other.countArea;
        }

        bool operator!=(const InternalIteratorPixelsOfCC& other) const {
            return !(*this == other);
        }

        int operator*() const {
            if (isFlatZones) {
                return *subIter;  
            } else {
                return *(iter->begin());
            }
        }
    };

    class IteratorPixelsOfCC { //Classe Externa para o Iterador
    private:
        NodeCT<CNPsType>* instance;
        int area;

    public:
        IteratorPixelsOfCC(NodeCT<CNPsType>* obj, int _area) : instance(obj), area(_area) {}

        InternalIteratorPixelsOfCC begin() { return InternalIteratorPixelsOfCC(instance, 0); }
        InternalIteratorPixelsOfCC end() { return InternalIteratorPixelsOfCC(instance, area); }
        int front() const{ return instance->getRepresentativeCNPs(); }
    };
    IteratorPixelsOfCC getPixelsOfCC() {
        return IteratorPixelsOfCC(this, this->areaCC);
    }



//============= Iterator para iterar os pixels compactos (CNPs) ==============//

class InternalIteratorCNPs {
    private:
        using FlatzoneIterator = std::list<std::list<int>>::iterator;
        using CNPsIterator = std::list<int>::iterator;
        
        FlatzoneIterator flatzoneIt, flatzoneEnd;
        CNPsIterator cnpsIt;
    
        void advance() {
            if (flatzoneIt != flatzoneEnd && cnpsIt == flatzoneIt->end()) {
                ++flatzoneIt;
                if (flatzoneIt != flatzoneEnd) {
                    cnpsIt = flatzoneIt->begin();
                }
            }
        }
       
        
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = int;
        using difference_type = std::ptrdiff_t;
        using pointer = int*;
        using reference = int&;
        
        InternalIteratorCNPs(FlatzoneIterator flatzoneIt, FlatzoneIterator flatzoneEnd)
            : flatzoneIt(flatzoneIt), flatzoneEnd(flatzoneEnd) {
            if (flatzoneIt != flatzoneEnd) {
                cnpsIt = flatzoneIt->begin();
                advance();
            }
        }
        
        reference operator*() {
            assert(flatzoneIt != flatzoneEnd && "Tentando acessar um iterador inválido");
            return *cnpsIt;
        }
        
        InternalIteratorCNPs& operator++() {
            ++cnpsIt;
            advance();
            return *this;
        }
        
        bool operator==(const InternalIteratorCNPs& other) const {
            return flatzoneIt == other.flatzoneIt && (flatzoneIt == flatzoneEnd || cnpsIt == other.cnpsIt);
        }
        
        bool operator!=(const InternalIteratorCNPs& other) const {
            return !(*this == other);
        }
    };
    
    class IteratorCNPs {
        private:
            std::list<std::list<int>>* cnpsByFlatzone;

        public:
            explicit IteratorCNPs(NodeCT* node) : cnpsByFlatzone(&node->cnps) {}
            explicit IteratorCNPs(std::list<std::list<int>>* cnpsByFlatzone) : cnpsByFlatzone(cnpsByFlatzone) {}
        
            InternalIteratorCNPs begin() { 
                return InternalIteratorCNPs(cnpsByFlatzone->begin(), cnpsByFlatzone->end()); 
            }
        
            InternalIteratorCNPs end() { 
                return InternalIteratorCNPs(cnpsByFlatzone->end(), cnpsByFlatzone->end()); 
            }
        
            int front() const { return cnpsByFlatzone->front().front(); }
            int back() const { return cnpsByFlatzone->back().back(); }
    };
    
    
    // Método getCNPs() quando CNPsType == FlatZones
    template<typename T = CNPsType>
    std::enable_if_t<std::is_same<T, FlatZones>::value, IteratorCNPs> getCNPs() {
        return IteratorCNPs(this);
    }

    // Método getCNPs() quando CNPsType == Pixels
    template<typename T = CNPsType>
    std::enable_if_t<std::is_same<T, Pixels>::value, std::list<int>&> getCNPs() {
        return cnps;
    }







    
};

#include "NodeCT.tpp"


#endif
