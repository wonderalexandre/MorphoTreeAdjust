#include <list>
#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <utility>
#include <functional> 
#include <optional>
#include <stdexcept>
#include "../include/Common.hpp"

#ifndef NODECT_H
#define NODECT_H

template <typename CNPsType>
class ComponentTree;  //Forward declaration


template <typename CNPsType>
class NodeCT : public std::enable_shared_from_this<NodeCT<CNPsType>> {
    private:
	int index=-1; 
    int threshold2; //for maxtree: maximal threshold, same that "level"
	int threshold1;  //for maxtree: minimal threshold
	long int areaCC;
    //long int numCNPs = -1;
	
	NodeCTPtr<CNPsType> parent;
	std::list<NodeCTPtr<CNPsType>> children;
    CNPsType cnps; //pixels of the proper part 

public:
	
    NodeCT();
    NodeCT(int index, NodeCTPtr<CNPsType> parent, int threshold1, int threshold2);
	~NodeCT() {
        parent = nullptr;
    }
    
    
    ///Métodos disponíveis SOMENTE para `FlatZones`
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    CNPsType moveCNPsByFlatZone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    CNPsType& getCNPsByFlatZone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    FlatZone& getFlatZone(int idFlatZone);

    //template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    int getNumFlatzone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzone(FlatZone&& flatZone, ComponentTreeFZPtr tree = nullptr, int capacity = -1);
    

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzones(CNPsType&& flatZones, ComponentTreeFZPtr tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsToConnectedFlatzone(FlatZone&& flatZone, ComponentTreeFZPtr tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void removeFlatzone(int idFlatZone);

    ///Métodos disponíveis SOMENTE para `Pixels`
    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, Pixels>::value, int> = 0>
    void addCNPs(int p);

    void setArea(long int area);
	long int getArea() const;
    void addChild(NodeCTPtr<CNPsType> child);
    int getIndex() const;
	int getThreshold1() const;
	int getThreshold2() const;
    int getNumCNPs() ;
	int getLevel() const;
	void setLevel(int level);
	bool isChild(NodeCTPtr<CNPsType> node) const;
    bool isLeaf() const;
	NodeCTPtr<CNPsType> getParent();
	void setParent(NodeCTPtr<CNPsType> parent);
	std::list<NodeCTPtr<CNPsType>>& getChildren();
	int getNumSiblings() const;
    int computerNumDescendants();
    


//============= Iterator para iterar os nodes do caminho até o root==============//
class InternalIteratorNodesOfPathToRoot {
    private:
        NodeCTPtr<CNPsType> currentNode;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCTPtr<CNPsType>;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCTPtr<CNPsType>*;
        using reference = NodeCTPtr<CNPsType>&;
    
        InternalIteratorNodesOfPathToRoot(NodeCTPtr<CNPsType> obj) : currentNode(obj) {}
    
        InternalIteratorNodesOfPathToRoot& operator++() {
            if (currentNode) {
                currentNode = currentNode->getParent();  // Retorna outro shared_ptr
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
            return currentNode;
        }
    };
    
    class IteratorNodesOfPathToRoot {
    private:
        NodeCTPtr<CNPsType> instance;
    
    public:
        explicit IteratorNodesOfPathToRoot(NodeCTPtr<CNPsType> obj) : instance(obj) {}
    
        InternalIteratorNodesOfPathToRoot begin() const { return InternalIteratorNodesOfPathToRoot(instance); }
        InternalIteratorNodesOfPathToRoot end() const { return InternalIteratorNodesOfPathToRoot(nullptr); }
    };
    
    // Chamador usa this como shared_ptr:
    IteratorNodesOfPathToRoot getNodesOfPathToRoot() {
        return IteratorNodesOfPathToRoot(this->shared_from_this());
    }
    

/////////
// **Iterador para coletar ramos em pós-ordem**
// Classe do Iterador Pós-Ordem por Ramos
class InternalIteratorBranchPostOrderTraversal {
    private:
        std::stack<NodeCTPtr<CNPsType>> processingStack;
        std::stack<NodeCTPtr<CNPsType>> postOrderStack;
        std::list<std::list<NodeCTPtr<CNPsType>>> branches;
        typename std::list<std::list<NodeCTPtr<CNPsType>>>::iterator branchIterator;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::list<NodeCTPtr<CNPsType>>;
        using pointer = std::list<NodeCTPtr<CNPsType>>*;
        using reference = std::list<NodeCTPtr<CNPsType>>&;
    
        InternalIteratorBranchPostOrderTraversal(NodeCTPtr<CNPsType> root) {
            if (!root) return;
    
            std::stack<NodeCTPtr<CNPsType>> tempStack;
            tempStack.push(root);
    
            while (!tempStack.empty()) {
                auto current = tempStack.top();
                tempStack.pop();
                postOrderStack.push(current);
    
                for (auto& child : current->getChildren()) {
                    tempStack.push(child);
                }
            }
    
            std::list<NodeCTPtr<CNPsType>> currentBranch;
            while (!postOrderStack.empty()) {
                auto node = postOrderStack.top();
                postOrderStack.pop();
    
                if (!currentBranch.empty()) {
                    auto lastNode = currentBranch.back();
                    if (lastNode->getParent() && lastNode->getParent()->getChildren().back() != lastNode) {
                        branches.push_back(currentBranch);
                        currentBranch.clear();
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
        NodeCTPtr<CNPsType> root;
    public:
        explicit IteratorBranchPostOrderTraversal(NodeCTPtr<CNPsType> root) : root(root) {}
    
        InternalIteratorBranchPostOrderTraversal begin() { return InternalIteratorBranchPostOrderTraversal(root); }
        InternalIteratorBranchPostOrderTraversal end() { return InternalIteratorBranchPostOrderTraversal(nullptr); }
    };
    
    // Método da classe
    IteratorBranchPostOrderTraversal getIteratorBranchPostOrderTraversal() {
        return IteratorBranchPostOrderTraversal(this->shared_from_this());
    }
    


//============= Iterator para iterar os nodes de um percuso em pos-ordem ==============//
	class InternalIteratorPostOrderTraversal {
    private:
        std::stack<NodeCTPtr<CNPsType>> nodeStack;
        std::stack<NodeCTPtr<CNPsType>> outputStack;
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCTPtr<CNPsType>;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCTPtr<CNPsType>*;
        using reference = NodeCTPtr<CNPsType>&; 

        InternalIteratorPostOrderTraversal(NodeCTPtr<CNPsType> root) {
            if (root) {
                nodeStack.push(root);
                while (!nodeStack.empty()) {
                    NodeCTPtr<CNPsType> current = nodeStack.top();nodeStack.pop();
                    outputStack.push(current);
                    for (NodeCTPtr<CNPsType> child : current->getChildren()) {
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
        NodeCTPtr<CNPsType> root;
    public:
        explicit IteratorPostOrderTraversal(NodeCTPtr<CNPsType> root) : root(root) {}

        InternalIteratorPostOrderTraversal begin() { return InternalIteratorPostOrderTraversal(root); }
        InternalIteratorPostOrderTraversal end() { return InternalIteratorPostOrderTraversal(nullptr); }
    };

    IteratorPostOrderTraversal getIteratorPostOrderTraversal() { return IteratorPostOrderTraversal(this->shared_from_this()); }



//============= Iterator para iterar os nodes de um percuso em largura ==============//
    class InternalIteratorBreadthFirstTraversal {
    private:
        std::queue<NodeCTPtr<CNPsType>> nodeQueue;

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCTPtr<CNPsType>;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCTPtr<CNPsType>*;
        using reference = NodeCTPtr<CNPsType>&; // Retorna ponteiro!

        InternalIteratorBreadthFirstTraversal(NodeCTPtr<CNPsType> root) {
            if (root) {
                nodeQueue.push(root);
            }
        }

        InternalIteratorBreadthFirstTraversal& operator++() {
            if (!nodeQueue.empty()) {
                NodeCTPtr<CNPsType> current = nodeQueue.front();
                nodeQueue.pop();
                for (NodeCTPtr<CNPsType> child : current->getChildren()) {
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
    NodeCTPtr<CNPsType> root;

    public:
        explicit IteratorBreadthFirstTraversal(NodeCTPtr<CNPsType> root) : root(root) {}

        InternalIteratorBreadthFirstTraversal begin() { return InternalIteratorBreadthFirstTraversal(root); }
        InternalIteratorBreadthFirstTraversal end() { return InternalIteratorBreadthFirstTraversal(nullptr); }
    };

    // Método para expor o iterador na classe NodeCT
    IteratorBreadthFirstTraversal getIteratorBreadthFirstTraversal() { 
        return IteratorBreadthFirstTraversal(this->shared_from_this()); 
    }


//////////////////    


//============= Iterator para iterar os pixels compactos (CNPs) ==============//
    class InternalIteratorCNPs {
    private:
        using FlatzoneIterator = FlatZones::iterator;
        using CNPsIterator = std::list<int>::iterator;
        
        FlatzoneIterator flatzoneIt, flatzoneEnd;
        CNPsIterator cnpsIt;
        
        void advance() {
            while (flatzoneIt != flatzoneEnd && cnpsIt == flatzoneIt->second.end()) {
                ++flatzoneIt;
                if (flatzoneIt != flatzoneEnd) {
                    cnpsIt = flatzoneIt->second.begin();
                }
            }
        }
        
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = int;
        using difference_type = std::ptrdiff_t;
        using pointer = int*;
        using reference = int&;
        
        InternalIteratorCNPs(FlatzoneIterator flatzoneBegin, FlatzoneIterator flatzoneEnd) : flatzoneIt(flatzoneBegin), flatzoneEnd(flatzoneEnd) {
            if (flatzoneIt != flatzoneEnd) {
                cnpsIt = flatzoneIt->second.begin();
                advance();
            }
        }
        
        reference operator*() {
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
        FlatZones* cnpsByFlatzone;

    public:
        explicit IteratorCNPs(NodeCTPtr<CNPsType> node) : cnpsByFlatzone(&node->cnps) {}
        explicit IteratorCNPs(FlatZones* cnpsByFlatzone) 
            : cnpsByFlatzone(cnpsByFlatzone) {}

        InternalIteratorCNPs begin() { 
            return InternalIteratorCNPs(cnpsByFlatzone->begin(), cnpsByFlatzone->end()); 
        }

        InternalIteratorCNPs end() { 
            return InternalIteratorCNPs(cnpsByFlatzone->end(), cnpsByFlatzone->end()); 
        }

        int front() const { return cnpsByFlatzone->begin()->second.front(); }

    };
   
    
    
    // Método getCNPs() quando CNPsType == FlatZones
    template<typename T = CNPsType>
    std::enable_if_t<std::is_same<T, FlatZones>::value, IteratorCNPs> getCNPs() {
        return IteratorCNPs(this->shared_from_this());
    }

    // Método getCNPs() quando CNPsType == Pixels
    template<typename T = CNPsType>
    std::enable_if_t<std::is_same<T, Pixels>::value, std::list<int>&> getCNPs() {
        return cnps;
    }




//============= Iterator para iterar os pixels de um CC==============//

// ================== Iterador para Pixels ==================
template <typename T = CNPsType, typename std::enable_if_t<std::is_same_v<T, Pixels>, int> = 0>
class InternalIteratorPixelsOfCC_Pixels{
    private:
        NodeCTPtr<CNPsType> currentNode;
        std::stack<NodeCTPtr<CNPsType>> s;

        std::list<int>::iterator iter;
        int countArea;
        using iterator_category = std::input_iterator_tag;
        using value_type = int; 
    public:
        InternalIteratorPixelsOfCC_Pixels(NodeCTPtr<CNPsType> obj, int area)  {
            this->currentNode = obj;
            this->countArea =area;
            this->iter = this->currentNode->cnps.begin();
            for (NodeCTPtr<CNPsType> child: this->currentNode->getChildren()){
                s.push(child);
            }	
        }
        InternalIteratorPixelsOfCC_Pixels& operator++() { 
            this->iter++; 
            if(this->iter == this->currentNode->cnps.end()){
                if(!s.empty()){
                    this->currentNode = s.top(); s.pop();
                    this->iter = this->currentNode->cnps.begin();
                    for (NodeCTPtr<CNPsType> child: currentNode->getChildren()){
                        s.push(child);
                    }
                }
            }
            this->countArea++;
            return *this; 
        }
        bool operator==(InternalIteratorPixelsOfCC_Pixels other) const { 
            return this->countArea == other.countArea; 
        }
        bool operator!=(InternalIteratorPixelsOfCC_Pixels other) const { 
            return !(*this == other);
        }
        int operator*() const { 
            return (*this->iter); 
        }  
};
// ================== Iterador para FlatZones ==================
template <typename T = CNPsType, typename std::enable_if_t<std::is_same_v<T, FlatZones>, int> = 0>
class InternalIteratorPixelsOfCC_FlatZones {
    private:
    NodeCTPtr<CNPsType> currentNode;
        std::stack<NodeCTPtr<CNPsType>> s;
        typename std::unordered_map<int, std::list<int>>::iterator iterMap;  
        typename std::unordered_map<int, std::list<int>>::iterator endIterMap;  
        typename std::list<int>::iterator iterList;
        typename std::list<int>::iterator endIterList;
        int countArea;

        using iterator_category = std::input_iterator_tag;
        using value_type = int;

    public:
     InternalIteratorPixelsOfCC_FlatZones(NodeCTPtr<CNPsType> obj, int area) : currentNode(obj), countArea(area) {
            //if (area > 0) {
                iterMap = this->currentNode->cnps.begin();
                endIterMap = this->currentNode->cnps.end();

                iterList = iterMap->second.begin(); 
                endIterList = iterMap->second.end(); 
                
                for (NodeCTPtr<CNPsType> child : this->currentNode->getChildren()) {
                    s.push(child);
                }
           // }

            
        }

        InternalIteratorPixelsOfCC_FlatZones& operator++() {
                ++iterList;  // Ir para o próximo elemento do mapa
                if(iterList == endIterList){
                    ++iterMap;
                    if(iterMap != endIterMap){
                        iterList = iterMap->second.begin(); 
                        endIterList = iterMap->second.end(); 
                    }
                }
            
            // Se acabar os elementos, ir para o próximo nó na árvore
            if (iterMap == endIterMap && !s.empty()) {
                this->currentNode = s.top(); s.pop();

                iterMap = this->currentNode->cnps.begin();
                endIterMap = this->currentNode->cnps.end();
                iterList = iterMap->second.begin(); 
                endIterList = iterMap->second.end(); 
                
                for (NodeCTPtr<CNPsType> child : this->currentNode->getChildren()) {
                    s.push(child);
                }
            }

            this->countArea++;
            return *this;
        }

        bool operator==(const InternalIteratorPixelsOfCC_FlatZones& other) const {
            return this->countArea == other.countArea;
        }

        bool operator!=(const InternalIteratorPixelsOfCC_FlatZones& other) const {
            return !(*this == other);
        }

        int operator*() const {
            return *iterList;  // Retorna um pixel dentro da `std::list<int>`
        }
    };

    class IteratorPixelsOfCC {
    private:
        NodeCTPtr<CNPsType> instance;
        int area;
    public:
        explicit IteratorPixelsOfCC(NodeCTPtr<CNPsType> obj, int _area) : instance(obj), area(_area) {}

        auto begin() {
            if constexpr (std::is_same_v<CNPsType, Pixels>) {
                return InternalIteratorPixelsOfCC_Pixels<CNPsType>(instance, 0);
            } else {
                return InternalIteratorPixelsOfCC_FlatZones<CNPsType>(instance, 0);
            }
        }

        auto end() {
            if constexpr (std::is_same_v<CNPsType, Pixels>) {
                return InternalIteratorPixelsOfCC_Pixels<CNPsType>(instance, area);
            } else {
                return InternalIteratorPixelsOfCC_FlatZones<CNPsType>(instance, area);
            }
        }

        int front() const {
            return instance->getRepresentativeCNPs();
        }
    };
    IteratorPixelsOfCC getPixelsOfCC() {
        return IteratorPixelsOfCC(this->shared_from_this(), this->areaCC);
    }




    
};


#include "../include/NodeCT.tpp"


#endif
