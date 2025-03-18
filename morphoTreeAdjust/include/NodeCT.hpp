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
    FlatZone& getFlatZone(int idFlatZone);

    //template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    int getNumFlatzone();

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzone(FlatZone&& flatZone, ComponentTree<CNPsType>* tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsOfDisjointFlatzones(CNPsType&& flatZones, ComponentTree<CNPsType>* tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void addCNPsToConnectedFlatzone(FlatZone&& flatZone, ComponentTree<CNPsType>* tree);

    template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    void removeFlatzone(int idFlatZone);

    //template<typename T = CNPsType, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int> = 0>
    //int getFlatZoneID(int pixel);

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
	NodeCT<CNPsType>* getParent();
	void setParent(NodeCT<CNPsType>* parent);
	std::list<NodeCT<CNPsType>*>& getChildren();
	int getNumSiblings() const;
    int getRepresentativeCNPs() const;
    int computerNumDescendants();


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

/////////
// **Iterador para coletar ramos em pós-ordem**
// Classe do Iterador Pós-Ordem por Ramos
class InternalIteratorBranchPostOrderTraversal {
    private:
        std::stack<NodeCT*> processingStack;
        std::stack<NodeCT*> postOrderStack;
        std::list<std::list<NodeCT*>> branches; // Lista de ramos
        //std::list<std::list<NodeCT*>>::iterator branchIterator;
        typename std::list<typename std::list<NodeCT<CNPsType>*>>::iterator branchIterator;
    
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::list<NodeCT*>;
        using pointer = std::list<NodeCT*>*;
        using reference = std::list<NodeCT*>&;
    
        // Construtor: Faz percurso pós-ordem e coleta ramos
        InternalIteratorBranchPostOrderTraversal(NodeCT* root) {
            if (!root) return;
    
            // Passo 1: Gerar percurso pós-ordem
            std::stack<NodeCT*> tempStack;
            tempStack.push(root);
    
            while (!tempStack.empty()) {
                NodeCT* current = tempStack.top();
                tempStack.pop();
                postOrderStack.push(current);
    
                for (NodeCT* child : current->getChildren()) {
                    tempStack.push(child);
                }
            }
    
            // Passo 2: Construir ramos baseando-se no percurso pós-ordem
            std::list<NodeCT*> currentBranch;
            while (!postOrderStack.empty()) {
                NodeCT* node = postOrderStack.top();
                postOrderStack.pop();
    
                if (!currentBranch.empty()) {
                    NodeCT* lastNode = currentBranch.back();
                    if (lastNode->parent && lastNode->parent->getChildren().back() != lastNode) {
                        branches.push_back(currentBranch);
                        currentBranch.clear();
                    }
                }
                currentBranch.push_back(node);
            }
            if (!currentBranch.empty()) {
                branches.push_back(currentBranch);
            }
    
            // Iniciar o iterador
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
    
    // Classe Externa para o Iterador
    class IteratorBranchPostOrderTraversal {
    private:
        NodeCT* root;
    public:
        explicit IteratorBranchPostOrderTraversal(NodeCT* root) : root(root) {}
    
        InternalIteratorBranchPostOrderTraversal begin() { return InternalIteratorBranchPostOrderTraversal(root); }
        InternalIteratorBranchPostOrderTraversal end() { return InternalIteratorBranchPostOrderTraversal(nullptr); }
    };
    
    // Método auxiliar para obter o iterador
    IteratorBranchPostOrderTraversal getIteratorBranchPostOrderTraversal() {
        return IteratorBranchPostOrderTraversal(this);
    }


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


//============= Iterator para iterar os pixels compactos (CNPs) ==============//
    class InternalIteratorCNPs {
    private:
        using FlatzoneIterator = std::unordered_map<int, std::list<int>>::iterator;
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
        
        InternalIteratorCNPs(FlatzoneIterator flatzoneBegin, FlatzoneIterator flatzoneEnd)
            : flatzoneIt(flatzoneBegin), flatzoneEnd(flatzoneEnd) {
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
        std::unordered_map<int, std::list<int>>* cnpsByFlatzone;

    public:
        explicit IteratorCNPs(NodeCT* node) : cnpsByFlatzone(&node->cnps) {}
        explicit IteratorCNPs(std::unordered_map<int, std::list<int>>* cnpsByFlatzone) 
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
        return IteratorCNPs(this);
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
        NodeCT<CNPsType>* currentNode;
        std::stack<NodeCT<CNPsType>*> s;

        std::list<int>::iterator iter;
        int countArea;
        using iterator_category = std::input_iterator_tag;
        using value_type = int; 
    public:
        InternalIteratorPixelsOfCC_Pixels(NodeCT<CNPsType>* obj, int area)  {
            this->currentNode = obj;
            this->countArea =area;
            this->iter = this->currentNode->cnps.begin();
            for (NodeCT<CNPsType>* child: this->currentNode->getChildren()){
                s.push(child);
            }	
        }
        InternalIteratorPixelsOfCC_Pixels& operator++() { 
            this->iter++; 
            if(this->iter == this->currentNode->cnps.end()){
                if(!s.empty()){
                    this->currentNode = s.top(); s.pop();
                    this->iter = this->currentNode->cnps.begin();
                    for (NodeCT<CNPsType>* child: currentNode->getChildren()){
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
        NodeCT<CNPsType>* currentNode;
        std::stack<NodeCT<CNPsType>*> s;
        typename std::unordered_map<int, std::list<int>>::iterator iterMap;  
        typename std::unordered_map<int, std::list<int>>::iterator endIterMap;  
        typename std::list<int>::iterator iterList;
        typename std::list<int>::iterator endIterList;
        int countArea;

        using iterator_category = std::input_iterator_tag;
        using value_type = int;

    public:
     InternalIteratorPixelsOfCC_FlatZones(NodeCT<CNPsType>* obj, int area) : currentNode(obj), countArea(area) {
            //if (area > 0) {
                iterMap = this->currentNode->cnps.begin();
                endIterMap = this->currentNode->cnps.end();

                iterList = iterMap->second.begin(); 
                endIterList = iterMap->second.end(); 
                
                for (NodeCT<CNPsType>* child : this->currentNode->getChildren()) {
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
                
                for (NodeCT<CNPsType>* child : this->currentNode->getChildren()) {
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
        NodeCT<CNPsType>* instance;
        int area;
    public:
        explicit IteratorPixelsOfCC(NodeCT<CNPsType>* obj, int _area) : instance(obj), area(_area) {}

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
        return IteratorPixelsOfCC(this, this->areaCC);
    }





    
};


#include "../include/NodeCT.tpp"


#endif
