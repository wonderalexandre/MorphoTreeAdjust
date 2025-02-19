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

class ComponentTree;  //Forward declaration

struct ListRefHash {
    size_t operator()(const std::reference_wrapper<std::list<int>>& ref) const {
        return reinterpret_cast<size_t>(&ref.get());  // Usa o endereço como hash
    }
};

struct ListRefEqual {
    bool operator()(const std::reference_wrapper<std::list<int>>& lhs, 
                    const std::reference_wrapper<std::list<int>>& rhs) const {
        return &lhs.get() == &rhs.get();  // Compara pelos endereços
    }
};


class NodeCT {
private:
	int index=-1; 
    int threshold2; //for maxtree: maximal threshold, same that "level"
	int threshold1;  //for maxtree: minimal threshold
	long int areaCC;
	
	NodeCT* parent;
	std::list<std::list<int>> cnpsByFlatzone; //pixels of the proper part by flatzone
    std::list<NodeCT*> children;


public:
	
    NodeCT();
    NodeCT(int index, NodeCT* parent, int threshold1, int threshold2);
	~NodeCT() {
        parent = nullptr;
    }
    
    
    std::list<std::list<int>> moveCNPsByFlatZone(); 
    std::list<std::list<int>> getCopyCNPsByFlatZone();
	
	
    int getNumFlatzone();
    //void setCNPsByFlatZone(std::list<std::list<int>>&& cnpsByFlatzone, ComponentTree* tree);
    
    void addCNPsOfDisjointFlatzone(std::list<int>&& flatZone, ComponentTree* tree);
    void addCNPsOfDisjointFlatzones(std::list<std::list<int>>&& flatZones, ComponentTree* tree);
    void addCNPsToConnectedFlatzone(std::list<int>&& flatZone, ComponentTree* tree);
    int getCNP(size_t index);
    void removeFlatzone(std::list<int>& flatzone);
    

    void setArea(long int area);
	long int getArea() const;
    void addChild(NodeCT* child);
    int getIndex() const;
	int getThreshold1() const;
	int getThreshold2() const;
    int getNumCNPs() const;
	int getLevel() const;
	void setLevel(int level);
	bool isChild(NodeCT* node);
    bool isLeaf();
	void setNumDescendants(int num);
	NodeCT* getParent();
	void setParent(NodeCT* parent);
	
	std::list<NodeCT*>& getChildren();
	int getNumSiblings();

    std::vector<int> getCNPsToVector() {
        std::vector<int> cnpsVector;
        for (const auto& region : cnpsByFlatzone) {
            cnpsVector.insert(cnpsVector.end(), region.begin(), region.end());
        }
        return cnpsVector;
    }

	std::vector<NodeCT*> getChildrenToVector(){
		std::vector<NodeCT*> childrenVector(children.begin(), children.end());
		return childrenVector;
	}


///////////////////////////////////////////////////
	/*
    class InternalIteratorPixelsOfCC{
		private:
			NodeCT *currentNode;
			std::stack<NodeCT*> s;
			std::list<int>::iterator iter;
			int countArea;
			using iterator_category = std::input_iterator_tag;
            using value_type = int; 
		public:
			InternalIteratorPixelsOfCC(NodeCT *obj, int area)  {
				this->currentNode = obj;
				this->countArea =area;
				this->iter = this->currentNode->cnpsByFlatzone.begin();
				for (NodeCT *child: this->currentNode->getChildren()){
					s.push(child);
				}	
			}
			InternalIteratorPixelsOfCC& operator++() { 
			    this->iter++; 
				if(this->iter == this->currentNode->cnpsByFlatzone.end()){
					if(!s.empty()){
            			this->currentNode = s.top(); s.pop();
						this->iter = this->currentNode->cnpsByFlatzone.begin();
						for (NodeCT *child: currentNode->getChildren()){
                		    s.push(child);
						}
					}
				}
				this->countArea++;
				return *this; 
            }
            bool operator==(InternalIteratorPixelsOfCC other) const { 
                return this->countArea == other.countArea; 
            }
            bool operator!=(InternalIteratorPixelsOfCC other) const { 
                return !(*this == other);
            }
            int operator*() const { 
                return (*this->iter); 
            }  
    };
    */
    class InternalIteratorPixelsOfCC {
    private:
        NodeCT* currentNode;
        std::stack<NodeCT*> s;
        std::list<std::list<int>>::iterator flatzoneIter;  // ✅ Iterador do vetor
        std::list<int>::iterator iter;  // ✅ Iterador da lista
        int countArea;

        using iterator_category = std::input_iterator_tag;
        using value_type = int;

    public:
        InternalIteratorPixelsOfCC(NodeCT* obj, int area) {
            this->currentNode = obj;
            this->countArea = area;

            // Verifica se há flatzones antes de acessar
            if (!this->currentNode->cnpsByFlatzone.empty()) {
                this->flatzoneIter = this->currentNode->cnpsByFlatzone.begin();
                this->iter = flatzoneIter->begin();
            }

            for (NodeCT* child : this->currentNode->getChildren()) {
                s.push(child);
            }
        }

        InternalIteratorPixelsOfCC& operator++() {
            if (this->flatzoneIter != this->currentNode->cnpsByFlatzone.end()) {
                ++this->iter;

                // ✅ Se terminamos a lista atual, passamos para a próxima flatzone
                while (this->iter == this->flatzoneIter->end()) {
                    ++this->flatzoneIter;
                    if (this->flatzoneIter == this->currentNode->cnpsByFlatzone.end()) {
                        break;
                    }
                    this->iter = this->flatzoneIter->begin();
                }
            }

            if (this->flatzoneIter == this->currentNode->cnpsByFlatzone.end() && !s.empty()) {
                this->currentNode = s.top();
                s.pop();

                if (!this->currentNode->cnpsByFlatzone.empty()) {
                    this->flatzoneIter = this->currentNode->cnpsByFlatzone.begin();
                    this->iter = this->flatzoneIter->begin();
                }

                for (NodeCT* child : this->currentNode->getChildren()) {
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
            return *this->iter;
        }
    };

	class IteratorPixelsOfCC{
		private:
			NodeCT *instance;
			int area;
		public:
			IteratorPixelsOfCC(NodeCT *obj, int _area): instance(obj), area(_area) {}
			InternalIteratorPixelsOfCC begin(){ return InternalIteratorPixelsOfCC(instance, 0); }
            InternalIteratorPixelsOfCC end(){ return InternalIteratorPixelsOfCC(instance, area); }
	};	
	IteratorPixelsOfCC getPixelsOfCC(){
	    return IteratorPixelsOfCC(this, this->areaCC);
	}



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
                explicit IteratorCNPs(NodeCT* node) : cnpsByFlatzone(&node->cnpsByFlatzone) {}
            
                InternalIteratorCNPs begin() { 
                    return InternalIteratorCNPs(cnpsByFlatzone->begin(), cnpsByFlatzone->end()); 
                }
            
                InternalIteratorCNPs end() { 
                    return InternalIteratorCNPs(cnpsByFlatzone->end(), cnpsByFlatzone->end()); 
                }
            
                int front() const { return cnpsByFlatzone->front().front(); }
                int back() const { return cnpsByFlatzone->back().back(); }
        };
        
        IteratorCNPs getCNPs() { return IteratorCNPs(this); }
        

};

#endif