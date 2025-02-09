#include <list>
#include <vector>
#include <stack>
#include <queue>
#include <iterator>
#include <utility>

#ifndef NODECT_H
#define NODECT_H

class NodeCT {
private:
	int index=-1; 
    int threshold2; //for maxtree: maximal threshold, same that "level"
	int threshold1;  //for maxtree: minimal threshold
	long int areaCC;
	
	NodeCT* parent;
	std::list<int> cnps; //pixels of the proper part
    std::list<NodeCT*> children;
    
public:
	
    NodeCT();
    NodeCT(int index, NodeCT* parent, int threshold1, int threshold2);
	~NodeCT() {
        parent = nullptr;
    }
    void addCNPs(int p);
	void setCNPs(std::list<int> cnps);
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
	void setNumDescendants(int num);
	NodeCT* getParent();
	void setParent(NodeCT* parent);
	std::list<int>& getCNPs();
	std::list<int> getCNPsCopy();
	void removeCNPs(std::list<int> cnps);
	std::list<NodeCT*>& getChildren();
	int getNumSiblings();

	std::vector<int> getCNPsToVector(){
		std::vector<int> cnpsVector(cnps.begin(), cnps.end());
		return cnpsVector;
	}
	std::vector<NodeCT*> getChildrenToVector(){
		std::vector<NodeCT*> childrenVector(children.begin(), children.end());
		return childrenVector;
	}


   // This function is used by unordered_set to compare elements of Test.
    bool operator==(const NodeCT *node) const;
	bool operator>(const NodeCT *node) const;
	bool operator>=(const NodeCT *node) const;
	bool operator<(const NodeCT *node) const;
	bool operator<=(const NodeCT *node) const;


	class NodeHashFunction {
		public:
			// id is returned as hash function
			size_t operator()(NodeCT* node) const{
				return node->index;
			}
	};

///////////////////////////////////////////////////
	
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
				this->iter = this->currentNode->cnps.begin();
				for (NodeCT *child: this->currentNode->getChildren()){
					s.push(child);
				}	
			}
			InternalIteratorPixelsOfCC& operator++() { 
			    this->iter++; 
				if(this->iter == this->currentNode->cnps.end()){
					if(!s.empty()){
            			this->currentNode = s.top(); s.pop();
						this->iter = this->currentNode->cnps.begin();
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

};

#endif