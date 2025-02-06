#include <list>
#include <vector>
#include <stack>
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
	
	class InternalIteratorNodesOfPathToRoot {
    private:
        NodeCT* currentNode;
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeCT*;
        using difference_type = std::ptrdiff_t;
        using pointer = NodeCT*;
        using reference = NodeCT*; // Retorna ponteiro!

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

};

#endif