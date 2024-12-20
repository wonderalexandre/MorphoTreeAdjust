#include <list>
#include <stack>
#include <iterator>
#include <utility>

#ifndef NODECT_H
#define NODECT_H

class NodeCT {
private:
	int index; 
    int level;
	long int areaCC;
	
	NodeCT* parent;
	std::list<int> cnps;
    std::list<NodeCT*> children;
    
public:
	
    NodeCT();
    NodeCT(int index, NodeCT* parent, int level);
    void addCNPs(int p);
	void setCNPs(std::list<int> cnps);
	void setArea(long int area);
	long int getArea();
    void addChild(NodeCT* child);
	int getIndex();
	int getLevel();
	bool isChild(NodeCT* node);
	void setNumDescendants(int num);
	NodeCT* getParent();
	void setParent(NodeCT* parent);
	std::list<int>& getCNPs();
	std::list<int> getCNPsCopy();
	std::list<NodeCT*>& getChildren();
	int getNumSiblings();

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
	class InternalIteratorNodesOfPathToRoot{
		private:
			NodeCT *currentNode;
			int index;
			using iterator_category = std::input_iterator_tag;
            using value_type = NodeCT; 
		public:
			InternalIteratorNodesOfPathToRoot(NodeCT *obj, int index)  {
				this->currentNode = obj;
				this->index = index;	
			}
			InternalIteratorNodesOfPathToRoot& operator++() { 
				this->index = this->currentNode->index;
				if(this->currentNode != nullptr){
					this->currentNode = this->currentNode->parent;
				}
				return *this; 
			}
			bool operator==(InternalIteratorNodesOfPathToRoot other) const { 
                return this->index == other.index; 
            }
            bool operator!=(InternalIteratorNodesOfPathToRoot other) const { 
                return !(*this == other);
            }
            NodeCT* operator*()  { 
                return (this->currentNode); 
            }  
	};
	class IteratorNodesOfPathToRoot{
		private:
			NodeCT *instance;
		public:
			IteratorNodesOfPathToRoot(NodeCT *obj): instance(obj){}
			InternalIteratorNodesOfPathToRoot begin(){ return InternalIteratorNodesOfPathToRoot(instance, instance->index); }
            InternalIteratorNodesOfPathToRoot end(){ return InternalIteratorNodesOfPathToRoot(instance, 0); }
	};
	IteratorNodesOfPathToRoot getNodesOfPathToRoot(){
	    IteratorNodesOfPathToRoot iter(this);
    	return iter;
	}

	
};

#endif