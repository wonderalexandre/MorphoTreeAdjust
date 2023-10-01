#include "../include/NodeCT.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <list>

NodeCT::NodeCT(){}

NodeCT::NodeCT(int index,  NodeCT* parent, int level) {
		this->index = index;
        this->parent = parent;
        this->level = level;
}

void NodeCT::addCNPs(int p) {
    this->cnps.push_back(p);
}

void NodeCT::addChild(NodeCT* child) {
	this->children.push_back(child);
} 

bool NodeCT::isChild(NodeCT* child){
    auto it = std::find(this->children.begin(), this->children.end(), child);
    return it != children.end();
}

int NodeCT::getIndex(){ return this->index; }

int NodeCT::getLevel(){ return this->level; }

NodeCT* NodeCT::getParent(){  return this->parent; }

void NodeCT::setParent(NodeCT* parent){ this->parent = parent; }

std::list<int>& NodeCT::getCNPs(){  return this->cnps; }

std::list<NodeCT*>& NodeCT::getChildren(){  return this->children; }

int NodeCT::getNumSiblings() {
    if(this->parent != nullptr)
		return this->parent->getChildren().size();
	else
		return 0;
}

bool NodeCT::operator==(const NodeCT* node) const{
    return (this->index == node->index);
}
bool NodeCT::operator>(const NodeCT *node) const{
    return (this->index > node->index);
}
bool NodeCT::operator>=(const NodeCT *node) const{
    return (this->index >= node->index);
}
bool NodeCT::operator<(const NodeCT *node) const{
    return (this->index < node->index);
}
bool NodeCT::operator<=(const NodeCT *node) const{
    return (this->index <= node->index);
}
