#include "../include/NodeCT.hpp"
#include "../include/AdjacencyRelation.hpp"

#include <list>
#include <iostream>
#include <algorithm>

NodeCT::NodeCT(): children(), cnps(){}

NodeCT::NodeCT(int index,  NodeCT* parent, int threshold1, int threshold2) {
		this->index = index;
        this->parent = parent;
        this->threshold2 = threshold2;
        this->threshold1 = threshold1;
        this->areaCC = 0;
}

void NodeCT::addCNPs(int p) {
    this->cnps.push_back(p);
}

void NodeCT::setCNPs(std::list<int> cnps){
    this->cnps = cnps;
}

void NodeCT::addChild(NodeCT* child) {
	this->children.push_back(child);
} 

bool NodeCT::isChild(NodeCT* child){
    auto it = std::find(this->children.begin(), this->children.end(), child);
    return it != children.end();
}

int NodeCT::getIndex(){ return this->index; }

int NodeCT::getThreshold1(){ return this->threshold1; }

int NodeCT::getThreshold2(){ return this->threshold1; }

int NodeCT::getLevel(){ return this->threshold2; }

void NodeCT::setLevel(int level){ this->threshold2 = level; }


void NodeCT::setArea(long int area){this->areaCC = area;}

long int NodeCT::getArea(){return this->areaCC;} 

NodeCT* NodeCT::getParent(){  return this->parent; }

void NodeCT::setParent(NodeCT* parent){ this->parent = parent; }

std::list<int>& NodeCT::getCNPs(){  return this->cnps; }

std::list<int> NodeCT::getCNPsCopy(){ 
    std::list<int> cnpsCopy;
    std::copy( this->cnps.begin(), this->cnps.end(), std::back_inserter(cnpsCopy) );
    return cnpsCopy; 
}

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
