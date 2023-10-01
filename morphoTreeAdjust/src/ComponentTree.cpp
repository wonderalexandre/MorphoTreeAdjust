#include <list>
#include <vector>
#include <stack>
#include <unordered_map>
#include <iostream>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"



int* ComponentTree::countingSort(int* img){
	int n = this->numRows * this->numCols;
	int maxvalue = img[0];
	for (int i = 1; i < n; i++)
		if(maxvalue < img[i]) maxvalue = img[i];
			
	std::vector<int> counter(maxvalue + 1, 0); 
	int *orderedPixels = new int[n];
		
	if(this->isMaxtree()){
		for (int i = 0; i < n; i++)
			counter[img[i]]++;

		for (int i = 1; i < maxvalue; i++) 
			counter[i] += counter[i - 1];
		counter[maxvalue] += counter[maxvalue-1];
		
		for (int i = n - 1; i >= 0; --i)
			orderedPixels[--counter[img[i]]] = i;	

	}else{
		for (int i = 0; i < n; i++)
			counter[maxvalue - img[i]]++;

		for (int i = 1; i < maxvalue; i++) 
			counter[i] += counter[i - 1];
		counter[maxvalue] += counter[maxvalue-1];

		for (int i = n - 1; i >= 0; --i)
			orderedPixels[--counter[maxvalue - img[i]]] = i;
	}
	
	return orderedPixels;
}

int ComponentTree::findRoot(int *zPar, int x) {
	if (zPar[x] == x)
		return x;
	else {
		zPar[x] = findRoot(zPar, zPar[x]);
		return zPar[x];
	}
}

int* ComponentTree::createTreeByUnionFind(int* orderedPixels, int* img) {
	const int n = this->numRows * this->numCols;
	int *zPar = new int[n];
	int *parent = new int[n];
		
	for (int p = 0; p < n; p++) {
		zPar[p] =  -1;
	}
		
	for(int i=n-1; i >= 0; i--){
		int p = orderedPixels[i];
		parent[p] = p;
		zPar[p] = p;
		for (int n : this->adj->getAdjPixels(p)) {
			if(zPar[n] != -1){
				int r = this->findRoot(zPar, n);
				if(p != r){
					parent[r] = p;
					zPar[r] = p;
				}
			}
		}
	}
			
	// canonizacao da arvore
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		int q = parent[p];
				
		if(img[parent[q]] == img[q]){
			parent[p] = parent[q];
		}
	}
		
	delete[] zPar;
	return parent;		
}

void ComponentTree::reconstruction(NodeCT* node, int* imgOut){
	for (int p : node->getCNPs()){
		imgOut[p] = node->getLevel();
	}
	for(NodeCT* child: node->getChildren()){
		reconstruction(child, imgOut);
	}
}

ComponentTree::ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation){
	this->numRows = numRows;
	this->numCols = numCols;
	this->maxtreeTreeType = isMaxtree;
	this->adj = new AdjacencyRelation(numRows, numCols, 1.5);	
 }

ComponentTree::ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
	: ComponentTree(numRows, numCols, isMaxtree, radiusOfAdjacencyRelation){
	
	build(img);
} 

ComponentTree::~ComponentTree(){
	delete this->adj;  
	std::stack<NodeCT*> s;
	s.push(this->root);
	while(!s.empty()){
    	NodeCT* node = s.top(); s.pop();
		for (NodeCT *child: node->getChildren()){
			s.push(child);
		}
		delete node;
	}
 }

void ComponentTree::build(int* img){
 
	int n = this->numRows * this->numCols;
	this->nodes = std::vector<NodeCT*>(n);
	
	int* orderedPixels = countingSort(img);
	int* parent = createTreeByUnionFind(orderedPixels, img);

	this->numNodes = 0;
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		if (p == parent[p]) { //representante do node raiz
			this->root = nodes[p] = new NodeCT(this->numNodes++, nullptr, img[p]);
			nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			nodes[p] = new NodeCT(this->numNodes++, nodes[parent[p]], img[p]);
			nodes[p]->addCNPs(p);
			nodes[parent[p]]->addChild(nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			nodes[parent[p]]->addCNPs(p);
			nodes[p] = nodes[parent[p]];
		}
	}


	delete[] parent;
	delete[] orderedPixels;
}


	
NodeCT* ComponentTree::getRoot() {
	return this->root;
}

void ComponentTree::setRoot(NodeCT* n){
	this->root = n;
}

 NodeCT* ComponentTree::getSC(int p) {
	return this->nodes[p];
}

void ComponentTree::setSC(int p, NodeCT* n){
	this->nodes[p] = n;
}

AdjacencyRelation* ComponentTree::getAdjacencyRelation(){
	return this->adj;
}

bool ComponentTree::isMaxtree(){
	return this->maxtreeTreeType;
}


int ComponentTree::getNumNodes(){
	return this->numNodes;
}

int ComponentTree::getNumRowsOfImage(){
	return this->numRows;
}

int ComponentTree::getNumColsOfImage(){
	return this->numCols;
}

bool ComponentTree::prunning(NodeCT* node){
	if(node != this->root){
		NodeCT* parent = node->getParent();
		parent->getChildren().remove(node);
		node->setParent(nullptr);
		std::stack<NodeCT*> s;
		s.push(node);
		while(!s.empty()){
    		NodeCT* child = s.top(); s.pop();	
			this->numNodes--;
			parent->getCNPs().splice(parent->getCNPs().end(), child->getCNPs());
			for(int p: child->getCNPs()){
				this->nodes[p] = parent;
			}
			for(NodeCT* n: child->getChildren()){
				s.push(n);
			}
		}
		delete node;
		
		return true;
	}
	return false;
	
		
}

std::list<NodeCT*> ComponentTree::getLeaves(){
	std::list<NodeCT*> leaves;
	std::stack<NodeCT*> s;
	s.push(this->root);
	while(!s.empty()){
    	NodeCT* node = s.top(); s.pop();
		if(node->getChildren().empty()){
			leaves.push_back(node);
		}else{
			for (NodeCT *child: node->getChildren()){
				s.push(child);
			}
		}
	}
	return leaves;
}

int* ComponentTree::reconstructionImage(){
	int n = this->numRows * this->numCols;
	int* img = new int[n];
	this->reconstruction(this->root, img);
	return img;
}
	