#include <list>
#include <array>
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
	this->adj = new AdjacencyRelation(numRows, numCols, radiusOfAdjacencyRelation);	
	this->nodes = new NodeCT*[numRows*numCols]();
 }

ComponentTree::ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
	: ComponentTree(numRows, numCols, isMaxtree, radiusOfAdjacencyRelation){
	
	build(img);
} 


 ComponentTree::~ComponentTree() {
    delete this->adj;  // Libera a estrutura auxiliar (supondo que seja alocada dinamicamente)
    if (this->root) {
        std::stack<NodeCT*> s;
        s.push(this->root);
        while (!s.empty()) {
            NodeCT* node = s.top();s.pop();
            for (NodeCT* child : node->getChildren()) {
                s.push(child);
            }
            delete node; 
			node = nullptr;
        }
    }
	delete[] nodes;
    nodes = nullptr; 
}


void ComponentTree::build(int* img){
	int n = this->numRows * this->numCols;
	
	int* orderedPixels = countingSort(img);
	int* parent = createTreeByUnionFind(orderedPixels, img);
	int countInitialized = 0;
	this->numNodes = 0;
	for (int i = 0; i < n; i++) {
		int p = orderedPixels[i];
		if (p == parent[p]) { //representante do node raiz
			int threshold1 = this->isMaxtree()? 0 : 255;
			int threshold2 = img[p];
			this->root = this->nodes[p] = new NodeCT(this->numNodes++, nullptr, threshold1, threshold2);
			this->nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
			this->nodes[p] = new NodeCT(this->numNodes++, this->nodes[parent[p]], threshold1, threshold2);
			this->nodes[p]->addCNPs(p);
			this->nodes[parent[p]]->addChild(nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			this->nodes[parent[p]]->addCNPs(p);
			this->nodes[p] = nodes[parent[p]];
		}else{
			std::cerr << "\n\n\n\nOps...falhou geral\n\n\n\n";
			break;
		}

		if (this->nodes[p] != nullptr) {
            countInitialized++;
        }

	}
	if(countInitialized != n)
		std::cerr << "DEBUG: Total de nós inicializados: " << countInitialized  << " / " << n << std::endl;


	//computer area
	computerArea(this->root);

	delete[] parent;
	delete[] orderedPixels;
}

void ComponentTree::computerArea(NodeCT* node){
	long int area = node->getCNPs().size();
	for(NodeCT* child: node->getChildren()){
		computerArea(child);
		area += child->getArea();
	}
	node->setArea(area);
}

	
NodeCT* ComponentTree::getRoot() {
	return this->root;
}

void ComponentTree::setRoot(NodeCT* n){
	this->root = n;
}

 /*NodeCT* ComponentTree::getSC(int p) {
	return this->nodes[p];
}*/
NodeCT* ComponentTree::getSC(int p) {
    if (p < 0 || p >= (this->numRows * this->numCols)) { 
        std::cerr << "Error function getSC: Índice " << p 
                  << " fora dos limites! (numNodes=" << this->numNodes
                  << ", matriz max=" << (this->numRows * this->numCols) << ")"
                  << std::endl;
        return nullptr; 
    }
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

void ComponentTree::prunning(NodeCT* node){
	if(node == nullptr){
		std::cout << "node is nullptr" << std::endl;
	}

	if(node != this->root){
		NodeCT* parent = node->getParent();
		parent->getChildren().remove(node);
        node->setParent(nullptr);

		std::stack<NodeCT*> s;
		s.push(node);
		while(!s.empty()){
    		NodeCT* child = s.top(); s.pop();	
			this->numNodes--;
			for(NodeCT* n: child->getChildren()){
				s.push(n);
			}
			
			for(int p: child->getCNPs()){
				this->nodes[p] = parent;
			}
			parent->getCNPs().splice(parent->getCNPs().end(), child->getCNPs());
			
			delete child;
			child = nullptr;
		}
		node = nullptr;
		
	}

	
		
}



std::vector<NodeCT*> ComponentTree::getNodesThreshold(int areaThreshold){
	std::vector<NodeCT*> lista;
	std::stack<NodeCT*> pilha;
	pilha.push(this->getRoot());

	while(!pilha.empty()) {
	    NodeCT* node = pilha.top(); pilha.pop();
	    if(node->getArea() > areaThreshold) {
			for(NodeCT* child: node->getChildren()) {
    			pilha.push(child);
    	    }
	    }
	    else {
			lista.push_back(node);
	    }
	}
	return lista;
}

std::vector<NodeCT*> ComponentTree::getLeaves(){
    std::vector<NodeCT*> leaves;
    std::stack<NodeCT*> s;
    s.push(this->root);
    
    while (!s.empty()) {
        NodeCT* node = s.top(); s.pop();
        if (node->getChildren().empty()) {
            leaves.push_back(node);  // Adiciona folha ao vetor
        } else {
            for (NodeCT* child : node->getChildren()) {
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
	