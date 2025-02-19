#include <list>
#include <array>
#include <queue>
#include <vector>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <utility>

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"

//#define NDEBUG  // Remove os asserts do código
#include <cassert>


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
	this->nodes = new NodeCT*[numRows*numCols](); // Inicializa com `nullptr`
	this->pixelToFlatzone = std::vector<std::list<int>*>(numRows * numCols, nullptr);
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

	//delete[] pixelToFlatzone;
	//pixelToFlatzone = nullptr; 
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
			//this->nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
			this->nodes[p] = new NodeCT(this->numNodes++, this->nodes[parent[p]], threshold1, threshold2);
			//this->nodes[p]->addCNPs(p);
			this->nodes[parent[p]]->addChild(nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			//this->nodes[parent[p]]->addCNPs(p);
			this->nodes[p] = nodes[parent[p]];
		}else{
			std::cerr << "\n\n\n\nOps...falhou geral\n\n\n\n";
			break;
		}

		if (this->nodes[p] != nullptr) {
            countInitialized++;
        }

	}
	if(countInitialized != n){
		std::cerr << "❌ ERRO FATAL: Mapeamento SC falhou, nem todos os pixels foram mapeados: " << countInitialized  << " / " << n << std::endl;
        exit(0);
	}

	this->assignCNPs();
	
	computerArea(this->root); //computer area

	delete[] parent;
	delete[] orderedPixels;
}

void ComponentTree::assignCNPs() {
    int n = this->numRows * this->numCols;

    std::vector<bool> visited(n, false); 
    for (int p = 0; p < n; p++) {
        if (visited[p]) continue; 

		if (this->nodes[p] == nullptr) {
			std::cerr << "❌ ERRO FATAL: Pixel " << p << " não foi mapeado para nenhum nó!" << std::endl;
			exit(0);
		}
		

        NodeCT* node = this->nodes[p];
        std::list<int> flatZone;
        std::queue<int> queue;
        queue.push(p);
        visited[p] = true;

        while (!queue.empty()) {
            int q_p = queue.front(); queue.pop();
            flatZone.push_back(q_p);
            for (int np : this->adj->getAdjPixels(q_p)) {
				if (!visited[np] && this->nodes[np] == node) {
                    visited[np] = true;
                    queue.push(np);
                }
            }
        }
		assert(flatZone.size() > 0 && "ERRO: Existem flatzones vazias!");
		node->addCNPsOfDisjointFlatzone(std::move(flatZone), this);
	}
	
	assert([this]() -> bool {
		for (int p = 0; p < this->numRows * this->numCols; p++) {
			if (!this->pixelToFlatzone[p] || this->pixelToFlatzone[p]->empty()) {
				std::cerr << "ERRO: O pixel " << p 
						  << " pertence ao nó ID " << this->nodes[p]->getIndex()
						  << " mas não está corretamente mapeado para uma flatzone válida!" << std::endl;
				return false; // Retorna `false` para falhar no assert
			}
		}
		return true; // Retorna `true` se tudo estiver correto
	}());
		
    
}



void ComponentTree::computerArea(NodeCT* node){
	long int area = node->getNumCNPs();
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

 
NodeCT* ComponentTree::getSC(int p) {
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getSC está acessando index fora dos limites");
    return this->nodes[p];
}


void ComponentTree::setSC(int p, NodeCT* n){
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método setSC está acessando index fora dos limites");
	assert(n != nullptr && "Erro: o método setSC recebeu um ponteiro nulo");

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


void ComponentTree::prunning(NodeCT* node) {
    assert(node != nullptr && "Erro: node is nullptr");
	assert(node->parent != nullptr && "Erro: node is root");

    if (node != this->root) {
        NodeCT* parent = node->getParent();
        parent->getChildren().remove(node); 
        node->setParent(nullptr); 

        std::stack<NodeCT*> s;
        s.push(node);

        std::list<int> cnpsCC;

        while (!s.empty()) {
            NodeCT* child = s.top(); s.pop();
            this->numNodes--;
            for (NodeCT* n : child->getChildren()) {
                s.push(n);
            }

            for (std::list<int>& flatZone : child->moveCNPsByFlatZone()) {
                cnpsCC.splice(cnpsCC.end(), flatZone);  // Mover flatzones para `cnpsCC`
            }

            delete child;  
        }

        parent->addCNPsToConnectedFlatzone(std::move(cnpsCC), this);
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
	

std::list<int>& ComponentTree::getFlatzoneRef(int p) {
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getFlatzoneRef está acessando index fora dos limites");
    assert(pixelToFlatzone[p] != nullptr && "Erro: o método getFlatzoneRef tentou acessar flatzone nula!");
	assert(!pixelToFlatzone[p]->empty() && "Erro: flatzone vazia!");
    return *pixelToFlatzone[p];
}

std::list<int>* ComponentTree::getFlatzonePointer(int p){
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getFlatzonePointer está acessando index fora dos limites");
    assert(pixelToFlatzone[p] != nullptr && "Erro: o método getFlatzonePointer tentou acessar flatzone nula!");
	assert(!pixelToFlatzone[p]->empty() && "Erro: flatzone vazia!");
    return pixelToFlatzone[p];
}

void ComponentTree::updatePixelToFlatzone(int p, std::list<int>* newFlatzone) {
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método updatePixelToFlatzone está acessando index fora dos limites");
    assert(newFlatzone != nullptr && "Erro: o método updatePixelToFlatzone recebeu um ponteiro nulo para uma flatzone!");
	assert(!newFlatzone->empty() && "Erro: flatzone vazia!");
	
    pixelToFlatzone[p] = newFlatzone;  
}

