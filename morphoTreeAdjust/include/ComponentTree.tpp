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
#include "../include/FlatZonesGraph.hpp"


template <typename CNPsType>
int* ComponentTree<CNPsType>::countingSort(int* img){
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

template <typename CNPsType>
int ComponentTree<CNPsType>::findRoot(int *zPar, int x) {
	if (zPar[x] == x)
		return x;
	else {
		zPar[x] = findRoot(zPar, zPar[x]);
		return zPar[x];
	}
}

template <typename CNPsType>
int* ComponentTree<CNPsType>::createTreeByUnionFind(int* orderedPixels, int* img) {
	int *zPar = new int[numPixels];
	int *parent = new int[numPixels];
		
	for (int p = 0; p < numPixels; p++) {
		zPar[p] =  -1;
	}
		
	for(int i=numPixels-1; i >= 0; i--){
		int p = orderedPixels[i];
		parent[p] = p;
		zPar[p] = p;
		for (int n : this->adj.getAdjPixels(p)) {
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
	for (int i = 0; i < numPixels; i++) {
		int p = orderedPixels[i];
		int q = parent[p];
				
		if(img[parent[q]] == img[q]){
			parent[p] = parent[q];
		}
	}
		
	delete[] zPar;
	return parent;		
}


template <typename CNPsType>
void ComponentTree<CNPsType>::reconstruction(NodeCTPtr<CNPsType> node, int* imgOut) {
    assert(node != nullptr && "Erro: Nó inválido passado para reconstrução!");
    assert(imgOut != nullptr && "Erro: Ponteiro de saída da imagem é nulo!");

    for (int p : node->getCNPs()) {
        imgOut[p] = node->getLevel();
    }

    for (NodeCTPtr<CNPsType> child : node->getChildren()) {
        reconstruction(child, imgOut);
    }
}

template <typename CNPsType>
ComponentTree<CNPsType>::ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation):
    numRows(numRows), numCols(numCols), 
    maxtreeTreeType(isMaxtree), numPixels(numRows*numCols), 
    adj(numRows, numCols, radiusOfAdjacencyRelation)
{   
    pixelToNode.resize(numPixels, nullptr);
}

template <typename CNPsType>
ComponentTree<CNPsType>::ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
    : ComponentTree<CNPsType>(numRows, numCols, isMaxtree, radiusOfAdjacencyRelation) {
    build(img);
}


template <typename CNPsType>
void ComponentTree<CNPsType>::build(int* img){ 
	int* orderedPixels = countingSort(img);
	int* parent = createTreeByUnionFind(orderedPixels, img);
	this->numNodes = 0;
	for (int i = 0; i < numPixels; i++) {
		int p = orderedPixels[i];
		if (p == parent[p]) { //representante do node raiz
			int threshold1 = this->isMaxtree()? 0 : 255;
			int threshold2 = img[p];
			this->root = this->pixelToNode[p] = std::make_shared< NodeCT<CNPsType> >(this->numNodes++, nullptr, threshold1, threshold2);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
			this->pixelToNode[p] = std::make_shared< NodeCT<CNPsType> >(this->numNodes++, this->pixelToNode[parent[p]], threshold1, threshold2);
			this->pixelToNode[parent[p]]->addChild(pixelToNode[p]);
		}
		else if (img[p] == img[parent[p]]) {
			this->pixelToNode[p] = pixelToNode[parent[p]];
		}else{
			std::cerr << "\n\n\n\nOps...falhou geral\n\n\n\n";
			break;
		}
	}
	
    this->assignCNPs(img);
	computerArea(this->root); //computer area

	delete[] parent;
	delete[] orderedPixels;
}

template <>
inline void ComponentTreeP::assignCNPs(int* img) {
    for (int p = 0; p < numPixels; p++) {
        pixelToNode[p]->addCNPs(p);
    }
}

/*
Esse método faz atribuição dos pixels as suas respectivas flatzones e depois constroí-se um grafo de flatzones.
O menor pixel (o primeiro a ser visitado durante a construção da flatzone) é considerado o id da flatzone.
Essa propriedade será mantida no grafo ao realizar operações de fusão de flatzones
*/
template <>
inline void ComponentTreeFZ::assignCNPs(int* img) {
    this->flatzoneGraph = FlatZonesGraph::createInstance(img, this->numRows, this->numCols, this->adj);
    for (FlatZone& flatZone : this->flatzoneGraph->getFlatzones()) {
        this->pixelToNode[flatZone.front()]->addCNPsOfDisjointFlatzone(std::move(flatZone));
    }


    assert([&]() {
        for(NodeFZPtr node: this->root->getIteratorBreadthFirstTraversal()){
            for (auto& [id, flatZone] : node->getCNPsByFlatZone()) {
                int minPixel = *std::min_element(flatZone.begin(), flatZone.end());
                if (flatZone.front() != minPixel) {
                    return false;
                }
            }
        }
        return true;
    }() && "ERRO: Alguma flatzone não possui como ID o menor pixel!");

}

template <typename CNPsType>
void ComponentTree<CNPsType>::computerArea(NodeCTPtr<CNPsType> node){
	long int area = node->getNumCNPs();
	for(NodeCTPtr<CNPsType> child: node->getChildren()){
		computerArea(child);
		area += child->getArea();
	}
	node->setArea(area);
}

template <typename CNPsType>
NodeCTPtr<CNPsType> ComponentTree<CNPsType>::getRoot() {
	return this->root;
}

template <typename CNPsType>
void ComponentTree<CNPsType>::setRoot(NodeCTPtr<CNPsType> n){
	this->root = n;
}

template <typename CNPsType>
NodeCTPtr<CNPsType> ComponentTree<CNPsType>::getSC(int p) {
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getSC está acessando index fora dos limites");
    return this->pixelToNode[p];
}

template <typename CNPsType>
void ComponentTree<CNPsType>::setSC(int p, NodeCTPtr<CNPsType> n){
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método setSC está acessando index fora dos limites");
	assert(n != nullptr && "Erro: o método setSC recebeu um ponteiro nulo");

	this->pixelToNode[p] = n;
}

template <typename CNPsType>
AdjacencyRelation& ComponentTree<CNPsType>::getAdjacencyRelation(){
	return this->adj;
}

template <typename CNPsType>
bool ComponentTree<CNPsType>::isMaxtree(){
	return this->maxtreeTreeType;
}

template <typename CNPsType>
int ComponentTree<CNPsType>::getNumNodes(){
	return this->numNodes;
}

template <typename CNPsType>
int ComponentTree<CNPsType>::getNumRowsOfImage(){
	return this->numRows;
}

template <typename CNPsType>
int ComponentTree<CNPsType>::getNumColsOfImage(){
	return this->numCols;
}

template <>
inline void ComponentTreeP::prunning(NodePPtr node) {
    assert(node != nullptr && "Erro: node is nullptr");
	assert(node->getParent() != nullptr && "Erro: node is root");

    if (node != this->root) {
        NodePPtr parent = node->getParent();
        parent->getChildren().remove(node); 
        node->setParent(nullptr); 

        std::stack<NodePPtr> s;
        s.push(node);
        while (!s.empty()) {
            NodePPtr child = s.top(); s.pop();
            this->numNodes--;
            for (NodePPtr n : child->getChildren()) {
                s.push(n);
            }
            for (int p : child->getCNPs()) {
                parent->addCNPs(p);
                this->setSC(p, parent);
            }
            
            //delete child;  
            child = nullptr;  
        }

    }

}

template <>
inline void ComponentTreeFZ::prunning(NodeFZPtr rootSubtree) {
    
    assert(rootSubtree != nullptr && "Erro: node is nullptr");
    assert(rootSubtree->getParent() != nullptr && "Erro: node é a raiz");
    std::list<FlatZone> flatZoneList;
    if (rootSubtree != this->root) {
        NodeFZPtr parent = rootSubtree->getParent();
        parent->getChildren().remove(rootSubtree);
        rootSubtree->setParent(nullptr);

        std::queue<NodeFZPtr> queue;
        queue.push(rootSubtree);

        // Coletar todas as flatzones que serão fundidas
        while (!queue.empty()) {
            NodeFZPtr node = queue.front(); queue.pop();
            this->numNodes--;

            for (NodeFZPtr child : node->getChildren()) {
                queue.push(child);
            }
            
            for (auto& [id, flatzone] : node->getCNPsByFlatZone()) {
                flatZoneList.push_back(flatzone);
            }
            
            if(node != rootSubtree){
                node = nullptr;
            }
        }
        if(flatZoneList.size() > 1){
            //unifica todas as flatzones em uma unica flatzone na raiz da subtree e atualiza o grafo
            FlatZone unifiedFlatzone;
            flatzoneGraph->updateGraphAfterPruning(flatZoneList, unifiedFlatzone, rootSubtree, this->shared_from_this());
            parent->addCNPsToConnectedFlatzone(std::move(unifiedFlatzone), this->shared_from_this());
        }else{
            parent->addCNPsToConnectedFlatzone(std::move(flatZoneList.front()), this->shared_from_this());
        }
        rootSubtree = nullptr;
    }
}



template <>
inline void ComponentTreeP::mergeWithParent(NodePPtr node){
    if(node->getParent() != nullptr){
        NodePPtr parent = node->getParent();
        std::list<NodePPtr>& childrenParent = parent->getChildren();
        childrenParent.remove( node );			
        this->numNodes--;

        for( int p: node->getCNPs() ) {				
            parent->addCNPs(p);
            this->pixelToNode[p] = parent;	
        }

        for(NodePPtr child : node->getChildren()) {							
            childrenParent.push_back(child);				
            child->setParent(parent);			
        }	
        
        node = nullptr;
    }
}

template <>
inline void ComponentTreeFZ::mergeWithParent(NodeFZPtr node){
    if(node->getParent() != nullptr){
        NodeFZPtr parent = node->getParent();
        std::list<NodeFZPtr>& childrenParent = parent->getChildren();
        childrenParent.remove( node );			
        this->numNodes--;

        for(auto& [id, flatzone]: node->getCNPsByFlatZone()){
            if(flatzoneGraph->isAdjacent(id, parent)){
                parent->addCNPsToConnectedFlatzone(std::move(flatzone), this->shared_from_this());
            }else{
                parent->addCNPsOfDisjointFlatzone(std::move(flatzone), this->shared_from_this());
            }
        }
        
        for(NodeFZPtr child : node->getChildren()) {							
            childrenParent.push_back(child);				
            child->setParent(parent);			
        }			

        node = nullptr;
    }
}




template <typename CNPsType>
std::vector<NodeCTPtr<CNPsType>> ComponentTree<CNPsType>::getNodesThreshold(int areaThreshold){
	std::vector<NodeCTPtr<CNPsType>> lista;
	std::queue<NodeCTPtr<CNPsType>> queue;
	queue.push(this->getRoot());
    
    int sumArea = 0; //somente para uso estatistico
    int numFlatZones=0; //somente para uso estatistico

	while(!queue.empty()) {
	    NodeCTPtr<CNPsType> node = queue.front(); queue.pop();
	    if(node->getArea() > areaThreshold) {
			for(NodeCTPtr<CNPsType> child: node->getChildren()) {
    			queue.push(child);
    	    }
	    }
	    else {
            if(PRINT_LOG){ //somente para uso estatistico
                sumArea += node->getArea(); 
                numFlatZones += node->getNumFlatzone(); 
            }
			lista.push_back(node);
	    }
	}
    if(PRINT_LOG){
        int areaImage = this->getNumColsOfImage() * this->getNumRowsOfImage();
        std::cout << "\tArea threshold: " << areaThreshold 
          << ", #Nodes: " << lista.size() 
          << ", #FlatZones: " << numFlatZones
          << ", #InputTreeNodes: " << this->getNumNodes()
          << ", |Pruning Area|: " << sumArea 
          << " (" << std::fixed << std::setprecision(2) 
          << (static_cast<double>(sumArea) / areaImage) * 100.0 << "% of the image area)" 
          << std::endl;
    }
	return lista;
}

template <typename CNPsType>
std::vector<NodeCTPtr<CNPsType>> ComponentTree<CNPsType>::getLeaves(){
    std::vector<NodeCTPtr<CNPsType>> leaves;
    std::stack<NodeCTPtr<CNPsType>> s;
    s.push(this->root);
    
    while (!s.empty()) {
        NodeCTPtr<CNPsType> node = s.top(); s.pop();
        if (node->getChildren().empty()) {
            leaves.push_back(node);  // Adiciona folha ao vetor
        } else {
            for (NodeCTPtr<CNPsType> child : node->getChildren()) {
                s.push(child);
            }
        }
    }
    return leaves;
}

template <typename CNPsType>
int* ComponentTree<CNPsType>::reconstructionImage(){
	int n = this->numRows * this->numCols;
	int* img = new int[n];
	this->reconstruction(this->root, img);
	return img;
}



template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZone& ComponentTreeFZ::getFlatzoneByID(int idFlatZone) {
    return pixelToNode[idFlatZone]->getFlatZone(idFlatZone);
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
ListOfAdjacentFlatzones& ComponentTreeFZ::getListOfAdjacentFlatzones(){
    return this->flatzoneGraph->getGraph();
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
std::unique_ptr<FlatZonesGraph>& ComponentTreeFZ::getFlatZonesGraph(){
    return this->flatzoneGraph;
}

