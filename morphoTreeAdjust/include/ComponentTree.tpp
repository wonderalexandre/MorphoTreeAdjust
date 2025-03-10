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
void ComponentTree<CNPsType>::reconstruction(NodeCT<CNPsType>* node, int* imgOut) {
    assert(node != nullptr && "Erro: Nó inválido passado para reconstrução!");
    assert(imgOut != nullptr && "Erro: Ponteiro de saída da imagem é nulo!");

    for (int p : node->getCNPs()) {
        imgOut[p] = node->getLevel();
    }

    for (NodeCT<CNPsType>* child : node->getChildren()) {
        reconstruction(child, imgOut);
    }
}

template <typename CNPsType>
ComponentTree<CNPsType>::ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation){   
    this->numRows = numRows;
    this->numCols = numCols;
    this->maxtreeTreeType = isMaxtree;
    this->numPixels = numRows*numCols;
    this->adj = new AdjacencyRelation(numRows, numCols, radiusOfAdjacencyRelation); 
    this->pixelToNode = new NodeCT<CNPsType>*[numRows * numCols]();
}

template <typename CNPsType>
ComponentTree<CNPsType>::ComponentTree(int* img, int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation)
    : ComponentTree<CNPsType>(numRows, numCols, isMaxtree, radiusOfAdjacencyRelation) {
    build(img);
}


template <typename CNPsType>
ComponentTree<CNPsType>::~ComponentTree() {
    delete this->adj;  

    if (this->root) {
        std::stack<NodeCT<CNPsType>*> s;
        s.push(this->root);
        while (!s.empty()) {
            NodeCT<CNPsType>* node = s.top();
            s.pop();
            for (NodeCT<CNPsType>* child : node->getChildren()) {
                s.push(child);
            }
            delete node; 
            node = nullptr;
        }
    }

    delete[] pixelToNode;
    pixelToNode = nullptr;

    if constexpr (std::is_same_v<CNPsType, FlatZones>) {
        for (int i = 0; i < numPixels; i++) {
                delete flatzoneGraph[i];  
            }
            delete[] flatzoneGraph;  
            flatzoneGraph = nullptr; 
        }

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
			this->root = this->pixelToNode[p] = new NodeCT<CNPsType>(this->numNodes++, nullptr, threshold1, threshold2);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
			this->pixelToNode[p] = new NodeCT<CNPsType>(this->numNodes++, this->pixelToNode[parent[p]], threshold1, threshold2);
			this->pixelToNode[parent[p]]->addChild(pixelToNode[p]);
		}
		else if (img[p] == img[parent[p]]) {
			this->pixelToNode[p] = pixelToNode[parent[p]];
		}else{
			std::cerr << "\n\n\n\nOps...falhou geral\n\n\n\n";
			break;
		}
	}
	
    this->assignCNPs();
	computerArea(this->root); //computer area

	delete[] parent;
	delete[] orderedPixels;
}

template <>
inline void ComponentTreeP::assignCNPs() {
    for (int p = 0; p < numPixels; p++) {
        pixelToNode[p]->addCNPs(p);
    }
}

template <>
inline void ComponentTreeFZ::assignCNPs() {
    int* pixelToFlatzone = new int[numPixels];
    bool* visited = new bool[numPixels]();
    bool* isCountor = new bool[numPixels]();
    flatzoneGraph = new AdjacentFlatzones*[numPixels]();

    for (int p = 0; p < numPixels; p++) {
        if (visited[p]) continue; 
        assert(this->pixelToNode[p] != nullptr && "Falha no mapeamento SC");
            
        NodeFZ* node = this->pixelToNode[p];
        std::list<int> flatZone;
        std::queue<int> queue;
        queue.push(p);
        visited[p] = true;
        this->flatzoneGraph[p] = new AdjacentFlatzones();

        while (!queue.empty()) {
            int q = queue.front(); queue.pop();
            flatZone.push_back(q);
            
            for (int nq : this->adj->getAdjPixels(q)) {
                if (!visited[nq] && this->pixelToNode[nq] == node) {
                    visited[nq] = true;
                    queue.push(nq);
                }
                else if (pixelToNode[nq] != node){
                    isCountor[nq] = true;
                    isCountor[q] = true;
                    pixelToFlatzone[q] = p; //id da flatzone
                }
            }

        }
        assert(flatZone.size() > 0 && "ERRO: Existem flatzones vazias!");
        node->addCNPsOfDisjointFlatzone(std::move(flatZone), this);
    }

    //build graph
    for (int p = 0; p < numPixels; p++) {
        if(!isCountor[p]) continue;
        int flatZoneID_P = pixelToFlatzone[p];
        AdjacentFlatzones& setP = *flatzoneGraph[flatZoneID_P];
        for (int np : this->adj->getAdjPixels(p)) {
            if(!isCountor[np]) continue;
            int flatZoneID_NP = pixelToFlatzone[np];
            if (flatZoneID_NP != flatZoneID_P) {
                AdjacentFlatzones& setNP = *flatzoneGraph[flatZoneID_NP];
                setP.insert(flatZoneID_NP);
                setNP.insert(flatZoneID_P);
            }
        }
    }
   

    delete[] pixelToFlatzone;
    delete[] isCountor;
    delete[] visited;
}


template <typename CNPsType>
void ComponentTree<CNPsType>::computerArea(NodeCT<CNPsType>* node){
	long int area = node->getNumCNPs();
	for(NodeCT<CNPsType>* child: node->getChildren()){
		computerArea(child);
		area += child->getArea();
	}
	node->setArea(area);
}

template <typename CNPsType>
NodeCT<CNPsType>* ComponentTree<CNPsType>::getRoot() {
	return this->root;
}

template <typename CNPsType>
void ComponentTree<CNPsType>::setRoot(NodeCT<CNPsType>* n){
	this->root = n;
}

template <typename CNPsType>
NodeCT<CNPsType>* ComponentTree<CNPsType>::getSC(int p) {
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getSC está acessando index fora dos limites");
    return this->pixelToNode[p];
}

template <typename CNPsType>
void ComponentTree<CNPsType>::setSC(int p, NodeCT<CNPsType>* n){
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método setSC está acessando index fora dos limites");
	assert(n != nullptr && "Erro: o método setSC recebeu um ponteiro nulo");

	this->pixelToNode[p] = n;
}

template <typename CNPsType>
AdjacencyRelation* ComponentTree<CNPsType>::getAdjacencyRelation(){
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
inline void ComponentTreeP::prunning(NodeP* node) {
    assert(node != nullptr && "Erro: node is nullptr");
	assert(node->getParent() != nullptr && "Erro: node is root");

    if (node != this->root) {
        NodeP* parent = node->getParent();
        parent->getChildren().remove(node); 
        node->setParent(nullptr); 

        std::stack<NodeP*> s;
        s.push(node);
        while (!s.empty()) {
            NodeP* child = s.top(); s.pop();
            this->numNodes--;
            for (NodeP* n : child->getChildren()) {
                s.push(n);
            }
            for (int p : child->getCNPs()) {
                parent->addCNPs(p);
                this->setSC(p, parent);
            }
            
            delete child;  
            child = nullptr;  
        }

    }

}

template <>
inline void ComponentTreeFZ::prunning(NodeFZ* rootSubtree) {
    assert(rootSubtree != nullptr && "Erro: node is nullptr");
    assert(rootSubtree->getParent() != nullptr && "Erro: node é a raiz");
    std::list<FlatZoneNode> flatZoneNodeList;
    std::list<NodeFZ*> toRemove;
    if (rootSubtree != this->root) {
        NodeFZ* parent = rootSubtree->getParent();
        parent->getChildren().remove(rootSubtree);
        rootSubtree->setParent(nullptr);

        std::queue<NodeFZ*> queue;
        queue.push(rootSubtree);

        // Coletar todas as flatzones que serão fundidas
        while (!queue.empty()) {
            NodeFZ* node = queue.front(); queue.pop();
            this->numNodes--;

            for (NodeFZ* child : node->getChildren()) {
                queue.push(child);
            }

            toRemove.push_back(node);
            for (auto& [id, flatzone] : node->getCNPsByFlatZone()) {
                flatZoneNodeList.emplace_back(node, flatzone);
            }
        }
        
        //unifica todas as flatzones em uma unica flatzone na raiz da subtree e atualiza o grafo
        FlatZone unifiedFlatzone;
        updateGraphAfterPruning(flatZoneNodeList, unifiedFlatzone, rootSubtree);
        
        //a flatzone unificada será novamente unificada as flatzones conexas a ela presente no parent e atuliza o grafo
        parent->addCNPsToConnectedFlatzone(std::move(unifiedFlatzone), this);

        for(NodeFZ* node: toRemove){
            delete node;
            node = nullptr; 
        }
    }
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::updateGraphAfterPruning(std::list<FlatZoneNode>& flatZoneNodeList, FlatZone& unifiedFlatzone, NodeFZ* node) {
    assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");
    int unifiedFlatzoneID = node->getRepresentativeCNPs();
    
    for(FlatZoneNode& pair: flatZoneNodeList)   {
        NodeFZ* child = pair.node;
        FlatZone& flatZone = *pair.flatzone;
        //atualiza o grafo e coleta os cnps em uma unica flatzone
        int flatZoneID = getIdFlatZone(flatZone);
        std::unordered_set<int> neighborsCopy = *(flatzoneGraph[flatZoneID]);
        if(flatZoneID == unifiedFlatzoneID){
            for (int neighborID : neighborsCopy) {
                //vizinhos que serão unificados
                if( (this->maxtreeTreeType && node->getLevel() <= getSC(neighborID)->getLevel()) || (!this->maxtreeTreeType && node->getLevel() >= getSC(neighborID)->getLevel())){
                    this->flatzoneGraph[neighborID]->erase(unifiedFlatzoneID);
                    this->flatzoneGraph[unifiedFlatzoneID]->erase(neighborID);
                }
            }
            unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); 
        }else{
            for (int neighborID : neighborsCopy) {
                if(this->flatzoneGraph[neighborID] ){
                    //vizinhos que NÃO serão unificados
                    if( (this->maxtreeTreeType && node->getLevel() >= getSC(neighborID)->getLevel()) || (!this->maxtreeTreeType && node->getLevel() <= getSC(neighborID)->getLevel())){
                        this->flatzoneGraph[neighborID]->erase(flatZoneID);
                        this->flatzoneGraph[flatZoneID]->erase(neighborID);
                            
                        this->flatzoneGraph[neighborID]->insert(unifiedFlatzoneID);
                        this->flatzoneGraph[unifiedFlatzoneID]->insert(neighborID);
                    }
                }
            }
            unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); // Adicionar pixels ao cnpsCC
            
            delete this->flatzoneGraph[flatZoneID];
            this->flatzoneGraph[flatZoneID] = nullptr;
        }
    }
    
    assert(!unifiedFlatzone.empty() && "ERRO: unifiedFlatzone está vazio após a fusão!");
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::updateGraph(std::list<FlatZoneNode>& flatZoneNodeList,  FlatZone& unifiedFlatzone, NodeFZ* node) {
    assert(!flatZoneNodeList.empty() && "ERRO: Lista de FlatZoneNode está vazia!");

    int unifiedFlatzoneID = -1;
    for(FlatZoneNode& pair: flatZoneNodeList)   {
        NodeFZ* child = pair.node;
        if(child == node){
            FlatZone* flatZone = pair.flatzone;
            unifiedFlatzoneID = flatZone->front();
            unifiedFlatzone.splice(unifiedFlatzone.end(), *flatZone); 
        }
    }
    
    std::unordered_set<int>* unifiedFlatzoneSet = this->flatzoneGraph[unifiedFlatzoneID];
    

    for(FlatZoneNode& pair: flatZoneNodeList)   {
        NodeFZ* child = pair.node;
        FlatZone& flatZone = *pair.flatzone;
        //atualiza o grafo e coleta os cnps em uma unica flatzone

        int flatZoneID = getIdFlatZone(flatZone);

        std::unordered_set<int> neighborsCopy = *(flatzoneGraph[flatZoneID]);
        


        if(flatZoneID == unifiedFlatzoneID){
            for (int neighborID : neighborsCopy) {
                //vizinhos que serão unificados
                if( (this->maxtreeTreeType && node->getLevel() <= getSC(neighborID)->getLevel()) || (!this->maxtreeTreeType && node->getLevel() >= getSC(neighborID)->getLevel())){
                    this->flatzoneGraph[neighborID]->erase(unifiedFlatzoneID);
                    this->flatzoneGraph[unifiedFlatzoneID]->erase(neighborID);
                }
            }
        }else{
            for (int neighborID : neighborsCopy) {

                if(this->flatzoneGraph[neighborID] ){

                    //vizinhos que NÃO serão unificados
                    if( (this->maxtreeTreeType && node->getLevel() <= getSC(neighborID)->getLevel()) || (!this->maxtreeTreeType && node->getLevel() >= getSC(neighborID)->getLevel())){
                        this->flatzoneGraph[neighborID]->erase(flatZoneID);
                        this->flatzoneGraph[flatZoneID]->erase(neighborID);
                        
                        this->flatzoneGraph[unifiedFlatzoneID]->insert(neighborID);
                        this->flatzoneGraph[neighborID]->insert(unifiedFlatzoneID);
                        
                            
                    }
                }
            }
            delete this->flatzoneGraph[flatZoneID];
            this->flatzoneGraph[flatZoneID] = nullptr;
            unifiedFlatzone.splice(unifiedFlatzone.end(), flatZone); // Adicionar pixels ao cnpsCC
        }
    }
    
    assert(!unifiedFlatzone.empty() && "ERRO: unifiedFlatzone está vazio após a fusão!");
}




template <typename CNPsType>
std::vector<NodeCT<CNPsType>*> ComponentTree<CNPsType>::getNodesThreshold(int areaThreshold){
	std::vector<NodeCT<CNPsType>*> lista;
	std::stack<NodeCT<CNPsType>*> pilha;
	pilha.push(this->getRoot());
    int sumArea = 0;
	while(!pilha.empty()) {
	    NodeCT<CNPsType>* node = pilha.top(); pilha.pop();
	    if(node->getArea() > areaThreshold) {
			for(NodeCT<CNPsType>* child: node->getChildren()) {
    			pilha.push(child);
    	    }
	    }
	    else {
            sumArea += node->getArea();
			lista.push_back(node);
	    }
	}
    if(PRINT_LOG){
        int areaImage = this->getNumColsOfImage() * this->getNumRowsOfImage();
        std::cout << "\tArea threshold: " << areaThreshold 
          << ", #Nodes: " << lista.size() 
          << ", |Pruning Area|: " << sumArea 
          << " (" << std::fixed << std::setprecision(2) 
          << (static_cast<double>(sumArea) / areaImage) * 100.0 << "% of the image area)" 
          << std::endl;
    }
	return lista;
}

template <typename CNPsType>
std::vector<NodeCT<CNPsType>*> ComponentTree<CNPsType>::getLeaves(){
    std::vector<NodeCT<CNPsType>*> leaves;
    std::stack<NodeCT<CNPsType>*> s;
    s.push(this->root);
    
    while (!s.empty()) {
        NodeCT<CNPsType>* node = s.top(); s.pop();
        if (node->getChildren().empty()) {
            leaves.push_back(node);  // Adiciona folha ao vetor
        } else {
            for (NodeCT<CNPsType>* child : node->getChildren()) {
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
int ComponentTreeFZ::getIdFlatZone(const FlatZone& fz){
    return fz.front();
}



template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZone& ComponentTreeFZ::getFlatzoneByID(int idFlatZone) {
    return pixelToNode[idFlatZone]->getFlatZone(idFlatZone);
}
