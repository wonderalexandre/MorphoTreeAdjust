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
ComponentTree<CNPsType>::ComponentTree(int numRows, int numCols, bool isMaxtree, double radiusOfAdjacencyRelation) {
    this->numRows = numRows;
    this->numCols = numCols;
    this->maxtreeTreeType = isMaxtree;
    this->adj = new AdjacencyRelation(numRows, numCols, radiusOfAdjacencyRelation);
    this->nodes = new NodeCT<CNPsType>*[numRows * numCols](); // Inicializa com `nullptr`
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

    delete[] nodes;
    nodes = nullptr;

    if constexpr (std::is_same_v<CNPsType, FlatZones>) {
        for (auto& ptr : pixelToFlatzone) {
            ptr = nullptr;  
        }
        pixelToFlatzone.clear();
        flatzoneGraph.clear();
    }
}

template <typename CNPsType>
void ComponentTree<CNPsType>::build(int* img){
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
			this->root = this->nodes[p] = new NodeCT<CNPsType>(this->numNodes++, nullptr, threshold1, threshold2);
            //if constexpr (std::is_same<CNPsType, Pixels>::value) this->nodes[p]->addCNPs(p);
		}
		else if (img[p] != img[parent[p]]) { //representante de um node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
			this->nodes[p] = new NodeCT<CNPsType>(this->numNodes++, this->nodes[parent[p]], threshold1, threshold2);
			//if constexpr (std::is_same<CNPsType, Pixels>::value) this->nodes[p]->addCNPs(p);
			this->nodes[parent[p]]->addChild(nodes[p]);
		}
		else if (img[p] == img[parent[p]]) {
			//if constexpr (std::is_same<CNPsType, Pixels>::value) this->nodes[parent[p]]->addCNPs(p);
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

template <>
inline void ComponentTreeP::assignCNPs() {
    int n = this->numRows * this->numCols;
    for (int p = 0; p < n; p++) {
        nodes[p]->addCNPs(p);
    }
}

template <>
inline void ComponentTreeFZ::assignCNPs() {
    this->pixelToFlatzone = std::vector<std::list<int>*>(numRows * numCols, nullptr);
    int n = this->numRows * this->numCols;
        std::vector<bool> visited(n, false); 
        for (int p = 0; p < n; p++) {
            if (visited[p]) continue; 
            if (this->nodes[p] == nullptr) {
                std::cerr << "❌ ERRO FATAL: Pixel " << p << " não foi mapeado para nenhum nó!" << std::endl;
                exit(0);
            }
            

            NodeFZ* node = this->nodes[p];
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
        this->buildFlatzoneGraph();
    
    
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::buildFlatzoneGraph() {
    
    int n = this->numRows * this->numCols;

    for (int p = 0; p < n; p++) {
        FlatZoneRef flatzoneP = this->getFlatzoneRef(p);

        // Se a flatzone ainda não foi registrada no grafo, cria uma entrada para ela
        if (this->flatzoneGraph.find(flatzoneP) == this->flatzoneGraph.end()) {
            this->flatzoneGraph[flatzoneP] = std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual>();
        }

        // Itera sobre os pixels vizinhos para encontrar vizinhos de flatzones
        for (int np : this->adj->getAdjPixels(p)) {
            FlatZoneRef flatzoneNP = this->getFlatzoneRef(np);

            if (&flatzoneNP.get() != &flatzoneP.get()) {
                // Se a flatzone vizinha ainda não foi registrada, cria uma entrada para ela
                if (this->flatzoneGraph.find(flatzoneNP) == this->flatzoneGraph.end()) {
                    this->flatzoneGraph[flatzoneNP] = std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual>();
                }

                // Conectar ambas no grafo
                this->flatzoneGraph[flatzoneP].insert(flatzoneNP);
                this->flatzoneGraph[flatzoneNP].insert(flatzoneP);
            }
        }
    }
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
    return this->nodes[p];
}

template <typename CNPsType>
void ComponentTree<CNPsType>::setSC(int p, NodeCT<CNPsType>* n){
	assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método setSC está acessando index fora dos limites");
	assert(n != nullptr && "Erro: o método setSC recebeu um ponteiro nulo");

	this->nodes[p] = n;
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

        std::list<int> cnpsCC;

        while (!s.empty()) {
            NodeP* child = s.top(); s.pop();
            this->numNodes--;
            for (NodeP* n : child->getChildren()) {
                s.push(n);
            }

            cnpsCC.splice(cnpsCC.end(), child->getCNPs());  // Mover flatzones para `cnpsCC`
            
            delete child;  
        }

        parent->getCNPs().splice(parent->getCNPs().end(), cnpsCC);
    }

}




template <>
inline void ComponentTreeFZ::prunning(NodeFZ* node) {
    assert(node != nullptr && "Erro: node is nullptr");
    assert(node->getParent() != nullptr && "Erro: node é a raiz");

    if (node != this->root) {
        NodeFZ* parent = node->getParent();
        parent->getChildren().remove(node);
        node->setParent(nullptr);

        std::stack<NodeFZ*> s;
        s.push(node);

        std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual> flatzonesToMergeSet;
        std::vector<NodeFZ*> toDelete;  // Lista de nós a serem deletados

        // Coletar todas as flatzones que serão fundidas
        while (!s.empty()) {
            NodeFZ* child = s.top(); s.pop();
            this->numNodes--;

            for (NodeFZ* n : child->getChildren()) {
                s.push(n);
            }

            for (FlatZone& flatZone : child->getCNPsByFlatZone()) {
                flatzonesToMergeSet.insert(std::ref(flatZone));
            }

            toDelete.push_back(child);  // Guardamos o nó para deletar depois
        }

        assert(!flatzonesToMergeSet.empty() && "ERRO: Nenhuma flatzone foi coletada para fusão!");

        if (flatzonesToMergeSet.size() == 1) {
            FlatZone& flatzoneToMerge = flatzonesToMergeSet.begin()->get();
            parent->addCNPsToConnectedFlatzone(std::move(flatzoneToMerge), this);
        } else {			
            FlatZone cnpsCC;  // Lista final de CNPs fundidos

            // Remover do grafo as flatzones que serão fundidas
            std::unordered_set<FlatZoneRef, ListRefHash, ListRefEqual> neighborsToTransfer;
            for (const auto& flatzoneRef : flatzonesToMergeSet) {
                FlatZone& flatzoneToMerge = flatzoneRef.get();

                // Coletar vizinhos que não fazem parte da fusão
                for (const auto& neighborRef : this->flatzoneGraph[flatzoneRef]) {
                    if (flatzonesToMergeSet.find(neighborRef) == flatzonesToMergeSet.end()) {
                        this->flatzoneGraph[neighborRef].erase(flatzoneRef);
                        neighborsToTransfer.insert(neighborRef);
                    }
                }

                // Remover do grafo
                this->flatzoneGraph.erase(flatzoneRef);

                // Adicionar pixels ao cnpsCC
                cnpsCC.splice(cnpsCC.end(), flatzoneToMerge);
            }

            assert(!cnpsCC.empty() && "ERRO: cnpsCC está vazio após a fusão!");

            // Adicionar `cnpsCC` ao grafo com suas conexões
            FlatZoneRef cnpsCCRef = cnpsCC;
            this->flatzoneGraph[cnpsCCRef] = neighborsToTransfer;

            for (const auto& neighborRef : neighborsToTransfer) {
                this->flatzoneGraph[neighborRef].insert(cnpsCCRef);
            }

            //Chamar `addCNPsToConnectedFlatzone` para finalizar a fusão
            parent->addCNPsToConnectedFlatzone(std::move(cnpsCC), this);
        }

        // Deletar os nós podados
        for (NodeFZ* child : toDelete) {
            delete child;
        }
    }
}

template <typename CNPsType>
std::vector<NodeCT<CNPsType>*> ComponentTree<CNPsType>::getNodesThreshold(int areaThreshold){
	std::vector<NodeCT<CNPsType>*> lista;
	std::stack<NodeCT<CNPsType>*> pilha;
	pilha.push(this->getRoot());

	while(!pilha.empty()) {
	    NodeCT<CNPsType>* node = pilha.top(); pilha.pop();
	    if(node->getArea() > areaThreshold) {
			for(NodeCT<CNPsType>* child: node->getChildren()) {
    			pilha.push(child);
    	    }
	    }
	    else {
			lista.push_back(node);
	    }
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
FlatZone& ComponentTreeFZ::getFlatzoneRef(int p) {
    assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getFlatzoneRef está acessando index fora dos limites");
    assert(pixelToFlatzone[p] != nullptr && "Erro: o método getFlatzoneRef tentou acessar flatzone nula!");
    assert(!pixelToFlatzone[p]->empty() && "Erro: flatzone vazia!");

    return *pixelToFlatzone[p];
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZone* ComponentTreeFZ::getFlatzonePointer(int p) {
    assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método getFlatzonePointer está acessando index fora dos limites");
    assert(pixelToFlatzone[p] != nullptr && "Erro: o método getFlatzonePointer tentou acessar flatzone nula!");
    assert(!pixelToFlatzone[p]->empty() && "Erro: flatzone vazia!");

    return pixelToFlatzone[p];
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::updatePixelToFlatzone(int p, std::list<int>* newFlatzone) {
    assert(p >= 0 && p < (this->numRows * this->numCols) && "Error: o método updatePixelToFlatzone está acessando index fora dos limites");
    assert(newFlatzone != nullptr && "Erro: o método updatePixelToFlatzone recebeu um ponteiro nulo para uma flatzone!");
    assert(!newFlatzone->empty() && "Erro: flatzone vazia!");

    pixelToFlatzone[p] = newFlatzone;  
}

