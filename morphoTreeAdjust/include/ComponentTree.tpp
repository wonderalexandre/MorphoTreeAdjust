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
#include <algorithm> 

#include "../include/NodeCT.hpp"
#include "../include/ComponentTree.hpp"
#include "../include/AdjacencyRelation.hpp"
#include "../include/FlatZonesGraph.hpp"


template <typename CNPsType>
int* ComponentTree<CNPsType>::countingSort(ImageUInt8Ptr imgPtr){
	int n = this->numRows * this->numCols;
    auto img = imgPtr->rawData();
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
int* ComponentTree<CNPsType>::createTreeByUnionFind(int* orderedPixels, ImageUInt8Ptr imgPtr) {
	int numPixels = numRows*numCols;
    int *zPar = new int[numPixels];
	int *parent = new int[numPixels];
	auto img = imgPtr->rawData();
    
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
void ComponentTree<CNPsType>::reconstruction(NodeCTPtr<CNPsType> node, uint8_t* imgOut) {
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
ComponentTree<CNPsType>::ComponentTree(ImageUInt8Ptr img, bool isMaxtree, AdjacencyRelationPtr adj) : numRows(img->getNumRows()), numCols(img->getNumCols()), maxtreeTreeType(isMaxtree), adj(adj){   
    pixelToNode.resize(numRows*numCols, nullptr);
    build(img);
}

template <typename CNPsType>
ComponentTree<CNPsType>::ComponentTree(ImageUInt8Ptr img, bool isMaxtree, AdjacencyRelationPtr adj, std::shared_ptr<FlatZonesGraph> graph) : numRows(img->getNumRows()), numCols(img->getNumCols()), maxtreeTreeType(isMaxtree), adj(adj), flatzoneGraph(graph){   
    pixelToNode.resize(numRows*numCols, nullptr);
    build(img);
}


template <typename CNPsType>
void ComponentTree<CNPsType>::build(ImageUInt8Ptr imgPtr){ 
	int* orderedPixels = countingSort(imgPtr);
	int* parent = createTreeByUnionFind(orderedPixels, imgPtr);
    auto img = imgPtr->rawData();
	this->numNodes = 0;
    int numPixels = numRows*numCols;
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
	
    this->assignCNPs();
	computerArea(this->root); //computer area

	delete[] parent;
	delete[] orderedPixels;
}

template <>
inline void ComponentTreeP::assignCNPs() {
    int numPixels = numRows*numCols;
    for (int p = 0; p < numPixels; p++) {
        pixelToNode[p]->addCNPs(p);
    }
}

/*
Esse método faz atribuição dos pixels as suas respectivas flatzones.
O menor pixel (o primeiro a ser visitado durante a construção da flatzone) é considerado o id da flatzone.
Essa propriedade será mantida no grafo ao realizar operações de fusão de flatzones
*/
template <>
inline void ComponentTreeFZ::assignCNPs() {
    
    u_int32_t *numFlatzones = new u_int32_t[this->numNodes]();
    for (FlatZone& flatZone : this->flatzoneGraph->getFlatzones()) {
        numFlatzones[this->pixelToNode[flatZone.front()]->getIndex()] += 1;        
    }
    for (FlatZone flatZone : this->flatzoneGraph->getFlatzones()) {
        this->pixelToNode[flatZone.front()]->addCNPsOfDisjointFlatzone(std::move(flatZone), nullptr, numFlatzones[this->pixelToNode[flatZone.front()]->getIndex()]);
    }
    delete[] numFlatzones;
    
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
AdjacencyRelationPtr ComponentTree<CNPsType>::getAdjacencyRelation(){
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
    FlatZone flatzonesPruned;

    if (rootSubtree != this->root) {
        NodeFZPtr parent = rootSubtree->getParent();
        parent->getChildren().remove(rootSubtree);
        rootSubtree->setParent(nullptr);

        std::queue<NodeFZPtr> queue;
        queue.push(rootSubtree);
        int idMinimum = rootSubtree->getCNPsByFlatZone().begin()->first; // pega um idFlatzone qualquer
        // Coletar todas as flatzones que serão fundidas
        while (!queue.empty()) {
            NodeFZPtr node = queue.front(); queue.pop();
            this->numNodes--;

            for (NodeFZPtr child : node->getChildren()) {
                queue.push(child);
            }
            
            auto& cnps = node->getCNPsByFlatZone();
            for (auto it = cnps.begin(); it != cnps.end(); ) {
                for (const int& p : it->second) {
                    this->setSC(p, parent);
                }
                if(it->first < idMinimum) {
                    idMinimum = it->first; 
                    flatzonesPruned.splice(flatzonesPruned.begin(), it->second); 
                }else{
                    flatzonesPruned.splice(flatzonesPruned.end(), it->second); 
                }
                it = cnps.erase(it); // remove e avança
            }
            
        }

    
        int idFlatzonePruned = flatzonesPruned.front();
        int idFlatzonePrunedRep = flatzoneGraph->findRepresentative(idFlatzonePruned);// pega o id da flatzone
        auto& cnps = parent->getCNPsByFlatZone();
        for (auto it = cnps.begin(); it != cnps.end(); ) {
            int idFlatzone = it->first;
            int idFlatzoneRep = flatzoneGraph->findRepresentative(idFlatzone);
            if (idFlatzoneRep == idFlatzonePrunedRep) {
                if (idMinimum > idFlatzone) {
                    flatzonesPruned.splice(flatzonesPruned.begin(), it->second); // flatzone contem o id da flatzone
                    idMinimum = idFlatzone; 
                } else {
                    flatzonesPruned.splice(flatzonesPruned.end(), it->second); // flatzonePruned contem o id da flatzone
                }
                
                it = cnps.erase(it); // remove e avança
            } else {
                ++it; // apenas avança se não remover
            }
        }

        // Adiciona a flatzone pruned ao parent
        cnps[flatzonesPruned.front()] = std::move(flatzonesPruned);

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
inline void ComponentTreeFZ::mergeWithParent(NodeFZPtr node) {
    if (!node || node->getParent() == nullptr) return;

    NodeFZPtr parent = node->getParent();
    std::list<NodeFZPtr>& childrenParent = parent->getChildren();
    childrenParent.remove(node);
    this->numNodes--;

    for (NodeFZPtr child : node->getChildren()) {
        childrenParent.push_back(child);
        child->setParent(parent);
    }

    auto& cnps = node->getCNPsByFlatZone();
    auto& cnpsParent = parent->getCNPsByFlatZone();
    for (auto it = cnps.begin(); it != cnps.end(); ) {
        FlatZone& fz = it->second;
        for(int pixel: fz){
            this->pixelToNode[pixel] = parent;	
        }
        int idFZ = it->first;
        int idFZRep = flatzoneGraph->findRepresentative(idFZ);
        auto itFound = cnpsParent.find(idFZRep);
        if (itFound != cnpsParent.end()) { // Se idFZRep está em parent, então adicionar no final os pixels de fz nessa flatzone do parent
            FlatZone& targetParentFZ = itFound->second;
            targetParentFZ.splice(targetParentFZ.end(), fz);

            for (auto itParent = cnpsParent.begin(); itParent != cnpsParent.end(); ) {
                int idFZParent = itParent->first;
                int idFRepParent = flatzoneGraph->findRepresentative(idFZParent);
                if(idFZParent != idFRepParent) { //se idFZParent tem um novo represetante em parent
                    FlatZone& fzParent = itParent->second;
                    targetParentFZ.splice(targetParentFZ.end(), fzParent);
                    itParent = cnpsParent.erase(itParent); // remove duplicata
                }else{
                    ++itParent;
                }
            }
        }
        else{ //Se idFZRep não está em parent, então idFZRep é representando de flatzones que estão em parent
            for (auto itParent = cnpsParent.begin(); itParent != cnpsParent.end(); ) {
                int idFZParent = itParent->first;
                int idFRepParent = flatzoneGraph->findRepresentative(idFZParent);
                if(idFZRep == idFRepParent) { //o representando de FZParent é fz, então adicionar no final os pixels de FZParent em fz e atualizar as chaves de cnps
                    FlatZone& fzParent = itParent->second;
                    fz.splice(fz.end(), fzParent);
                    itParent = cnpsParent.erase(itParent); // remove chave antiga
                }else{
                    ++itParent;
                }
            }
            cnpsParent[idFZRep] = std::move(fz);      
        }
        it = cnps.erase(it);
    }

}



template <>
inline void ComponentTreeP::mergeWithParent(FlatZone* flatzone){
    int idFlatzone = flatzone->front();
    NodePPtr node = this->pixelToNode[idFlatzone];
    if(node->getNumCNPs() == static_cast<int>(flatzone->size())) {
        this->mergeWithParent(node);
    }
    else{
        NodePPtr parent = node->getParent();
        for( int p: *flatzone) {				
            parent->addCNPs(p);
            this->pixelToNode[p] = parent;	
        }

    }
}

template <>
inline void ComponentTreeFZ::mergeWithParent(FlatZone* fz) {
    int idFZ = fz->front();
    NodeFZPtr node = this->pixelToNode[idFZ];
    if (node->getNumFlatzone() == 1) {
        // Caso trivial: apenas um flat zone — usa merge tradicional
        this->mergeWithParent(node);
    } else {
        node->removeFlatzone(idFZ);
        int idFZRep = flatzoneGraph->findRepresentative(idFZ);
        NodeFZPtr parent = node->getParent();
        for(int pixel: *fz){
            this->pixelToNode[pixel] = parent;	
        }
        auto& cnpsParent = parent->getCNPsByFlatZone();

        auto itFound = cnpsParent.find(idFZRep);
        if (itFound != cnpsParent.end()) {
            // Se já existe o representante, funde em seu final
            FlatZone& targetParentFZ = itFound->second;
            targetParentFZ.splice(targetParentFZ.end(), *fz);

            // Remove duplicatas: todas as flat zones do mesmo grupo, exceto a chave do representante
            for (auto itParent = cnpsParent.begin(); itParent != cnpsParent.end(); ) {
                int idFZParent = itParent->first;
                int idFRepParent = flatzoneGraph->findRepresentative(idFZParent);
                if (idFZParent != idFZRep && idFRepParent == idFZRep) {
                    FlatZone& fzParent = itParent->second;
                    targetParentFZ.splice(targetParentFZ.end(), fzParent);
                    itParent = cnpsParent.erase(itParent);
                } else {
                    ++itParent;
                }
            }
        } else {
            // Não existe entrada para o grupo, então funde todas as flat zones do grupo (exceto fz)
            for (auto itParent = cnpsParent.begin(); itParent != cnpsParent.end(); ) {
                int idFZParent = itParent->first;
                int idFRepParent = flatzoneGraph->findRepresentative(idFZParent);
                if (idFZRep == idFRepParent) {
                    FlatZone& fzParent = itParent->second;
                    fz->splice(fz->end(), fzParent);
                    itParent = cnpsParent.erase(itParent);
                } else {
                    ++itParent;
                }
            }
            // Agora registra tudo sob a chave do representante
            cnpsParent[idFZRep] = std::move(*fz);
        }
    }
}


template <typename CNPsType>
std::vector<NodeCTPtr<CNPsType>> ComponentTree<CNPsType>::getNodesThreshold(int areaThreshold){
	std::vector<NodeCTPtr<CNPsType>> lista;
	std::queue<NodeCTPtr<CNPsType>> queue;
	queue.push(this->getRoot());
    
    int sumArea = 0; //somente para uso estatistico
    int numFlatZones=0; //somente para uso estatistico
    int numNodes = 0; //somente para uso estatistico
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
                numFlatZones += (node->computerNumFlatzoneDescendants() + node->getNumFlatzone()); 
                numNodes += node->computerNumDescendants() + 1;
            }
			lista.push_back(node);
	    }
	}
    if(PRINT_LOG){
        int areaImage = this->getNumColsOfImage() * this->getNumRowsOfImage();
        std::cout << "\tArea threshold: " << areaThreshold 
          << ", #Nodes: " << numNodes
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
ImageUInt8Ptr ComponentTree<CNPsType>::reconstructionImage(){
	ImageUInt8Ptr imgPtr = ImageUInt8::create(this->numRows, this->numCols);
	this->reconstruction(this->root, imgPtr->rawData());
	return imgPtr;
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
FlatZone& ComponentTreeFZ::getFlatzoneByID(int idFlatZone) {
    return pixelToNode[idFlatZone]->getFlatZone(idFlatZone);
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
std::shared_ptr<FlatZonesGraph>& ComponentTreeFZ::getFlatZonesGraph(){
    return this->flatzoneGraph;
}

