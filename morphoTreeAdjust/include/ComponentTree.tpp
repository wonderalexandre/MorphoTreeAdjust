#include <list>
#include <array>
#include <vector>
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
NodeCT<CNPsType> ComponentTree<CNPsType>::proxy(NodeId id) const {
    if (id < 0) 
        return NodeCT<CNPsType>(); // handle vazio
    else
        return NodeCT<CNPsType>(const_cast<ComponentTree<CNPsType>*>(this), id);
}


template<typename CNPsType>
NodeCT<CNPsType> ComponentTree<CNPsType>::createNode(int repNode, NodeCT<CNPsType> parent, int threshold1, int threshold2) {
    NodeId id = makeNode(repNode, parent ? parent.getIndex() : -1, threshold1, threshold2);
    return proxy(id);
}

template <typename CNPsType>
NodeId ComponentTree<CNPsType>::makeNode(int repNode, NodeId parentId, int threshold1, int threshold2){
    // Aloca ID contíguo
    NodeId id = this->arena.allocate(repNode, threshold1, threshold2);

    // Encadeia no pai (por ID) se houver
    if (parentId >= 0) {
        addChildById(parentId, id);
    }
    
    //contador de nós
    this->numNodes++;

    
    return id;
}

template<typename CNPsType>
inline void ComponentTree<CNPsType>::setParentById(NodeId nodeId, NodeId parentId) {
    if (parentId == arena.parentId[nodeId]) return;
    if (parentId == -1) {
        if (arena.parentId[nodeId] != -1) 
            removeChildById(arena.parentId[nodeId], parentId, false);
    } else {
        addChildById(parentId, nodeId);
    }        
}

template<typename CNPsType>
void ComponentTree<CNPsType>::addChildById(int parentId, int childId) {
    if (parentId < 0 || childId < 0) return;

    // Se o filho já tem pai, desconecta antes de ler P/C
    if (arena.parentId[childId] != -1) {
        removeChildById(arena.parentId[childId], childId, false);
    }

    arena.parentId[childId]      = parentId;
    arena.prevSiblingId[childId] = arena.lastChildId[parentId];
    arena.nextSiblingId[childId] = -1;

    if (arena.firstChildId[parentId] == -1) {
        arena.firstChildId[parentId] = arena.lastChildId[parentId] = childId;
    } else {
        arena.nextSiblingId[arena.lastChildId[parentId]] = childId;
        arena.lastChildId[parentId] = childId;
    }
    ++arena.childCount[parentId];

    
    //assert(ComponentTree::validateStructure(this) && "ComponentTree topology invariant failed after addChildById");
    
}

// Remove um filho 'childId' da lista encadeada de filhos do pai 'parentId'.
template<typename CNPsType>
inline void ComponentTree<CNPsType>::removeChildById(int parentId, int childId, bool release) {
    if (parentId < 0 || childId < 0) return;
    if (arena.parentId[childId] != parentId) return;

    const int prev = arena.prevSiblingId[childId];
    const int next = arena.nextSiblingId[childId];

    if (prev == -1) arena.firstChildId[parentId] = next;
    else            arena.nextSiblingId[prev] = next;

    if (next == -1) arena.lastChildId[parentId] = prev;
    else            arena.prevSiblingId[next] = prev;

    if (arena.childCount[parentId] > 0) 
        --arena.childCount[parentId];

    arena.parentId[childId] = -1;
    arena.prevSiblingId[childId] = -1;
    arena.nextSiblingId[childId] = -1;
    if(release){
        releaseNode(childId);
        //assert(ComponentTree::validateStructure(this) && "ComponentTree topology invariant failed after removeChildById");
    }
}


// Move todos os filhos de 'fromId' para o fim da lista de filhos de 'toId'.
template<typename CNPsType>
inline void ComponentTree<CNPsType>::spliceChildrenById(int toId, int fromId) {
    if (toId < 0 || fromId < 0 || toId == fromId) return;

    NodeId firstFrom = arena.firstChildId[fromId];
    if (firstFrom == -1) return; // nada para mover

    // 1) todos os filhos de 'fromId' passam a ter pai 'toId'
    for (int c = arena.firstChildId[fromId]; c != -1; c = arena.nextSiblingId[c]) {
        arena.parentId[c] = toId;
    }

    // 2) concatena a lista de filhos de 'fromId' no fim da lista de 'toId'
    if (arena.firstChildId[toId] == -1) {
        // 'toId' não tinha filhos — vira exatamente a lista de 'fromId'
        arena.firstChildId[toId] = arena.firstChildId[fromId];
        arena.lastChildId[toId]  = arena.lastChildId[fromId];
        // o primeiro filho já tem prevSiblingId == -1 porque era o primeiro de 'fromId'
    } else {
        // 'toId' já tinha filhos — encadeia no final
        arena.nextSiblingId[ arena.lastChildId[toId] ] = arena.firstChildId[fromId];
        arena.prevSiblingId[ arena.firstChildId[fromId] ] = arena.lastChildId[toId];
        arena.lastChildId[toId] = arena.lastChildId[fromId];
    }

    // 3) atualiza contadores
    arena.childCount[toId] += arena.childCount[fromId];

    // 4) zera a lista de 'fromId'
    arena.firstChildId[fromId] = -1;
    arena.lastChildId[fromId]  = -1;
    arena.childCount[fromId]   = 0;

   // assert(ComponentTree::validateStructure(this) && "ComponentTree topology invariant failed after spliceChildrenById");
}

template <typename CNPsType>
std::vector<int> ComponentTree<CNPsType>::countingSort(ImageUInt8Ptr imgPtr){
	int n = this->numRows * this->numCols;
    auto img = imgPtr->rawData();
	int maxvalue = img[0];
	for (int i = 1; i < n; i++)
		if(maxvalue < img[i]) maxvalue = img[i];
			
	std::vector<uint32_t> counter(maxvalue + 1, 0);  
	std::vector<int> orderedPixels(n);
    	
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

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
std::vector<int> ComponentTreeFZ::countingSort() {
    int numFZ = flatzoneGraph->getNumFlatZones();
    auto img = flatzoneGraph->getImage()->rawData();
    
    // Passo 1: encontrar o valor máximo de nível de cinza entre as flat zones
    uint8_t maxvalue = img[ pixelView.indexToPixel[0] ];
    for (int i = 1; i < numFZ; ++i) {
        int pixelID = pixelView.indexToPixel[i]; //flatzoneGraph->getPixelFromIndex(i);
        uint8_t val = img[pixelID];
        if (val > maxvalue) maxvalue = val;
    }
    std::vector<uint32_t> counter(maxvalue + 1);  
	std::vector<int> orderedFlatzones(numFZ);
    
    if (this->isMaxtree()) {
        // Passo 2: contar número de flat zones por nível de cinza
        for (int i = 0; i < numFZ; ++i) {
            int pixelID = pixelView.indexToPixel[i];//flatzoneGraph->getPixelFromIndex(i);
            uint8_t gray = img[pixelID];
            counter[gray]++;
        }

        // Passo 3: soma acumulada
        for (int i = 1; i < maxvalue; ++i)
            counter[i] += counter[i - 1];
        counter[maxvalue] += counter[maxvalue - 1];

        // Passo 4: ordenação dos representantes das flat zones
        for (int i = numFZ - 1; i >= 0; --i) {
            int pixelID = pixelView.indexToPixel[i];//flatzoneGraph->getPixelFromIndex(i);
            uint8_t gray = img[pixelID];
            orderedFlatzones[--counter[gray]] = pixelID;
        }

    } else { // min-tree
        // Passo 2: contar número de flat zones por (maxvalue - gray)
        for (int i = 0; i < numFZ; ++i) {
            int pixelID = pixelView.indexToPixel[i];//flatzoneGraph->getPixelFromIndex(i);
            uint8_t gray = img[pixelID];
            counter[maxvalue - gray]++;
        }

        // Passo 3: soma acumulada
        for (int i = 1; i < maxvalue; ++i)
            counter[i] += counter[i - 1];
        counter[maxvalue] += counter[maxvalue - 1];

        // Passo 4: ordenação dos representantes das flat zones
        for (int i = numFZ - 1; i >= 0; --i) {
            int pixelID = pixelView.indexToPixel[i];//flatzoneGraph->getPixelFromIndex(i);
            uint8_t gray = img[pixelID];
            orderedFlatzones[--counter[maxvalue - gray]] = pixelID;
        }
    }

    return orderedFlatzones;
}




template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::createTreeByUnionFind(std::vector<int>& orderedPixelFlatzones) {
    int numFZ = flatzoneGraph->getNumFlatZones();
    //parent.resize(numFZ);
    std::vector<int> zPar(numFZ, -1); 
    std::vector<int> parent(numFZ, -1); 
    auto findRoot = [&](int p) {
        while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
        return p;
    };

    //criando a arvore já canônizada
	for(int i=numFZ-1; i >= 0; i--){
		int p = orderedPixelFlatzones[i];
        int idxP = pixelView.pixelToIndex[p];
        zPar[idxP] = idxP;
        parent[idxP] = idxP; 
		for(int q: flatzoneGraph->getAdjacentFlatzonesFromPixel(p)){
            int idxQ = pixelView.pixelToIndex[q];
			if(zPar[idxQ] != -1){
				int idxR = findRoot(idxQ);
				if(idxP != idxR){
					parent[idxR] = idxP;
					zPar[idxR] = idxP;
				}
			}
		}
	}

    
    //int numFZ = flatzoneGraph->getNumFlatZones();
    auto img = flatzoneGraph->getImage()->rawData();
	for (int i = 0; i < numFZ; i++) {
		int p = orderedPixelFlatzones[i];
        int idxP = pixelView.pixelToIndex[p];// flatzoneGraph->getIndexFromPixel(p);
        int idxPParent = parent[idxP];
        int pParent = pixelView.indexToPixel[idxPParent];//flatzoneGraph->getPixelFromIndex(idxPParent);
        
		if (idxP == idxPParent) { //representante do node raiz
			int threshold1 = this->isMaxtree()? 0 : 255;
			int threshold2 = img[p];
            this->pixelToNodeId[p] = this->root = makeNode(p, -1, threshold1, threshold2);
		}
		else if (img[p] != img[pParent]) { //representante de um node
			int threshold1 = this->isMaxtree()? img[pParent]+1 : img[pParent]-1;
			int threshold2 = img[p];
            this->pixelToNodeId[p] = makeNode(p, this->pixelToNodeId[pParent], threshold1, threshold2);
		}
		else if (img[p] == img[pParent]) {
			this->pixelToNodeId[p] = pixelToNodeId[pParent];
		}else{
			std::cerr << "\n\n\n\nOps...falhou geral\n\n\n\n";
			break;
		}
	}
			
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, Pixels>::value, int>>
void ComponentTreeP::createTreeByUnionFind(std::vector<int>& orderedPixels, ImageUInt8Ptr imgPtr) {
	int numPixels = numRows*numCols;
    std::vector<int> zPar(numPixels, -1); 
    std::vector<int> parent(numPixels, -1); 
    auto findRoot = [&](int p) {
        while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
        return p;
    };
	auto img = imgPtr->rawData();
    
    //criando a arvore
	for(int i=numPixels-1; i >= 0; i--){
		int p = orderedPixels[i];
		parent[p] = p;
		zPar[p] = p;
		for (int q : this->adj->getNeighborPixels(p)) {
			if(zPar[q] != -1){
				int r = findRoot(q);
				if(p != r){
					parent[r] = p;
					zPar[r] = p;
				}
			}
		}
	}
    
	// canonizacao da arvore
    int numNodes = 0;
	for (int i = 0; i < numPixels; i++) {
		int p = orderedPixels[i];
		int q = parent[p];
				
		if(img[parent[q]] == img[q]){
			parent[p] = parent[q];
		}
        
        if (parent[p] == p || img[parent[p]] != img[p]) 
            ++numNodes;  // p é representante
    }

    reserveNodes(numNodes); 

    for (int i = 0; i < numPixels; i++) {
		int p = orderedPixels[i];
		//int q = parent[p];
		if (p == parent[p]) { //representante do node raiz
			int threshold1 = this->isMaxtree()? 0 : 255;
			int threshold2 = img[p];
            this->pixelToNodeId[p] = this->root = makeNode(p, -1, threshold1, threshold2);
		}
		else if (img[p] != img[parent[p]]) { //representante de um novo node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
            this->pixelToNodeId[p] = makeNode(p, this->pixelToNodeId[parent[p]], threshold1, threshold2);
		}
		else if (img[p] == img[parent[p]]) { 
			this->pixelToNodeId[p] = pixelToNodeId[parent[p]];
		}   
	}

    //colentando os pixels
    pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numNodes);
    pixelView = pixelBuffer->view();
    int indice = 0;
    for (int i = 0; i < numPixels; i++) {
		int p = orderedPixels[i];
        
        if (p == parent[p] || img[p] != img[parent[p]]){ //representante do set
            pixelView.indexToPixel[indice] = p;
            pixelView.pixelToIndex[p] = indice;
            pixelView.sizeSets[indice] = 1;
            pixelView.pixelsNext[p] = p;
            indice++;
        }
        else if (img[p] == img[parent[p]]) { 
            pixelView.pixelsNext[p] = pixelView.pixelsNext[parent[p]];
            pixelView.pixelsNext[parent[p]] = p;
            int idx = pixelView.pixelToIndex[ parent[p] ];
            pixelView.sizeSets[idx]++;
		}
    } 
      
    
    assert((indice == numNodes) && "Erro na contagem de sets");
    	
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::createTreeByUnionFind(std::vector<int>& orderedPixels, ImageUInt8Ptr imgPtr) {
    int numPixels = numRows * numCols;
    
    
    std::vector<int> zPar(numPixels, -1);
    std::vector<int> parent(numPixels, -1);
    auto findRoot = [&](int p) {
        while (zPar[p] != p) { zPar[p] = zPar[zPar[p]]; p = zPar[p]; }
        return p;
    };
    auto img = imgPtr->rawData();

    // --- construção da árvore (como já estava) ---
    for (int i = numPixels - 1; i >= 0; i--) {
        int p = orderedPixels[i];
        parent[p] = p;
        zPar[p] = p;
        for (int q : this->adj->getNeighborPixels(p)) {
            if (zPar[q] != -1) {
                int r = findRoot(q);
                if (p != r) {
                    parent[r] = p;
                    zPar[r] = p;
                }
            }
        }
    }
    
    int numNodes = 0;
    // --- canonização (como já estava) ---
    for (int i = 0; i < numPixels; i++) {
        int p = orderedPixels[i];
        int q = parent[p];
        if (img[parent[q]] == img[q]) {
            parent[p] = parent[q];
        }
        if(p == parent[p] || img[p] != img[parent[p]])
            ++numNodes;
    
    }

    reserveNodes(numNodes); 

    for (int i = 0; i < numPixels; i++) {
        int p = orderedPixels[i];
    
		if (p == parent[p]) { //representante do node raiz
			int threshold1 = this->isMaxtree()? 0 : 255;
			int threshold2 = img[p];
            this->pixelToNodeId[p] = this->root = makeNode(p, -1, threshold1, threshold2);
		}
		else if (img[p] != img[parent[p]]) { //representante de um novo node
			int threshold1 = this->isMaxtree()? img[parent[p]]+1 : img[parent[p]]-1;
			int threshold2 = img[p];
			this->pixelToNodeId[p] = makeNode(p, this->pixelToNodeId[parent[p]], threshold1, threshold2);
		}
		else if (img[p] == img[parent[p]]) { 
			this->pixelToNodeId[p] = pixelToNodeId[parent[p]];
		}
    }



    // --- aloca no teto: numFZs <= numPixels ---
    pixelBuffer = std::make_shared<PixelSetManager>(numPixels, numPixels);
    pixelView   = pixelBuffer->view();

    // Reuso de zPar como visited: marcar com -2 quando visitado
    auto isVisited  = [&](int p) { return zPar[p] == -2; };
    auto setVisited = [&](int p) { zPar[p] = -2; };
    FastQueue<int> queue(numPixels/3);
    int nextIdxFZ = 0; // quantas FZs realmente criadas

    for (int p0 = 0; p0 < numPixels; ++p0) {
        if (isVisited(p0)) continue;

        // platô de p0 (após canonização)
        int repNode = (parent[p0] == p0 || img[parent[p0]] != img[p0]) ? p0 : parent[p0];
        int L = img[p0];

        // head da FZ é p0 (o menor pixel da FZ, pela ordem do laço)
        int headFZ = p0;

        // cria idx para esta FZ
        int idx = nextIdxFZ++;
        pixelView.indexToPixel[idx] = headFZ;  // head
        pixelView.sizeSets[idx]     = 0;       // vamos contar dentro da BFS
        pixelView.pixelsNext[headFZ] = headFZ; // ciclo unitário inicial (será expandido no loop)

        // BFS da FZ
        setVisited(p0);
        queue.push(p0);

        while (!queue.empty()) {
            int p = queue.pop();

            // insere p após o head
            pixelView.pixelsNext[p]      = pixelView.pixelsNext[headFZ];
            pixelView.pixelsNext[headFZ] = p;
            pixelView.pixelToIndex[p]    = idx;
            pixelView.sizeSets[idx]++;

            for (int q : this->adj->getNeighborPixels(p)) {
                if (!isVisited(q) && img[q] == L) {
                    int rq = (parent[q] == q || img[parent[q]] != img[q]) ? q : parent[q];
                    if (rq == repNode) {
                        setVisited(q);
                        queue.push(q);
                    }
                }
            }
        }

        // registra a FZ no nó correspondente ao representante do platô
        arena.repCNPs[this->pixelToNodeId[repNode]].push_back(headFZ);
    }

    // encolhe vetores de sets para o número real de FZs encontradas
    pixelBuffer->shrinkToNumSets(nextIdxFZ);
	pixelView = pixelBuffer->view();  // reobter spans após o resize
}


template <typename CNPsType>
ComponentTree<CNPsType>::ComponentTree(ImageUInt8Ptr img, bool isMaxtree, AdjacencyRelationPtr adj) : numRows(img->getNumRows()), numCols(img->getNumCols()), maxtreeTreeType(isMaxtree), adj(adj), numNodes(0){   
    this->pixelToNodeId.resize(numRows * numCols, -1);
    build(img);
    
}

template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
ComponentTreeFZ::ComponentTree(std::shared_ptr<FlatZonesGraph> graph, bool isMaxtree) :  maxtreeTreeType(isMaxtree), numNodes(0), flatzoneGraph(graph){   
    ImageUInt8Ptr img = flatzoneGraph->getImage();
    this->numRows = img->getNumRows();
    this->numCols = img->getNumCols();
    this->adj = flatzoneGraph->getAdjacencyRelation();
    this->pixelToNodeId.resize(numRows * numCols, -1);
    this->pixelBuffer = flatzoneGraph->getPixelSetManager();
    this->pixelView = pixelBuffer->view();
    build();    
}

    
template <typename CNPsType>
void ComponentTree<CNPsType>::build(){ 
    
    if(flatzoneGraph)
        reserveNodes(flatzoneGraph->getNumFlatZones());

    std::vector<int> orderedFlatzones = countingSort();
	createTreeByUnionFind(orderedFlatzones);
    
    //atribui os representantes das FZ nos respectivos nodes
    for(int rep : pixelBuffer->getFlatzoneRepresentatives()) {
        const NodeId nid = this->pixelToNodeId[rep];
        arena.repCNPs[nid].push_back(rep);
    }

	computerArea(this->root); //computer area
    
}


template <typename CNPsType>
void ComponentTree<CNPsType>::build(ImageUInt8Ptr imgPtr){ 
	std::vector<int> orderedPixels = countingSort(imgPtr);
    createTreeByUnionFind(orderedPixels, imgPtr);
	computerArea(this->root); //computer area
}




template <typename CNPsType>
void ComponentTree<CNPsType>::computerArea(NodeId id){
    int32_t area = getNumCNPsById(id);
    for (NodeId c : arena.children(id)) {
        computerArea(c);
        area += arena.areaCC[c]; // já calculado na recursão
    }
    arena.areaCC[id] = area;
}

template<>
template<typename T, typename std::enable_if_t<std::is_same<T, Pixels>::value, int>>
void ComponentTreeP::prunning(NodeId nodeId){
    assert(nodeId && "node is invalid");
    assert(getParentById(nodeId) && "node is root");

    const int parentId = arena.parentId[nodeId]; 
    if(parentId >= 0){ 
        // 1) desconecta 'node' do pai
        removeChildById(parentId, nodeId, false);

        // 2) BFS ID-first na subárvore para redirecionar UF/SC e contabilizar remoções
        FastQueue<NodeId> q;
        q.push(nodeId);

        const int parentRep = arena.repNode[parentId];
        while (!q.empty()) {
            NodeId curId = q.pop();

            // enfileira filhos por IDs
            for (NodeId c : arena.children(curId)) {
                q.push(c);
            }

            // une representantes no UF e atualiza pixel->node para o representante
            int repCur = arena.repNode[curId];
            setSCById(repCur, parentId);
            pixelBuffer->mergeSetsByRep(parentRep, repCur);

            // release no node removido
            releaseNode(curId);
        }
    }
}

template<>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
void ComponentTreeFZ::prunning(NodeId nodeId) {
    assert(nodeId && "node is invalid");
    assert(getParentById(nodeId) && "node is root");

    int parentId   = arena.parentId[nodeId];
    if(parentId >=0 ){

        removeChildById(parentId, nodeId, false);

        FastQueue<NodeId> q;
        q.push(nodeId);
        while (!q.empty()) {
            NodeId curId = q.pop();

            // 1) transfere CNPs (todas FZs representadas nesse nó) para o pai
            for (int rep : arena.repCNPs[curId]) {
                arena.repCNPs[parentId].push_back(rep);
                setSCById(rep, parentId);
            }
            // 2) une representante de nó no UF para manter coerência com pixel->node de reps
            //parent[ arena.repNode[curId] ] = arena.repNode[parentId];
            setSCById(arena.repNode[curId], parentId);

            // 3) fila dos filhos por IDs
            for (NodeId c : arena.children(curId)) {
                q.push(c);
            }

            releaseNode(curId);
        }
    }
}




template<>
template<typename T, typename std::enable_if_t<std::is_same<T, Pixels>::value, int>>
void ComponentTreeP::mergeWithParent(NodeP node)
{
    if (!node || !node.getParent()) return;

    const int nodeId   = node.getIndex();
    const int parentId = arena.parentId[nodeId];

    // 1) tira 'node' da lista de filhos do pai
    removeChildById(parentId, nodeId, false);

    // 2) move os filhos de 'node' para o pai (preserva ordem) — tudo por IDs
    spliceChildrenById(parentId, nodeId);

    // 3) une representantes e atualiza pixel->node para o representante do nó colapsado
    const int parentRep = arena.repNode[parentId];
    const int nodeRep   = arena.repNode[nodeId];
    setSCById(nodeRep, parentId);
    pixelBuffer->mergeSetsByRep(parentRep, nodeRep);

    // 4) marca o nó como desconectado
    releaseNode(nodeId);
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, Pixels>::value, int>>
void ComponentTreeP::mergeWithParent(std::vector<int>& flatzone){
    int idFlatzone = flatzone.front();
    NodeP node = proxy(this->pixelToNodeId[idFlatzone]);
    if(getNumCNPsById(node) == static_cast<int>(flatzone.size())) {
        this->mergeWithParent(node);
    }
    else{
        //TODO: pensar em como otimizar esse caso
        //winner ganha os pixels de flatzone e o loser perde esses pixels
/*
        NodeP parent = node.getParent();
        int repWinner = parent.getRepNode(); //representante do pai
        int repLoser  = node.getRepNode();   //representante do filho

        //1. Recupera índices dos representantes
        int idxRootWinner = pixelView.pixelToIndex[repWinner];
        int idxRootLoser  = pixelView.pixelToIndex[repLoser];

        //2. Atualiza a quantidade de cnps
        pixelView.sizeSets[idxRootWinner] += flatzone->size();
        pixelView.sizeSets[idxRootLoser] -= flatzone->size();

        for( int p: *flatzone) {				
            parent.addRepCNPs(p);
            this->pixelToNode[p] = parent;	


            // 3. Splice O(1) das listas circulares (pixels)
            int nextWinner = pixelView.pixelsNext[repWinner];
            int nextLoser  = pixelView.pixelsNext[repLoser];
            pixelView.pixelsNext[repWinner] = nextLoser;
            pixelView.pixelsNext[repLoser]  = nextWinner;

            // 4. Invalida slot perdedor
            pixelView.sizeSets[idxRootLoser]  = 0;
            pixelView.indexToPixel[idxRootLoser] = -1;

            // 5. Redireciona lookups pelo antigo rep pixel
            pixelView.pixelToIndex[repLoser] = idxRootWinner;
        }
        */
    }
}




template <typename CNPsType>
std::vector<NodeId> ComponentTree<CNPsType>::getLeaves(){
    std::vector<NodeId> leaves;
    FastQueue<NodeId> s;
    s.push(this->root);

    while (!s.empty()) {
        NodeId id = s.pop();
        if (arena.childCount[id] == 0) {
            leaves.push_back(id);
        } else {
            for(NodeId c: arena.children(id)){
                s.push(c);
            }
        }
    }
    return leaves;
}


template <>
template<typename T, typename std::enable_if_t<std::is_same<T, FlatZones>::value, int>>
std::shared_ptr<FlatZonesGraph>& ComponentTreeFZ::getFlatZonesGraph(){
    return this->flatzoneGraph;
}


template <typename CNPsType>
ImageUInt8Ptr ComponentTree<CNPsType>::reconstructionImage(){
    ImageUInt8Ptr imgPtr = ImageUInt8::create(this->numRows, this->numCols);
    this->reconstruction(this->root, imgPtr->rawData());
    return imgPtr;
}


// FlatZones
template <>
inline void ComponentTreeFZ::reconstruction(NodeId id, uint8_t* imgOut) {
    assert(imgOut && "Erro: Ponteiro de saída da imagem é nulo!");
    for (int p : pixelBuffer->getPixelsBySet(arena.repCNPs[id])) {
        imgOut[p] = static_cast<uint8_t>(arena.threshold2[id]);
    }
    for (int c : arena.children(id)) {
        reconstruction(c, imgOut);
    }
}

// Pixels
template <>
inline void ComponentTreeP::reconstruction(NodeId id, uint8_t* imgOut) {
    assert(imgOut && "Erro: Ponteiro de saída da imagem é nulo!");
    for (int p : pixelBuffer->getPixelsBySet(arena.repNode[id])) {
        imgOut[p] = static_cast<uint8_t>(arena.threshold2[id]);
    }
    for (int c : arena.children(id)) {
        reconstruction(c, imgOut);
    }
}



// FlatZones
template <>
inline auto ComponentTreeFZ::getCNPsById(NodeId id) const {
    return pixelBuffer->getPixelsBySet(arena.repCNPs[id]);
}

// Pixels
template <>
inline auto ComponentTreeP::getCNPsById(NodeId id) const {
    return pixelBuffer->getPixelsBySet(arena.repNode[id]);
}

// FlatZones
template <>
inline int ComponentTreeFZ::getNumFlatzoneById(NodeId id) const {
    return static_cast<int>(arena.repCNPs[id].size());
}

// Pixels
template <>
inline int ComponentTreeP::getNumFlatzoneById(NodeId id) const {
    (void)id;
    return 1;
}

// FlatZones
template <>
inline int ComponentTreeFZ::getNumCNPsById(NodeId id) const {
    return this->pixelBuffer->numPixelsInSets(arena.repCNPs[id]);
}

// Pixels
template <>
inline int ComponentTreeP::getNumCNPsById(NodeId id) const {
    return this->pixelBuffer->numPixelsInSet(arena.repNode[id]);
}





