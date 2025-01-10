#include "../include/ComponentTreeGeneralAdjustment.hpp"
#include "../include/NodeCT.hpp"
#include <unordered_set>
#include <list>
#include <iostream>

#include <queue>
#include <vector>
#include <set>
#include <functional> 
#include <algorithm>



void ComponentTreeGeneralAdjustment::adjustMaxTree(ComponentTree &maxtree, std::list<int> flatZoneP, int newGrayLevel) {
	
	AdjacencyRelation* adj = maxtree.getAdjacencyRelation();
	NodeCT* nodeContainsP = maxtree.getSC( flatZoneP.front() ); //Nó contendo Xa
	
	int a = nodeContainsP->getThreshold2();
	int b = newGrayLevel;
	
	NodeCT* nodeP = new NodeCT(); //lines 1 and 2
	nodeP->setCNPs(flatZoneP); //line 3
	nodeP->setLevel(b); //line 5
	nodeP->setParent( nodeContainsP ); //line 6
	
	nodeContainsP->removeCNPs(flatZoneP); //line 4: ρ(Xa) := ρ(Xa) \ P

	int v = -1; //line 7 (inicializacao)

	//lines between 1 and 6
	std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> neighborNodesAdded; // Conjunto para evitar duplicatas nos nós adicionados
	auto comparator = [](NodeCT* a, NodeCT* b) { return a->getLevel() > b->getLevel(); }; // Definição de comparador para ordenar a fila de prioridade pelo nível (getLevel)
	std::priority_queue<NodeCT*, std::vector<NodeCT*>, decltype(comparator)> neighborNodesOpen(comparator);
	std::priority_queue<NodeCT*, std::vector<NodeCT*>, decltype(comparator)> neighborNodesClose(comparator);
	for (int p : flatZoneP) { 
		for (int q : adj->getAdjPixels(p)) {
			if (maxtree.getSC(p)->getLevel() < maxtree.getSC(q)->getLevel()) { // Isso garante que estamos lidando com {q_i} corretamente
				if (neighborNodesAdded.find(maxtree.getSC(q)) == neighborNodesAdded.end()) { // Verificar se o nó contendo q já foi adicionado a fila 
					neighborNodesAdded.insert(maxtree.getSC(q)); // Adicionar o nó q ao conjunto para evitar duplicatas
					if (maxtree.getSC(q)->getThreshold2() <= b) {
						neighborNodesClose.push(maxtree.getSC(q));
					} else {
						neighborNodesOpen.push(maxtree.getSC(q));
						if(v < maxtree.getSC(q)->getThreshold1()){ //construindo v (line 7)
							v = maxtree.getSC(q)->getThreshold1(); 
						}
					}
				}
			}
		}
	}

	int w = -1; //line 13

	//lines between 14 and 35
	while (!neighborNodesOpen.empty()) {
		NodeCT* node = neighborNodesOpen.top();
		neighborNodesOpen.pop();
		std::cout << "neighborNodesOpen: index: " << node->getIndex() << ", level: " << node->getLevel() << std::endl;
	}

	std::cout << "\n\n";

	//lines between 14 and 35
	while (!neighborNodesClose.empty()) {
		NodeCT* node = neighborNodesClose.top();
		neighborNodesClose.pop();
		std::cout << "neighborNodesClose: index: " << node->getIndex() << ", level: " << node->getLevel() << std::endl;
	}

}
