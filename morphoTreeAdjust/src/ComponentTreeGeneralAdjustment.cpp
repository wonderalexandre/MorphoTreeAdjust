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
	
	//lines between 7 and 12
    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> neighborNodesSet1; //line 2: open
	std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> neighborNodesSet2;  //line3: close
	for (int p : flatZoneP) {
	    for (int q : adj->getAdjPixels(p)) {
			if(maxtree.getSC(p)->getLevel() < maxtree.getSC(q)->getLevel()){  //elements of {q_i}
				if (maxtree.getSC(q)->getThreshold2() <= b) { 
					neighborNodesSet2.insert( maxtree.getSC(q) );
				}else{
					neighborNodesSet1.insert( maxtree.getSC(q) );
				}
			}
	    }
	}
	auto comparator = [](NodeCT* a, NodeCT* b) { return a->getLevel() > b->getLevel(); };
    std::priority_queue<NodeCT*, std::vector<NodeCT*>, decltype(comparator)> neighborNodes1(comparator);
	for (NodeCT* node : neighborNodesSet1) {
        neighborNodes1.push(node);
    }

	std::priority_queue<NodeCT*, std::vector<NodeCT*>, decltype(comparator)> neighborNodes2(comparator);
	for (NodeCT* node : neighborNodesSet2) {
        neighborNodes2.push(node);
    }

	int w = -1; //line 13

	//lines between 14 and 35
	while (!neighborNodes1.empty()) {
		NodeCT* node = neighborNodes1.top();
		neighborNodes1.pop();
		std::cout << "neighborNodes1: index: " << node->getIndex() << ", level: " << node->getLevel() << std::endl;
	}

	std::cout << "\n\n";

	//lines between 14 and 35
	while (!neighborNodes1.empty()) {
		NodeCT* node = neighborNodes1.top();
		neighborNodes1.pop();
		std::cout << "neighborNodes2: index: " << node->getIndex() << ", level: " << node->getLevel() << std::endl;
	}

}
