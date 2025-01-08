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
    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> neighborNodesSet; 
	neighborNodesSet.insert(nodeP);
	for (int p : flatZoneP) {
	    for (int q : adj->getAdjPixels(p)) {
			if (a < maxtree.getSC(q)->getThreshold2() && b >= maxtree.getSC(q)->getThreshold2()) { 
				neighborNodesSet.insert( maxtree.getSC(q) );
		    }
	    }
	}
	auto comparator = [](NodeCT* a, NodeCT* b) { return a->getLevel() > b->getLevel(); };
    std::priority_queue<NodeCT*, std::vector<NodeCT*>, decltype(comparator)> neighborNodes(comparator);
	for (NodeCT* node : neighborNodesSet) {
        neighborNodes.push(node);
    }


	int w = -1; //line 13

	//lines between 14 and 35
	while (!neighborNodes.empty()) {
		NodeCT* node = neighborNodes.top();
		neighborNodes.pop();
		std::cout << "neighborNodes: index: " << node->getIndex() << ", level: " << node->getLevel() << std::endl;
	}

}
