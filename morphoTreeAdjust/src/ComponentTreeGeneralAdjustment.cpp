#include "../include/ComponentTreeGeneralAdjustment.hpp"
#include "../include/NodeCT.hpp"
#include <unordered_set>
#include <list>
#include <iostream>
 



void ComponentTreeGeneralAdjustment::adjustMaxTree(ComponentTree &maxtree, std::list<int> flatZoneP, int newGrayLevel) {
	
	AdjacencyRelation* adj = maxtree.getAdjacencyRelation();
	NodeCT* nodeContainsP = maxtree.getSC( flatZoneP.front() );
	
	int a = nodeContainsP->getLevel();
	int b = newGrayLevel;
	
	NodeCT* nodeP = new NodeCT();
	nodeP->setLevel(b);
	nodeP->setCNPs(flatZoneP);
	nodeP->setParent( nodeContainsP );
	
    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> neighborNodes;
	for (int p : flatZoneP) {
	    for (int q : adj->getAdjPixels(p)) {
			if (a > maxtree.getSC(q)->getLevel() && b <= maxtree.getSC(q)->getLevel()) { 
				neighborNodes.insert( maxtree.getSC(q) );
		    }
	    }
	}
	

}
