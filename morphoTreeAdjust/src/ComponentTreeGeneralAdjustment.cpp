#include "../include/ComponentTreeGeneralAdjustment.hpp"
#include "../include/NodeCT.hpp"
#include <unordered_set>
#include <list>
#include <iostream>
 



void ComponentTreeGeneralAdjustment::adjustMaxTree(ComponentTree &maxtree, std::list<int> flatZoneP, int newGrayLevel) {
	
	AdjacencyRelation* adj = maxtree.getAdjacencyRelation();
	NodeCT* nodeContainsP = maxtree.getSC( flatZoneP.front() ); //Nó contendo Xa
	
	int a = nodeContainsP->getLevel();
	int b = newGrayLevel;
	
	NodeCT* nodeP = new NodeCT(); //lines 1 and 2
	nodeP->setCNPs(flatZoneP); //line 3
	nodeP->setLevel(b); //line 5
	nodeP->setParent( nodeContainsP ); //line 6
	
	nodeContainsP->removeCNPs(flatZoneP); //line 4: ρ(Xa) := ρ(Xa) \ P
	
	//lines between 7 and 12
    std::unordered_set<NodeCT*, NodeCT::NodeHashFunction> neighborNodes; 
	neighborNodes.insert(nodeP );
	for (int p : flatZoneP) {
	    for (int q : adj->getAdjPixels(p)) {
			if (a > maxtree.getSC(q)->getLevel() && b <= maxtree.getSC(q)->getLevel()) { 
				neighborNodes.insert( maxtree.getSC(q) );
		    }
	    }
	}
	

}
