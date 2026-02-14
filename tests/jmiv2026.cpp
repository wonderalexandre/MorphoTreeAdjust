
#include <iostream>
#include <list>
#include <chrono>
#include <filesystem>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentBySubtree.hpp"
#include "../morphoTreeAdjust/include/Common.hpp"
#include "../morphoTreeAdjust/include/FlatZonesGraph.hpp"

#include "../morphoTreeAdjust/include/AttributeComputer.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"
#include "../tests/Tests.hpp"

namespace fs = std::filesystem;
const bool ENABLE_PRINT = false;
constexpr bool ENABLE_LOG = false;

static inline long long elapsedMs(const Stopwatch& sw) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(sw.elapsed()).count();
}

enum class GraphKind {
    OnDemand,
    Eager,
    Naive
};

static const char* graphKindToString(GraphKind kind) {
    switch (kind) {
        case GraphKind::OnDemand: return "on_demand";
        case GraphKind::Eager: return "eager";
        case GraphKind::Naive: return "naive";
    }
    return "on_demand";
}

static bool parseGraphKind(const std::string& value, GraphKind* out) {
    if (value == "on_demand" || value == "ondemand" || value == "default") {
        *out = GraphKind::OnDemand;
        return true;
    }
    if (value == "eager") {
        *out = GraphKind::Eager;
        return true;
    }
    if (value == "naive") {
        *out = GraphKind::Naive;
        return true;
    }
    return false;
}

static void printUsage(const char* argv0) {
    std::cerr << "Use: " << argv0 << " <image> [graph]\n"
              << "  graph: on_demand|eager|naive (default: on_demand)\n";
}


struct PauseTimers {
    Stopwatch* sw;
    Stopwatch* swAll;
    bool swWasRunning;
    bool swAllWasRunning;

    PauseTimers(Stopwatch* sw, Stopwatch* swAll)
        : sw(sw),
          swAll(swAll),
          swWasRunning(sw && sw->running()),
          swAllWasRunning(swAll && swAll->running()) {
        if (swWasRunning) {
            sw->pause();
        }
        if (swAllWasRunning) {
            swAll->pause();
        }
    }

    ~PauseTimers() {
        if (swWasRunning) {
            sw->resume();
        }
        if (swAllWasRunning) {
            swAll->resume();
        }
    }
};

template <typename CNPsType, typename GraphT>
std::vector<NodeId> getNodesThreshold(ComponentTree<CNPsType, GraphT>* tree,
                                      std::span<float> attribute,
                                      int threshold,
                                      Stopwatch* sw = nullptr,
                                      Stopwatch* swAll = nullptr) {
    std::vector<NodeId> lista;
    FastQueue<NodeId> queue;
    queue.push(tree->getRootById());

    int numFlatZones = 0;
    int numNodes = 0;
    while (!queue.empty()) {
        NodeId id = queue.pop();
        if (attribute[id] > threshold) {
            for(NodeId c: tree->getChildrenById(id)) { 
                queue.push(c);
            }
        } else {
            if (ENABLE_LOG) {
                PauseTimers pause(sw, swAll);
                numFlatZones += (tree->computerNumFlatzoneDescendants(id) + tree->getNumFlatzoneById(id));
                numNodes     += tree->computerNumDescendants(id) + 1;
            }
            lista.push_back(id);
        }
    }
    if (ENABLE_LOG) {
        PauseTimers pause(sw, swAll);
        std::cout << "\tThreshold: " << threshold
                << ", #Nodes: " << numNodes
                << ", #FlatZones: " << numFlatZones
                << ", #InputTreeNodes: " << tree->getNumNodes()
                << std::endl;
    }
    return lista;
}

ImageUInt8Ptr computerCASF_naive(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    Stopwatch sw;
    Stopwatch swAll;
    swAll.start();
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            sw.start();
            swAll.start();
        }
        ComponentTreePPtr maxtreePtr = std::make_shared<ComponentTreeP>(imgOut, true, adj);
        
        ComponentTreeP* maxtree = maxtreePtr.get();
        
        AreaComputerP computerAttrMax( maxtree );
        std::vector<float> attribute = computerAttrMax.compute();
        
        for(NodeId node: getNodesThreshold(maxtree, attribute, threshold, &sw, &swAll)) {
	        maxtree->prunning(node);    
	    }
        
        imgOut = maxtree->reconstructionImage();
        if(ENABLE_PRINT){
            sw.pause();
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (build/prunning/rec maxtree) naive: " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }

        ComponentTreePPtr mintreePtr = std::make_shared<ComponentTreeP>(imgOut, false, adj);
        ComponentTreeP* mintree = mintreePtr.get();
        AreaComputerP computerAttrMin ( mintree );
        attribute = computerAttrMin.compute();
        
	    for(NodeId node: getNodesThreshold(mintree, attribute, threshold, &sw, &swAll)) {
	        mintree->prunning(node);    
	    }

        imgOut = mintree->reconstructionImage();  
        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (build/prunning/rec mintree) naive: " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }
  	    
	}
    
    return imgOut;
}


template <typename GraphT = DefaultFlatZonesGraph>
ImageUInt8Ptr computerCASF_hybrid(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds, size_t cutoffPointHybrid){
    Stopwatch sw;
    Stopwatch swAll;
    if(ENABLE_PRINT){
        sw.start();
    }

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    for(size_t i=0; i < cutoffPointHybrid; i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            sw.start();
            swAll.start();
        }
        
        ComponentTreePPtr maxtreePtr = std::make_shared<ComponentTreeP>(imgOut, true, adj);
        ComponentTreeP* maxtree = maxtreePtr.get();
        AreaComputerP computerAttrMax ( maxtree );
        std::vector<float> attribute = computerAttrMax.compute();
        
	    for(NodeId node: getNodesThreshold(maxtree, attribute, threshold, &sw, &swAll)) {
	        maxtree->prunning(node);    
	    }
        
        imgOut = maxtree->reconstructionImage();
        if(ENABLE_PRINT){
            sw.pause();
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (update mintree and pruning maxtree) hybrid: " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }
	    ComponentTreePPtr mintreePtr = std::make_shared<ComponentTreeP>(imgOut, false, adj);
        ComponentTreeP* mintree = mintreePtr.get();

        AreaComputerP computerAttrMin ( mintree );
        attribute = computerAttrMin.compute();
        
	    for(NodeId node: getNodesThreshold(mintree, attribute, threshold, &sw, &swAll)) {
	        mintree->prunning(node);    
	    }

        imgOut = mintree->reconstructionImage();  
        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (update maxtree and pruning mintree) hybrid: " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }

    }

    sw.start();
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(imgOut, adj);
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    ComponentTreeAdjustmentByLeaf<DefaultAttributeComputerT<GraphT>, GraphT> adjust(mintree, maxtree);
    AreaComputerFZT<GraphT> computerAttrMax ( maxtree );
    AreaComputerFZT<GraphT> computerAttrMin ( mintree );
    std::vector<float> attributeMax = computerAttrMax.compute();
    std::vector<float> attributeMin = computerAttrMin.compute();
    adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

    if(ENABLE_PRINT){
        sw.pause();
       // std::cout << "\tTime (build trees): " << elapsedMs(sw) << " ms\n";
    }
    
    for(size_t i=cutoffPointHybrid; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            swAll.start();
            sw.start();
        }
        auto nodesToPruning = getNodesThreshold(maxtree, attributeMax, threshold, &sw, &swAll);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        
        if(ENABLE_PRINT){
            sw.pause();
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (update mintree and pruning maxtree) hybrid: " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }
        nodesToPruning = getNodesThreshold(mintree, attributeMin, threshold, &sw, &swAll);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning); 

        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (update maxtree and pruning mintree) hybrid: " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }

	}
    imgOut = mintree->reconstructionImage();
    return imgOut;
}


template <typename GraphT = DefaultFlatZonesGraph>
ImageUInt8Ptr computerCASF(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    Stopwatch sw;
    Stopwatch swAll;
    if(ENABLE_PRINT){
        sw.start();
    }

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(img, adj);
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();
    
    ComponentTreeAdjustmentByLeaf<DefaultAttributeComputerT<GraphT>, GraphT> adjust(mintree, maxtree);
    AreaComputerFZT<GraphT> computerAttrMax ( maxtree );
    AreaComputerFZT<GraphT> computerAttrMin ( mintree );
    std::vector<float> attributeMax = computerAttrMax.compute();
    std::vector<float> attributeMin = computerAttrMin.compute();
    adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

    if(ENABLE_PRINT){
        sw.pause();
        {
            PauseTimers pause(nullptr, &swAll);
            std::cout << "\n\tTime (build trees): " << elapsedMs(sw) << " ms\n";
        }
    }
    
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            swAll.start();
            sw.start();
        }
        auto nodesToPruning = getNodesThreshold(maxtree, attributeMax, threshold, &sw, &swAll);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        sw.pause();

        if(ENABLE_PRINT){
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (update mintree and pruning maxtree): " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }
        nodesToPruning = getNodesThreshold(mintree, attributeMin, threshold, &sw, &swAll);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning); 
        
        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (update maxtree and pruning mintree): " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }
	}
    auto imgOut = mintree->reconstructionImage();
    return imgOut;
}


template <typename GraphT = DefaultFlatZonesGraph>
ImageUInt8Ptr computerCASF_hybridSubtree(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds, size_t cutoffPointHybrid){
    Stopwatch sw;
    Stopwatch swAll;
    if(ENABLE_PRINT){
        sw.start();
    }
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();

    for(size_t i=0; i < cutoffPointHybrid; i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            sw.start();
            swAll.start();
        }
        
        ComponentTreePPtr maxtreePtr = std::make_shared<ComponentTreeP>(imgOut, true, adj);
        ComponentTreeP* maxtree = maxtreePtr.get();
        AreaComputerP computerAttrMax ( maxtree );
        std::vector<float> attributeMax = computerAttrMax.compute();
	    for(NodeId node: getNodesThreshold(maxtree, attributeMax, threshold, &sw, &swAll)) {
	        maxtree->prunning(node);
	    }
        
        imgOut = maxtree->reconstructionImage();
        if(ENABLE_PRINT){
            sw.pause();
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (update mintree and pruning maxtree) hybrid: " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }
	    ComponentTreePPtr mintreePtr = std::make_shared<ComponentTreeP>(imgOut, false, adj);
        ComponentTreeP* mintree = mintreePtr.get();

        AreaComputerP computerAttrMin ( mintree );
        std::vector<float> attributeMin = computerAttrMin.compute();
	    for(NodeId node: getNodesThreshold(mintree, attributeMin, threshold, &sw, &swAll)) {
	        mintree->prunning(node);    
	    }

        imgOut = mintree->reconstructionImage();  
        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (update maxtree and pruning mintree) hybrid: " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }

    }

    sw.start();
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(imgOut, adj);

    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<GraphT>, GraphT> adjust(mintree, maxtree);
    AreaComputerFZT<GraphT> computerAttrMax ( maxtree );
    AreaComputerFZT<GraphT> computerAttrMin ( mintree );
    std::vector<float> attributeMax = computerAttrMax.compute();
    std::vector<float> attributeMin = computerAttrMin.compute();
    adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);
    if(ENABLE_PRINT){
        sw.pause();
        //std::cout << "\tTime (build trees): " << elapsedMs(sw) << " ms\n";
    }
    
    for(size_t i=cutoffPointHybrid; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            swAll.start();
            sw.start();
        }
        auto nodesToPruning = getNodesThreshold(maxtree, attributeMax, threshold, &sw, &swAll);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        
        if(ENABLE_PRINT){
            sw.pause();
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (update mintree and pruning maxtree) hybrid: " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }
        nodesToPruning = getNodesThreshold(mintree, attributeMin, threshold, &sw, &swAll);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning); 

        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (update maxtree and pruning mintree) hybrid: " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }

	}
    imgOut = mintree->reconstructionImage();
    return imgOut;
}



template <typename GraphT = DefaultFlatZonesGraph>
ImageUInt8Ptr computerCASF_subtree(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    Stopwatch sw;
    Stopwatch swAll;
    if(ENABLE_PRINT){
        sw.start();
    }
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    sw.start();
    std::shared_ptr<GraphT> graph = std::make_shared<GraphT>(img, adj);
    ComponentTreeFZPtr<GraphT> maxtreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, true);
    ComponentTreeFZPtr<GraphT> mintreePtr = std::make_shared<ComponentTreeFZ<GraphT>>(graph, false);
    
    ComponentTreeFZ<GraphT>* maxtree = maxtreePtr.get();
    ComponentTreeFZ<GraphT>* mintree = mintreePtr.get();

    ComponentTreeAdjustmentBySubtree<DefaultAttributeComputerT<GraphT>, GraphT> adjust(mintree, maxtree);
    AreaComputerFZT<GraphT> computerAttrMax ( maxtree );
    AreaComputerFZT<GraphT> computerAttrMin ( mintree );
    std::vector<float> attributeMax = computerAttrMax.compute();
    std::vector<float> attributeMin = computerAttrMin.compute();
    adjust.setAttributeComputer(computerAttrMin, computerAttrMax, attributeMin, attributeMax);

    if(ENABLE_PRINT){
        sw.pause();
        {
            PauseTimers pause(nullptr, &swAll);
            std::cout << "\n\tTime (build trees): " << elapsedMs(sw) << " ms\n";
        }
    }
    
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            swAll.start();
            sw.start();
        }
        auto nodesToPruning = getNodesThreshold(maxtree, attributeMax, threshold, &sw, &swAll);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        sw.pause();
        
        if(ENABLE_PRINT){
            {
                PauseTimers pause(nullptr, &swAll);
                std::cout << "\t- Time (update mintree and pruning maxtree): " << elapsedMs(sw) << " ms\n";
            }
            sw.start();
        }
        nodesToPruning = getNodesThreshold(mintree, attributeMin, threshold, &sw, &swAll);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning); 
        if(ENABLE_PRINT){
            sw.pause();
            swAll.pause();
            std::cout << "\t- Time (update maxtree and pruning mintree): " << elapsedMs(sw) << " ms\n";
            std::cout << "Time: " << elapsedMs(swAll) << " ms\n\n";
        }

	}
    auto imgOut = mintree->reconstructionImage();
    return imgOut;
}


int main(int argc, char* argv[]) {
    std::string filename;
    GraphKind graphKind = GraphKind::Eager;
    if (argc < 2 || argc > 3) {
        printUsage(argv[0]);
        return 1;
    }
    filename = argv[1];
    if (argc == 3) {
        if (!parseGraphKind(argv[2], &graphKind)) {
            std::cerr << "Erro: tipo de grafo invalido: " << argv[2] << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }
    
    std::cout << "Image: " << filename << std::endl;
    std::cout << "Graph: " << graphKindToString(graphKind) << std::endl;
    int numCols, numRows, nchannels;
    uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    
    if (!data) {
        std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
        return 1;
    }   
    ImageUInt8Ptr img = ImageUInt8::fromRaw(data, numRows, numCols);
    std::cout << "Resolution: " << numCols << "x" << numRows << std::endl;
    
    //std::vector<int> thresholds = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
    //std::vector<int> thresholds = {111, 222, 333, 445, 556, 667, 778, 889, 1000, 1111, 1222, 1333, 1445, 1556, 1667, 1778, 1889, 2000, 2111, 2222, 2333, 2445, 2556, 2667, 2778, 2889, 3000, 3111, 3222, 3333, 3445, 3556, 3667, 3778, 3889, 4000, 4111, 4222, 4333, 4445, 4556, 4667, 4778, 4889, 5000, 5111, 5222, 5333, 5445, 5556, 5667};
    std::vector<int> thresholds = {50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500};

    

    double radioAdj = 1.5;
    
    Stopwatch sw;
    sw.start();
    auto imgOut1 = computerCASF_naive(img, radioAdj, thresholds);
    sw.pause();
    std::cout << "Time naive approach: " << elapsedMs(sw) << " ms\n";

    ImageUInt8Ptr imgOut2;
    ImageUInt8Ptr imgOut3;
    ImageUInt8Ptr imgOut4;
    ImageUInt8Ptr imgOut5;
    switch (graphKind) {
        case GraphKind::OnDemand:
            sw.start();
            imgOut2 = computerCASF<DefaultFlatZonesGraph>(img, radioAdj, thresholds);
            sw.pause();
            std::cout << "Time our approach:   " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut3 = computerCASF_hybrid<DefaultFlatZonesGraph>(img, radioAdj, thresholds, 2);
            sw.pause();
            std::cout << "Time our approach_hybrid: " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut4 = computerCASF_subtree<DefaultFlatZonesGraph>(img, radioAdj, thresholds);
            sw.pause();
            std::cout << "Time our approach (subtree): " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut5 = computerCASF_hybridSubtree<DefaultFlatZonesGraph>(img, radioAdj, thresholds, 1);
            sw.pause();
            std::cout << "Time our approach_hybrid (subtree): " << elapsedMs(sw) << " ms\n";
            break;
        case GraphKind::Eager:
            sw.start();
            imgOut2 = computerCASF<FlatZonesGraphFullEdges>(img, radioAdj, thresholds);
            sw.pause();
            std::cout << "Time our approach:   " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut3 = computerCASF_hybrid<FlatZonesGraphFullEdges>(img, radioAdj, thresholds, 2);
            sw.pause();
            std::cout << "Time our approach_hybrid: " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut4 = computerCASF_subtree<FlatZonesGraphFullEdges>(img, radioAdj, thresholds);
            sw.pause();
            std::cout << "Time our approach (subtree): " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut5 = computerCASF_hybridSubtree<FlatZonesGraphFullEdges>(img, radioAdj, thresholds, 1);
            sw.pause();
            std::cout << "Time our approach_hybrid (subtree): " << elapsedMs(sw) << " ms\n";
            break;
        case GraphKind::Naive:
            sw.start();
            imgOut2 = computerCASF<FlatZonesGraphOnDemandEdgesByPixel>(img, radioAdj, thresholds);
            sw.pause();
            std::cout << "Time our approach:   " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut3 = computerCASF_hybrid<FlatZonesGraphOnDemandEdgesByPixel>(img, radioAdj, thresholds, 2);
            sw.pause();
            std::cout << "Time our approach_hybrid: " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut4 = computerCASF_subtree<FlatZonesGraphOnDemandEdgesByPixel>(img, radioAdj, thresholds);
            sw.pause();
            std::cout << "Time our approach (subtree): " << elapsedMs(sw) << " ms\n";

            sw.start();
            imgOut5 = computerCASF_hybridSubtree<FlatZonesGraphOnDemandEdgesByPixel>(img, radioAdj, thresholds, 1);
            sw.pause();
            std::cout << "Time our approach_hybrid (subtree): " << elapsedMs(sw) << " ms\n";
            break;
    }
    
    std::cout << "The images (naive/our) are equals: " << (imgOut1->isEqual(imgOut2)? "True":"False") << "\n";
    std::cout << "The images (naive/our_hybrid) are equals: " << (imgOut1->isEqual(imgOut3)? "True":"False") << "\n";
    std::cout << "The images (naive/our_subtree) are equals: " << (imgOut1->isEqual(imgOut4)? "True":"False") << "\n";
    std::cout << "The images (naive/our_hybridSubtree) are equals: " << (imgOut1->isEqual(imgOut5)? "True":"False") << "\n\n";

    /*
    data = new unsigned char[n];./
    for (int i = 0; i < numCols * numRows; i++) {
        data[i] = static_cast<int>(imgOut3[i]);  // Converte de `unsigned char` para `int`
    }
    stbi_write_png(("/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/out_our_" + entry.path().filename().string()).c_str(), numCols, numRows, 1, data, 0);
    delete[] data;
    */

    return 0;
}
