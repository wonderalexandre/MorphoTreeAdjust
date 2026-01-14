
#include <iostream>
#include <list>
#include <chrono>
#include <filesystem>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentByLeaf.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustmentBySubtree.hpp"
#include "../morphoTreeAdjust/include/FlatZonesGraph.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"
#include "../tests/Tests.hpp"

namespace fs = std::filesystem;
const bool ENABLE_PRINT = true;


ImageUInt8Ptr computerCASF_naive(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    start_all = std::chrono::high_resolution_clock::now();
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            
            start = std::chrono::high_resolution_clock::now();
        }
        //auto startTmp1 = std::chrono::high_resolution_clock::now();
        ComponentTreePPtr maxtreePtr = std::make_shared<ComponentTreeP>(imgOut, true, adj);
        ComponentTreeP* maxtree = maxtreePtr.get();
        //auto timeTmp1 = std::chrono::high_resolution_clock::now() - startTmp1;
        

	    for(NodeId node: ComponentTreeP::getNodesThreshold(maxtree, threshold)) {
	        maxtree->prunning(node);    
	    }
        
	    imgOut = maxtree->reconstructionImage();
        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec maxtree) naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
	    //auto startTmp2 = std::chrono::high_resolution_clock::now();
        ComponentTreePPtr mintreePtr = std::make_shared<ComponentTreeP>(imgOut, false, adj);
        ComponentTreeP* mintree = mintreePtr.get();
        //auto timeTmp2 = std::chrono::high_resolution_clock::now() - startTmp2;
       // if(i==0) std::cout << "\t- Time build trees: " << std::chrono::duration_cast<std::chrono::milliseconds>(timeTmp1+timeTmp2).count() << " ms\n";

	    for(NodeId node: ComponentTreeP::getNodesThreshold(mintree, threshold)) {
	        mintree->prunning(node);    
	    }

	    imgOut = mintree->reconstructionImage();  
        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec mintree) naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }
        
	}
    
    return imgOut;
}


ImageUInt8Ptr computerCASF_hybrid(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds, size_t cutoffPointHybrid){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    if(ENABLE_PRINT){
        start = std::chrono::high_resolution_clock::now();
    }

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    
    for(size_t i=0; i < cutoffPointHybrid; i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start = std::chrono::high_resolution_clock::now();
            start_all = std::chrono::high_resolution_clock::now();
        }
        
        ComponentTreePPtr maxtreePtr = std::make_shared<ComponentTreeP>(imgOut, true, adj);
        ComponentTreeP* maxtree = maxtreePtr.get();
	    for(NodeId node: ComponentTreeP::getNodesThreshold(maxtree, threshold)) {
	        maxtree->prunning(node);    
	    }
        
	    imgOut = maxtree->reconstructionImage();
        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            //std::cout << "\t- Time (build/prunning/rec maxtree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            std::cout << "\t- Time (update mintree and pruning maxtree) hybrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
	    ComponentTreePPtr mintreePtr = std::make_shared<ComponentTreeP>(imgOut, false, adj);
        ComponentTreeP* mintree = mintreePtr.get();
        
	    for(NodeId node: ComponentTreeP::getNodesThreshold(mintree, threshold)) {
	        mintree->prunning(node);    
	    }

	    imgOut = mintree->reconstructionImage();  
        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update maxtree and pruning mintree) hybrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }

    }

    start = std::chrono::high_resolution_clock::now();
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(imgOut, adj);
    ComponentTreeFZPtr maxtreePtr = std::make_shared<ComponentTreeFZ>(graph, true);
    ComponentTreeFZPtr mintreePtr = std::make_shared<ComponentTreeFZ>(graph, false);

    ComponentTreeFZ* mintree = mintreePtr.get();
    ComponentTreeFZ* maxtree = maxtreePtr.get();

    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);
    
    if(ENABLE_PRINT){
        end = std::chrono::high_resolution_clock::now();
        //std::cout << "\tTime (build trees): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    }
    
    for(size_t i=cutoffPointHybrid; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start_all = std::chrono::high_resolution_clock::now();
            start = std::chrono::high_resolution_clock::now();
        }
        auto nodesToPruning = ComponentTreeFZ::getNodesThreshold(maxtree, threshold);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        
        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update mintree and pruning maxtree) hybrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
        
        nodesToPruning = ComponentTreeFZ::getNodesThreshold(mintree, threshold);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning); 

        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update maxtree and pruning mintree) hybrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }

	}
    imgOut = mintree->reconstructionImage();
    return imgOut;
}


ImageUInt8Ptr computerCASF(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    if(ENABLE_PRINT){
        start = std::chrono::high_resolution_clock::now();
    }

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    
    //auto startTmp1 = std::chrono::high_resolution_clock::now();
    std::shared_ptr<FlatZonesGraph> graph = std::make_shared<FlatZonesGraph>(img, adj);
    //std::cout << "\t- Time build graph: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTmp1).count() << " ms\n";

    //auto startTmp2 = std::chrono::high_resolution_clock::now();
    ComponentTreeFZPtr maxtreePtr = std::make_shared<ComponentTreeFZ>(graph, true);
    ComponentTreeFZPtr mintreePtr = std::make_shared<ComponentTreeFZ>(graph, false);
    ComponentTreeFZ* mintree = mintreePtr.get();
    ComponentTreeFZ* maxtree = maxtreePtr.get();
    //std::cout << "\t- Time build trees: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTmp2).count() << " ms\n";
    //std::cout << "\t- Time build trees + graph: " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTmp1).count() << " ms\n";

    ComponentTreeAdjustmentByLeaf adjust(mintree, maxtree);
    
    if(ENABLE_PRINT){
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\tTime (build trees): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    }
    
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        
        if(ENABLE_PRINT){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start_all = std::chrono::high_resolution_clock::now();
            start = std::chrono::high_resolution_clock::now();
        }
        auto nodesToPruning = ComponentTreeFZ::getNodesThreshold(maxtree, threshold);
        adjust.adjustMinTree(mintree, maxtree, nodesToPruning);
        end = std::chrono::high_resolution_clock::now();
        
        if(ENABLE_PRINT){
            std::cout << "\t- Time (update mintree and pruning maxtree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
        nodesToPruning = ComponentTreeFZ::getNodesThreshold(mintree, threshold);
        adjust.adjustMaxTree(maxtree, mintree, nodesToPruning); 
        
        if(ENABLE_PRINT){
            end = std::chrono::high_resolution_clock::now();
            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update maxtree and pruning mintree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }

       
	}
    auto imgOut = mintree->reconstructionImage();
    return imgOut;
}




int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Use: " << argv[0] << " <image>\n";
        return 1;
    }

    std::string filename = argv[1];
    std::cout << "Image: " << filename << std::endl;
    
    std::vector<int> thresholds = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
    //std::vector<int> thresholds = {50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500};
    
    ImageUInt8Ptr img = nullptr;
    if(!std::filesystem::exists(filename)) {
        img = getSimpleImage();
        std::cout << "Using default image: Wonder" << std::endl;
    }else{
        int numCols, numRows, nchannels;
        uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
        
        if (!data) {
            std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
            return 1;
        }
        std::cout << "Resolution: " << numCols << "x" << numRows << std::endl;
        img = ImageUInt8::fromRaw(data, numRows, numCols);
        
    }
    
    double radioAdj = 1.5;
    
    auto start = std::chrono::high_resolution_clock::now();
    auto imgOut1 = computerCASF_naive(img, radioAdj, thresholds);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time naive approach: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n\n";

    start = std::chrono::high_resolution_clock::now();
    auto imgOut2 = computerCASF(img, radioAdj, thresholds);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time our approach:   " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n\n";

    start = std::chrono::high_resolution_clock::now();
    int cutoffPointHybrid = 1;
    auto imgOut3 = computerCASF_hybrid(img, radioAdj, thresholds, cutoffPointHybrid);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time our approach_hybrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n\n";
    
    std::cout << "The images (naive/our) are equals: " << (imgOut1->isEqual(imgOut2)? "True":"False") << "\n";
    std::cout << "The images (naive/our_hybrid) are equals: " << (imgOut1->isEqual(imgOut3)? "True":"False") << "\n\n";
    
    /*
    data = new unsigned char[n];
    for (int i = 0; i < numCols * numRows; i++) {
        data[i] = static_cast<int>(imgOut3[i]);  // Converte de `unsigned char` para `int`
    }
    stbi_write_png(("/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/out_our_" + entry.path().filename().string()).c_str(), numCols, numRows, 1, data, 0);
    delete[] data;
    */

    return 0;
}
