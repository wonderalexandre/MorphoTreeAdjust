
#include <iostream>
#include <list>
#include <chrono>
#include <filesystem>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustment.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"
#include "../tests/Tests.hpp"

namespace fs = std::filesystem;


/*
import numpy as np
def imprime(n):
	X = np.linspace(0, 200, n+1)
	s = "teste.push_back( std::vector<int>{"
	for x in X:
		if x != 0:
			s += str(int(x))
		if x != 0 and x != 200:
			s += ", "
	return s+"});"

for i in range(20):
	print(imprime(i))

*/

ImageUInt8Ptr computerCASF_naive(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    start_all = std::chrono::high_resolution_clock::now();
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(PRINT_LOG){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            
            start = std::chrono::high_resolution_clock::now();
        }
        
        ComponentTreePPtr maxtree = std::make_shared<ComponentTreeP>(imgOut, true, adj);
	    for(NodePPtr node: maxtree->getNodesThreshold(threshold)) {
	        maxtree->prunning(node);    
	    }
        
	    imgOut = maxtree->reconstructionImage();
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec maxtree) naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
	    ComponentTreePPtr mintree = std::make_shared<ComponentTreeP>(imgOut, false, adj);
	    for(NodePPtr node: mintree->getNodesThreshold(threshold)) {
	        mintree->prunning(node);    
	    }

	    imgOut = mintree->reconstructionImage();  
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec mintree) naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }
        //delete maxtree;
        //delete mintree;      	    
	}
    
    return imgOut;
}


ImageUInt8Ptr computerCASF_hibrid(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    if(PRINT_LOG){
        start = std::chrono::high_resolution_clock::now();
    }

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    ImageUInt8Ptr imgOut = img->clone();
    size_t cutPoint_Hibrid = 3;
    for(size_t i=0; i < cutPoint_Hibrid; i++) {
        int threshold = thresholds[i];
        if(PRINT_LOG){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start = std::chrono::high_resolution_clock::now();
        }
        
        ComponentTreePPtr maxtree = std::make_shared<ComponentTreeP>(imgOut, true, adj);
	    for(NodePPtr node: maxtree->getNodesThreshold(threshold)) {
	        maxtree->prunning(node);    
	    }
        
	    imgOut = maxtree->reconstructionImage();
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec maxtree) hibrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
	    ComponentTreePPtr mintree = std::make_shared<ComponentTreeP>(imgOut, false, adj);
	    for(NodePPtr node: mintree->getNodesThreshold(threshold)) {
	        mintree->prunning(node);    
	    }

	    imgOut = mintree->reconstructionImage();  
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec mintree) hibrid: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }

    }


    std::unique_ptr<FlatZonesGraph> maxTreeFZGraph = std::make_unique<FlatZonesGraph>(imgOut, adj);
    std::unique_ptr<FlatZonesGraph> minTreeFZGraph = std::make_unique<FlatZonesGraph>(*maxTreeFZGraph);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(imgOut, true, adj, std::move(maxTreeFZGraph));
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(imgOut, false, adj, std::move(minTreeFZGraph));
    
    ComponentTreeAdjustment adjust(mintree, maxtree);
    
    if(PRINT_LOG){
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\tTime (build trees): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    }
    
    for(size_t i=cutPoint_Hibrid; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(PRINT_LOG){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start_all = std::chrono::high_resolution_clock::now();
            start = std::chrono::high_resolution_clock::now();
        }
        adjust.adjustMinTree(mintree, maxtree, maxtree->getNodesThreshold(threshold));
        end = std::chrono::high_resolution_clock::now();
        
        if(PRINT_LOG){
            std::cout << "\t- Time (update maxtree and pruning mintree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
        adjust.adjustMaxTree(maxtree, mintree, mintree->getNodesThreshold(threshold)); 
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update mintree and pruning maxtree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }

	}
    imgOut = mintree->reconstructionImage();
    //delete maxtree;
    //delete mintree;
    return imgOut;
}


ImageUInt8Ptr computerCASF(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    if(PRINT_LOG){
        start = std::chrono::high_resolution_clock::now();
    }

    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::unique_ptr<FlatZonesGraph> maxTreeFZGraph = std::make_unique<FlatZonesGraph>(img, adj);
    std::unique_ptr<FlatZonesGraph> minTreeFZGraph = std::make_unique<FlatZonesGraph>(*maxTreeFZGraph);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, std::move(maxTreeFZGraph));
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, std::move(minTreeFZGraph));
    
    ComponentTreeAdjustment adjust(mintree, maxtree);
    
    if(PRINT_LOG){
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\tTime (build trees): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    }
    
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(PRINT_LOG){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start_all = std::chrono::high_resolution_clock::now();
            start = std::chrono::high_resolution_clock::now();
        }
        adjust.adjustMinTree(mintree, maxtree, maxtree->getNodesThreshold(threshold));
        end = std::chrono::high_resolution_clock::now();
        
        if(PRINT_LOG){
            std::cout << "\t- Time (update maxtree and pruning mintree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
        adjust.adjustMaxTree(maxtree, mintree, mintree->getNodesThreshold(threshold)); 
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update mintree and pruning maxtree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }

	}
    auto imgOut = mintree->reconstructionImage();
    //delete maxtree;
    //delete mintree;
    return imgOut;
}



ImageUInt8Ptr computerCASF_subtree(ImageUInt8Ptr img, double radioAdj, const std::vector<int>& thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    if(PRINT_LOG){
        start = std::chrono::high_resolution_clock::now();
    }
    AdjacencyRelationPtr adj =std::make_shared<AdjacencyRelation>(img->getNumRows(), img->getNumCols(), radioAdj);
    std::unique_ptr<FlatZonesGraph> maxTreeFZGraph = std::make_unique<FlatZonesGraph>(img, adj);
    std::unique_ptr<FlatZonesGraph> minTreeFZGraph = std::make_unique<FlatZonesGraph>(*maxTreeFZGraph);

    ComponentTreeFZPtr maxtree = std::make_shared<ComponentTreeFZ>(img, true, adj, std::move(maxTreeFZGraph));
    ComponentTreeFZPtr mintree = std::make_shared<ComponentTreeFZ>(img, false, adj, std::move(minTreeFZGraph));
    
    ComponentTreeAdjustment adjust(mintree, maxtree);
    
    if(PRINT_LOG){
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\tTime (build trees): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    }
    
    for(size_t i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(PRINT_LOG){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start_all = std::chrono::high_resolution_clock::now();
            start = std::chrono::high_resolution_clock::now();
        }
        adjust.adjustMinTree2(mintree, maxtree, maxtree->getNodesThreshold(threshold));
        end = std::chrono::high_resolution_clock::now();
        
        if(PRINT_LOG){
            std::cout << "\t- Time (update maxtree and pruning mintree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
        adjust.adjustMaxTree2(maxtree, mintree, mintree->getNodesThreshold(threshold)); 
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            end_all = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (update mintree and pruning maxtree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }

	}
    auto imgOut = mintree->reconstructionImage();
    //delete maxtree;
    //delete mintree;
    return imgOut;
}


int main(int argc, char* argv[]) {
   /* if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <diretorio_imagens>\n";
        return 1;
    }*/

    std::string directoryPath = "/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/dat/dataset-icdar2024_occluded";//argv[1];

    if (!fs::is_directory(directoryPath)) {
        std::cerr << "Caminho inválido ou não é um diretório.\n";
        return 1;
    }

    //std::vector<int> thresholds{50, 176, 303, 430, 557, 684, 811, 938, 1065, 1192, 1319, 1446, 1573, 1700, 1826, 1953, 2080, 2207, 2334, 2461, 2588, 2715, 2842, 2969, 3096, 3223, 3350, 3476, 3603, 3730, 3857, 3984, 4111, 4238, 4365, 4492, 4619, 4746, 4873, 5000, 5400, 5911, 6422, 6933, 7444, 7955, 8466, 8977, 9488, 10000};
    std::vector<int> thresholds = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500};
    
   /*
    std::vector<int> thresholds = {10,   30,   50,   70,   90,  110,  130,  150,  170,  190,  210, //linear
        230,  250,  270,  290,  310,  330,  350,  370,  390,  410,  430,
        450,  470,  490,  510,  530,  550,  570,  590,  610,  630,  650,
        670,  690,  710,  730,  750,  770,  790,  810,  830,  850,  870,
        890,  910,  930,  950,  970, 1000};



    std::vector<int> thresholds = {10,   19,   31,   44,   58,   73,   89,  105,  122,  139,  157, //power=1.2
        174,  192,  211,  230,  249,  268,  287,  307,  327,  347,  368,
        388,  409,  430,  451,  472,  494,  515,  537,  559,  581,  603,
        626,  648,  671,  693,  716,  739,  762,  786,  809,  832,  856,
        880,  903,  927,  951,  975, 1000};
   
      

    std::vector<int> thresholds = {10,   14,   21,   29,   39,   50,   62,   74,   88,  102,  116, //power=1.4
        132,  148,  164,  181,  198,  216,  234,  253,  272,  292,  312,
        332,  353,  374,  395,  417,  439,  462,  485,  508,  531,  555,
        579,  603,  628,  652,  678,  703,  729,  755,  781,  807,  834,
        861,  888,  916,  943,  971, 1000};


        std::vector<int> thresholds = {10,    19,    41,    75,   119,   174,   237,   310,   392,
         482,   581,   688,   803,   926,  1057,  1196,  1342,  1496,
        1657,  1825,  2000,  2183,  2373,  2570,  2774,  2985,  3202,
        3427,  3658,  3896,  4140,  4391,  4649,  4913,  5184,  5461,
        5745,  6035,  6331,  6634,  6942,  7258,  7579,  7906,  8240,
        8580,  8926,  9278,  9636, 10000};
*/

    std::cout << "Thresholds: {";
    for(size_t i=0; i < thresholds.size(); i++)
        std::cout << thresholds[i] << (i+1 < thresholds.size()? ", ": "}");
    std::cout << "\n\n";

    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.path().extension() != ".png") continue;
        

        std::string filename = entry.path().string();
        std::cout << "Image: " << filename << std::endl;

        int numCols, numRows, nchannels;
        uint8_t* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
        
        if (!data) {
            std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
            return 1;
        }
        std::cout << "Resolution: " << numCols << "x" << numRows << std::endl;
        ImageUInt8Ptr img = ImageUInt8::fromRaw(data, numRows, numCols);
        
        double radioAdj = 1.5;
        
        auto start = std::chrono::high_resolution_clock::now();
        auto imgOut1 = computerCASF_naive(img, radioAdj, thresholds);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Times" << std::endl;
        std::cout << "\tnaive approach: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
        
        start = std::chrono::high_resolution_clock::now();
        auto imgOut2 = computerCASF(img, radioAdj, thresholds);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\tour approach:   " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

        start = std::chrono::high_resolution_clock::now();
        auto imgOut3 = computerCASF_subtree(img,radioAdj, thresholds);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\tour approach (subtree): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

        
        std::cout << "The images (naive/our) are equals: " << (imgOut1->isEqual(imgOut2)? "True":"False") << "\n";
        std::cout << "The images (our/our2) are equals: " << (imgOut1->isEqual(imgOut3)? "True":"False") << "\n\n";
        
        /*
        data = new unsigned char[n];
        for (int i = 0; i < numCols * numRows; i++) {
            data[i] = static_cast<int>(imgOut3[i]);  // Converte de `unsigned char` para `int`
        }
        stbi_write_png(("/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/out_our_" + entry.path().filename().string()).c_str(), numCols, numRows, 1, data, 0);
        delete[] data;
        */
        break;
    }
    return 0;
}
