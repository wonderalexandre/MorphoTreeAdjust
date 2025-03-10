
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

int* computerCASF_naive(int* img, int numRows, int numCols, double radioAdj, std::vector<int> thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    int* imgOut = new int[numRows*numCols];
    for(int i=0; i < numRows*numCols; i++)
        imgOut[i] = img[i];

    for(int i=0; i < thresholds.size(); i++) {
        int threshold = thresholds[i];
        if(PRINT_LOG){
            std::cout << "Opening/Closing: " << (i+1) << " \t\tthreshold:" << threshold << std::endl;
            start_all = std::chrono::high_resolution_clock::now();
            start = std::chrono::high_resolution_clock::now();
        }
        ComponentTreeP* maxtree = new ComponentTreeP(imgOut, numRows, numCols, true, radioAdj);
	    for(NodeP* node: maxtree->getNodesThreshold(threshold)) {
	        maxtree->prunning(node);    
	    }
        delete[] imgOut;
	    imgOut = maxtree->reconstructionImage();
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec maxtree) naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
            start = std::chrono::high_resolution_clock::now();
        }
	    ComponentTreeP* mintree = new ComponentTreeP(imgOut, numRows, numCols, false, radioAdj);
	    for(NodeP* node: mintree->getNodesThreshold(threshold)) {
	        mintree->prunning(node);    
	    }

        delete[] imgOut;
	    imgOut = mintree->reconstructionImage();  
        if(PRINT_LOG){
            end = std::chrono::high_resolution_clock::now();
            std::cout << "\t- Time (build/prunning/rec mintree) naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

            auto end_all = std::chrono::high_resolution_clock::now();
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_all - start_all).count() << " ms\n\n";
        }
        delete maxtree;
        delete mintree;      	    
	}
    
    return imgOut;
}


int* computerCASF(int* img, int numRows, int numCols, double radioAdj, std::vector<int> thresholds){
    std::chrono::high_resolution_clock::time_point start, start_all, end, end_all;
    if(PRINT_LOG){
        start = std::chrono::high_resolution_clock::now();
    }

    ComponentTreeFZ* maxtree = new ComponentTreeFZ(img, numRows, numCols, true, radioAdj);
    ComponentTreeFZ* mintree = new ComponentTreeFZ(img, numRows, numCols, false, radioAdj);
    
    ComponentTreeAdjustment adjust(mintree, maxtree);
    
    if(PRINT_LOG){
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\n\tTime (build trees): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
    }
    
    for(int i=0; i < thresholds.size(); i++) {
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
    int* imgOut = mintree->reconstructionImage();
    delete maxtree;
    delete mintree;
    return imgOut;
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Uso: " << argv[0] << " <diretorio_imagens>\n";
        return 1;
    }

    std::string directoryPath = argv[1];

    if (!fs::is_directory(directoryPath)) {
        std::cerr << "Caminho inválido ou não é um diretório.\n";
        return 1;
    }

    std::vector<int> thresholds{50, 176, 303, 430, 557, 684, 811, 938, 1065, 1192, 1319, 1446, 1573, 1700, 1826, 1953, 2080, 2207, 2334, 2461, 2588, 2715, 2842, 2969, 3096, 3223, 3350, 3476, 3603, 3730, 3857, 3984, 4111, 4238, 4365, 4492, 4619, 4746, 4873, 5000, 5400, 5911, 6422, 6933, 7444, 7955, 8466, 8977, 9488, 10000};
    std::cout << "Thresholds: {";
    for(int i=0; i < thresholds.size(); i++)
        std::cout << thresholds[i] << (i+1 < thresholds.size()? ", ": "}");
    std::cout << "\n\n";

    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.path().extension() != ".png") continue;
        
        std::string filename = entry.path().string();
        std::cout << "Image: " << filename << std::endl;

        int numCols, numRows, nchannels;
        unsigned char* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
        
        if (!data) {
            std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
            return 1;
        }

        std::cout << "Resolution: " << numCols << "x" << numRows << std::endl;

        int* img = new int[numCols * numRows];
        for (int i = 0; i < numCols * numRows; i++) {
            img[i] = static_cast<int>(data[i]);  // Converte de `unsigned char` para `int`
        }

        // Liberar a memória da imagem carregada
        stbi_image_free(data);

        int n = numRows * numCols;
        double radioAdj = 1.5;
        
        auto start = std::chrono::high_resolution_clock::now();
        int* imgOut1 = computerCASF_naive(img, numRows, numCols, radioAdj, thresholds);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Times" << std::endl;
        std::cout << "\tnaive approach: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
        
        start = std::chrono::high_resolution_clock::now();
        int* imgOut2 = computerCASF(img, numRows, numCols, radioAdj, thresholds);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "\tour approach:   " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

        
        std::cout << "The images are equals: " << (isEquals(imgOut1, imgOut2, n)? "True":"False") << "\n\n";
        delete[] imgOut1;
        delete[] imgOut2;
        delete[] img;

        //stbi_write_png(("/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/out_naive_" + filename).c_str(), numCols, numRows, 1, imgOut1, numCols * sizeof(int));
        //stbi_write_png(("/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/build/out_our_" + filename).c_str(), numCols, numRows, 1, imgOut2, numCols * sizeof(int));
    }
    return 0;
}
