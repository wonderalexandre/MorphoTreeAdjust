
#include <iostream>
#include <list>
#include <chrono>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustment.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"
#include "../tests/Tests.hpp"

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
    std::string filename;
    if (argc > 2) {
        std::cerr << "Use: " << argv[0] << " <file>\n";
        return 1;  
    }else if (argc == 1){
        filename = "/Users/wonderalexandre/GitHub/MorphoTreeAdjust/tests/dat/lena.png";
    }else{
        filename = argv[1];  
    }

    
    //std::string filename = argv[1];  
    std::cout << "Image: " << filename << std::endl;

    int numCols, numRows, nchannels;
    unsigned char* data = stbi_load(filename.c_str(), &numCols, &numRows, &nchannels, 1);
    
    if (!data) {
        std::cerr << "Erro: Não foi possível carregar a imagem " << filename << std::endl;
        return 1;
    }

    std::cout << "Loaded image:" << numCols << "x" << numRows << std::endl;

    int* img = new int[numCols * numRows];
    for (int i = 0; i < numCols * numRows; i++) {
        img[i] = static_cast<int>(data[i]);  // Converte de `unsigned char` para `int`
    }

    // Liberar a memória da imagem carregada
    stbi_image_free(data);

    int n = numRows * numCols;
    double radioAdj = 1.5;
    

    std::list< std::vector<int> > teste;
    
    teste.push_back( std::vector<int>{2500, 5000, 7500, 10000});
    teste.push_back( std::vector<int>{2000, 4000, 6000, 8000, 10000});
    teste.push_back( std::vector<int>{1666, 3333, 5000, 6666, 8333, 10000});
    teste.push_back( std::vector<int>{1428, 2857, 4285, 5714, 7142, 8571, 10000});
    teste.push_back( std::vector<int>{1250, 2500, 3750, 5000, 6250, 7500, 8750, 10000});
    teste.push_back( std::vector<int>{1111, 2222, 3333, 4444, 5555, 6666, 7777, 8888, 10000});
    teste.push_back( std::vector<int>{1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000});
    teste.push_back( std::vector<int>{909, 1818, 2727, 3636, 4545, 5454, 6363, 7272, 8181, 9090, 10000});
    teste.push_back( std::vector<int>{833, 1666, 2500, 3333, 4166, 5000, 5833, 6666, 7500, 8333, 9166, 10000});
    teste.push_back( std::vector<int>{769, 1538, 2307, 3076, 3846, 4615, 5384, 6153, 6923, 7692, 8461, 9230, 10000});
    teste.push_back( std::vector<int>{714, 1428, 2142, 2857, 3571, 4285, 5000, 5714, 6428, 7142, 7857, 8571, 9285, 10000});
    teste.push_back( std::vector<int>{666, 1333, 2000, 2666, 3333, 4000, 4666, 5333, 6000, 6666, 7333, 8000, 8666, 9333, 10000});
    teste.push_back( std::vector<int>{625, 1250, 1875, 2500, 3125, 3750, 4375, 5000, 5625, 6250, 6875, 7500, 8125, 8750, 9375, 10000});
    teste.push_back( std::vector<int>{588, 1176, 1764, 2352, 2941, 3529, 4117, 4705, 5294, 5882, 6470, 7058, 7647, 8235, 8823, 9411, 10000});
    teste.push_back( std::vector<int>{555, 1111, 1666, 2222, 2777, 3333, 3888, 4444, 5000, 5555, 6111, 6666, 7222, 7777, 8333, 8888, 9444, 10000});
       
    for (const std::vector<int>& thresholds : teste) {
        std::cout << "\n\n#" << thresholds.size() << "   => ";
	    for(int threshold: thresholds)
		    std::cout << threshold << ", ";
	    std::cout << "\n";

        auto start = std::chrono::high_resolution_clock::now();
        int* imgOut1 = computerCASF_naive(img, numRows, numCols, radioAdj, thresholds);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Time, naive approach: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

        start = std::chrono::high_resolution_clock::now();
        int* imgOut2 = computerCASF(img, numRows, numCols, radioAdj, thresholds);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Time, our approach: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

       
        std::cout << "The images are equals: " << (isEquals(imgOut1, imgOut2, n)? "True":"False");
        delete[] imgOut1;
        delete[] imgOut2;
    }
    delete[] img;
    return 0;
}
