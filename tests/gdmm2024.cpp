
#include <iostream>
#include <list>
#include <chrono>

#include "../morphoTreeAdjust/include/NodeCT.hpp"
#include "../morphoTreeAdjust/include/ComponentTree.hpp"
#include "../morphoTreeAdjust/include/AdjacencyRelation.hpp"
#include "../morphoTreeAdjust/include/ComponentTreeAdjustment.hpp"

#include "./external/stb/stb_image.h"
#include "./external/stb/stb_image_write.h"


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
    long ti = 0;
    
    int* imgOut = new int[numRows*numCols];
    for(int i=0; i < numRows*numCols; i++)
        imgOut[i] = img[i];

    for(int threshold: thresholds) {

        ComponentTree maxtree(imgOut, numRows, numCols, true, radioAdj);
	    for(NodeCT* node: maxtree.getNodesThreshold(threshold)) {
	        maxtree.prunning(node);    
	    }
	    imgOut = maxtree.reconstructionImage();


	    ComponentTree mintree(imgOut, numRows, numCols, false, radioAdj);
	    for(NodeCT* node: mintree.getNodesThreshold(threshold)) {
	        mintree.prunning(node);    
	    }
	    imgOut = mintree.reconstructionImage();        	    
	}
    

    long tf = 0;
    long time = ti - tf;
    return imgOut;
}

int* computerCASF(int* img, int numRows, int numCols, double radioAdj, std::vector<int> thresholds){
    ComponentTree maxtree(img, numRows, numCols, true, radioAdj);
    ComponentTree mintree(img, numRows, numCols, false, radioAdj);
    ComponentTreeAdjustment adjust;
    for(int threshold: thresholds) {
		adjust.adjustMaxTree(maxtree, mintree, mintree.getNodesThreshold(threshold));
		adjust.adjustMinTree(mintree, maxtree, maxtree.getNodesThreshold(threshold)); 
	}
    int* imgOut = mintree.reconstructionImage();
    return imgOut;
}

 bool isEquals(int* imgOut1, int* imgOut2, int size){
    int equals = 0;
    for(int p=0; p < size; p++){
        if(imgOut1[p] != imgOut2[p]){
            equals++;
            
        }
    }
    std::cout << "\tDiff:" << equals <<"\t";
    return equals == 0;
 }

int main()
{

    std::cout << "DGMM...\n";

    int* img=new int[255]{
        122, 127, 166, 201, 152,  96,  54,  44,  40,  41,  42,  43,  44,
        44,  37, 133, 143, 213, 246, 236, 196, 137,  85,  55,  43,  44,
        45,  35,  40,  42, 133, 168, 231, 242, 246, 246, 228, 172, 111,
        74,  76,  80,  54,  52,  41, 147, 215, 222, 199, 220, 235, 244,
       237, 205, 172, 181, 186, 106,  57,  47, 164, 235, 224, 149, 168,
       208, 231, 244, 248, 246, 246, 230, 133,  58,  62, 140, 224, 237,
       161, 128, 149, 180, 227, 245, 248, 247, 243, 189, 103,  94, 134,
       211, 240, 181, 109, 105, 120, 168, 223, 240, 241, 246, 237, 176,
       110, 117, 188, 244, 210, 111,  74,  86, 144, 215, 230, 219, 227,
       232, 212, 133,  66, 159, 242, 238, 149,  75,  78, 163, 238, 212,
       172, 198, 219, 175, 111,  75, 144, 231, 244, 171,  81, 113, 212,
       222, 149, 108, 115, 137, 118,  99,  78, 139, 222, 245, 185, 115,
       176, 229, 176,  85,  62,  79,  95,  98, 107,  48, 102, 199, 241,
       220, 171, 220, 208, 125,  47,  45,  73,  90,  98, 104,  41,  72,
       171, 240, 242, 233, 226, 149,  65,  39,  60,  97, 104, 106, 112,
        54,  68, 140, 228, 238, 236, 194, 100,  44,  48,  85, 100, 104,
       107, 122,  54,  54,  94, 181, 222, 214, 141,  67,  40,  72,  99,
       105, 106, 109, 123,  54,  48,  59,  95, 145, 158,  84,  52,  60,
        96, 110, 115, 116, 110, 113,  49,  45,  44,  48,  71,  89,  49,
        47,  71,  95, 162, 156, 119, 122, 111};
    

    int numRows=17;
    int numCols=15;
    int n = numRows * numCols;
    double radioAdj = 1.5;
    
    // read image
    //int nchannels;
    //unsigned char *data = stbi_load("/Users/wonderalexandre/lena.png", &numCols, &numRows, &nchannels, 1);
    //stbi_image_free(data);


    std::list< std::vector<int> > teste;
    teste.push_back( std::vector<int>{20, 50});
    teste.push_back( std::vector<int>{100, 200});
    teste.push_back( std::vector<int>{66, 133, 200});
    teste.push_back( std::vector<int>{50, 100, 150, 200});
    teste.push_back( std::vector<int>{40, 80, 120, 160, 200});
    teste.push_back( std::vector<int>{33, 66, 100, 133, 166, 200});
    teste.push_back( std::vector<int>{28, 57, 85, 114, 142, 171, 200});
    teste.push_back( std::vector<int>{25, 50, 75, 100, 125, 150, 175, 200});
    teste.push_back( std::vector<int>{22, 44, 66, 88, 111, 133, 155, 177, 200});
    teste.push_back( std::vector<int>{20, 40, 60, 80, 100, 120, 140, 160, 180, 200});
    teste.push_back( std::vector<int>{18, 36, 54, 72, 90, 109, 127, 145, 163, 181, 200});
    teste.push_back( std::vector<int>{16, 33, 50, 66, 83, 100, 116, 133, 150, 166, 183, 200});
    teste.push_back( std::vector<int>{15, 30, 46, 61, 76, 92, 107, 123, 138, 153, 169, 184, 200});
    teste.push_back( std::vector<int>{14, 28, 42, 57, 71, 85, 100, 114, 128, 142, 157, 171, 185, 200});
    teste.push_back( std::vector<int>{13, 26, 40, 53, 66, 80, 93, 106, 120, 133, 146, 160, 173, 186, 200});
    teste.push_back( std::vector<int>{12, 25, 37, 50, 62, 75, 87, 100, 112, 125, 137, 150, 162, 175, 187, 200});
    teste.push_back( std::vector<int>{11, 23, 35, 47, 58, 70, 82, 94, 105, 117, 129, 141, 152, 164, 176, 188, 200});
    teste.push_back( std::vector<int>{11, 22, 33, 44, 55, 66, 77, 88, 100, 111, 122, 133, 144, 155, 166, 177, 188, 200});
    teste.push_back( std::vector<int>{10, 21, 31, 42, 52, 63, 73, 84, 94, 105, 115, 126, 136, 147, 157, 168, 178, 189, 200});
    
    for(std::vector thresholds: teste) {
        std::cout << "\n\n#" << thresholds.size() << "   => ";
	    for(int threshold: thresholds)
		    std::cout << threshold << ", ";
	    std::cout << "\n";

        auto start = std::chrono::high_resolution_clock::now();
        int* imgOut1 = computerCASF_naive(img, numRows, numCols, radioAdj, thresholds);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "Tempo naive: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

        start = std::chrono::high_resolution_clock::now();
        int* imgOut2 = computerCASF(img, numRows, numCols, radioAdj, thresholds);
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Tempo our approach: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";

        std::cout << "Sao iguais:" << isEquals(imgOut1, imgOut2, numCols*numRows);

    }
    
    return 0;
}
