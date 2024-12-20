
#include "../include/CASF.hpp"


int* CASF::doFilter(int* img, int numRows, int numCols, int thresholds[]){
    ComponentTree maxtree(img, numRows, numCols, true, 1.5);
    ComponentTree mintree(img, numRows, numCols, false, 1.5);

    return mintree.reconstructionImage();
}
