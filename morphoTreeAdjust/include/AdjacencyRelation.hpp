
#include <list>
#include <vector>
#include "../include/Common.hpp"

#ifndef ADJACENCY_H
#define ADJACENCY_H

class AdjacencyRelation {
private:
    int id;
    
    int row;
    int col;    
    int numCols;
    int numRows;
    int radius;
    int n;
 
    int *offsetRow;
    int *offsetCol;
    
      

public:

    AdjacencyRelation(int numCols, int numRows, double radius);
    ~AdjacencyRelation();
    int nextValid();
    int getSize();
    AdjacencyRelation& getAdjPixels(int row, int col);
    AdjacencyRelation& getAdjPixels(int index);
    bool isAdjacent(int p, int q);
    bool isAdjacent(int px, int py, int qx, int qy);
    double getRadius();
    class IteratorAdjacency { 
    private:
        int index;
        AdjacencyRelation* instance; 

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = int;

        IteratorAdjacency(AdjacencyRelation* obj, int id) : instance(obj), index(id) { }

        AdjacencyRelation* getInstance() { return instance; } 

        IteratorAdjacency& operator++() { 
            this->index = instance->nextValid();  
            return *this; 
        }

        bool operator==(const IteratorAdjacency& other) const { 
            return index == other.index; 
        }
        bool operator!=(const IteratorAdjacency& other) const { 
            return !(*this == other);
        }

        int operator*() const { 
            return (instance->row + instance->offsetRow[index]) * instance->numCols + (instance->col + instance->offsetCol[index]); 
        }
    };
    IteratorAdjacency begin();
    IteratorAdjacency end();	 
};

#endif
