#ifndef BOUNDARY_CHANGE_H
#define BOUNDARY_CHANGE_H

#include "Bundle.h"

/* 
 * Stores information about changes near the processor domain boundaries 
 * Records:
 *  1) key of object which count has changed
 *  2) which volume element the change occurred
 *  3) how much it changed (eg. +/-1) 
 */
struct BoundaryChange
{
    int64 objKey;
    int pointIndex;
    int change;
    
    BoundaryChange(int64, int, int); // Constructor for object
};

#endif