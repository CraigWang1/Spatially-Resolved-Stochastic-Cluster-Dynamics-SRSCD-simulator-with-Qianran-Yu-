#include "BoundaryChange.h"

BoundaryChange::BoundaryChange() : objKey(0), pointIndex(0), change(0)
{
}

BoundaryChange::BoundaryChange(
               int64 oKey,
               int i,
               int ch)
{
    objKey = oKey;
    pointIndex = i;
    change = ch;
}