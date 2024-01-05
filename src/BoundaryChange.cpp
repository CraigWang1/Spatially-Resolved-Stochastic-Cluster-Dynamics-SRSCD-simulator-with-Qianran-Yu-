#include "BoundaryChange.h"

BoundaryChange::BoundaryChange(
               int64 oKey,
               int i,
               int ch)
{
    objKey = oKey;
    pointIndex = i;
    change = ch;
}