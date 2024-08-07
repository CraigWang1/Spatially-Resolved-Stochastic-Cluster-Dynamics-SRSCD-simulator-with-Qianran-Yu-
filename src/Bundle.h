// Bundle.h -- struct Bundle
#ifndef BUNDLE_H
#define BUNDLE_H

#include"OneLine.h"
struct Bundle {
    OneLine* lines[POINTS];  /* pointers that point to line */
    Bundle(const Object* const, unordered_map<int64, Object*>&, unordered_map<int64, Object*>&, unordered_map<Object*, Bundle*>&);/* constructor */
    ~Bundle();               /* destructor */
};

#endif
