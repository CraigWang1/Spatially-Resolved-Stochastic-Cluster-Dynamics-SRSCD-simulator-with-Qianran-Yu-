// Bundle.h -- struct Bundle
#ifndef BUNDLE_H
#define BUNDLE_H

#include"OneLine.h"
#include"robin_hood.h"
struct Bundle {
    OneLine* lines[POINTS];  /* pointers that point to line */
    Bundle(const Object* const, robin_hood::unordered_flat_map<int64, Object*>&, robin_hood::unordered_flat_map<int64, Object*>&);/* constructor */
    ~Bundle();               /* destructor */
};

#endif
