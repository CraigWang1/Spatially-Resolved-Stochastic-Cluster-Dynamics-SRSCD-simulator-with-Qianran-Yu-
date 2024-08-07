// CascadeDamage.h -- class CascadeDamage
#ifndef CASCADE_DAMAGE_H
#define CASCADE_DAMAGE_H

#include<vector>
#include<cstdlib>
#include<cmath>
#include"rvgs.h"
#include"constants.h"

using namespace std;
class CascadeDamage {
private:
    vector<int*> damage;
public:
    CascadeDamage() {};/* constructor */
    ~CascadeDamage(); /* destructor */
    /* generate cascade damage functions */
    void generateNeutronDamage(const double, int&);
    void generateIonDamage(const double, int&);
    const int getDamage(const int, const int) const;
    void cleanDamage();
    int size();
};

#endif