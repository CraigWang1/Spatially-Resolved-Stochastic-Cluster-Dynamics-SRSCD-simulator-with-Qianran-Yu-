//Damage.h -- Damage class to hold damage
#ifndef DAMAGE_H
#define DAMAGE_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include "constants.h"
#include "Object.h"
#include "robin_hood.h"
using namespace std;

class Damage {
private:
    // private data member:
    double DPA_RATE[POINTS];
    double NRT[POINTS];
    double damage[POINTS][CHANNELS];
    double totalIonRate;
    // private functions:
    void readFile();
    void computeDamageZero(const int);
    void computeDamageOne(const int);
    void computeDamageTwo(const int, robin_hood::unordered_flat_map<int64, Object*>&);
    void computeDamageOther(const int, const int);
public:
    Damage(robin_hood::unordered_flat_map<int64, Object*>&); // constructor
    Reaction selectDamage(const int, long double&);
    void display(const int) const; // display damage in one element
    const double getTotalDamage(const int) const;
    double getDpaRate(const int);
    double getDamageTwo(const int);
    void updateDamageTwo(const int, robin_hood::unordered_flat_map<int64, Object*>&);
    double getTotalIonRate();
};

#endif