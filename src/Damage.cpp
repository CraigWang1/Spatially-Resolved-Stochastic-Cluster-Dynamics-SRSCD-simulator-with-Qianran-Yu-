#include "Damage.h"
// Damage.cpp --  implementations of the damage class

Damage::Damage()
{
    int index;
    for (index = 0; index < POINTS; ++index) {
        totalRate[index] = 0.0;
    }
    readFile();
    //for(index=0; index<1; ++index){
    for (index = 0; index < POINTS; ++index) {
        computeDamageZero(index);
        if (CHANNELS > 1) {
            computeDamageOne(index);
        }
        if (CHANNELS > 2) {
            computeDamageTwo(index);
        }
        if (CHANNELS > 3) {
            int iindex;
            for (iindex = 3; iindex < CHANNELS; ++iindex) {
                computeDamageOther(index, iindex);
            }
        }
    }
}


Reaction Damage::selectDamage(const int n, long double & randRate)
{
    int index = 0;
    long double tempRate = randRate;

    if (totalRate[n] < randRate) {
        randRate -= totalRate[n];
        return NONE;
        /* this block will never be excuted */
    }
    while (index < CHANNELS) {
        if (damage[n][index] >= tempRate) {
            if (index == 0) {
                return PARTICLE;
            }
            else if (index == 1) {
                return HE;
            }
            else if (index == 2) {
                return H;
            }
        }
        else {
            tempRate -= damage[n][index];
            ++index;
        }
    }
    return NONE;
}

void Damage::display(const int count) const
{
    cout << "Damage information in element " << count << endl;
    for (int i = 0; i < CHANNELS; i++) {
        cout << "Damage[" << i << "]: " << damage[count][i] << endl;
    }
    cout << "Total Damage Rate: " << totalRate[count] << endl;
}

const double Damage::getTotalDamage(const int n) const
{
    return totalRate[n];
}

/* private method implementation */
void Damage::readFile()
{
    int index = 0;
    ifstream fd;
    fd.open("damage.txt");
    string oneLine;
    stringstream lineHold;
    while (getline(fd, oneLine)&& index < POINTS) {
        lineHold.str(oneLine);
        lineHold >> DPA_RATE[index] >> NRT[index];
        lineHold.clear();
        ++index;
    }
    fd.close();
}

void Damage::computeDamageZero(const int n)
{
    if (!IRRADIATION_ON)
    {
        damage[n][0] = 0.0;
        return;
    }

    if(n != 0){
        if( NRT[n] == 0.0 ){
            
            damage[n][0] = 0.0;
        }else{
            damage[n][0] = (DPA_RATE[n] * DENSITY*VOLUME / NRT[n]);
        }
        
    }else{
        damage[n][0] = 0.0;
    }

    //damage[n][0] = DPA_RATE[n] * DENSITY*VOLUME / NRT[n];
    totalRate[n] += damage[n][0];
}

void Damage::computeDamageOne(const int n)
{
    damage[n][1]= 0.;
    // damage[n][1]= RATIO_HE*1.0e-06*DPA_RATE[n]*DENSITY*VOLUME;
    totalRate[n] += damage[n][1];
}

void Damage::computeDamageTwo(const int n)
{
    if (!HYDROGEN_ON)
    {
        damage[n][2] = 0.0;
        return;
    }

    double concentration_H = 1.34e+4;
    if (n == 0) {
        damage[n][2] = concentration_H * FLUX_H * SURFACE_VOLUME; // rate (num hydrogen insertions/s)
        // cout << "damage2: " << damage[n][2] << endl;
    }
    else {
        damage[n][2] = 0.0;
    }
    //damage[n][2] = RATIO_H*1.0e-06*DPA_RATE[n] * DENSITY*VOLUME;
    totalRate[n] += damage[n][2];
}

void Damage::computeDamageOther(const int n, const int m)
{
    damage[n][m] = 0.0;
    totalRate[n] += damage[n][m];
}

double Damage::getDpaRate(const int n){
    return DPA_RATE[n];
}

double Damage::getDamageTwo(const int n)
{
    return damage[n][2];
}