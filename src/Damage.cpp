#include "Damage.h"
// Damage.cpp --  implementations of the damage class

Damage::Damage(unordered_map<int64, Object*>& allObjects)
{
    int index;
    readFile();
    totalIonRate = 0.0;
    for (index = 0; index < POINTS; ++index) {
        computeDamageZero(index);
        if (CHANNELS > 1) {
            computeDamageOne(index);
        }
        if (CHANNELS > 2) {
            computeDamageTwo(index, allObjects);
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

    double totalRate = getTotalDamage(n);
    if (totalRate < randRate) {
        randRate -= totalRate;
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
    cout << "Total Damage Rate: " << getTotalDamage(count) << endl;
}

const double Damage::getTotalDamage(const int n) const
{
    double totalRate = 0.0;
    for (int i = 0; i < CHANNELS; i++)
    {
        totalRate += damage[n][i];
    }
    return totalRate;
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

    totalIonRate += damage[n][0];
    //damage[n][0] = DPA_RATE[n] * DENSITY*VOLUME / NRT[n];
}

void Damage::computeDamageOne(const int n)
{
    damage[n][1]= 0.;
    // damage[n][1]= RATIO_HE*1.0e-06*DPA_RATE[n]*DENSITY*VOLUME;
}

void Damage::computeDamageTwo(const int n, unordered_map<int64, Object*>& allObjects)
{
    if (!HYDROGEN_ON)
    {
        damage[n][2] = 0.0;
        return;
    }

    double reflectionCoeff = -0.074 * log(H_DEPOSITION_ENERGY) + 0.96; // Data regression from Ogorodnikova 2015
    double maxSurfaceConc = pow(DENSITY, 2.0/3.0);   // insertion rate formula form Zhenhou Wang 2020
    double surfaceConc = 0.0;

    int64 HKey = 1;
    if (allObjects.find(HKey) != allObjects.end())
        surfaceConc = allObjects[HKey]->getNumber(0) / DIVIDING_AREA;  // [cm^-2] concentration

    double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;

    if (n == 0 && surfaceConc < maxSurfaceConc) {
        damage[n][2] = FLUX_H * (1 - reflectionCoeff) * (1 - surfaceSaturationFraction) * (1 - H_DIRECT_IMPLANTATION_FRACTION) * DIVIDING_AREA;
    }
    else if (n == 1) {
        damage[n][2] = FLUX_H * (1 - reflectionCoeff) * H_DIRECT_IMPLANTATION_FRACTION * DIVIDING_AREA;
    }
    else {
        damage[n][2] = 0.0;
    }
}

void Damage::computeDamageOther(const int n, const int m)
{
    damage[n][m] = 0.0;
}

double Damage::getDpaRate(const int n){
    return DPA_RATE[n];
}

double Damage::getDamageTwo(const int n)
{
    return damage[n][2];
}

void Damage::updateDamageTwo(const int n, unordered_map<int64, Object*>& allObjects)
{
    computeDamageTwo(n, allObjects);
}

double Damage::getTotalIonRate()
{
    return totalIonRate;
}
