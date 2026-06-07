#include "Damage.h"
// Damage.cpp --  implementations of the damage class

Damage::Damage(unordered_map<int64, Object*>& allObjects)
{
    int index;
    readFile();
    totalIonRate = 0.0;

    for (int n = 0; n < POINTS; n++) {
        double elementLen = length(n);
        if (n == 0) {
            positions[n] = elementLen / 2.;
        }
        else {
            positions[n] = positions[n-1] + length(n-1)/2. + elementLen/2.;
        }
    }

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
            damage[n][0] = (DPA_RATE[n] * DENSITY*volumeAtIndex(n) / NRT[n]);
        }
        
    }else{
        damage[n][0] = 0.0;
    }

    totalIonRate += damage[n][0];
    //damage[n][0] = DPA_RATE[n] * DENSITY*volumeAtIndex(n) / NRT[n];
}

void Damage::computeDamageOne(const int n)
{
    damage[n][1]= 0.;
    // damage[n][1]= RATIO_HE*1.0e-06*DPA_RATE[n]*DENSITY*volumeAtIndex(n);
}

void Damage::computeDamageTwo(const int n, unordered_map<int64, Object*>& allObjects)
{
    if (!HYDROGEN_ON)
    {
        damage[n][2] = 0.0;
        return;
    }

    // Valid for ion energies of 50-200eV
    double reflectionCoeff = -0.074 * log(H_DEPOSITION_ENERGY) + 0.96; // Data regressions from Ogorodnikova 2015
    double meanRange = (0.393651 * sqrt(H_DEPOSITION_ENERGY) - 0.000805616) * NM_TO_CM;
    double penetratingRate = FLUX_H * (1 - reflectionCoeff) * DIVIDING_AREA;

    if (positions[n] <= meanRange && n < POINTS-1 && positions[n+1] >= meanRange)
    {
        damage[n][2] = penetratingRate * (positions[n+1] - meanRange) / (positions[n+1] - positions[n]);
    }
    else if (positions[n] >= meanRange && n >= 1 && positions[n-1] <= meanRange)
    {
        damage[n][2] = penetratingRate * (meanRange - positions[n-1]) / (positions[n] - positions[n-1]);
    }
    else
    {
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
