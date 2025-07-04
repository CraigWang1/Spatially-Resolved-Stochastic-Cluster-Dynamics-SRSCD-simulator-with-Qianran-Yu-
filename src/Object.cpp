#include "Object.h"
#include<cmath>
// Object.cpp -- implementations of Object class

/* public method implementations */
Object::Object(
               const int64 key,
               const int count,
               const int n) :oKey(key), totalNumber(0), bindSH(0.0)
{
    setAttributes(key);
    setProperties(count, n);
}

Object::Object(const int * attr,
               const int count,
               const int n): totalNumber(0)
{
    for (int i = 0; i < LEVELS; i++) {
        attributes[i] = attr[i];
    }
    oKey = 0;
    setKey();
    setProperties(count, n);
}

Object::Object(const int64 key, const int *number):oKey(key), totalNumber(0)
{
    setAttributes(key);
    dimensionality = setDimensionality();
    computeThermalProperties();
    computeR1R1e();
    computeSinks();
    setNumber();
    for (int i = 0; i < POINTS; i++) {
        addNumber(i, number[i]);
    }
    
}


void Object::addNumber(const int count, const int n)
{
    number[count] += n;
    totalNumber += n;
}

void Object::reduceNumber(const int count)
{
    --number[count];
    --totalNumber;
}

int Object::signof(const int64 key) const
{
    return (key < 0) ? -1 : 1;
}

double Object::zero(const int defect)
{
    return (abs(defect) > 1) ? 1.0 : 0.0;
}

int64 Object::getKey() const
{
    return oKey;
}

double Object::getDiff() const
{
    return diffusivity;
}

int Object::getNumber(const int count) const
{
    return number[count];
}

int Object::getTotalNumber() const
{
    return totalNumber;
}

int Object::getAttri(const int index) const
{
    return attributes[index];
}

double Object::getSinkDislocation() const
{
    return sinkStrengthDislocation;
}

double Object::getSinkGrainBndry() const
{
    return sinkStrengthGrainBndry;
}

long double Object::getBind(const int index) const
{
    return bind[index];
}

double Object::getR1() const
{
    return r1;
}

double Object::getR1e() const
{
    return r1e;
}

int Object::getDim() const
{
    return dimensionality;
}

double Object::getBindSH() const
{
	return bindSH;
}

void Object::getThreeNumber(const int count, int* objectN) const
{
    /*
     * objectN[0] = object number in this element
     * objectN[1] = object number in the previous element
     * objectN[2] = object number in the next element
     */
    objectN[0] = number[count];
    if (count == 0) {
        /* when this is the surface element */
        objectN[1] = 0; /* Question, ask Jaime */
        objectN[2] = number[count + 1];
    }
    else if (count == POINTS - 1) {
        /* when this is the last element*/
        objectN[1] = number[count - 1];
        objectN[2] = 0; /* Question, ask Jaime */
    }
    else {
        objectN[1] = number[count - 1];
        objectN[2] = number[count + 1];
    }
}

void Object::display() const
{
    cout << "Information of Object " << oKey << ": " << endl;
    cout << "Attributes: ";
    for (int i = 0; i < LEVELS; ++i) {
        cout << attributes[i] << "    ";
    }
    cout << endl;
    cout << "number: ";
    for (int i = 0; i < POINTS; ++i) {
        cout << number[i] << "    ";
    }
    cout << endl;
    cout << "total number: " << totalNumber << endl;
    cout << "dimensionality: " << dimensionality << endl;
    cout << "diffusivity: " << diffusivity << endl;
    cout << "bind term: ";
    for (int i = 0; i < LEVELS; ++i) {
        cout << bind[i] << "    ";
    }
    cout << endl;
    cout << "r1 = " << r1 << endl;
    cout << "r1e = " << r1e << endl;
    cout << "sink strength: " << sinkStrengthDislocation <<  " " << sinkStrengthGrainBndry << endl;
    cout << endl;
}

/* private method impementations */
void Object::setKey()
{
    oKey = attrToKey(attributes);
}

void Object::setAttributes(const int64 key)
{
    int64 tempKey = abs(key);
    for (int i = 0; i < LEVELS; i++) {
        attributes[i] = double(tempKey) / pow(10.0, (double)EXP10*(LEVELS - 1 - i));
        tempKey -= ((int64)attributes[i])*((int64)pow(10.0, (double)EXP10*(LEVELS - 1 - i)));
    }
    attributes[0] *= signof(key);

    // For now, assume we are not working with He, so assume the digit is for H
    attributes[2] += attributes[1] * 1000;
    attributes[1] = 0;
}

void Object::setNumber()
{
    for (int i = 0; i < POINTS; i++) {
        number[i] = 0;
        totalNumber += number[i];
    }
}

int Object::setDimensionality()
{
    return attributes[0] > 4 ? 1 : 3;
}

void Object::computeR1R1e()
{
    int ndef = attributes[0];
    if (ndef < 0) {       // V cluster (assume to be spherical cluster)
        r1 = r1e = pow(3.0*fabs((double)ndef)*avol/4.0/PI, 1.0/3.0);
    }
    else if (ndef > 0) {  // SIA cluster (assumed to form circular loop)
        r1 = r1e = sqrt((double)ndef*avol/jumped/PI);
    }
    else if (ndef == 0 && attributes[2] > 1) { // H cluster, assumed to be square 2D platelet (Jie Hou 2018)
        r1 = r1e = sqrt((double)attributes[2])*ALATT/2.0;   // half of the side length of the square
    }
    else if (ndef == 0 && attributes[2] == 1) { // 1H
        r1 = r1e = ALATT * sqrt(3.0)/4.0;   // tetrahedral interstitial site radius
    }
    else { // shouldn't get here
        r1 = r1e = jumped;
    }
}

void Object::computeDiffCoeff()
{
    const double fi = 0.9, fv = 0.7; // Diffusion correlationm factors.
    const double gi = 0.5, gv = 0.125; // Geometric factor for diffusion.
    double prefactor = 0, energy_m = 0;
    int check_all = 0;
    int check_He = 0;
    int check_H = 0;
    // int check_C = 0;
    int i;
    
    for (i = 1; i < LEVELS; i++) {
        check_all |= attributes[i];
        if (i >= 2) check_He |= attributes[i];
        if (i >= 3) check_H |= attributes[i];
    }
    /* check_all: 0 when this is a pure defect cluster.
     *  check_He:  0 when this is a cluster with He.
     *  check_H:   0 when this is a cluster with H.
     */
    
    // Pure defect clusters:
    if (!check_all) {
        
        if (attributes[0] > 0) { // SIAs
            if (abs(attributes[0]) == 1) { // 1I
                prefactor = 8.744e-4;
                energy_m = 0.009;
            }else if (abs(attributes[0]) == 2) { // 1I
                prefactor = 7.97e-4;
                energy_m = 0.024;
            }else if (abs(attributes[0]) == 3) { // 1I
                prefactor = 3.92e-4;
                energy_m = 0.033;
            }
            else if(abs(attributes[0]) > 3) { // >1I
                prefactor = gi*jumped*jumped*fi*NU0*pow(fabs(attributes[0]), -0.5);
                energy_m = 0.013;
            }
        }
        else if (attributes[0] < 0) { // Vacancies.
            if (abs(attributes[0]) == 1) { // 1V
                
                prefactor = 1.77e-2;
                energy_m = 1.29;
            }
            else if (abs(attributes[0]) == 2) { // >1V
                
                prefactor = 2.91e-5;
                energy_m = 1.66;
            }else if (abs(attributes[0]) > 2) { // >1V
                
                prefactor = gv*jumped*jumped*fv*NU0*pow(0.001, fabs(attributes[0]) - 1.0);
                energy_m = 1.66;
            }
        }
    }
    else if (!check_He) {
        
        if (attributes[0] != 0) { // SIA-He and V-He clusters immobile.
            prefactor = 0.0;
        }
        else if (attributes[1] == 1) { // He1.
            prefactor = gi*jumped*jumped*NU0;
            energy_m = 0.01;
        }
        else  if (attributes[1] == 2) { // He2.
            prefactor = gi*jumped*jumped*NU0*0.01;
            energy_m = 0.03;
        }
        else  if (attributes[1] == 3) { // He3.
            prefactor = gi*jumped*jumped*NU0*0.01;
            energy_m = 0.05;
        }
        else  if (attributes[1] > 3) { // He>3.
            prefactor = gi*jumped*jumped*NU0*0.01;
            energy_m = 0.06;
        }
        else prefactor = 0.0;
    }
    else if (!check_H) {
        if (attributes[0] != 0) { // SIA-H and V-H clusters immobile.
            prefactor = 0.0;
        }
        else if (attributes[2] == 1) { // H1
            /*
            //number below is from ppt Liu(2015) at 300K D_H=10e-9 m^2/s = 10e-5 cm^2/s
            prefactor = 3.8e-3;
            energy_m = 0.41;
             */
            prefactor = 1.58e-3;
            energy_m = 0.25;
        }
        else
            prefactor = 0.0;
    }
    /* All data from [CS Becquart et al., J Nucl Mater 403 (2010) 75] */
    diffusivity = prefactor*exp(-energy_m / KB / TEMPERATURE);
}

void Object::computeBindTerm()
{
    long double energy_d[LEVELS] = { 0.0 };
    long double energy_b = 0.0;
    double attfreq = 1.0;
    double efi = 9.96, emi = 0.013; // Ab initio migration and formation energies of V and SIA in pure W.
    double efv = 3.23, emv = 1.66;
    double eb2v = -0.1, eb2i = 2.12, eb2he = 1.03;
    double efhe = 4.0, emhe = 0.01;
    double emh = H_MIGRATION_ENERGY;
    int check_all = 0, check_He = 0, check_H = 0;
    int i;
    
    for (i = 0; i < LEVELS; i++) {
        bind[i] = 0.0;
        if (i >= 1) check_all |= attributes[i];
        if (i >= 2) check_He |= attributes[i];
        if (i >= 3) check_H |= attributes[i];
    }
    /**
     * bind energy positive(eg. 5 eV) means easy to get together, hard to dissociate, when dissociate, absorb 5eV energy
     * bind energy negative(eg. -5 eV) means hard to get together, easy to dissociate, when dissociate, release 5eV energy
     **/
    // Pure defect clusters:
    if (!check_all) {
        if (attributes[0]>0) { // SIAs
            if (abs(attributes[0]) == 1) { // 1I
                attfreq = 0.0;
            }
            else if (abs(attributes[0]) == 2) { // 2I
                energy_b = 2.12;
            }
            else if (abs(attributes[0]) == 3) { // 3I
                energy_b = 3.02;
            }
            else if (abs(attributes[0]) == 4) { // 4I
                energy_b = 3.60;
            }
            else if (abs(attributes[0]) == 5) { // 5I
                energy_b = 3.98;
            }
            else if (abs(attributes[0]) == 6) { // 6I
                energy_b = 4.27;
            }
            else if (abs(attributes[0]) == 7) { // 7I
                energy_b = 5.39;
            }
            else if (abs(attributes[0])>7) // > 7I
                energy_b = efi + (eb2i - efi)*(pow(fabs((double)attributes[0]), 0.6666667) - pow((fabs((double)attributes[0]) - 1.0), 0.6666667)) / 0.5847;
        }
        else if (attributes[0]<0) { // Vacancies.
            if (abs(attributes[0]) == 1) { // 1V
                attfreq = 0.0;
            }
            else if (abs(attributes[0]) == 2) { // 2V
                energy_b = eb2v;
            }
            else if (abs(attributes[0]) == 3) { // 3V
                energy_b = 0.04;
            }
            else if (abs(attributes[0]) == 4) { // 4V
                energy_b = 0.64;
            }
            else if (abs(attributes[0]) == 5) { // 5V
                energy_b = 0.72;
            }
            else if (abs(attributes[0]) == 6) { // 6V
                energy_b = 0.89;
            }
            else if (abs(attributes[0]) == 7) { // 7V
                energy_b = 0.72;
            }
            else if (abs(attributes[0]) == 8) { // 8V
                energy_b = 0.88;
            }
            else if (abs(attributes[0])>8) // > 8V
                energy_b = efv + (eb2v - efv)*(pow(fabs((double)attributes[0]), 0.6666667) - pow((fabs((double)attributes[0]) - 1.0), 0.6666667)) / 0.5874;
        }
        energy_d[0] = energy_b;
        bind[0] = attfreq*exp(-energy_d[0] / KB / TEMPERATURE);
    }
    // He-defect clusters:
    else if (!check_He) {
        
        if (attributes[0]<0) { // He-V clusters:
            double ratio = fabs(((double)attributes[1]) / ((double)attributes[0]));
            printf("%dV - %dHe\n", abs(attributes[0]), abs(attributes[1]));
            if (abs(attributes[0]) == 1 && abs(attributes[1]) == 1) { // 1V-1He
                energy_d[0] = 4.6;
                energy_d[1] = energy_d[0];
            }
            else {
                energy_d[0] = 2.4 + 3.5*log10(ratio) + 1.7*log10(ratio)*log10(ratio);
                // binding energy of V to cluster.
                energy_d[1] = 4.6 - 1.1*log10(ratio) - 0.3*log10(ratio)*log10(ratio);
                // binding energy of He to cluster.
            }
            bind[0] = attfreq*exp(-energy_d[0] / KB / TEMPERATURE);
            bind[1] = attfreq*exp(-energy_d[1] / KB / TEMPERATURE);
        }
        else if (attributes[0]>0) // He-SIA clusters.
            attfreq = 0.0; // No dissociation between He and SIA clusters.
        else if (attributes[0] == 0) { // pure He clusters.
            if (attributes[1] == 1) { // He1.
                attfreq = 0;
            }
            else if (attributes[1] == 2) { // He2.
                energy_b = 1.03;
            }
            else if (attributes[1] == 3) { // He3.
                energy_b = 1.36;
            }
            else if (attributes[1] == 4) { // He4.
                energy_b = 1.52;
            }
            else { // He>4.
                energy_b = efhe + (eb2he - efhe)*(pow(fabs((double)attributes[0]), 0.6666667) - pow((fabs((double)attributes[0]) - 1.0), 0.6666667)) / 0.5874;
            }
            energy_d[1] = energy_b + emhe;
            bind[1] = attfreq*exp(-energy_d[1] / KB / TEMPERATURE);
        }
    }
    
    // H-defect clusters:
    else if (!check_H){
        if (attributes[0]<0) { // H-V clusters:
            double ratio = fabs( ((double) attributes[2])/((double) attributes[0]) );
            //printf("%dV - %dH\n", abs(attr[0]), abs(attr[2]));
            /*this part is for binding energy of mV-nH that try to dissociate a V from Li Xiaochun(2015)*/
            if (ratio==1 /*&& attributes[0]> -5 && attributes[2]< 15*/) { // V-H.
                energy_b = 1.60;
            } else if (ratio==2) { // V-H2
                energy_b = 2.10;
            } else if (ratio==3) { // V-H3
                energy_b = 2.40;
            } else if (ratio==4) { // V-H4
                energy_b = 3.90;
            } else if (ratio==5) { // V-H5
                energy_b = 4.35;
            } else if (ratio==6) { // V-H6
                energy_b = 5.80;
            } else if (ratio==7) { // V-H7
                energy_b = 7.0;
            } else if(ratio == 8){ // V-H8
                energy_b = 8.25;
            } else if(ratio == 9){ // V-H9
                energy_b = 9.80;
            } else if(ratio == 10){ // V-H10
                energy_b = 11.25;
            } else{
                
                energy_b = 1.91 + 0.0974 * ratio * ratio; //extrapolation
                // if happened form this, dissociate as soon as possible
            }
            energy_d[0] = energy_b + emv;

            /*this part is for binding energy of mV-nH+1H based on Ogorodnikova (2015)*/
            int numV = abs(attributes[0]);
            int numH = attributes[2];
            int monovacancyMaxH = 12;

            double radius = pow(3.0*fabs((double)numV)*avol/4.0/PI, 1.0/3.0);
            double surfArea = 4*PI*pow(radius, 2);
            double surfHDensity = numH / surfArea;

            double monovacancyRadius = pow(3.0*avol/4.0/PI, 1.0/3.0);
            double monovacancySurfArea = 4*PI*pow(monovacancyRadius, 2);
            double maxSurfHDensity = monovacancyMaxH / monovacancySurfArea;

            double maxBindE;
            if (numV == 1)  // Data from Mason 2023
                maxBindE = 1.23;
            else if (numV == 2)
                maxBindE = 1.55;
            else if (numV == 3)
                maxBindE = 1.82;
            else if (numV == 4)
                maxBindE = 1.91;
            else if (numV == 5)
                maxBindE = 1.97;
            else if (numV == 6)
                maxBindE = 1.9;
            else if (numV == 8)
                maxBindE = 1.68;
            else if (numV == 10)
                maxBindE = 1.78;
            else if (numV == 12)
                maxBindE = 1.69;
            else
                maxBindE = 1.86;

            // Regression so that when H = 1, energy_b = maxBindE and when we reach maxSurfHDensity, energy_b = 0
            energy_b = -maxBindE / pow(fabs(maxSurfHDensity - 1.0/surfArea), 1.1) * pow(fabs(surfHDensity - 1.0/surfArea), 1.1) + maxBindE;
            energy_d[2] = energy_b + emh;
            bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
            bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
        }
        else if (attributes[0]>0){ // H-SIA clusters.
            double ratio = fabs( ((double) attributes[2])/((double) attributes[0]) );

            /*this part is for binding energy of mSIA-nH that try to dissociate a SIA from Daniel Mason(2023)*/
            energy_b = 7.34 - 1.50*1.0/(1+exp(ratio/0.812));
            energy_d[0] = energy_b + emi;

            /*this part is for binding energy of mSIA-nH that try to dissociate a H from Daniel Mason(2023)*/
            energy_b = -0.147 * pow(ratio, 0.652) + 0.547;
            energy_d[2] = energy_b + emh;

            bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
            bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
        }
        else if (attributes[0]== 0) { // nH clusters, this is binding energy of nH cluster dissociating 1 H from Qin(2015)
            if (attributes[2]==1) { // H
                energy_b = 0;
            }   /* Data from Hou 2018 for planar H platelet clusters */
            else {
                energy_b = 0.38 - 0.45*exp(-attributes[2]/12.04);
            }
            // else if (attributes[2] == 2) { // 2H
            //     // energy_b = 0.04;
            //     energy_b = 1;
            // }
            // else if (attributes[2] == 3) {
            //     // energy_b = -0.03;
            //     energy_b = 2;
            // }
            // else if (attributes[2] == 4) {
            //     energy_b = 0.32;
            // }
            // else if (attributes[2] == 5) {
            //     energy_b = 0.07;
            // }
            // else if (attributes[2] == 6) {
            //     energy_b = 0.30;
            // }
            // else if (attributes[2] == 7) {
            //     energy_b = 0.37;
            // }
            // else if (attributes[2] == 8) {
            //     energy_b = 0.28;
            // }
            // else if (attributes[2] == 9) {
            //     energy_b = 0.33;
            // }
            // else if (attributes[2] == 10) {
            //     energy_b = 0.35;
            // }
            // else if (attributes[2] == 11) {
            //     energy_b = 0.24;
            // }
            // else if (attributes[2] == 12) {
            //     energy_b = 0.26;
            // }
            // else if (attributes[2] == 13) {
            //     energy_b = 0.35;
            // }
            // else if (attributes[2] == 14) {
            //     energy_b = 0.26;
            // }
            // else if (attributes[2] == 15) {
            //     energy_b = 0.27;
            // }
            // else if (attributes[2] == 16) {
            //     energy_b = 0.39;
            // }
            // else {
            //     energy_b = 0.32; // extrapolation
            // }



            // else {  // nH
                // energy_b = 0.38 - 0.45*exp(-(double)attributes[2]/12.04);
            // }
            // else if (attributes[2] == 2) { // 2H
            //     energy_b = 0.01;
            // }
            // else if (attributes[2] == 3) { // 3H
            //     energy_b = -0.01;
            // }
            // else if (attributes[2] == 4) {
            //     energy_b = -0.05;
            // }
            // else if (attributes[2] == 5) {
            //     energy_b = -0.06;
            // }
            // else if (attributes[2] == 6) {
            //     energy_b = -0.01;
            // }
            // else if (attributes[2] == 7) {
            //     energy_b = 0.05;
            // }
            // else if (attributes[2] == 8) {
            //     energy_b = 0.03;
            // }
            // else {
            //     energy_b = 0;
            // }

            energy_d[2] = energy_b + emh;
            if (attributes[2] >= 100000)
                energy_d[2] = 0.0;

            bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
            bind[0] = 0; //because there's no V/SIA in cluster
        }
    }
}

void Object::computeSinks()
{
    /* The total sink strength for all defects are stored in the array s.
     [0] for vacancies; [1] for SIAs; */
    double Zdv = 1.0, Zdi = 1.1;
    // double Zsv = 1.0, Zsi = 1.1;
    double Sd;
    /* 1. Dislocation sink strength: */
    Sd = DISLOCATION;
    
    /* 2. Grain boundary sink strength: */
    
    if (attributes[0] <= 0) // Vacancies or H objects
    {
        sinkStrengthDislocation = Sd * Zdv;
        sinkStrengthGrainBndry = 6*sqrt(sinkStrengthDislocation)/GRAIN_SIZE;
    }
    else // SIA objects
    {
        sinkStrengthDislocation = Sd * Zdi;
        sinkStrengthGrainBndry = 6*sqrt(sinkStrengthDislocation)/GRAIN_SIZE;
    }
}

void Object::computeThermalProperties()
{
    /* Properties that depend on temperature */
    computeDiffCoeff();
    computeBindTerm();
}

void Object::setProperties(const int count, const int n)
{
    setNumber();
    addNumber(count, n);
    dimensionality = setDimensionality();
    computeThermalProperties();
    computeR1R1e();
    computeSinks();
}

int64 attrToKey(const int * const attr)
{
    int64 key = 0;
    int sign = (attr[0] < 0) ? -1 : 1;
    for (int i = 0; i < LEVELS; ++i) {
        key += labs(attr[i])*((int64)pow(10.0, (double)EXP10*(LEVELS - 1 - i)));
    }
    key *= sign;
    return key;
}
