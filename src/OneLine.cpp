// OneLine.cpp -- implementations of class OneLine
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include "OneLine.h"
#include "Bundle.h"
using namespace std;

/* public function implementations */
/* Helper class to calculate reaction rates of a species in one spatial element */
OneLine::OneLine(
                 const Object* const hostObject,
                 const int count,
                 unordered_map<int64, Object*>& mobileObjects,
                 unordered_map<int64, Object*>& allObjects,
                 unordered_map<Object*, Bundle*>& linePool) :totalRate(0.0)
{
    setOneLine(hostObject, count, mobileObjects, allObjects, linePool);
}

OneLine::OneLine() : diffRToF(0.0), diffRToB(0.0), sinkRDislocation(0.0), sinkRGrainBndry(0.0), SAVR(0.0), recombRER(0.0), recombRLH(0.0), totalRate(0.0)
{
    for (int i = 0; i < LEVELS; i++)
        dissociationR[i] = 0.0;
}

Reaction OneLine::selectReaction(
                                 const Object* const hostObject,
                                 int64& theOtherKey,
                                 long double& randRate)
{
    int index = 0;
    long double tempRate = randRate;
    std::unordered_map<int64, long double>::iterator iter = secondR.begin();
    if (totalRate < tempRate) {
        randRate -= totalRate;
        return NONE;
    } // the reaction is not positioned in this line
    if(diffRToF >= tempRate){
        return DIFFUSETOF;
    }
    else {
        tempRate -= diffRToF;
    }
    if (diffRToB >= tempRate) {
        return DIFFUSETOB;
    }
    else {
        tempRate -= diffRToB;
    }
    if (sinkRDislocation >= tempRate) {
        return SINKDISLOCATION;
    }
    else {
        tempRate -= sinkRDislocation;
    }
    if (sinkRGrainBndry >= tempRate) {
        return SINKGRAINBNDRY;
    }
    else {
        tempRate -= sinkRGrainBndry;
    }
    if (SAVR >= tempRate) {
        return SAV;
    }
    else {
        tempRate -= SAVR;
    }
    if (recombRER >= tempRate) {
        return RECOMBER;
    }
    else {
        tempRate -= recombRER;
    }
    if (recombRLH >= tempRate) {
        return RECOMBLH;
    }
    else {
        tempRate -= recombRLH;
    }
    while (index < LEVELS) {
        if (dissociationR[index] >= tempRate) {
            /* generate the other cluster key for monomer! */
            int attribute = hostObject->getAttri(index);
            theOtherKey = hostObject->signof(attribute)*((int64)pow(10.0, (double)EXP10*(LEVELS - index - 1)));
            return DISSOCIATION;
        }
        else {
            tempRate -= dissociationR[index];
            ++index;
        }
    }
    while (iter != secondR.end()) {
        if (iter->second >= tempRate) {
            theOtherKey = iter->first;
            return COMBINATION;
        }
        else {
            tempRate -= iter->second;
            ++iter;
        }
    }
    return ERROR;
}

void OneLine::addReaction(
                          const Object* const hostObject,
                          const Object* const newObject,
                          unordered_map<int64, Object*>& allObjects,
                          unordered_map<Object*, Bundle*>& linePool,
                          const int count)
{
    double rate = computeCombReaction(hostObject, newObject, allObjects, linePool, count);
    std::pair<int64, double> oneReaction(newObject->getKey(), rate);
    secondR.insert(oneReaction);
}

void OneLine::removeReaction(const int64 deleteKey)
{
    secondR.erase(deleteKey);
}

void OneLine::updateReaction(
                             Object const * const hostObject,
                             Object const * const mobileObject,
                             unordered_map<int64, Object*>& allObjects,
                             unordered_map<Object*, Bundle*>& linePool,
                             const int n)
{
    double rate = computeCombReaction(hostObject, mobileObject, allObjects, linePool, n);
    secondR[mobileObject->getKey()] = rate;
}

void OneLine::updateLine(
                         const Object* const hostObject,
                         const int count,
                         unordered_map<int64, Object*>& mobileObjects,
                         unordered_map<int64, Object*>& allObjects,
                         unordered_map<Object*, Bundle*>& linePool)
{
    secondR.clear();
    setOneLine(hostObject, count, mobileObjects, allObjects, linePool);
}

void OneLine::updateDiff(
                        const Object* const hostObject, 
                        const int count,
                        unordered_map<int64, Object*>& allObjects)
{
    computeDiffReaction(hostObject, count, allObjects);
}

const long double OneLine::computeTotalRate()
{
    int i;
    unordered_map<int64,long double>::iterator iter;
    totalRate = 0.0;
    totalRate += diffRToF; /* add one diffusion rate */
    totalRate += diffRToB; /* add another diffusion rate*/
    totalRate += sinkRDislocation;    /* add dislocation sink rate */
    totalRate += sinkRGrainBndry;     /* add grain boundary sink rate */
    totalRate += SAVR;     /* add super abundant vacancy rate */
    totalRate += recombRER; /* add one recombination rate */
    totalRate += recombRLH; /* add another recombination rate */
    for (i = 0; i < LEVELS; i++) {
        totalRate += dissociationR[i];
    }
    for (iter = secondR.begin(); iter != secondR.end(); ++iter) {
        totalRate += iter->second;
    }
    return totalRate;
}

void OneLine::display(Object const * const hostObject)
{
    ofstream fs;
    fs.open("lines.txt", ios::app);
    fs << "Line for Oject" << hostObject->getKey() << ":    ";
    fs << "(diff)" << diffRToF << ", " << diffRToB << "    " << "(sink)" << sinkRDislocation << ", " << sinkRGrainBndry << "    ";
    for (int i = 0; i < LEVELS; ++i) {
        fs <<"(diss)"<< dissociationR[i] << "    ";
    }
    unordered_map<int64, long double>::iterator iter;
    for (iter = secondR.begin(); iter != secondR.end(); ++iter) {
        fs << "(" << iter->first << ")" << iter->second << "    ";
    }
    fs << "(SAV)" << SAVR << "    ";
    fs << "(recomb)" << recombRER << ", " << recombRLH;
    fs << endl;
    fs.close();
}

/* private function implementations */
void OneLine::setOneLine(
                         const Object* const hostObject,
                         const int count,
                         unordered_map<int64, Object*>& mobileObjects,
                         unordered_map<int64, Object*>& allObjects,
                         unordered_map<Object*, Bundle*>& linePool)
{
    computeDiffReaction(hostObject, count, allObjects);
    computeSinkReaction(hostObject, count);
    for (int index = 0; index < LEVELS; index++) {
        dissociationR[index] = computeDissReaction(hostObject, allObjects, linePool, index, count);
    }
    unordered_map<int64, Object*>::iterator iter;
    for (iter = mobileObjects.begin(); iter != mobileObjects.end(); ++iter) {
        double rate = computeCombReaction(hostObject, iter->second, allObjects, linePool, count);
        std::pair<int64, double> oneReaction(iter->first, rate);
        secondR.insert(oneReaction);
    }
    computeSAVReaction(hostObject, count);
    computeRecombReaction(hostObject, count, allObjects);
    computeTotalRate();
}

void OneLine::computeDiffReaction(const Object* const hostObject, const int count, unordered_map<int64, Object*>& allObjects)
{
    if (!DIFF_ON)
    {
        diffRToF = 0.0;
        diffRToB = 0.0;
        return;
    }

    double prefactor = 0.0;
    int objectN[3];   
    hostObject->getThreeNumber(count, objectN);
    double concentration = 0;
    double frontConcentration = 0;
    double backConcentration = 0;
    if (count == SURFACE_INDEX)
    {
        concentration = 0;  /* surface layer corresponds to adsorbed layer of hydrogen at surface */
        frontConcentration = H_SATURATION_CONCENTRATION;
        backConcentration = objectN[2] / SUBSURFACE_VOLUME;
    }
    else if (count == SUBSURFACE_INDEX)
    {
        concentration = objectN[0] / SUBSURFACE_VOLUME;
        frontConcentration = 0;
        backConcentration = objectN[2] / FIRST_BULK_VOLUME;
    }
    else if (count == FIRST_BULK_INDEX)
    {
        concentration = objectN[0] / FIRST_BULK_VOLUME;
        frontConcentration = objectN[1] / SUBSURFACE_VOLUME;
        backConcentration = objectN[2] / VOLUME;
    }
    else if (count == FIRST_BULK_INDEX + 1)
    {
        concentration = objectN[0] / VOLUME;
        frontConcentration = objectN[1] / FIRST_BULK_VOLUME;
        backConcentration = objectN[2] / VOLUME;
    }
    else
    {
        concentration = objectN[0] / VOLUME;
        frontConcentration = objectN[1] / VOLUME;
        backConcentration = objectN[2] / VOLUME;
    }
    const int64 SIAKey = 1000000;

    // Account for special cases from 2020 Zhenhou Wang for hydrogen moving between surface and bulk
    if((count == SUBSURFACE_INDEX && hostObject->getKey() == 1) ||
        (count == SURFACE_INDEX && hostObject->getKey() == 1))
    {
        double surfaceConc = 0.0; 
        int64 HKey = 1;
        if (allObjects.find(HKey) != allObjects.end())
            surfaceConc = allObjects[HKey]->getNumber(0) / DIVIDING_AREA;  // [cm^-2] concentration

        double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0);
        double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;

        // special case for 1H diffusion from Subsurface to Surface 
        if (count == SUBSURFACE_INDEX)
        {
            double jumpingDist = maxSurfaceConc / 6 / DENSITY;
            double freq = NU0 * exp(-H_MIGRATION_ENERGY / KB / TEMPERATURE);
        
            prefactor = freq * jumpingDist * (1 - surfaceSaturationFraction) * DIVIDING_AREA;
            diffRToF = prefactor * concentration;
            if (surfaceConc >= maxSurfaceConc)
                diffRToF = 0.0;
            diffRToB = 0.0;

            // Normal Diffusion to the first bulk element
            if (concentration > backConcentration) // allow 1H to diffuse into neigbouring volume element if both vol elems have 1H
            {
                double lengthb = (SUBSURFACE_THICKNESS + FIRST_BULK_THICKNESS)/2.0 * NM_TO_CM; // first element to second element distance (20nm) 
                /* if diffusable */
                prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthb;
                diffRToB = prefactor*(concentration - backConcentration);
            }
            return;
        }
        // special case for 1H diffusion from Surface to Subsurface
        else if (count == 0)
        {
            double absorbE;
            // double desorbE = -0.00195416 * exp(5.87242*surfaceSaturationFraction) + 1.48996;            // Ajmalghan 2019
            // double desorbE = 2.0*(0.9 - 0.2*surfaceSaturationFraction - 0.7*pow(surfaceSaturationFraction, 12));
            // double desorbE = 1.023 + 0.584/(1.0 + exp(7.38e-16 * surfaceConc - 2.85));
            // double desorbE = 2.0*(0.525 + 0.591*(1.0/(1.0+exp( (surfaceSaturationFraction-0.247)/0.0692 )))); // from Hodille 2020
            // if (TEMP_INCREASE_RATE == 0)  // no thermal desorption, assume atmosphere environment so use experiment data
                absorbE = 1.10 + 0.939*(1.0/(1.0+exp( (surfaceSaturationFraction-0.232)/0.0683 )));  // from Hodille 2020
            // else                          // doing thermal desorption, assume vacuum environment so use DFT data
                // absorbE = -3.6592e-8 * exp(16.9129*surfaceSaturationFraction) + 1.71738;             // Ajmalghan 2019
                // absorbE = desorbE/2. + HEAT_OF_SOLUTION + H_MIGRATION_ENERGY + 0.02;   // Add 0.02 from Tajuki Oda 2023
            double freq = NU0 * exp(-absorbE / KB / TEMPERATURE);            
            prefactor = freq * surfaceConc * DIVIDING_AREA;
            diffRToB = prefactor;
            diffRToF = 0.0;
            return;
        }
    }
    // SIA moving between surface and bulk (assume it has the same equations as hydrogen)
    else if (count == SUBSURFACE_INDEX && hostObject->getKey() == SIAKey)
    {
        double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0);
        double jumpingDist = maxSurfaceConc / 6 / DENSITY;
        double siaMigrationEnergy = 0.009; // [eV]
        double freq = NU0 * exp(-siaMigrationEnergy / KB / TEMPERATURE);
    
        prefactor = freq * jumpingDist * DIVIDING_AREA;
        diffRToF = prefactor * concentration;
        diffRToB = 0.0;

        // Normal Diffusion to the first bulk element
        if (concentration > backConcentration)
        {
            double lengthb = (SUBSURFACE_THICKNESS + FIRST_BULK_THICKNESS)/2.0 * NM_TO_CM; // first element to second element distance (20nm) 
            /* if diffusable */
            prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthb;
            diffRToB = prefactor*(concentration - backConcentration);
        }
        return;
    }
    else // avoid unnecessary calculation if we are processing special cases
    {
        /* by having two diffusion rates, this rate will never be less than 0 */
        /* length measured in cm */
        double lengthf = 0.0, lengthb = 0.0;
        if (count == SURFACE_INDEX)
        {   
            lengthf = (ELEMENT_THICKNESS + SUBSURFACE_THICKNESS) / 2. * NM_TO_CM; // this value doesn't matter b/c can't diffuse out into vacuum 
            lengthb = (SUBSURFACE_THICKNESS) / 2. * NM_TO_CM; /* surface to first element distance */
        }
        else if (count == SUBSURFACE_INDEX)
        {
            lengthf =  SUBSURFACE_THICKNESS / 2. * NM_TO_CM; /* surface to first element distance */
            lengthb = (SUBSURFACE_THICKNESS + FIRST_BULK_THICKNESS) / 2. * NM_TO_CM;  
        }
        else if (count == FIRST_BULK_INDEX)
        {
            lengthf = (SUBSURFACE_THICKNESS + FIRST_BULK_THICKNESS) / 2. * NM_TO_CM;
            lengthb = (FIRST_BULK_THICKNESS + ELEMENT_THICKNESS) / 2. * NM_TO_CM;
        }
        else if (count == FIRST_BULK_INDEX + 1)
        {
            lengthf = (FIRST_BULK_THICKNESS + ELEMENT_THICKNESS) / 2. * NM_TO_CM;
            lengthb = ELEMENT_THICKNESS * NM_TO_CM;
        }
        else
        {
            lengthf = lengthb = ELEMENT_THICKNESS * NM_TO_CM; /* other element distances */
        }

        /* 
        * 1. compute diffusion rate to the front element 
        * Diffusion goes from area of higher concentration to lower concentration
        * Object not allowed to diffuse out through the front
        * Allow neighbouring elements with 1H each to diffuse into each other like real life.
        */
        if ((concentration > frontConcentration)
             && count != SURFACE_INDEX && count != SUBSURFACE_INDEX) 
        {
            /* if diffusable, surface objects diffusing into vacuum is considered */
            prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthf;
            diffRToF = prefactor*(concentration - frontConcentration);
        }
        else {
            diffRToF = 0.0;
        }
        /* 2. compute diffusion rate to the back element
        * Object are allowed to diffuse out through the back (assume infinite W sample)
        */
        if (concentration > backConcentration && count != SURFACE_INDEX) {
            /* if diffusable */
            prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthb;
            diffRToB = prefactor*(concentration - backConcentration);
        }
        else {
            diffRToB = 0.0;
        }
        if (TEMP_INCREASE_RATE > 0 && count == POINTS - 1) { // if doing TDS, don't let stuff escape thru the back, so we can count it as it emerges from surface
            diffRToB = 0.0;
        }
        if (hostObject->getAttri(0) > 0 && hostObject->getAttri(2) == 0 && count == POINTS - 1) {
            diffRToB = 0.0; // Don't let SIA diffuse out of box so it can recombine with vacancies at the edge
        }
    }
}

void OneLine::computeSinkReaction(const Object* const hostObject, const int count)
{
    if (!SINK_ON || count == SURFACE_INDEX || count == SUBSURFACE_INDEX)
    {
        sinkRDislocation = 0.0;
        sinkRGrainBndry = 0.0;
        return;
    }

    sinkRDislocation = hostObject->getNumber(count)*hostObject->getDiff()*hostObject->getSinkDislocation();
    sinkRGrainBndry = hostObject->getNumber(count)*hostObject->getDiff()*hostObject->getSinkGrainBndry();
}

long double OneLine::computeBaseDissReaction(
                                  const Object* const hostObject,
                                  const int index,
                                  const int count) const
{
    if (!DISS_ON || count == SURFACE_INDEX || count == SUBSURFACE_INDEX)
    {
        return 0.0;
    }


    // All types of monomer (I1, V1, H1, He1 object, etc.) should have 0 dissociation rate
    int elementNum = 0;
    for (int i = 0; i < LEVELS; i++)
    {  
        elementNum += abs(hostObject->getAttri(i));
    }

    if (hostObject->getAttri(index) != 0 && elementNum > 1) {
        return jumped / (jumped + hostObject->getR1e()) * 4.0 * PI * pow(hostObject->getR1e(), 2) / pow(ALATT, 2) * NU0 * hostObject->getBind(index) * hostObject->getNumber(count);
    }
    return 0.0;
}

long double OneLine::computeDissReaction(
                                  const Object* const hostObject,
                                  unordered_map<int64, Object*>& allObjects,
                                  unordered_map<Object*, Bundle*>& linePool,
                                  const int index,
                                  const int count) const
{
    /* Compute the net diss reaction rate for this object and comb reaction rate for its predecessor conjugate for approximation speedup, this updates the predecessor's comb rate as well 
       For now it is only enabled for VH clusters, because chopping off a diss pathway for those doesn't seem to matter much */   
    long double baseDissRate = computeBaseDissReaction(hostObject, index, count);
    return baseDissRate;
    long double netDissRate = baseDissRate;
    int attrIndexH = 2;
    int64 HKey = 1;
    if (index == attrIndexH && 
        ((hostObject->getAttri(0) <= -1 && hostObject->getAttri(2) >= 7) ||
         (hostObject->getAttri(0) >= 20 && hostObject->getAttri(2) >= 41)) )
    {
        int predAttr[LEVELS]; // predecessor attributes
        for (int level = 0; level < LEVELS; level++)
        {
            predAttr[level] = hostObject->getAttri(level);
        }
        predAttr[attrIndexH]--;  // for an H dissociation
        int64 predKey = attrToKey(predAttr);

        if (allObjects.find(predKey) != allObjects.end() && allObjects[predKey]->getNumber(count) > 0 &&
            allObjects.find(HKey) != allObjects.end() && allObjects[HKey]->getNumber(count) > 0)
        {
            Object* predObj = allObjects[predKey];
            Object* HObj = allObjects[HKey];
            OneLine* predLine = linePool[predObj]->lines[count];
            long double baseCombRate = computeBaseCombReaction(predObj, HObj, count);
            if (baseCombRate > baseDissRate)
            {
                if (predLine != nullptr)
                    predLine->setCombReaction(HKey, baseCombRate - baseDissRate);
                netDissRate = 0.0;
            }
            else if (baseDissRate > baseCombRate)
            {
                if (predLine != nullptr)
                    predLine->setCombReaction(HKey, 0.0);
                netDissRate = baseDissRate - baseCombRate;
            }
        }
        else
        {
            Object* predObj;
            Object* HObj;
            if (allObjects.find(predKey) == allObjects.end() || allObjects[predKey]->getNumber(count) <= 0)
            {
                predObj = new Object(predKey, count);
            }
            else
            {
                predObj = allObjects[predKey];
            }
            if (allObjects.find(HKey) == allObjects.end() || allObjects[HKey]->getNumber(count) <= 0)
            {
                HObj = new Object(HKey, count);
            }
            else
            {
                HObj = allObjects[HKey];
            }

            long double baseCombRate = computeBaseCombReaction(predObj, HObj, count);
            if (baseCombRate > baseDissRate)
            {
                netDissRate = 0.0;
            }

            if (allObjects.find(predKey) == allObjects.end() || allObjects[predKey]->getNumber(count) <= 0)
            {
                delete predObj;
            }   
            if (allObjects.find(HKey) == allObjects.end() || allObjects[HKey]->getNumber(count) <= 0)
            {
                delete HObj;
            }
        }
    }

    return netDissRate;
}

long double OneLine::computeBaseCombReaction(
                                    const Object* const hostObject,
                                    const Object* const mobileObject,
                                    const int count) const
{    
    if (!COMB_ON || count == SURFACE_INDEX || count == SUBSURFACE_INDEX)
    {
        return 0.0;
    }

    double concentration;
    double r12;
    double dimensionTerm;
    double volume;
    double adjustmentFactor = 1;
    if (count == SURFACE_INDEX || count == SUBSURFACE_INDEX)
        volume = SUBSURFACE_VOLUME; // shouldn't get here
    else if (count == FIRST_BULK_INDEX)
        volume = FIRST_BULK_VOLUME;
    else
        volume = VOLUME;
    
    if (hostObject->getKey() != mobileObject->getKey()) {
        concentration = hostObject->getNumber(count)*mobileObject->getNumber(count) / volume;
        
    }else {
        concentration = hostObject->getNumber(count)*(hostObject->getNumber(count) - 1) / volume;
        
    }

    // H+H-->2H
    // Disable H clustering for now to reduce complexity
    if (hostObject->getKey() == 1 && mobileObject->getKey() == 1){
        return 0.0;
    }

    // Disable multiples of 1V-12H + 1H because vacancy can store max 12H, save sim time
    if (hostObject->getAttri(0) < 0 && 
        hostObject->getAttri(2) == abs(hostObject->getAttri(0))*12 &&
        mobileObject->getAttri(0) == 0 && mobileObject->getAttri(2) > 0)
    {
        return 0.0;
    }

    // Harder for SIA to recombine with V-H cluster if H/V ratio is high
    if (hostObject->getAttri(0) < 0 &&
        hostObject->getAttri(2) > 0 &&
        mobileObject->getAttri(0) > 0 &&
        mobileObject->getAttri(2) == 0)
    {
        // Assume that when we reach the max H surface density inside vacancy cluster, SIA can't recombine with VH cluster
        int numH = hostObject->getAttri(2);
        int monovacancyMaxH = 12;

        double radius = hostObject->getR1();
        double surfArea = 4*PI*pow(radius, 2);
        double surfHDensity = numH / surfArea;

        double monovacancyRadius = pow(3.0*avol/4.0/PI, 1.0/3.0);
        double monovacancySurfArea = 4*PI*pow(monovacancyRadius, 2);
        double maxSurfHDensity = monovacancyMaxH / monovacancySurfArea;

        adjustmentFactor = 1 - surfHDensity / maxSurfHDensity;

        if (adjustmentFactor < 0)
            adjustmentFactor = 0;
    }

    /*
    if(hostObject->getKey() == 1 && mobileObject->getKey() == 2){
        return 0.0;
    }
    if(hostObject->getKey() == 2 && mobileObject->getKey() == 1){
        return 0.0;
    }
    */

    r12 = hostObject->getR1() + mobileObject->getR1();
    dimensionTerm = computeDimensionTerm(r12, hostObject, mobileObject, count);
    return 4.0*PI*concentration*r12*dimensionTerm*adjustmentFactor;
}

long double OneLine::computeCombReaction(
                                    const Object* const hostObject,
                                    const Object* const mobileObject,
                                    unordered_map<int64, Object*>& allObjects,
                                    unordered_map<Object*, Bundle*>& linePool,
                                    const int count) const
{
    /* Compute the net comb/diss reaction rate for approximation speedup, this updates the product's diss rate as well 
       For now it is only enabled for VH clusters, because chopping off a diss pathway for those doesn't seem to matter much */
    long double baseCombRate = computeBaseCombReaction(hostObject, mobileObject, count);
    return baseCombRate;
    long double netCombRate = baseCombRate;
    int attrIndexH = 2;
    if (((hostObject->getAttri(0) <= -1 && hostObject->getAttri(2) >= 6) ||
        (hostObject->getAttri(0) >= 20 && hostObject->getAttri(2) >= 40)) && 
        mobileObject->getAttri(0) == 0 && mobileObject->getAttri(2) == 1)
    {
        int prodAttr[LEVELS];
        for (int level = 0; level < LEVELS; level++)
        {
            prodAttr[level] = hostObject->getAttri(level) + mobileObject->getAttri(level);
        }
        int64 prodKey = attrToKey(prodAttr);

        if (allObjects.find(prodKey) != allObjects.end() && allObjects[prodKey]->getNumber(count) > 0)
        {
            Object* prodObj = allObjects[prodKey];
            OneLine* prodLine = linePool[prodObj]->lines[count];
            long double baseDissRate = computeBaseDissReaction(prodObj, attrIndexH, count);
            if (baseCombRate > baseDissRate)
            {
                if (prodLine != nullptr)
                    prodLine->setDissReaction(attrIndexH, 0.0);
                netCombRate = baseCombRate - baseDissRate;
            }
            else if (baseDissRate > baseCombRate)
            {
                if (prodLine != nullptr)
                    prodLine->setDissReaction(attrIndexH, baseDissRate - baseCombRate);
                netCombRate = 0.0;
            }
        }
        else
        {
            /* Product doesn't exit */
            Object* tempProdObj = new Object(prodKey, count);
            long double baseDissRate = computeBaseDissReaction(tempProdObj, attrIndexH, count);
            if (baseDissRate > baseCombRate)
            {
                netCombRate = 0.0;
            }
            delete tempProdObj;
        }
    }

    return netCombRate;
}

void OneLine::computeSAVReaction(
                                 const Object* const hostObject,
                                 const int count)
{
    /*
     * Allow overpressurized VH cluster to eject W atom to create another vacancy.
     * And allow excess 1H to eject W atom when H is oversaturated.
     */
    SAVR = 0;
    double volume = VOLUME;
    if (count == SUBSURFACE_INDEX)
        volume = SUBSURFACE_VOLUME;
    else if (count == FIRST_BULK_INDEX)
        volume = FIRST_BULK_VOLUME;

    if (!SAV_ON || count == SURFACE_INDEX || count == SUBSURFACE_INDEX)
    {
        return;
    }

    // If we have a mV-nH object, or nH object
    if (hostObject->getAttri(0) <= 0 && hostObject->getAttri(2) > 0 && hostObject->getNumber(count) > 0)
    {
        int numHPerCluster = hostObject->getAttri(2);
        int numVacancies = abs(hostObject->getAttri(0));
        double clusterThresholdH = 4.0*numVacancies;   // Qianran Yu 2020, did linear fit from graph of excess sav energies
        if (numHPerCluster > clusterThresholdH)
        {
            // 1H is SAV candidate only if dissolved H concentration is oversaturated
            if (hostObject->getKey() == 1) 
            {
                double maxSolubilityConc = DENSITY * exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE);
                double maxNumH = maxSolubilityConc * volume;
                int extraH = ceil(hostObject->getNumber(count) - maxNumH);
                if (extraH > 0)
                    SAVR = NU0 * exp(-SAV_ENERGY/KB/TEMPERATURE) * extraH;
            }
            // Overpressurized VH cluster is always SAV candidate
            else
                SAVR = NU0 * exp(-SAV_ENERGY/KB/TEMPERATURE) * hostObject->getNumber(count);
        }
    }
}

void OneLine::computeRecombReaction(
                                    const Object* const hostObject,
                                    const int count,
                                    unordered_map<int64, Object*>& allObjects)
{
    recombRER = 0.0;
    recombRLH = 0.0;

    // only H can recombine at surface and leave surface
    if (!RECOMB_ON || count != SURFACE_INDEX || hostObject->getKey() != 1)
    {
        return;
    }

    double surfaceConc = 0.0;
    int64 HKey = 1;
    int numH = allObjects[HKey]->getNumber(0);
    if (allObjects.find(HKey) != allObjects.end())
        surfaceConc = numH / DIVIDING_AREA;  // [cm^-2] concentration

    double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0);
    double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;

    // Calculate ER recomb rate
    if (numH >= 1 && HYDROGEN_ON)   // incident atom collides with adsorbed atom, so hydrogen must be on for this to work
    {
        double crossSectionERRecomb = 1.7e-17; // [cm^2] cross-section of ER recombination from Zhenhou Wang 2020
        recombRER = FLUX_H * crossSectionERRecomb * surfaceConc * DIVIDING_AREA; 
    }
    else
        recombRER = 0.0;
    
    // Calculate LH recomb rate
    if (numH >= 2)                  // two atoms combine at the surface to form H2 and desorb
    {
        double desorbE;
        // if (TEMP_INCREASE_RATE == 0)  // assume atmospheric environment, use experiment desorption energy barrier
            // desorbE = 2.0*(0.525 + 0.591*(1.0/(1.0+exp( (surfaceSaturationFraction-0.247)/0.0692 )))); // from Hodille 2020
        // else                          // doing thermal desorption, assume vacuum environment, use DFT desorption energy barrier
            // desorbE = -0.00195416 * exp(5.87242*surfaceSaturationFraction) + 1.48996;            // Ajmalghan 2019
            // desorbE = 1.40259 - 0.00881176*exp(5.45029*surfaceSaturationFraction - 1.22515);
            // desorbE = 0.019+1.453/(1.0+exp((surfaceSaturationFraction-1.000)/0.111));
            // desorbE = 2.0*(0.9 - 0.2*surfaceSaturationFraction - 0.7*pow(surfaceSaturationFraction, 12));
            // desorbE = 1.023 + 0.584/(1.0 + exp(7.38e-16 * surfaceConc - 2.85));
            desorbE = 1.029 + 0.700/(1.0+exp((surfaceSaturationFraction-0.475)/0.151));
        double desorptionR = NU0 * pow(ALATT, 2); // [cm^2 s^-1]
        recombRLH = desorptionR * exp(-desorbE / KB / TEMPERATURE) * surfaceConc * surfaceConc * DIVIDING_AREA;
    }
    else
        recombRLH = 0.0;
}

double OneLine::computeDimensionTerm(
                                     const double r12,
                                     const Object* const hostObject,
                                     const Object* const mobileObject,
                                     const int count) const
{
    double term = 0.0;
    double hostDiff = hostObject->getDiff(), mobileDiff = mobileObject->getDiff();
    int hostDim = hostObject->getDim(), mobileDim = mobileObject->getDim();
    int hostN = hostObject->getNumber(count), mobileN = mobileObject->getNumber(count);
    // int dimsum = hostDim + mobileDim; 
    int dimsum = 6;
    double alpha_a = -log(PI*PI*pow(r12, 3.0) / VOLUME / hostN);
    double alpha_b = -log(PI*PI*pow(r12, 3.0) / VOLUME / mobileN);
    switch (dimsum) {
        case 6: // 3D + 3D
            term = hostDiff + mobileDiff;
            break;
        case 4: // 3D + 1D
            if (hostDim == 1 && mobileDim == 3)
                term = hostObject->getDiff()*(mobileN / VOLUME)*(2.0*PI*pow(r12, 3.0)) + mobileDiff;
            else
                term = mobileDiff*(mobileN / VOLUME)*(2.0*PI*pow(r12, 3.0)) + hostDiff;
            break;
        case 2: // 1D + 1D
            term = hostDiff / alpha_b + mobileDiff / alpha_a;
            break;
    }
    return term;
}

void OneLine::setDissReaction(const int index, long double rate)
{
    dissociationR[index] = rate;
}

void OneLine::setCombReaction(const int64 mobileObjectKey, long double rate)
{
    secondR[mobileObjectKey] = rate;
}