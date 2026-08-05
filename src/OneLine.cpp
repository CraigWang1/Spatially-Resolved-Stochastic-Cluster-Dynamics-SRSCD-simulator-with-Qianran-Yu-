// OneLine.cpp -- implementations of class OneLine
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include "OneLine.h"
#include "Bundle.h"
#include "rvgs.h"
using namespace std;

/* public function implementations */
/* Helper class to calculate reaction rates of a species in one spatial element */
OneLine::OneLine(
                 const Object* const hostObject,
                 const int count,
                 unordered_map<int64, Object*>& mobileObjects,
                 unordered_map<int64, Object*>& allObjects) :totalRate(0.0)
{
    setOneLine(hostObject, count, mobileObjects, allObjects);
}

OneLine::OneLine() : diffRToF(0.0), diffRToB(0.0), sinkRDislocationScrew(0.0), sinkRDislocationEdge(0.0), sinkRGrainBndry(0.0), SAVR(0.0), recombRER(0.0), recombRLH(0.0), totalRate(0.0)
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
    if (sinkRDislocationScrew >= tempRate) {
        return SINKDISLOCATIONSCREW;
    }
    else {
        tempRate -= sinkRDislocationScrew;
    }
    if (sinkRDislocationEdge >= tempRate) {
        return SINKDISLOCATIONEDGE;
    }
    else {
        tempRate -= sinkRDislocationEdge;
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
                          const int count)
{
    double rate = computeCombReaction(hostObject, newObject, allObjects, count);
    if (rate > 0)
    {
        std::pair<int64, double> oneReaction(newObject->getKey(), rate);
        secondR.insert(oneReaction);
    }
    else
    {
        secondR.erase(newObject->getKey());
    }
}

void OneLine::removeReaction(const int64 deleteKey)
{
    secondR.erase(deleteKey);
}

void OneLine::updateReaction(
                             Object const * const hostObject,
                             Object const * const mobileObject,
                             unordered_map<int64, Object*>& allObjects,
                             const int n)
{
    double rate = computeCombReaction(hostObject, mobileObject, allObjects, n);
    if (rate > 0)
        secondR[mobileObject->getKey()] = rate;
    else
        secondR.erase(mobileObject->getKey());
}

void OneLine::updateLine(
                         const Object* const hostObject,
                         const int count,
                         unordered_map<int64, Object*>& mobileObjects,
                         unordered_map<int64, Object*>& allObjects)
{
    secondR.clear();
    setOneLine(hostObject, count, mobileObjects, allObjects);
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
    totalRate += sinkRDislocationScrew;    /* add dislocation sink rate */
    totalRate += sinkRDislocationEdge;
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
    fs << "(diff)" << diffRToF << ", " << diffRToB << "    " << "(sink)" << sinkRDislocationScrew << ", " << sinkRDislocationEdge << ", " << sinkRGrainBndry << "    ";
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
                         unordered_map<int64, Object*>& allObjects)
{
    computeDiffReaction(hostObject, count, allObjects);
    computeSinkReaction(hostObject, count);
    for (int index = 0; index < LEVELS; index++) {
        dissociationR[index] = computeDissReaction(hostObject, allObjects, index, count);
    }
    unordered_map<int64, Object*>::iterator iter;
    for (iter = mobileObjects.begin(); iter != mobileObjects.end(); ++iter) {
        double rate = computeCombReaction(hostObject, iter->second, allObjects, count);
        if (rate > 0)
        {
            std::pair<int64, double> oneReaction(iter->first, rate);
            secondR.insert(oneReaction);
        }
        else
        {
            secondR.erase(iter->second->getKey());
        }
    }
    computeSAVReaction(hostObject, count);
    computeRecombReaction(hostObject, count, allObjects);
    computeTotalRate();
}

void OneLine::computeDiffReaction(const Object* const hostObject, const int count, unordered_map<int64, Object*>& allObjects)
{
    diffRToF = 0.0;
    diffRToB = 0.0;

    if (!DIFF_ON
        || (UNIFORM_FREE_H_ON && hostObject->getKey() == 1))
    {
        // For uniform concentration of free hydrogen mode,
        // don't need to simulate diffusion because there is assuemd to be
        // uniform, constant concentration
        return;
    }

    long double prefactor = 0.0;
    int objectN[3];   
    hostObject->getThreeNumber(count, objectN);
    double concentration = 0;
    double frontConcentration = 0;
    double backConcentration = 0;
    double vol = volumeAtIndex(count);
    double volf = volumeAtIndex(count-1);
    double volb = volumeAtIndex(count+1);
    if (vol != 0)
        concentration = objectN[0] / vol;
    if (volf != 0)
        frontConcentration = objectN[1] / volf;
    if (volb != 0)
        backConcentration = objectN[2] / volb;
    if (count == SURFACE_INDEX)
    {
        concentration = 0;  /* surface layer corresponds to adsorbed layer of hydrogen at surface, which doesn't follow standard diffusion (see below) */
        frontConcentration = 0;
    }
    else if (count == SUBSURFACE_INDEX)
    {
        frontConcentration = 0;
    }
    else if (count == BACK_SUBSURFACE_INDEX && BACK_DESORB)
    {
        backConcentration = 0;
    }
    else if (count == BACK_SURFACE_INDEX && BACK_DESORB)
    {
        concentration = 0;
        backConcentration = 0;
    }

    /* length measured in cm */
    double distf = lengthf(count);
    double distb = lengthb(count);

    // Account for special cases from 2020 Zhenhou Wang for hydrogen moving between surface and bulk
    if ((
        count == SURFACE_INDEX
        || count == SUBSURFACE_INDEX
        || (count == BACK_SUBSURFACE_INDEX && BACK_DESORB)
        || (count == BACK_SURFACE_INDEX && BACK_DESORB)
        ) 
        && hostObject->getKey() == 1)
    {
        double surfaceConc = 0.0; 
        int64 HKey = 1;
        if (allObjects.find(HKey) != allObjects.end())
        {
            if (count == SURFACE_INDEX || count == SUBSURFACE_INDEX)
                surfaceConc = allObjects[HKey]->getNumber(SURFACE_INDEX) / DIVIDING_AREA;  // [cm^-2] concentration
            else
                surfaceConc = allObjects[HKey]->getNumber(BACK_SURFACE_INDEX) / DIVIDING_AREA;  // [cm^-2] concentration
        }

        double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0); // [110] surface
        double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;

        // special case for 1H diffusion from Subsurface to Surface 
        if (count == SUBSURFACE_INDEX || count == BACK_SUBSURFACE_INDEX)
        {
            double jumpingDist = maxSurfaceConc / 6 / DENSITY;
            double freq = NU0 * hostObject->getExpMig();
            long double diffRToSurf = 0;
            diffRToF = diffRToB = 0;
        
            prefactor = freq * jumpingDist * (1 - surfaceSaturationFraction) * DIVIDING_AREA;
            diffRToSurf = prefactor * concentration;
            if (surfaceConc >= maxSurfaceConc)
                diffRToSurf = 0.0;

            if (count == SUBSURFACE_INDEX)
            {
                diffRToF = diffRToSurf;
                if (concentration > backConcentration)
                {
                    /* if diffusable */
                    prefactor = hostObject->getDiff() * DIVIDING_AREA / distb;
                    diffRToB = prefactor*(concentration - backConcentration);
                }
            }
            else if (count == BACK_SUBSURFACE_INDEX)
            {
                diffRToB = diffRToSurf;
                if (concentration > frontConcentration)
                {
                    prefactor = hostObject->getDiff() * DIVIDING_AREA / distf;
                    diffRToF = prefactor*(concentration - frontConcentration);
                }
            }
            return;
        }
        // special case for 1H diffusion from Surface to Subsurface
        else if (count == SURFACE_INDEX || count == BACK_SURFACE_INDEX)
        {
            diffRToF = diffRToB = 0;
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
            double freq = NU0 * exp(-absorbE / (KB * TEMPERATURE));            
            prefactor = freq * surfaceConc * DIVIDING_AREA;

            if (count == SURFACE_INDEX)
                diffRToB = prefactor;
            else if (count == BACK_SURFACE_INDEX)
                diffRToF = prefactor;

            return;
        }
    }
    else // avoid unnecessary calculation if we are processing special cases
    {
        /* 
         * If the number of objects is small, 
         * transition from net diffusion rate
         * to a (non-net) hop rate to capture 
         * realistic local concentration spikes
         */
        const int hopThres = 10;
        bool hopFront = HOP_ON && (objectN[0] < hopThres) && (objectN[1] < hopThres);
        bool hopBack = HOP_ON && (objectN[0] < hopThres) && (objectN[2] < hopThres); 

        diffRToF = diffRToB = 0;

        /* 
        * 1. compute diffusion rate to the front element 
        * Diffusion goes from area of higher concentration to lower concentration
        * Object not allowed to diffuse out through the front
        */
        if ((concentration > frontConcentration || hopFront)
            && count != SURFACE_INDEX
            && (count != BACK_SURFACE_INDEX || !BACK_DESORB)  
            && (count != SUBSURFACE_INDEX || (hostObject->getAttri(0) != 0 && hostObject->getAttri(2) == 0))) 
        {
            prefactor = hostObject->getDiff() * DIVIDING_AREA / distf;
            if (hopFront)
            {
                diffRToF = prefactor*concentration;
            }
            else
            {
                diffRToF = prefactor*(concentration - frontConcentration);
            }
        }

        /* 
         * 2. compute diffusion rate to the back element
         */
        if ((concentration > backConcentration || hopBack) 
            && count != SURFACE_INDEX 
            && count != BACK_SURFACE_INDEX
            && (count != BACK_SUBSURFACE_INDEX || !BACK_DESORB || (hostObject->getAttri(0) != 0 && hostObject->getAttri(2) == 0))) 
        {
            prefactor = hostObject->getDiff() * DIVIDING_AREA / distb;
            if (hopBack)
            {
                diffRToB = prefactor*concentration;
            }
            else
            {
                diffRToB = prefactor*(concentration - backConcentration);
            }
        }
    }
}

void OneLine::computeSinkReaction(const Object* const hostObject, const int count)
{
    if (!SINK_ON 
        || count == SURFACE_INDEX 
        || count == SUBSURFACE_INDEX
        || (count == BACK_SUBSURFACE_INDEX && BACK_DESORB)
        || (count == BACK_SURFACE_INDEX && BACK_DESORB)
        || (UNIFORM_FREE_H_ON && hostObject->getKey() == 1)) // Assume that if you're in UNIFORM_FREE_H_ON mode, that since your temperature was high enough such that you assumed the H could penetrate to the back of the sample relatively quickly, that the temperature is high enough that the amount of H trapped in intrinsic traps (dislocations, grain boundaries) is negligible. If this is not the case, comment out this line.
    {
        sinkRDislocationScrew = 0.0;
        sinkRDislocationEdge = 0.0;
        sinkRGrainBndry = 0.0;
        return;
    }

    double hostNumber = hostObject->getNumber(count);
    if (UNIFORM_FREE_H_ON && hostObject->getKey() == 1)
        hostNumber = UNIFORM_H_CONCENTRATION*volumeAtIndex(count);

    sinkRDislocationScrew = hostNumber*hostObject->getDiff()*hostObject->getSinkDislocation() * (1 - EDGE_DISLOCATION_FRAC);
    sinkRDislocationEdge = hostNumber*hostObject->getDiff()*hostObject->getSinkDislocation() * EDGE_DISLOCATION_FRAC;
    sinkRGrainBndry = hostNumber*hostObject->getDiff()*hostObject->getSinkGrainBndry();
}

long double OneLine::computeBaseDissReaction(
                                  const Object* const hostObject,
                                  const int index,
                                  const int count) const
{
    if (!DISS_ON 
        || count == SURFACE_INDEX 
        || count == SUBSURFACE_INDEX
        || (count == BACK_SUBSURFACE_INDEX && BACK_DESORB)
        || (count == BACK_SURFACE_INDEX && BACK_DESORB))
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
            OneLine* predLine = predObj->lines[count];
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
    if (hostObject->getNumber(count) * mobileObject->getNumber(count) == 0)
        return 0;

    if (!COMB_ON 
        || count == SURFACE_INDEX 
        || count == SUBSURFACE_INDEX
        || (count == BACK_SUBSURFACE_INDEX && BACK_DESORB)
        || (count == BACK_SURFACE_INDEX && BACK_DESORB))
    {
        return 0.0;
    }

    double concentration;
    double r12;
    double dimensionTerm;
    double volume = volumeAtIndex(count);
    double adjustmentFactor = 1;
    
    double hostNumber = hostObject->getNumber(count);
    double mobileNumber = mobileObject->getNumber(count);
    if (UNIFORM_FREE_H_ON && hostObject->getKey() == 1)
        hostNumber = UNIFORM_H_CONCENTRATION*volume;
    else if (UNIFORM_FREE_H_ON && mobileObject->getKey() == 1)
        mobileNumber = UNIFORM_H_CONCENTRATION*volume;

    if (hostObject->getKey() != mobileObject->getKey()) {
        concentration = hostNumber*mobileNumber / volume;
    }else {
        concentration = hostNumber*(hostNumber - 1) / volume;
    }

    // H+H-->2H
    // Disable H clustering for now to reduce complexity
    if (hostObject->getKey() == 1 && mobileObject->getKey() == 1 && !NH_CLUSTERING){
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

    // if ((hostObject->getAttri(0) >= 3 && mobileObject->getAttri(0) > 0)
    //     || (hostObject->getAttri(0) > 0 && mobileObject->getAttri(0) >= 3))
    //     return 0;

    // if (mobileObject->getKey() == 2000000 ||
        // (mobileObject->getKey() == 2000000 && hostObject->getKey() == -1000000)){
    // if (mobileObject->getKey() == 2000000){
    //     return 8*PI*mobileObject->getDiff()*pow(r12, 2.0)*mobileObject->getNumber(count)/volume*pow(hostObject->getNumber(count)/volume, 4.0/3.0)*volume;
    // }

    // if (hostObject->getKey() == 2000000 && mobileObject->getKey() == -1000000){
    //     return 8*PI*hostObject->getDiff()*pow(r12, 2.0)*hostObject->getNumber(count)/volume*pow(mobileObject->getNumber(count)/volume, 4.0/3.0)*volume;
    // }

    
    // 1D + 1D
    if (hostObject->getDim() == 1 && mobileObject->getDim() == 1)
    {
        return 8*r12*(hostObject->getDiff()+mobileObject->getDiff())*(r12*pow(mobileObject->getNumber(count)/volume, 1.0/3)+1.0/(2*log(pow(mobileObject->getNumber(count)/volume, -1.0/3)/2.0/r12)))*concentration;
    }

    // 1D + 3D
    else if (hostObject->getDim() == 1 && mobileObject->getDim() == 3)
    {
        double rate_3D_to_0D = 4*PI*r12*mobileObject->getDiff()*concentration;
        double rate_1D_to_0D = 8*PI*hostObject->getDiff()*pow(r12, 2)*hostObject->getNumber(count)/volume*pow(mobileObject->getNumber(count)/volume, 4.0/3)*volume;
        
        return rate_3D_to_0D + rate_1D_to_0D;
    }

    // 3D + 1D
    else if (hostObject->getDim() == 3 && mobileObject->getDim() == 1)
    {
        double rate_3D_to_0D = 4*PI*r12*hostObject->getDiff()*concentration;
        double rate_1D_to_0D = 8*PI*mobileObject->getDiff()*pow(r12, 2)*mobileObject->getNumber(count)/volume*pow(hostObject->getNumber(count)/volume, 4.0/3)*volume;
        
        return rate_3D_to_0D + rate_1D_to_0D;
    }

    // 3D + 3D
    else
    {
        return 4.0*PI*concentration*r12*dimensionTerm;
    }


    /*
    // Number of SIA in SIA cluster for it to travel in 1D only (no rotations)
    int numSIAfor1D = 3;

    if (hostObject->getAttri(0) >= numSIAfor1D && mobileObject->getAttri(0) >= numSIAfor1D)
        return 0;   // Assume 1D-1D collision negligibly happens

    // 1D + immobile object (rate formula from Sicong He 2025)
    if (mobileObject->getAttri(0) >= numSIAfor1D && mobileObject->getAttri(2) == 0)
        return 8*PI*mobileObject->getDiff()*pow(r12, 2.0)*mobileObject->getNumber(count)/volume*pow(hostObject->getNumber(count)/volume, 4.0/3.0)*volume;

    if (hostObject->getAttri(0) >= numSIAfor1D && hostObject->getAttri(2) == 0)
    {
        // cout << hostObject->getKey() << " " << mobileObject->getKey() << endl;
        // cout << 8*PI*hostObject->getDiff()*pow(r12, 2.0)*hostObject->getNumber(count)/volume*pow(mobileObject->getNumber(count)/volume, 4.0/3.0)*volume << endl;
        // cout << endl;
        return 8*PI*hostObject->getDiff()*pow(r12, 2.0)*hostObject->getNumber(count)/volume*pow(mobileObject->getNumber(count)/volume, 4.0/3.0)*volume;
    }

    // Otherwise it's a 3D+3D reaction
    return 4.0*PI*concentration*r12*dimensionTerm;
    */
}

long double OneLine::computeCombReaction(
                                    const Object* const hostObject,
                                    const Object* const mobileObject,
                                    unordered_map<int64, Object*>& allObjects,
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
            OneLine* prodLine = prodObj->lines[count];
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

    if (!SAV_ON 
        || count == SURFACE_INDEX 
        || count == SUBSURFACE_INDEX
        || (count == BACK_SUBSURFACE_INDEX && BACK_DESORB)
        || (count == BACK_SURFACE_INDEX && BACK_DESORB))
    {
        return;
    }

    // If we have a mV-nH object, or nH object
    if (hostObject->getAttri(0) <= 0 && hostObject->getAttri(2) > 0 && hostObject->getNumber(count) > 0)
    {
        int numHPerCluster = hostObject->getAttri(2);
        int numVacancies = abs(hostObject->getAttri(0));
        double clusterThresholdH;
        if (numVacancies == 0)
        {
            clusterThresholdH = 0;
        }
        else if (numVacancies <= 7)
        {
            int savHThres[8] = {0, 9, 14, 17, 22, 29, 34, 36}; // index = #vac, value = numH that will trigger sav
            clusterThresholdH = savHThres[numVacancies];
        }
        else
        {
            clusterThresholdH = 4.75*numVacancies + 4; // Qianran Yu 2020, did linear fit from graph of excess sav energies
        }
        // Allow overpressured HV clusters and nH clusters to be SAV candidates
        if (numHPerCluster >= clusterThresholdH && numVacancies > 0)
        {
            SAVR = NU0 * exp(-SAV_ENERGY/KB/TEMPERATURE) * hostObject->getNumber(count);
        }
        else if (numHPerCluster >= 1 && numVacancies == 0)
        {
            // double coeff = 0;
            // if (TEMPERATURE < 383)   // custom fitting based on experiments at 383K (simmonds 2017) and 823K (nobuta 2022)
            //     coeff = 0.007;
            // else if (TEMPERATURE > 823)
            //     coeff = 0.0015;
            // else
            //     coeff = 0.007 + (TEMPERATURE-383.0)/(823.0-383.0) * (0.0015-0.007);   // linear interpolation
            
            // SAVR = coeff * hostObject->getNumber(count);
            double hostNumber = hostObject->getNumber(count);
            if (UNIFORM_FREE_H_ON)
                hostNumber = UNIFORM_H_CONCENTRATION*volumeAtIndex(count);
            // SAVR = 1.88 * exp(-0.097/(KB*TEMPERATURE)) * hostNumber;
            SAVR = 0.087*hostNumber;
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
    int64 HKey = 1;

    // only H can recombine at surface and leave surface
    if (!RECOMB_ON 
        || (count != SURFACE_INDEX && (count != BACK_SURFACE_INDEX || !BACK_DESORB)) 
        || hostObject->getKey() != HKey
        || allObjects.find(HKey) == allObjects.end()
        || UNIFORM_FREE_H_ON)
    {
        return;
    }

    double surfaceConc = 0.0;
    int numH = allObjects[HKey]->getNumber(count);  
    surfaceConc = numH / DIVIDING_AREA;  // [cm^-2] concentration

    double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0);  // [110] surface
    double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;

    // Calculate ER recomb rate, only on plasma facing surface
    if (numH >= 1 && HYDROGEN_ON && count == SURFACE_INDEX)   // incident atom collides with adsorbed atom, so hydrogen must be on for this to work
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
            // desorbE = 2.0*(0.525 + 0.591*(1.0/(1.0+exp( (surfaceSaturationFraction-0.247)/0.0692 )))); // from Hodille 2020
            // desorbE = -0.00213989*exp(5.78271*surfaceSaturationFraction) + 1.4965;                  // Ajmalghan 2019 NEB
            // desorbE = -0.00195416 * exp(5.87242*surfaceSaturationFraction) + 1.48996;            // Ajmalghan 2019
            // desorbE = 1.40259 - 0.00881176*exp(5.45029*surfaceSaturationFraction - 1.22515);
            // desorbE = 0.019+1.453/(1.0+exp((surfaceSaturationFraction-1.000)/0.111));
            // desorbE = 2.0*(0.9 - 0.2*surfaceSaturationFraction - 0.7*pow(surfaceSaturationFraction, 12));
            // desorbE = 1.023 + 0.584/(1.0 + exp(7.38e-16 * surfaceConc - 2.85));
            // desorbE = 1.029 + 0.700/(1.0+exp((surfaceSaturationFraction-0.475)/0.151));
            desorbE = 0.8 + 1.4/(1.0+exp((surfaceSaturationFraction-0.3)/0.2));
            // desorbE = -0.58707 * surfaceSaturationFraction + 1.54351;
            // desorbE = 0; // Assume H immediately desorbs at the surface
        double desorptionR = NU0 / maxSurfaceConc; // [cm^2 s^-1]
        recombRLH = desorptionR * exp(-desorbE / (KB * TEMPERATURE)) * surfaceConc * surfaceConc * DIVIDING_AREA;
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
    double volume = volumeAtIndex(count);
    int dimsum = 6;
    double alpha_a = -log(PI*PI*pow(r12, 3.0) / volume / hostN);
    double alpha_b = -log(PI*PI*pow(r12, 3.0) / volume / mobileN);
    switch (dimsum) {
        case 6: // 3D + 3D
            term = hostDiff + mobileDiff;
            break;
        case 4: // 3D + 1D
            if (hostDim == 1 && mobileDim == 3)
                term = hostObject->getDiff()*(mobileN / volume)*(2.0*PI*pow(r12, 3.0)) + mobileDiff;
            else
                term = mobileDiff*(mobileN / volume)*(2.0*PI*pow(r12, 3.0)) + hostDiff;
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
    if (rate > 0)
        secondR[mobileObjectKey] = rate;
    else
        secondR.erase(mobileObjectKey);
}
