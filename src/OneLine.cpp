// OneLine.cpp -- implementations of class OneLine
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include "OneLine.h"
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

OneLine::OneLine() : diffRToF(0.0), diffRToB(0.0), sinkR(0.0), SAVR(0.0), recombRER(0.0), recombRLH(0.0), totalRate(0.0)
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
    if (sinkR >= tempRate) {
        return SINK;
    }
    else {
        tempRate -= sinkR;
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
                          const int count)
{
    double rate = computeCombReaction(hostObject, newObject, count);
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
                             const int n)
{
    double rate = computeCombReaction(hostObject, mobileObject, n);
    secondR[mobileObject->getKey()] = rate;
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
    totalRate += sinkR;    /* add sink rate */
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
    fs << "(diff)" << diffRToF << ", " << diffRToB << "    " << "(sink)" << sinkR << "    ";
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
        dissociationR[index] = computeDissReaction(hostObject, index, count);
    }
    unordered_map<int64, Object*>::iterator iter;
    for (iter = mobileObjects.begin(); iter != mobileObjects.end(); ++iter) {
        double rate = computeCombReaction(hostObject, iter->second, count);
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
    double frontConcentration = objectN[1] / VOLUME;
    double backConcentration = objectN[2] / VOLUME;
    if (count == SURFACE_INDEX)
    {
        concentration = objectN[0] / SURFACE_VOLUME;  /* First mesh element has slightly larger volume of normal mesh element size + surface layer size */
        frontConcentration = H_SATURATION_CONCENTRATION;
    }
    else if (count == 1)
    {
        concentration = objectN[0] / VOLUME;
        frontConcentration = objectN[1] / SURFACE_VOLUME;
    }
    else if (count == HUGE_INDEX - 1)
    {
        concentration = objectN[0] / VOLUME;
        backConcentration = objectN[2] / HUGE_VOLUME;
    }
    else if (count == HUGE_INDEX)
    {  
        concentration = objectN[0] / HUGE_VOLUME;
        backConcentration = objectN[2] / SURFACE_VOLUME;
    }
    else if (count == BACKSURFACE_INDEX)
    {
        concentration = objectN[0] / SURFACE_VOLUME;
        frontConcentration = objectN[1] / HUGE_VOLUME;
        backConcentration = 0.0;
    }
    else
    {
        concentration = objectN[0] / VOLUME;
    }

    // Account for special cases from 2020 Zhenhou Wang for surface and back surface
    if((hostObject -> getKey() == 1)
        && (count == SURFACE_INDEX || count == 1 || count == HUGE_INDEX || count == BACKSURFACE_INDEX))
    {
        double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0);
        double surfaceConc = 0.0; 
        int64 HKey = 1;
        if (allObjects.find(HKey) != allObjects.end())
            surfaceConc = allObjects[HKey]->getNumber(0) / DIVIDING_AREA;  // [cm^-2] concentration
        double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;

        // special case for 1H diffusion from Bulk to Surface 
        if (count == 1)
        {
            double jumpingDist = maxSurfaceConc / 6.0 / DENSITY;
            double freq = NU0 * exp(-H_MIGRATION_ENERGY / KB / TEMPERATURE);
        
            prefactor = freq * jumpingDist * (1 - surfaceSaturationFraction) * DIVIDING_AREA;
            diffRToF = prefactor * concentration;
            diffRToB = 0.0;

            // Normal Diffusion to the second bulk element
            if (concentration > backConcentration)
            {
                double lengthb = ELEMENT_THICKNESS * NM_TO_CM; // first element to second element distance (20nm) 
                /* if diffusable */
                prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthb;
                diffRToB = prefactor*(concentration - backConcentration);
            }
            return;
        }
        // special case for 1H diffusion from Surface to Bulk
        else if (count == SURFACE_INDEX)
        {
            double absorbE = 1.10 + 0.939*( 1.0 / ( 1.0 + exp( (surfaceSaturationFraction - 0.232)/0.0683 ) ) );  // absorption energy barrier into bulk from Hodille 2020
            double freq = NU0 * exp(-absorbE / KB / TEMPERATURE);
            diffRToB = freq * surfaceConc * DIVIDING_AREA;
            diffRToF = 0.0;
            return;
        }
        else if (count == HUGE_INDEX)
        {
            double jumpingDist = maxSurfaceConc / 6.0 / DENSITY;
            double freq = NU0 * exp(-H_MIGRATION_ENERGY / KB / TEMPERATURE);
        
            prefactor = freq * jumpingDist * (1 - surfaceSaturationFraction) * DIVIDING_AREA;
            diffRToB = prefactor * concentration;
            diffRToF = 0.0;

            // Normal Diffusion to the regular bulk element before the Huge element
            if (concentration > backConcentration)
            {
                double lengthf = (ELEMENT_THICKNESS / 2.0 + HUGE_THICKNESS / 2.0) * NM_TO_CM; // first element to second element distance (20nm) 
                /* if diffusable */
                prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthf;
                diffRToF = prefactor*(concentration - backConcentration);
            }
            return;
        }
        else if (count == BACKSURFACE_INDEX)
        {
            double absorbE = 1.10 + 0.939*( 1.0 / ( 1.0 + exp( (surfaceSaturationFraction - 0.232)/0.0683 ) ) );  // absorption energy barrier into bulk from Hodille 2020
            double freq = NU0 * exp(-absorbE / KB / TEMPERATURE);
            diffRToF = freq * surfaceConc * DIVIDING_AREA;
            diffRToB = 0.0;
            return;
        }
    }
    else // avoid unnecessary calculation if we are processing special cases
    {
        /* by having two diffusion rates, this rate will never be less than 0 */
        /* length measured in cm */
        double lengthf = 0.0, lengthb = 0.0;
        if(count == SURFACE_INDEX){
            lengthf = (ELEMENT_THICKNESS + SURFACE_THICKNESS) / 2. * NM_TO_CM; // thickness of W surface is 0.544nm, this length is centroid to vacuum 
            lengthb =  (ELEMENT_THICKNESS + SURFACE_THICKNESS / 2.) * NM_TO_CM; /* surface to first element distance */
        }else if(count == 1){
            lengthf =  (ELEMENT_THICKNESS + SURFACE_THICKNESS / 2.) * NM_TO_CM; /* surface to first element distance */
            lengthb = ELEMENT_THICKNESS * NM_TO_CM; // first element to second element distance (20nm) 
        }else{
            lengthf = lengthb = ELEMENT_THICKNESS * NM_TO_CM; /* other element distances */
        }

        /* 
        * 1. compute diffusion rate to the front element 
        * Diffusion goes from area of higher concentration to lower concentration
        * Object not allowed to diffuse out through the front
        */
        if (concentration > frontConcentration && count != 0) {
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
        if (concentration > backConcentration) {
            /* if diffusable */
            prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthb;
            diffRToB = prefactor*(concentration - backConcentration);
        }
        else {
            diffRToB = 0.0;
        }
    }
}

void OneLine::computeSinkReaction(const Object* const hostObject, const int count)
{
    if (!SINK_ON || count == SURFACE_INDEX || count == HUGE_INDEX || count == BACKSURFACE_INDEX)
    {
        sinkR = 0.0;
        return;
    }

    sinkR = hostObject->getNumber(count)*hostObject->getDiff()*hostObject->getSink();
}

long double OneLine::computeDissReaction(
                                  const Object* const hostObject,
                                  const int index,
                                  const int count) const
{
    if (!DISS_ON || count == SURFACE_INDEX || count == HUGE_INDEX || count == BACKSURFACE_INDEX)
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
        int attr[LEVELS] = { 0 };
        attr[index] = hostObject->signof(hostObject->getAttri(index));
        Object tempObject(attr, count);
        return 4.0 * PI * hostObject->getR1e() / avol * tempObject.getDiff() * hostObject->getBind(index) * hostObject->getNumber(count);
    }
    return 0.0;
}

long double OneLine::computeCombReaction(
                                    const Object* const hostObject,
                                    const Object* const mobileObject,
                                    const int count) const
{    
    if (!COMB_ON || count == SURFACE_INDEX || count == HUGE_INDEX || count == BACKSURFACE_INDEX)
    {
        return 0.0;
    }

    double concentration;
    double r12;
    double dimensionTerm;
    double volume;
    if(count == SURFACE_INDEX || count == BACKSURFACE_INDEX){
        volume = SURFACE_VOLUME;
        // volume on surface layer = volume/20nm * 0.54nm
        
    }
    else if (count == HUGE_INDEX)
    {
        volume = HUGE_VOLUME;
    }
    else{
        volume = VOLUME;
    }
    
    if (hostObject->getKey() != mobileObject->getKey()) {
        concentration = hostObject->getNumber(count)*mobileObject->getNumber(count) / volume;
        
    }else {
        concentration = hostObject->getNumber(count)*(hostObject->getNumber(count) - 1) / volume;
        
    }

    // H+H-->2H
    // Disable H clustering for now b/c not sure how it works with SAV
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
    return 4.0*PI*concentration*r12*dimensionTerm;
}

void OneLine::computeSAVReaction(
                                 const Object* const hostObject,
                                 const int count)
{
    SAVR = 0;

    // TODO: no SAV for huge as well
    if (!SAV_ON || count == SURFACE_INDEX || count == HUGE_INDEX || count == BACKSURFACE_INDEX)
    {
        return;
    }

    // If we have a mV-nH object, or nH object
    if (hostObject->getAttri(0) <= 0 && hostObject->getAttri(2) > 0 && hostObject->getNumber(count) > 0)
    {
        int numH = hostObject->getAttri(2);
        int numVacancies = abs(hostObject->getAttri(0));
        // double thresholdH = 4.75 + 4*numVacancies;
        double thresholdH = 4*numVacancies;
        if (numH > thresholdH)
        {
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
    if (!RECOMB_ON || count != SURFACE_INDEX || count != BACKSURFACE_INDEX|| hostObject->getKey() != 1)
    {
        return;
    }

    double maxSurfaceConc = 6.9 * pow(DENSITY, 2.0/3.0);  // Hodille 2020
    double surfaceConc = 0.0;
    int64 HKey = 1;
    int numH = allObjects[HKey]->getNumber(0);
    if (allObjects.find(HKey) != allObjects.end())
        surfaceConc = numH / DIVIDING_AREA;  // [cm^-2] concentration

    // Calculate ER recomb rate
    if (count == SURFACE_INDEX && numH >= 1)
    {
        double crossSectionERRecomb = 1.7e-17; // [cm^2] cross-section of ER recombination from Zhenhou Wang 2020
        recombRER = FLUX_H * crossSectionERRecomb * surfaceConc * DIVIDING_AREA; 
    }
    else
        recombRER = 0.0;
    
    // Calculate LH recomb rate
    if (numH >= 2)
    {
        double surfaceSaturationFraction = surfaceConc / maxSurfaceConc;
        double desorptionR = NU0 * pow(ALATT, 2.0); // [cm^2 s^-1]
        double desorbE = 0.525 + 0.591*( 1.0 / ( 1.0 + exp( (surfaceSaturationFraction - 0.247)/0.0692 ) ) );  // Energy barrier to desorb from surface (Hodille 2020)
        recombRLH = desorptionR * exp(-2.0 * desorbE / KB / TEMPERATURE) * surfaceConc * surfaceConc * DIVIDING_AREA;
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
    double volume;
    if(count == SURFACE_INDEX || count == BACKSURFACE_INDEX)
    {
        volume = SURFACE_VOLUME; 
    }
    else if (count == HUGE_INDEX)
    {
        volume = HUGE_VOLUME;
    }
    else{
        volume = VOLUME;
    }

    double term = 0.0;
    double hostDiff = hostObject->getDiff(), mobileDiff = mobileObject->getDiff();
    int hostDim = hostObject->getDim(), mobileDim = mobileObject->getDim();
    int hostN = hostObject->getNumber(count), mobileN = mobileObject->getNumber(count);
    // int dimsum = hostDim + mobileDim; 
    int dimsum = 6;
    double alpha_a = -log(PI*PI*pow(r12, 3.0) / volume / hostN);
    double alpha_b = -log(PI*PI*pow(r12, 3.0) / volume / mobileN);
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