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
                 unordered_map<int64, Object*>& mobileObjects) :totalRate(0.0)
{
    setOneLine(hostObject, count, mobileObjects);
    
}

Reaction OneLine::selectReaction(
                                 const Object* const hostObject,
                                 int64& theOtherKey,
                                 double& randRate)
{
    int index = 0;
    double tempRate = randRate;
    std::unordered_map<int64, long double>::iterator iter = secondR.begin();
    if (totalRate < tempRate) {
        randRate -= totalRate;
        return NONE;
    } // the reaction is not positioned in this line
    if(diffRToF > tempRate){
        return DIFFUSETOF;
    }
    else {
        tempRate -= diffRToF;
    }
    if (diffRToB > tempRate) {
        return DIFFUSETOB;
    }
    else {
        tempRate -= diffRToB;
    }
    if (sinkR > tempRate) {
        return SINK;
    }
    else {
        tempRate -= sinkR;
    }
    while (index < LEVELS) {
        if (dissociationR[index] > tempRate) {
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
        if (iter->second > tempRate) {
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
                         unordered_map<int64, Object*>& mobileObjects)
{
    secondR.clear();
    setOneLine(hostObject, count, mobileObjects);
}

const long double OneLine::computeTotalRate()
{
    int i;
    unordered_map<int64,long double>::iterator iter;
    totalRate = 0.0;
    totalRate += diffRToF; /* add one diffusion rate */
    totalRate += diffRToB; /* add another diffusion rate*/
    totalRate += sinkR; /* add sink rate */
    for (i = 0; i < LEVELS; i++) {
        totalRate += dissociationR[i];
    }
    for (iter = secondR.begin(); iter != secondR.end(); ++iter) {
        totalRate += iter->second;
    }
    return totalRate;
}

const double OneLine::getDiffRateF() const
{
    return diffRToF;
}

const double OneLine::getDiffRateB() const
{
    return diffRToB;
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
    fs << endl;
    fs.close();
}

/* private function implementations */
void OneLine::setOneLine(
                         const Object* const hostObject,
                         const int count,
                         unordered_map<int64, Object*>& mobileObjects)
{
    computeDiffReaction(hostObject, count);
    computeSinkReaction(hostObject, count);
    for (int index = 0; index < LEVELS; index++) {
        computeDissReaction(hostObject, index, count);
    }
    unordered_map<int64, Object*>::iterator iter;
    for (iter = mobileObjects.begin(); iter != mobileObjects.end(); ++iter) {
        double rate = computeCombReaction(hostObject, iter->second, count);
        std::pair<int64, double> oneReaction(iter->first, rate);
        secondR.insert(oneReaction);
    }
    computeTotalRate();
}

void OneLine::computeDiffReaction(const Object* const hostObject, const int count)
{
    if (!DIFF_ON)
    {
        diffRToF = 0.0;
        diffRToB = 0.0;
        return;
    }

    /* by having two diffusion rates, this rate will never be less than 0 */
    /* length measured in cm */
	double lengthf = 0.0, lengthb = 0.0;
	if(count == 0){
		// lengthf = 2.74e-8; // thickness of W surface is 0.544nm, this length is centroid to vacuum 
		lengthf = 1.0e-6 + SURFACE_THICKNESS / 2. * 1e-7;
        // lengthf = 2.0e-6 + SURFACE_THICKNESS / 2. * 1e-7;    /* when H is oversaturated, act as if there is another spatial element to the left that it can diffuse out into */
        // lengthf = SURFACE_THICKNESS / 2. * 1e-7;
        lengthb = 2.0e-6 + SURFACE_THICKNESS / 2. * 1e-7; /* surface to first element distance */
    }else if(count == 1){
		lengthf = 2.0e-6 + SURFACE_THICKNESS / 2. * 1e-7;
		lengthb = 2.0e-6; // first element to second element distance (20nm) 
	}else{
		lengthf = lengthb = 2.0e-6; /* other element distances */
	}

    double prefactor = 0.0;
    int objectN[3];   
    hostObject->getThreeNumber(count, objectN);
    double concentration = 0;
    double frontConcentration = 0;
    double backConcentration = objectN[2] / VOLUME;
    if (count == 0)
    {
        concentration = objectN[0] / FRONTMOST_VOLUME;
        frontConcentration = H_SATURATION_CONCENTRATION;
    }
    else if (count == 1)
    {
        concentration = objectN[0] / VOLUME;
        frontConcentration = objectN[1] / FRONTMOST_VOLUME;
    }
    else
    {
        concentration = objectN[0] / VOLUME;
        frontConcentration = objectN[1] / VOLUME;
    }

    /* 
     * 1. compute diffusion rate to the front element 
     * Diffusion goes from area of higher concentration to lower concentration
     */
    if (concentration > frontConcentration) {
        /* if diffusable, surface objects diffusing into vacuum is considered */
        prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthf;
        diffRToF = prefactor*(concentration - frontConcentration);
        //diffRToF = 0.0;
        
        int64 HKey = 1;
        if(count == 0 && hostObject->getKey() != HKey){
            diffRToF = 0.0;
        } // only allow H to diffuse out when it's oversaturated
    }
    else {
        diffRToF = 0.0;
    }
    /* 2. compute diffusion rate to the back element*/
    if (concentration > backConcentration) {
        /* if diffusable */
        prefactor = hostObject->getDiff() * DIVIDING_AREA / lengthb;
        diffRToB = prefactor*(concentration - backConcentration);
    }
    else {
        diffRToB = 0.0;
    }


    if (count == POINTS - 1) {
        //objects at bottom is not allowed to diffuse into vacuum
        diffRToB = 0.0;
    }
}

void OneLine::computeSinkReaction(const Object* const hostObject, const int count)
{
    if (!SINK_ON)
    {
        sinkR = 0.0;
        return;
    }

    sinkR = hostObject->getNumber(count)*hostObject->getDiff()*hostObject->getSink();
}

void OneLine::computeDissReaction(
                                  const Object* const hostObject,
                                  const int index,
                                  const int count)
{
    if (!DISS_ON)
    {
        dissociationR[index] = 0.0;
        return;
    }

    if (hostObject->getAttri(index) != 0) {
        int attr[LEVELS] = { 0 };
        attr[index] = hostObject->signof(hostObject->getAttri(index));
        Object tempObject(attr, count);
        if(count == 0 && hostObject->getKey() == 2){
            dissociationR[index] = 4.0 * PI * hostObject->getR1e() / avol * tempObject.getDiff() * hostObject->getBindSH() * hostObject->getNumber(count);
            
        }else{
            dissociationR[index] = 4.0 * PI * hostObject->getR1e() / avol * tempObject.getDiff() * hostObject->getBind(index) * hostObject->getNumber(count);
            
        }
    }
    else {
        dissociationR[index] = 0.0;
        
    }
}

double OneLine::computeCombReaction(
                                    const Object* const hostObject,
                                    const Object* const mobileObject,
                                    const int count)
{    
    if (!COMB_ON)
    {
        return 0.0;
    }

    double concentration;
    double r12;
    double dimensionTerm;
    double volume;
    if(count == 0){
        volume = FRONTMOST_VOLUME;
        // thin surface element + first bulk element are merged
    }else{
        volume = VOLUME;
    }
    
    if (hostObject->getKey() != mobileObject->getKey()) {
        concentration = hostObject->getNumber(count)*mobileObject->getNumber(count) / volume;
        
    }else {
        concentration = hostObject->getNumber(count)*(hostObject->getNumber(count) - 1) / volume;
    }

    // H+H-->2H
    if(count != 0 && hostObject->getKey() == 1 && mobileObject->getKey() == 1){
        return 0.0;
    }

    // Since H+H->2H can only happen on the thin surface element,
    // temporarily "unmerge" the first bulk and surface element to process this.
    // (First bulk and surface element are merged otherwise for more stable diffusion processing)
    if (count == 0 && hostObject->getKey() == 1 && mobileObject->getKey() == 1){
        // Assume uniform concentration between among surface and first bulk element
        double ratio = SURFACE_THICKNESS / ELEMENT_THICKNESS;
        int numSurfaceH = round(hostObject->getNumber(count) * ratio);
        concentration = numSurfaceH*(numSurfaceH - 1) / SURFACE_VOLUME;
    }
    
    if(hostObject->getKey() == 1 && mobileObject->getKey() == 2){
        return 0.0;
    }
    if(hostObject->getKey() == 2 && mobileObject->getKey() == 1){
        return 0.0;
    }
    r12 = hostObject->getR1() + mobileObject->getR1();
    dimensionTerm = computeDimensionTerm(r12, hostObject, mobileObject, count);
    return 4.0*PI*concentration*r12*dimensionTerm;
}

double OneLine::computeDimensionTerm(
                                     const double r12,
                                     const Object* const hostObject,
                                     const Object* const mobileObject,
                                     const int count)
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

