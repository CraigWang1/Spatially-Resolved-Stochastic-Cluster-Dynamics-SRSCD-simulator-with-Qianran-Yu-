// OneLine.h -- class OneLine
#pragma once
#include<iostream>
#include<unordered_map>
#include"Object.h"

class OneLine {
private:
    // private data member
    int sinks[LEVELS+1][POINTS];
    long double diffRToF;            // diffusion rate from jth element to (j-1)th element
    long double diffRToB;            // diffusion rate from jth element to (j+1)th latter one
    long double sinkR;                 // go to sink reaction rate
    long double dissociationR[LEVELS];
    long double SAVR;                // super abundant vacancy reaction rate

    // Prefactors to make calculating rates easier.
    // Usually they are the rates associated with just 1 instance of the hostObject
    long double diffRPrefactor;
    long double sinkRPrefactor;
    long double dissociationRPrefactor[LEVELS];
    long double SAVRPrefactor;
    long double combRPrefactor;

    // dissociation reaction rate
    long double totalRate;             // totalRate of the line;
    std::unordered_map<int64, long double> secondR;  // second reaction rate
    std::unordered_map<int64, long double> secondRPrefactor; // prefactors for comb rate
    
    // private functions
    void computePrefactors(const Object* const, const int, unordered_map<int64, Object*>&);
    void computeSinkRPrefactor(const Object* const);
    void computeDissRPrefactor(const Object* const, const int, const int);
    void computeSAVRPrefactor(const Object* const, const int);
    void computeCombRPrefactor(const Object* const, const Object* const, const int);

    void setOneLine(const Object* const, const int, unordered_map<int64, Object*>&); // build one line
    void computeDiffReaction(const Object* const, const int);  // compute diffusion rate
    void computeSinkReaction(const Object* const, const int);         // compute absorption reaction rate
    void computeDissReaction(const Object* const, const int, const int); // compute dissociation reaction rate
    void computeSAVReaction(const Object* const, const int);   // compute super-abundant-vacancy rate
    double computeCombReaction(const Object* const, const Object* const, const int);  // compute 2nd order reaction rate
    double computeDimensionTerm(const double, const Object* const, const Object* const, const int);
    
public:
    OneLine(const Object* const, const int, unordered_map<int64, Object*>&);
    Reaction selectReaction(const Object* const,int64&, long double&);
    void addReaction(const Object* const, const Object* const, const int); /* add one reaction when a new object is been created */
    void removeReaction(const int64); /* delete one reaction when an object is deleted */
    void updateReaction(const Object* const, const Object* const, const int);/* when the number of another object has changed, one rate in this line should be changed */
    void updateLine(const Object* const, const int, unordered_map<int64, Object*>&); /* when number of this object has changed, rates in this line should be updated */
    const long double computeTotalRate();
    const long double getDiffRateF() const;
    const long double getDiffRateB() const;
    void display(const Object* const);
};
