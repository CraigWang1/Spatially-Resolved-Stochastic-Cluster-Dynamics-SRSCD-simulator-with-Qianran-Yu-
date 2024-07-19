// OneLine.h -- class OneLine
#pragma once
#include<iostream>
#include<unordered_map>
#include"Object.h"

class OneLine {
private:
    // private data member
    long double diffRToF;            // diffusion rate from jth element to (j-1)th element
    long double diffRToB;            // diffusion rate from jth element to (j+1)th latter one
    long double sinkR;                 // go to sink reaction rate
    long double dissociationR[LEVELS];
    long double SAVR;                // super abundant vacancy reaction rate
    long double recombRToF;          // surface recombination rate at front surface of material (i.e. for H+H = H2 and leaves material at surface)
    long double recombRToB;          // surface recombination rate at back surface of material
    // dissociation reaction rate
    long double totalRate;             // totalRate of the line;
    std::unordered_map<int64, long double> secondR;  // second reaction rate
    
    // private functions
    void setOneLine(const Object* const, const int, unordered_map<int64, Object*>&); // build one line
    void computeDiffReaction(const Object* const, const int);  // compute diffusion rate
    void computeSinkReaction(const Object* const, const int);         // compute absorption reaction rate
    void computeSAVReaction(const Object* const, const int);   // compute super-abundant-vacancy rate
    void computeRecombReaction(const Object* const, const int);
    double computeDimensionTerm(const double, const Object* const, const Object* const, const int) const;
    
public:
    OneLine(const Object* const, const int, unordered_map<int64, Object*>&);
    OneLine();  /* this constructor is used when the user only wants to calculate up to a couple rate and not all the rates for a temporary purpose */
    Reaction selectReaction(const Object* const,int64&, long double&);
    void addReaction(const Object* const, const Object* const, const int); /* add one reaction when a new object is been created */
    void removeReaction(const int64); /* delete one reaction when an object is deleted */
    void updateReaction(const Object* const, const Object* const, const int);/* when the number of another object has changed, one rate in this line should be changed */
    void updateLine(const Object* const, const int, unordered_map<int64, Object*>&); /* when number of this object has changed, rates in this line should be updated */
    void updateDiff(const Object* const, const int);
    const long double computeTotalRate();
    long double computeDissReaction(const Object* const, const int, const int) const; // compute dissociation reaction rate
    long double computeCombReaction(const Object* const, const Object* const, const int) const;  // compute 2nd order reaction rate
    void setDissReaction(const int, long double);    // set this rate from an external source if they wish to correct it
    void setCombReaction(const int64, long double);  // set this rate from an external source if they wish to correct it (i.e. net comb/diss rate)
    void display(const Object* const);
};
