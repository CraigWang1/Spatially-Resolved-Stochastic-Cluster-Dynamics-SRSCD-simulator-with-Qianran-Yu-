// SCDWrapper -- SCDWrapper that manipulate the whole bulk
#pragma once
#include"Damage.h"
#include"Bundle.h"
#include"cpdf.h"
#include"rvgs.h"
#include"CascadeDamage.h"
#include"BoundaryChange.h"
#include"constants.h"
// #include"gnuplot_i.h"
#include <string>
#include <iomanip>
#include <cassert>
#include <random>

class SCDWrapper {
private:
    /* private data member */
    unordered_map<int64, Object*> allObjects;  // map that store all object
    unordered_map<int64, Object*> mobileObjects;  // map that store mobile object
    unordered_map<int64, Object*> HObjects;    // map that stores nH objects
    unordered_map<Object*, Bundle*> linePool;
    unordered_map<int64, int> surface;
    unordered_map<int64, int> bottom;
    unordered_map<int, double> formationE; 
    unordered_map<int, double> HSaturationLimit; // look up table for HSaturationLimit
    vector<BoundaryChange> leftBoundaryChangeQ, rightBoundaryChangeQ;

    Damage damage;
    Cpdf cpdf;
    int sinks[LEVELS+1][POINTS];
    long double sinkDissRate[2][POINTS];
    // dissociation rate of V/H from dislocations
    int reactions[8][POINTS];
    int startIndex, endIndex; // the indices of which points this processor is responsible for
    
    /* hold reactions, 1st dimension is reaction type, second dimension is element, value is total number of this reaction */
    long double matrixRate[POINTS];    // total rate in every element(point)
    long double bulkRate;  // total rate in the whole bulk;
    long double domainRate; // total rate in the volume elements that this processor is responsible for
    long double noneRate;  // the total rate of null event occuring in the domain that this processor is responsible for
    double totalDpa;
    bool lastElemSaturated;
    enum InsertStyle {INTERSTITIAL, SUBSTITUTIONAL};
    ofstream fs1, fs2, fs3, fs4, fs5, fs6;
    fstream fs;
    fstream fv;
    // GnuplotS gs, gr; /* plot species.out and reaction */
    // GnuplotS gd1, gd2; /* damage graph 1 and damage graph 2*/
    // GnuplotS gh1, gh2, gh3; /* H deposition graph 1,2,3 */
    // GnuplotS gv;
    
    /* private functions */
    /* set sinks function */
    void setSinks();
    /* function to compute dissociation reaction rate */
    void computeSinkDissRate(const int, const int);
    /* a way to get key from attribute */
    int64 attrToKey(const int* const); /* same with Object::setKey() */
    /* decide insertion mode function */
    int64 atomProperty(InsertStyle, const int);
    /* moniter map functions */
    void updateObjectInMap(Object *, const int);
    /* Updates the object's reaction rates in that nth mesh element */
    void addNewObjectToMap(Object*);
    /* add to map by object pointer */
    void addToObjectMap(const int64, const int, const int number=1);
    /* general-use function to add a number of instances to an object in the nth mesh element */
    void reduceFromObjectMap(const int64, const int, const int number=1);
    /* general-use function to remove one instance of an object in the nth mesh element */
    void removeObjectFromMap(const int64); 
    /* remove one object from map */
    void addReactionToOther(const Object* const);
    /* Impact of one new mobile object to other existing objects */
    void updateRateToOther(const Object* const, const int);
    /* when number of this object changes, rates related to this object change also */
    void removeRateToOther(const int64);
    /* when one mobile object has been removed,
     ** rates of this object with other extisting objects should also be removed
     */
    void updateCombDissRatePair(int combinedObjectAttr[], const int);
    void updateNetCombDissRate(const Object* const hostObject, const int);
    void updateSinks(const int, const int*); /* only for restart use */
    /* process event functions */
    void processDiffEvent(Object*, const int, const char);     /* process diffusion reactionObject */
    void processSinkEvent(Object*, const int);     /* process absorption reaction */
    void processDissoEvent(Object*, const int, const int64, fstream& ); /* process dissociation event */
    void processCombEvent(Object*, const int, const int64, fstream& );  /* process combination reaction */
    void processSAVEvent(Object*, const int);      /* process super-abundant-vacancy reaction */
    void processRecombEvent(Object*, const int, bool);   /* process surface recombination event: 1H+1H forms H2 and leaves material surface */
    void processSinkDissEvent(const int, const int); /* process dissociation from sink event */
    /* get insertion functions */
    void getElectronInsertion(const int);
    void getNeutronInsertion(const int);
    void getIonInsertion(const int, const double, fstream&);
    void getParticleInsertion(const int, const double, fstream&);    // deal with damage[0]
    void getHeInsertion(const int);  // deal with damage[1]
    void getHInsertion(const int, const double, fstream&);   // deal with damage[2]
    /* write file funcitons*/
    void writeSinkFile(const Object* const, const long int n, const double);
    /* sinks.out only writes when sink reaction happens, now this function is not "writing things" but only updating sinks[][] */
    void writeSpeciesFile(const double, const long int, const int);
    void writeClusterFile(const double, const long int);
    /* distinguish reaction type function */
    int countDefectNumber(const int, string);
    /* this function counts the number of defects(V,SIA,H...) in every element, returns total number of this defect in the bulk */
    void countRatioDistribution(double&);
    /* count ratio of H to V in the bulk and at the same time decide whether to supplement H to system */
    void sizeDistribution(); /* get size distribution */
    void writeReaction(); /* take down reactions for drawing */
    double getHSaturationConcentration() const;

public:  
    SCDWrapper();  // constructor: for start ;
    //SCDWrapper();  // constructor: for restart;
    void computeMatrixRate(const int n);  // computes total rate in element n
    void updateMatrixRate(const int n, const Reaction reaction); // computes total rate in elements that are affected by reaction
    void computeBulkRate();
    void computeDomainRate();
    long double getBulkRate();
    long double getDomainRate();
    Object* selectReaction(int64&, Reaction&, int&);  // select reaction
    Object* selectDomainReaction(int64&, Reaction&, int&);
    void processEvent(const Reaction, Object*, const int, const int64, const double, const double);    // deal with reactions
    ~SCDWrapper();          /* destructor to delete everything newed */
    // get series functions that allow direct manipulation on private data member
    unordered_map<int64, Object*>* getAllObjects();
    unordered_map<int64, Object*>* getMobileObjects();
    unordered_map<Object*, Bundle*>* getLinePool();
    void examineRate(); /* computes matrix rate in all points*/
    void examineDomainRate(); /* computes matrix rate in points that this processor is responsible for */
    /* output file functions */
    void displayDamage();
    void displayAllObject();
    void writeFile(const double, const long int, const int);
    /* testing functions which includes display/output file functions that might be used when debugging */
    void writeSinkFile(const double, const long int, const int);
    void display() const;   // display rateMatrix within in one element
    friend void restart(long int&, double&, SCDWrapper*);
    /* test functions, don't need in the future */
    void test(const int); /* count how many vacancies now in the bulk */
    
    /* draw picture functions */
    void drawSpeciesAndReactions(double&t);
    void drawDamage(double&); /* this function draw all diagrams in damage process */
    void drawHD(double&); /* this function draw all diagrams in H deposition process */
    void writeVacancy();
    double getTotalDpa();

    /* Additional parallel processing functions */
    void setDomain(int, int);
    void fillNoneReaction(const double&);
    void clearNoneReaction();
    vector<BoundaryChange>* getLeftBoundaryChangeQ();
    vector<BoundaryChange>* getRightBoundaryChangeQ();
    void clearBoundaryChangeQs();
    void implementBoundaryChanges(vector<BoundaryChange>&);
};
