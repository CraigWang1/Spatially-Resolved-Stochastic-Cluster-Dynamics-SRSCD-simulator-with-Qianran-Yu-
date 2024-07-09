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
#include <algorithm>
#include <list>
#include <omp.h>

struct SectorBound {
    int leftIdx;
    int rightIdx;  // inclusive
};

class SCDWrapper {
private:
    /* private data member */
    unordered_map<int64, Object*> allObjects;  // map that store all object
    unordered_map<int64, Object*> mobileObjects;  // map that store mobile object
    unordered_map<Object*, Bundle*> linePool;
    unordered_map<int64, int> surface;
    unordered_map<int64, int> bottom;
    unordered_map<int, double> formationE; 
    unordered_map<int, double> HSaturationLimit; // look up table for HSaturationLimit

    Damage damage;
    Cpdf cpdf;
    int sinks[LEVELS+1][POINTS];
    long double sinkDissRate[2][POINTS];
    // dissociation rate of V/H from dislocations
    int reactions[8][POINTS];

    /* hold reactions, 1st dimension is reaction type, second dimension is element, value is total number of this reaction */
    long double matrixRate[POINTS];    // total rate in every element(point)
    int matrixNumReactions[POINTS];    // number of nonzero reactions in each spatial element
    long double bulkRate;  // total rate in the whole bulk;
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

    /* Parameters for running in parallel */
    int startIdx, endIdx;
    vector<SectorBound> sectorBounds; // the start and end indices (inclusive) of each sector in the spatial map
    long double sectorRates[3]; // the total reaction rates in each sector
    list<BoundaryChange*> leftBoundaryChanges;  // record changes in mobile species counts at the boundaries and ghost region of the domain this processor is in charge of, so can later communicate to other processors
    list<BoundaryChange*> rightBoundaryChanges; 

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
    void reduceFromObjectMap(const int64, const int);
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
    void updateSinks(const int, const int*); /* only for restart use */
    /* process event functions */
    void processDiffEvent(Object*, const int, const char);     /* process diffusion reactionObject */
    void processSinkEvent(Object*, const int);     /* process absorption reaction */
    void processDissoEvent(Object*, const int, const int64, fstream& ); /* process dissociation event */
    void processCombEvent(Object*, const int, const int64, fstream& );  /* process combination reaction */
    void processSAVEvent(Object*, const int);      /* process super-abundant-vacancy reaction */
    void processSinkDissEvent(const int, const int); /* process dissociation from sink event */
    /* get insertion functions */
    void getElectronInsertion(const int);
    void getNeutronInsertion(const int);
    void getIonInsertion(const int, fstream&);
    void getParticleInsertion(const int, fstream&);    // deal with damage[0]
    void getHeInsertion(const int);  // deal with damage[1]
    void getHInsertion(const int, fstream&);   // deal with damage[2]
    /* write file funcitons*/
    void writeSinkFile(const Object* const, const long int n, const double);
    /* sinks.out only writes when sink reaction happens, now this function is not "writing things" but only updating sinks[][] */
    void writeSpeciesFile(const double, const long int);
    void writeClusterFile(const double, const long int);
    /* distinguish reaction type function */
    int countDefectNumber(const int, string);
    /* this function counts the number of defects(V,SIA,H...) in every element, returns total number of this defect in the bulk */
    void countRatioDistribution(double&);
    /* count ratio of H to V in the bulk and at the same time decide whether to supplement H to system */
    void sizeDistribution(); /* get size distribution */
    void writeReaction(); /* take down reactions for drawing */

    double getHSaturationLimit(int m); /*calculate the saturation limit for H with correction terms based on Max vacancy cluter V number*/
    double getVMonomerBindingE(int m); /*helper funciton for HSatE*/
    double getVFreeFormationE(int m); /*helper funciton for HSatE*/
   
    int getMaxVNum(); /* retrieve the biggest number of V among all VmHn objects*/

public:  
    SCDWrapper();  // constructor: for start ;
    //SCDWrapper();  // constructor: for restart;
    void computeMatrixRate(const int n);  // computes total rate in element n
    void updateMatrixRate(const int n, const Reaction reaction); // computes total rate in elements that are affected by reaction
    void computeBulkRate();
    long double getBulkRate();
    Object* selectReaction(int64&, Reaction&, int&);  // select reaction
    void processEvent(const Reaction, Object*, const int, const int64, const double);    // deal with reactions
    ~SCDWrapper();          /* destructor to delete everything newed */
    // get series functions that allow direct manipulation on private data member
    unordered_map<int64, Object*>* getAllObjects();
    unordered_map<int64, Object*>* getMobileObjects();
    unordered_map<Object*, Bundle*>* getLinePool();
    void examineRate(); /* computes matrix rate in all points*/
    /* output file functions */
    void displayDamage();
    void displayAllObject();
    void writeFile(const double, const long int);
    /* testing functions which includes display/output file functions that might be used when debugging */
    void writeSinkFile(const double, const long int);
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

    /* Parallel functions */
    void setDomain(int, int);
    void computeSectorRate(int); /* Compute the total rate in the specified sector for this processor */
    int computeSectorNumReactions(int); /* Compute the total number of reactions in the specified sector for this processor */
    long double getSectorRate(int); /* Get the total rate in the specified sector */ 
    long double findMaxAvgSectorRate();  /* Finds the average reaction rate in each sector, and returns the largest one */
    Object* selectSectorReaction(int64&, Reaction&, int&, int); /* Select a reaction inside this specified sector */
    int getStartIdx();  /* Returns the first spatial element's index that this processor is responsible for */
    int getEndIdx();    /* Returns the last spatial element's index that this processor is responsible for */
    void processBoundaryChanges(list<BoundaryChange*>*); /* Updates boundary configuration based on changes from neighbouring processors */
    list<BoundaryChange*>* getLeftBoundaryChanges();     /* Returns the changes on the left boundary of this processor's domain */
    list<BoundaryChange*>* getRightBoundaryChanges();    /* Returns the changes on the right boundary of this processor's domain */
    void combineProcessors(const vector<SCDWrapper*> &);   /* Usually for the master srscd object to combine information from each thread */
};
