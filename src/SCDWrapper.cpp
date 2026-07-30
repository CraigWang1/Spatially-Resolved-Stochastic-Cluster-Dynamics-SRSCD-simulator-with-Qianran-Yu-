#include<ctime>
#include<cstdlib>
#include "SCDWrapper.h"
// SCDWrapper -- implementaions of class SCDWrapper

static double fluenceH = 0.0;
static double in_time = 0.0;
static double doseIon[POINTS] = {0.0};
static int annilV = 0;
static int start = 1;
// static int plotTime = 1;
// static int plotTime1 = 1;
// static int plotTime2 = 1;
static int event = 0;
static int sinkV = 0;
// static int sinkH = 0;
static int generationNumber = 0;
static int dissV = 0; //this only counts number of events
static int dissH = 0; //this only counts number of events
/* public funciton */
SCDWrapper::SCDWrapper():allObjects(), engine(rd()), distribution(0.0L, 1.0L), damage(allObjects), cpdf(), startIndex(0), endIndex(0), matrixRateTree(POINTS), totalDpa(0)
{
    formationE[1] = V_FORM_E; 

    /* initialize sink numbers */
    setSinks();

    for (int i = 0; i < POINTS; ++i) {
        computeMatrixRate(i);
        totalVacInElement[i] = 0;
    } /* initialized matrix rate in every element */
    computeBulkRate();  /* initialized total rate in the bulk */
    
    /* initialize reactions[][] */
    for(int i=0; i<8; i++){
        for(int j=0; j<POINTS; j++){
            reactions[i][j] = 0;
        }
    }

    lastElemSaturated = false;
    numHDesorbed = 0;

    /* Set-up damage parameters */
    int maxClusterSize = 500;
    double vacExp = 1.6;    // inverse power law exponents, Sicong He 2025
    double siaExp = 1.8;
    vacClusterSizeCDF.resize(maxClusterSize);   // starts counting from cluster size 1
    siaClusterSizeCDF.resize(maxClusterSize);
    for (int size = 1; size <= maxClusterSize; size++)
    {
        int idx = size - 1;
        vacClusterSizeCDF[idx] = 1.0 / pow(size, vacExp);
        siaClusterSizeCDF[idx] = 1.0 / pow(size, siaExp);
        if (size > 1)
        {
            vacClusterSizeCDF[idx] += vacClusterSizeCDF[idx - 1];
            siaClusterSizeCDF[idx] += siaClusterSizeCDF[idx - 1];
        }
    }
    // Normalize cdf for max value to be = 1
    double maxVacCDFValue = vacClusterSizeCDF.back();
    double maxSiaCDFValue = siaClusterSizeCDF.back();
    for (int i = 0; i < maxClusterSize; i++)
    {
        vacClusterSizeCDF[i] /= maxVacCDFValue;
        siaClusterSizeCDF[i] /= maxSiaCDFValue;
    }

    selectReactionFile.open("selectReaction.txt", ios::app);
    processEventFile.open("Reactions.txt", ios::app);
    desorbedFile.open("Desorbed.txt", ios::out);
    /* set format of gs --- species */
    /*
    gs.set_title("Species");
    gs.set_xlabel("#(SIA/V)");
    gs.set_ylabel("#(H)");
    gs.set_xautoscale();
    gs.set_yautoscale();
    gs.cmd("set xtics 1\n");
    gs.cmd("set ytics 1\n");
    */
    /* set format of gs --- reactions */
    /*
    gr.set_title("Reactions");
    gr.set_xlabel("Element Number");
    gr.set_ylabel("Reaction Type");
    gr.set_yrange(-1,8);
    gr.set_xrange(-1,POINTS);
    gr.cmd("set xtics 1\n");
    gr.cmd("set ytics(\"DF\" 0,\"DB\" 1, \"SINK\" 2, \"DIS\" 3, \"COM\" 4, \"ION\" 5, \"H\" 6, \"SAV\" 7)\n");
    */
    /* set format of gd1 -- size distribution */
    /*
    gd1.set_title("Size Distribution");
    gd1.set_xlabel("size(nm)");
    gd1.set_ylabel("C(#/cm^3)");
    gd1.set_xrange(-6, 15);
    gd1.set_yautoscale();
    gd1.cmd("set xtics 1\n");
    */
    /* set format of gd2 ---depth distribution */
    /*
    gd2.set_title("Depth Distribution");
    gd2.set_xlabel("Depth(nm)");
    gd2.set_ylabel("C(#/cm^3)");
    gd2.set_yautoscale();
    gd2.cmd("set xtics 72\n");
    */
    /* set format of gh1 ---depth distribution */
    /*
    gh1.set_title("Depth Distribution");
    gh1.set_xlabel("Depth(nm)");
    gh1.set_ylabel("C(#/cm^3)");
    gh1.set_ylogscale();
    //gh1.set_yautoscale();
    gh1.cmd("set xtics 72\n");
    */
    /* set format of gh2 --  H to V ratio */
    /*
    gh2.set_title("H to V ratio");
    gh2.set_xlabel("time(s)");
    gh2.set_ylabel("ratio");
    gh2.cmd("set grid\n");
    gh2.cmd("set ytics 0.5\n");
    gh2.set_yautoscale();
    gh2.set_xautoscale();
    //gh2.cmd("set xtics 250\n");
    */
}

void SCDWrapper::computeMatrixRate(const int n)
{
    /* Qianran 0925 */
    //cout << "Element " << n + 1 << endl;
    matrixRate[n] = 0.0;

    unordered_map<int64, Object*>& objectsInThisElement = objectsInElement[n];
    if (objectsInThisElement.size() > 0)
    {
        unordered_map<int64, Object*>::iterator iter;
        for (iter = objectsInThisElement.begin(); iter != objectsInThisElement.end(); ++iter) {
            Object* tempObject = iter->second;
            OneLine* tempLine = tempObject->lines[n];
            if (tempLine != nullptr) {
                matrixRate[n] += tempLine->computeTotalRate();
                //tempLine->display(tempObject);/* Qianran 0925 */
            }
        }
    }
    matrixRate[n] += damage.getTotalDamage(n);
    matrixRate[n] += sinkDissRateDislocationScrew[n][0];
    matrixRate[n] += sinkDissRateDislocationScrew[n][1];
    matrixRate[n] += sinkDissRateDislocationEdge[n][0];
    matrixRate[n] += sinkDissRateDislocationEdge[n][1];
    matrixRate[n] += sinkDissRateGrainBndry[n][0];
    matrixRate[n] += sinkDissRateGrainBndry[n][1];

    matrixRateTree.set_val(n, matrixRate[n]);
}

void SCDWrapper::updateMatrixRate(const int n, const Reaction reaction)
{
    /* Qianran 0925 */
    int affectedStart = n - 1;
    int affectedEnd = n + 1;
    if (reaction == DIFFUSETOF)
        affectedStart -= 1;
    else if (reaction == DIFFUSETOB)
        affectedEnd += 1;
    if (affectedStart < 0)
        affectedStart = 0;
    if (affectedEnd >= POINTS)
        affectedEnd = POINTS - 1;

    for (int i = affectedStart; i <= affectedEnd; i++) {
        computeMatrixRate(i);
    }
}

void SCDWrapper::computeBulkRate()
{
    // bulkRate = 0.0;
    // for (int i = 0; i < POINTS; ++i) {
    //     bulkRate += matrixRate[i];
    // }
    bulkRate = matrixRateTree.sum(0, POINTS - 1);
}

void SCDWrapper::computeDomainRate()
{
    // domainRate = 0.0;
    // for (int i = startIndex; i <= endIndex; i++)
    // {
    //     domainRate += matrixRate[i];
    // }
    domainRate = matrixRateTree.sum(startIndex, endIndex);
}

void SCDWrapper::examineDomainRate()
{
    for (int i = startIndex; i <= endIndex; i++)
    {
        computeMatrixRate(i);
    }
    computeDomainRate();
}

long double SCDWrapper::getBulkRate()
{
    return bulkRate;
}

long double SCDWrapper::getDomainRate()
{
    return domainRate;
}

Object* SCDWrapper::selectDomainReaction(
                                   int64& theOtherKey,
                                   Reaction& reaction,
                                   int& count)
{
    int pointIndex;
    
    // Generate a random long double
    long double randomNum = distribution(engine);
    
    long double randRate = randomNum * (domainRate+noneRate);
    long double tempRandRate = randRate;
    Object* tempObject = nullptr;
    OneLine* tempLine;
    //fs << "BulkRate = " << bulkRate << "RandRate = " << randRate << endl;

    // Select which spatial element our reaction is in
    pointIndex = matrixRateTree.first_prefix_at_least_from(startIndex, tempRandRate);

    // If the reaction isn't in our domain (if running in parallel)
    if (pointIndex > endIndex || tempRandRate == 0)
    {
        reaction = NONE;
        return tempObject;
    }

    // Remove rates of other spatial elements
    if (pointIndex > startIndex)
        tempRandRate -= matrixRateTree.sum(startIndex, pointIndex - 1);

    // Select the reaction inside of our spatial element
    reaction = NONE;
    unordered_map<int64, Object*>& objectsInThisElement = objectsInElement[pointIndex];
    if (objectsInThisElement.size() > 0)
    {
        unordered_map<int64, Object*>::iterator iter = objectsInThisElement.begin();
        while (reaction == NONE && iter != objectsInThisElement.end()) {
            tempObject = iter->second;
            tempLine = tempObject->lines[pointIndex];
            if (tempLine != nullptr) {
                reaction = tempLine->selectReaction(tempObject, theOtherKey, tempRandRate);
            }
            ++iter;
        }
    }
    if (reaction == NONE) {
        reaction = damage.selectDamage(pointIndex, tempRandRate);
    }
    if (reaction == NONE){
        if (sinkDissRateDislocationScrew[pointIndex][0] >= tempRandRate){
            reaction = DISSVDISLOCATIONSCREW;
        }else{
            tempRandRate -= sinkDissRateDislocationScrew[pointIndex][0];
        }
    }
    if (reaction == NONE){
        if (sinkDissRateDislocationScrew[pointIndex][1] >= tempRandRate){
            reaction = DISSHDISLOCATIONSCREW;
        }else{
            tempRandRate -= sinkDissRateDislocationScrew[pointIndex][1];
        }
    }
    if (reaction == NONE){
        if (sinkDissRateDislocationEdge[pointIndex][0] >= tempRandRate){
            reaction = DISSVDISLOCATIONEDGE;
        }else{
            tempRandRate -= sinkDissRateDislocationEdge[pointIndex][0];
        }
    }
    if (reaction == NONE){
        if (sinkDissRateDislocationEdge[pointIndex][1] >= tempRandRate){
            reaction = DISSHDISLOCATIONEDGE;
        }else{
            tempRandRate -= sinkDissRateDislocationEdge[pointIndex][1];
        }
    }
    if (reaction == NONE){
        if (sinkDissRateGrainBndry[pointIndex][0] >= tempRandRate){
            reaction = DISSVGRAINBNDRY;
        }else{
            tempRandRate -= sinkDissRateGrainBndry[pointIndex][0];
        }
    }
    if (reaction == NONE){
        if (sinkDissRateGrainBndry[pointIndex][1] >= tempRandRate){
            reaction = DISSHGRAINBNDRY;
        }else{
            tempRandRate -= sinkDissRateGrainBndry[pointIndex][1];
        }
    }
    
    count = pointIndex;
    if (LOG_REACTIONS)
        selectReactionFile << "Element = " << pointIndex + 1 <<", "<< "Reaction = " << reaction << endl << endl;
    return tempObject;
}

void SCDWrapper::processEvent(
                              const Reaction reaction,
                              Object* hostObject,
                              const int n,
                              const int64 theOtherKey,
                              const double time,
                              const double dt
                              )
{
    if (reaction == NONE)
        return;

    ++event;
    switch (reaction) {
        case DIFFUSETOF:
            processDiffEvent(hostObject, n, 'f');
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() <<"  diffuses from element "<< n << " to element " << n-1 << endl;
            break;
        case DIFFUSETOB:
            processDiffEvent(hostObject, n, 'b');
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() <<"  diffuses from element "<< n << " to element " << n+1 << endl;
            break;
        case SINKDISLOCATIONSCREW:
            processSinkEvent(hostObject, n);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() <<"  in element "<< n << " goes to screw dislocation sink."<< endl;
            writeSinkFile(hostObject, n, time, reaction); /* this step also updated sinkDissRate */
            break;
        case SINKDISLOCATIONEDGE:
            processSinkEvent(hostObject, n);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() <<"  in element "<< n << " goes to edge dislocation sink."<< endl;
            writeSinkFile(hostObject, n, time, reaction); /* this step also updated sinkDissRate */
            break;
        case SINKGRAINBNDRY:
            processSinkEvent(hostObject, n);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() <<"  in element "<< n << " goes to grain boundary sink."<< endl;
            writeSinkFile(hostObject, n, time, reaction); /* this step also updated sinkDissRate */
            break;
        case DISSOCIATION:
            processDissoEvent(hostObject, n, theOtherKey, processEventFile);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() <<"  experiences a dissociation." << endl;
            break;
        case COMBINATION:
            processCombEvent(hostObject, n, theOtherKey, processEventFile);
            break;
        case SAV:
            processSAVEvent(hostObject, n);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() << "  in element " << n << " experiences SAV (ejects interstitial)" << endl;
            break;
        case RECOMBER:
            processRecombEvent(hostObject, n, true, time);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() << "  in element " << n << " experiences an ER recombination to form H2 and leave the front of the material" << endl;
            break;
        case RECOMBLH:
            processRecombEvent(hostObject, n, false, time);
            if (LOG_REACTIONS)
                processEventFile << hostObject->getKey() << "  in element " << n << " experiences a LH recombination to form H2 and leave the back of the material" << endl;
            break;
        case PARTICLE:
            getParticleInsertion(n, dt, processEventFile);
            break;
        case HE:
            getHeInsertion(n);
            break;
        case H:
            getHInsertion(n, dt, processEventFile);
            break;
        case DISSVDISLOCATIONSCREW:
        case DISSVDISLOCATIONEDGE:
        case DISSVGRAINBNDRY:
            processSinkDissEvent(0, n, reaction);
            dissV++;
            //cout << "dissV = " << dissV << endl;
            break;
        case DISSHDISLOCATIONSCREW:
        case DISSHDISLOCATIONEDGE:
        case DISSHGRAINBNDRY:
            processSinkDissEvent(1, n, reaction);
            dissH++;
            //cout << "dissH = " << dissH << endl;
            break;
        default:
            break;
    }

    // Keep track of affected reaction rates
    removeDestroyedObjects();
    updateMatrixRate(n, reaction);
    computeDomainRate();
}

SCDWrapper::~SCDWrapper()
{
    unordered_map<int64, Object*>::iterator iter1;
    /* clean contents of object* */
    for (iter1 = allObjects.begin(); iter1 != allObjects.end(); ++iter1) {
        delete iter1->second;
    }
    for (iter1 = HObjects.begin(); iter1 != HObjects.end(); ++iter1) {
        delete iter1->second;
    }

    selectReactionFile.close();
    processEventFile.close();
}

unordered_map<int64, Object*>* SCDWrapper::getAllObjects()
{
    return &allObjects;
}

unordered_map<int64, Object*>* SCDWrapper::getMobileObjects()
{
    return &mobileObjects;
}

void SCDWrapper::examineRate()
{
    for (int i = 0; i < POINTS; ++i) {
        computeMatrixRate(i);
    }
}

void SCDWrapper::writeSinkFile(const Object * const hostObject, const long int n, const double time, Reaction reaction)
{
    int i;
    for (i = 0; i < LEVELS; i++) {
        int number = hostObject->getAttri(i);
        if (i == 0 && number<0) {
            if (reaction == SINKDISLOCATIONSCREW)
                sinksDislocationScrew[n][i] += abs(number);
            else if (reaction == SINKDISLOCATIONEDGE)
                sinksDislocationEdge[n][i] += abs(number);
            else
                sinksGrainBndry[n][i] += abs(number);
            computeSinkDissRate(i, n);
        }
        else {
            if (reaction == SINKDISLOCATIONSCREW)
                sinksDislocationScrew[n][i + 1] += number;
            else if (reaction == SINKDISLOCATIONEDGE)
                sinksDislocationEdge[n][i + 1] += number;
            else
                sinksGrainBndry[n][i + 1] += number;
            if(i == 2){
                computeSinkDissRate(1, n);
            }
        }
    }
}

void SCDWrapper::writeSpeciesFile(const double time, const long int n, const int threadID)
{
    if (threadID == 0)
    {
        if (HYDROGEN_ON)
            fluenceH = FLUX_H * time; // Only take simulation progress parameters from thread id 0
        else
            fluenceH = 0;
    }

    ofstream fo;
    unordered_map<int64, Object*>::iterator iter;
    fo.open(std::string("species") + std::to_string(threadID) + std::string(".txt"), ios::out);
    fo << "step = " << n << endl;
    fo << "time = " << time << endl;
    fo << "fluenceH = " << fluenceH << endl;
    fo << "startIndex = " << startIndex << endl;
    fo << "endIndex = " << endIndex << endl;
    
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter) {
        fo << "object " << iter->second->getKey() << "    ";
        for(int i = 0; i < POINTS; i++){
            fo << (iter->second->getNumber(i)) << "    ";
            //totalDose += doseIon[i];
        }
        fo<<endl;
    }
    fo.close();
    /* write species.out file */
    
}

void SCDWrapper::writeClusterFile(const double time, const long int n)
{
    int i;
    ofstream fc;
    fc.open("clusters.out", ios::app);
    unordered_map<int64, Object*>::iterator iter;
    int sia[POINTS] = { 0 }, siac[POINTS] = { 0 }, siah[POINTS] = { 0 }, v[POINTS] = { 0 }, vc[POINTS] = { 0 }, vh[POINTS] = { 0 };
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter) {
        int attr0 = iter->second->getAttri(0);
        int attr2 = iter->second->getAttri(2);
        for (i = 0; i < POINTS; ++i) {
            int number = iter->second->getNumber(i);
            if (attr0 > 0) {
                if (attr0 == 1) {
                    sia[i] += number;
                }
                else {
                    siac[i] += number;
                    if (attr2 != 0) {
                        siah[i] += number;
                    }
                }
            }
            if (attr0 < 0) {
                if (attr0 == -1) {
                    v[i] += number;
                }
                else {
                    vc[i] += number;
                    if (attr2 != 0) {
                        vh[i] += number;
                    }
                }
            }
        }
    }
    fc << "Aggregate time = " << time << "  step = " << n << endl;
    for (i = 0; i < POINTS; ++i) {
        fc << "Element " << i + 1 << "    ";
        fc << fluenceH <<"  "<< doseIon[i]<<"   " << sia[i] / VOLUME << "    " << siac[i] / VOLUME << "    " << siah[i] / VOLUME << "    " << v[i] / VOLUME
        << "    " << vc[i] / VOLUME << "       " << vh[i] / VOLUME;
        fc << endl;
    }
    fc.close();
}

void SCDWrapper::writeSinkFile(const double time, const long int n, const int threadID)
{
    int i, j;
    ofstream outFile;
    outFile.open("sink" + std::to_string(threadID) + ".txt", ios::out);
    outFile << "startIndex = " << startIndex << endl;
    outFile << "endIndex = " << endIndex << endl;
    //fs << "Aggregate time = " << time << "  step = " << n << endl;
    for (i = 0; i < POINTS; ++i) {
        for (j = 0; j < LEVELS + 1; ++j) {
            outFile << sinksDislocationScrew[i][j]<< "    ";
        }
        for (j = 0; j < LEVELS + 1; ++j) {
            outFile << sinksDislocationEdge[i][j]<< "    ";
        }
        for (j = 0; j < LEVELS + 1; ++j) {
            outFile << sinksGrainBndry[i][j]<< "    ";
        }
        outFile << endl;
    }
    outFile.close();
}

void SCDWrapper::displayDamage(){
    for(int i = 0; i < POINTS; i++){
        damage.display(i);
    }
}

void SCDWrapper::displayAllObject(){
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter) {
        std::cout << iter->second->getKey() << "    ";
        for(int i = 0; i < POINTS; i++){
            std::cout << (iter->second->getNumber(i)) << "    " ;
        }
        std::cout<<endl<<endl;
    }
}

void SCDWrapper::writeFile(const double time, const long int n, const int threadID)
{
    writeSpeciesFile(time, n, threadID);
    writeSinkFile(time, n, threadID);
    //writeClusterFile(time, n);
}

void SCDWrapper::setSinks()
{
    int i, j;
    for (i = 0; i < POINTS; ++i){
        for (j = 0; j < LEVELS + 1; ++j) {
            sinksDislocationScrew[i][j] = 0;
            sinksDislocationEdge[i][j] = 0;
            sinksGrainBndry[i][j] = 0;
        }
    }
    for (i = 0; i < POINTS; ++i){
        for (j = 0; j < 2; ++j) {
            sinkDissRateDislocationScrew[i][j] = 0;
            sinkDissRateDislocationEdge[i][j] = 0;
            sinkDissRateGrainBndry[i][j] = 0;
        }
    }
    for (i = 0; i < POINTS; ++i){
        computeSinkDissRate(0, i);    // vacancy diss
        computeSinkDissRate(1, i);    // H diss
    }
}

void SCDWrapper::computeSinkDissRate(const int type, const int point)
{
    if (point == SURFACE_INDEX 
        || point == SUBSURFACE_INDEX
        || (point == BACK_SUBSURFACE_INDEX && BACK_DESORB)
        || (point == BACK_SURFACE_INDEX && BACK_DESORB))
    {
        sinkDissRateDislocationScrew[point][type] = 0;
        sinkDissRateDislocationEdge[point][type] = 0;
        sinkDissRateGrainBndry[point][type] = 0;
        return;
    }

    double ebHDislocationScrew = 0.55, ebHDislocationEdge = 0.89, ebHGrainBndry = 0.86; //binding and migration energy of hydrogen

    // vacancy emission (neglect vac emission from sinks)
    if(type == 0)
    {
        // Dislocations and grain boundaries are always a source of vacancy emission
        // vacVolTerm = max(1.0-totalVacInElement[point]*avol/volume, 0.);                // so that the mesh element doesn't become 100% vac
        // sinkDissRateDislocationScrew[point][type] = 2.0*PI*volume*DISLOCATION*(1-EDGE_DISLOCATION_FRAC)/b*NU0*exp(-(ebVDislocation+vacMigrationEnergy)/KB/TEMPERATURE)*vacVolTerm;
        // sinkDissRateDislocationEdge[point][type] = 2.0*PI*volume*DISLOCATION*EDGE_DISLOCATION_FRAC/b*NU0*exp(-(ebVDislocation+vacMigrationEnergy)/KB/TEMPERATURE)*vacVolTerm;
        // sinkDissRateGrainBndry[point][type] = 6.0*volume/GRAIN_SIZE/b/b*NU0*exp(-(ebVGrainBndry+vacMigrationEnergy)/KB/TEMPERATURE)*vacVolTerm;
    }
    // hydrogen emission
    else if(type == 1){
        if (sinksDislocationScrew[point][3] > 0)
            sinkDissRateDislocationScrew[point][type] = NU0*exp(-(ebHDislocationScrew+H_MIGRATION_ENERGY)/KB/TEMPERATURE)*sinksDislocationScrew[point][3];
        else
            sinkDissRateDislocationScrew[point][type] = 0;

        if (sinksDislocationEdge[point][3] > 0)
            sinkDissRateDislocationEdge[point][type] = NU0*exp(-(ebHDislocationEdge+H_MIGRATION_ENERGY)/KB/TEMPERATURE)*sinksDislocationEdge[point][3];
        else
            sinkDissRateDislocationEdge[point][type] = 0;

        if (sinksGrainBndry[point][3] > 0)
            sinkDissRateGrainBndry[point][type] = NU0*exp(-(ebHGrainBndry+H_MIGRATION_ENERGY)/KB/TEMPERATURE)*sinksGrainBndry[point][3];
        else
            sinkDissRateGrainBndry[point][type] = 0;
    }
}

int64 SCDWrapper::atomProperty(SCDWrapper::InsertStyle mode, const int n)
{
    if (mode == SUBSTITUTIONAL) {
        return (-1) * ((int64)pow(10.0, (double)EXP10 * (LEVELS - 1)) + (int64)pow(10.0, (double)EXP10 * (LEVELS - n)));
    }
    else if (mode == INTERSTITIAL) {
        return (int64)pow(10.0, (double)EXP10 * (LEVELS - n));
    }
    return 0; /* If return 0, means the code is wrong */
}

void SCDWrapper::addNewObjectToMap(Object* newObject)
{
    if (newObject->getKey() == 0) {
        delete newObject;
    }/* if this object is invalid, just delete this object */
    else {
        pair<int64, Object*> newNode(newObject->getKey(), newObject);
        allObjects.insert(newNode); /* add to all object */
        if (newObject->getDiff() > 0) {
            mobileObjects.insert(newNode);
        } /* add to mobile object if necessary */
        if (newObject->getAttri(0) == 0 && newObject->getAttri(2) > 0) {
            HObjects.insert(newNode);
        } /* Keep track of nH objects */
    }/* if this object is valid, add it to map */
}

void SCDWrapper::addToObjectMap(const int64 key, const int n, const int number)
{
    // If the object is nothing (eg. product of 1V + 1SIA comb), don't process it
    // For uniform free H concentration mode, don't touch the existing 1H, it's needed there for calculating rates
    if (key == 0 
        || number == 0
        || (UNIFORM_FREE_H_ON && key == 1))
    {
        return;
    }

    Object* anObject;

    bool objExists = allObjects.find(key) != allObjects.end();
    if (objExists)
    {
        anObject = allObjects[key];
    }

    if ( (objExists && anObject->getNumber(n) <= 0 && number < 0) 
        || (!objExists && number <= 0) )
    {
        return;
    }

    /* If the object exists, add to it. Otherwise create the object. */
    if (objExists) 
    {
        /* object found! then number of instances in this element increases by number*/
        anObject->addNumber(n, number);
    }
    else
    {
        /* object wasn't found! build new object and insert it into map */
        anObject = new Object(key, n, number);
        addNewObjectToMap(anObject);
    }
    if (anObject->getAttri(0) < 0)
    {
        int VPerCluster = abs(anObject->getAttri(0));
        totalVacInElement[n] += VPerCluster * number;
    }

    // If this object exists in this spatial element, account for it
    if (anObject->getNumber(n) > 0)
    {
        objectsInElement[n][key] = anObject;
    }
    else
    {
        objectsInElement[n].erase(key);
    }

    updateObjectInMap(anObject, n);

    affectedObjects.insert(key);
}

void SCDWrapper::reduceFromObjectMap(const int64 key, const int n, const int number)
{
    /* Count is mesh element index, number is number to add to object population */
    addToObjectMap(key, n, -number);
}

void SCDWrapper::updateObjectInMap(Object * hostObject, const int count)
{
    // Update the OneLine associated with this object in this count
    OneLine* tempLine = hostObject->lines[count];
    double diffusivity = hostObject->getDiff();
    int number = hostObject->getNumber(count);
    if (tempLine != nullptr) {
        if (number > 0) {
            tempLine->updateLine(hostObject, count, mobileObjects, allObjects);
        }
        else {
            delete tempLine;
            hostObject->lines[count] = nullptr;
        }
    }
    else {
        if (number > 0) {
            tempLine = new OneLine(hostObject, count, mobileObjects, allObjects);
            hostObject->lines[count] = tempLine;
        }
    }

    // Update the OneLines of other objects impacted by this mobile object
    if (diffusivity > 0) {
        updateRateToOther(hostObject, count);

        // Update diffusion rates of this object in neighbouring elements
        if((count-1) >= 0){
            OneLine* tempLine = hostObject->lines[count - 1];
            if(tempLine != nullptr){
                tempLine->updateDiff(hostObject, count - 1, allObjects);
            }
        }
        if((count + 1) < POINTS){
            OneLine* tempLine = hostObject->lines[count + 1];
            if(tempLine != nullptr){
                tempLine->updateDiff(hostObject, count + 1, allObjects);
            }
        }
    }

    // If the sink diss rates depend on free species count, uncomment the below lines
    // if (hostObject->getKey() == -1000000)
    //     computeSinkDissRate(0, count);
    // else if (hostObject->getKey() == 1)
    //     computeSinkDissRate(1, count);
}

void SCDWrapper::updateRateToOther(Object const * const mobileObject, const int count)
{
    unordered_map<int64, Object*>& objectsInThisElement = objectsInElement[count];
    
    if (objectsInThisElement.size() > 0)
    {
        unordered_map<int64, Object*>::iterator iter;
        for (iter = objectsInThisElement.begin(); iter != objectsInThisElement.end(); ++iter)
        {
            Object* hostObject = iter->second;
            OneLine* tempLine = hostObject->lines[count];
            if (tempLine != nullptr) {
                /* If both are mobile objects, mobileObject will have already recorded the combination rate, no need to record it again */
                if (hostObject->getDiff() > 0 && mobileObject->getKey() != hostObject->getKey()) {
                    tempLine->setCombReaction(mobileObject->getKey(), 0.0);
                }
                else {
                    tempLine->updateReaction(hostObject, mobileObject, allObjects, count);
                }
            }
        }
    }
}

void SCDWrapper::removeDestroyedObjects()
{
    unordered_set<int64>::iterator iter;
    for (iter = affectedObjects.begin(); iter != affectedObjects.end(); ++iter)
    {
        int64 key = *iter;
        if (allObjects[key]->getTotalNumber() <= 0)
        {
            removeObjectFromMap(key);
        }
    }
    affectedObjects.clear();
}

void SCDWrapper::removeObjectFromMap(const int64 deleteKey)
{
    Object* deleteObject = allObjects[deleteKey];
    double diffusivity = deleteObject->getDiff();
    bool is_nH = deleteObject->getAttri(0) == 0 && deleteObject->getAttri(2) > 0;  // if it's an nH object
    delete deleteObject;  /* delete the content of this object */
    allObjects.erase(deleteKey); /* delete this object from map allObjects */
    if (diffusivity > 0) {
        mobileObjects.erase(deleteKey); /* delete this object from map mobileObjects */
    }
    if (is_nH) {
        HObjects.erase(deleteKey);
    }
}

void SCDWrapper::updateSinks(const int point, const int* number){
    // Type refers to screw dislocation, edge dislocation, or grain boundary
    // Level refers to vacancy, sia, He, or hydrogen
    for (int type = 0; type < NUM_SINKS; type++)
    {
        for (int level = 0; level < LEVELS+1; level++)
        {
            if (type == 0) // dislocation
            {
                sinksDislocationScrew[point][level] = number[type*(LEVELS+1)+level];
            }
            else if (type == 1)
            {
                sinksDislocationEdge[point][level] = number[type*(LEVELS+1)+level];
            }
            else // grain boundary
            {
                sinksGrainBndry[point][level] = number[type*(LEVELS+1)+level];
            }
        }
    }
    computeSinkDissRate(0, point);  // vacancy diss
    computeSinkDissRate(1, point);  // H diss
}

/* private function */
void SCDWrapper::processDiffEvent(Object* hostObject, const int n, const char signal)
{
    int64 key = hostObject->getKey();
    reduceFromObjectMap(key, n);

    if (signal == 'f') {
        ++reactions[0][n];
        if(n != 0){ /* when not surface */
            /* diffuse to the previous element */
            addToObjectMap(key, n - 1);            
        }else{
            //surface diffuse to vacuum
            if(surface.find(key) != surface.end()){
                ++surface[key];
                
            }else{
                pair<int64, int> newNode(key, 1);
                surface.insert(newNode); /* add to all object */
            } 
        }
    }
    else if(signal == 'b'){
        ++reactions[1][n];
        if ((n + 1) != POINTS) {
            /* diffuse to the latter element */
            addToObjectMap(key, n + 1);
        }else{
            // bottom diffuses to vacuum
            if(bottom.find(key) != bottom.end()){
                ++bottom[key];
            }else{
                pair<int64, int> newNode(key, 1);
                bottom.insert(newNode); /* add to all object */
            }
        }
    }
}

void SCDWrapper::processSinkEvent(Object * hostObject, const int n)
{
    ++reactions[2][n];
    reduceFromObjectMap(hostObject->getKey(), n);
} /* if sink rate != 0, diffusivity of this object is not 0 */

void SCDWrapper::processDissoEvent(
                                   Object * hostObject,
                                   const int n,
                                   const int64 monomerKey,
                                   fstream& fs)
{
    ++reactions[3][n];
    int64 HKey = 1;
    int number = 1;
    int theOtherAttr[LEVELS] = { 0 };   /* this holds the attribute of the other cluster(product) */
    
    /* deal with the host object */
    reduceFromObjectMap(hostObject->getKey(), n);

    /* deal with monomer */
    addToObjectMap(monomerKey, n);

    /* generate the other cluster */
    Object* monomer = allObjects[monomerKey];
    for (int i = 0; i < LEVELS; i++) {
        theOtherAttr[i] = hostObject->getAttri(i) - monomer->getAttri(i);
    }/* now I have the attribute for the other cluster */
    int64 theOtherKey = attrToKey(theOtherAttr);
    
    if(hostObject->getAttri(0)== -1 && hostObject->getAttri(2) >0 && monomerKey == -1000000){
        theOtherKey = HKey;
        number = hostObject->getAttri(2);
        /* V1-Hm -> Hm + V */
        /* change above reaction to V1-Hm -> m * H +V */
        
    }
    else if (hostObject->getAttri(0) == 1 && hostObject->getAttri(2) > 0 && monomerKey == 1000000){
        theOtherKey = HKey;
        number = hostObject->getAttri(2);
        /* SIA1-Hm -> Hm + SIA */
        /* change above reaction to SIA1-Hm -> m*H + SIA */
    }
    
    addToObjectMap(theOtherKey, n, number);

    /*
    if(n==0){
        fs1 <<"Dissociation: " << hostObject->getKey() << " -> " << theOtherKey << " + "<<monomerKey<<endl;
    }
    if(n==1){
        fs3 <<"Dissociation: " << hostObject->getKey() << " -> " << theOtherKey << " + "<<monomerKey<<endl;
    }if(n==2){
        fs5 <<"Dissociation: " << hostObject->getKey() << " -> " << theOtherKey << " + "<<monomerKey<<endl;
    }
     */
    if (LOG_REACTIONS)
        fs <<"Dissociation: " << hostObject->getKey() << " -> " << theOtherKey << " + "<<monomerKey<<" in Element "<<n<<endl;
}

void SCDWrapper::processCombEvent(
                                  Object * hostObject,
                                  const int n,
                                  const int64 theOtherKey,
                                  fstream& fs)
{
    ++reactions[4][n];
    /* 1. find the other reactant */
    Object* theOtherObject = allObjects[theOtherKey];
    int productAttr[LEVELS] = { 0 };  /* attributes of the product object */
    int64 productKey = 0;
    int64 SIAKey = 1000000;
    int64 HKey = 1;
    int number = 1; /*the number of this material added */
    for (int i = 0; i < LEVELS; ++i) {/* 2. get the attributes of the product */
        productAttr[i] = hostObject->getAttri(i) + theOtherObject->getAttri(i);
    } /* now I have the attribute of the product object */
    productKey = attrToKey(productAttr);
    
    /* update reactant */
    reduceFromObjectMap(hostObject->getKey(), n);
    reduceFromObjectMap(theOtherObject->getKey(), n);

    if(hostObject->getAttri(0)<0 && hostObject->getAttri(2)>0 && productAttr[0] > 0 && productAttr[2] > 0){
        /* Vn-Hm + xxx -> SIAp-Hq (n, m are not zero) */
        /* change above reaction to Vn-Hm + xxx -> SIAp + q*H */
        productKey = HKey;
        number = productAttr[2];
        addToObjectMap(productKey, n, number);

        productKey = productAttr[0] * SIAKey;
        number = 1;
        addToObjectMap(productKey, n, number);
    }else if(hostObject->getAttri(0)<0 && hostObject->getAttri(2)>0 && theOtherObject->getAttri(0)>0 && theOtherObject->getAttri(2)==0 && abs(hostObject->getAttri(0)) == theOtherObject->getAttri(0)){
        /* Vn-Hm + SIAn -> Hm (n, m are not zero) */
        /* change above reaction to Vn-Hm + SIAn -> m * H */
        productKey = HKey;
        number = hostObject->getAttri(2);
        addToObjectMap(productKey, n, number);
    }
    else if(hostObject->getAttri(0)<0 && hostObject->getAttri(2)>0 && theOtherObject->getAttri(0)>0 && theOtherObject->getAttri(2)==0 && productAttr[0] < 0 && productAttr[2] > 12*abs(productAttr[0])){
        /* Vn-Hm + SIAx -> V(n-x)-Hm (n, m, and x are not zero) */
        /* change above reaction to Vn-Hm + SIAx -> V(n-x)-H(max) + (excess)*H */
        int maxH = abs(productAttr[0]) * 12;  // 1 vacancy stores max 12 H
        int excessH = productAttr[2] - maxH;
        productAttr[2] = maxH;
        productKey = attrToKey(productAttr);

        // Create the new V(n-x)-H(max) object (eject the excess H above vacancy containment limit)
        addToObjectMap(productKey, n, number);

        // Create (excess)*H object
        addToObjectMap(HKey, n, excessH);
    }
    else{
        addToObjectMap(productKey, n, number);
    }

    int attrZeroA = hostObject->getAttri(0);
    int attrZeroB = theOtherObject -> getAttri(0);
    if(attrZeroA * attrZeroB < 0){
        /* vacancy, interstitial anniles */
        annilV += abs(attrZeroA) < abs(attrZeroB) ? abs(attrZeroA) : abs(attrZeroB);
    }
}

void SCDWrapper::processSAVEvent(Object* hostObject, const int n)
{
    /* 
     * Superabundant vacancy mechanism.
     * Eject an interstitial (which increases vacancy by 1)
     */
    // Eject interstitial
    int64 SIAKey = (int64)pow(10.0, (double)EXP10 * (LEVELS - 1)); /* Key for SIA. */
    addToObjectMap(SIAKey, n);

    // Generate 1 vacancy
    int productAttr[LEVELS] = { 0 };
    for (int i = 0; i < LEVELS; i++)
        productAttr[i] = hostObject->getAttri(i);

    productAttr[0] -= 1;

    int64 productKey = attrToKey(productAttr);

    addToObjectMap(productKey, n);

    reduceFromObjectMap(hostObject->getKey(), n);
}

void SCDWrapper::processRecombEvent(Object* hostObject, const int n, bool ER, double time)
{
    if (n != SURFACE_INDEX && (n != BACK_SURFACE_INDEX || !BACK_DESORB))
        cerr << "Recomb Error" << endl;

    /* 
     * Recombination: Two 1H instances combine to form H2 molecule, which
     * leaves the material through either the front or the back of the material 
     */

    // The implementation is easy: just remove one 1H object for ER Recomb and two 1H objects for LH Recomb
    if (ER)
        reduceFromObjectMap(hostObject->getKey(), n, 1);
    else
    {
        reduceFromObjectMap(hostObject->getKey(), n, 2);
        numHDesorbed += 2;
        if (abs(TEMP_INCREASE_RATE) > 0) // when doing thermal desorption
            writeDesorbedFile(time);
    }
}

void SCDWrapper::processSinkDissEvent(const int type, const int point, Reaction reaction)
{
    // dissV event
    int64 productKey = 0;
    if(type == 0)
    {
        if (reaction == DISSVDISLOCATIONSCREW)
            sinksDislocationScrew[point][0]--;
        else if (reaction == DISSVDISLOCATIONEDGE)
            sinksDislocationEdge[point][0]--;
        else
            sinksGrainBndry[point][0]--;
        productKey = -1000000; //1V
    }
    // dissH event
    else if(type ==1)
    {
        if (reaction == DISSHDISLOCATIONSCREW)
            sinksDislocationScrew[point][3]--;
        else if (reaction == DISSHDISLOCATIONEDGE)
            sinksDislocationEdge[point][3]--;
        else
            sinksGrainBndry[point][3]--;
        productKey = 1;
    }
    addToObjectMap(productKey, point);
    computeSinkDissRate(type, point);
}

void SCDWrapper::getElectronInsertion(const int n)
{
    /* Insert Frenkel pairs: */
    int64 SIAKey = (int64)pow(10.0, (double)EXP10 * (LEVELS - 1)); /* Key for SIA. */
    int64 vacancyKey = (-1) * SIAKey; /* Key for vacancy */
    addToObjectMap(SIAKey, n);
    addToObjectMap(vacancyKey, n);
}

void SCDWrapper::getNeutronInsertion(const int count)
{
    // log10 of Epka(eV) values vs. cdf of RB19J sample holder, from SPECTER
    // vector<double> log10Epkas = {-5.0000, -3.0000, -2.0000, -1.0000, 0.0000, 0.6990, 1.0000, 1.3010, 1.4771, 1.6021, 1.6990, 1.7782, 1.8451, 1.9031, 1.9542, 2.0000, 2.3010, 2.4771, 2.6021, 2.6990, 2.7782, 2.8451, 2.9031, 2.9542, 3.0000, 3.1761, 3.3010, 3.3979, 3.4771, 3.5441, 3.6021, 3.6532, 3.6990, 3.7404, 3.7782, 3.8129, 3.8451, 3.8751, 3.9031, 3.9294, 3.9542, 3.9777, 4.0000, 4.1761, 4.3010, 4.3979, 4.4771, 4.5441, 4.6021, 4.6532, 4.6990, 4.7404, 4.7782, 4.8129, 4.8451, 4.8751, 4.9031, 4.9294, 4.9542, 4.9777, 5.0000, 5.0792, 5.1461, 5.2041, 5.2553, 5.3010, 5.3424};
    // vector<double> EpkaCDF = {0.000000000E+00, 4.487262029E-02, 4.926943359E-01, 7.778088348E-01, 8.218171672E-01, 8.468413160E-01, 8.718453651E-01, 8.836338497E-01, 8.915471086E-01, 8.982785042E-01, 9.041154623E-01, 9.081313859E-01, 9.113171912E-01, 9.139804038E-01, 9.162064476E-01, 9.301958513E-01, 9.366760003E-01, 9.408507117E-01, 9.440003375E-01, 9.466886748E-01, 9.490785313E-01, 9.510312189E-01, 9.527849193E-01, 9.543577222E-01, 9.607142580E-01, 9.658215963E-01, 9.699953027E-01, 9.733409008E-01, 9.760161732E-01, 9.781909627E-01, 9.799707928E-01, 9.815335458E-01, 9.829093715E-01, 9.840580704E-01, 9.850871760E-01, 9.860221143E-01, 9.868329370E-01, 9.875701946E-01, 9.882468516E-01, 9.888702444E-01, 9.894376594E-01, 9.899579406E-01, 9.936944782E-01, 9.959506715E-01, 9.974239808E-01, 9.983597232E-01, 9.989551773E-01, 9.993204897E-01, 9.995423906E-01, 9.996847971E-01, 9.997787834E-01, 9.998402383E-01, 9.998843270E-01, 9.999157931E-01, 9.999364757E-01, 9.999517214E-01, 9.999630375E-01, 9.999717155E-01, 9.999781103E-01, 9.999829483E-01, 9.999935810E-01, 9.999972784E-01, 9.999987848E-01, 9.999994797E-01, 9.999998196E-01, 9.999999595E-01, 1.000000000E+00};

    vector<double> Epkas = {7.50E-06,1.02E-05,1.07E-05,1.12E-05,1.18E-05,1.23E-05,1.29E-05,1.35E-05,1.41E-05,1.48E-05,1.55E-05,1.62E-05,1.70E-05,1.78E-05,1.86E-05,1.95E-05,2.04E-05,2.14E-05,2.24E-05,2.34E-05,2.46E-05,2.57E-05,2.69E-05,2.82E-05,2.95E-05,3.09E-05,3.24E-05,3.39E-05,3.55E-05,3.72E-05,3.89E-05,4.07E-05,4.27E-05,4.47E-05,4.68E-05,4.90E-05,5.13E-05,5.37E-05,5.62E-05,5.89E-05,6.17E-05,6.46E-05,6.76E-05,7.08E-05,7.42E-05,7.76E-05,8.13E-05,8.51E-05,8.91E-05,9.34E-05,9.77E-05,1.02E-04,1.07E-04,1.12E-04,1.18E-04,1.23E-04,1.29E-04,1.35E-04,1.41E-04,1.48E-04,1.55E-04,1.62E-04,1.70E-04,1.78E-04,1.86E-04,1.95E-04,2.04E-04,2.14E-04,2.24E-04,2.34E-04,2.46E-04,2.57E-04,2.69E-04,2.82E-04,2.95E-04,3.09E-04,3.24E-04,3.39E-04,3.55E-04,3.72E-04,3.89E-04,4.07E-04,4.27E-04,4.47E-04,4.68E-04,4.90E-04,5.13E-04,5.37E-04,5.62E-04,5.89E-04,6.17E-04,6.46E-04,6.76E-04,7.08E-04,7.42E-04,7.76E-04,8.13E-04,8.51E-04,8.91E-04,9.34E-04,9.77E-04,1.02E-03,1.07E-03,1.12E-03,1.18E-03,1.23E-03,1.29E-03,1.35E-03,1.41E-03,1.48E-03,1.55E-03,1.62E-03,1.70E-03,1.78E-03,1.86E-03,1.95E-03,2.04E-03,2.14E-03,2.24E-03,2.34E-03,2.46E-03,2.57E-03,2.69E-03,2.82E-03,2.95E-03,3.09E-03,3.24E-03,3.39E-03,3.55E-03,3.72E-03,3.89E-03,4.07E-03,4.27E-03,4.47E-03,4.68E-03,4.90E-03,5.13E-03,5.37E-03,5.62E-03,5.89E-03,6.17E-03,6.46E-03,6.76E-03,7.08E-03,7.42E-03,7.76E-03,8.13E-03,8.51E-03,8.91E-03,9.34E-03,9.77E-03,1.02E-02,1.07E-02,1.12E-02,1.18E-02,1.23E-02,1.29E-02,1.35E-02,1.41E-02,1.48E-02,1.55E-02,1.62E-02,1.70E-02,1.78E-02,1.86E-02,1.95E-02,2.04E-02,2.14E-02,2.24E-02,2.34E-02,2.46E-02,2.57E-02,2.69E-02,2.82E-02,2.95E-02,3.09E-02,3.24E-02,3.39E-02,3.55E-02,3.72E-02,3.89E-02,4.07E-02,4.27E-02,4.47E-02,4.68E-02,4.90E-02,5.13E-02,5.37E-02,5.62E-02,5.89E-02,6.17E-02,6.46E-02,6.76E-02,7.08E-02,7.42E-02,7.76E-02,8.13E-02,8.51E-02,8.91E-02,9.34E-02,9.77E-02,1.02E-01,1.07E-01,1.12E-01,1.18E-01,1.23E-01,1.29E-01,1.35E-01,1.41E-01,1.48E-01,1.55E-01,1.62E-01,1.70E-01,1.78E-01,1.86E-01,1.95E-01,2.04E-01,2.14E-01,2.24E-01,2.34E-01,2.46E-01,2.57E-01,2.69E-01,2.82E-01,2.95E-01,3.09E-01,3.24E-01,3.39E-01,3.55E-01,3.72E-01,3.89E-01,4.07E-01,4.27E-01,4.47E-01,4.68E-01,4.90E-01,5.13E-01,5.37E-01,5.63E-01,5.88E-01,6.13E-01,6.38E-01,6.63E-01,6.88E-01,7.13E-01,7.38E-01,7.63E-01,7.88E-01,8.13E-01,8.38E-01,8.63E-01,8.88E-01,9.13E-01,9.38E-01,9.63E-01,9.88E-01,1.01E+00,1.04E+00,1.06E+00,1.09E+00,1.11E+00,1.14E+00,1.16E+00,1.19E+00,1.21E+00,1.24E+00,1.26E+00,1.29E+00,1.31E+00,1.34E+00,1.36E+00,1.39E+00,1.41E+00,1.44E+00,1.46E+00,1.49E+00,1.51E+00,1.54E+00,1.56E+00,1.59E+00,1.61E+00,1.64E+00,1.66E+00,1.69E+00,1.71E+00,1.74E+00,1.76E+00,1.79E+00,1.81E+00,1.84E+00,1.86E+00,1.89E+00,1.91E+00,1.94E+00,1.96E+00,1.99E+00,2.01E+00,2.04E+00,2.06E+00,2.09E+00,2.11E+00,2.14E+00,2.16E+00,2.19E+00,2.21E+00,2.24E+00,2.26E+00,2.29E+00,2.31E+00,2.34E+00,2.36E+00,2.39E+00,2.41E+00,2.44E+00,2.46E+00,2.49E+00,2.51E+00,2.54E+00,2.56E+00,2.59E+00,2.61E+00,2.64E+00,2.66E+00,2.69E+00,2.71E+00,2.74E+00,2.76E+00,2.79E+00,2.81E+00,2.84E+00,2.86E+00,2.89E+00,2.91E+00,2.94E+00,2.96E+00,2.99E+00,3.01E+00,3.04E+00,3.06E+00,3.09E+00,3.11E+00,3.14E+00,3.16E+00,3.19E+00,3.21E+00,3.24E+00,3.26E+00,3.29E+00,3.31E+00,3.34E+00,3.36E+00,3.39E+00,3.41E+00,3.44E+00,3.46E+00,3.49E+00,3.51E+00,3.54E+00,3.56E+00,3.59E+00,3.61E+00,3.64E+00,3.66E+00,3.69E+00,3.71E+00,3.74E+00,3.76E+00,3.79E+00,3.81E+00,3.84E+00,3.86E+00,3.89E+00,3.91E+00,3.94E+00,3.96E+00,3.99E+00,4.01E+00,4.04E+00,4.06E+00,4.09E+00,4.11E+00,4.14E+00,4.16E+00,4.19E+00,4.21E+00,4.24E+00,4.26E+00,4.29E+00,4.31E+00,4.34E+00,4.36E+00,4.39E+00,4.41E+00,4.44E+00,4.46E+00,4.49E+00,4.51E+00,4.54E+00,4.56E+00,4.59E+00,4.61E+00,4.64E+00,4.66E+00,4.69E+00,4.71E+00,4.74E+00,4.76E+00,4.79E+00,4.81E+00,4.84E+00,4.86E+00,4.89E+00,4.91E+00,4.94E+00,4.96E+00,4.99E+00,5.01E+00,5.04E+00,5.06E+00,5.09E+00,5.11E+00,5.14E+00,5.16E+00,5.19E+00,5.21E+00,5.24E+00,5.26E+00,5.29E+00,5.31E+00,5.34E+00,5.36E+00,5.39E+00,5.41E+00,5.44E+00,5.46E+00,5.49E+00,5.51E+00,5.54E+00,5.56E+00,5.59E+00,5.61E+00,5.64E+00,5.66E+00,5.69E+00,5.71E+00,5.74E+00,5.76E+00,5.79E+00,5.81E+00,5.84E+00,5.86E+00,5.89E+00,5.91E+00,5.94E+00,5.96E+00,5.99E+00,6.01E+00,6.04E+00,6.06E+00,6.09E+00,6.11E+00,6.14E+00,6.16E+00,6.19E+00,6.21E+00,6.24E+00,6.26E+00,6.29E+00,6.31E+00,6.34E+00,6.36E+00,6.39E+00,6.41E+00,6.44E+00,6.46E+00,6.49E+00,6.51E+00,6.54E+00,6.56E+00,6.59E+00,6.61E+00,6.64E+00,6.66E+00,6.69E+00,6.71E+00,6.74E+00,6.76E+00,6.79E+00,6.81E+00,6.84E+00,6.86E+00,6.89E+00,6.91E+00,6.94E+00,6.96E+00,6.99E+00,7.01E+00,7.04E+00,7.06E+00,7.09E+00,7.11E+00,7.14E+00,7.16E+00,7.19E+00,7.21E+00,7.24E+00,7.26E+00,7.29E+00,7.31E+00,7.34E+00,7.36E+00,7.39E+00,7.41E+00,7.44E+00,7.46E+00,7.49E+00,7.51E+00,7.54E+00,7.56E+00,7.59E+00,7.61E+00,7.64E+00,7.66E+00,7.69E+00,7.71E+00,7.74E+00,7.76E+00,7.79E+00,7.81E+00,7.84E+00,7.86E+00,7.89E+00,7.91E+00,7.94E+00,7.96E+00,7.99E+00,8.01E+00,8.04E+00,8.06E+00,8.09E+00,8.11E+00,8.14E+00,8.16E+00,8.19E+00,8.21E+00,8.24E+00,8.26E+00,8.29E+00,8.31E+00,8.34E+00,8.36E+00,8.39E+00,8.41E+00,8.44E+00,8.46E+00,8.49E+00,8.51E+00,8.54E+00,8.56E+00,8.59E+00,8.61E+00,8.64E+00,8.66E+00,8.69E+00,8.71E+00,8.74E+00,8.76E+00,8.79E+00,8.81E+00,8.84E+00,8.86E+00,8.89E+00,8.91E+00,8.94E+00,8.96E+00,8.99E+00,9.01E+00,9.04E+00,9.06E+00,9.09E+00,9.11E+00,9.14E+00,9.16E+00,9.19E+00,9.21E+00,9.24E+00,9.26E+00,9.29E+00,9.31E+00,9.34E+00,9.36E+00,9.39E+00,9.41E+00,9.44E+00,9.46E+00,9.49E+00,9.51E+00,9.54E+00,9.56E+00,9.59E+00,9.61E+00,9.64E+00,9.66E+00,9.69E+00,9.71E+00,9.74E+00,9.76E+00,9.79E+00,9.81E+00,9.84E+00,9.86E+00,9.89E+00,9.91E+00,9.94E+00,9.96E+00,9.99E+00,1.02E+01,1.07E+01,1.12E+01,1.18E+01,1.23E+01,1.29E+01,1.35E+01,1.41E+01,1.48E+01,1.55E+01,1.62E+01,1.70E+01,1.78E+01,1.86E+01,1.95E+01,2.04E+01,2.14E+01,2.24E+01,2.34E+01,2.46E+01,2.57E+01,2.69E+01,2.82E+01,2.95E+01,3.09E+01,3.24E+01,3.39E+01,3.55E+01,3.72E+01,3.89E+01,4.07E+01,4.27E+01,4.47E+01,4.68E+01,4.90E+01,5.13E+01,5.37E+01,5.62E+01,5.89E+01,6.17E+01,6.46E+01,6.76E+01,7.08E+01,7.42E+01,7.76E+01,8.13E+01,8.51E+01,8.91E+01,9.34E+01,9.77E+01,1.02E+02,1.07E+02,1.12E+02,1.18E+02,1.23E+02,1.29E+02,1.35E+02,1.41E+02,1.48E+02,1.55E+02,1.62E+02,1.70E+02,1.78E+02,1.86E+02,1.95E+02,2.04E+02,2.14E+02,2.24E+02,2.34E+02,2.46E+02,2.57E+02,2.69E+02,2.82E+02,2.95E+02,3.09E+02,3.24E+02,3.39E+02,3.55E+02,3.72E+02,3.89E+02,4.07E+02,4.27E+02,4.47E+02,4.68E+02,4.90E+02,5.13E+02,5.37E+02,5.62E+02,5.89E+02,6.17E+02,6.46E+02,6.76E+02,7.08E+02,7.42E+02,7.76E+02,8.13E+02,8.51E+02,8.91E+02,9.34E+02,9.77E+02,1.02E+03,1.07E+03,1.12E+03,1.18E+03,1.23E+03,1.29E+03,1.35E+03,1.41E+03,1.48E+03,1.55E+03,1.62E+03,1.70E+03,1.78E+03,1.86E+03,1.95E+03,2.04E+03,2.14E+03,2.24E+03,2.34E+03,2.46E+03,2.57E+03,2.69E+03,2.82E+03,2.95E+03,3.09E+03,3.24E+03,3.39E+03,3.55E+03,3.72E+03,3.89E+03,4.07E+03,4.27E+03,4.47E+03,4.68E+03,4.90E+03,5.13E+03,5.37E+03,5.62E+03,5.89E+03,6.17E+03,6.46E+03,6.76E+03,7.08E+03,7.42E+03,7.76E+03,8.13E+03,8.51E+03,8.91E+03,9.34E+03,9.77E+03,1.02E+04,1.07E+04,1.12E+04,1.18E+04,1.23E+04,1.29E+04,1.35E+04,1.41E+04,1.48E+04,1.55E+04,1.62E+04,1.70E+04,1.78E+04,1.86E+04,1.95E+04,2.04E+04,2.14E+04,2.24E+04,2.34E+04,2.46E+04,2.57E+04,2.69E+04,2.82E+04,2.95E+04,3.09E+04,3.24E+04,3.39E+04,3.55E+04,3.72E+04,3.89E+04,4.07E+04,4.27E+04,4.47E+04,4.68E+04,4.90E+04,5.13E+04,5.37E+04,5.62E+04,5.89E+04,6.17E+04,6.46E+04,6.76E+04,7.08E+04,7.42E+04,7.76E+04,8.13E+04,8.51E+04,8.91E+04,9.34E+04,9.77E+04,1.02E+05,1.07E+05,1.12E+05,1.18E+05,1.23E+05,1.29E+05,1.35E+05,1.41E+05,1.48E+05,1.55E+05,1.62E+05,1.70E+05,1.78E+05,1.86E+05,1.95E+05,2.04E+05,2.14E+05,2.24E+05,2.34E+05,2.46E+05};
    vector<double> EpkaCDF = {0.0000E+00,2.0310E-05,2.1259E-05,2.2266E-05,2.3304E-05,2.4398E-05,2.5555E-05,2.6750E-05,2.8007E-05,2.9333E-05,3.0720E-05,3.2165E-05,3.3675E-05,3.5257E-05,3.6915E-05,3.8652E-05,4.0468E-05,4.2372E-05,4.4371E-05,4.6466E-05,4.8659E-05,5.0953E-05,5.3358E-05,5.5872E-05,5.8509E-05,6.1269E-05,6.4161E-05,6.7189E-05,7.0358E-05,7.3675E-05,7.7149E-05,8.0787E-05,8.4599E-05,8.8588E-05,9.2765E-05,9.7141E-05,1.0172E-04,1.0652E-04,1.1154E-04,1.1680E-04,1.2231E-04,1.2807E-04,1.3411E-04,1.4044E-04,1.4706E-04,1.5399E-04,1.6125E-04,1.6886E-04,1.7682E-04,1.8516E-04,1.9389E-04,2.0303E-04,2.1261E-04,2.2264E-04,2.3313E-04,2.4411E-04,2.5563E-04,2.6768E-04,2.8029E-04,2.9351E-04,3.0735E-04,3.2184E-04,3.3700E-04,3.5289E-04,3.6953E-04,3.8695E-04,4.0519E-04,4.2429E-04,4.4429E-04,4.6522E-04,4.8716E-04,5.1013E-04,5.3416E-04,5.5934E-04,5.8572E-04,6.1333E-04,6.4223E-04,6.7251E-04,7.0421E-04,7.3740E-04,7.7216E-04,8.0855E-04,8.4666E-04,8.8657E-04,9.2835E-04,9.7211E-04,1.0179E-03,1.0659E-03,1.1161E-03,1.1687E-03,1.2238E-03,1.2815E-03,1.3419E-03,1.4052E-03,1.4714E-03,1.5407E-03,1.6134E-03,1.6894E-03,1.7690E-03,1.8524E-03,1.9397E-03,2.0311E-03,2.1269E-03,2.2272E-03,2.3321E-03,2.4419E-03,2.5571E-03,2.6776E-03,2.8037E-03,2.9360E-03,3.0743E-03,3.2192E-03,3.3709E-03,3.5298E-03,3.6962E-03,3.8704E-03,4.0527E-03,4.2438E-03,4.4438E-03,4.6531E-03,4.8725E-03,5.1022E-03,5.3425E-03,5.5943E-03,5.8581E-03,6.1342E-03,6.4232E-03,6.7260E-03,7.0430E-03,7.3749E-03,7.7225E-03,8.0864E-03,8.4676E-03,8.8666E-03,9.2845E-03,9.7221E-03,1.0180E-02,1.0660E-02,1.1162E-02,1.1686E-02,1.2234E-02,1.2806E-02,1.3404E-02,1.4028E-02,1.4680E-02,1.5361E-02,1.6072E-02,1.6814E-02,1.7589E-02,1.8397E-02,1.9239E-02,2.0118E-02,2.1035E-02,2.1992E-02,2.2989E-02,2.4029E-02,2.5114E-02,2.6245E-02,2.7423E-02,2.8654E-02,2.9935E-02,3.1272E-02,3.2665E-02,3.4118E-02,3.5634E-02,3.7215E-02,3.8864E-02,4.0586E-02,4.2383E-02,4.4256E-02,4.6213E-02,4.8254E-02,5.0383E-02,5.2606E-02,5.4929E-02,5.7356E-02,5.9891E-02,6.2541E-02,6.5309E-02,6.8201E-02,7.1223E-02,7.4380E-02,7.7682E-02,8.1134E-02,8.4745E-02,8.8522E-02,9.2472E-02,9.6602E-02,1.0092E-01,1.0544E-01,1.1017E-01,1.1511E-01,1.2029E-01,1.2570E-01,1.3137E-01,1.3729E-01,1.4349E-01,1.4998E-01,1.5676E-01,1.6380E-01,1.7115E-01,1.7883E-01,1.8688E-01,1.9529E-01,2.0409E-01,2.1328E-01,2.2293E-01,2.3300E-01,2.4353E-01,2.5458E-01,2.6612E-01,2.7820E-01,2.9083E-01,3.0406E-01,3.1789E-01,3.3237E-01,3.4751E-01,3.6336E-01,3.7995E-01,3.9729E-01,4.1546E-01,4.3446E-01,4.5432E-01,4.7510E-01,4.9684E-01,5.1955E-01,5.4329E-01,5.6809E-01,5.9392E-01,6.2065E-01,6.4731E-01,6.5986E-01,6.6371E-01,6.6683E-01,6.6895E-01,6.7074E-01,6.7252E-01,6.7435E-01,6.7610E-01,6.7769E-01,6.7907E-01,6.8043E-01,6.8176E-01,6.8307E-01,6.8437E-01,6.8566E-01,6.8694E-01,6.8821E-01,6.8947E-01,6.9072E-01,6.9197E-01,6.9320E-01,6.9443E-01,6.9566E-01,6.9687E-01,6.9808E-01,6.9908E-01,6.9993E-01,7.0073E-01,7.0151E-01,7.0229E-01,7.0306E-01,7.0383E-01,7.0459E-01,7.0535E-01,7.0611E-01,7.0686E-01,7.0761E-01,7.0836E-01,7.0911E-01,7.0985E-01,7.1059E-01,7.1133E-01,7.1206E-01,7.1280E-01,7.1353E-01,7.1426E-01,7.1498E-01,7.1571E-01,7.1643E-01,7.1716E-01,7.1788E-01,7.1860E-01,7.1932E-01,7.2004E-01,7.2076E-01,7.2147E-01,7.2219E-01,7.2290E-01,7.2362E-01,7.2433E-01,7.2504E-01,7.2575E-01,7.2646E-01,7.2717E-01,7.2788E-01,7.2859E-01,7.2930E-01,7.3000E-01,7.3071E-01,7.3141E-01,7.3212E-01,7.3282E-01,7.3352E-01,7.3421E-01,7.3489E-01,7.3558E-01,7.3626E-01,7.3694E-01,7.3762E-01,7.3830E-01,7.3898E-01,7.3966E-01,7.4034E-01,7.4102E-01,7.4169E-01,7.4228E-01,7.4282E-01,7.4335E-01,7.4389E-01,7.4442E-01,7.4495E-01,7.4548E-01,7.4600E-01,7.4653E-01,7.4706E-01,7.4758E-01,7.4811E-01,7.4864E-01,7.4916E-01,7.4968E-01,7.5021E-01,7.5073E-01,7.5125E-01,7.5178E-01,7.5230E-01,7.5282E-01,7.5334E-01,7.5386E-01,7.5439E-01,7.5491E-01,7.5543E-01,7.5595E-01,7.5646E-01,7.5698E-01,7.5750E-01,7.5802E-01,7.5853E-01,7.5905E-01,7.5957E-01,7.6008E-01,7.6059E-01,7.6109E-01,7.6160E-01,7.6210E-01,7.6260E-01,7.6310E-01,7.6360E-01,7.6410E-01,7.6460E-01,7.6511E-01,7.6561E-01,7.6611E-01,7.6660E-01,7.6710E-01,7.6760E-01,7.6810E-01,7.6859E-01,7.6909E-01,7.6958E-01,7.7008E-01,7.7057E-01,7.7106E-01,7.7155E-01,7.7203E-01,7.7246E-01,7.7281E-01,7.7313E-01,7.7343E-01,7.7374E-01,7.7404E-01,7.7434E-01,7.7463E-01,7.7493E-01,7.7522E-01,7.7552E-01,7.7581E-01,7.7610E-01,7.7639E-01,7.7669E-01,7.7698E-01,7.7727E-01,7.7756E-01,7.7785E-01,7.7814E-01,7.7842E-01,7.7871E-01,7.7900E-01,7.7929E-01,7.7958E-01,7.7986E-01,7.8015E-01,7.8042E-01,7.8066E-01,7.8089E-01,7.8112E-01,7.8135E-01,7.8158E-01,7.8181E-01,7.8204E-01,7.8227E-01,7.8249E-01,7.8272E-01,7.8295E-01,7.8317E-01,7.8340E-01,7.8363E-01,7.8385E-01,7.8408E-01,7.8431E-01,7.8453E-01,7.8476E-01,7.8498E-01,7.8521E-01,7.8543E-01,7.8566E-01,7.8588E-01,7.8610E-01,7.8633E-01,7.8655E-01,7.8677E-01,7.8699E-01,7.8721E-01,7.8740E-01,7.8757E-01,7.8774E-01,7.8789E-01,7.8805E-01,7.8821E-01,7.8837E-01,7.8852E-01,7.8867E-01,7.8881E-01,7.8895E-01,7.8909E-01,7.8923E-01,7.8937E-01,7.8951E-01,7.8964E-01,7.8978E-01,7.8992E-01,7.9006E-01,7.9019E-01,7.9033E-01,7.9046E-01,7.9060E-01,7.9074E-01,7.9087E-01,7.9101E-01,7.9114E-01,7.9127E-01,7.9141E-01,7.9154E-01,7.9167E-01,7.9180E-01,7.9193E-01,7.9206E-01,7.9220E-01,7.9233E-01,7.9246E-01,7.9259E-01,7.9272E-01,7.9285E-01,7.9298E-01,7.9311E-01,7.9324E-01,7.9337E-01,7.9350E-01,7.9363E-01,7.9376E-01,7.9389E-01,7.9401E-01,7.9414E-01,7.9427E-01,7.9440E-01,7.9453E-01,7.9466E-01,7.9478E-01,7.9491E-01,7.9503E-01,7.9516E-01,7.9529E-01,7.9541E-01,7.9554E-01,7.9566E-01,7.9579E-01,7.9591E-01,7.9604E-01,7.9616E-01,7.9629E-01,7.9641E-01,7.9654E-01,7.9666E-01,7.9678E-01,7.9691E-01,7.9703E-01,7.9715E-01,7.9728E-01,7.9740E-01,7.9753E-01,7.9765E-01,7.9777E-01,7.9789E-01,7.9802E-01,7.9814E-01,7.9826E-01,7.9839E-01,7.9851E-01,7.9863E-01,7.9875E-01,7.9888E-01,7.9900E-01,7.9912E-01,7.9924E-01,7.9936E-01,7.9948E-01,7.9960E-01,7.9972E-01,7.9984E-01,7.9996E-01,8.0009E-01,8.0021E-01,8.0033E-01,8.0045E-01,8.0057E-01,8.0069E-01,8.0081E-01,8.0093E-01,8.0105E-01,8.0117E-01,8.0129E-01,8.0141E-01,8.0153E-01,8.0165E-01,8.0177E-01,8.0188E-01,8.0200E-01,8.0212E-01,8.0224E-01,8.0235E-01,8.0247E-01,8.0259E-01,8.0270E-01,8.0282E-01,8.0294E-01,8.0305E-01,8.0317E-01,8.0329E-01,8.0340E-01,8.0352E-01,8.0363E-01,8.0375E-01,8.0387E-01,8.0398E-01,8.0410E-01,8.0421E-01,8.0433E-01,8.0444E-01,8.0456E-01,8.0467E-01,8.0479E-01,8.0490E-01,8.0502E-01,8.0513E-01,8.0524E-01,8.0536E-01,8.0547E-01,8.0559E-01,8.0570E-01,8.0582E-01,8.0593E-01,8.0604E-01,8.0616E-01,8.0627E-01,8.0638E-01,8.0650E-01,8.0661E-01,8.0672E-01,8.0684E-01,8.0695E-01,8.0706E-01,8.0718E-01,8.0729E-01,8.0739E-01,8.0750E-01,8.0761E-01,8.0772E-01,8.0783E-01,8.0794E-01,8.0805E-01,8.0816E-01,8.0827E-01,8.0837E-01,8.0848E-01,8.0859E-01,8.0870E-01,8.0881E-01,8.0892E-01,8.0903E-01,8.0913E-01,8.0924E-01,8.0935E-01,8.0946E-01,8.0957E-01,8.0967E-01,8.0978E-01,8.1181E-01,8.1387E-01,8.1597E-01,8.1809E-01,8.2027E-01,8.2252E-01,8.2485E-01,8.2725E-01,8.2965E-01,8.3196E-01,8.3421E-01,8.3652E-01,8.3873E-01,8.4100E-01,8.4331E-01,8.4563E-01,8.4773E-01,8.4983E-01,8.5187E-01,8.5372E-01,8.5553E-01,8.5730E-01,8.5904E-01,8.6079E-01,8.6245E-01,8.6406E-01,8.6558E-01,8.6711E-01,8.6864E-01,8.7015E-01,8.7163E-01,8.7311E-01,8.7456E-01,8.7596E-01,8.7733E-01,8.7872E-01,8.8009E-01,8.8144E-01,8.8276E-01,8.8408E-01,8.8540E-01,8.8666E-01,8.8790E-01,8.8913E-01,8.9034E-01,8.9155E-01,8.9273E-01,8.9391E-01,8.9507E-01,8.9622E-01,8.9737E-01,8.9851E-01,8.9987E-01,9.0112E-01,9.0233E-01,9.0352E-01,9.0470E-01,9.0591E-01,9.0716E-01,9.0840E-01,9.0955E-01,9.1068E-01,9.1179E-01,9.1288E-01,9.1395E-01,9.1501E-01,9.1607E-01,9.1711E-01,9.1814E-01,9.1917E-01,9.2018E-01,9.2118E-01,9.2218E-01,9.2317E-01,9.2416E-01,9.2514E-01,9.2611E-01,9.2709E-01,9.2807E-01,9.2904E-01,9.3001E-01,9.3098E-01,9.3196E-01,9.3293E-01,9.3391E-01,9.3489E-01,9.3587E-01,9.3685E-01,9.3784E-01,9.3882E-01,9.3981E-01,9.4080E-01,9.4180E-01,9.4279E-01,9.4379E-01,9.4479E-01,9.4579E-01,9.4679E-01,9.4779E-01,9.4879E-01,9.4979E-01,9.5080E-01,9.5180E-01,9.5279E-01,9.5379E-01,9.5479E-01,9.5579E-01,9.5678E-01,9.5778E-01,9.5877E-01,9.5976E-01,9.6075E-01,9.6173E-01,9.6271E-01,9.6368E-01,9.6465E-01,9.6562E-01,9.6657E-01,9.6752E-01,9.6846E-01,9.6939E-01,9.7032E-01,9.7123E-01,9.7213E-01,9.7302E-01,9.7389E-01,9.7475E-01,9.7560E-01,9.7644E-01,9.7725E-01,9.7806E-01,9.7884E-01,9.7961E-01,9.8037E-01,9.8110E-01,9.8182E-01,9.8253E-01,9.8321E-01,9.8388E-01,9.8452E-01,9.8515E-01,9.8577E-01,9.8636E-01,9.8695E-01,9.8751E-01,9.8806E-01,9.8859E-01,9.8911E-01,9.8962E-01,9.9011E-01,9.9059E-01,9.9105E-01,9.9151E-01,9.9195E-01,9.9238E-01,9.9281E-01,9.9322E-01,9.9362E-01,9.9401E-01,9.9439E-01,9.9476E-01,9.9511E-01,9.9546E-01,9.9579E-01,9.9611E-01,9.9642E-01,9.9671E-01,9.9700E-01,9.9726E-01,9.9752E-01,9.9775E-01,9.9797E-01,9.9818E-01,9.9837E-01,9.9855E-01,9.9871E-01,9.9886E-01,9.9900E-01,9.9912E-01,9.9923E-01,9.9933E-01,9.9942E-01,9.9950E-01,9.9957E-01,9.9963E-01,9.9968E-01,9.9973E-01,9.9977E-01,9.9980E-01,9.9983E-01,9.9986E-01,9.9988E-01,9.9990E-01,9.9992E-01,9.9993E-01,9.9994E-01,9.9995E-01,9.9996E-01,9.9997E-01,9.9998E-01,9.9998E-01,9.9998E-01,9.9999E-01,9.9999E-01,9.9999E-01,9.9999E-01,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00,1.0000E+00};

    double randomNum = distribution(engine);
    size_t interpIdx = lower_bound(EpkaCDF.begin(), EpkaCDF.end(), randomNum) - EpkaCDF.begin();
    if (interpIdx >= EpkaCDF.size())
    {
        interpIdx = EpkaCDF.size() - 1;
    }

    double Epka;
    if (interpIdx == 0)
    {
        Epka = Epkas[0];
    }
    else // interpolate
    {
        Epka = Epkas[interpIdx-1] + (randomNum-EpkaCDF[interpIdx-1])/(EpkaCDF[interpIdx]-EpkaCDF[interpIdx-1]) * (Epkas[interpIdx]-Epkas[interpIdx-1]);
    }

    if (Epka < 620)
        return;

    double numFP;
    // if (Epka <= 48000)  // Sicong He 2025
    //     numFP = 3.81 * pow(Epka/1000.0, 0.62);
    // else
    //     numFP = 0.50 * pow(Epka/1000.0, 1.15);

    if (Epka <= 43000)   // Qianran Yu 2020
        numFP = 1.15e-2 * pow(Epka, 0.74);
    else
        numFP = 1.89e-5*pow(Epka, 1.34);

    // Fraction of vacs and interstitials that are clustered (from Byggmaster et al)
    // vector<double> pkaEnergyTable = {40,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,600,700,800,900,1000,2000,5000,10000,20000,50000,100000,200000,300000,500000,1000000,2000000}; // eV
    // vector<double> fracSIAinCluster = {0,0,0,0,0,0,0,0,0,0.006,0.01,0.078,0.103,0.218,0.298,0.32,0.308,0.305,0.414,0.425,0.458,0.445,0.533,0.565,0.611,0.674,0.795,0.841,0.884,0.894,0.9,0.893};
    // vector<double> fracVacinCluster = {0,0,0,0,0,0,0,0,0,0.006,0,0.022,0.084,0.138,0.162,0.168,0.124,0.169,0.178,0.191,0.208,0.175,0.179,0.181,0.164,0.284,0.397,0.505,0.588,0.624,0.652,0.629};

    double fcli;  // fraction clustered sias and vac
    double fclv;

    // interpIdx = lower_bound(pkaEnergyTable.begin(), pkaEnergyTable.end(), Epka) - pkaEnergyTable.begin();
    // if (interpIdx == 0)
    // {
    //     fcli = fracSIAinCluster[0];
    //     fclv = fracVacinCluster[0];
    // }
    // else if (interpIdx == pkaEnergyTable.size())
    // {
    //     fcli = fracSIAinCluster.back();
    //     fclv = fracVacinCluster.back();
    // }
    // else
    // {
    //     fcli = fracSIAinCluster[interpIdx-1] + (Epka - pkaEnergyTable[interpIdx-1]) / (pkaEnergyTable[interpIdx] - pkaEnergyTable[interpIdx-1]) * (fracSIAinCluster[interpIdx] - fracSIAinCluster[interpIdx-1]);
    //     fclv = fracVacinCluster[interpIdx-1] + (Epka - pkaEnergyTable[interpIdx-1]) / (pkaEnergyTable[interpIdx] - pkaEnergyTable[interpIdx-1]) * (fracVacinCluster[interpIdx] - fracVacinCluster[interpIdx-1]);
    // }

    fcli = 0.0185*pow(Epka, 0.326);
    fclv = 0.625 - 1.750e-4*TEMPERATURE;

    // Convert decimal number of frenkel pairs created into integer through sampling from Poisson distribution
    int n;
    int lower = int(numFP);
    int upper = lower + 1;
    randomNum = distribution(engine);
    if (lower + randomNum < numFP)
        n = upper;
    else
        n = lower;

    // Number of vac and interstitials that are clustered
    int ncli = Binomial(n, fcli);
    int nclv = Binomial(n, fclv);
    // int ncli = round(n * fcli);
    // int nclv = round(n * fclv);

    // Generate vacancy clusters and monovacancies by sampling cluster size distribution
    int nv = 0;

    for (nv = 0; nv < n-nclv; nv++)
    {
        int productAttr[3] = {-1, 0, 0}; // vac1
        int64 key = attrToKey(productAttr);
        addToObjectMap(key, count);
    }

    while (nv < n)
    {
        randomNum = distribution(engine);
        int clusterNum = 1 + (lower_bound(vacClusterSizeCDF.begin(), vacClusterSizeCDF.end(), randomNum) - vacClusterSizeCDF.begin());

        if (clusterNum == 1)
            continue;

        if (nv + clusterNum > n)
            clusterNum = n - nv;

        int productAttr[3] = {0};
        productAttr[0] = -clusterNum;
        int64 key = attrToKey(productAttr);
        addToObjectMap(key, count);

        nv += clusterNum;
    }

    // Generate SIA clusters and single sia
    int nsia = 0;

    for (nsia = 0; nsia < n-ncli; nsia++)
    {
        int productAttr[3] = {1, 0, 0}; // sia1
        int64 key = attrToKey(productAttr);
        addToObjectMap(key, count);
    }

    while (nsia < n)
    {
        randomNum = distribution(engine);
        int clusterNum = 1 + (lower_bound(siaClusterSizeCDF.begin(), siaClusterSizeCDF.end(), randomNum) - siaClusterSizeCDF.begin());

        if (clusterNum == 1)
            continue;

        if (nsia + clusterNum > n)
            clusterNum = n - nsia;

        int productAttr[3] = {0};
        productAttr[0] = clusterNum;
        int64 key = attrToKey(productAttr);
        addToObjectMap(key, count);

        nsia += clusterNum;
    }
}

// void SCDWrapper::getNeutronInsertion(const int count)
// {
//     // log10 of Epka(eV) values vs. cdf of RB19J sample holder, from SPECTER
//     // vector<double> log10Epkas = {-5.0000, -3.0000, -2.0000, -1.0000, 0.0000, 0.6990, 1.0000, 1.3010, 1.4771, 1.6021, 1.6990, 1.7782, 1.8451, 1.9031, 1.9542, 2.0000, 2.3010, 2.4771, 2.6021, 2.6990, 2.7782, 2.8451, 2.9031, 2.9542, 3.0000, 3.1761, 3.3010, 3.3979, 3.4771, 3.5441, 3.6021, 3.6532, 3.6990, 3.7404, 3.7782, 3.8129, 3.8451, 3.8751, 3.9031, 3.9294, 3.9542, 3.9777, 4.0000, 4.1761, 4.3010, 4.3979, 4.4771, 4.5441, 4.6021, 4.6532, 4.6990, 4.7404, 4.7782, 4.8129, 4.8451, 4.8751, 4.9031, 4.9294, 4.9542, 4.9777, 5.0000, 5.0792, 5.1461, 5.2041, 5.2553, 5.3010, 5.3424};
//     // vector<double> EpkaCDF = {0.000E+00, 4.487E-02, 4.928E-01, 7.778E-01, 8.218E-01, 8.468E-01, 8.719E-01, 8.836E-01, 8.916E-01, 8.983E-01, 9.041E-01, 9.081E-01, 9.114E-01, 9.140E-01, 9.163E-01, 9.303E-01, 9.367E-01, 9.409E-01, 9.440E-01, 9.467E-01, 9.491E-01, 9.511E-01, 9.529E-01, 9.544E-01, 9.607E-01, 9.658E-01, 9.701E-01, 9.734E-01, 9.761E-01, 9.782E-01, 9.800E-01, 9.816E-01, 9.829E-01, 9.841E-01, 9.851E-01, 9.860E-01, 9.868E-01, 9.876E-01, 9.882E-01, 9.889E-01, 9.894E-01, 9.899E-01, 9.938E-01, 9.960E-01, 9.975E-01, 9.984E-01, 9.990E-01, 9.994E-01, 9.996E-01, 9.997E-01, 9.998E-01, 9.999E-01, 9.999E-01, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00, 1.000E+00};

//     vector<double> log10Epkas = {-5.0000, -3.0000, -2.0000, -1.0000, 0.0000, 0.6990, 1.0000, 1.3010, 1.4771, 1.6021, 1.6990, 1.7782, 1.8451, 1.9031, 1.9542, 2.0000, 2.3010, 2.4771, 2.6021, 2.6990, 2.7782, 2.8451, 2.9031, 2.9542, 3.0000, 3.1761, 3.3010, 3.3979, 3.4771, 3.5441, 3.6021, 3.6532, 3.6990, 3.7404, 3.7782, 3.8129, 3.8451, 3.8751, 3.9031, 3.9294, 3.9542, 3.9777, 4.0000, 4.1761, 4.3010, 4.3979, 4.4771, 4.5441, 4.6021, 4.6532, 4.6990, 4.7404, 4.7782, 4.8129, 4.8451, 4.8751, 4.9031, 4.9294, 4.9542, 4.9777, 5.0000, 5.0792, 5.1461, 5.2041, 5.2553, 5.3010, 5.3424};
//     vector<double> EpkaCDF = {0.000000000E+00, 4.487262029E-02, 4.926943359E-01, 7.778088348E-01, 8.218171672E-01, 8.468413160E-01, 8.718453651E-01, 8.836338497E-01, 8.915471086E-01, 8.982785042E-01, 9.041154623E-01, 9.081313859E-01, 9.113171912E-01, 9.139804038E-01, 9.162064476E-01, 9.301958513E-01, 9.366760003E-01, 9.408507117E-01, 9.440003375E-01, 9.466886748E-01, 9.490785313E-01, 9.510312189E-01, 9.527849193E-01, 9.543577222E-01, 9.607142580E-01, 9.658215963E-01, 9.699953027E-01, 9.733409008E-01, 9.760161732E-01, 9.781909627E-01, 9.799707928E-01, 9.815335458E-01, 9.829093715E-01, 9.840580704E-01, 9.850871760E-01, 9.860221143E-01, 9.868329370E-01, 9.875701946E-01, 9.882468516E-01, 9.888702444E-01, 9.894376594E-01, 9.899579406E-01, 9.936944782E-01, 9.959506715E-01, 9.974239808E-01, 9.983597232E-01, 9.989551773E-01, 9.993204897E-01, 9.995423906E-01, 9.996847971E-01, 9.997787834E-01, 9.998402383E-01, 9.998843270E-01, 9.999157931E-01, 9.999364757E-01, 9.999517214E-01, 9.999630375E-01, 9.999717155E-01, 9.999781103E-01, 9.999829483E-01, 9.999935810E-01, 9.999972784E-01, 9.999987848E-01, 9.999994797E-01, 9.999998196E-01, 9.999999595E-01, 1.000000000E+00};


//     double randomNum = distribution(engine);
//     int interpIdx = lower_bound(EpkaCDF.begin(), EpkaCDF.end(), randomNum) - EpkaCDF.begin();
    
//     double log10Epka;
//     if (interpIdx == 0)
//     {
//         log10Epka = log10Epkas[0];
//     }
//     else // interpolate
//     {
//         log10Epka = log10Epkas[interpIdx-1] + (randomNum-EpkaCDF[interpIdx-1])/(EpkaCDF[interpIdx]-EpkaCDF[interpIdx-1]) * (log10Epkas[interpIdx]-log10Epkas[interpIdx-1]);
//     }

//     double Epka = pow(10., log10Epka); // eV

//     if (Epka < TDE)
//         return;

//     double numFP;
//     if (Epka < 43000)  // Setyawan 2015
//         numFP = 1.15e-2 * pow(Epka, 0.74);
//     else
//         numFP = 1.89e-5 * pow(Epka, 1.34);

//     double fcli = clamp(0.0185 * pow(Epka, 0.326), 0., 1.);     // fraction of generated interstitials that are clustered
//     double fclv = clamp(0.625 - 1.75e-4 * TEMPERATURE, 0., 1.); // fraction of generated vacancies that are clustered

//     // Convert decimal number of frenkel pairs created into integer through sampling from Poisson distribution
//     int n = Poisson(numFP);
//     int numClusteredVac = Binomial(n, fclv);
//     int numClusteredSia = Binomial(n, fcli);

//     // Vac and SIA cluster sizes (starting at size 2) from A.E. sand 2014
//     vector<double> vacClusterSizeCDF = {2.8528E-01, 3.4358E-01, 3.8455E-01, 4.1752E-01, 4.4561E-01, 4.7032E-01, 4.9251E-01, 5.1271E-01, 5.3131E-01, 5.4857E-01, 5.6469E-01, 5.7982E-01, 5.9409E-01, 6.0760E-01, 6.2042E-01, 6.3262E-01, 6.4426E-01, 6.5540E-01, 6.6606E-01, 6.7629E-01, 6.8612E-01, 6.9558E-01, 7.0469E-01, 7.1348E-01, 7.2195E-01, 7.3015E-01, 7.3807E-01, 7.4573E-01, 7.5315E-01, 7.6034E-01, 7.6731E-01, 7.7407E-01, 7.8063E-01, 7.8700E-01, 7.9319E-01, 7.9921E-01, 8.0505E-01, 8.1074E-01, 8.1627E-01, 8.2165E-01, 8.2688E-01, 8.3198E-01, 8.3695E-01, 8.4179E-01, 8.4650E-01, 8.5109E-01, 8.5556E-01, 8.5992E-01, 8.6417E-01, 8.6832E-01, 8.7236E-01, 8.7630E-01, 8.8014E-01, 8.8389E-01, 8.8755E-01, 8.9111E-01, 8.9459E-01, 8.9798E-01, 9.0129E-01, 9.0452E-01, 9.0767E-01, 9.1074E-01, 9.1374E-01, 9.1667E-01, 9.1952E-01, 9.2230E-01, 9.2501E-01, 9.2766E-01, 9.3024E-01, 9.3275E-01, 9.3520E-01, 9.3759E-01, 9.3992E-01, 9.4219E-01, 9.4440E-01, 9.4656E-01, 9.4866E-01, 9.5070E-01, 9.5269E-01, 9.5463E-01, 9.5652E-01, 9.5835E-01, 9.6013E-01, 9.6187E-01, 9.6356E-01, 9.6520E-01, 9.6679E-01, 9.6834E-01, 9.6984E-01, 9.7130E-01, 9.7271E-01, 9.7408E-01, 9.7541E-01, 9.7670E-01, 9.7794E-01, 9.7915E-01, 9.8032E-01, 9.8144E-01, 9.8253E-01, 9.8358E-01, 9.8459E-01, 9.8557E-01, 9.8651E-01, 9.8741E-01, 9.8828E-01, 9.8912E-01, 9.8992E-01, 9.9068E-01, 9.9141E-01, 9.9211E-01, 9.9278E-01, 9.9341E-01, 9.9402E-01, 9.9459E-01, 9.9513E-01, 9.9564E-01, 9.9612E-01, 9.9657E-01, 9.9699E-01, 9.9738E-01, 9.9774E-01, 9.9808E-01, 9.9838E-01, 9.9866E-01, 9.9891E-01, 9.9914E-01, 9.9934E-01, 9.9951E-01, 9.9965E-01, 9.9977E-01, 9.9987E-01, 9.9994E-01, 9.9998E-01, 1.0000E+00};
//     vector<double> siaClusterSizeCDF = {1.1313E-01, 1.4542E-01, 1.7255E-01, 1.9680E-01, 2.1904E-01, 2.3976E-01, 2.5925E-01, 2.7772E-01, 2.9532E-01, 3.1215E-01, 3.2829E-01, 3.4383E-01, 3.5882E-01, 3.7329E-01, 3.8731E-01, 4.0089E-01, 4.1407E-01, 4.2687E-01, 4.3932E-01, 4.5144E-01, 4.6325E-01, 4.7475E-01, 4.8598E-01, 4.9694E-01, 5.0764E-01, 5.1809E-01, 5.2831E-01, 5.3830E-01, 5.4808E-01, 5.5765E-01, 5.6701E-01, 5.7619E-01, 5.8517E-01, 5.9397E-01, 6.0260E-01, 6.1106E-01, 6.1935E-01, 6.2749E-01, 6.3547E-01, 6.4329E-01, 6.5097E-01, 6.5851E-01, 6.6591E-01, 6.7317E-01, 6.8030E-01, 6.8730E-01, 6.9417E-01, 7.0092E-01, 7.0755E-01, 7.1406E-01, 7.2046E-01, 7.2674E-01, 7.3292E-01, 7.3898E-01, 7.4494E-01, 7.5079E-01, 7.5655E-01, 7.6220E-01, 7.6775E-01, 7.7321E-01, 7.7857E-01, 7.8384E-01, 7.8902E-01, 7.9411E-01, 7.9911E-01, 8.0403E-01, 8.0886E-01, 8.1360E-01, 8.1826E-01, 8.2284E-01, 8.2734E-01, 8.3176E-01, 8.3611E-01, 8.4038E-01, 8.4457E-01, 8.4869E-01, 8.5273E-01, 8.5670E-01, 8.6060E-01, 8.6443E-01, 8.6819E-01, 8.7189E-01, 8.7551E-01, 8.7907E-01, 8.8256E-01, 8.8599E-01, 8.8936E-01, 8.9266E-01, 8.9589E-01, 8.9907E-01, 9.0218E-01, 9.0524E-01, 9.0823E-01, 9.1117E-01, 9.1405E-01, 9.1687E-01, 9.1963E-01, 9.2234E-01, 9.2499E-01, 9.2759E-01, 9.3013E-01, 9.3262E-01, 9.3505E-01, 9.3743E-01, 9.3976E-01, 9.4204E-01, 9.4427E-01, 9.4644E-01, 9.4857E-01, 9.5065E-01, 9.5267E-01, 9.5465E-01, 9.5658E-01, 9.5846E-01, 9.6030E-01, 9.6209E-01, 9.6383E-01, 9.6552E-01, 9.6717E-01, 9.6878E-01, 9.7034E-01, 9.7186E-01, 9.7333E-01, 9.7476E-01, 9.7614E-01, 9.7748E-01, 9.7878E-01, 9.8004E-01, 9.8126E-01, 9.8243E-01, 9.8357E-01, 9.8466E-01, 9.8571E-01, 9.8673E-01, 9.8770E-01, 9.8863E-01, 9.8953E-01, 9.9038E-01, 9.9120E-01, 9.9198E-01, 9.9272E-01, 9.9342E-01, 9.9409E-01, 9.9472E-01, 9.9531E-01, 9.9586E-01, 9.9638E-01, 9.9687E-01, 9.9731E-01, 9.9773E-01, 9.9810E-01, 9.9845E-01, 9.9875E-01, 9.9903E-01, 9.9927E-01, 9.9947E-01, 9.9964E-01, 9.9978E-01, 9.9988E-01, 9.9996E-01, 9.9999E-01, 1.0000E+00};

//     // Generate vacancy clusters and monovacancies by sampling cluster size distribution
//     int nv = 0;
//     while (nv < numClusteredVac)
//     {
//         randomNum = distribution(engine);
//         int clusterNum = 2 + (lower_bound(vacClusterSizeCDF.begin(), vacClusterSizeCDF.end(), randomNum) - vacClusterSizeCDF.begin());

//         if (nv + clusterNum > numClusteredVac)
//             clusterNum = numClusteredVac - nv;

//         int productAttr[3] = {0};
//         productAttr[0] = -clusterNum;
//         int64 key = attrToKey(productAttr);
//         addToObjectMap(key, count);

//         nv += clusterNum;
//     }

//     int64 vacKey = -1000000;
//     addToObjectMap(vacKey, count, n - numClusteredVac);

//     // Generate SIA clusters and single sia
//     int nsia = 0;
//     while (nsia < numClusteredSia)
//     {
//         randomNum = distribution(engine);
//         int clusterNum = 2 + (lower_bound(siaClusterSizeCDF.begin(), siaClusterSizeCDF.end(), randomNum) - siaClusterSizeCDF.begin());

//         if (nsia + clusterNum > numClusteredSia)
//             clusterNum = numClusteredSia - nsia;

//         int productAttr[3] = {0};
//         productAttr[0] = clusterNum;
//         int64 key = attrToKey(productAttr);
//         addToObjectMap(key, count);

//         nsia += clusterNum;
//     }

//     int64 siaKey = 1000000;
//     addToObjectMap(siaKey, count, n - numClusteredSia);
// }

void SCDWrapper::getIonInsertion(const int n, const double dt, fstream& fs)
{
    ++reactions[5][n];
    double ionEnergy = 0.0;
    double totalEnergy = (double)Poisson(AVG_ION_EN[n]);
    int i, j, ndef = 0;
    int64 clusterKey;
    //fstream fk;
    //fk.open("V_insertion.txt", ios::app);
    doseIon[n] += damage.getDpaRate(n) / damage.getTotalIonRate();
    totalDpa += doseIon[n];
    in_time += dt;
    //fk << damage.getDpaRate(n) << " * " << dt << "  " << in_time << "   ";
    CascadeDamage damage;
    while (ionEnergy < totalEnergy) 
    {
        double pkaEnergy = 0.0;
        while (pkaEnergy < 620.0) 
        {
            // Limit to produce a stable Frenkel pair in keV (from Troev et al (2011)).
            double xi = ((double)rand() / RAND_MAX) * cpdf.getMaxPossibility(n);
            ionEnergy += pkaEnergy;
            pkaEnergy = cpdf.samplePkaEnergy(xi, n);
        }
        damage.generateIonDamage(pkaEnergy, ndef);
        for (i = 0;i < damage.size();++i) 
        {
            int sign = (i == 0) ? -1 : 1;
            for (j = 0; j < ndef; ++j) 
            {
                const int number = damage.getDamage(i, j);
                if (number != 0) 
                {
                    clusterKey = (int64)sign*(j + 1)*(pow(10.0, (double)EXP10*(LEVELS - 1)));
                    /*
                    if(clusterKey > 1000){
                        if (allObjects.find(clusterKey) == allObjects.end()) {
                            Object* newObject = new Object(clusterKey, n, number);
                            writeSInkFile(newObject, n, 0.0);
                        } // we don't have this cluster
                        else {
                            Object* tempObject = allObjects[clusterKey];
                            writeSInkFile(newObject, n, 0.0);
                        }// we have this cluster
                        continue;
                    }//park SIAs
                    */
                    /* generate cluster key */
                    addToObjectMap(clusterKey, n, number);
                    if (LOG_REACTIONS)
                        fs <<"ion insertion in element "<< n <<", gain " << number <<" " << clusterKey <<endl;
                }
            }
        }
        damage.cleanDamage();
        ionEnergy += pkaEnergy;
    }
    //fk<<dt<<"   "<<nov<<endl;
    //fk.close();
}

void SCDWrapper::getParticleInsertion(const int n, const double dt, fstream& fs)
{
#ifdef ELECTRON
    getElectronInsertion(n);
#elif defined (NEUTRON)
    getNeutronInsertion(n);
#elif defined (ION)
    getIonInsertion(n, dt, fs);
#endif
}

void SCDWrapper::getHeInsertion(const int n)
{
    int channel = 1;
    int64 clusterKey = atomProperty(INTERSTITIAL, channel + 1);
    addToObjectMap(clusterKey, n);
}

void SCDWrapper::getHInsertion(const int n, const double dt, fstream& fs)
{
    ++reactions[6][n];
    int channel = 2;
    int64 clusterKey = atomProperty(INTERSTITIAL, channel + 1);
    addToObjectMap(clusterKey, n);
    if (LOG_REACTIONS)
        fs << "H insertion: get 1 " << clusterKey <<" in element "<< n <<endl;
}

void restart(long int & iStep, double & advTime, SCDWrapper *srscd)
{
    int64 objectKey;
    int numberSinks[NUM_SINKS*(LEVELS+1)] = { 0 };  // LEVELS + 1 to separate vacancies and sia, *2 because separate dislocations and grain boundaries
    int number[POINTS] = { 0 };
    int step = 0;
    string skip;
    string oneLine;
    stringstream lineHold;
    /* update species informtion */
    ifstream ofile("restart.txt");
    if (ofile.good()) {
        while (getline(ofile, oneLine)) {
            step++;
            lineHold.str(oneLine);
            if (step == 1) {
                /* update iStep */
                /* line 1 */
                lineHold >> skip >> skip >> iStep;
            }
            else if (step == 2){
                /* update advTime */
                /* line 2 */
                lineHold >> skip >> skip >> advTime;
            }
            else if(step == 3){
                /* line 3 */
                lineHold >> skip >> skip >> fluenceH;
            }else{
                lineHold >> skip >> objectKey;
                for (int i = 0; i < POINTS; i++) {
                    lineHold >> number[i];
                    srscd->addToObjectMap(objectKey, i, number[i]);
                }
            }
            lineHold.clear();
        }
    }
    ofile.close();

    /* update sink numbers */
    ifstream file("sink.txt");
    if (file.good()) {
        for(int j = 0; j < POINTS; j++){
            if (getline(file, oneLine)) {
                lineHold.str(oneLine);
                for (int i = 0; i < NUM_SINKS*(LEVELS+1); i++) {
                    lineHold >> numberSinks[i];
                }
                srscd->updateSinks(j,numberSinks);
                lineHold.clear();
            }
        }
    }
    srscd->writeSinkFile(0, 0, 0);
    file.close();

    srscd->clearBoundaryChangeQs();
}

int SCDWrapper::countDefectNumber(const int count, string type){
    fstream vd, id, hd;
    fstream v1d, v2d, v3d;
    /**
     * vd : vacancy number vs depth file\
     * id : SIA number vs depth file\
     * hd : hydrogen number vs depth file\
     **/
    if(type == "V"){
        vd.open("vd.txt", ios::out);
        /*
        v1d.open("v1d.txt", ios::out);
        v2d.open("v2d.txt", ios::out);
        v3d.open("v3d.txt", ios::out);
        */
        fv.open("vt.txt", ios::app);
        
    }else if(type == "SIA"){
        //id.open("id.txt", ios::out);
    }else{
        hd.open("hd.txt", ios::out);
    }
    double v = 0.0;
    /* this file work after damage*/
    
    int ndef[POINTS] = {0}; /* number of this object in every element*/
    int tndef = 0; /*number of this kind of defect in total */
    unordered_map<int64, Object*>::iterator iter;
    for(int i=0; i<POINTS; i++){
        double volume = volumeAtIndex(i);
        
        for (iter = allObjects.begin(); iter != allObjects.end(); ++iter){
            Object* thisObject = iter -> second;
            int totalNumber  = thisObject -> getNumber(i);
            int attribute = thisObject -> getAttri(count);
            if(count == 0){
                if(attribute>0 && type=="SIA"){
                    ndef[i] += totalNumber * abs(attribute);
                }
                else if (attribute < 0 && type == "V")
                {
                    ndef[i] += totalNumber * abs(attribute);
                }
            }else{
                ndef[i] += totalNumber * abs(attribute);
            }
        }
        tndef += ndef[i];
        if(type == "V"){
            vd << 10+20*i << "      "<<ndef[i]/volume<<endl;
            v += ndef[i]/volume;
         }else if(type == "SIA"){
            //id<< 18+36*(i-1) << "     "<<ndef[i]/volume<<endl;
         }else{
            hd<< 10+20*(i) << "     " << ndef[i]/volume<<endl;
        }
    }
    
    if(type == "V"){
        fv << v << endl;
        vd.close();
        //v1d.close();
        //v2d.close();
        //v3d.close();
        fv.close();
        
    }else if(type == "SIA"){
        //id.close();
    }else{
        hd.close();
    }
    return tndef;
}

void SCDWrapper::sizeDistribution(){
    unordered_map<int, int> sizeD;
    /* preset size from -6 to 6 */
    for(int preSize = -6; preSize<= 6; preSize++){
        std::pair<int, int> oneSize(preSize, 0);
        sizeD.insert(oneSize);
    }
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter){
        Object* thisObject = iter -> second;
        int totalNumber  = thisObject -> getTotalNumber();
        int attribute = thisObject -> getAttri(0);
        if(sizeD.find(attribute) != sizeD.end()){
           /* find this size */
            sizeD[attribute] += totalNumber;
        }else{
            std::pair<int, int> oneSize(attribute, totalNumber);
            sizeD.insert(oneSize);
        }
    }
    /* write to a file */
    fstream sd;
    sd.open("sd.txt", ios::out);
    unordered_map<int, int>::iterator iterr;
    for (iterr = sizeD.begin(); iterr != sizeD.end(); ++iterr){
        int n = iterr -> first;
        if(n > 0){
            //SIA cluster
            sd << iterr->first << "  " <<pow(ATOMICVOLUME*n/PI/BURGER, 0.5)<<"  "<< iterr->second/VOLUME<<endl;
            
        }else{
            //V cluster
            sd << iterr->first << "  " <<pow(3*ATOMICVOLUME*abs(n)/4/PI/BURGER, 0.333)<<"   "<< iterr->second/VOLUME<<endl;
        }
    }
    sd.close();
}

void SCDWrapper::writeReaction(){
    fstream fr;
    fr.open("reaction.txt", ios::out);
    for(int i=0; i<8; i++){
        for(int j=0; j<POINTS; j++){
            fr << i <<"    "<< j <<"    "<<reactions[i][j]<<endl;
        }
    }
    fr.close();
}


/* test functions */
void SCDWrapper::countRatioDistribution(double& t){
    // double sW = DENSITY * (VOLUME / 20) * SURFACE_THICKNESS;
    /*number of surface tungsten */
    int sH = 0; /* number of surface hydrogen*/
    unordered_map<int64, Object*>::iterator iter;
    fstream fo;
    fstream vha; // vacancy hydrogen number of all objects
    fstream vhc; // vacancy hydrogen number of only (V-H)clusters
    //fo.open("object.txt", ios::out);
    //vha.open("vha.txt", ios::app);
    //vhc.open("vhc.txt", ios::app);
    double ah = 0.0, av = 0.0, ch = 0.0, cv = 0.0;
    // all H, all V, cluster H, cluster V
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter){
        Object* thisObject = iter -> second;
        int attrTwo = thisObject->getAttri(2);
        int attrZero = thisObject->getAttri(0);
        int totalNumber = thisObject->getTotalNumber();
        
        if(attrZero < 0){
            av += abs(attrZero) * totalNumber;
        }
        ah += attrTwo * totalNumber;
        if(attrZero < 0 && attrTwo >0){
            /* if have vacancy */
            /* a v-h cluster */
            cv += abs(attrZero) * totalNumber;
            ch += attrTwo * totalNumber;
        }
        int surfaceNumber = thisObject->getNumber(0); /* get the number on surface */
        sH += attrTwo * surfaceNumber;  /* get total number of H on surface */
        //fo << attrZero << " " << attrTwo << "   "<< totalNumber << endl;
    }
    /*
    if(av == 0.0){
        vha << t << "   -1" << endl;
    }else{
        vha << t<< "    "<< ah/av << endl;
    }
    if(cv == 0.0){
        vhc << t << "   -1" << endl;
    }else{
        vhc << t << "   " << ch/cv <<endl;
    }
    */
    //fo.close();
    //vha.close();
    //vhc.close();
   
}

void SCDWrapper::test(const int inputV){
    ofstream fo;
    fo.open("SAV.txt", ios::app);
    if(start){
        fo<<"DoseH"<<"  "<<"InputVacancy"<<"    "<<"NowV"<<"    "<<"AnnilV"<<"  "<<"CreateV"<<" "<<"NowH"<<"    "<<"V/H"<<endl;
        start = 0;
    }
    int nowV = countDefectNumber(0, "V");
    int hydrogen = countDefectNumber(2, "H");
    fo << fluenceH << "    " << inputV/VOLUME << " " << nowV/VOLUME << " " << sinkV/VOLUME << " " << annilV/VOLUME << "  " << (nowV+sinkV+annilV-inputV)/VOLUME << "   " << hydrogen/VOLUME << " ";
    if(hydrogen == 0){
        fo << "---" << endl;
    }else{
        fo << double(hydrogen)/double(nowV) << endl;
    }
    fo.close();
}

void SCDWrapper::drawSpeciesAndReactions( double&t ){
    /* draw species */
    countRatioDistribution(t);
    char buffer[20];
    char v[20];
    char* __attribute__((unused)) result = gcvt(t, 8, buffer);
    result = gcvt(VOLUME, 1, v);
    std::ostringstream cmdstr;
    std::ostringstream cmdstr1;
    //gs.cmd("set key bmargin left horizontal Right noreverse enhanced title \"Td = 573K \" box\n");
    //gs.cmd("set grid \n");
    cmdstr<<"plot \"object.txt\" using 1:2:3 with labels title \" t = "<<buffer<<" s,   V= " << v << " cm^3 \" font \"New-Roman, 20, Bold \" textcolor lt 7 \n";
    //gs.cmd(cmdstr.str());
    
    /* draw reactions */
    //writeReaction();
    //gr.cmd("set key bmargin left horizontal Right noreverse enhanced title \"Td = 300K \" box\n");
    //gr.cmd("set grid \n");
    // gr.cmd("set view map\n");
    //gr.cmd("set cbrange[0:500]\n");
    //gr.cmd("set cbtics add ('>1000' 1000) \n");
    //gr.cmd("set pm3d map\n");
    //gr.cmd("set palette model CMY rgbformulae 7, 5, 15\n");
    //gr.cmd("set palette defined(0 \"green\", 1 \"blue\", 2 \"red\" , 3 \"orange\", 4 \"yellow\", 5 \"grey\")\n");
    //cmdstr1<<"splot \"reaction.txt\" using 2:1:3 with points palette ps 2 pt 3 title \" t = "<<buffer<<" s,   V= " << v << " cm^3 \" \n";
    //gr.cmd("set pm3d map\n");
    //gr.cmd(cmdstr1.str());
    /*
    if(plotTime1 == 1){
        cout<<"press a char"<<endl;
        getchar();
        plotTime1++;
    }
    */
    //sleep(0.02);
    
}

void SCDWrapper::drawDamage(double& t){
    /* draw size distribution */
    char buffer[20];
    char* __attribute__((unused)) result = gcvt(t,8,buffer);
    double dpa = getTotalDpa();
    std::ostringstream cmdstr;
    std::ostringstream cmdstr1;
    std::ostringstream cmdstr2;
    sizeDistribution();
    cmdstr2 <<"set key bmargin left horizontal Right noreverse enhanced title \"T = "<< TEMPERATURE <<"K \" box\n";
    //gd1.cmd(cmdstr2.str());
    cmdstr<<"plot \"sd.txt\" using 1:2 with histeps title \" t = "<<buffer<<" \" lw 3\n";
    //gd1.cmd(cmdstr.str());
    /* draw concentration-depth */
    fv.open("vt.txt", ios::app);
    fv << dpa << "    ";
    fv.close();
    countDefectNumber(0, "V");
    countDefectNumber(0, "SIA");
    //gd2.cmd(cmdstr2.str());
    //cmdstr1<<"plot \"vd.txt\" using 1:2 with lp title \" t = "<<buffer<<" V \" lw 3\n";
    //cmdstr2<<"replot \"id.txt\" using 1:2 with lp title \" t = "<<buffer<<" SIA \" lw 3\n";
    //gd2.cmd(cmdstr1.str());
    //gd2.cmd(cmdstr2.str());
    //cmdstr1<<"plot \"vd.txt\" using 1:2 with lp title \" t = "<<buffer<<"       V \" lw 3, \"id.txt\" using 1:2 with lp title \" SIA \" lw 3\n";
    cmdstr1<<"plot \"vd.txt\" using 1:2 with lp title \" dpa = " << dpa << "       V \" lw 3, \"id.txt\" using 1:2 with lp title \" SIA \" lw 3\n";
    //gd2.cmd(cmdstr1.str());
    //sleep(0.02);
    //gv.cmd("plot \"vt.txt\" using 1:2 w lp lw 3 title \"V\" \n");
    /* wait to adjust the screen locations */
    /*
    if(plotTime == 1){
        cout<<"press a char"<<endl;
        getchar();
        plotTime++;
    }
    */
}

void SCDWrapper::drawHD(double& t){
    /* draw concentration-depth */
    char buffer[20];
    char* __attribute__((unused)) result = gcvt(t,8,buffer);
    std::ostringstream cmdstr1;
    countDefectNumber(0, "V");
    countDefectNumber(0, "SIA");
    countDefectNumber(2, "H");
    // system("python3 plot_hd.py");
    // gh1.cmd("set key bmargin left horizontal Right noreverse enhanced title \"Td= 573K \" box\n");
    // cmdstr1<<"plot \"vd.txt\" using 1:2 with lp title \" t = "<<buffer<<" V \" lw 3, \"id.txt\" using 1:2 with lp title \"SIA\" lw 3, \"hd.txt\" using 1:2 with lp title \"H\" lw 3\n";
    // gh1.cmd(cmdstr1.str());
    
   
    /* draw H to V ratio */
    //gh2.cmd("plot \"vha.txt\" using 1:2 w l title \" whole H/V \" lw 3, \"vhc.txt\" using 1:2 w l title \" cluster H/V \" lw 3 \n");
    
    // if(plotTime2==1){
    //     getchar();
    //     ++plotTime2;
    // }
    //sleep(0.02);
}

void SCDWrapper::writeVacancy(){
    double dpa = getTotalDpa();
    fstream fv;
    fv.open("V.txt", ios::app);
    fv<<"At dpa = " <<dpa<<", "<<generationNumber<<" vacancies are generated, "<<annilV<<" consumed from combination, "<<sinkV<<" goes to sink, "<< generationNumber - annilV - sinkV<<" left in bulk"<<endl;
    fv.close();
}

double SCDWrapper::getTotalDpa(){
    return totalDpa;
}

double SCDWrapper::getHSaturationConcentration() const
{
    return DENSITY * exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE);
    double concentration = H_SATURATION_CONCENTRATION;
    bool dimer = false;
    bool vhpair = false;

    for (unordered_map<int64, Object*>::const_iterator iter = allObjects.begin(); iter != allObjects.end() && !dimer && !vhpair; iter ++)
    {
        Object* tempObj = iter->second;
        // if HH objects is present

        if (tempObj->getKey() == 2 && !dimer)
        {
            concentration += DENSITY * 8 * exp(-(2 * HEAT_OF_SOLUTION - HH_BIND_E) / KB / TEMPERATURE);
            dimer = true;
        }
        // if mV-nH objects are present
        if (tempObj->getAttri(0) <= -1 && tempObj->getAttri(2) >= 1 && !vhpair)
        {
            concentration += DENSITY * 8 * exp(-(HEAT_OF_SOLUTION + V_FORM_E - VH_BIND_E) / KB / TEMPERATURE);
            vhpair = true;
        }
    }
    return concentration;
    // 8 is the coordination number in a BCC lattice
}

void SCDWrapper::setDomain(int start, int end)
{
    startIndex = start;
    endIndex = end;

    if (UNIFORM_FREE_H_ON)
    {
        // To simulate a constant concentration of free 1H in each element,
        // put a placeholder of 1 count of 1H in each mesh element as a placeholder
        // to make it easier to calculate rates. There is not actually exactly 1H in each element.
        int64 HKey = 1;
        Object* HObj;
        if (allObjects.find(HKey) != allObjects.end())
        {
            HObj = allObjects[HKey];
        }
        else
        {
            HObj = new Object(HKey, 0, 0);
            addNewObjectToMap(HObj);
        }

        for (int n = startIndex; n <= endIndex; n++)
        {
            if (n == SURFACE_INDEX 
                || n == SUBSURFACE_INDEX
                || (n == BACK_SUBSURFACE_INDEX && BACK_DESORB)
                || (n == BACK_SURFACE_INDEX && BACK_DESORB))
            {
                // 0 H in these mesh elements are needed bc no relevant reactions in uniform free H mode would take place here
                HObj->addNumber(n, -HObj->getNumber(n));
                updateObjectInMap(HObj, n);
                objectsInElement[n].erase(HKey);
            }
            else
            {
                // place 1 count of 1H in mesh element
                HObj->addNumber(n, 1 - HObj->getNumber(n));
                updateObjectInMap(HObj, n);
                objectsInElement[n][HKey] = HObj;
            }
        }
    }
}

void SCDWrapper::fillNoneReaction(long double maxDomainRate)
{
    /*
     * Fill the remainder of the domain rate with NONE reaction
     * so that each processor can move at the same time step in parallel
     * (Dunn 2016)
     */
    noneRate = maxDomainRate - domainRate;
}

void SCDWrapper::clearNoneReaction()
{
    noneRate = 0.0;
}

vector<BoundaryChange>* SCDWrapper::getLeftBoundaryChangeQ()
{
    return &leftBoundaryChangeQ;
}

vector<BoundaryChange>* SCDWrapper::getRightBoundaryChangeQ()
{
    return &rightBoundaryChangeQ;
}

void SCDWrapper::clearBoundaryChangeQs()
{
    leftBoundaryChangeQ.clear();
    rightBoundaryChangeQ.clear();
}

void SCDWrapper::implementBoundaryChanges(vector<BoundaryChange>& boundaryChanges)
{
    bool updateFront = false, updateBack = false;
    for (BoundaryChange& bc: boundaryChanges)
    {
        if (allObjects.find(bc.objKey) != allObjects.end()) {
            /* object found! then number of instances in this element increase */
            Object* anObject = allObjects[bc.objKey];
            anObject->addNumber(bc.pointIndex, bc.change); // negative change accounts for decrease
            // cout << "object: " << bc->objKey << " add " << bc->change << " at " << bc->pointIndex << endl;
            updateObjectInMap(anObject, bc.pointIndex);
        }
        else {
            /* object didn't find! build new object and insert it into map */
            Object* anObject = new Object(bc.objKey, bc.pointIndex, bc.change);
            addNewObjectToMap(anObject);
        }

        if (bc.pointIndex == startIndex || bc.pointIndex == startIndex - 1)
            updateFront = true;
        else
            updateBack = true;
    }

    if (updateFront)
        updateMatrixRate(startIndex, Reaction::NONE);
    if (updateBack)
        updateMatrixRate(endIndex, Reaction::NONE);
    computeDomainRate();
}

int SCDWrapper::getStartIndex()
{
    return startIndex;
}

int SCDWrapper::getEndIndex()
{
    return endIndex;
}

vector<BoundaryChange> SCDWrapper::getSpatialElement(int n)
{
    /* Return all of the object counts (not including sinks) at this spatial element */
    vector<BoundaryChange> objects;
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter)
    {
        Object* obj = iter->second;
        if (obj->getNumber(n) > 0)
        {
            objects.push_back(BoundaryChange(obj->getKey(), n, obj->getNumber(n)));
        }
    }
    return objects;
}

void SCDWrapper::getSink(int n, int* output)
{
    for (int type = 0; type < NUM_SINKS; type++)
    {
        for (int level = 0; level < LEVELS+1; level++)  // levels + 1 because vacancy and sia each have their own levels, so vac + sia + helium + H = 4 levels = 3 + 1
        {
            if (type == 0)
                output[type*(LEVELS+1)+level] = sinksDislocationScrew[n][level];
            else if (type == 1)
                output[type*(LEVELS+1)+level] = sinksDislocationEdge[n][level];
            else
                output[type*(LEVELS+1)+level] = sinksGrainBndry[n][level];
        }
    }
}

void SCDWrapper::addSpatialElement(int newGhostIndex, vector<BoundaryChange> newGhostObjects, int newBoundaryIndex, int* newBoundarySinks)
{
    /* Clear out spatial element to use as our new ghost index */
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter)
    {
        Object* obj = iter->second;
        if (obj->getNumber(newGhostIndex) > 0)
        {
            reduceFromObjectMap(obj->getKey(), newGhostIndex, obj->getNumber(newGhostIndex));
        } 
    }

    /* Put in new objects into our new ghost index */
    for (BoundaryChange& bc: newGhostObjects)
    {
        addToObjectMap(bc.objKey, bc.pointIndex, bc.change);
    }

    /* Put in new sink counts into our new boundary index */
    for (int type = 0; type < NUM_SINKS; type++)
    {
        for (int level = 0; level < LEVELS+1; level++)
        {
            if (type == 0)
                sinksDislocationScrew[newBoundaryIndex][level] = newBoundarySinks[type*(LEVELS+1)+level];
            else if (type == 1)
                sinksDislocationEdge[newBoundaryIndex][level] = newBoundarySinks[type*(LEVELS+1)+level];
            else
                sinksGrainBndry[newBoundaryIndex][level] = newBoundarySinks[type*(LEVELS+1)+level];
            if (level == 0)
                computeSinkDissRate(0, newBoundaryIndex);  // V Diss rate
            else if (level == 3)
                computeSinkDissRate(1, newBoundaryIndex);  // H Diss rate
        }
    }

    clearBoundaryChangeQs();  // from the add/reduceToObjectMap calls
}

void SCDWrapper::recalculateAllRates()
{
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter)
    {
        Object* object = iter->second;
        object->computeThermalProperties();
    }

    for (int point = 0; point < POINTS; point++)
    {
        for (iter = allObjects.begin(); iter != allObjects.end(); ++iter)
        {
            Object* object = iter->second;
            if (object->getNumber(point) > 0)
                updateObjectInMap(object, point);
        }

        // Omit damage recalculation for now, assume it doesn't change with temperature

        computeSinkDissRate(0, point);
        computeSinkDissRate(1, point);
    }

    examineDomainRate();
}

void SCDWrapper::writeDesorbedFile(double time)
{
    desorbedFile << time << " " << numHDesorbed << endl;
}
