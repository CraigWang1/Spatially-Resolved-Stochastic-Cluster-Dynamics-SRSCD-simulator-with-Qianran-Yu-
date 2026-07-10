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
    if (pointIndex > endIndex)
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

    double ebHDislocationScrew = 0.55, ebHDislocationEdge = 0.89, ebHGrainBndry = 0.85; //binding and migration energy of hydrogen

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
    updateObjectInMap(anObject, n);

    // If this object exists in this spatial element, account for it
    if (anObject->getNumber(n) > 0)
    {
        objectsInElement[n][key] = anObject;
    }
    else
    {
        objectsInElement[n].erase(key);
    }

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
    vector<double> log10Epkas = {-5.0000, -3.0000, -2.0000, -1.0000, 0.0000, 0.6990, 1.0000, 1.3010, 1.4771, 1.6021, 1.6990, 1.7782, 1.8451, 1.9031, 1.9542, 2.0000, 2.3010, 2.4771, 2.6021, 2.6990, 2.7782, 2.8451, 2.9031, 2.9542, 3.0000, 3.1761, 3.3010, 3.3979, 3.4771, 3.5441, 3.6021, 3.6532, 3.6990, 3.7404, 3.7782, 3.8129, 3.8451, 3.8751, 3.9031, 3.9294, 3.9542, 3.9777, 4.0000, 4.1761, 4.3010, 4.3979, 4.4771, 4.5441, 4.6021, 4.6532, 4.6990, 4.7404, 4.7782, 4.8129, 4.8451, 4.8751, 4.9031, 4.9294, 4.9542, 4.9777, 5.0000, 5.0792, 5.1461, 5.2041, 5.2553, 5.3010, 5.3424};
    vector<double> EpkaCDF = {0.000000000E+00, 4.487262029E-02, 4.926943359E-01, 7.778088348E-01, 8.218171672E-01, 8.468413160E-01, 8.718453651E-01, 8.836338497E-01, 8.915471086E-01, 8.982785042E-01, 9.041154623E-01, 9.081313859E-01, 9.113171912E-01, 9.139804038E-01, 9.162064476E-01, 9.301958513E-01, 9.366760003E-01, 9.408507117E-01, 9.440003375E-01, 9.466886748E-01, 9.490785313E-01, 9.510312189E-01, 9.527849193E-01, 9.543577222E-01, 9.607142580E-01, 9.658215963E-01, 9.699953027E-01, 9.733409008E-01, 9.760161732E-01, 9.781909627E-01, 9.799707928E-01, 9.815335458E-01, 9.829093715E-01, 9.840580704E-01, 9.850871760E-01, 9.860221143E-01, 9.868329370E-01, 9.875701946E-01, 9.882468516E-01, 9.888702444E-01, 9.894376594E-01, 9.899579406E-01, 9.936944782E-01, 9.959506715E-01, 9.974239808E-01, 9.983597232E-01, 9.989551773E-01, 9.993204897E-01, 9.995423906E-01, 9.996847971E-01, 9.997787834E-01, 9.998402383E-01, 9.998843270E-01, 9.999157931E-01, 9.999364757E-01, 9.999517214E-01, 9.999630375E-01, 9.999717155E-01, 9.999781103E-01, 9.999829483E-01, 9.999935810E-01, 9.999972784E-01, 9.999987848E-01, 9.999994797E-01, 9.999998196E-01, 9.999999595E-01, 1.000000000E+00};

    double randomNum = distribution(engine);
    size_t interpIdx = lower_bound(EpkaCDF.begin(), EpkaCDF.end(), randomNum) - EpkaCDF.begin();
    if (interpIdx >= EpkaCDF.size())
    {
        interpIdx = EpkaCDF.size() - 1;
    }

    double log10Epka;
    if (interpIdx == 0)
    {
        log10Epka = log10Epkas[0];
    }
    else // interpolate
    {
        log10Epka = log10Epkas[interpIdx-1] + (randomNum-EpkaCDF[interpIdx-1])/(EpkaCDF[interpIdx]-EpkaCDF[interpIdx-1]) * (log10Epkas[interpIdx]-log10Epkas[interpIdx-1]);
    }

    double Epka = pow(10., log10Epka); // eV

    double numFP;
    if (Epka <= 48000)  // Sicong He 2025
        numFP = 3.81 * pow(Epka/1000.0, 0.62);
    else
        numFP = 0.50 * pow(Epka/1000.0, 1.15);

    // Convert decimal number of frenkel pairs created into integer through sampling from Poisson distribution
    int n;
    int lower = int(numFP);
    int upper = lower + 1;
    randomNum = distribution(engine);
    if (lower + randomNum < numFP)
        n = upper;
    else
        n = lower;

    // Generate vacancy clusters and monovacancies by sampling cluster size distribution
    int nv = 0;
    while (nv < n)
    {
        randomNum = distribution(engine);
        int clusterNum = 1 + (lower_bound(vacClusterSizeCDF.begin(), vacClusterSizeCDF.end(), randomNum) - vacClusterSizeCDF.begin());

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
    while (nsia < n)
    {
        randomNum = distribution(engine);
        int clusterNum = 1 + (lower_bound(siaClusterSizeCDF.begin(), siaClusterSizeCDF.end(), randomNum) - siaClusterSizeCDF.begin());

        if (nsia + clusterNum > n)
            clusterNum = n - nsia;

        int productAttr[3] = {0};
        productAttr[0] = clusterNum;
        int64 key = attrToKey(productAttr);
        addToObjectMap(key, count);

        nsia += clusterNum;
    }
}

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
