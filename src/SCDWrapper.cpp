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
SCDWrapper::SCDWrapper():allObjects(), damage(allObjects), cpdf(), startIndex(0), endIndex(0), totalDpa(0)
{
    formationE[1] = V_FORM_E; 

    for (int i = 0; i < POINTS; ++i) {
        computeMatrixRate(i);
    } /* initialized matrix rate in every element */
    computeBulkRate();  /* initialized total rate in the bulk */
    
    /* initialize reactions[][] */
    for(int i=0; i<8; i++){
        for(int j=0; j<POINTS; j++){
            reactions[i][j] = 0;
        }
    }
    /* initialize sink numbers */
    setSinks();

    lastElemSaturated = false;
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
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter) {
        int64 tempKey = iter->first;
        Object* tempObject = iter->second;
        int totalNumber = tempObject->getTotalNumber();
        while(totalNumber == 0) {
            iter++;
            removeObjectFromMap(tempKey);
            if (iter != allObjects.end()) {
                tempObject = iter->second;
                tempKey = iter->first;
                totalNumber = tempObject->getTotalNumber();
            }
            else {
                break;
            }
        }
        if (iter == allObjects.end()) {
            break;
        }
        Bundle* tempBundle = linePool[tempObject];
        OneLine* tempLine = tempBundle->lines[n];
        if (tempLine != nullptr) {
            matrixRate[n] += tempLine->computeTotalRate();
            //tempLine->display(tempObject);/* Qianran 0925 */
        }
    }
    matrixRate[n] += damage.getTotalDamage(n);
    matrixRate[n] += sinkDissRate[0][n];
    matrixRate[n] += sinkDissRate[1][n];
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
        matrixRate[i] = 0;
        matrixRate[i] += damage.getTotalDamage(i);
        matrixRate[i] += sinkDissRate[0][i];
        matrixRate[i] += sinkDissRate[1][i];
    }

    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter) {
        int64 tempKey = iter->first;
        Object* tempObject = iter->second;
        int totalNumber = tempObject->getTotalNumber();
        while(totalNumber == 0) {
            iter++;
            removeObjectFromMap(tempKey);
            if (iter != allObjects.end()) {
                tempObject = iter->second;
                tempKey = iter->first;
                totalNumber = tempObject->getTotalNumber();
            }
            else {
                break;
            }
        }
        if (iter == allObjects.end()) {
            break;
        }

        for (int i = affectedStart; i <= affectedEnd; i++) {
            Bundle* tempBundle = linePool[tempObject];
            OneLine* tempLine = tempBundle->lines[i];
            if (tempLine != nullptr) {
                matrixRate[i] += tempLine->computeTotalRate();
                //tempLine->display(tempObject);/* Qianran 0925 */
            }
        }

    }
}

void SCDWrapper::computeBulkRate()
{
    bulkRate = 0.0;
    for (int i = 0; i < POINTS; ++i) {
        bulkRate += matrixRate[i];
    }
}

void SCDWrapper::computeDomainRate()
{
    domainRate = 0.0;
    for (int i = startIndex; i <= endIndex; i++)
    {
        domainRate += matrixRate[i];
    }
    domainRate += noneRate;
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

Object* SCDWrapper::selectReaction(
                                   int64& theOtherKey,
                                   Reaction& reaction,
                                   int& count)
{
    ofstream fs;
    fs.open("selectReaction.txt", ios::app);
    int pointIndex;

    // Generate a random number
    // Create a random device to seed the random number engine
    std::random_device rd;

    // Initialize the random number engine with the seed
    std::default_random_engine engine(rd());

    // Create a uniform real distribution that produces values in the range [0.0, 1.0)
    std::uniform_real_distribution<long double> distribution(0.0L, 1.0L);

    // Generate a random long double
    long double randomNum = distribution(engine);

    long double randRate = randomNum * bulkRate;
    // long double randRate = ((double)rand() / (RAND_MAX+1.0))*bulkRate; // add 1 to make the random rate [0, 1)
    long double tempRandRate = randRate;
    Object* tempObject = nullptr;
    Bundle* tempBundle;
    OneLine* tempLine;
    //fs << "BulkRate = " << bulkRate << "RandRate = " << randRate << endl;
    for (pointIndex = 0; pointIndex < POINTS; ++pointIndex) {
        if (matrixRate[pointIndex] < tempRandRate) {
            tempRandRate -= matrixRate[pointIndex];
            continue;
        }/* if the event is not positioned in this element, move on to the next one */
        else {
            reaction = NONE;
            unordered_map<int64, Object*>::iterator iter = allObjects.begin();
            while (reaction == NONE && iter != allObjects.end()) {
                tempObject = iter->second;
                tempBundle = linePool[tempObject];
                tempLine = tempBundle->lines[pointIndex];
                if (tempLine != nullptr) {
                    reaction = tempLine->selectReaction(tempObject, theOtherKey, tempRandRate);
                }
                ++iter;
            }
            if (reaction == NONE) {
                reaction = damage.selectDamage(pointIndex, tempRandRate);
            }
            if (reaction == NONE){
                if(sinkDissRate[0][pointIndex] >= tempRandRate){
                    reaction = DISSV;
                }else{
                    reaction = DISSH;
                }
                
            }
            
            count = pointIndex;
            if (LOG_REACTIONS)
                fs << "Element = " << pointIndex + 1 <<", "<< "Reaction = " << reaction << endl << endl;
            fs.close();
            return tempObject;
        }
    }
    // Sometimes there is a case where randRate == bulkRate, which might result in a 
    // tiny bit of leftover number (eg. doing 9e14-9e14 gives 0.0001 when it should give 0)
    // causing the program to think that no event was selected
    cout << "returned nullptr" << endl;
    cout << randRate << " " << bulkRate << " " << tempRandRate << endl;
    return tempObject;
}

Object* SCDWrapper::selectDomainReaction(
                                   int64& theOtherKey,
                                   Reaction& reaction,
                                   int& count)
{
    ofstream fs;
    fs.open("selectReaction.txt", ios::app);
    int pointIndex;

    // Generate a random number
    // Create a random device to seed the random number engine
    std::random_device rd;

    // Initialize the random number engine with the seed
    std::default_random_engine engine(rd());

    // Create a uniform real distribution that produces values in the range [0.0, 1.0)
    std::uniform_real_distribution<long double> distribution(0.0L, 1.0L);

    // Generate a random long double
    long double randomNum = distribution(engine);

    long double randRate = randomNum * domainRate;
    // long double randRate = ((double)rand() / (RAND_MAX+1.0))*bulkRate; // add 1 to make the random rate [0, 1)
    long double tempRandRate = randRate;
    Object* tempObject = nullptr;
    Bundle* tempBundle;
    OneLine* tempLine;
    //fs << "BulkRate = " << bulkRate << "RandRate = " << randRate << endl;
    for (pointIndex = startIndex; pointIndex <= endIndex; ++pointIndex) {
        if (matrixRate[pointIndex] < tempRandRate) {
            tempRandRate -= matrixRate[pointIndex];
            continue;
        }/* if the event is not positioned in this element, move on to the next one */
        else {
            reaction = NONE;
            unordered_map<int64, Object*>::iterator iter = allObjects.begin();
            while (reaction == NONE && iter != allObjects.end()) {
                tempObject = iter->second;
                tempBundle = linePool[tempObject];
                tempLine = tempBundle->lines[pointIndex];
                if (tempLine != nullptr) {
                    reaction = tempLine->selectReaction(tempObject, theOtherKey, tempRandRate);
                }
                ++iter;
            }
            if (reaction == NONE) {
                reaction = damage.selectDamage(pointIndex, tempRandRate);
            }
            if (reaction == NONE){
                if(sinkDissRate[0][pointIndex] >= tempRandRate){
                    reaction = DISSV;
                }else{
                    tempRandRate -= sinkDissRate[0][pointIndex];
                    if (sinkDissRate[1][pointIndex] >= tempRandRate){
                        reaction = DISSH;
                    }
                }
            }
            
            count = pointIndex;
            if (LOG_REACTIONS)
                fs << "Element = " << pointIndex + 1 <<", "<< "Reaction = " << reaction << endl << endl;
            fs.close();
            return tempObject;
        }
    }
    // Sometimes there is a case where randRate == bulkRate, which might result in a 
    // tiny bit of leftover number (eg. doing 9e14-9e14 gives 0.0001 when it should give 0)
    // causing the program to think that no event was selected
    reaction = NONE;
    // cout << "returned nullptr" << endl;
    // cout << randRate << " " << domainRate << " " << tempRandRate << endl;
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
    fs.open("Reactions.txt", ios::app);
    switch (reaction) {
        case DIFFUSETOF:
            processDiffEvent(hostObject, n, 'f');
            if (LOG_REACTIONS)
                fs << hostObject->getKey() <<"  diffuses from element "<< n << " to element " << n-1 << endl;
            break;
        case DIFFUSETOB:
            processDiffEvent(hostObject, n, 'b');
            if (LOG_REACTIONS)
                fs << hostObject->getKey() <<"  diffuses from element "<< n << " to element " << n+1 << endl;
            break;
        case SINK:
            processSinkEvent(hostObject, n);
            if (LOG_REACTIONS)
                fs << hostObject->getKey() <<"  in element "<< n << " goes to sink."<< endl;
            writeSinkFile(hostObject, n, time); /* this step also updated sinkDissRate */
            break;
        case DISSOCIATION:
            processDissoEvent(hostObject, n, theOtherKey, fs);
            if (LOG_REACTIONS)
                fs << hostObject->getKey() <<"  experiences a dissociation." << endl;
            break;
        case COMBINATION:
            processCombEvent(hostObject, n, theOtherKey, fs);
            break;
        case SAV:
            processSAVEvent(hostObject, n);
            if (LOG_REACTIONS)
                fs << hostObject->getKey() << "  in element " << n << " experiences SAV (ejects interstitial)" << endl;
            break;
        case RECOMBER:
            processRecombEvent(hostObject, n, true);
            if (LOG_REACTIONS)
                fs << hostObject->getKey() << "  in element " << n << " experiences an ER recombination to form H2 and leave the front of the material" << endl;
            break;
        case RECOMBLH:
            processRecombEvent(hostObject, n, false);
            if (LOG_REACTIONS)
                fs << hostObject->getKey() << "  in element " << n << " experiences a LH recombination to form H2 and leave the back of the material" << endl;
            break;
        case PARTICLE:
            getParticleInsertion(n, dt, fs);
            break;
        case HE:
            getHeInsertion(n);
            break;
        case H:
            getHInsertion(n, dt, fs);
            break;
        case DISSV:
            processSinkDissEvent(0, n);
            dissV++;
            //cout << "dissV = " << dissV << endl;
            break;
        case DISSH:
            processSinkDissEvent(1, n);
            dissH++;
            //cout << "dissH = " << dissH << endl;
            break;
        default:
            break;
    }

    // Keep track of affected reaction rates
    updateMatrixRate(n, reaction);
    computeBulkRate();

    fs.close();
}

SCDWrapper::~SCDWrapper()
{
    /* clean linePool*/
    unordered_map<Object*, Bundle*>::iterator iter;
    unordered_map<int64, Object*>::iterator iter1;
    for (iter = linePool.begin(); iter != linePool.end(); ++iter) {
        delete iter->second;
    }
    /* clean contents of object* */
    for (iter1 = allObjects.begin(); iter1 != allObjects.end(); ++iter1) {
        delete iter1->second;
    }
    for (iter1 = HObjects.begin(); iter1 != HObjects.end(); ++iter1) {
        delete iter1->second;
    }
}

unordered_map<int64, Object*>* SCDWrapper::getAllObjects()
{
    return &allObjects;
}

unordered_map<int64, Object*>* SCDWrapper::getMobileObjects()
{
    return &mobileObjects;
}

unordered_map<Object*, Bundle*>* SCDWrapper::getLinePool()
{
    return &linePool;
}

void SCDWrapper::examineRate()
{
    fstream fs;
    //fs.open("Rate.txt",ios::app);
    for (int i = 0; i < POINTS; ++i) {
        computeMatrixRate(i);
        //fs << matrixRate[i] <<"        ";
    }
    // computeBulkRate();
    //fs << "BulkRate" << bulkRate << endl<<endl<<endl;
    //fs.close();
}

void SCDWrapper::writeSinkFile(const Object * const hostObject, const long int n, const double time)
{
    int i;
    for (i = 0; i < LEVELS; i++) {
        int number = hostObject->getAttri(i);
        if (i == 0 && number<0) {
            sinks[i][n] += abs(number);
            computeSinkDissRate(i, n);
        }
        else {
            sinks[i + 1][n] += number;
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
        fluenceH = FLUX_H * time; // Only take simulation progress parameters from thread id 0
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
    ofstream fs;
    fs.open("sink" + std::to_string(threadID) + ".txt", ios::out);
    fs << "startIndex = " << startIndex << endl;
    fs << "endIndex = " << endIndex << endl;
    //fs << "Aggregate time = " << time << "  step = " << n << endl;
    for (i = 0; i < POINTS; ++i) {
        for (j = 0; j < LEVELS + 1; ++j) {
            fs << sinks[j][i]<< "    ";
        }
        fs << endl;
    }
    fs.close();
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
    for (i = 0; i < LEVELS + 1; ++i){
        for (j = 0; j < POINTS; ++j) {
            sinks[i][j] = 0;
        }
    }
    for (i = 0; i < 2; ++i){
        for (j = 0; j < POINTS; ++j) {
            sinkDissRate[i][j] = 0;
        }
    }
}

void SCDWrapper::computeSinkDissRate(const int type, const int point)
{
    double b = 2.8e-8; //burger's vector 2.8e-8 cm
    // double nu = 1e+13; // attempt frequency s^-1 averaged with dt
    double ebH = 0.6, emH = 0.25; //binding and migration energy of hydrogen
    double ebV = 1.0, emV = 1.29; //binding and migration energy of vacancy
    double prefactor = 1.58e-3;
    double prefactorv = 177e-4;
    // vacancy dissociation
    if(type == 0){
        sinkDissRate[type][point] = (2*PI*DISLOCATION*VOLUME/ALATT/ALATT/b)*sinks[0][point]*prefactorv*exp(-(ebV+emV)/KB/TEMPERATURE);
        //sinkDissRate[type][point] = 0.0;
    }
    // hydrogen dissocaiton
    if(type == 1){
        sinkDissRate[type][point] = (2*PI*DISLOCATION*VOLUME/ALATT/ALATT/b)*sinks[3][point]*prefactor*exp(-(ebH+emH)/KB/TEMPERATURE);
    }
}

int64 SCDWrapper::attrToKey(const int * const attr)
{
    int64 key = 0;
    int sign = (attr[0] < 0) ? -1 : 1;
    for (int i = 0; i < LEVELS; ++i) {
        key += labs(attr[i])*((int64)pow(10.0, (double)EXP10*(LEVELS - 1 - i)));
    }
    key *= sign;
    return key;
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
            addReactionToOther(newObject); /* add to other objects' lines */
        } /* add to mobile object if necessary */
        if (newObject->getAttri(0) == 0 && newObject->getAttri(2) > 0) {
            HObjects.insert(newNode);
        } /* Keep track of nH objects */
        Bundle* newBundle = new Bundle(newObject, mobileObjects, allObjects);
        pair<Object*, Bundle*> bundle(newObject, newBundle);
        linePool.insert(bundle); /* add to line pool */
    }/* if this object is valid, add it to map */
}

void SCDWrapper::addToObjectMap(const int64 key, const int n, const int number)
{
    /* If the object exists, add to it. Otherwise create the object. */
    Object* anObject;
    if (allObjects.find(key) != allObjects.end()) 
    {
        /* object found! then number of instances in this element increases by number*/
        anObject = allObjects[key];
        anObject->addNumber(n, number);
        updateObjectInMap(anObject, n);
    }
    else if (number > 0)
    {
        /* object wasn't found! build new object and insert it into map */
        anObject = new Object(key, n, number);
        addNewObjectToMap(anObject);
    }
    else
    {
        return;
    }

    bool leftBoundary = (
        (n == startIndex || n == startIndex - 1) && n != 0
    );

    bool rightBoundary = (
        (n == endIndex || n == endIndex + 1) && n != POINTS - 1
    );

    if (leftBoundary)
        leftBoundaryChangeQ.push_back(
            BoundaryChange(key, n, number));
    else if (rightBoundary)  
        rightBoundaryChangeQ.push_back(
            BoundaryChange(key, n, number)); 
}

void SCDWrapper::reduceFromObjectMap(const int64 key, const int n, const int number)
{
    /* Count is mesh element index, number is number to add to object population */
    addToObjectMap(key, n, -number);
}

void SCDWrapper::updateObjectInMap(Object * hostObject, const int count)
{
    // Update the OneLine associated with this object in this count
    OneLine* tempLine = linePool[hostObject]->lines[count];
    double diffusivity = hostObject->getDiff();
    int number = hostObject->getNumber(count);
    if (tempLine != nullptr) {
        if (number > 0) {
            tempLine->updateLine(hostObject, count, mobileObjects, allObjects);
        }
        else {
            delete tempLine;
            linePool[hostObject]->lines[count] = nullptr;
        }
    }
    else {
        if (number > 0) {
            tempLine = new OneLine(hostObject, count, mobileObjects, allObjects);
            linePool[hostObject]->lines[count] = tempLine;
        }
    }

    // Update the OneLines of other objects impacted by this mobile object
    if (diffusivity > 0) {
        updateRateToOther(hostObject, count);

        // Update diffusion rates of this object in neighbouring elements
        if((count-1) >= 0){
            OneLine* tempLine = linePool[hostObject]->lines[count - 1];
            if(tempLine != nullptr){
                tempLine->updateDiff(hostObject, count - 1, allObjects);
            }
        }
        if((count + 1) < POINTS){
            OneLine* tempLine = linePool[hostObject]->lines[count + 1];
            if(tempLine != nullptr){
                tempLine->updateDiff(hostObject, count + 1, allObjects);
            }
        }
    }

    if (hostObject->getKey() == 1) // if is 1H object
    {
        damage.updateDamageTwo(count, allObjects);
    }

    // updateNetCombDissRate(hostObject, count);
}

void SCDWrapper::addReactionToOther(Object const * const mobileObject)
{
    unordered_map<Object*, Bundle*>::iterator iter;
    for (iter = linePool.begin(); iter != linePool.end(); ++iter) {
        Object* hostObject = iter->first;
        Bundle* tempBundle = iter->second;
        for (int i = 0; i < POINTS; ++i) {
            OneLine* tempLine = tempBundle->lines[i];
            if (tempLine != nullptr) {
                tempLine->addReaction(hostObject, mobileObject, i);
            }
        }
    }
}

void SCDWrapper::updateRateToOther(Object const * const mobileObject, const int count)
{
    unordered_map<Object*, Bundle*>::iterator iter;
    for (iter = linePool.begin(); iter != linePool.end(); ++iter) {
        Object* hostObject = iter->first;
        Bundle* tempBundle = iter->second;
        OneLine* tempLine = tempBundle->lines[count];
        if (tempLine != nullptr) {
            /* If both are mobile objects, mobileObject will have already recorded the combination rate, no need to record it again */
            if (hostObject->getDiff() > 0 && mobileObject->getKey() != hostObject->getKey()) {
                tempLine->setCombReaction(mobileObject->getKey(), 0.0);
            }
            else {
                tempLine->updateReaction(hostObject, mobileObject, count);
            }
        }
    }
}

void SCDWrapper::removeObjectFromMap(const int64 deleteKey)
{
    Object* deleteObject = allObjects[deleteKey];
    Bundle* tempBundle = linePool[deleteObject];
    double diffusivity = deleteObject->getDiff();
    bool is_nH = deleteObject->getAttri(0) == 0 && deleteObject->getAttri(2) > 0;  // if it's an nH object
    delete tempBundle; /* delete bundle */
    linePool.erase(deleteObject); /* remove this object from map linePool */
    delete deleteObject;  /* delete the content of this object */
    allObjects.erase(deleteKey); /* delete this object from map allObjects */
    if (diffusivity > 0) {
        mobileObjects.erase(deleteKey); /* delete this object from map mobileObjects */
        removeRateToOther(deleteKey);
    }
    if (is_nH) {
        HObjects.erase(deleteKey);
    }
}

void SCDWrapper::removeRateToOther(const int64 deleteKey)
{
    unordered_map<int64, Object*>::iterator iter;
    for (iter = allObjects.begin(); iter != allObjects.end(); ++iter) {
        Object* tempObject = iter->second;
        Bundle* tempBundle = linePool[tempObject];
        for (int i = 0; i < POINTS; ++i) {
            OneLine* tempLine = tempBundle->lines[i];
            if (tempLine != nullptr) {
                tempLine->removeReaction(deleteKey);
            }
        }
    }
}

void SCDWrapper::updateCombDissRatePair(int combinedObjectAttr[], const int count)
{
    /*
     * Find and update the net combination rate of a species that combines and its conjugate product species that dissociates
     * For now, only do this calculation for H + SIA/V cluster, because that is the 
     * most common reaction.
     *
     * combinedObject - the Object that results after a combination event
     */

    int attrIndex = 2;   /* For now just focus on H monomer */
    if (combinedObjectAttr[attrIndex] <= 0)
    {
        // cout << "updateCombDissRatePair is only meant to be used with the product of cluster + H reaction, this function is being used wrong" << endl;
        return;
    }

    int64 HKey = 1;
    int prevObjectAttr[LEVELS] = {0}; /* The Object that combines with H to form the combinedObject */
    for (int i = 0; i < LEVELS; i++)
        prevObjectAttr[i] = combinedObjectAttr[i];
    prevObjectAttr[attrIndex] -= 1;  /* The object from before it combined with H */

    int64 prevObjectKey = attrToKey(prevObjectAttr);
    int64 combinedObjectKey = attrToKey(combinedObjectAttr);

    bool foundHObj = allObjects.find(HKey) != allObjects.end() && allObjects[HKey]->getNumber(count) > 0;
    bool foundPrevObj = allObjects.find(prevObjectKey) != allObjects.end() && allObjects[prevObjectKey]->getNumber(count) > 0;
    bool foundCombSpecies = foundHObj && foundPrevObj;
    bool foundCombResult = allObjects.find(combinedObjectKey) != allObjects.end() && allObjects[combinedObjectKey]->getNumber(count) > 0;
    Object *hostObject, *mobileObject, *combinedObject;
    OneLine *hostLine, *combinedLine;

    if (!foundHObj && !foundPrevObj && !foundCombResult)
        return;  /* All species don't exist, nothing to compare */

    if (foundHObj)
        mobileObject = allObjects[HKey];
    else
        mobileObject = new Object(HKey, count);

    if (foundPrevObj)
    {
        hostObject = allObjects[prevObjectKey];
        hostLine = linePool[hostObject]->lines[count];
    }
    else
    {
        hostObject = new Object(prevObjectKey, count);
        hostLine = new OneLine();
    }

    if (foundCombResult)
    {
        combinedObject = allObjects[combinedObjectKey];
        combinedLine = linePool[combinedObject]->lines[count];
    }
    else
    {
        combinedObject = new Object(combinedObjectKey, count);
        combinedLine = new OneLine();
    }

    if (foundCombSpecies && foundCombResult)
    {
        long double baseCombRate = 0.0;
        long double baseDissRate = 0.0;
        if (hostLine != nullptr)
            baseCombRate = hostLine->computeCombReaction(hostObject, mobileObject, count);
        if (combinedLine != nullptr)
            baseDissRate = combinedLine->computeDissReaction(combinedObject, attrIndex, count);
        if (baseCombRate > baseDissRate)
        {
            if (hostLine != nullptr && foundCombSpecies)
                hostLine->setCombReaction(mobileObject->getKey(), baseCombRate - baseDissRate);
            if (combinedLine != nullptr && foundCombResult)
                combinedLine->setDissReaction(attrIndex, 0.0);
        }
        else  /* diss >= comb */
        {
            if (hostLine != nullptr && foundCombSpecies)
                hostLine->setCombReaction(mobileObject->getKey(), 0.0);
            if (combinedLine != nullptr && foundCombResult)
                combinedLine->setDissReaction(attrIndex, baseDissRate - baseCombRate);
        }
    }

    if (!foundHObj)
        delete mobileObject;
    if (!foundPrevObj)
    {
        delete hostObject;
        delete hostLine;
    }
    if (!foundCombResult)
    {
        delete combinedObject;
        delete combinedLine;
    }
}

void SCDWrapper::updateNetCombDissRate(const Object* const hostObject, const int count)
{
    /* 
     * Find the net combination rate (comb - diss, or vice versa if diss >= comb).
     * For now, only do this calculation for H + SIA/V cluster, because that is the 
     * most common reaction.
     */

    /* Update this comb rate and the product's diss rate */
    int64 HKey = 1;
    int attrIndex = 2; /* for now, only focus on H + cluster comb and diss reactions */
    int originalAttr[LEVELS] = {0};
    int productAttr[LEVELS] = {0};

    if (hostObject->getKey() != HKey)
    {
        for (int i = 0; i < LEVELS; i++)
        {
            originalAttr[i] = hostObject->getAttri(i);
            productAttr[i] = hostObject->getAttri(i);
        }
        productAttr[attrIndex]++;   /* Product of this cluster + H comb reaction */
        updateCombDissRatePair(productAttr, count);  /* Update this comb rate, and next object's diss rate */
        updateCombDissRatePair(originalAttr, count); /* Update this diss rate, and prev object's comb rate */ 
    }
    else /* hostObject is 1H object */
    {
        for (unordered_map<int64, Object*>::iterator iter = allObjects.begin(); iter != allObjects.end(); ++iter)
        {
            Object* otherObject = iter->second;
            for (int i = 0; i < LEVELS; i++)
                productAttr[i] = hostObject->getAttri(i) + otherObject->getAttri(i);
            updateCombDissRatePair(productAttr, count); /* Update this comb rate, and next object's diss rate. No need to update this object's diss rate because a change in the 1H species count doesn't affect this object's diss rate directly. */
        }
    }
}

void SCDWrapper::updateSinks(const int point, const int* number){
    for(int i = 0; i< LEVELS+1; i++){
        sinks[i][point] = number[i];
        if(i==0){
            computeSinkDissRate(i,point);
        }else if(i==3){
            computeSinkDissRate(1,point);
        }
    }
}

/* private function */
void SCDWrapper::processDiffEvent(Object* hostObject, const int n, const char signal)
{
    int64 key = hostObject->getKey();
    
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
    reduceFromObjectMap(key, n);
}

void SCDWrapper::processSinkEvent(Object * hostObject, const int n)
{
    ++reactions[2][n];
    reduceFromObjectMap(hostObject->getKey(), n);
    /*
    if(key == -1000000){
        //case 1V
        sinks[0][n]++;
        computeSinkDissRate(0, n);
    }else if(key == 1){
        // case 1H
        sinks[3][n]++;
        computeSinkDissRate(1, n);
    }
    */
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
    /* 1) deal with monomer */
    addToObjectMap(monomerKey, n);
    /* 2) generate the other cluster */
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

    /* 3) deal with the host object */
    reduceFromObjectMap(hostObject->getKey(), n);
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
    else if(hostObject->getAttri(0)<0 && hostObject->getAttri(2)>0 && theOtherObject->getAttri(0)>0 && theOtherObject->getAttri(2)==0 && productAttr[0] < 0 && productAttr[2] > 4.75 + 4*abs(productAttr[0])){
        /* Vn-Hm + SIAx -> V(n-x)-Hm (n, m, and x are not zero) */
        /* change above reaction to Vn-Hm + SIAx -> V(n-x)-H(max) + (excess)*H */
        int maxH = abs(productAttr[0]) * 4;
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
    /* update reactant */
    reduceFromObjectMap(hostObject->getKey(), n);
    reduceFromObjectMap(theOtherObject->getKey(), n);

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
    const int64 HKey = 1;
    double HConcentration = 0;
    unordered_map<int64, Object*>::iterator iter;
    double volume = VOLUME;
    if (n == 0)
        volume = SURFACE_VOLUME;  // first mesh element is slightly bigger because of added surface layer
    for (iter = HObjects.begin(); iter != HObjects.end(); ++iter)
    {
        Object* HObject = iter->second;
        HConcentration += HObject->getNumber(n) * HObject->getAttri(2) / volume;
    }

    double saturation_concentration = getHSaturationConcentration();

    // If it's not saturated, SAV doesn't need to create more vacancies
    if (HConcentration < saturation_concentration)
    {
        return;
    }

    // Eject interstitial
    int64 SIAKey = (int64)pow(10.0, (double)EXP10 * (LEVELS - 1)); /* Key for SIA. */
    addToObjectMap(SIAKey, n);

    // Generate 1 vacancy
    int productAttr[LEVELS] = { 0 };
    for (int i = 0; i < LEVELS; i++)
        productAttr[i] = hostObject->getAttri(i);

    productAttr[0] -= 1;

    if (productAttr[0] < 0 && productAttr[2] > abs(productAttr[0])*4)
    {
        /* If the product has too many H per vacancy, eject the rest of the H */
        int excessH = productAttr[2] - abs(productAttr[0])*4;
        productAttr[2] -= excessH;
        addToObjectMap(HKey, n, excessH);
    }

    int64 productKey = attrToKey(productAttr);

    addToObjectMap(productKey, n);

    reduceFromObjectMap(hostObject->getKey(), n);
}

void SCDWrapper::processRecombEvent(Object* hostObject, const int n, bool ER)
{
    if (n != 0)
        cerr << "Recomb Error" << endl;
    /* 
     * Recombination: Two 1H instances combine to form H2 molecule, which
     * leaves the material through either the front or the back of the material 
     */

    // The implementation is easy: just remove one 1H object for ER Recomb and two 1H objects for LH Recomb
    if (ER)
        reduceFromObjectMap(hostObject->getKey(), n, 1);
    else
        reduceFromObjectMap(hostObject->getKey(), n, 2);
}

void SCDWrapper::processSinkDissEvent(const int type, const int point)
{
    // dissV event
    int64 productKey = 0;
    if(type == 0){
        sinks[0][point]--;
        productKey = -1000000; //1V
    }else if(type ==1){
        //dissV event
        sinks[3][point]--;
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

void SCDWrapper::getNeutronInsertion(const int n)
{
    double neutronEnergy = 0.0;
    double totalEnergy = (double)Poisson(AVG_NEUTRON_EN);
    CascadeDamage damage;
    int i, j, ndef = 0;
    int64 clusterKey;
    while (neutronEnergy < totalEnergy) {
        double pkaEnergy = 0.0;
        while (pkaEnergy < 0.62) {
            // Limit to produce a stable Frenkel pair in keV (from Troev et al (2011)).
            double xi = (rand() / RAND_MAX) * cpdf.getMaxPossibility(n);
            neutronEnergy += pkaEnergy;
            pkaEnergy = cpdf.samplePkaEnergy(xi, n);
        }
        damage.generateNeutronDamage(pkaEnergy, ndef);
        for (i = 0; i < damage.size(); ++i) {
            int sign = (i == 0) ? -1 : 1;
            for (j = 0; j < ndef; ++j) {
                const int number = damage.getDamage(i, j);
                if (number != 0) {
                    clusterKey = (int64)sign*(j + 1)*(pow(10.0, (double)EXP10*(LEVELS - 1)));
                    /* generate cluster key */
                    addToObjectMap(clusterKey, n, number);
                }
            }
        }
        damage.cleanDamage();
        neutronEnergy += pkaEnergy;
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
    fluenceH += FLUX_H*dt; /* 4.00e+16 is H flux */

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
    int numberSinks[LEVELS+1] = { 0 };
    int number[POINTS] = { 0 };
    int step = 0;
    string skip;
    string oneLine;
    stringstream lineHold;
    /* update sink numbers */
    ifstream file("sink.txt");
    if (file.good()) {
        for(int j = 0; j < POINTS; j++){
            while (getline(file, oneLine)) {
                lineHold.str(oneLine);
                for (int i = 0; i < LEVELS+1; i++) {
                    lineHold >> numberSinks[i];
                }
                srscd->updateSinks(j,numberSinks);
                lineHold.clear();
            }
        }
    }
    file.close();
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
        double volume;
        if(i==0){
            volume = SURFACE_VOLUME; /* volume at front is thin surface layer */
        }else{
            volume = VOLUME;
        }
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
    return H_SATURATION_CONCENTRATION;
    double concentration = H_SATURATION_CONCENTRATION;
    bool dimer = false;
    bool vhpair = false;

    for (unordered_map<int64, Object*>::const_iterator iter = allObjects.begin(); iter != allObjects.end() && !dimer && !vhpair; iter ++)
    {
        Object* tempObj = iter->second;
        // if HH objects is present

        if (tempObj->getKey() == 2 && !dimer)
        {
            concentration += DENSITY * 8 * exp(-(2 * H_FORM_E - HH_BIND_E) / KB / TEMPERATURE);
            dimer = true;
        }
        // if mV-nH objects are present
        if (tempObj->getAttri(0) <= -1 && tempObj->getAttri(2) >= 1 && !vhpair)
        {
            concentration += DENSITY * 8 * exp(-(H_FORM_E + V_FORM_E - VH_BIND_E) / KB / TEMPERATURE);
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
}

void SCDWrapper::fillNoneReaction(const double& maxDomainRate)
{
    /*
     * Fill the remainder of the domain rate with NONE reaction
     * so that each processor can move at the same time step in parallel
     * (Dunn 2016)
     */
    noneRate = maxDomainRate - domainRate;
    computeDomainRate();
}

void SCDWrapper::clearNoneReaction()
{
    noneRate = 0.0;
    computeDomainRate();
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
    for (int level = 0; level < LEVELS+1; level++)  // levels + 1 because vacancy and sia each have their own levels, so vac + sia + helium + H = 4 levels = 3 + 1
    {
        output[level] = sinks[level][n];
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

    /* Put it new sink counts into our new boundary index */
    for (int level = 0; level < LEVELS+1; level++)
    {
        sinks[level][newBoundaryIndex] = newBoundarySinks[level];
        if (level == 0)
            computeSinkDissRate(0, newBoundaryIndex);  // V Diss rate
        else if (level == 3)
            computeSinkDissRate(1, newBoundaryIndex);  // H Diss rate
    }

    clearBoundaryChangeQs();  // from the add/reduceToObjectMap calls
}
