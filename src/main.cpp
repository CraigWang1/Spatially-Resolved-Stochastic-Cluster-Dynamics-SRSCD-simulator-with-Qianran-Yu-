// SRSCD -- Spatially Resolved Stochastic Cluster Dynamics
// Q Yu, January 2017.

#include<mpi.h>
#include<ctime>
#include<unistd.h>
#include<time.h>
#include <cassert>
#include"SCDWrapper.h"

int main(int argc, char** argv) 
{
    int threadID, numThreads;
    int rootThreadID = 0;
    int tag = 1;
    int lengthTag = 2;
    int dataTag = 3;
    MPI_Status status;
    MPI_Request requests[4];
    int intPerMessage = 3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);

    int rightNeighbor = threadID + 1;
    int leftNeighbor = threadID - 1;

    SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd */
    int64 theOtherKey = 0;
    Object* hostObject = nullptr;
    Reaction reaction = ERROR;
    vector<BoundaryChange> receivedBoundaryChanges;
    int pointIndex = -1;
    long int iStep = 0;
    double random;
    double advTime = 0.0;
    double dt = 0.0;
    double system_dt = 0.0;
    double accTime = 0.0;
    double dpa = 0.0;
    double progress = 0.0; /* progress of 50 = 50% done with simulation */
    double eta_min = 0.0;
    double prev_progress = 0.0;
    double prev_eta_min = 0.0;
    double write_time = 0.2;
    double write_increment = 0.2;
    long double avgDomainRate = 0.0;
    unsigned int transferTurn = 0;
    int barWidth, pos; /* progress bar parameters */
    int numDigits = 7;
    int magnitude = 0;
    bool done = false;
    fstream st;
    /* check whether to restart*/
    restart(iStep, advTime, srscd);
    while (advTime > write_time)
        write_time += write_increment;
    srscd->examineRate();
    /* check ended */
    srand(time(0));
    if (threadID == 0)
    {
        srscd->displayAllObject();
        srscd->drawSpeciesAndReactions(advTime);
    }

    // Assign volume elements to this processor
    double indexIncrement = (double) POINTS / numThreads;
    int startIndex = threadID * indexIncrement;
    int endIndex = (threadID + 1) * indexIncrement - 1;
    srscd->setDomain(startIndex, endIndex);
    srscd->clearNoneReaction();
    srscd->examineDomainRate();

    double prev_time = MPI_Wtime();
    cout << H_SATURATION_CONCENTRATION * VOLUME << endl;
    
    while(!done)
    {
        long double localDomainRate = srscd->getDomainRate();
        long double maxDomainRate;

        // Each MPI process sends its rate to reduction, every thread collects result
        MPI_Allreduce(&localDomainRate, &maxDomainRate, 1, MPI_LONG_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        srscd->fillNoneReaction(maxDomainRate);
        hostObject = srscd->selectDomainReaction(theOtherKey, reaction, pointIndex);/* choose an event */

        /*
        string reactions[] = {"diffF", "diffB", "sink", "diss", "comb", "sav", "recombER", "recombLH", "none", "particle", "HE", "H", "dissV", "dissH", "error"};        
        cout << reactions[reaction] << " " << pointIndex << " " << theOtherKey;
        if (hostObject != nullptr)
            cout << " " << hostObject->getKey();
        cout << endl;
        */

        srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime, accTime); /* process event */

        if(reaction == Reaction::H || reaction == Reaction::PARTICLE)
        {
            accTime = 0.0;
        }

        // Update average domain rate (running average)
        avgDomainRate += 0.0001*(localDomainRate - avgDomainRate);

        // Even processors send to right, odd processors receive from left
        vector<BoundaryChange>* leftBoundaryChanges = srscd->getLeftBoundaryChangeQ();
        vector<BoundaryChange>* rightBoundaryChanges = srscd->getRightBoundaryChangeQ();
        int numToSendLeft = leftBoundaryChanges->size();
        int numToSendRight = rightBoundaryChanges->size();
        int numToRecvLeft = 0, numToRecvRight = 0;

        if (leftNeighbor >= 0)
        {
            MPI_Isend(&numToSendLeft, 1, MPI_INT, leftNeighbor, lengthTag, MPI_COMM_WORLD, &requests[0]);
            if (numToSendLeft > 0)
            {
                long int data[numToSendLeft*intPerMessage];
                for (size_t i = 0; i < leftBoundaryChanges->size(); i++)
                {
                    BoundaryChange& bc = leftBoundaryChanges->at(i);
                    data[intPerMessage*i] = bc.objKey;
                    data[intPerMessage*i+1] = bc.pointIndex;
                    data[intPerMessage*i+2] = bc.change; 
                }
                MPI_Isend(data, numToSendLeft*intPerMessage, MPI_LONG, leftNeighbor, dataTag, MPI_COMM_WORLD, &requests[1]);
            }
        }
        if (rightNeighbor < numThreads)
        {
            MPI_Isend(&numToSendRight, 1, MPI_INT, rightNeighbor, lengthTag, MPI_COMM_WORLD, &requests[2]);
            if (numToSendRight > 0)
            {
                long int data[numToSendRight*intPerMessage];
                for (size_t i = 0; i < rightBoundaryChanges->size(); i++)
                {
                    BoundaryChange& bc = rightBoundaryChanges->at(i);
                    data[intPerMessage*i] = bc.objKey;
                    data[intPerMessage*i+1] = bc.pointIndex;
                    data[intPerMessage*i+2] = bc.change; 
                }
                MPI_Isend(data, numToSendRight*intPerMessage, MPI_LONG, rightNeighbor, dataTag, MPI_COMM_WORLD, &requests[3]);
            }
        }
        if (leftNeighbor >= 0)
        {
            MPI_Recv(&numToRecvLeft, 1, MPI_INT, leftNeighbor, lengthTag, MPI_COMM_WORLD, &status);
            if (numToRecvLeft > 0)
            {
                long int data[intPerMessage*numToRecvLeft];
                MPI_Recv(data, intPerMessage*numToRecvLeft, MPI_LONG, leftNeighbor, dataTag, MPI_COMM_WORLD, &status);
                for (int i = 0; i < numToRecvLeft; i++)
                {
                    BoundaryChange bc;
                    bc.objKey = data[intPerMessage*i];
                    bc.pointIndex = data[intPerMessage*i+1];
                    bc.change = data[intPerMessage*i+2];
                    receivedBoundaryChanges.push_back(bc);
                }
            }
        }
        if (rightNeighbor < numThreads)
        {
            MPI_Recv(&numToRecvRight, 1, MPI_INT, rightNeighbor, lengthTag, MPI_COMM_WORLD, &status);
            if (numToRecvRight > 0)
            {
                long int data[intPerMessage*numToRecvRight];
                MPI_Recv(data, intPerMessage*numToRecvRight, MPI_LONG, rightNeighbor, dataTag, MPI_COMM_WORLD, &status);
                for (int i = 0; i < numToRecvRight; i++)
                {
                    BoundaryChange bc;
                    bc.objKey = data[intPerMessage*i];
                    bc.pointIndex = data[intPerMessage*i+1];
                    bc.change = data[intPerMessage*i+2];
                    receivedBoundaryChanges.push_back(bc);
                }
            }
        }

        srscd->clearBoundaryChangeQs();

        if (receivedBoundaryChanges.size() > 0)
        {
            srscd->implementBoundaryChanges(receivedBoundaryChanges);
            receivedBoundaryChanges.clear();
        }

        if (IRRADIATION_ON)
        {
            double localDpa = srscd->getTotalDpa();
            MPI_Reduce(&localDpa, &dpa, 1, MPI_DOUBLE, MPI_SUM, rootThreadID, MPI_COMM_WORLD);
        }

        if (threadID == rootThreadID)
        {
            do {
                random = (double)rand() / RAND_MAX;
            } while (random == 0);
            dt = (-1) / maxDomainRate*log(random);
            accTime += dt;
            advTime += dt;

            if (IRRADIATION_ON)
            {
                done = (dpa >= TOTAL_DPA);
            }
            else
            {
                done = (advTime >= TOTAL_TIME);
            }
            if (done)
            {
                cout << "Finished. Bye" << endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }

        ++iStep;
        if(iStep%PSTEPS == 0)
        {
            srscd->writeFile(advTime, iStep, threadID);

            if (threadID == rootThreadID)
            {
                system_dt = (MPI_Wtime() - prev_time);
                prev_time = MPI_Wtime();

                // Chose between dpa or time to calculate progress
                if (IRRADIATION_ON)
                {
                    progress = (dpa / TOTAL_DPA) * 100.;
                }
                else
                {
                    progress = (advTime / TOTAL_TIME) * 100.;
                }

                // Initialize the first eta estimate, or update the previous one
                if (prev_progress != 0 && progress - prev_progress != 0)
                {
                    eta_min = (100 - progress) / ((progress - prev_progress) / system_dt) / 60.;
                    
                    if (prev_eta_min != 0)
                    {
                        eta_min = prev_eta_min + 0.1 * (eta_min - prev_eta_min); // moving average

                        // Print progress bar
                        cout << "[";
                        barWidth = 70;
                        pos = barWidth * (progress/100.);
                        for (int i = 0; i < barWidth; i++)
                        {
                            if (i < pos) cout << "=";
                            else if (i == pos) cout << ">";
                            else cout << " ";
                        }

                        // Print a set amount of digits
                        magnitude = 0;
                        while ((int)(advTime / pow(10, magnitude)))
                        {
                            magnitude++;
                        }
                        cout << "] " << std::fixed << std::setprecision(2) << progress << "%";
                        cout << "   eta: " << (long int) round(eta_min) << " min";
                        cout << "   time: " << std::fixed << std::setprecision(numDigits - magnitude) << advTime << " s            \r";
                        cout.flush();
                    }
                }

                prev_eta_min = eta_min;
                prev_progress = progress;
                st.close();
            }

            // Transfer spatial elements between processors to increase parallel efficiency
            if (iStep%(10*PSTEPS) == 0)
            {
                if (transferTurn % 2 == 0)
                {
                    // Even processors transfer spatial element with left processor if necessary
                    if (threadID % 2 == 0 && leftNeighbor >= 0)
                    {
                        long double otherAvgDomainRate = 0;
                        MPI_Send(&avgDomainRate, 1, MPI_LONG_DOUBLE, leftNeighbor, tag, MPI_COMM_WORLD); // Don't give spatial element if we only have one remaining
                        MPI_Recv(&otherAvgDomainRate, 1, MPI_LONG_DOUBLE, leftNeighbor, tag, MPI_COMM_WORLD, &status);

                        bool canGiveSpatialElement = srscd->getStartIndex() < srscd->getEndIndex();
                        bool otherCanGiveMeSpatialElement;
                        MPI_Send(&canGiveSpatialElement, 1, MPI_CXX_BOOL, leftNeighbor, tag, MPI_COMM_WORLD);
                        MPI_Recv(&otherCanGiveMeSpatialElement, 1, MPI_CXX_BOOL, leftNeighbor, tag, MPI_COMM_WORLD, &status);

                        if (canGiveSpatialElement &&
                            avgDomainRate > otherAvgDomainRate)
                        {
                            /* If we can give a spatial element to the left neighbor */
                            /* Give up a mesh element by transferring information of the spatial element just before the boundary, which is the other processor's new ghost region */
                            vector<BoundaryChange> transferChanges = srscd->getSpatialElement(srscd->getStartIndex()+1);
                            for (BoundaryChange& bc: transferChanges)
                            {
                                long int data[3] = {bc.objKey, bc.pointIndex, bc.change};
                                MPI_Send(data, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD);
                            }
                            long int data[3] = {};
                            MPI_Send(data, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD); // Send end of transmission message                        }
                            int sinks[LEVELS+1];
                            srscd->getSink(srscd->getStartIndex(), sinks);

                            /* Transfer sinks of the boundary element that we are going to give */
                            MPI_Send(sinks, LEVELS+1, MPI_INT, leftNeighbor, tag, MPI_COMM_WORLD);

                            /* Set up this new domain */
                            srscd->setDomain(srscd->getStartIndex()+1, srscd->getEndIndex());
                            srscd->examineDomainRate();
                        }
                        else if (otherCanGiveMeSpatialElement &&
                            otherAvgDomainRate > avgDomainRate)
                        {
                            /* If we can receive a mesh element from the left neighbor */
                            /* Receive a mesh element by receiving information of the spatial element one past the other processor's boundary (this new ghost region) */
                            long int recv[3];
                            MPI_Recv(recv, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD, &status);
                            vector<BoundaryChange> receivedTransferChanges;

                            while (recv[0] != 0) // while it's a valid message
                            {
                                receivedTransferChanges.push_back(BoundaryChange(recv[0], recv[1], recv[2]));
                                MPI_Recv(recv, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD, &status);
                            }

                            /* Receive sinks of the new boundary element the other processor gave us */
                            int newBoundarySinks[LEVELS+1];
                            MPI_Recv(newBoundarySinks, LEVELS+1, MPI_INT, leftNeighbor, tag, MPI_COMM_WORLD, &status);

                            /* Set up this new domain */
                            srscd->addSpatialElement(srscd->getStartIndex()-2, receivedTransferChanges, srscd->getStartIndex()-1, newBoundarySinks);
                            srscd->setDomain(srscd->getStartIndex()-1, srscd->getEndIndex());
                            srscd->examineDomainRate();
                        }
                    }
                    else if (threadID % 2 == 1 && rightNeighbor < numThreads)
                    {
                        long double otherAvgDomainRate = 0;
                        MPI_Recv(&otherAvgDomainRate, 1, MPI_LONG_DOUBLE, rightNeighbor, tag, MPI_COMM_WORLD, &status);
                        MPI_Send(&avgDomainRate, 1, MPI_LONG_DOUBLE, rightNeighbor, tag, MPI_COMM_WORLD); // Don't give spatial element if we only have one remaining

                        bool canGiveSpatialElement = srscd->getStartIndex() < srscd->getEndIndex();
                        bool otherCanGiveMeSpatialElement;
                        MPI_Recv(&otherCanGiveMeSpatialElement, 1, MPI_CXX_BOOL, rightNeighbor, tag, MPI_COMM_WORLD, &status);
                        MPI_Send(&canGiveSpatialElement, 1, MPI_CXX_BOOL, rightNeighbor, tag, MPI_COMM_WORLD);

                        if (otherCanGiveMeSpatialElement &&
                            avgDomainRate < otherAvgDomainRate)
                        {
                            /* If we can receive a mesh element from the right neighbor */
                            /* Receive a mesh element by receiving information of the spatial element one past the other processor's boundary (this new ghost region) */
                            long int recv[3];
                            MPI_Recv(recv, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD, &status);
                            vector<BoundaryChange> receivedTransferChanges;

                            while (recv[0] != 0) // while it's a valid message
                            {
                                receivedTransferChanges.push_back(BoundaryChange(recv[0], recv[1], recv[2]));
                                MPI_Recv(recv, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD, &status);
                            }

                            /* Receive sinks of the new boundary element the other processor gave us */
                            int newBoundarySinks[LEVELS+1];
                            MPI_Recv(newBoundarySinks, LEVELS+1, MPI_INT, rightNeighbor, tag, MPI_COMM_WORLD, &status);

                            /* Set up this new domain */
                            srscd->addSpatialElement(srscd->getEndIndex()+2, receivedTransferChanges, srscd->getEndIndex()+1, newBoundarySinks);
                            srscd->setDomain(srscd->getStartIndex(), srscd->getEndIndex()+1);
                            srscd->examineDomainRate();
                        }
                        else if (canGiveSpatialElement && 
                            avgDomainRate > otherAvgDomainRate)
                        {
                            /* If we can give a spatial element to the right neighbor */
                            /* Give up a mesh element by transferring information of the spatial element just before the boundary, which is the other processor's new ghost region */
                            vector<BoundaryChange> transferChanges = srscd->getSpatialElement(srscd->getEndIndex()-1);
                            for (BoundaryChange& bc: transferChanges)
                            {
                                long int data[3] = {bc.objKey, bc.pointIndex, bc.change};
                                MPI_Send(data, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD);
                            }
                            long int data[3] = {};
                            MPI_Send(data, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD); // Send end of transmission message                        }
                            int sinks[LEVELS+1];
                            srscd->getSink(srscd->getEndIndex(), sinks);

                            /* Transfer sinks of the boundary element that we are going to give */
                            MPI_Send(sinks, LEVELS+1, MPI_INT, rightNeighbor, tag, MPI_COMM_WORLD);

                            /* Set up this new domain */
                            srscd->setDomain(srscd->getStartIndex(), srscd->getEndIndex()-1);
                            srscd->examineDomainRate();
                        }
                    }
                }
                else
                {
                    // Even processors transfer spatial element with right processor if necessary
                    if (threadID % 2 == 0 && rightNeighbor < numThreads)
                    {
                        long double otherAvgDomainRate = 0;
                        MPI_Send(&avgDomainRate, 1, MPI_LONG_DOUBLE, rightNeighbor, tag, MPI_COMM_WORLD); // Don't give spatial element if we only have one remaining
                        MPI_Recv(&otherAvgDomainRate, 1, MPI_LONG_DOUBLE, rightNeighbor, tag, MPI_COMM_WORLD, &status);

                        bool canGiveSpatialElement = srscd->getStartIndex() < srscd->getEndIndex();
                        bool otherCanGiveMeSpatialElement;
                        MPI_Send(&canGiveSpatialElement, 1, MPI_CXX_BOOL, rightNeighbor, tag, MPI_COMM_WORLD);
                        MPI_Recv(&otherCanGiveMeSpatialElement, 1, MPI_CXX_BOOL, rightNeighbor, tag, MPI_COMM_WORLD, &status);

                        if (canGiveSpatialElement &&
                            avgDomainRate > otherAvgDomainRate)
                        {
                            /* If we can give a spatial element to the right neighbor */
                            /* Give up a mesh element by transferring information of the spatial element just before the boundary, which is the other processor's new ghost region */
                            vector<BoundaryChange> transferChanges = srscd->getSpatialElement(srscd->getEndIndex()-1);
                            for (BoundaryChange& bc: transferChanges)
                            {
                                long int data[3] = {bc.objKey, bc.pointIndex, bc.change};
                                MPI_Send(data, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD);
                            }
                            long int data[3] = {};
                            MPI_Send(data, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD); // Send end of transmission message                        }
                            int sinks[LEVELS+1];
                            srscd->getSink(srscd->getEndIndex(), sinks);

                            /* Transfer sinks of the boundary element that we are going to give */
                            MPI_Send(sinks, LEVELS+1, MPI_INT, rightNeighbor, tag, MPI_COMM_WORLD);

                            /* Set up this new domain */
                            srscd->setDomain(srscd->getStartIndex(), srscd->getEndIndex()-1);
                            srscd->examineDomainRate();
                        }
                        else if (otherCanGiveMeSpatialElement &&
                            otherAvgDomainRate > avgDomainRate)
                        {
                            /* If we can receive a mesh element from the right neighbor */
                            /* Receive a mesh element by receiving information of the spatial element one past the other processor's boundary (this new ghost region) */
                            long int recv[3];
                            MPI_Recv(recv, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD, &status);
                            vector<BoundaryChange> receivedTransferChanges;

                            while (recv[0] != 0) // while it's a valid message
                            {
                                receivedTransferChanges.push_back(BoundaryChange(recv[0], recv[1], recv[2]));
                                MPI_Recv(recv, 3, MPI_LONG, rightNeighbor, tag, MPI_COMM_WORLD, &status);
                            }

                            /* Receive sinks of the new boundary element the other processor gave us */
                            int newBoundarySinks[LEVELS+1];
                            MPI_Recv(newBoundarySinks, LEVELS+1, MPI_INT, rightNeighbor, tag, MPI_COMM_WORLD, &status);

                            /* Set up this new domain */
                            srscd->addSpatialElement(srscd->getEndIndex()+2, receivedTransferChanges, srscd->getEndIndex()+1, newBoundarySinks);
                            srscd->setDomain(srscd->getStartIndex(), srscd->getEndIndex()+1);
                            srscd->examineDomainRate();
                        }
                    }
                    else if (threadID % 2 == 1 && leftNeighbor >= 0)
                    {
                        long double otherAvgDomainRate = 0;
                        MPI_Recv(&otherAvgDomainRate, 1, MPI_LONG_DOUBLE, leftNeighbor, tag, MPI_COMM_WORLD, &status);
                        MPI_Send(&avgDomainRate, 1, MPI_LONG_DOUBLE, leftNeighbor, tag, MPI_COMM_WORLD); // Don't give spatial element if we only have one remaining

                        bool canGiveSpatialElement = srscd->getStartIndex() < srscd->getEndIndex();
                        bool otherCanGiveMeSpatialElement;
                        MPI_Recv(&otherCanGiveMeSpatialElement, 1, MPI_CXX_BOOL, leftNeighbor, tag, MPI_COMM_WORLD, &status);
                        MPI_Send(&canGiveSpatialElement, 1, MPI_CXX_BOOL, leftNeighbor, tag, MPI_COMM_WORLD);

                        if (otherCanGiveMeSpatialElement &&
                            avgDomainRate < otherAvgDomainRate)
                        {
                            /* If we can receive a mesh element from the left neighbor */
                            /* Receive a mesh element by receiving information of the spatial element one past the other processor's boundary (this new ghost region) */
                            long int recv[3];
                            MPI_Recv(recv, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD, &status);
                            vector<BoundaryChange> receivedTransferChanges;

                            while (recv[0] != 0) // while it's a valid message
                            {
                                receivedTransferChanges.push_back(BoundaryChange(recv[0], recv[1], recv[2]));
                                MPI_Recv(recv, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD, &status);
                            }

                            /* Receive sinks of the new boundary element the other processor gave us */
                            int newBoundarySinks[LEVELS+1];
                            MPI_Recv(newBoundarySinks, LEVELS+1, MPI_INT, leftNeighbor, tag, MPI_COMM_WORLD, &status);

                            /* Set up this new domain */
                            srscd->addSpatialElement(srscd->getStartIndex()-2, receivedTransferChanges, srscd->getStartIndex()-1, newBoundarySinks);
                            srscd->setDomain(srscd->getStartIndex()-1, srscd->getEndIndex());
                            srscd->examineDomainRate();
                        }
                        else if (canGiveSpatialElement && 
                            avgDomainRate > otherAvgDomainRate)
                        {
                            /* If we can give a spatial element to the left neighbor */
                            /* Give up a mesh element by transferring information of the spatial element just before the boundary, which is the other processor's new ghost region */
                            vector<BoundaryChange> transferChanges = srscd->getSpatialElement(srscd->getStartIndex()+1);
                            for (BoundaryChange& bc: transferChanges)
                            {
                                long int data[3] = {bc.objKey, bc.pointIndex, bc.change};
                                MPI_Send(data, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD);
                            }
                            long int data[3] = {};
                            MPI_Send(data, 3, MPI_LONG, leftNeighbor, tag, MPI_COMM_WORLD); // Send end of transmission message                        }
                            int sinks[LEVELS+1];
                            srscd->getSink(srscd->getEndIndex(), sinks);

                            /* Transfer sinks of the boundary element that we are going to give */
                            MPI_Send(sinks, LEVELS+1, MPI_INT, leftNeighbor, tag, MPI_COMM_WORLD);

                            /* Set up this new domain */
                            srscd->setDomain(srscd->getStartIndex()+1, srscd->getEndIndex());
                            srscd->examineDomainRate();
                        }
                    }
                }
                transferTurn++;
            }
        }
    }
    srscd->drawSpeciesAndReactions(advTime);
    srscd->drawDamage(advTime);
    srscd->writeVacancy();
    srscd->writeSinkFile(advTime, iStep, threadID);
    cout<<"dpa = "<<dpa<<endl;
    cout << "Finished, Bye" << endl;
    return 0;
}
