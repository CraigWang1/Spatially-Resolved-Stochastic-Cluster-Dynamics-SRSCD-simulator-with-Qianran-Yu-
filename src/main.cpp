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

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);

    int rightNeighbor = threadID + 1;
    int leftNeighbor = threadID - 1;

    SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd */
    int64 theOtherKey = 0;
    Object* hostObject = nullptr;
    Reaction reaction = ERROR;
    int pointIndex = -1;
    int rootThreadID = 0;
    long int iStep = 0;
    double random;
    double advTime = 0.0;
    double dt = 0.0;
    double system_dt = 0.0;
    double accTime = 0.0;
    long double bulkRate = 0.0; /* total rate of the whole bulk */
    double dpa = 0.0;
    double progress = 0.0; /* progress of 50 = 50% done with simulation */
    double eta_min = 0.0;
    double prev_progress = 0.0;
    double prev_eta_min = 0.0;
    double write_time = 0.2;
    double write_increment = 0.2;
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
    double splitIndex = 0;
    int startIndex = threadID * indexIncrement;
    int endIndex = (threadID + 1) * indexIncrement - 1;
    srscd->setDomain(startIndex, endIndex);
    srscd->examineDomainRate();

    double prev_time = MPI_Wtime();
    cout << H_SATURATION_CONCENTRATION * VOLUME << endl;
    
    while(!done)
    {
        srscd->clearNoneReaction();
        long double localDomainRate = srscd->getDomainRate();
        long double maxDomainRate;

        // Each MPI process sends its rate to reduction, every thread collects result
        MPI_Allreduce(&localDomainRate, &maxDomainRate, 1, MPI_LONG_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        srscd->fillNoneReaction(maxDomainRate);
        hostObject = srscd->selectDomainReaction(theOtherKey, reaction, pointIndex);/* choose an event */
        srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime, accTime); /* process event */
        
        if(reaction == Reaction::H || reaction == Reaction::PARTICLE)
        {
            accTime = 0.0;
        }

        // Even processors send to right, odd processors receive from left
        vector<BoundaryChange>* leftBoundaryChanges = srscd->getLeftBoundaryChangeQ();
        vector<BoundaryChange>* rightBoundaryChanges = srscd->getRightBoundaryChangeQ();
        if (threadID % 2 == 0 && rightNeighbor < numThreads)
        {
            for (size_t i = 0; i < leftBoundaryChanges->size(); i++)
            {
                BoundaryChange& bc = leftBoundaryChanges->at(i);
                long int data = {bc.objKey, bc.pointIndex, bc.change};
                MPI_SEND()
            }
            MPI_Send()
        }
        
        bulkRate = srscd->getBulkRate(); /* calculate the bulk rate */
        ++iStep;
        do {
            random = (double)rand() / RAND_MAX;
        } while (random == 0);
        dt = (-1) / bulkRate*log(random);
        accTime += dt;
        advTime += dt;
        dpa = srscd->getTotalDpa();

        if(iStep%PSTEPS == 0)
        {
            system_dt = (MPI_Wtime() - prev_time);
            prev_time = MPI_Wtime();
            st.open("st.txt", ios::app);
            st << (float)prev_time << "  "<< dpa << endl;

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
            srscd->writeFile(advTime, iStep);
            st.close();
        }

        if (IRRADIATION_ON)
        {
            done = (dpa >= TOTAL_DPA);
        }
        else
        {
            done = (advTime >= TOTAL_TIME);
        }
    }
    srscd->drawSpeciesAndReactions(advTime);
    srscd->drawDamage(advTime);
    srscd->writeVacancy();
    srscd->writeSinkFile(advTime, iStep);
    cout<<"dpa = "<<dpa<<endl;
    cout << "Finished, Bye" << endl;
    return 0;
}
