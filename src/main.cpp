// SRSCD -- Spatially Resolved Stochastic Cluster Dynamics
// Q Yu, January 2017.

#include<mpi.h>
#include<ctime>
#include<unistd.h>
#include<time.h>
#include <cassert>
#include"SCDWrapper.h"

double TEMPERATURE = 864;  // [K], this is extern so all files have access to this
const double startingTemp = TEMPERATURE;

int main(int argc, char** argv) 
{
    int threadID, numThreads;
    int rootThreadID = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &threadID);

    SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd */
    int64 theOtherKey = 0;
    Object* hostObject = nullptr;
    Reaction reaction = ERROR;
    int pointIndex = -1;
    long int iStep = 0;
    double random;
    double advTime = 0.0;
    double prevRecalculateTDSRatesTime = 0;
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
    int barWidth, pos; /* progress bar parameters */
    int numDigits = 7;
    int magnitude = 0;
    bool done = false;
    fstream st;
    /* check whether to restart*/
    restart(iStep, advTime, srscd);
    srscd->writeFile(advTime, iStep, threadID);
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
    
    while(!done)
    {
        if (abs(TEMP_INCREASE_RATE) > 0 && advTime - prevRecalculateTDSRatesTime > 1) // if doing thermal desorption, recalculate thermal properties once per second (not all the time) to reduce computational burden (stepwise temperature increase)
        {
            TEMPERATURE = startingTemp + advTime * TEMP_INCREASE_RATE;
            srscd->recalculateAllRates();
            prevRecalculateTDSRatesTime = advTime;
        }
        long double localDomainRate = srscd->getDomainRate();
        long double maxDomainRate = localDomainRate;   // parallelism not implemented fully yet

        if (threadID == rootThreadID)
        {
            do {
                random = (double)rand() / RAND_MAX;
            } while (random == 0);
            dt = (-1) / maxDomainRate*log(random);
            accTime += dt;
            advTime += dt;

            // if (IRRADIATION_ON)
            // {
                // done = (dpa >= TOTAL_DPA);
            // }
            // else
            // {
                done = (advTime >= TOTAL_TIME);
            // }
            if (done)
            {
                MPI_Bcast(&done, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
            }
        }

        srscd->fillNoneReaction(maxDomainRate);
        hostObject = srscd->selectDomainReaction(theOtherKey, reaction, pointIndex);/* choose an event */

        /*
        string reactions[] = {"diffF", "diffB", "sinkDislocScrew", "sinkDislocEdge", "sinkGrain", "diss", "comb", "sav", "recombER", "recombLH", "none", "particle", "HE", "H", "dissVDislocScrew", "dissVDislocEdge", "dissVGrain", "dissHDislocScrew", "dissHDislocEdge", "dissHGrain", "error"};
        if (hostObject != nullptr)
            cout << hostObject->getKey() << " ";
        cout << reactions[reaction] << " at pt " << pointIndex << " with other obj " << theOtherKey;
        cout << endl;
        */

        srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime, accTime); /* process event */

        if(reaction == Reaction::H || reaction == Reaction::PARTICLE)
        {
            accTime = 0.0;
        }

        if (IRRADIATION_ON)
        {
            double localDpa = srscd->getTotalDpa();
            MPI_Reduce(&localDpa, &dpa, 1, MPI_DOUBLE, MPI_SUM, rootThreadID, MPI_COMM_WORLD);
        }

        ++iStep;
        if(iStep%PSTEPS == 0 || done)
        {
            srscd->writeFile(advTime, iStep, threadID);

            if (threadID == rootThreadID)
            {
                system_dt = (MPI_Wtime() - prev_time);
                prev_time = MPI_Wtime();

                // Chose between dpa or time to calculate progress
                // if (IRRADIATION_ON)
                // {
                    // progress = (dpa / TOTAL_DPA) * 100.;
                // }
                // else
                // {
                    progress = (advTime / TOTAL_TIME) * 100.;
                // }

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
            if (done)
            {
                cout << "\n" << endl;
            }
        }
    }
    srscd->writeFile(advTime, iStep, threadID);
    srscd->drawSpeciesAndReactions(advTime);
    srscd->drawDamage(advTime);
    srscd->writeVacancy();

    MPI_Barrier(MPI_COMM_WORLD);  // Make sure all other threads are done printing if applicable
    if (threadID == rootThreadID)
    {
        cout << "dpa = " << dpa << endl;
        cout << "Finished, Bye" << endl;
    }
    MPI_Finalize();

    return 0;
}
