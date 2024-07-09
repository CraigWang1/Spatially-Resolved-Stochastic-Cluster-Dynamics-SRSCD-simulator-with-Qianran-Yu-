// SRSCD -- Spatially Resolved Stochastic Cluster Dynamics
// Q Yu, January 2017.

#include<ctime>
#include<unistd.h>
#include<time.h>
#include <cassert>
#include <omp.h>
#include"SCDWrapper.h"

int main() 
{
    SCDWrapper* master_srscd = new SCDWrapper(); /* establish spatially resolved scd */
    int64 theOtherKey = 0;
    Object* hostObject = nullptr;
    Reaction reaction = ERROR;
    int pointIndex = -1;
    long int iStep = 0;
    double random;
    double advTime = 0.0;
    double dt = 0.0;
    double system_dt = 0.0;
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
    double splitIdx = 0; /* running in parallel parameters */
    // const int numThreads = omp_get_max_threads();
    const int numThreads = 4;
    omp_set_num_threads(numThreads);
    double idxIncrement = (double) POINTS / numThreads;
    fstream st;
    /* check whether to restart*/
    restart(iStep, advTime, master_srscd);
    while (advTime > write_time)
        write_time += write_increment;
    /* check ended */
    srand(time(0));
    master_srscd->displayAllObject();
    master_srscd->drawSpeciesAndReactions(advTime);
    double prev_time = omp_get_wtime();
    vector<list<BoundaryChange*>*> leftBoundaryChanges(numThreads, nullptr);
    vector<list<BoundaryChange*>*> rightBoundaryChanges(numThreads, nullptr);
    vector<SCDWrapper*> processors(numThreads);
    long double globalMaxAvgRate = 0;

    cout << H_SATURATION_CONCENTRATION * VOLUME << endl;

    // Run in parallel with multiple processors
    #pragma omp parallel private(theOtherKey, reaction, pointIndex, hostObject, advTime)
    {
        vector<BoundaryChange*> neighborChanges;
        const int thread_id = omp_get_thread_num();

        SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd for single processor */
        processors[thread_id] = srscd;
        restart(iStep, advTime, srscd);

        // Split the volume elements among each processor
        for (int i = 0; i < numThreads; i++)
        {
            if (i == thread_id)
            {
                int startIndex = round(splitIdx);
                splitIdx += idxIncrement;
                int endIndex = round(splitIdx) - 1;
                srscd->setDomain(startIndex, endIndex);
            }
            #pragma omp barrier
        }

        srscd->examineRate();

        while(!done)
        {
            long double localMaxAvgRate = srscd->findMaxAvgSectorRate();
            globalMaxAvgRate = max(globalMaxAvgRate, localMaxAvgRate);
            #pragma omp barrier

            #pragma omp single
            {
                if (globalMaxAvgRate > 0)
                    dt = NUM_REACTIONS_PER_CYCLE / globalMaxAvgRate;
                else
                    dt = TOTAL_TIME - advTime;
                dt = min(dt, TOTAL_TIME - advTime);
                globalMaxAvgRate = 0.0;
            }
            #pragma omp barrier

            for (int sector = 0; sector < NUM_SECTORS; sector++)
            {
                double local_advTime = 0;
                while (local_advTime < dt)
                {
                    long double sectorRate = srscd->getSectorRate(sector); /* calculate the bulk rate */
                    do {
                        random = (double)rand() / RAND_MAX;
                    } while (random == 0);
                    double local_dt = (-1) / sectorRate * log(random);

                    if (local_advTime + local_dt > dt)
                        break;

                    hostObject = srscd->selectSectorReaction(theOtherKey, reaction, pointIndex, sector);
                    srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime);
                    local_advTime += local_dt;
                }

                // Communicate boundary changes
                if (sector == 0)
                    leftBoundaryChanges[thread_id] = srscd->getLeftBoundaryChanges();
                else if (sector == NUM_SECTORS - 1)
                    rightBoundaryChanges[thread_id] = srscd->getRightBoundaryChanges();
                #pragma omp barrier

                if (thread_id > 0)
                    srscd->processBoundaryChanges(rightBoundaryChanges[thread_id-1]);
                if (thread_id < numThreads-1)
                    srscd->processBoundaryChanges(leftBoundaryChanges[thread_id+1]);
            }

            advTime += dt;
            dpa = srscd->getTotalDpa();

            #pragma omp single
            {
                ++iStep;
                if(iStep%PSTEPS == 0)
                {
                    master_srscd->combineProcessors(processors);
                    system_dt = omp_get_wtime() - prev_time;
                    prev_time = omp_get_wtime();
                    st.open("st.txt", ios::app);
                    st << prev_time << "  "<< dpa << endl;

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
                    master_srscd->writeFile(advTime, iStep);
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
            #pragma omp barrier
        }
    }
    master_srscd->drawSpeciesAndReactions(advTime);
    master_srscd->drawDamage(advTime);
    master_srscd->writeVacancy();
    master_srscd->writeSinkFile(advTime, iStep);
    cout<<"dpa = "<<dpa<<endl;
    cout << "Finished, Bye" << endl;
    return 0;
}
