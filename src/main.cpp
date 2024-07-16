// SRSCD -- Spatially Resolved Stochastic Cluster Dynamics
// Q Yu, January 2017.

#include<ctime>
#include<unistd.h>
#include<time.h>
#include <cassert>
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
    double accTime = 0.0;
    long double bulkRate = 0.0; /* total rate of the whole bulk */
    long double maxDomainRate = 0.0;
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
    double split_index = 0;
    omp_set_num_threads(4);
    const int num_threads = omp_get_max_threads();
    double index_increment = (double) POINTS / num_threads;
    fstream st;
    /* check whether to restart*/
    restart(iStep, advTime, master_srscd);
    while (advTime > write_time)
        write_time += write_increment;
    master_srscd->examineRate();
    /* check ended */
    srand(time(0));
    master_srscd->displayAllObject();
    master_srscd->drawSpeciesAndReactions(advTime);
    double prev_time = omp_get_wtime();

    vector<BoundaryChange*>* boundaryChanges[num_threads];
    vector<SCDWrapper*> processors(num_threads);

    cout << H_SATURATION_CONCENTRATION * VOLUME << endl;

    // Run in parallel with multiple processors
    #pragma omp parallel private(theOtherKey, reaction, pointIndex, hostObject, advTime, accTime)
    {
        const int thread_id = omp_get_thread_num();
        const int left_neighbor = thread_id - 1;
        const int right_neighbor = thread_id + 1;

        SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd for single processor */
        restart(iStep, advTime, srscd);

        // Split the volume elements among each processor
        for (int i = 0; i < num_threads; i++)
        {
            if (i == thread_id)
            {
                int startIndex = round(split_index);
                split_index += index_increment;
                int endIndex = round(split_index) - 1;
                srscd->setDomain(startIndex, endIndex);
            }
            #pragma omp barrier
        }

        srscd->examineDomainRate();

        while(!done)
        {
            // Find the greatest domain rate, that all the other processors will adopt (Dunn 2016)
            srscd->clearNoneReaction();
            double localDomainRate = srscd->getDomainRate();

            #pragma omp critical
            {
                if (localDomainRate > maxDomainRate)
                {
                    maxDomainRate = localDomainRate;
                }
            }

            // Wait for all threads to find the global max domain rate
            #pragma omp barrier

            // Calculate global dt
            #pragma omp master
            {
                do {
                    random = (double)rand() / RAND_MAX;
                } while (random == 0);
                dt = (-1) / maxDomainRate*log(random);
            }

            srscd->fillNoneReaction(maxDomainRate); 
            hostObject = srscd->selectDomainReaction(theOtherKey, reaction, pointIndex);/* choose an event */
            srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime, accTime); /* process event */

            if(reaction == Reaction::H || reaction == Reaction::PARTICLE)
            {
                accTime = 0.0;
            }

            // Communicate boundary events.
            boundaryChanges[thread_id] = srscd->getTxBoundaryChangeQueue();

            // Wait for all threads to finish uploading boundary events
            #pragma omp barrier 

            // Have processor read the boundary changes from its neighbors
            if (left_neighbor >= 0 && 
                boundaryChanges[left_neighbor]->size() > 0)
            {
                srscd->implementBoundaryChanges(boundaryChanges[left_neighbor]);
            }
            if (right_neighbor < num_threads &&
                boundaryChanges[right_neighbor]->size() > 0)
            {
                srscd->implementBoundaryChanges(boundaryChanges[right_neighbor]);
            }

            if (iStep%PSTEPS == 0)
            {
                processors[thread_id] = srscd;
                #pragma omp barrier

                #pragma omp master
                {
                    // Keep track of all simulation elements between processors
                    master_srscd->combineVolumeElements(processors);

                    double time = omp_get_wtime();
                    double system_dt = (time - prev_time);
                    prev_time = time;
                    st.open("st.txt", ios::app);
                    st << (float)prev_time/CLOCKS_PER_SEC << "  "<< dpa << endl;

                    /*
                    cout << "\nt = " << advTime << endl;
                    cout <<"iStep = "<< iStep << endl;
                    cout<<"dt= " << dt <<endl;
                    cout<<"BulkRate = "<<bulkRate<<endl;
                    */

                    // Chose between dpa or time to calculate progress
                    if (dpa != 0)
                    {
                        progress = (dpa / TOTAL_DPA) * 100.;
                    }
                    else
                    {
                        progress = (advTime / TOTAL_TIME) * 100.;
                    }

                    // Initialize the first eta estimate, or update the previous one
                    if (prev_progress != 0)
                    {
                        eta_min = (100 - progress) / ((progress - prev_progress) / system_dt) / 60.;
                        
                        if (prev_eta_min != 0)
                        {
                            eta_min = prev_eta_min + 0.1 * (eta_min - prev_eta_min);

                            // Print progress bar
                            cout << "[";
                            int barWidth = 70;
                            int pos = barWidth * (progress/100.);
                            for (int i = 0; i < barWidth; i++)
                            {
                                if (i < pos) cout << "=";
                                else if (i == pos) cout << ">";
                                else cout << " ";
                            }

                            // Print a set amount of digits
                            int numDigits = 7;
                            int magnitude = 0;
                            while ((int)(advTime / pow(10, magnitude)))
                            {
                                magnitude++;
                            }
                            int numSigFigs = 2;
                            cout << "] " << std::fixed << std::setprecision(2) << progress << "%";
                            cout << "   eta: " << (long int) round(eta_min) << " min";
                            cout << "   time: " << std::fixed << std::setprecision(numDigits - magnitude) << advTime << " s            \r";
                            cout.flush();
                        }
                    }

                    prev_eta_min = eta_min;
                    prev_progress = progress;

                    // Keep track of the volume elements on all processors combined

                    // master_srscd->drawSpeciesAndReactions(advTime);
                    // master_srscd->drawDamage(advTime);
                    master_srscd->writeFile(advTime, iStep);
                    //master_srscd->writeVacancy();
                    //master_srscd->drawHD(advTime);
                    // master_srscd->countDefectNumber(2, "H");

                    st.close();

                    if (advTime / write_time >= 1)
                    {
                        // master_srscd->writeFile(advTime, iStep);
                        write_time += write_increment;
                    }
                }
            }

            #pragma omp barrier

            srscd->clearTxBoundaryChangeQueue();

            #pragma omp master
            {
                ++iStep;
                dpa = 0;
                maxDomainRate = 0.0;
            }

            #pragma omp barrier

            #pragma omp atomic
            dpa += srscd->getDomainDpa();

            accTime += dt;
            advTime += dt;
        }

        // Keep track of the combined simulation volume for logging
        processors[thread_id] = srscd;
    }
    master_srscd->drawSpeciesAndReactions(advTime);
    master_srscd->drawDamage(advTime);
    master_srscd->writeVacancy();
    master_srscd->writeSinkFile(advTime, iStep);
    cout<<"dpa = "<<dpa<<endl;
    cout << "Finished, Bye" << endl;
    return 0;
}