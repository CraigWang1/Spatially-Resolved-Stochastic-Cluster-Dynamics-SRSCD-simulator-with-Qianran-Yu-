// SRSCD -- Spatially Resolved Stochastic Cluster Dynamics
// Q Yu, January 2017.

#include<ctime>
#include<unistd.h>
#include<time.h>
#include"SCDWrapper.h"

int main() {
    //fork();
    cout << H_SATURATION_CONCENTRATION * VOLUME << endl;
    SCDWrapper* srscd = new SCDWrapper(); /* establish spatially resolved scd */
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
    int inputH = 0;
    double totalDPA = 0.2;
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
    srscd->getAndExamineRate();
    /* check ended */
    srand(time(0));
    /*display damage*/
    // srscd->displayDamage();
    srscd->displayAllObject();
    srscd->drawSpeciesAndReactions(advTime);
    clock_t prev_time = clock();
    //srscd->drawHD(advTime);
    while(!done)
    {
        hostObject = srscd->selectReaction(theOtherKey, reaction, pointIndex);/* choose an event */
        srscd->processEvent(reaction, hostObject, pointIndex, theOtherKey, advTime, accTime); /* process event */
        if(reaction == 6)
        {
            accTime = 0.0;
        }
        
        if(iStep%PSTEPS == 0)
        {
            system_dt = (clock() - prev_time) / (double)CLOCKS_PER_SEC;
            prev_time = clock();
            st.open("st.txt", ios::app);
            st << (float)prev_time/CLOCKS_PER_SEC << "  "<< dpa << endl;
            
            /*
            cout << "\nt = " << advTime << endl;
            cout <<"iStep = "<< iStep << endl;
            cout<<"dt= " << dt <<endl;
            cout<<"BulkRate = "<<bulkRate<<endl;
            */

            // Chose between dpa or time to calculate progress
            if (IRRADIATION_ON)
            {
                progress = (dpa / totalDPA) * 100.;
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
             
            // srscd->drawSpeciesAndReactions(advTime);
            // srscd->drawDamage(advTime);
            srscd->writeFile(advTime, iStep);
            // srscd->writeVacancy();
            // srscd->drawHD(advTime);
            st.close();
        }
        // if(iStep % 5 == 0)
        if(iStep % PSTEPS == 0)
        {
            // srscd->drawHD(advTime);
            // srscd->countDefectNumber(2, "H");
        }

        if (advTime / write_time >= 1)
        {
            // srscd->writeFile(advTime, iStep);
            write_time += write_increment;
        }

        bulkRate = srscd->getAndExamineRate(); /* calculate the bulk rate */
        ++iStep;
        do {
            random = (double)rand() / RAND_MAX;
        } while (random == 0);
        dt = (-1) / bulkRate*log(random);
        accTime += dt;
        advTime += dt;
        dpa = srscd->getTotalDpa();

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
