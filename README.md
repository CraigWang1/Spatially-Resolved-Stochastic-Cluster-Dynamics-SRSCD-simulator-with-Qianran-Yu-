# Spatially-Resolved-Stochastic-Cluster-Dynamics-SRSCD-simulator-with-Qianran-Yu-
The SRSCD simulation package contains a computer code written in C++. It is an alternative mean field rate theory model that dynamically simulates microstructure (as clusters) evolutions of multi-species systems in solid materials in 1-D space. We have used the code to study Zr-H and W-H system.

****Developed by****

Qianran Yu

Research projects on:
1. Zr-Hydride nucleation and growth under LWR fission condition
2. Hydrogen retention in heavy ion irradiated W materials

Please send an email to the following address for more information:

Qianran Yu (yuqianran0709@gmail.com)

Jaime Marian (jmarian@g.ucla.edu)

Craig Wang (craigwang@g.ucla.edu)

****Introduction****

The SRSCD is a stochastic variant of the mean field rate theory method that is typically developed to simulate microstructure evolutions during permeation of long-term migrating elements, such as hydrogen, into metallic materials. Please see the following journal papers for detail:

[1] Qianran Yu, Michael Reyes, Nachiket Shah and Jaime Marian, "Kinetic Model of Incipient Hydride Formation in Zr Clad under Dynamic Oxide Growth Conditions", Materials 13(5), 1088 (2020). (link: https://www.mdpi.com/1996-1944/13/5/1088)

[2] Qianran Yu, Michael J. Simmonds, Russ. Doerner, George R. Tynan, Li Yang, Brian D. Wirth and Jaime Marian, “Understanding hydrogen retention in damaged tungsten using experimentally-guided models of complex multispecies evolution”, Nuclear Fusion 60, 096003 (2020). (link: https://iopscience.iop.org/article/10.1088/1741-4326/ab9b3c)

****How to use****

Before use: please install gcc/g++ for the newest version (c++11 or newer)

Set up parameters, check equations in relevant function, copy/create input files in "src" folder, and then simply type "make". An
executable file named "scdexe" will be generated. Run the simulations by using command "./scdexe".

****Physical mechanisms****

1. (0th) External particle insertion :
   - Once a non-hydrogen particle insertion event is selected, a set of random pka energies are selected based on cpdf (obtained from SRIM) until the sum of the energies reaches total incident energy that is expended on lattice damage. Then the code performs collision cascade damage by generating certain number and/or size of defect clusters with statistical variations. Such statistical variation as well as the fraction of SIA/V clusters being generated comes from MD studies.  
2. (1st) Monomer dissociation:
   - single SIA/V atom departs from a cluster.
3. (1st) Defect absorption by dislocations (or sinks)
   - Defects come and get trapped into sinks (such as dislocations, grain boundaries, precipitates etc.)
4. (2nd) Binary combination:
   - Two clusters combine and become a larger cluster, or get annhilated. At least one of the reactants need to be mobile.
5. Long term migration (diffusion):
   - Clusters' long-term migration by Fick's law. 
6. Hydrogen desorption from surface:
   - Hydrogen adsorbing onto the surface of the material and leaving.
7. Super Abundant Vacancy (SAV):
   - Overpressurized VH clusters, or oversaturated H atoms, can eject a tungsten atom to create another vacancy.

****Program structure****

A. In the "src" folder:
- Object.cpp/Object.h: store information of species.
- constant.h: store parameters. 
- Damage.cpp/Damage.h: store external particle insertion rates.
- CascadeDamage.cpp/CascadeDamage.h: functions that process cascade damage.
- cpdf.cpp/cpdf.h: sample pka energies using cpdf function.
- rvgs.cpp/rvgs.h: store statistical functions.
- gnuplot_i.h: a library for plotting figures by gnuplot.
- OneLine.cpp/OneLine.h: compuate 1st, 2nd order reaction and diffusion reaction rates. This class handles the information of one species in one spatial element 
- Bundle.cpp/Bundle.h: a class that links information of one species in all spatial elements together.
- SCDWrapper.cpp/SCDWrapper.h: the class that handles the whole rate matrix. It includes ways to select and process events and functions to update rates. There are also output functions. 
- main.cpp: main function. Also is where Temperature parameter is set.
- makefile

B. example_input.zip:
This folder includes the input files to run the example of 3.4 MeV Cu ion irradiation on tungsten materials

****Input files****

The users need to manually create input files to run specific cases. These input files are:
- damage.txt: see "damage.cpp/damage()"
- cpdf#.txt : the number of cpdf files should be equal to the declared number of points (constant.h), see "cpdf.cpp/cpdf()"
- restart.txt : this file is needed when you want to continue the simulation from a previous result, see "SCDWrapper.cpp/restart()"
- sink.txt: store sink absorption information from previous result, see "SCDWrapper.cpp/restart()"

****Output files****

After running the simulation, the following files will appear in the src folder:
- species0.txt: How much of each species is stored at each spatial element. Each object is represented by a 9 digit code. The first 3 digits represent the number of W interstitials (if positive) or vacancies (if negative). The next 3 digits represent the number of He. The final 3 digits represent the number of H in the cluster. If you wish to restart from a checkpoint for the next simulation, delete the "startIndex = ..." and "endIndex = ..." lines inside the file and rename it to "restart.txt".
- sink0.txt: How much of each species is stored in sinks (eg. grain boundaries, dislocations) at each element. The first column is number of W vacancies stored in dislocations, the second column is number of W interstitials stored in dislocations, the third column is number of He stored in dislocations, and the fourth column is number of H stored in dislocations. The 5th, 6th, 7th, and 8th columns follow the same pattern, but they are stored in grain boundaries instead. Each row represents the next spatial element. If you wish to restart from a checkpoint for the next simulation, delete the "startIndex = ..." and "endIndex = ..." lines inside the file and rename it to "sink.txt".
- selectReaction.txt: Log of which reactions were selected at each simulation step.
- Reactions.txt: Log of reaction performed at each time step.

****Additional remarks****

The current code is currently under development to simulate hydrogen deposition into W and thermal desorption processes. It can also be used to simulate irradiation damage irradiation damage, but that is also unfinished. Each case includes different/additional functions. The user is welcomed to contact Qianran or Craig for specific version of the SRSCD code.
