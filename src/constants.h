//
//  constants.h
//  TestSRSCD
//
//  Created by fx on 2/16/18.
//  Copyright © 2018 fx. All rights reserved.
//
#include <cmath>

#pragma once
/**
 * Definition of constants in the code.
 * All definition are in CAP letters.
 */
// Universal constants:

#define PI 3.1415926
#define KB 8.617e-05       // [ev/K] Boltzmann's constant.

// Material properties:

#define DENSITY 6.30705e+22   // [atoms/cm^3] Atomic density for W.
#define ALATT 3.165e-08     // [cm] Lattice parameter for W.
#define ATOMICVOLUME 0.0158 //[nm^3] Atomic volume of W.
#define BURGER 0.28 //[nm] burger's vector of W.
#define DISLOCATION 1e+8   // [cm^-2] Dislocation density.
//#define DISLOCATION 0.0
#define ODS_R 2.5e-07       // [cm] ODS-particle radius.
//#define ODS_DENSITY 2.6e+17 // [cm^-3] ODS-particle density.
#define ODS_DENSITY 0.0
//#define GRAIN_SIZE 0.01    // [cm] Typical grain size.
#define GRAIN_SIZE 0.0002 //[cm]
#define FOIL_THICKNESS 0.0002 //[cm] Foil thickness from UCSD (2000nm)
#define ELEMENT_THICKNESS 20.0 // [nm] thickness of bulk spatial element
#define SURFACE_THICKNESS 0.544 //[nm] thickness of surface (conrresponds to two monolayers of tungsten)
#define NU0 6.1e+12           // [Hz] Attempt frequency.
#define C_DENSITY 10        // [appm] C-atom density
#define GAMMA 1.0           // Fraction of surface emission.
#define TDE 90              // [eV] Threshold displacement energy for W.
#define HEAT_OF_SOLUTION 1.04   // [eV] Heat of solution of H in W.
#define V1_FORMATION_ENERGY 3.23 // [eV] formation energy of V1
#define SAV_ENERGY 0.95          // [eV] energy that controls SAV rate
#define H_DEPOSITION_ENERGY 200.0 // [eV] energy of an H atom in the incident flux
#define H_FORM_E 1.04           // [eV] Heat of solution of H in W (formation energy of H)
#define V_FORM_E 3.23           // [eV] Formation Energy of 1 Vacancy
#define VH_BIND_E 1.17486       // [eV] Binding Energy of V-H
#define HH_BIND_E 0.02          // [eV] Binding Energy of H-H
// Run parameters:
#define ION               // Irradiation type.
#define TOTAL_TIME 125   // [s] Total simulated time.
#define TOTAL_DPA 2        // Total DPA damage to reach.
#define VOLUME 1.0e-17    // [cm^3] System volume.
#define DIVIDING_AREA 5.0e-12  // [cm^2] Area of the surface in between volume elements.
#define TEMPERATURE 300.0  // [K] System temperature.
//#define RATIO_HE 1.1       // [appm/dpa] He-to-dpa ratio.
#define RATIO_HE 0       // [appm/dpa] He-to-dpa ratio.
#define RATIO_H 0
#define FLUX_H 4.0e16    // [1/(cm^2*s)]
//#define DPA_RATE 0       //When only H exposure is available. no self-damage at all
//#define DPA_RATE 3.55e-6   // [dpa/s] Damage rate.
#define CHANNELS 3         // Irradiation channels used (1:W, 2:He, 3:H,...). the number of different particle insertion(irradiation) process.
#define PSTEPS 30000 // Print data every so many.
#define TSTEPS 50000 // Run these many steps.
#define LEVELS 3
#define EXP10 3 // number of 0s assigned to each group in object's key in hash table
#define POINTS 300       // number of elements: one surface(Point 0), other bulk elements(NO.1,2,3,4,...,99)
// Auxiliary definitions:
enum Reaction { DIFFUSETOF, DIFFUSETOB, SINK, DISSOCIATION, COMBINATION, SAV, NONE, PARTICLE, HE, H, DISSV, DISSH, ERROR};

/*
 ** Reaction Types:
 ** DIFFUSION: diffusion reaction
 ** SINK: absorption
 ** DISSOCIATION: dissociation reaction(one monomer splited out)
 ** COMBINATION: 2nd order reaction, 2 reactants goes to 1 product
 ** SAV: Super Abundant Vacancy mechanism
 ** NONE: SIGNAL-> no reaction been selected in the rateMatrix
 ** PARTICLE: damage[0], ion/neutron insertion
 ** HE: damage[1], He insersion
 ** H: damage[2], H insersion
 ** DISSV : dissociation of vacancy from dislocation
 ** DISSH : dissociation of hydrogren from dislocation
 ** END: SIGNAL-> non damage been selected, rate index in this element goes to the end
 */

const double AVG_ION_EN[POINTS] = {22342, 21202, 16767, 15427, 19102, 23477, 27934, 25724, 21525, 23409, 26167, 30148, 28405, 28622, 33293, 25380, 28119, 35640, 41737, 36054, 38115, 40023, 43367, 40388, 43523, 46082, 31469, 37833, 40849, 39033, 39739, 37562, 35608, 39525, 35038, 33527, 31842, 28732, 28872, 31789, 26140, 25478, 23078, 18575, 16248, 16883, 13112, 11587, 10769, 9441, 7145, 5443, 4832, 4088, 3359, 2229, 2682, 2621, 1144, 727, 630, 408, 514, 163, 207, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
const double H_SATURATION_CONCENTRATION = DENSITY * exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE);
const double SURFACE_VOLUME = VOLUME + (VOLUME / 20. * SURFACE_THICKNESS); // surface element is the thin surface layer + one bulk element


const bool LOG_REACTIONS = false;

/* Configure which features to run */
const bool IRRADIATION_ON = false;
const bool HYDROGEN_ON = true;
const bool COMB_ON = true; // combination reaction
const bool SINK_ON = false; // sink reaction
const bool DISS_ON = true; // dissociation reaction
const bool DIFF_ON = true; // diffusion reaction
const bool SAV_ON = true;  // super abundant vacancy reaction

/* Unit conversion factors */
const double NM_TO_CM = 1.0e-7;

//#define AVG_ION_EN 1.71e+6 // (from TRIM) Average ion energy (in eV) expended on damage from 5MeV Cu.
#define AVG_NEUTRON_EN 40.6 // (from SPECTER) Total damage energy in keV produced by a neutron in ITER.
#define RESTART            // Do restart. */
/* #define RATE_DUMP          // Dump rate spectrum every PSTEPS. */
/* #define DEBUG              // Check cascade damage. */
