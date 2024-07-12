import math

POINTS = 100                             # num spatial elements in the simulation (1 surface + 100 bulk)
VOLUME = 1e-17                           # volume of a spatial element [cm^3]
SURFACE_THICKNESS = 0.544                # [nm]
SURFACE_VOLUME = VOLUME / 20 * SURFACE_THICKNESS + VOLUME # [cm^3]
DENSITY = 6.30705e+22                      # [atoms/cm^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             # [ev/K] Boltzmann's constant.
TEMPERATURE = 300

Z = 8                                    # coordination number in a BCC lattice
H_FORM_E = 1.04                          # [eV] Heat of solution of H in W (formation energy of H)
V_FORM_E = 3.23                            # [eV] Formation Energy of 1 Vacancy
VH_BIND_E = 1.17486                           # [eV] Binding Energy of V-H
HH_BIND_E = 0.02                           # [eV] Binding Energy of H-H

monomer_c = DENSITY * math.exp(-H_FORM_E / KB / TEMPERATURE)
dimer_c = DENSITY * Z * math.exp(-(2 * H_FORM_E - HH_BIND_E) / KB / TEMPERATURE)
vhpair_c = DENSITY * Z * math.exp(-(H_FORM_E + V_FORM_E - VH_BIND_E) / KB / TEMPERATURE)
total_c = monomer_c + dimer_c + vhpair_c

with open('restart.txt', 'w') as f:
    f.write("step = 0\n")
    f.write("time = 0.0\n")
    f.write("fluenceH = 0.0\n")

    num_surface = math.floor(SURFACE_VOLUME * total_c)
    f.write(f"object 1    {num_surface}")

    num_bulk = math.floor(VOLUME * total_c)
    for i in range(100):
        f.write(f"    {num_bulk}")
    f.write("\n")