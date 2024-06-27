import math

POINTS = 100                             # num spatial elements in the simulation (1 surface + 100 bulk)
VOLUME = 1e-17                           # volume of a spatial element [cm^3]
SURFACE_THICKNESS = 0.544                # [nm]
SURFACE_VOLUME = VOLUME / 20 * SURFACE_THICKNESS + VOLUME # [cm^3]
DENSITY = 6.30705e+22                      # [atoms/cm^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             # [ev/K] Boltzmann's constant.
TEMPERATURE = 1000
H_SATURATION_CONCENTRATION = DENSITY * math.exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE)


with open('restart.txt', 'w') as f:
    f.write("step = 0\n")
    f.write("time = 0.0\n")
    f.write("fluenceH = 0.0\n")

    num_surface = math.floor(SURFACE_VOLUME * H_SATURATION_CONCENTRATION)
    f.write(f"object 1    {num_surface}")

    num_bulk = math.floor(VOLUME * H_SATURATION_CONCENTRATION)
    for i in range(100):
        f.write(f"    {num_bulk}")
    f.write("\n")