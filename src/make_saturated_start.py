import math
import random
import numpy as np

# Change data files list, times list, and flux for custom use case
POINTS = 832                            # num spatial elements in the simulation (1 surface + 100 bulk)
FIRST_EXP_INDEX = 739
DIVIDING_AREA = 0.4583e-12                    # [cm]
FIRST_BULK_THICKNESS = 6.77                 # [nm]
ELEMENT_THICKNESS = 6.77                   # [nm]
TOTAL_TIME = 7692.3         # [s]
BACK_DESORB = False

NM_TO_CM = 1e-7
NM_TO_UM = 1e-3
CM_TO_UM = 1e4
VOLUME = DIVIDING_AREA * ELEMENT_THICKNESS * NM_TO_CM  
SUBSURFACE_THICKNESS = 0.544                # [nm]
SURFACE_INDEX = 0
SUBSURFACE_INDEX = 1
FIRST_BULK_INDEX = 2
BACK_SUBSURFACE_INDEX = POINTS - 2
BACK_SURFACE_INDEX = POINTS - 1
EXP_LENGTH_MULT = 1.1


DENSITY = 6.30705e+22                      # [atoms/cm^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             # [ev/K] Boltzmann's constant.
TEMPERATURE = 300
H_SATURATION_CONCENTRATION = DENSITY * math.exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE) / DENSITY * 100
dpi = 100

TOTAL_NUM_H = 4248

def volumeAtIndex(i):
    """
    Returns the volume (cm^3) of volume element i
    """
    return DIVIDING_AREA * length(i)

def length(i):
    """
    Returns the length (cm) of volume element i
    """
    if i == SURFACE_INDEX or (i == BACK_SURFACE_INDEX and BACK_DESORB):
        return 0
    if i == SUBSURFACE_INDEX or (i == BACK_SUBSURFACE_INDEX and BACK_DESORB):
        return SUBSURFACE_THICKNESS * NM_TO_CM
    if i == FIRST_BULK_INDEX:
        return FIRST_BULK_THICKNESS * NM_TO_CM
    if i < FIRST_EXP_INDEX:
        return ELEMENT_THICKNESS * NM_TO_CM
    return ELEMENT_THICKNESS * EXP_LENGTH_MULT ** (i - FIRST_EXP_INDEX) * NM_TO_CM

positions = []
for i in range(POINTS):
    if i == 0:
        positions.append(0)  # surface element
    else:
        positions.append( positions[-1] + (length(i) + length(i-1))/2 * CM_TO_UM )

positions = np.array(positions)
sample_depth_um = positions[-1] + length(POINTS-1)/2*CM_TO_UM
# print(sample_depth_um)

obj_keys = [
    -1000000,
    -1000001,
    -1000002,
    -1000003,
    -1000004,
    -2000005,
]
num_H_per_key = [
    0,
    1,
    2,
    3,
    4,
    5,
]
num_H_cdf = [0.25, 0.59, 0.79, 0.89, 0.96, 1]
num_H = 0

counts = {}

while num_H < TOTAL_NUM_H:
    # Choose which mesh element to insert the cluster to
    point_index = 0
    while point_index == 0 or point_index == 1:
        depth = random.random() * sample_depth_um
        point_index = np.argmin(np.abs(positions - depth))
    
    # Choose the type of cluster this is
    random_num = random.random()
    cdf_idx = np.searchsorted(num_H_cdf, random_num)
    obj_key = obj_keys[cdf_idx]

    if obj_key not in counts:
        counts[obj_key] = np.zeros(POINTS)

    counts[obj_key][point_index] += 1
    num_H += num_H_per_key[cdf_idx]

with open('restart.txt', 'w') as f:
    f.write("step = 0\n")
    f.write("time = 0.0\n")
    f.write("fluenceH = 0.0\n")

    for obj_key in counts:
        f.write(f"object {str(obj_key)}")
        for i in range(len(counts[obj_key])):
            f.write(f"    {str(int(counts[obj_key][i]))}")
        f.write("\n")
