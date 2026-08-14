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

TOTAL_NUM_H = 23593

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

obj_keys = [-1000000, -1000001, -1000002, -1000003, -1000004, -1000005, -1000006, -2000003, -2000004, -2000005, -2000006, -2000007, -2000008, -3000006, -3000007, -3000008, -3000009, -3000010, -4000009, -4000010, -4000011, -4000012, -5000012, -5000013, -5000014, -6000016, -7000015, -7000017, -9000021, -10000026, -11000023, -12000027, -15000029, -17000032, -20000038, -20000039, -32000050, -58000079, -75000097, -218000183, -225000193, -228000198, -347000244, -639000348, -680000362, -2302000690, -3136000820, -6258001156]  
num_H_per_key = [0, 1, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 8, 6, 7, 8, 9, 10, 9, 10, 11, 12, 12, 13, 14, 16, 15, 17, 21, 26, 23, 27, 29, 32, 38, 39, 50, 79, 97, 183, 193, 198, 244, 348, 362, 690, 820, 1156]  
num_H_cdf = [0.12389380530973451, 0.25073746312684364, 0.40117994100294985, 0.504424778761062, 0.5545722713864307, 0.5575221238938054, 0.5604719764011801, 0.5663716814159293, 0.5840707964601771, 0.6637168141592922, 0.7728613569321535, 0.7905604719764013, 0.793510324483776, 0.7964601769911507, 0.8171091445427731, 0.8613569321533925, 0.873156342182891, 0.8761061946902656, 0.8879056047197642, 0.8997050147492627, 0.9144542772861358, 0.9174041297935105, 0.9233038348082597, 0.9262536873156344, 0.9292035398230091, 0.9321533923303837, 0.9351032448377584, 0.9410029498525077, 0.9439528023598823, 0.946902654867257, 0.9498525073746317, 0.9528023598820063, 0.955752212389381, 0.9587020648967557, 0.9616519174041304, 0.964601769911505, 0.9675516224188797, 0.9705014749262544, 0.973451327433629, 0.9764011799410037, 0.9793510324483784, 0.9823008849557531, 0.9852507374631277, 0.9882005899705024, 0.9911504424778771, 0.9941002949852518, 0.9970501474926264, 1.000000000000001] 
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
