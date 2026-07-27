"""
Script to read clusters from one case and distribute those clusters
into a sample of different cross sectional area but same length.

Assumes the clusters are uniformly distributed along the length
of the sample.
"""


import os, re, cv2, time
import math
import random
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from tqdm import tqdm
from textwrap import wrap
from math import floor, ceil
from scipy.signal import butter, filtfilt
from make_speciesfile import combine_species_files


ORIGINAL_POINTS = 252
ORIGINAL_DIVIDING_AREA = 1.0e-10

# Change data files list, times list, and flux for custom use case
# Downsized sample's parameters:
POINTS = 3695                            # num spatial elements in the simulation (1 surface + 100 bulk)
FIRST_EXP_INDEX = 3800
DIVIDING_AREA = 0.4583e-12                    # [cm]
FIRST_BULK_THICKNESS = 6.77                 # [nm]
ELEMENT_THICKNESS = 6.77                   # [nm]
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


combine_species_files()

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

def vacancies_per_cluster(obj_key):
	"""
	obj_key: string of the object's key
	Assumes obj_key is in form of -xxx000mmm, xxx is num vacancies, mmm is num H
	"""
	if int(obj_key) < 0:
		return abs(int(obj_key[:len(obj_key) - 6]))
	else:
		return 0

def interstitials_per_cluster(obj_key):
	"""
	obj_key: string of the object's key
	Assumes obj_key is in form of xxx000mmm, xxx is num sia, mmm is num H
	"""
	if int(obj_key) >= 1000000:
		return abs(int(obj_key[:len(obj_key) - 6]))
	else:
		return 0

with open("species.txt") as f:
	vacancy_c = np.zeros(ORIGINAL_POINTS-1) # don't include vac and sia trapped on the surface element
	sia_c = np.zeros(ORIGINAL_POINTS-1)
	f.readline() #step
	time = float(f.readline().split()[2]) #time
	f.readline() #fluenceH

	vac_cluster_counts = {}
	sia_cluster_counts = {}

	for line_hold in f:
		line_hold = line_hold.split()
		obj_key = int(line_hold[1])

		v_per_cluster = vacancies_per_cluster(line_hold[1])
		if v_per_cluster > 0:
			num_clusters = np.array(line_hold[3:]).astype(float)
			vacancy_c += num_clusters * v_per_cluster
			if np.sum(num_clusters) > 0:
				vac_cluster_counts[obj_key] = np.sum(num_clusters)

		sia_per_cluster = interstitials_per_cluster(line_hold[1])
		if sia_per_cluster > 0:
			num_clusters = np.array(line_hold[3:]).astype(float)
			sia_c += num_clusters * sia_per_cluster
			if np.sum(num_clusters) > 0:
				sia_cluster_counts[obj_key] = np.sum(num_clusters)

print('Vacancies:', np.sum(vacancy_c))
print('SIAs:', np.sum(sia_c))

print()

num_downsampled_vac = round(np.sum(vacancy_c)*DIVIDING_AREA/ORIGINAL_DIVIDING_AREA)
num_downsampled_sia = round(np.sum(sia_c)*DIVIDING_AREA/ORIGINAL_DIVIDING_AREA)

print('Downsampled vacancies:', num_downsampled_vac)
print('Downsampled sias:', num_downsampled_sia)

print()

vac_cluster_counts = dict(sorted(vac_cluster_counts.items(), reverse=True))
vac_cluster_keys = np.array(list(vac_cluster_counts.keys()))
vac_cluster_probs = np.array(list(vac_cluster_counts.values())) / np.sum(list(vac_cluster_counts.values()))
vac_cluster_cdf = np.cumsum(vac_cluster_probs)


sia_cluster_counts = dict(sorted(sia_cluster_counts.items()))
sia_cluster_keys = np.array(list(sia_cluster_counts.keys()))
sia_cluster_probs = np.array(list(sia_cluster_counts.values())) / np.sum(list(sia_cluster_counts.values()))
sia_cluster_cdf = np.cumsum(sia_cluster_probs)

# Sample from cluster distributions to get our desired number of vac and sia for downsized sample
counts = {}
nv = 0

while nv < num_downsampled_vac:
    # Choose which mesh element to insert the cluster to
    point_index = 0
    while point_index == 0 or point_index == 1:
        depth = random.random() * sample_depth_um
        point_index = np.argmin(np.abs(positions - depth))
    
    # Choose the type of cluster this is
    random_num = random.random()
    cdf_idx = np.searchsorted(vac_cluster_cdf, random_num)
    obj_key = vac_cluster_keys[cdf_idx]

    num_vac_in_cluster = abs(int(obj_key / 1e6))
    if num_vac_in_cluster + nv > num_downsampled_vac:
    	num_vac_in_cluster = num_downsampled_vac - nv
    	obj_key = -1000000 * num_vac_in_cluster

    if obj_key not in counts:
        counts[obj_key] = np.zeros(POINTS)

    counts[obj_key][point_index] += 1
    nv += num_vac_in_cluster


nsia = 0

while nsia < num_downsampled_sia:
    # Choose which mesh element to insert the cluster to
    point_index = 0
    while point_index == 0 or point_index == 1:
        depth = random.random() * sample_depth_um
        point_index = np.argmin(np.abs(positions - depth))
    
    # Choose the type of cluster this is
    random_num = random.random()
    cdf_idx = np.searchsorted(sia_cluster_cdf, random_num)
    obj_key = sia_cluster_keys[cdf_idx]

    num_sia_in_cluster = abs(int(obj_key / 1e6))
    if num_sia_in_cluster + nsia > num_downsampled_sia:
    	num_sia_in_cluster = num_downsampled_sia - nsia
    	obj_key = 1000000 * num_sia_in_cluster

    if obj_key not in counts:
        counts[obj_key] = np.zeros(POINTS)

    counts[obj_key][point_index] += 1
    nsia += num_sia_in_cluster

print(list(vac_cluster_keys))
print(list(vac_cluster_probs))

# with open('restart.txt', 'w') as f:
#     f.write("step = 0\n")
#     f.write("time = 0.0\n")
#     f.write("fluenceH = 0.0\n")

#     for obj_key in counts:
#         f.write(f"object {str(obj_key)}")
#         for i in range(len(counts[obj_key])):
#             f.write(f"    {str(int(counts[obj_key][i]))}")
#         f.write("\n")
