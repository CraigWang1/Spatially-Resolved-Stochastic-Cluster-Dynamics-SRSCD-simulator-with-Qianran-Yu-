import os, re, cv2, time
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import addcopyfighandler
from tqdm import tqdm
from textwrap import wrap
from math import floor, ceil
from scipy.signal import butter, filtfilt
from make_speciesfile import combine_species_files

combine_species_files()

plt.rcParams.update({'font.size': 14})

# Change data files list, times list, and flux for custom use case
POINTS = 252                            # num spatial elements in the simulation (1 surface + 100 bulk)
FIRST_EXP_INDEX = 700
DIVIDING_AREA = 0.5141e-12                    # [cm]
FIRST_BULK_THICKNESS = 7.17                 # [nm]
ELEMENT_THICKNESS = 7.17                   # [nm]
TOTAL_TIME = 10000         # [s]
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

def getConcentration(x, t):
	""" 
	Returns concentration (# hydrogen / m^3) 
	given input of position (x, in meters) and time (t, in seconds) (analytical solution)
	https://www.desmos.com/calculator/crbfp7vgdh
	"""
	prefactor = flux * l  / D
	term1 = D * t / l ** 2
	term2 = (3 * (x - l) ** 2 - l ** 2) / (6 * l ** 2)
	term3 = 0
	for n in range(1, 101):
		term3 += ((-1) ** n / n ** 2) * math.exp(-D * n ** 2 * math.pi ** 2 * t / l ** 2) * math.cos(n * math.pi * (x - l) / l)
	term3 *= -2 / math.pi ** 2
	return prefactor * (term1 + term2 + term3)
	
def hydrogen_per_cluster(obj_key):
	"""
	obj_key: string of the object's key
	Assumes obj_key is in form of xxx, xxx is num H
	"""
	if len(obj_key) >= 6:
		return abs(int(obj_key[len(obj_key)-6:]))
	else:
		return int(obj_key)

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

# out = cv2.VideoWriter('output.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 40.0, (640,480))

# plt.figure(figsize=(640/dpi, 480/dpi), dpi=dpi)

fluences = [
	"1e24",
	# "5e23",
	# "1e24"
]

# Plot Experiment
for i in range(len(fluences)):
	with open(f"/home/craig/research/experiment_retention_383K/fluence_{fluences[i]}.txt") as f:
		# f.readline() # column titles
		experiment_positions = []
		concentrations = []
		for line_hold in f:
			line_hold = line_hold.split(", ")
			depth_micrometer = float(line_hold[0])
			# if depth_micrometer > 2:
				# break
			concentration_at = float(line_hold[1]) # at. % units
			experiment_positions.append(depth_micrometer)
			concentrations.append(concentration_at)
		# plt.plot(positions, concentrations, label=f"Experiment {fluences[i]}")
		# plt.plot(experiment_positions, concentrations, label="Experiment", color='r')

# Plot Simulation
# with open("/home/craig/Downloads/Spatially-Resolved-Stochastic-Cluster-Dynamics-SRSCD-simulator-with-Qianran-Yu-/src/species.txt") as f:
with open("species.txt") as f:
	positions = []
	for i in range(POINTS):
		if i == 0:
			positions.append(0)  # surface element
		else:
			positions.append( positions[-1] + (length(i) + length(i-1))/2 * CM_TO_UM )

	print('Sample length:', round(positions[-1] + length(POINTS-1)/2*CM_TO_UM, 2), 'um')

	trapped_hydrogen_c = np.zeros(POINTS)
	free_hydrogen_c = np.zeros(POINTS)
	vacancy_c = np.zeros(POINTS)
	sia_c = np.zeros(POINTS)
	plot_h = False
	plot_v = False
	f.readline() #step
	time = float(f.readline().split()[2]) #time
	print("time (s):", round(time, 2))
	f.readline() #fluenceH
	for line_hold in f:
		line_hold = line_hold.split()
		obj_key = int(line_hold[1])
		h_per_cluster = hydrogen_per_cluster(line_hold[1])
		# if h_per_cluster > 0:    # for plotting all H
		if obj_key > 0 and obj_key < 1000000:  # free H
			free_hydrogen_c += np.array(line_hold[2:]).astype(float) * h_per_cluster
			plot_h = True
		elif (obj_key < 0 or obj_key > 1000000) and h_per_cluster > 0: # trapped H
			trapped_hydrogen_c += np.array(line_hold[2:]).astype(float) * h_per_cluster
			plot_h = True

		v_per_cluster = vacancies_per_cluster(line_hold[1])
		if v_per_cluster > 0:
			vacancy_c += np.array(line_hold[2:]).astype(float) * v_per_cluster
			plot_v = True

		sia_per_cluster = interstitials_per_cluster(line_hold[1])
		if sia_per_cluster > 0:
			sia_c += np.array(line_hold[2:]).astype(float) * sia_per_cluster

with open("sink0.txt") as f:
	f.readline()
	f.readline()
	numH = []
	for line_hold in f:
		line_hold = line_hold.split()
		numH.append(int(line_hold[3]) + int(line_hold[7]) + int(line_hold[11]))
	trapped_hydrogen_c += np.array(numH).astype(float)
	# print(sum(numH)/np.sum(trapped_hydrogen_c))
print("Retained fluence [m^-2]:", np.sum(trapped_hydrogen_c/DIVIDING_AREA*1e4))
print("Projected fluence [m^-2]:", np.sum(trapped_hydrogen_c/DIVIDING_AREA*1e4*TOTAL_TIME/time))
print('Vacancies:', np.sum(vacancy_c) - vacancy_c[0])
print('SIAs:', np.sum(sia_c) - sia_c[0])
# print('Num Free H:', np.sum(free_hydrogen_c)-free_hydrogen_c[0])
for i in range(len(trapped_hydrogen_c)):
	if i != 0:
		trapped_hydrogen_c[i] /= volumeAtIndex(i)
		free_hydrogen_c[i] /= volumeAtIndex(i)

all_hydrogen_c = free_hydrogen_c + trapped_hydrogen_c

for i in range(len(trapped_hydrogen_c)):
	trapped_hydrogen_c[i] = trapped_hydrogen_c[i] / (DENSITY) * 100
	free_hydrogen_c[i] = free_hydrogen_c[i] / (DENSITY) * 100

# Apply Butterworth filter with filtfilt for zero phase shift
def lowpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
	sos = scipy.signal.butter(poles, cutoff, 'lowpass', fs=sample_rate, output='sos')
	filtered_data = scipy.signal.sosfiltfilt(sos, data)
	return filtered_data	

cutoff = 5  # Cutoff frequency
fs = 1 / (positions[5] - positions[4])  # Sampling frequency

# Create a 5-pole low-pass filter with an 80 Hz cutoff
# b, a = scipy.signal.butter(5, 2.5, fs=fs)

# Apply the filter using Gustafsson's method
# smoothed_hydrogen_c = scipy.signal.filtfilt(b, a, trapped_hydrogen_c[2:], method="gust")
smoothed_hydrogen_c = scipy.signal.savgol_filter(trapped_hydrogen_c[2:], 50, 1)

concentrations = [c for c in concentrations]

retained_experiment_fluence = 0  # arbitrary units
retained_sim_fluence = 0
for i in range(len(experiment_positions)-1):
	# if i == 0:
		# print(concentrations[i] * (experiment_positions[i+1]-experiment_positions[i]))
	retained_experiment_fluence += concentrations[i] * (experiment_positions[i+1]-experiment_positions[i])
for i in range(len(positions)-1):
	retained_sim_fluence += trapped_hydrogen_c[i] * (positions[i+1]-positions[i])
# print(retained_experiment_fluence)
print()
# print("Sim retained vs. experiment retained: "+str(retained_sim_fluence/retained_experiment_fluence))
if plot_h:
	# plt.plot(positions[2:], free_hydrogen_c[2:], label="Free Hydrogen Concentration", marker='^', linestyle='-', markersize=0)
	plt.plot(positions[2:], trapped_hydrogen_c[2:], label="Simulation", alpha=0.3, marker='^')
	plt.plot(positions[2:], smoothed_hydrogen_c[:], label="Simulation Filtered", color='blue', marker='^', markersize=0)
	# plt.plot(positions[2:], all_hydrogen_c[2:], label="Hydrogen Concentration")
# if plot_v:
	# indices_to_delete = [i for i in range(len(vacancy_c)) if vacancy_c[i] == 0]		
	# positions_vacancy = np.delete(positions, indices_to_delete)
	# nonzero_vacancy_c = np.delete(vacancy_c, indices_to_delete)
	# plt.plot(positions[:upto], vacancy_c[:upto], color='r', label="Vacancy Concentration")
# print("Summed retained concentration: "+str(np.sum(trapped_hydrogen_c)))
# plt.axhline(y=H_SATURATION_CONCENTRATION, color='black', linestyle='--', label="Free Hydrogen Saturation Limit")
plt.yscale('log')
plt.ylim(2*10**-3, 10**0)
plt.xlim(0, 5)
plt.plot(experiment_positions, concentrations, label="Experiment", color='r')
plt.legend()
plt.title("Trapped Hydrogen Concentration Vs. Depth\n $T = 383K, Fluence = 1 \cdot 10^{24}$ $[m^{-2}]$")
plt.xlabel("Depth $[\mu m]$")
plt.ylabel("Trapped Hydrogen Concentration $[at. \%]$")
plt.show()
