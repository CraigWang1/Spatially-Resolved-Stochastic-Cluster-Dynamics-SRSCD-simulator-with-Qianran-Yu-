import os, re, cv2, time
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from tqdm import tqdm
from textwrap import wrap
from math import floor, ceil
from scipy.signal import butter, filtfilt

# matplotlib.style.use('ggplot')


# Change data files list, times list, and flux for custom use case
POINTS = 301                             # num spatial elements in the simulation (1 surface + 100 bulk)
VOLUME = 1e-17                           # volume of a spatial element [cm^3]
SURFACE_THICKNESS = 0.544                # [nm]
SURFACE_VOLUME = VOLUME / 20 * SURFACE_THICKNESS + VOLUME # [cm^3]
DENSITY = 6.30705e+22                      # [atoms/cm^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             # [ev/K] Boltzmann's constant.
TEMPERATURE = 300
H_SATURATION_CONCENTRATION = DENSITY * math.exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE) / DENSITY * 100
dpi = 100

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

# out = cv2.VideoWriter('output.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 40.0, (640,480))

# plt.figure(figsize=(640/dpi, 480/dpi), dpi=dpi)

fluences = [
	"5e22",
	# "5e23",
	# "1e24"
]

# Plot Experiment
for i in range(len(fluences)):
	with open(f"/home/craig/research/experiment_retention_300K/fluence_{fluences[i]}.txt") as f:
		f.readline() # column titles
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
with open("restart_store16.txt") as f:
	positions = [0.010272 + i*.020 for i in range(POINTS)]  #micrometer
	trapped_hydrogen_c = np.zeros(POINTS)
	free_hydrogen_c = np.zeros(POINTS)
	vacancy_c = np.zeros(POINTS)
	plot_h = False
	plot_v = False
	f.readline() #step
	time = float(f.readline().split()[2]) #time
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
	trapped_hydrogen_c[0] *= VOLUME / SURFACE_VOLUME
	free_hydrogen_c[0] *= VOLUME / SURFACE_VOLUME
	vacancy_c[0] *= VOLUME / SURFACE_VOLUME

	trapped_hydrogen_c /= VOLUME
	free_hydrogen_c /= VOLUME
	vacancy_c /= VOLUME

	all_hydrogen_c = free_hydrogen_c + trapped_hydrogen_c

	for i in range(len(trapped_hydrogen_c)):
		trapped_hydrogen_c[i] = trapped_hydrogen_c[i] / (DENSITY + trapped_hydrogen_c[i]) * 100

	# Apply Butterworth filter with filtfilt for zero phase shift
	def lowpass(data: np.ndarray, cutoff: float, sample_rate: float, poles: int = 5):
		sos = scipy.signal.butter(poles, cutoff, 'lowpass', fs=sample_rate, output='sos')
		filtered_data = scipy.signal.sosfiltfilt(sos, data)
		return filtered_data	

	cutoff = 5  # Cutoff frequency
	fs = 1 / (positions[1] - positions[0])  # Sampling frequency

	# Create a 5-pole low-pass filter with an 80 Hz cutoff
	b, a = scipy.signal.butter(5, 1.25/3, fs=fs)

	# Apply the filter using Gustafsson's method
	avg_hydrogen_c = scipy.signal.filtfilt(b, a, trapped_hydrogen_c, method="gust")

concentrations = [c*time/26315 for c in concentrations]
plt.plot(experiment_positions, concentrations, label="Experiment", color='r', linewidth=0, marker="^")

upto = POINTS
if plot_h:
	# plt.plot(positions[:upto], free_hydrogen_c[:upto], label="Free Hydrogen Concentration", marker='^', linestyle='-', markersize=0)
	plt.plot(positions[:upto], trapped_hydrogen_c[:upto], label="Simulation", alpha=0.3)
	plt.plot(positions[:upto], avg_hydrogen_c[:upto], label="Simulation Filtered", color='blue', marker='^', linestyle='-', markersize=0)
	# plt.plot(positions[:upto], all_hydrogen_c[:upto], label="Hydrogen Concentration")
# if plot_v:
	# indices_to_delete = [i for i in range(len(vacancy_c)) if vacancy_c[i] == 0]		
	# positions_vacancy = np.delete(positions, indices_to_delete)
	# nonzero_vacancy_c = np.delete(vacancy_c, indices_to_delete)
	# plt.plot(positions[:upto], vacancy_c[:upto], color='r', label="Vacancy Concentration")
print(np.sum(trapped_hydrogen_c))
plt.axhline(y=H_SATURATION_CONCENTRATION, color='black', linestyle='--', label="Free Hydrogen Saturation Limit")
plt.legend()
plt.title("Trapped Hydrogen Concentration Vs. Depth\n $T = 300K, Fluence = 5 \cdot 10^{22}$ $[m^{-2}]$")
plt.xlabel("Depth $[\mu m]$")
plt.ylabel("Trapped Hydrogen Concentration $[at. \%]$")
plt.show()