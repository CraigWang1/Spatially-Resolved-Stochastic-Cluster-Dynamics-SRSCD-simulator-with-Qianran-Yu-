import os, re, cv2, time
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from textwrap import wrap

# matplotlib.style.use('ggplot')


# Change data files list, times list, and flux for custom use case
POINTS = 100                             # num spatial elements in the simulation (1 surface + 100 bulk)
VOLUME = 1e-17                           # volume of a spatial element [cm^3]
SURFACE_THICKNESS = 0.544                # [nm]
SURFACE_VOLUME = VOLUME / 20 * SURFACE_THICKNESS + VOLUME # [cm^3]
DENSITY = 6.30705e+22                      # [atoms/cm^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             # [ev/K] Boltzmann's constant.
TEMPERATURE = 1000
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

positions_sim = [10.272 + i*20 for i in range(100)]

plt.figure(figsize=(640/dpi, 480/dpi), dpi=dpi)

fluences = [
	"5e22",
	# "5e23",
	# "1e24"
]

for i in range(len(fluences)):
	with open(f"/home/craig/research/experiment_retention_300K/fluence_{fluences[i]}.txt") as f:
		f.readline() # column titles
		positions = []
		concentrations = []
		for line_hold in f:
			line_hold = line_hold.split(", ")
			depth_micrometer = float(line_hold[0])
			concentration_at = float(line_hold[1]) # at. % units
			positions.append(depth_micrometer)
			concentrations.append(concentration_at)
		# plt.plot(positions, concentrations, label=f"Experiment {fluences[i]}")
		plt.plot(positions, concentrations, label="Experiment")

plt.axhline(y=H_SATURATION_CONCENTRATION, color='black', linestyle='--', label="Free Hydrogen Saturation Limit")

plt.legend()
plt.title("Trapped Hydrogen Concentration Vs. Depth\n $T = 300K, Fluence = 5 \cdot 10^{22}$ $[m^{-2}]$")
plt.xlabel("Depth $[\mu m]$")
plt.ylabel("Concentration $[at. \%]$")
plt.show()