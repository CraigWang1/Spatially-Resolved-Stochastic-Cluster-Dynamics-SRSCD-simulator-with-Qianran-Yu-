import os, re, cv2, time
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from textwrap import wrap
from make_speciesfile import combine_species_files


# Change data files list, times list, and flux for custom use case
POINTS = 302                             # num spatial elements in the simulation (1 surface + 100 bulk)
VOLUME = 1e-17                           # volume of a spatial element [cm^3]
SURFACE_THICKNESS = 0.544                # [nm]
HUGE_THICKNESS = 6000                  # [nm]
SURFACE_VOLUME = VOLUME / 20 * SURFACE_THICKNESS  # [cm^3]
HUGE_VOLUME = VOLUME / 20 * HUGE_THICKNESS        # [cm^3]
DENSITY = 6.30705e+22                      # [atoms/cm^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             # [ev/K] Boltzmann's constant.
TEMPERATURE = 300
H_SATURATION_CONCENTRATION = DENSITY * math.exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE);

combine_species_files(POINTS)

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

positions = [10.272 + i*20 for i in range(POINTS)]
# positions[-1] += HUGE_THICKNESS
dpi = 100
count = 0

# for filename in tqdm(sorted(os.listdir("1700K AllH"), key=lambda x:float(re.findall("(\d+)",x)[0]))):
count += 1
for i in range(1):
	with open(f"species.txt") as f:
		trapped_hydrogen_c = np.zeros(POINTS)
		free_hydrogen_c = np.zeros(POINTS)
		vacancy_c = np.zeros(POINTS)
		plot_h = False
		plot_v = False
		f.readline() #step
		t = float(f.readline().split()[2]) #time
		f.readline() #fluenceH
		for line_hold in f:
			line_hold = line_hold.split()
			obj_key = int(line_hold[1])
			h_per_cluster = hydrogen_per_cluster(line_hold[1])
			# if h_per_cluster > 0:    # for plotting all H
			if obj_key > 0 and obj_key < 1000000:  # free H
				free_hydrogen_c += np.array(line_hold[2:]).astype(float) * h_per_cluster
				plot_h = True
			elif obj_key < 0 or (obj_key > 1000000) and h_per_cluster > 0: # trapped H
				trapped_hydrogen_c += np.array(line_hold[2:]).astype(float) * h_per_cluster
				plot_h = True

			v_per_cluster = vacancies_per_cluster(line_hold[1])
			if v_per_cluster > 0:
				vacancy_c += np.array(line_hold[2:]).astype(float) * v_per_cluster
				plot_v = True
		trapped_hydrogen_c[0] *= VOLUME / SURFACE_VOLUME
		free_hydrogen_c[0] *= VOLUME / SURFACE_VOLUME
		vacancy_c[0] *= VOLUME / SURFACE_VOLUME

		trapped_hydrogen_c[-1] *= VOLUME / HUGE_VOLUME
		free_hydrogen_c[-1] *= VOLUME / HUGE_VOLUME
		vacancy_c[-1] *= VOLUME / HUGE_VOLUME

		trapped_hydrogen_c /= VOLUME
		free_hydrogen_c /= VOLUME
		vacancy_c /= VOLUME

		all_hydrogen_c = free_hydrogen_c + trapped_hydrogen_c

		plt.figure(figsize=(640/dpi, 480/dpi), dpi=dpi)
		plt.axhline(y=H_SATURATION_CONCENTRATION, color='black', linestyle='--', label="Free Hydrogen Saturation Limit")
		upto = POINTS
		if plot_h:
			plt.plot(positions[1:upto], free_hydrogen_c[1:upto], label="Free Hydrogen Concentration", marker='^', linestyle='-', markersize=0)
			# plt.plot(positions[:upto], trapped_hydrogen_c[:upto], label="Trapped Hydrogen Concentration", color='darkgreen', marker='^', linestyle='-', markersize=0)
			# plt.plot(positions[:upto], all_hydrogen_c[:upto], label="Hydrogen Concentration")
		if plot_v:
			indices_to_delete = [i for i in range(len(vacancy_c)) if vacancy_c[i] == 0]		
			positions_vacancy = np.delete(positions, indices_to_delete)
			nonzero_vacancy_c = np.delete(vacancy_c, indices_to_delete)
			plt.plot(positions_vacancy[1:upto], nonzero_vacancy_c[1:upto], label="Nonzero Vacancy Concentration", color='r', linestyle='-', linewidth=0, marker='x')
			plt.plot(positions[1:upto], vacancy_c[1:upto], color='r')
		plt.title(f"Hydrogen Concentration vs. Depth\n$T = {TEMPERATURE} K, t = {round(t*1e6, 1)} \mu s$" + "$, Flux=4.0 \cdot 10^{23}$ $[cm^{-2}s^{-1}]$")
		plt.xlabel("Depth $[nm]$")
		plt.ylabel("Concentration $[cm^{-3}]$")

		# if not plot_v:
		# plt.ylim(0, 6e25)
		plt.legend()
		# plt.draw()
		# plt.tight_layout()
		# plt.savefig("fig.png")
		# plt.pause(5)
		# plt.close()
		plt.show()
		# img = cv2.imread("fig.png")
		# out.write(img)

# out.release()