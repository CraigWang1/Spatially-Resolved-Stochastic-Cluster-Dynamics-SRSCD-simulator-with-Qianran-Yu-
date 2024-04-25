import os, re, cv2
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from textwrap import wrap

# matplotlib.style.use('ggplot')


# Change data files list, times list, and flux for custom use case

flux = 1.6199111 * 10 ** 20              # [1/(m^2*s)]
l = 2000.544 * 10 ** -9                  # total length of sample [m]
D = 9.96964 * 10 ** -12                  # diffusivity of H (m^2/s)
POINTS = 101                             # num spatial elements in the simulation (1 surface + 100 bulk)
VOLUME = 1e-23                           # volume of a spatial element [m^3]
SURFACE_THICKNESS = 0.544                # [nm]
SURFACE_VOLUME = VOLUME / 20 * SURFACE_THICKNESS  # [m^3]
DENSITY = 6.30705e+28                      # [atoms/m^3] Atomic density for W.
HEAT_OF_SOLUTION = 1.04                    # [eV] Heat of solution of H in W.
KB = 8.617e-05                             #[ev/K] Boltzmann's constant.
TEMPERATURE = 1700
H_SATURATION_CONCENTRATION = DENSITY * math.exp(-HEAT_OF_SOLUTION/KB/TEMPERATURE);


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

out = cv2.VideoWriter('output.mp4', cv2.VideoWriter_fourcc(*'mp4v'), 20.0, (640,480))

positions = [10.272*10**-9 + i*20*10**-9 for i in range(100)]
positions.insert(0, 0)
dpi = 100
count = 0

for filename in tqdm(sorted(os.listdir(f"{TEMPERATURE}K"), key=lambda x:float(re.findall("(\d+)",x)[0]))):
	count += 1
	if count > 600:
		break
	with open(f"1700K New/{filename}") as f:
		hydrogen_c = np.zeros(POINTS)
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
			# if h_per_cluster > 0:
			if obj_key > 0 and obj_key < 1000000:
				hydrogen_c += np.array(line_hold[2:]).astype(float) * h_per_cluster
				plot_h = True

			v_per_cluster = vacancies_per_cluster(line_hold[1])
			if v_per_cluster > 0:
				vacancy_c += np.array(line_hold[2:]).astype(float) * v_per_cluster
				plot_v = True
		hydrogen_c[0] *= (20/0.544)
		vacancy_c[0] *= (20/0.544)
		hydrogen_c /= VOLUME
		vacancy_c /= VOLUME

		plt.figure(figsize=(640/dpi, 480/dpi), dpi=dpi)
		plt.axhline(y=H_SATURATION_CONCENTRATION, color='black', linestyle='--', label="Hydrogen Saturation Limit")
		if plot_h:
			plt.plot(positions, hydrogen_c, label="Free Hydrogen Concentration")
		# if plot_v:
			# plt.plot(positions, vacancy_c, label="Vacancy Concentration", color='r')
		plt.title(f"Free Hydrogen Concentration vs. Position\n$T = {TEMPERATURE} K, t = {round(time*1e6, 1)} \mu s$")
		plt.xlabel("Position $[m]$")
		plt.ylabel("Concentration $[m^{-3}]$")
		# if not plot_v:
		# plt.ylim(0, 6e25)
		plt.legend()
		# plt.tight_layout()
		plt.savefig("fig.png")
		plt.close()
		img = cv2.imread("fig.png")
		out.write(img)

out.release()