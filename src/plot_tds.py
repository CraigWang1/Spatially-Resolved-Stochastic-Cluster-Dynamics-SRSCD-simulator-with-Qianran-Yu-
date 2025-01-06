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
from make_speciesfile import combine_species_files

DIVIDING_AREA = 5.0e-16 # m^2

times = []
desorbed = []
with open("Desorbed.txt", "r") as f:
	for line_hold in f:
		line_hold = line_hold.split()
		time = float(line_hold[0])
		numH = int(line_hold[1])
		if len(desorbed) == 0:
			times.append(time)
			desorbed.append(numH)
		elif numH != desorbed[-1] and time != times[-1]:
			times.append(time)
			desorbed.append(numH)

experiment_temperatures = []
experiment_desorbed_flux = []
with open("/home/craig/research/experiment_retention_383K/tds.txt") as f:
	f.readline() # Header
	for line_hold in f:
		line_hold = line_hold.split(", ")
		experiment_temperatures.append(float(line_hold[0]))
		experiment_desorbed_flux.append(float(line_hold[1])*10**17)

total_fluence = 0
prev_time = 0
for i in range(len(experiment_temperatures)-1):
	temperature = experiment_temperatures[i]
	time = (temperature - 350) * 2
	flux = experiment_desorbed_flux[i]
	total_fluence += flux * (time - prev_time)
	prev_time = time

print(total_fluence)


temperatures = [350, 360, 370, 380, 390, 400, 410]
desorbed_flux = [0, 0, 0, 0, 0, 0, 0]
for i in range(len(times)-200):
	temperatures.append(350 + times[i] * 0.5) # 0.5 K/s heating
	dt = times[i+200] - times[i]
	dN = desorbed[i+200] - desorbed[i]
	desorbed_flux.append(dN/dt/DIVIDING_AREA)

temperatures.append(800)
desorbed_flux.append(0)
temperatures.append(1050)
desorbed_flux.append(0)

weight = 0.01
current = desorbed_flux[0]
avg = [current]
for i in range(1, len(temperatures)):
	current = weight * desorbed_flux[i] + (1-weight) * current
	avg.append(current)

plt.plot(temperatures, desorbed_flux, label="Simulation", color='b')
plt.plot(experiment_temperatures, experiment_desorbed_flux, color='r', label="Experiment")
plt.xlabel("Temperature $[K]$")
plt.ylabel("Desorption Flux $[D/m^{2}/s]$")
plt.title("Desorbed Flux Vs. Temperature")
plt.legend()
plt.show()