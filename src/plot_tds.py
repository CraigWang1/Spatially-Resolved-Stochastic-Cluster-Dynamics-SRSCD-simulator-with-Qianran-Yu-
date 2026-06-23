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

DIVIDING_AREA = 0.64e-16 # m^2

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
	time = (temperature - 300) * 2
	flux = experiment_desorbed_flux[i]
	total_fluence += flux * (time - prev_time)
	prev_time = time

print(total_fluence)

starting_temp = 300
temperatures = [starting_temp]
desorbed_flux = [0]
window_size = 50
for i in range(window_size, len(times)-window_size):
	avg_time = (times[i + window_size] + times[i - window_size])/2
	dt = times[i+window_size] - times[i-window_size]
	dN = desorbed[i+window_size] - desorbed[i-window_size]
	temperatures.append(starting_temp + avg_time * 0.5) # 0.5 K/s heating
	desorbed_flux.append(dN/dt/DIVIDING_AREA)
temperatures.insert(1, temperatures[0]+(temperatures[1]-temperatures[0])*0.7)
desorbed_flux.insert(1, 0)

temperatures.pop(0)
desorbed_flux.pop(0)

plt.plot(temperatures, desorbed_flux, label="Simulation", color='b')
plt.plot(experiment_temperatures, experiment_desorbed_flux, color='r', label="Experiment")
plt.xlabel("Temperature $[K]$")
plt.ylabel("Desorption Flux $[D/m^{2}/s]$")
plt.title("Desorbed Flux Vs. Temperature")
plt.legend()
plt.show()