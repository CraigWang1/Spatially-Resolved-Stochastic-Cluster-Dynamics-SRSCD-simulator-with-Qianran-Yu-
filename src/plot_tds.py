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

times = [0]
desorbed = [0]
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

temperatures = []
desorbed_flux = []
for i in range(len(times)-200):
	temperatures.append(350 + times[i] * 0.5) # 0.5 K/s heating
	dt = times[i+200] - times[i]
	dN = desorbed[i+200] - desorbed[i]
	desorbed_flux.append(dN/dt/DIVIDING_AREA)

fs = 1 / (times[1] - times[0])  # Sampling frequency

# Create a 5-pole low-pass filter with an 80 Hz cutoff
b, a = scipy.signal.butter(5, 0.01/5, fs=fs)

# Apply the filter using Gustafsson's method
avg = scipy.signal.filtfilt(b, a, desorbed_flux, method="gust")

plt.plot(temperatures, desorbed_flux, label="Simulation")
plt.plot(experiment_temperatures, experiment_desorbed_flux, color='r', label="Experiment")
plt.xlabel("Temperature $[K]$")
plt.ylabel("Desorption Flux $[D/m^{2}/s]$")
plt.title("Desorbed Flux Vs. Temperature")
plt.legend()
plt.show()