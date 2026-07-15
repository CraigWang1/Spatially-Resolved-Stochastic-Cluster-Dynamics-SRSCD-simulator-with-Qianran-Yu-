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
from scipy.signal import butter, filtfilt, savgol_coeffs
from scipy.interpolate import make_interp_spline
from make_speciesfile import combine_species_files

DIVIDING_AREA = 0.5141e-16 # m^2
STARTING_TEMP = 300   # K
TEMP_RISE_RATE = 0.5  # K/s
FWHM = 140            # Full width at half maximum in K (for smoothing later)

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

# Resample to get uniform spacing between points
times_uniform = np.linspace(min(times), max(times), 5000)
desorbed_uniform = np.interp(times_uniform, times, desorbed)
dt = times_uniform[1] - times_uniform[0]

times = times_uniform
desorbed = desorbed_uniform

desorbed_flux = np.gradient(desorbed, dt)/DIVIDING_AREA
temperatures = STARTING_TEMP + times*TEMP_RISE_RATE

desorbed_flux = np.insert(desorbed_flux, 0, 0)
temperatures = np.insert(temperatures, 0, temperatures[0]-5)

plt.plot(temperatures, desorbed_flux)
plt.show()

window_size = int( (FWHM*0.15 / (TEMP_RISE_RATE*dt)) )
print('window size:', window_size)

def zero_phase_ma(data, window_size):
    """
    Applies a zero-phase (lag-free) moving average filter to offline data.
    """
    # Create the moving average filter coefficients (b) and denominator (a)
    b = np.ones(window_size) / window_size
    a = 1.0
    
    # Use filtfilt to run the filter forwards and backwards to remove lag
    filtered_data = filtfilt(b, a, data)
    
    return filtered_data

# Apply the forward-backward zero-lag filter
smooth_flux = zero_phase_ma(desorbed_flux, window_size=window_size)

# plt.xlim([350, 1050])

plt.plot(temperatures, smooth_flux, label="Simulation", color='b')
plt.plot(experiment_temperatures, experiment_desorbed_flux, color='r', label="Experiment")
plt.xlabel("Temperature $[K]$")
plt.ylabel("Desorption Flux $[D/m^{2}/s]$")
plt.title("Desorbed Flux Vs. Temperature")
plt.legend()
plt.show()