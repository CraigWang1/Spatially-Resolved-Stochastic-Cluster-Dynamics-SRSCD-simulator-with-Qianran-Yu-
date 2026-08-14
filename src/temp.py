import matplotlib.pyplot as plt
import numpy as np
import addcopyfighandler

plt.rcParams.update({'font.size': 14})

"""
colors = ['#2b5c8f', '#d95f02']  # Clean, professional distinct palette

plt.rcParams.update({'font.size': 14})

plt.bar(['Experiment', 'Simulation'], [8.79e19, 8.97e19], width=0.5, color=colors, alpha=0.9)
plt.ylabel("Total Amount Retained $[D/m^{2}]$")
plt.title("Calibrated Retention for Experiment #1")
plt.savefig('fig.png', dpi=300, transparent=True)
plt.show()
"""

"""
surf_coverage = np.linspace(0, 1, 100)

custom_desorbE = 0.8 + 1.4 / (1+np.exp((surf_coverage-0.3)/0.2))
hodille_desorbE = 2.0*(0.525 + 0.591*(1.0/(1.0+np.exp( (surf_coverage-0.247)/0.0692 ))))
ajmalghan_desorbE = -0.00213989*np.exp(5.78271*surf_coverage) + 1.4965

plt.figure(figsize=(10, 6)) 
plt.plot(surf_coverage, custom_desorbE, label='Calibrated in This Work', linewidth=5)
plt.plot(surf_coverage, hodille_desorbE, label='Hodille 2020', alpha=0.5)
plt.plot(surf_coverage, ajmalghan_desorbE, label='Ajmalghan 2019', alpha=0.5)
plt.legend()
plt.title("Comparison of Calibrated Desorption Energy vs. Literature")
plt.ylabel("Desorption Energy (eV)")
plt.xlabel("Surface Site Coverage Fraction")
plt.savefig('fig.png', dpi=300, transparent=True)
plt.show()
"""

"""
temperatures = [383, 823]
H_induced_SAVR = [0.010, 0.079]

temps_interp = np.linspace(150, 1000, 300)
interp = 0.477 * np.exp(-0.128 / 8.617e-5 / temps_interp)

plt.figure(figsize=(10, 6))
plt.scatter(temperatures[0], H_induced_SAVR[0], label='Experiment #1 Fit', color='red', s=100)
plt.scatter(temperatures[1], H_induced_SAVR[1], label='Experiment #2 Fit', color='green', s=100)
plt.plot(temps_interp, interp, label='0.477*exp(-(0.128 eV)/kT)')
plt.xlabel("Temperature [K]")
plt.ylabel(r"H Induced Vacancy Generation Rate Per H [$s^{-1}$]")
plt.title("H Induced Vacancy Generation Rate Per H vs. Temperature")
plt.legend()
plt.tight_layout()
plt.savefig('fig.png', dpi=300, transparent=True)
plt.show()
"""


"""
colors = ['#2b5c8f', '#d95f02']  # Clean, professional distinct palette

plt.figure(figsize=(10, 6))
plt.bar(['Experiment', 'Simulation'], [5.36e20, 5.12e20], width=0.5, color=colors, alpha=0.9)
plt.ylabel("Total Amount Retained $[D/m^{2}]$")
plt.title("Calibrated H Retention for Neutron Damaged Experiment")
plt.tight_layout()
plt.savefig('fig.png', dpi=300, transparent=True)
plt.show()
"""

"""
# For 823K unirr case, before extending to 0.5mm
from collections import Counter

def hydrogen_per_cluster(obj_key):
	# obj_key: string of the object's key
	# Assumes obj_key is in form of xxx, xxx is num H

	if len(obj_key) >= 6:
		return abs(int(obj_key[len(obj_key)-6:]))
	else:
		return int(obj_key)

run1 = {-1000000: 6.0, -1000001: 14.0, -1000002: 11.0, -1000003: 8.0, -1000004: 2.0, -2000004: 3.0, -2000005: 4.0, -2000006: 10.0, -2000007: 1.0, -3000005: 1.0, -3000007: 1.0, -3000008: 3.0, -3000009: 2.0, -3000010: 2.0, -4000010: 1.0}
run2 = {-1000000: 14.0, -1000001: 9.0, -1000002: 11.0, -1000003: 7.0, -1000004: 5.0, -1000005: 1.0, -2000004: 2.0, -2000005: 7.0, -2000006: 6.0, -2000007: 3.0, -3000006: 1.0, -3000007: 1.0, -3000008: 4.0, -4000010: 1.0, -4000011: 1.0, -4000012: 1.0, -5000013: 1.0}
run3 = {-1000000: 4.0, -1000001: 9.0, -1000002: 7.0, -1000003: 5.0, -1000004: 1.0, -2000003: 1.0, -2000004: 3.0, -2000005: 6.0, -2000006: 3.0, -2000007: 1.0, -3000007: 1.0, -3000008: 2.0, -3000010: 1.0, -4000010: 1.0}
run4 = {-1000000: 8.0, -1000001: 11.0, -1000002: 16.0, -1000003: 5.0, -1000004: 1.0, -2000003: 1.0, -2000004: 1.0, -2000005: 7.0, -2000006: 6.0, -2000007: 1.0, -2000008: 1.0, -3000006: 1.0, -3000008: 1.0, -3000009: 1.0, -4000009: 1.0, -4000010: 1.0, -4000011: 1.0, -4000012: 1.0}
run5 = {-1000000: 8.0, -1000001: 8.0, -1000002: 6.0, -1000003: 12.0, -1000004: 1.0, -2000004: 3.0, -2000005: 10.0, -2000006: 9.0, -2000007: 1.0, -3000007: 1.0, -3000008: 1.0, -3000009: 1.0, -4000010: 1.0}
run6 = {-1000000: 8.0, -1000001: 5.0, -1000002: 6.0, -1000003: 5.0, -1000004: 6.0, -1000005: 2.0, -2000005: 4.0, -2000006: 4.0, -2000007: 2.0, -3000008: 3.0, -4000008: 2.0, -4000009: 1.0, -4000010: 2.0, -4000011: 1.0, -5000013: 1.0}
run7 = {-1000000: 8.0, -1000001: 12.0, -1000002: 11.0, -1000003: 4.0, -1000004: 1.0, -2000004: 1.0, -2000005: 4.0, -2000006: 4.0, -3000007: 1.0, -4000010: 2.0, -4000011: 1.0}
run8 = {-1000000: 7.0, -1000001: 15.0, -1000002: 12.0, -1000003: 8.0, -1000004: 1.0, -2000004: 2.0, -2000005: 4.0, -2000006: 6.0, -2000007: 2.0, -3000007: 1.0, -3000008: 1.0, -4000009: 2.0, -4000010: 1.0, -4000013: 1.0}
run9 = {-1000000: 18.0, -1000001: 15.0, -1000002: 12.0, -1000003: 6.0, -1000004: 1.0, -2000005: 4.0, -2000006: 5.0, -2000007: 1.0, -3000007: 3.0, -4000008: 1.0, -4000009: 1.0, -4000010: 1.0, -4000012: 2.0}
run10 = {-1000000: 5.0, -1000001: 8.0, -1000002: 13.0, -1000003: 4.0, -1000004: 4.0, -2000004: 1.0, -2000005: 1.0, -2000006: 7.0, -2000007: 2.0, -3000006: 1.0, -3000007: 3.0, -3000008: 3.0, -4000009: 1.0, -4000010: 1.0, -5000012: 1.0, -6000013: 1.0}

total_counts = dict(sorted(dict(Counter(run1) + Counter(run2) + Counter(run3) + Counter(run4) + Counter(run5) + Counter(run6) + Counter(run7) + Counter(run8) + Counter(run9) + Counter(run10)).items(), reverse=True))
total_keys = np.array(list(total_counts.keys()))
total_probs = np.array(list(total_counts.values())) / np.sum(list(total_counts.values()))
total_cdf = np.cumsum(total_probs)
num_H_per_key = [hydrogen_per_cluster(str(key)) for key in total_keys]

print(list(total_keys), '\n')
print(num_H_per_key, '\n')
print(list(total_cdf), '\n')
"""

# For W irr H retention case, before extending to 0.5mm
from collections import Counter

def hydrogen_per_cluster(obj_key):
	# obj_key: string of the object's key
	# Assumes obj_key is in form of xxx, xxx is num H
	
	if len(obj_key) >= 6:
		return abs(int(obj_key[len(obj_key)-6:]))
	else:
		return int(obj_key)

run1 = {-1000000: 10.0, -1000001: 9.0, -1000002: 14.0, -1000003: 5.0, -1000004: 1.0, -1000005: 1.0, -2000004: 1.0, -2000005: 2.0, -2000006: 6.0, -2000007: 2.0, -3000008: 4.0, -3000009: 2.0, -4000009: 1.0, -4000011: 1.0, -5000014: 1.0, -10000026: 1.0, -20000038: 1.0, -20000039: 1.0, -32000050: 1.0, -347000244: 1.0, -2302000690: 1.0, -3136000820: 1.0, -6258001156: 1.0}
run2 = {-1000000: 8.0, -1000001: 9.0, -1000002: 7.0, -1000003: 4.0, -1000004: 4.0, -2000003: 1.0, -2000004: 2.0, -2000005: 5.0, -2000006: 8.0, -2000007: 1.0, -2000008: 1.0, -3000006: 1.0, -3000007: 3.0, -3000008: 5.0, -4000010: 2.0, -5000012: 2.0, -9000021: 1.0, -11000023: 1.0, -15000029: 1.0, -75000097: 1.0}
run3 = {-1000000: 6.0, -1000001: 5.0, -1000002: 12.0, -1000003: 5.0, -1000004: 4.0, -2000005: 8.0, -2000006: 11.0, -2000007: 3.0, -3000007: 1.0, -3000008: 2.0, -3000010: 1.0, -4000009: 1.0, -4000011: 1.0, -5000013: 1.0, -7000017: 1.0, -680000362: 1.0}
run4 = {-1000000: 9.0, -1000001: 12.0, -1000002: 10.0, -1000003: 10.0, -1000004: 6.0, -1000006: 1.0, -2000004: 2.0, -2000005: 4.0, -2000006: 4.0, -3000007: 1.0, -3000008: 4.0, -3000009: 2.0, -4000009: 2.0, -4000010: 1.0, -4000011: 1.0, -4000012: 1.0, -7000015: 1.0, -58000079: 1.0, -218000183: 1.0, -639000348: 1.0}
run5 = {-1000000: 9.0, -1000001: 8.0, -1000002: 8.0, -1000003: 11.0, -1000004: 2.0, -2000003: 1.0, -2000004: 1.0, -2000005: 8.0, -2000006: 8.0, -3000007: 2.0, -4000010: 1.0, -4000011: 2.0, -6000016: 1.0, -7000017: 1.0, -12000027: 1.0, -17000032: 1.0, -225000193: 1.0, -228000198: 1.0}

total_counts = dict(sorted(dict(Counter(run1) + Counter(run2) + Counter(run3) + Counter(run4) + Counter(run5)).items(), reverse=True))
total_keys = np.array(list(total_counts.keys()))
total_probs = np.array(list(total_counts.values())) / np.sum(list(total_counts.values()))
total_cdf = np.cumsum(total_probs)
num_H_per_key = [hydrogen_per_cluster(str(key)) for key in total_keys]

print(list(total_keys), '\n')
print(num_H_per_key, '\n')
print(list(total_cdf), '\n')