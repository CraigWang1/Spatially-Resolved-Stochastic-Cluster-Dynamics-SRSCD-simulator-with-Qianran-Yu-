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
plt.show()
"""

"""
surf_coverage = np.linspace(0, 1, 100)

custom_desorbE = 0.8 + 1.4 / (1+np.exp((surf_coverage-0.3)/0.2))
hodille_desorbE = 2.0*(0.525 + 0.591*(1.0/(1.0+np.exp( (surf_coverage-0.247)/0.0692 ))))
ajmalghan_desorbE = -0.00213989*np.exp(5.78271*surf_coverage) + 1.4965

# plt.figure(figsize=(10, 6)) 
plt.plot(surf_coverage, custom_desorbE, label='Calibrated in This Work', linewidth=5)
plt.plot(surf_coverage, hodille_desorbE, label='Hodille 2020', alpha=0.5)
plt.plot(surf_coverage, ajmalghan_desorbE, label='Ajmalghan 2019', alpha=0.5)
plt.legend()
plt.title("Comparison of Calibrated Desorption Energy \nvs. Literature")
plt.ylabel("Desorption Energy (eV)")
plt.xlabel("Surface Site Coverage Fraction")
plt.show()
"""

temperatures = [383, 823]
H_induced_SAVR = [0.10, 0.48]

temps_interp = np.linspace(150, 1000, 300)
interp = 1.88 * np.exp(-0.097 / 8.617e-5 / temps_interp)

plt.scatter(temperatures, H_induced_SAVR, color='black')
plt.plot(temps_interp, interp, label='1.88*exp(-0.097/kT)')
plt.xlabel("Temperature [K]")
plt.ylabel(r"H Induced SAV Rate Per H [$s^{-1}$]")
plt.title("H Induced SAV Rate Per H vs. Temperature")
plt.legend()
plt.tight_layout()
plt.show()