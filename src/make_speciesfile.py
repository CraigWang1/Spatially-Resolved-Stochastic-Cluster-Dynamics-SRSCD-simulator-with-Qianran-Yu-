"""
Combine multiple speciesx.txt files from different processors into one species.txt file
"""

import numpy as np
from math import *

def combine_species_files():
	NUM_PROCS = 1
	partitions = []
	objCounts = {}

	for i in range(NUM_PROCS):
		with open(f"species{i}.txt") as f:
			if i == 0:
				step = int(f.readline().split()[2])
				time = float(f.readline().split()[2])
				fluenceH = float(f.readline().split()[2])
			else:
				for i in range(3):
					f.readline()
			startIdx = int(f.readline().split()[2])
			endIdx = int(f.readline().split()[2])
			print(startIdx, endIdx)
			for line in f:
				line = line.split()
				objKey = int(line[1])
				line = line[2:]
				counts = np.zeros(len(line), dtype=int)
				for j in range(startIdx, endIdx+1):
					counts[j] = int(line[j])
				if objKey in objCounts:
					objCounts[objKey] += counts
				elif np.sum(counts) > 0:
					objCounts[objKey] = counts
	print()
	f = open("species.txt", "w")
	f.write(f"step = {step}\n")
	f.write(f"time = {time}\n")
	f.write(f"fluenceH = {fluenceH}\n")

	for objKey in objCounts:
		f.write(f"object {objKey}")
		for count in objCounts[objKey]:
			f.write(f"    {count}")
		f.write("\n")

	f.close()

if __name__ == "__main__":
	combine_species_files()
