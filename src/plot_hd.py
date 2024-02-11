import matplotlib
import datetime
import matplotlib.pyplot as plt

def to_int(str):
	split = str.split('e+')
	if len(split) < 2:
		return 0
	num = float(split[0]) * 10 ** int(split[1])
	return num

position = []
concentration = []
with open("hd.txt") as f:
	for line in f.readlines():
		split = line.split()
		position.append(int(split[0]))
		concentration.append(to_int(split[1]))
	plt.title("H Concentration Vs. Position")
	plt.xlabel(r"Position $(nm)$")
	plt.ylabel(r"Concentration $(cm^{-3})$")
	if concentration[-1] > 1.991e20:
		plt.ylim(1.99333e20-0.003e20, 1.99333e20+0.003e20)
	plt.plot(position, concentration, label="Hydrogen Concentration")
	plt.plot([10, 1990], [1.99333e20, 1.99333e20], color='r', linestyle='-', label="Saturation Limit")
	plt.legend()
	save_location = f"images/{datetime.datetime.now().isoformat()}.jpg"
	plt.savefig(save_location)
	# plt.show()