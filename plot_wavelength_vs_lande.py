import matplotlib.pyplot as plt
import numpy as np
from pylab import *

filenames = ["all_G2_data.txt", "all_K0_data.txt", "all_K5_data.txt", "all_M2_data.txt"]

for item in filenames:
	wavelength = []
	lande = []
	with open(item) as file:
		for line in file:
			column = line.split()
			w = float(column[2])
			l = float(column[3])
			wavelength.append(w)
			lande.append(l)
	plt.plot(wavelength, lande, "ro")
	plt.axis([3500, 7000, -2, 4])
	#NEED TO ADD AXIS LABELS
	plt.show()
