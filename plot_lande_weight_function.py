import matplotlib.pyplot as plt
import numpy as np
import pyfits
import matplotlib.mlab as mlab
import powerlaw

reg_lande_g = []
with open("all_G2_data.txt") as template:
	for line in template:
		column = line.split()
		l = (float(column[3]))
		reg_lande_g.append(l)

lande_g1 = abs(np.array(reg_lande_g))

std = np.loadtxt("3long_p_std_data_hist.txt")

sorted_g = sorted(lande_g1, reverse = True)

cut_g = sorted_g[0:200:5]

pre = np.poly1d(np.polyfit(cut_g, std, 3))
#-0.01772 x + 0.1758 x - 0.5785 x + 1.417
p = np.array(pre)

a = p[0] / p[3]
b = p[1] / p[3]
c = p[2] / p[3]
d = 1

def weight(x):
	new_value = a*(x**3) + b*(x**2) + c*(x) + d
	return new_value

weighted_lande = []
for value in cut_g:
	new = weight(value)
	weighted_lande.append(new)

print weighted_lande

"""
plt.plot(cut_g, std)
plt.plot(cut_g, p(cut_g))
plt.gca().invert_xaxis()
plt.plot()
plt.xlabel("Lande g value")
plt.ylabel("Standard deviation in m/s")
plt.show()
"""
