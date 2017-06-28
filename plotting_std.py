import matplotlib.pyplot as plt
import numpy as np
import pyfits
import matplotlib.mlab as mlab

"""gaussian = np.loadtxt("gaussian_data_new.txt").transpose()
amp = abs(gaussian[0])
std = gaussian[2]
RV = gaussian[1]

right2 = []
left2 = []
for i in range(len(RV)):
	left = RV[i] - 3 * std[i]
	right = RV[i] + 3 * std[i]
	left2.append(left)
	right2.append(right)


amp2 = []
for value in amp:
	new = value / max(amp)
	amp2.append(new)

print amp2

sigma1 = abs(gaussian[1])
cut_sigma = sigma1[:741]

#index =  np.where(sigma == -43.18872095079048279)
#print index

mean = np.mean(cut_sigma)
sigma = np.std(cut_sigma)
print sigma

difference = cut_sigma[1] - min(cut_sigma)

if difference > 3 * sigma:
	print 3 * sigma
	print difference

lines_removed = np.arange(len(amp))
plt.plot(lines_removed, amp2)
plt.xlabel("Number of lines removed")
plt.ylabel("Flux of fit")
plt.show()"""

"""
RV_m = np.loadtxt("RV2000test.txt")

x = np.linspace(min(RV_m), max(RV_m))

mean = np.mean(RV_m)
sigma = np.std(RV_m)

#Creates a histogram of RV_m values based on probability and fits gaussian to it
plt.hist(RV_m, normed = 1, bins = 7)
plt.xlim((min(RV_m), max(RV_m)))
plt.title("RV Precision")
plt.xlabel("RV Precision in meters per second")
plt.ylabel("Probability")
plt.plot(x, mlab.normpdf(x, mean, sigma), linewidth = 4, color = 'r')
plt.show()
plt.close()

"""

std = np.loadtxt("blue_std_data_hist.txt")
rand = np.loadtxt("rand_std_data_hist.txt")
longer = np.loadtxt("3long_p_std_data_hist.txt")

std = np.delete(std, [0])
rand = np.delete(rand, [0])

std_short = std[0:40]
rand_short = rand[0:40]

print len(std_short), len(rand_short), len(longer)

lines_removed = np.arange(0, len(std_short)) * 5
#plt.xlim(0,200)
plt.plot(lines_removed, std_short, "blue", label = "Blue Removal")
plt.plot(lines_removed, rand_short, "red", label = "Random Removal")
plt.plot(lines_removed, longer, "green", label = "Removal of Max Lande-g Values")
plt.plot()
plt.xlabel("Number of lines removed")
plt.ylabel("Standard deviation in m/s")
plt.legend(loc = 4)
plt.show()


