import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
import pyfits, time
import matplotlib.mlab as mlab
import pylab as py

start = time.time()

#Reading in the binary mask with red/blue cutoffs and weight
mask = np.loadtxt("G2_masks.txt").transpose()
blue = mask[0]
red = mask[1]
weight = mask[2]

#Splitting file with line centroid, weight, VALD wavelength, and lande g into arrays
centroid = []
vald_wave =[]
reg_lande_g = []
with open("all_G2_data.txt") as template:
	for line in template:
		column = line.split()
		c = (float(column[0]))
		v = (float(column[2]))
		l = (float(column[3]))
		centroid.append(c)
		vald_wave.append(v)
		reg_lande_g.append(l)

#Making all values of lande g positive (negatives due to sign convention)
lande_g = abs(np.array(reg_lande_g))

#Making 2 histograms of lande g distribution, one with signs and one with abs. value
plt.hist(reg_lande_g)
plt.title("Lande g value occurences")
plt.xlabel("Lande g value")
plt.ylabel("Frequency")
plt.savefig("Lande_MC/Lande_histogram.png")
plt.close()
plt.hist(lande_g)
plt.title("Lande g value occurences")
plt.xlabel("Lande g value")
plt.ylabel("Frequency")
plt.savefig("Lande_MC/Lande_abs_histogram.png")
plt.close()

#Defining a function to fit the gaussian
def gaussian(x, amp, cen, wid, c):
    return amp*np.exp(-(x-cen)**2/(2*wid**2)) + c

#Reading in fits file for HARPS stellar spectra
hdu = pyfits.open("HARPS.2005-09-03T21_28_39.386_e2ds_A.fits")
h = hdu[0].header	
ss = hdu[0].shape
m = ss[0] #number of orders, 72
n = ss[1] #pixels per order, 4096
dim = hdu[0].header['NAXIS2']

#Creates array of zeros with dimensions m and n
xarr = np.arange(n)
wave = np.zeros((n,m),float)
string_1 = str("HIERARCH ESO DRS CAL TH COEFF LL")

#Appending in order coefficients of polynomial from fits header info
#D = index 0, C = index 1, B = index 2, C = index 3, repeat for all 72 orders)
all_coeffs = []
for i in range(0,72):
	coeffs = [] 
	index_0 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4))]
	index_1 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +1))]
	index_2 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +2))]
	index_3 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +3))]
	coeffs.append(index_0)
	coeffs.append(index_1)
	coeffs.append(index_2)
	coeffs.append(index_3)
	all_coeffs.append(coeffs)

#Pulls out signal to noise value for each order from fits header
SN = []
for i in range(0,72):
	signal_to_noise = hdu[0].header["HIERARCH ESO DRS SPE EXT SN%s " % (str(i))]
	SN.append(signal_to_noise)

#From polynomial, building solar spectra order by order
all_waves= []
for order in all_coeffs:
	waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
	all_waves.append(waves)

#Creating empty text files to which information will be added later
blank = []
std_hist = np.savetxt("Lande_MC/std_data_hist.txt", blank, delimiter = ",", newline = "\n")
mad = np.savetxt("Lande_MC/mad_data.txt", blank, delimiter = ",", newline = "\n")
std_MCavg = np.savetxt("Lande_MC/std_data_MCavg.txt", blank, delimiter = ",", newline = "\n")

#Data_t pulls data from fits header info and gives y axis values of stellar spectra
data_t = hdu[0].data

#Builds loop to remove 1000 lines from stellar spectra
count_lines = 0
while count_lines < 1001:
	start2 = time.time()

	#Builds loop to remove 5 lines at a time based on maximum lande g values
	count_del = 0
	while count_del < 5:
		max_index = np.argmax(lande_g)
		lande_g = np.delete(lande_g, [max_index])
		red = np.delete(red, [max_index])
		blue = np.delete(blue, [max_index])
		weight = np.delete(weight, [max_index])
		count_del += 1

	amp = []
	RV = []
	width = []

	#Building loop to do Monte Carlo simulation and find RV precision error
	count_MC = 0
	while count_MC < 2:

		#Collects values of RV and cc for each of the orders
		all_rv = []
		all_cc = []

		for j in range(0,1):
			#Gets x (wavelength) value for stellar spectra one order at a time
			x = all_waves[j]

			#data_t is 72 arrays each with 4096 values, giving y for each x
			#data gives one order of stellar spectra data at a time
			data = np.array(data_t[j])

			max_signal = max(data)
			noise = max_signal / SN[j]
			#order_noise = noise[j]
	
			#Generates Monte Carlo y value within gaussian distribution per x value
			np.random.seed()
			MC_data = data + np.random.normal(0, noise, data.shape)

			#Building template based on range of x wavelengths for one order
			x1 = np.arange(x[0]-100., x[0], x[1]-x[0])
			x2 = np.arange(x[-1], x[-1]+100, x[-1]-x[-2])
			xtem = np.hstack([x1, x, x2])
    		
    		#Assigns correct y values (flux) to template for each x (wavelength)
			ftem = np.zeros(xtem.size)
			for i in range(len(weight)):
				indices = np.where((xtem >= blue[i]) & (xtem <= red[i]))[0]
   				if indices.size>0:
    					ftem[indices] = weight[i]

    		#Creates cross correlation RV values using stellar spectrum and template
			rv, cc = pyasl.crosscorrRV(x, MC_data, xtem, ftem, -20., 20., 50./200, mode = "doppler")
			all_rv.append(rv) 
			all_cc.append(cc)

        	"""
        	fig = plt.figure()
        	ax1 = fig.add_subplot(111)
        	ax2 = fig.add_subplot(211)
        	ax3 = fig.add_subplot(311)
        	ax1.plot(x, MC_data)
        	plt.xlabel("Wavelength (angstroms)")
        	plt.ylabel("Flux")
       		plt.show()"""

	       	"""fig = plt.figure()
			ax1 = fig.add_subplot(111)
			ax1.plot(mean_index_rv, mean_index_cc)
			minind = np.argmin(mean_index_cc)
			plt.plot(mean_index_rv[minind], mean_index_cc[minind], 'ro')
			plt.xlabel("RV in km/s")
			plt.show()"""

		#Finds mean value of RV and cc from all RV curves of stacked orders
		all_rv_array = np.array(all_rv)
		mean_rv = np.mean(all_rv_array, axis = 0)
		all_cc_array = np.array(all_cc)
		mean_cc = np.mean(all_cc_array, axis = 0)

		"""fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.plot(mean_index_rv, mean_index_cc)
		minind = np.argmin(mean_index_cc)
		plt.plot(mean_index_rv[minind], mean_index_cc[minind], 'ro')
		plt.xlabel("RV in km/s")
		plt.show()"""

		#Fits a gaussian to RV cross correlation and appends info to empty arrays
		xt = np.array(mean_rv)
		yt = np.array(mean_cc)
		popt, pcov = curve_fit(gaussian, xt, yt, p0 = [-700000, 0, 10, 1700000])
		amp.append(popt[0])
		RV.append(popt[1])
		width.append(popt[2])

		#Finishes one round of Monte Carlo test and increases count
		print count_MC
		count_MC += 1

	#Turns the 100 RV values from MC from km/s to m/s
	RV_m = []
	for value in RV:
		meters_per_sec = value * 1000
		RV_m.append(meters_per_sec)

	#Computes statistical characteristics of RV info after the MC tests
	mean = np.mean(RV_m)
	variance = np.var(RV_m)
	sigma = np.std(RV_m)
	med = np.median(RV_m)
	mad = np.median(np.abs(np.array(RV_m - med)))

	#Computes mean of standard deviation of all MC gaussian fits (should decrease)
	mean_wid = np.mean(width)

	#Produces graph with x only in the range of the corresponding MC test RV_m values
	x = np.linspace(min(RV_m), max(RV_m))

	#Creates a histogram of RV_m values based on probability and fits gaussian to it
	plt.hist(RV_m, normed = 1)
	plt.xlim((min(RV_m), max(RV_m)))
	plt.title("RV Precision; Fit Results: mu = %s, std = %s" % (str(mean), str(sigma)))
	plt.xlabel("RV Precision in meters per second")
	plt.ylabel("Probability")
	plt.plot(x, mlab.normpdf(x, mean, sigma), linewidth = 4, color = 'r')
	plt.savefig("Lande_MC/RV_precision_2_%d.png" % count_lines)
	plt.close()

	#Appends STD(hist and avg) and MAD info to the three files created at the beginning
	with open("Lande_MC/std_data_hist.txt", "a") as file:
		file.write(str(sigma))
		file.write(str("\n"))
		file.close()
	with open("Lande_MC/mad_data.txt", "a") as file:
		file.write(str(mad))
		file.write(str("\n"))
		file.close()
	with open("Lande_MC/std_data_MCavg.txt", "a") as file:
		file.write(str(mean_wid))
		file.write(str("\n"))
		file.close()
	
	#Finishes loop to remove lines and prints count of lines removed
	print "Finished count %d" % count_lines
	count_lines += 5
