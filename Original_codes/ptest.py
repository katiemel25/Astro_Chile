import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
import pyfits, time
import matplotlib.mlab as mlab
import pylab as py
import multiprocessing as mp
from multiprocessing import Pool, freeze_support
import itertools

#Defining a function to fit the gaussian
def gaussian(x, amp, cen, wid, c):
    return amp*np.exp(-(x-cen)**2/(2*wid**2)) + c

#Creates function to sectionalize code for parallelization
def section(start, end):
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

	#Data_t pulls data from fits header info and gives y axis values of stellar spectra
	data_t = hdu[0].data
		
	RV = []
	amp = []
	width = []
		
	#Building loop to do Monte Carlo simulation and find RV precision error
	count_MC = start
	while count_MC < end:
			#Collects values of RV and cc for each of the orders
		all_rv = []		
		all_cc = []

		for j in range(0,72):
				#Gets x (wavelength) value for stellar spectra one order at a time
			x = all_waves[j]

			#data_t is 72 arrays each with 4096 values, giving y for each x
			#data gives one order of stellar spectra data at a time
			data = np.array(data_t[j])

			#Generates array with noise for each order, order_noise is for one order
			max_signal = max(data)
			noise = max_signal / SN
			order_noise = noise[j]

			#Generates Monte Carlo y value within gaussian distribution per x value
			MC_data = data + np.random.normal(0, np.array(order_noise), data.shape)

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

		#Finds mean value of RV and cc from all RV curves of stacked orders
		all_rv_array = np.array(all_rv)
		mean_rv = np.mean(all_rv_array, axis = 0)
		all_cc_array = np.array(all_cc)
		mean_cc = np.mean(all_cc_array, axis = 0)

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

	return RV, amp, width

#Splits job up into specific sections, to be fastest must have more sections than # of CPUs
x = []
for xval in range(0,2000,20):
	x.append(xval)

y = []
for yval in range(20,2020,20):
	y.append(yval)

arg_list = list(zip(x,y))

#Makes my 2 argument function into a one argument function in order to use .Pool
def multi_run_wrapper(args):
	return section(*args)

#Specifying number of cores based on cores available in computer and calling function by section
pool = mp.Pool(50)
results = pool.map(multi_run_wrapper, arg_list)

RV = np.array([])
amp = np.array([])
width = np.array([])
for array in np.array(results):
	RV = np.append(RV, np.array(array[0]))
	amp = np.append(amp, np.array(array[1]))
	width = np.append(width, np.array(array[2]))

	#Turns the 100 RV values from MC from km/s to m/s
RV_m = []
for value in RV:
	meters_per_sec = value * 1000
	RV_m.append(meters_per_sec)

RV_m_array = np.array(RV_m)

#Computes statistical characteristics of RV info after the MC tests
mean = np.mean(RV_m)
sigma = np.std(RV_m)

p_std_hist = np.savetxt("/home/katie/files/RV2.txt", RV_m)

"""y, binedges = np.histogram(RV_m_array) #, normed = True)
y = np.array(y)
binEdges = np.array(binedges)

print binEdges
print y

#x = np.linspace(min(RV_m), max(RV_m))

menStd = np.sqrt(y)
bincenters = 0.5 * (binEdges[1:]+binEdges[:-1])
plt.bar(bincenters, y, color='r', yerr=menStd)
plt.title("RV Precision; Fit Results: mu = %s, std = %s" % (str(mean), str(sigma)))
#plt.plot(x, mlab.normpdf(x, mean, sigma), linewidth = 4, color = 'r')
plt.xlabel("RV Precision in meters per second")
plt.ylabel("Distribution")
plt.show()

#plt.savefig("/home/files/ptest.png")
plt.close()

###################

x = np.linspace(min(RV_m), max(RV_m))

#Creates a histogram of RV_m values based on probability and fits gaussian to it
plt.hist(RV_m, normed = 1)
plt.xlim((min(RV_m), max(RV_m)))
plt.title("RV Precision; Fit Results: mu = %s, std = %s" % (str(mean), str(sigma)))
plt.xlabel("RV Precision in meters per second")
plt.ylabel("Probability")
plt.plot(x, mlab.normpdf(x, mean, sigma), linewidth = 4, color = 'r')
plt.show()
plt.savefig("/home/katie/files/ptest.png")
plt.close()"""

