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
from functools import partial

#Defining a function to fit the gaussian
def gaussian(x, amp, cen, wid, c):
    return amp*np.exp(-(x-cen)**2/(2*wid**2)) + c

#Creates function to sectionalize code for parallelization
def section(start, end):
	std_hist_array = []
	std_MCavg_array = []
	mad_array = []
	mean_array = []
	amp_array = []

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

	#Making each section of the multiprocessor start where it needs to
	count_del = 0
	while count_del < int(start):
		max_index = np.argmax(lande_g)
		lande_g = np.delete(lande_g, [max_index])
		red = np.delete(red, [max_index])
		blue = np.delete(blue, [max_index])
		weight = np.delete(weight, [max_index])
		count_del += 1

	#Builds loop to remove 1000 lines from stellar spectra
	count_lines = int(start)
	while count_lines < int(end):
		start2 = time.time()

		#Builds loop to remove given # of lines at a time based on maximum lande g values
		count_del = 0
		while count_del < 10:
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
		while count_MC < 5:
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
				#order_noise = noise[j] --> I think this needs to be not included

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

		#Computes mean of standard deviation of all MC gaussian fits (should decrease) in m/s
		mean_wid = np.mean(width) * 1000
		mean_amp = np.mean(amp)
		
		#Adds calculated values to empty lists above
		std_hist_array.append(sigma)
		std_MCavg_array.append(mean_wid)
		mad_array.append(mad)
		mean_array.append(mean)
		amp_array.append(mean_amp)

		#Just to check!
		print "STD from hist: " + str(sigma)
		print "STD from MC average: " + str(mean_wid)
		print "MAD: " + str(mad)
		print "RV value (mean): " + str(mean)
		print "Depth of ccf: " + str(mean_amp)

		#Finishes loop to remove lines and prints count of lines removed
		print "Finished count %d" % count_lines
		count_lines += 10

	#When func is called, returns 4 arrays with info that appear in each section of results
	return std_hist_array, std_MCavg_array, mad_array, mean_array, amp_array

#Splits job up into specific sections, to be fastest must have more sections than # of CPUs
x = []
for xval in range(0,60,20):
	x.append(xval)

y = []
for yval in range(20,80,20):
	y.append(yval)

arg_list = list(zip(x,y))
print "Arg List is:"
print arg_list

#Makes my 2 argument function into a one argument function in order to use .Pool
def multi_run_wrapper(args):
	return section(*args)

#Specifying number of cores based on cores available in computer and calling function by section
pool = mp.Pool(3)
results = pool.map(multi_run_wrapper, arg_list)
print results

final_hist = np.array([])
final_MCavg = np.array([])
final_mad = np.array([])
final_mean = np.array([])
final_amp = np.array([])
for array in np.array(results):
	final_hist = np.append(final_hist, np.array(array[0]))
	final_MCavg = np.append(final_MCavg, np.array(array[1]))
	final_mad = np.append(final_mad, np.array(array[2]))
	final_mean = np.append(final_mean, np.array(array[3]))
	final_amp = np.append(final_amp, np.array(array[4]))

p_std_hist = np.savetxt("/home/katie/files/2_p_std_data_hist.txt", final_hist)
p_mad = np.savetxt("/home/katie/files/2_p_mad_data.txt", final_mad)
p_std_MCavg = np.savetxt("/home/katie/files/2_p_std_data_MCavg.txt", final_MCavg)
p_std_RV = np.savetxt("/home/katie/files/2_p_RV_data.txt", final_mean)
p_amp = np.savetxt("/home/katie/files/2_p_amp_data.txt", final_amp)

###############################################################################################
###     results returns one large array with as many arrays as there are sections in arg_list.
###     The array for each section contains 4 arrays, each with a value for every line deleted.
###     I pull out the information from each individual statistic from each section in order.
###     Then I append this information to an empty array and write the array to a txt file.
###     Data is now ready to analyze.
###############################################################################################


