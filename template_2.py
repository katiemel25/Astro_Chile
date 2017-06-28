import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
import pyfits, time, itertools, random, os, sys, glob, re
import matplotlib.mlab as mlab
import pylab as py
from scipy import interpolate
import gc

start_time = time.time()

####################################################

#Sorts all file names in order by number (date)
import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

#sorts files into lists
ccf_original = sorted(glob.glob("Spectra/ccf_G2_a/*.fits"), key = numericalSort)
e2ds_original = sorted(glob.glob("Spectra/e2ds_A/*.fits"), key = numericalSort)

#cuts list to test short version if needed
ccf_short = ccf_original #[0:5]
e2ds_short = e2ds_original #[0:5]

#empty arrays for added info later
ccf = []
e2ds = []
all_SN = []
SN_each_file = []
bary_earth_each_file = []
RVC_each_file = []
max_berv_each_file = []
max_berv_wave_each_file = []
th_drift = []
pixels = 4096
orders = 72
c = 299792458

####################################################

#Reads files and extracts important info
for f in range(len(ccf_short)):
	read = pyfits.open(ccf_short[f])
	RVC = read[0].header["HIERARCH ESO DRS CCF RVC "] * 1000
	max_berv = read[0].header["HIERARCH ESO DRS BERVMX "] * 1000
	SN = []
	for i in range(0,72):
		signal_to_noise = read[0].header["HIERARCH ESO DRS SPE EXT SN%s " % (str(i))]
		SN.append(signal_to_noise)

	mean_SN = np.mean(SN)
	all_SN.append(mean_SN)

	#eliminates files with SN less than limit
	if mean_SN > 60:
		RVC_each_file.append(RVC)
		max_berv_each_file.append(max_berv)
		SN_each_file.append(mean_SN)
		ccf.append(ccf_short[f])

####################################################

#calculates shift corrections on file with greatest SN
max_SN_index = np.argmax(all_SN)

max_hdu = pyfits.open(ccf_short[max_SN_index])

xarr = np.arange(pixels)
wave = np.zeros((pixels,orders),float)

all_coeffs = []
for i in range(0,72):
	coeffs = [] 
	index_0 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4))]
	index_1 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +1))]
	index_2 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +2))]
	index_3 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +3))]
	coeffs.append(index_0)
	coeffs.append(index_1)
	coeffs.append(index_2)
	coeffs.append(index_3)
	all_coeffs.append(coeffs)

all_waves = []
for order in all_coeffs:
	waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
	all_waves.append(waves)

data_t = read[0].data

shifted = []
for j in range(0,72):
	x = all_waves[j]
	new_x = x / np.sqrt((1 + RVC_each_file[max_SN_index]/c) / (1 - RVC_each_file[max_SN_index]/c))
	shifted.append(new_x)

#important outputs: data_t (y) and shifted (x)

####################################################

#del ccf_short[max_SN_index]
#del e2ds_short[max_SN_index]

#corrected
all_new_y = []
all_x = []

#originals
unshifted_x = []
all_y = []

####################################################

#calculates shifted values for all files
for g in range(len(e2ds_short)):
	read = pyfits.open(e2ds_short[g])
	earth = read[0].header["HIERARCH ESO DRS BERV "]
	drift = read[0].header["HIERARCH ESO DRS DRIFT SPE RV "]

	SN = []
	for i in range(0,72):
		signal_to_noise = read[0].header["HIERARCH ESO DRS SPE EXT SN%s " % (str(i))]
		SN.append(signal_to_noise)

	mean_SN = np.mean(SN)

	xarr = np.arange(pixels)
	wave = np.zeros((pixels,orders),float)

	if mean_SN > 60:
		bary_earth_each_file.append(earth)
		e2ds.append(e2ds_short[g])
		th_drift.append(drift)

for g in range(len(e2ds)):
	read = pyfits.open(e2ds[g])
	all_coeffs = []
	for i in range(0,72):
		coeffs = [] 
		index_0 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4))]
		index_1 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +1))]
		index_2 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +2))]
		index_3 = read[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +3))]
		coeffs.append(index_0)
		coeffs.append(index_1)
		coeffs.append(index_2)
		coeffs.append(index_3)
		all_coeffs.append(coeffs)

	all_waves = []
	for order in all_coeffs:
		waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
		all_waves.append(waves)

	unshifted_x.append(np.array(all_waves))

	data_tb = read[0].data
	all_y.append(np.array(data_tb))

	shiftedb = []
	for j in range(0,72):
		x = all_waves[j]
		new_x = x / np.sqrt((1 + RVC_each_file[g]/c) / (1 - RVC_each_file[g]/c))
		shiftedb.append(new_x)

	all_x.append(np.array(shiftedb))

	new_y_spectrum = []
	for k in range(0,72):
		test_x = shifted[k]
		test_y = data_t[k]
		test_xb = shiftedb[k]
		test_yb = data_tb[k]
		interpolated = np.interp(test_x, test_xb, test_yb)
		new_y_spectrum.append(interpolated)

	all_new_y.append(new_y_spectrum)

#taking mean of fluxes to get combined spectrum
averaged = np.mean(all_new_y, axis = 0)
mean_template_flux = np.mean(np.array(averaged), axis = 1)

####################################################
###		x values for template in shifted
###		y values for template in averaged
####################################################

"""
for f in range(0,72):
	plt.plot(shifted[f], averaged[f])
	#plt.axvline(x = 6562.801)
	plt.xlabel("Wavelength for order %d" % f) 
	plt.ylabel("Flux")
	#plt.axvline(x = 4862.721)
	plt.savefig("Spectra/order_graphs/high_SN_order_%d.png" % f)
	plt.close()
"""

final_spectra = []
for h in range(len(ccf)):
	mean_flux = np.mean(np.array(all_y[h]), axis = 1)

	scaling = (mean_flux[int(60)]) / (mean_template_flux[int(60)])

	scaled_spectrum = []
	for i in range(0,72):
		divided = np.array(all_y[h])[i] / mean_flux[i]
		multiplied = np.array(divided) * mean_template_flux[i]
		final_spectrum = np.array(multiplied) * scaling
		scaled_spectrum.append(final_spectrum)

	final_spectra.append(scaled_spectrum)

####################################################

#RV calculations on each file individually
def gaussian(x, amp, cen, wid, c):
    return amp*np.exp(-(x-cen)**2/(2*wid**2)) + c

mask = np.loadtxt("G2_masks.txt").transpose()
blue = mask[0]
red = mask[1]
weight = mask[2]

amp = []
RV = []
width = []

for f in range(len(ccf)):
	x = unshifted_x[f][60]
	y = all_y[f][60]
	"""
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

	rv, cc = pyasl.crosscorrRV(x, y, xtem, ftem, -200., 200., 10./200, mode = "doppler", skipedge = 500)

	cc_max_index = np.argmax(cc)
	center = (rv[cc_max_index] + bary_earth_each_file[i] * 1000)
	print center"""

	all_rv = []
	all_cc = []

	for j in range(0,72):
		x = unshifted_x[f][j]
		data = all_y[f][j]

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

		rv, cc = pyasl.crosscorrRV(x, data, xtem, ftem, 0, 75, 50./200, mode = "doppler", skipedge = 500)
		#Bounds (0 and 75) must be wider to account for Max BERV)
		all_rv.append(rv) #+ bary_earth_each_file[f])
		all_cc.append(cc)
		#print "working"

	all_rv_array = np.array(all_rv)
	mean_rv = np.mean(all_rv_array, axis = 0)
	all_cc_array = np.array(all_cc)
	mean_cc = np.mean(all_cc_array, axis = 0)

	adjusted_mean_rv = []
	for value in mean_rv:
		new1 = value + bary_earth_each_file[f]
		new2 = new1 * 1000
		adjusted_mean_rv.append(new2)

	xt = np.array(adjusted_mean_rv)
	yt = np.array(mean_cc)
	popt, pcov = curve_fit(gaussian, xt, yt, p0 = [(min(yt)-max(yt)), xt[np.argmin(yt)], 10000, max(yt)])
	amp.append(popt[0])
	RV.append(popt[1])
	width.append(popt[2])

	plt.plot(xt, gaussian(xt, *popt), label = "Gaussian Fit")
	plt.plot(xt, yt, label = "Data")
	plt.xlim(popt[1]-20000, popt[1]+20000)
	plt.xlabel("RV in m/s")
	plt.ylabel("Flux")
	#plt.legend()
	plt.title("RV graph with Gaussian fit")
	plt.savefig("Spectra/spectra_gaussians/fit_index_%d.png" % f)
	plt.close()
	print "saved fig %d" % f

print RV

drift_adjusted_RV = []
for i in range(len(RV)):
	new = RV[i] - th_drift[i]
	drift_adjusted_RV.append(new)

print drift_adjusted_RV

#max_berv_wave = drift_adjusted_RV * np.sqrt((1 + max_berv_each_file[g]/c) / (1 - max_berv_each_file[g]/c))	
#max_berv_wave_each_file.append(max_berv_wave)

print RVC_each_file

####################################################

end_time = time.time() - start_time
print end_time

####################################################

