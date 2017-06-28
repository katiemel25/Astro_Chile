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

import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

ccf_original = sorted(glob.glob("Spectra/ccf_G2_a/*.fits"), key = numericalSort) 
e2ds_original = sorted(glob.glob("Spectra/e2ds_A/*.fits"), key = numericalSort)

ccf_read = []
e2ds_read = []

for file in ccf_original:
	read = pyfits.open(file, memmap = True)
	ccf_read.append(read)

for file in e2ds_original:
	read = pyfits.open(file, memmap = True)
	e2ds_read.append(read)

ccf = []
e2ds = []
SN_each_file = []
bary_earth = []

for f in range(len(ccf_read)):
	SN = []
	earth = ccf_read[f][0].header["HIERARCH ESO DRS BERV "]
	for i in range(0,72):
		signal_to_noise = ccf_read[f][0].header["HIERARCH ESO DRS SPE EXT SN%s " % (str(i))]
		SN.append(signal_to_noise)

	mean_SN = np.mean(SN)

	if mean_SN > 60:
		SN_each_file.append(mean_SN)
		bary_earth.append(earth)
		ccf.append(ccf_read[f])
		e2ds.append(e2ds_read[f])

max_SN_index = np.argmax(SN_each_file)

#####################################################

hdu1 = ccf[max_SN_index]
hdu2 = e2ds[max_SN_index]
RVC = hdu1[0].header["HIERARCH ESO DRS CCF RVC "] * 1000

pixels = hdu2[0].header["NAXIS1 "] #4096
orders = hdu2[0].header["NAXIS2 "] #72

xarr = np.arange(pixels)
wave = np.zeros((pixels,orders),float)

all_coeffs = []
for i in range(0,72):
	coeffs = [] 
	index_0 = hdu2[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4))]
	index_1 = hdu2[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +1))]
	index_2 = hdu2[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +2))]
	index_3 = hdu2[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +3))]
	coeffs.append(index_0)
	coeffs.append(index_1)
	coeffs.append(index_2)
	coeffs.append(index_3)
	all_coeffs.append(coeffs)

all_waves = []
for order in all_coeffs:
	waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
	all_waves.append(waves)

data_t = hdu2[0].data

c = 299792458 #speed of light

shifted = []
for j in range(0,72):
	x = all_waves[j]
	new_x = x / np.sqrt((1 + RVC/c) / (1 - RVC/c))
	shifted.append(new_x)


#plt.plot(shifted[67], data_t[67])
#plt.axvline(x = 6562.801)
#plt.axvline(x = 4862.721)
#plt.show()

#####################################################

del ccf[max_SN_index]
del e2ds[max_SN_index]

all_new_y = []

unshifted_x = []
all_x = []
all_y = []

#####################################################

for l in range(len(ccf)):
	hdu1b = ccf[l]
	hdu2b = e2ds[l]
	RVCb = hdu1b[0].header["HIERARCH ESO DRS CCF RVC "] * 1000

	pixelsb = hdu2b[0].header["NAXIS1 "] #4096
	ordersb = hdu2b[0].header["NAXIS2 "] #72

	xarrb = np.arange(pixelsb)
	waveb = np.zeros((pixelsb,ordersb),float)

	all_coeffsb = []
	for i in range(0,72):
		coeffs = [] 
		index_0 = hdu2b[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4))]
		index_1 = hdu2b[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +1))]
		index_2 = hdu2b[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +2))]
		index_3 = hdu2b[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +3))]
		coeffs.append(index_0)
		coeffs.append(index_1)
		coeffs.append(index_2)
		coeffs.append(index_3)
		all_coeffsb.append(coeffs)

	all_wavesb = []
	for order in all_coeffsb:
		waves = (order[3]*(xarrb**3)) + (order[2]*(xarrb**2)) + (order[1]*(xarrb)) + order[0]
		all_wavesb.append(waves)

	unshifted_x.append(np.array(all_wavesb))

	data_tb = hdu2b[0].data
	all_y.append(np.array(data_tb))

	shiftedb = []
	for j in range(0,72):
		x = all_wavesb[j]
		new_x = x / np.sqrt((1 + RVCb/c) / (1 - RVCb/c))
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

	print e2ds[l]

averaged = np.mean(all_new_y, axis = 0)

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

mean_template_flux = np.mean(np.array(averaged), axis = 1)

####################################################

final_spectra = []
for f in range(len(ccf)):
	mean_flux = np.mean(np.array(all_y[f]), axis = 1)

	scaling = (mean_flux[int(60)]) / (mean_template_flux[int(60)])

	scaled_spectrum = []
	for i in range(0,72):
		divided = np.array(all_y[f])[i] / mean_flux[i]
		multiplied = np.array(divided) * mean_template_flux[i]
		final_spectrum = np.array(multiplied) * scaling
		scaled_spectrum.append(final_spectrum)

	final_spectra.append(scaled_spectrum)

####################################################

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
	
	all_rv = []
	all_cc = []

	x = unshifted_x[f][60]
	y = all_y[f][60]


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

	rv, cc = pyasl.crosscorrRV(x, data, xtem, ftem, -200., 200., 50./200, mode = "doppler")

	cc_min_index = np.argmin(cc)
	center = rv[cc_min_index]

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

		rv, cc = pyasl.crosscorrRV(x, data, xtem, ftem, float(center-20), float(center-20), 50./200, mode = "doppler")
		all_rv.append(np.array(rv) + bary_earth[f])
		all_cc.append(cc)

	all_rv_array = np.array(all_rv)
	mean_rv = np.mean(all_rv_array, axis = 0)
	all_cc_array = np.array(all_cc)
	mean_cc = np.mean(all_cc_array, axis = 0)

	xt = np.array(mean_rv * 1000)
	yt = np.array(mean_cc * 1000)
	popt, pcov = curve_fit(gaussian, xt, yt, p0 = [(min(yt)-max(yt)), xt[np.argmin(yt)], 10000, max(yt)])
	amp.append(popt[0])
	RV.append(popt[1])
	width.append(popt[2])

	plt.plot(xt, gaussian(xt, *popt), label = "Gaussian Fit")
	plt.plot(xt, yt, label = "Data")
	plt.xlim(popt[1]-20000, popt[1]+20000)
	plt.xlabel("RV in m/s")
	plt.ylabel("Flux")
	plt.legend()
	plt.title("RV graph with Gaussian fit")
	plt.savefig("Spectra/spectra_gaussians/fit_index_%d.png" % f)
	plt.close()

print amp
print RV
print width

####################################################

end_time = time.time() - start_time
print end_time

####################################################

