import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
import pyfits, time, itertools, random, os, sys, glob, re
import matplotlib.mlab as mlab
import pylab as py
from scipy import interpolate

start_time = time.time()

import re
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

ccf_original = sorted(glob.glob("Spectra/ccf_G2_a/*.fits"), key = numericalSort) 
e2ds_original = sorted(glob.glob("Spectra/e2ds_A/*.fits"), key = numericalSort)

ccf = []
e2ds = []
SN_each_file = []

for file in ccf_original:
	SN = []
	for i in range(0,72):
		signal_to_noise = pyfits.open(file)[0].header["HIERARCH ESO DRS SPE EXT SN%s " % (str(i))]
		SN.append(signal_to_noise)

	mean_SN = np.mean(SN)

	if mean_SN > 60:
		SN_each_file.append(mean_SN)
		ccf.append(file)
		print file

for file in e2ds_original:
	SN = []
	for i in range(0,72):
		signal_to_noise = pyfits.open(file)[0].header["HIERARCH ESO DRS SPE EXT SN%s " % (str(i))]
		SN.append(signal_to_noise)

	mean_SN = np.mean(SN)

	if mean_SN > 60:
		e2ds.append(file)
		print file

ccf_write = np.savetxt("Spectra/ccf_file.txt", np.array(ccf))
e2ds_write = np.savetxt("Spectra/e2ds_file.txt", np.array(e2ds))
SN_write = np.savetxt("Spectra/SN_file.txt", np.array(SN_each_file))

end_time = time.time() - start_time
print end_time

