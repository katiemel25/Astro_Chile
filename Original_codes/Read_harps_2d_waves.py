import pyfits
import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
from scipy.optimize import curve_fit

#Opening file
hdu = pyfits.open("HARPS.2005-09-03T21_28_39.386_e2ds_A.fits")
h = hdu[0].header
ss = hdu[0].shape
m = ss[0] #72
n = ss[1] #4096

#print ss
#print n
#Lesson learned: need to refer back to the original hdu[0] to get functions to work

#Getting number of orders (there are 72)
dim = hdu[0].header['NAXIS2']
#print dim

#size = hdu[0].size
#print size

#Array of pixels
xarr = np.arange(n)
wave = np.zeros((n,m),float)
#print xarr
#print wave
#print wave.shape

string_1 = str("HIERARCH ESO DRS CAL TH COEFF LL")

all_coeffs = []

for i in range(0,72):
	coeffs = [] #print coeffs
	index_0 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4))]
	index_1 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +1))]
	index_2 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +2))]
	index_3 = hdu[0].header["HIERARCH ESO DRS CAL TH COEFF LL%s " % (str(i*4 +3))]
	coeffs.append(index_0)
	coeffs.append(index_1)
	coeffs.append(index_2)
	coeffs.append(index_3)
	all_coeffs.append(coeffs)

all_waves= []

for order in all_coeffs:
	waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
	all_waves.append(waves)
	plt.plot(waves)
	#plt.show()
	#plt.plot(all_waves)
	#plt.show()

#print np.array(all_waves)

#plt.xlabel("Pixel")
#plt.ylabel("Wavelength")
#plt.show()

all_rv = []
all_cc = []

for j in range(0,72):
	x = all_waves[j]
	mask = "G2_masks.txt"
	data = hdu[0].data
	data = data[j]

	x1 = np.arange(x[0]-100., x[0], x[1]-x[0])
	x2 = np.arange(x[-1], x[-1]+100, x[-1]-x[-2])
	xtem = np.hstack([x1, x, x2])
    
	c = 299792.458
    
	lines1, lines2, flux_l = np.genfromtxt(mask, dtype=None, unpack=True)
	nlines = len(flux_l)
    
	ftem = np.zeros(xtem.size)

	#Assigns correct fluxes to template
	for i in range(len(flux_l)):
		indices = np.where((xtem >= lines1[i]) & (xtem <= lines2[i]))[0]
   		if indices.size>0:
    			ftem[indices] = flux_l[i]

	#print len(x)
	#print len(data)
	#print len(xtem)
	#print len(ftem)

	#fig = plt.figure()
	#ax1 = fig.add_subplot(211)
	#ax2 = fig.add_subplot(212)
	#ax1.plot(x, data)
	#ax2.plot(xtem, ftem)
	#plt.show()

	rv, cc = pyasl.crosscorrRV(x, data, xtem, ftem, -20., 20., 50./200, mode = "doppler")
	all_rv.append(rv) #appends all rv values into an array for each order (320 data points)
	all_cc.append(cc) #appends all cc values into an array for each oder (320 data points)

#print len(rv)
#print len(cc)
#print len(all_rv)
#print len(all_cc)

all_rv_array = np.array(all_rv)
mean_index_rv = np.mean(all_rv_array, axis = 0)
all_cc_array = np.array(all_cc)
mean_index_cc = np.mean(all_cc_array, axis = 0)

#print mean_index_rv
#print mean_index_cc

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(mean_index_rv, mean_index_cc)

#Find the index of maximum cross-correlation function
minind = np.argmin(mean_index_cc)

plt.plot(mean_index_rv[minind], mean_index_cc[minind], 'ro')
plt.xlabel("RV in km/s")

plt.show()

#print mean_index_rv[minind]

x = np.array(mean_index_rv)
y = np.array(mean_index_cc)

#n = len(x)
#mean = sum(x*y)/n
#sigma = np.sqrt(sum(y*(x-mean)**2)/n)

def gaussian(x, amp, cen, wid, c):
    return 1 - (amp/(np.sqrt(2*np.pi)*wid))*np.exp(-(x-cen)**2/(2*wid**2)) + c

popt, pcov = curve_fit(gaussian, x, y, p0 = [700000, 0, 10, 1700000])
print popt[1], popt[2]
plt.plot(x, y, 'bo')
plt.plot(x, gaussian(x, *popt))
plt.show()
