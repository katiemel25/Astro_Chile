import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from PyAstronomy import pyasl
import pyfits, time

start = time.time()

mask = np.loadtxt("G2_masks.txt").transpose()
red = mask[0]
blue = mask[1]
weight = mask[2]

centroid = []
vald_wave =[]
lande_g = []

with open("all_G2_data.txt") as template:
	for line in template:
		column = line.split()
		c = (float(column[0]))
		v = (float(column[2]))
		l = abs(float(column[3]))
		centroid.append(c)
		vald_wave.append(v)
		lande_g.append(l)

def gaussian(x, amp, cen, wid, c):
    return amp*np.exp(-(x-cen)**2/(2*wid**2)) + c

print len(red)
print len(blue)
print len(weight)
print len(lande_g)

hdu = pyfits.open("HARPS.2005-09-03T21_28_39.386_e2ds_A.fits")
h = hdu[0].header	
ss = hdu[0].shape
m = ss[0] 
n = ss[1] 

dim = hdu[0].header['NAXIS2']

xarr = np.arange(n)
wave = np.zeros((n,m),float)

string_1 = str("HIERARCH ESO DRS CAL TH COEFF LL")

all_coeffs = []

mean = []
sigma = []
amp = []
area = []

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

all_waves= []

for order in all_coeffs:
	waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
	all_waves.append(waves)
	plt.plot(waves)	

count = 0
data_t = hdu[0].data

while count < 1000:
	start2 = time.time()
	max_index = np.argmax(lande_g)
	lande_g = np.delete(lande_g, [max_index])
	red = np.delete(red, [max_index])
	blue = np.delete(blue, [max_index])
	weight = np.delete(weight, [max_index])

	
	all_rv = []
	all_cc = []

	for j in range(0,72):
		x = all_waves[j]

		data = data_t[j]

		x1 = np.arange(x[0]-100., x[0], x[1]-x[0])
		x2 = np.arange(x[-1], x[-1]+100, x[-1]-x[-2])
		xtem = np.hstack([x1, x, x2])
    
		ftem = np.zeros(xtem.size)

	#Assigns correct fluxes to template
		for i in range(len(weight)):
			indices = np.where((xtem >= red[i]) & (xtem <= blue[i]))[0]
   			if indices.size>0:
    				ftem[indices] = weight[i]

		rv, cc = pyasl.crosscorrRV(x, data, xtem, ftem, -20., 20., 50./200, mode = "doppler")
		all_rv.append(rv) 
		all_cc.append(cc) 

	all_rv_array = np.array(all_rv)
	mean_index_rv = np.mean(all_rv_array, axis = 0)
	all_cc_array = np.array(all_cc)
	mean_index_cc = np.mean(all_cc_array, axis = 0)

	#fig = plt.figure()
	#ax1 = fig.add_subplot(111)
	#ax1.plot(mean_index_rv, mean_index_cc)

	#minind = np.argmin(mean_index_cc)

	#plt.plot(mean_index_rv[minind], mean_index_cc[minind], 'ro')
	#plt.xlabel("RV in km/s")

	#plt.show()

	xt = np.array(mean_index_rv)
	yt = np.array(mean_index_cc)

	#print xt, yt

	popt, pcov = curve_fit(gaussian, xt, yt, p0 = [-700000, 0, 10, 1700000])
	#print popt[1], popt[2]
	amp.append(popt[0])
	mean.append(popt[1])
	sigma.append(popt[2])


	offset = popt[3]
	area2 = np.trapz(yt, dx = 0.1)
	area.append(area2)

	print count, popt[1], popt[2]
	#print popt
	#print pcov
	#print sigma

	if count % 50 == 0:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(xt, yt, '.-')
		ynew = gaussian(xt, *popt)
		ax.plot(xt, ynew)
		fig.savefig('gaussian_%d.pdf' % count)
		plt.close(fig)


	if count % 10 == 0:
		np.savetxt("area.txt", area)
		print "Saved!"

	count += 1

	#print time.time() - start2

	#plt.plot(x, y, 'bo')
	#plt.plot(x, gaussian(x, *popt))
	#plt.show()

print mean
print sigma

print (time.time() - start)/60.

#oranges

