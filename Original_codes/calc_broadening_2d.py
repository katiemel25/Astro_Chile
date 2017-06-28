import pyfits, time
import numpy as np
from scipy.optimize import curve_fit
from scipy import interpolate
from PyAstronomy import pyasl
from astropy.io import ascii
from astropy.io import fits
import matplotlib
import matplotlib.backends
matplotlib.use('GTK')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#builds x and y
def values(h, j):# function that extracts the header wavelengths
    N = h['NAXIS' + str(j)];
    val = np.zeros(N);
    for i in range(0, N):
        #assigning values to the indices in the above created array of zeros
        val[i] = (i+1-float(h['CRPIX'+str(j)]))*float(h['CDELT'+str(j)])+float(h['CRVAL'+str(j)]);
    return val;


def compute_snr(x, data):
    ranges = [[5840.22, 5844.01], [5979.25, 5983.13], [6066.49, 6075.55], [6195.95, 6198.74], [6386.36, 6392.14], [5257.83, 5258.58], [5256.03, 5256.7]]
    sn = []
    sigma_e = []
    signal = []
    for r in range(len(ranges)):
        data_range = np.array([data[i] for i in range(len(data)) if (ranges[r][0] <= x[i] <= ranges[r][1]) == True])
        if len(data_range) != 0:
            m = np.mean(data_range)
            s = np.std(data_range)
            sn.append(m/s)
            sigma_e.append(s)
            signal.append(m)

    return np.mean(sn), np.mean(signal), np.mean(sigma_e)


def compute_vel(starname, x, data, xtem, ftem, make_plot_broadening = False):
    
    #calculating cross-correlation
    start = time.time()
    rv_temp, cc = pyasl.crosscorrRV(x, data, xtem, ftem, -40., 40., 50./200, mode = "doppler")

    return rv_temp, cc
    #print 'Finished cc for ' + starname + ' in ' + str(time.time() - start) + ' seconds'
    """
    #gives maximized cross-correlation function
    max_cc = max(cc)
    cc = cc/max_cc
    line_depth = 1. - cc
    
    ynew = line_depth
    xnew = rv_temp

    #computing Fourier Transform (meausure every possibility and
    #then give details about what the combined values consist of)

    sp = np.fft.fft(ynew)
    #returns transformed information on line depths in an array called "sp"
    freq = np.fft.fftfreq(len(xnew), d = (xnew[1] - xnew[0]))
    #frequency bins for the Fourier transformation parameters
    
    sp_final = np.sqrt(sp.real**2. + sp.imag**2.)
    #returns real and imaginary values of the Fourier transformed array
    #??? why is this necessary?

    i_pos = np.where(freq>0.0)
    #returns array of frequency bins only greater than zero
    sp_pos = sp_final[i_pos]
    #returns array of all fourier transformed line depths, indexed by bin number
    max_sp = max(sp_pos)
    #finding the maximum of all the line depths
    sp_pos = sp_pos/max_sp
    #??? does this relate all to the maximum value by kind of canceling it out?
    freq_pos = freq[i_pos]
    #returns array of fourier transformed frequencies, indexed by bin number too

    tck = interpolate.UnivariateSpline(freq_pos, sp_pos, s = 0, k = 4)
    #smoothes the data -- takes a noisy gaussian and adds moresamples for the x and y data
    zeros = tck.derivative().roots()

    new_vel = 0.66/zeros[0]
    
    
    if make_plot_broadening == True:
	   fig = plt.figure()
	   ax1 = fig.add_subplot(211)
	   ax1.plot(rv_temp, cc)
	   ax1.set_xlabel('Vel (km/s)')
	   ax1.set_ylabel('CC')
	
	   ax2 = fig.add_subplot(212)
	   ax2.plot(freq_pos, sp_pos, '.-')
	   ax2.plot(freq_pos, tck.__call__(freq_pos))
	   ax2.axvline(x = zeros[0], ls = '--')
	   ax2.set_xlabel('Frequency (s/km)')
	   ax2.set_ylabel('Power')
	   ax2.set_xscale('log')
	   ax2.set_yscale('log')
	
	   plt.tight_layout()
	   fig.savefig('./output/plots_broadening/' + starname + '.ps')
	   plt.close(fig)"""
    
    #return new_vel


def compute_error_broadening(sn):
    a = 4.414E-01
    b = 2.932E-02
    c = 3.2955E-03
    
    if sn < 150.:
	   return a*np.exp(-b*sn)
    else:
	   return c


def calc_broadening(f, mask, p, debug = False, file_debug = None, make_plot_broadening = False):
    hdu = pyfits.open(f)
    data = hdu[0].data
    header = hdu[0].header
    x = values(header, 1)
    #plots same stellar spectra as above
    
    x1 = np.arange(x[0]-100., x[0], x[1]-x[0])
    x2 = np.arange(x[-1], x[-1]+100, x[-1]-x[-2])
    xtem = np.hstack([x1, x, x2])
    
    c = 299792.458
    #speed of light
    
    lines1, lines2, flux_l = np.genfromtxt(mask, dtype=None, unpack=True)
    nlines = len(flux_l)
    
    ftem = np.zeros(xtem.size)
    
    for i in range(len(flux_l)):
	   indices = np.where((xtem >= lines1[i]) & (xtem <= lines2[i]))[0]
	   if indices.size>0:
	       ftem[indices] = flux_l[i]

    sn, mean_signal, mean_sigma = compute_snr(x, data)
	    
    broadening = compute_vel(f, x, data, xtem, ftem, make_plot_broadening)
    
    error_broadening = compute_error_broadening(sn)
    
    if debug == True:
	   file_debug.writelines('%d S/N=%f, broadening = %f +/- %f\n' % (p, sn, broadening, error_broadening))
    
    return broadening, error_broadening

"""
#Building the HARPS data plot
hdulist = fits.open("sun01_harps.fits")

data = hdulist[0].data #flux
hdu = hdulist[0].header

x = values(hdu, 1) #wavelength
print x

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x, data)
plt.show()"""

hdu = pyfits.open("HARPS.2005-09-03T21_28_39.386_e2ds_A.fits")
h = hdu[0].header
data = hdu[0].data
ss = hdu[0].shape
m = ss[0] #72
n = ss[1] #4096

dim = hdu[0].header['NAXIS2']

xarr = np.arange(n)
wave = np.zeros((n,m),float)

string_1 = str("HIERARCH ESO DRS CAL TH COEFF LL")

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

print all_coeffs

all_waves= []

for order in all_coeffs:
    waves = (order[3]*(xarr**3)) + (order[2]*(xarr**2)) + (order[1]*(xarr)) + order[0]
    all_waves.append(waves)
    plt.plot(waves)

#Plots both data and mask next to each other

mask = "G2_masks.txt"

for wave in all_waves:
        #hdu = pyfits.open(f)
        #data = hdu[0].data
        #header = hdu[0].header
        #x = values(wave, 1)
        #plots same stellar spectra as above
        x = wave

        #Extends template so that all data will fit    
        x1 = np.arange(x[0]-100., x[0], x[1]-x[0])
        x2 = np.arange(x[-1], x[-1]+100, x[-1]-x[-2])
        xtem = np.hstack([x1, x, x2])
    
        c = 299792.458
        #speed of light
    
        lines1, lines2, flux_l = np.genfromtxt(mask, dtype=None, unpack=True)
        nlines = len(flux_l)
    
        ftem = np.zeros(xtem.size)

        #Assigns correct fluxes to template
        for i in range(len(flux_l)):
            indices = np.where((xtem >= lines1[i]) & (xtem <= lines2[i]))[0]
            if indices.size>0:
                ftem[indices] = flux_l[i]

        """fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        ax1.plot(x, data)
        ax2.plot(xtem, ftem)"""
    

        rv, ccf = compute_vel("G2", x, data, xtem, ftem)
    
        """fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(rv, ccf)"""

        # Find the index of maximum cross-correlation function
    minind = np.argmin(ccf)
    plt.plot(rv[minind], ccf[minind], 'ro')
    plt.xlabel("RV in km/s")
    plt.show()


    
    """print("Cross-correlation function is maximized at dRV = ", rv[minind], " km/s for", mask)
    if rv[minind] > 0.0:
        print("  A red-shift with respect to the template")
    else:
        print("  A blue-shift with respect to the template")

    

    x = rv
    y = ccf

    n = len(x)
    mean = sum(x*y)/n
    sigma = sum(y*(x-mean)**2)/n

    def gaussian(x, a, b, sigma):
        val = a * np.exp(-(x-b)**2/sigma**2)
        return val

    popt, pcov = curve_fit(gaussian, x, y, p0 = [1, mean, sigma])

    plt.plot(x, gaussian(x, *popt), label = "Gaussian Fit")
    plt.plot(x, y, label = "Data")
    plt.legend
    plt.title("RV graph with Gaussian fit")
    plt.show()"""



"""
#DUMMY CODE FROM WEBSITE TO TEST

# Create the template
tw = np.linspace(5000,5010,1000)
tf = np.exp(-(tw-5004.0)**2/(2.*0.1**2))

# Create data, which are not that well sampled
dw = np.linspace(5000,5010,200)
df = np.exp(-(dw-5004.17)**2/(2.*0.1**2))

# Plot template and data
plt.title("Template (blue) and data (red)")
plt.plot(tw, tf, 'b.-')
plt.plot(dw, df, 'r.-')
plt.show()

# Carry out the cross-correlation.
# The RV-range is -30 - +30 km/s in steps of 0.6 km/s.
# The first and last 20 points of the data are skipped.
rv, cc = pyasl.crosscorrRV(dw, df, tw, tf, -30., 30., 30./50., skipedge=20)

# Find the index of maximum cross-correlation function
maxind = np.argmax(cc)

print("Cross-correlation function is maximized at dRV = ", rv[maxind], " km/s")
if rv[maxind] > 0.0:
  print("  A red-shift with respect to the template")
else:
  print("  A blue-shift with respect to the template")

plt.plot(rv, cc, 'bp-')
plt.plot(rv[maxind], cc[maxind], 'ro')
plt.show()"""