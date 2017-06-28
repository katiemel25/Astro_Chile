import numpy as np

#h alpha test

alpha = 656.281*10**(-9) #meters

h = 4.135667662*10**(-15) #eVs
c = 299792458 #meters/sec

h_alpha_energy = h*c/alpha
h_alpha_wavenumber = 1/(alpha*100)

#print h_alpha_energy
#print h_alpha_wavenumber

#test with neon

iron = 4045.8119*10**(-8) #centimeters

iron_wavenumber = 1/(iron)

print iron_wavenumber

