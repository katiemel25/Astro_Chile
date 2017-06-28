import numpy as np
from astropy.io import ascii
import itertools

red = []
blue = []
weight = []

#reading in red cutoffs, blue cutoffs, and line weights from G2 mask

with open('G2_masks.txt') as masks:
	for line in masks:
		column = line.split()
		r = float(column[0])
		b = float(column[1])
		w = float(column[2])
		red.append(r)
		blue.append(b)
		weight.append(w)

#adding information line by line into separate arrays

red_array = np.array(red)
blue_array = np.array(blue)
weight_array = np.array(weight)

#using the arrays from red and blue cutoffs to average and find expected line centroid

centroid = (red_array + blue_array) / 2
centroid_array = np.array(centroid)

#stacked the two columns together into a 2D array

G2_info = np.column_stack((centroid_array, weight_array))

#print G2_info

#pulling info on element, wavelength, and lande g factor from unzipped VALD data

data = np.loadtxt("VALD_full.txt", usecols=[2,8]).transpose()
all_elem_name = np.genfromtxt("VALD_full.txt", dtype=str, usecols=(0))
all_elem_type = np.genfromtxt("VALD_full.txt", dtype=str, usecols=(1))

wavelength = data[0]
VALD_column_2 = data[1]

#print wavelength

#making array of just the VALD wavelengths that are closest to my centroid values
def find_nearest(array, value):
	idx = (np.abs(array - value)).argmin()
	return array[idx]

needed_values = []
for value in centroid_array:
	selected = find_nearest(wavelength, value)
	needed_values.append(selected)

needed_values_array = np.array(needed_values)

#needed_values gives the matched VALD wavelengths

#need to match selected wavelength data to corresponding lande g factors
#then will make a file with centroids, weights, VALD wavelengths, and 

lande = []
elem_name = []
elem_type = []

for value in needed_values_array:
	i = 0
	for line in wavelength:
		i = i + 1
		if value == line:
			lande.append(VALD_column_2[i-1])
			elem_name.append(all_elem_name[i-1])
			elem_type.append(all_elem_type[i-1])
			break

def remove_duplicates(values):
    output = []
    seen = set()
    for value in values:
        # If value has not been encountered yet,
        # ... add it to both list and set.
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output

#need to write all of these into a stacked array and then a txt file

#all_info = np.column_stack((elem_name, elem_type))

all_info = zip(elem_name, elem_type)

print all_info

#np.savetxt("element_info.txt", all_info , fmt="%s")

all_info_no_duplicates = remove_duplicates(all_info)

print all_info_no_duplicates

#np.savetxt("element_info_duplicates_removed.txt", all_info_no_duplicates , fmt="%s")

h = 4.135667662*10**(-15) #eVs
c = 299792458 #m/s

centroid_energies = []
for value in centroid_array:
	energy = h*c/(value*10**(-10))
	centroid_energies.append(energy)

VALD_energies = []
for value in needed_values_array:
	energy_2 = h*c/(value*10**(-10))
	VALD_energies.append(energy_2)


np.savetxt("elements_w_VALD_energies_and_wavlengths.txt", np.column_stack((all_info, needed_values_array, VALD_energies, lande)), fmt='%s')



