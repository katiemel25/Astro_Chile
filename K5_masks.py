import numpy as np
from astropy.io import ascii

red = []
blue = []
weight = []

#reading in red cutoffs, blue cutoffs, and line weights from G2 mask

with open('K5_masks.txt') as masks:
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

for value in needed_values_array:
	i = 0
	for line in wavelength:
		if value == line:
			lande.append(VALD_column_2[i])
			break
		else:
			i = i + 1

lande_array = np.array(lande)

#need to write all of these into a stacked array and then a txt file

all_info = np.column_stack((centroid_array, weight_array, needed_values_array, lande_array))

print all_info

np.savetxt("all_K5_data.txt", all_info)