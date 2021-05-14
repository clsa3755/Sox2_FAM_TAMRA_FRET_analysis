#FRET Analysis for SOX2 bound to DNA
#Fluorophores used for example data are FAM and TAMRA excited at 490 nm and 560 nm respectively. Constants and Forster Radius calculated for this experiment are defined where applicable.
#To execute donor (FAM) and acceptor (TAMRA) emission data will need to be in separate csv files with wavelengths (nm) in the left most column, the length of Sox2 (nm) in the upper most row and the corresponding data between.
#There will need to be an additional csv file of the emission intensities by wavelength of just the donor(FAM) to later subtract out the background
#Tasks to execute:
##Read donor and accpetor emission spectra and output plots
##normalize donor emission spectra and construct new array to plot
##calculate FRET Effect, FRET Efficiency, distance between fluorophores, and DNA bend angle
##plot DNA bend angle over length of Sox2
##output DNA ben angle in degrees for each length of Sox2

#import the following:
import numpy as np
import os
import matplotlib.pyplot as plt
import math

#open donor emission spectra file as donor_e_file
#open acceptor emission spectra file as acceptor_e_file
#open donor only spectra data as donor_only_file
donor_e_file = os.path.abspath('sox2_fret_donor_emission.csv')
acceptor_e_file = os.path.abspath('sox2_fret_acceptor_emission.csv')
donor_only_file = os.path.abspath('FAM_sox2_donor_only.csv')

#******extract data needed******
#extract donor emission spectra data as donor_e_data
donor_e_data = np.genfromtxt(fname=donor_e_file,delimiter = ',', dtype='unicode')
#extract acceptor emission spectra data as acceptor_e_data
acceptor_e_data = np.genfromtxt(fname=acceptor_e_file,delimiter = ',', dtype='unicode')
#extract donor only emission data as donor_only_data
donor_only_data_str = np.genfromtxt(fname=donor_only_file,skip_header=1,delimiter = ',', dtype='unicode')


#extract the header of the donor data:varying length of SOX2 in nm as donor_e_sox2
#repeat for acceptor data as acceptor_e_sox2
donor_e_sox2=donor_e_data[0,1:]
acceptor_e_sox2=acceptor_e_data[0,1:]

#extract wavelengths for donor emission spectra as donor_e_wavelengths_str and convert to float(donor_e_wavelengths)
#repeat for acceptor wavelengths
donor_e_wavelengths_str = donor_e_data[1:,0]
donor_e_wavelengths = donor_e_wavelengths_str.astype(np.float)

acceptor_e_wavelengths_str = acceptor_e_data[1:,0]
acceptor_e_wavelengths = acceptor_e_wavelengths_str.astype(np.float)


#extract intensities at all SOX2 lengths for donor as donor_intensities_str and convert to float(donor_intensities)
#repeat for acceptor intensities
donor_intensities_str = donor_e_data[1:,1:]
donor_intensities = donor_intensities_str.astype(np.float)

acceptor_intensities_str = acceptor_e_data[1:,1:]
acceptor_intensities = acceptor_intensities_str.astype(np.float)

#convert donor only data to float
donor_only_data = donor_only_data_str.astype(np.float)

#extract intensities for each wavelength from the donor_only_data
donor_only_intensities = donor_only_data[:,1]

#use length function of donor_e_sox2 to determine number of columns
num_columns = len(donor_e_sox2)
#use length of donor_e_wavelengths to find number of rows
num_rows = len(donor_e_wavelengths)


#******Normalize donor emission spectra*******
#find all donor emission intensities at 520 nm which was determined to be the first peak
location_of_first_peak_wavelength = np.where(donor_e_wavelengths == 520)
print(location_of_first_peak_wavelength)
#use print output to find row number of intensities
#extract all intensities from that wavelength as peak1_donor_intensities
peak1_donor_intensities = donor_intensities[10,0:]

#in donor only emission spectra data, find intensity at 520 nm
location_of_donor_only_wavelength = np.where(donor_only_data == 520)
print(location_of_donor_only_wavelength)
#use print output to find row number of intensity
#extract that intensity value as donor_only_intensity
donor_only_intensity= donor_only_data[10,1]


#divide each intensity of the 520nm (peak1) donor excitation intensities by the donor only intensity at 520nm to get a corresponding multiple for each length of Sox2
#append this to a new array 'FAM_sox2_fit_multiple'
FAM_sox2_fit_multiple=[]
FAM_sox2_fit_intensities = []
for intensity in peak1_donor_intensities:
    multiple = intensity/donor_only_intensity
    FAM_sox2_fit_multiple.append(multiple)

#to get FAM fit arrays, run loop through the array of multiples for each corresponsing SOX2 length
#define 'multiples' as every multiple in array from 0 to number of columns
#run a second loop within that loop that goes through each of the donor only intensities at every wavelength
#define intensities as every donor only intensity from 0 to number of rows
#multiply the multiples and intensities together to get the fit donor intensities as fit_intensities
#append to array 'FAM_sox2_fit_intensities' should be 2D array of all intensities fit with each multiple for the corresponding length of sox2
for multiples in FAM_sox2_fit_multiple:
    mutiples = FAM_sox2_fit_multiple[0:num_columns]
    for intensities in donor_only_intensities:
        intensities = donor_only_intensities[0:num_rows]
        fit_intensities = intensities*multiples
    FAM_sox2_fit_intensities.append(fit_intensities)


#create new array for normalizing the donor intensity values in which donor intensities are in a 2D array separated by columns
#separate donor intensities by column through a range function
donor_values=[]
for j in range(0,num_columns):
    values_donor_columns = donor_intensities[:,j]
    donor_values.append(values_donor_columns)

#subtract the FAM_sox2_fit_intensities from the donor_values array to get the final array of 'normalized_intensities'
normalized_intensities = np.subtract(donor_values,FAM_sox2_fit_intensities)

#******begin calculations*******
#all calculations are representative of a corresponding length of SOX2 in nm (equal to array donor_e_sox2)
#*FRET Effect*
#extract intensities at 580 nm (second peak of donor emission spectra) for normalized donor values and acceptor values
#define as row_donor and row_acceptor
row_donor=np.where(donor_e_wavelengths == 580)
row_acceptor=np.where(acceptor_e_wavelengths == 580)
print(row_donor)
print(row_acceptor)
#use print output to extract intensities at 580 nm as FE_intensities_donor and FE_intensities_acceptor
FE_intensities_donor = normalized_intensities[:,70]
FE_intensities_acceptor = acceptor_intensities[10,:]
#divide FE_intensities_donor by FE_intensities_acceptor
Fret_Effect = np.divide(FE_intensities_donor,FE_intensities_acceptor)

#make array of zeroed FRET Effect by subtracting each value by the zeroth element of the Fret_Effect array
#define array as zeroed_FE
zeroed_FE = [Fret_Effect[n] - Fret_Effect[0] for n in range(len(Fret_Effect))]

#*E=Fret Efficiency*
#E = FE[εFAM490/εTAMRA560] + εTAMRA490/εTAMRA560
#εFAM490/εTAMRA560 and εTAMRA490/εTAMRA560 are constants and were calculated to be 0.543 and 0.121
#loop through array'Fret_Effect' and plug each 'effect' value into the equation
#append the results to an array of the corresponding Fret Efficiencies defined as 'E_array'
E_array = []
for effect in Fret_Effect:
    E = (effect/0.543) - 0.121
    E_array.append(E)

#make array of zeroed E values by subtracting each value by the zeroth element of E_array
#define array as zeroed_E
zeroed_E = [E_array[n] - E_array[0] for n in range(len(E_array))]


#*Rda (A) or 'r' = distance between donor fluorophore and acceptor fluorophore*
#r is equal to (((1/E-1)^(1/6))*Forster radius
#Forster radius = 60.5 and was calculated previously for this experiment
#loop through all 'E_value' in 'E_array' and plug each into the equation
#append resulting values of r to 'r_array'
r_array = []
for E_value in E_array:
        r = (((1/E_value)-1)**(1/6))*60.5
        r_array.append(r)

#*delta r by length of SOX2*
#calculate the change in r for each length of Sox2 by subtracting each element of the array from the zeroth element
#the resulting array is defined as delta_r
delta_r = [r_array[0] - r_array[n] for n in range(len(r_array))]

#*DNA bend angle in radians*
#the DNA bend angle is equal to arccos((r/2)/(Forster radius/2))*2
#calculate DNA bend angle for each length of sox2 by looping through every 'r_value' in 'r_array'
#append resulting bend_angle_rad values to bend_angle_rad_array
bend_angle_rad_array = []
for r_value in r_array:
    bend_angle_rad = np.arccos((r_value/2)/30.25)*2
    bend_angle_rad_array.append(bend_angle_rad)


#*DNA bend angle in degrees*
#covert all 'radian' values in 'bend_angle_rad_array' to degrees by looping through the array and using np.degrees
#append resulting 'bend_angle_deg' to 'bend_angle_deg_array'
#print array to get output of bend angles for each varying length of Sox2 bound to DNA
bend_angle_deg_array=[]
for radian in bend_angle_rad_array:
    bend_angle_deg = np.degrees(radian)
    bend_angle_deg_array.append(bend_angle_deg)
print(bend_angle_deg_array)

#Prepare subplots with 4 subplots in one column
fig,ax = plt.subplots(5,1, figsize = (10,20))
fig.subplots_adjust(hspace=.5)

#Run a for loop to graph all intensities for the donor emission spectra (defined as 'donor_y_axis') against all wavelengths(defined as 'donor_x_axis')
#Repeat for the acceptor emission spectra with intensities defined as 'acceptor_y_axis' and wavelengths defined as 'acceptor_x_axis'
for i in range(0,num_columns):
    donor_x_axis = donor_e_wavelengths
    donor_y_axis = donor_intensities[:,i]
    acceptor_x_axis = acceptor_e_wavelengths
    acceptor_y_axis = acceptor_intensities[:,i]
    ax[0].plot(donor_x_axis,donor_y_axis)
    ax[1].plot(acceptor_x_axis,acceptor_y_axis)
    ax[0].set_xlabel("Wavelength (nm)")
    ax[0].set_ylabel("Intensities (a.u.)")
    ax[0].title.set_text('ex 490 nm Donor Emission')
    ax[1].set_xlabel("Wavelength (nm)")
    ax[1].set_ylabel("Intensities (a.u.)")
    ax[1].title.set_text('ex 560 nm Acceptor Emission')

#create for loop to get y axes for the FAM fit emission spectra
#loop through the 2D array FAM_sox2_fit_intensities from 0 to number of columns and define each value in the array as 'intensity_arrays'
#define y axes as each of the intensity_arrays from 0 to number of rows
#the x axis should be the same for all donor emission spectra
#repeat this same for loop process to graph the normalized emission spectra using the variable 'normalized_y_axis' as 'normalized_arrays' from 0 to number of rows
for intensity_arrays in FAM_sox2_fit_intensities[0:num_columns]:
    FAM_sox2_fit_y_axis = intensity_arrays[0:num_rows]
    ax[2].plot(donor_x_axis,FAM_sox2_fit_y_axis)
    ax[2].set_xlabel("Wavelength (nm)")
    ax[2].set_ylabel("Intensities (a.u.)")
    ax[2].title.set_text('FAM SOX2 Fit')

for normalized_arrays in normalized_intensities[0:num_columns]:
    normalized_y_axis = normalized_arrays[0:num_rows]
    ax[3].plot(donor_x_axis,normalized_y_axis)
    ax[3].set_xlabel("Wavelength (nm)")
    ax[3].set_ylabel("Intensities (a.u.)")
    ax[3].title.set_text('NSTAM Extracted')

#plot bend angle of DNA vs length of Sox2
#'bend_angle_y_axis' as 'bend_angle_rad_array' without the zero value and 'bend_angle_x_axis' as 'donor_e_sox' without the zero value

bend_angle_x_axis = donor_e_sox2[1:]
bend_angle_y_axis = bend_angle_rad_array[1:]
ax[4].plot(bend_angle_x_axis,bend_angle_y_axis)
ax[4].set_xlabel("Sox2 (nm)")
ax[4].set_ylabel("DNA Bend Angle (radians)")

plt.savefig('Sox2_FAM_TAMRA.png',dpi=300)
