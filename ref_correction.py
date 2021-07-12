import os
import csv
import time
import numpy as np
import scipy
from scipy.optimize import curve_fit
import seaborn as sb
import matplotlib.pyplot as plt
from spectral import *
from project_library import *


start = time.time()                                                 ### To check time

location = 'D:/Study Table/soumik_ceres/data/'                      ### Location of Image

ss_txtname = 'ss-ceres-dawn.txt'                                    ### Solar flux data
ss_data= read_csv(ss_txtname, location=location, delimiter='\t', poprows=1)

location = 'D:/Study Table/soumik_ceres/Report 6th sem/Results/ernutet_survey/'
imagefile = location + str('artifact_corrected.hdr')

img = open_image(imagefile)                                         ### opening image
img_open = img.open_memmap(writeable = True)

wavelengths= img.bands.centers                                      ### wavelengths from metadata

Temperature = []

corrected_data = np.zeros((len(img_open), len(img_open[0]), (len(wavelengths) - 11)))

for row in range(len(img_open)):
    if (row+1)%10==0:
        print("currently on row", row+1, "of", str(len(img_open)))                  ### to keep track of time
    Temp_col = []

    for col in range(len(img_open[0])):
        img_band = img_open[row, col, :]

        wavelength_list, uncorr_reflectance = [], []
        for val in range(len(img_band)):
            wavelength_list.append(wavelengths[val])
            uncorr_reflectance.append(img_band[val])

        if any(i<=0 for i in uncorr_reflectance):
            Temp_col.append(0)
        else:
            wa_range1, uncorr_ref_range1 = wavelength_range(wavelength_list, uncorr_reflectance, 1.7, 2.5)
            m, c = st_line_fit(wa_range1, uncorr_ref_range1, fig=None)  # st line fit in no absorpton range

            ###### EMISSIVITY ###########
            emissivity = 1 - np.mean(uncorr_ref_range1)

            straight_line, thermal_list = [], []
            for n in range(len(wavelength_list)):
                straight_line.append(m[0] * wavelength_list[n] + c[0])
                thermal_list.append(uncorr_reflectance[n] - straight_line[-1])

            temp, ceres_thermal_emission, ceres_corrected = \
                temp_fitting(wavelength_list, uncorr_reflectance, straight_line, thermal_list, emissivity, ss_data)

            Temp_col.append(round(temp, 4))

            for th in range(len(wavelength_list)):
                corrected_data[row, col, th] = ceres_corrected[th]

    Temperature.append(Temp_col)

print("Iteration completed.")



save_location= 'D:/Study Table/soumik_ceres/Report 6th sem/Results/ernutet_survey/'

Temperature_file = save_location + str("Temperature Profile.csv")
with open(Temperature_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(Temperature)

envi.save_image(save_location+str("corrected_reflectance.hdr"), corrected_data, metadata = img.metadata, force = True)

print('Corrected HDR file saved')


print('Total time taken:', int((time.time()-start)/60), 'min', int((time.time()-start)%60), 'sec')

