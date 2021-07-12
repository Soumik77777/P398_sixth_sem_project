import os
import csv
import numpy as np
import statistics
import scipy
import scipy.constants as const
import time
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sb
from spectral import *
from project_library import *


start = time.time()

#### Accessing image
location = 'D:/Study Table/soumik_ceres/data/'                      ### Location of Image

imagefile = location + str('survey-coarser/m-VIR_IR_1B_1_487085891_ernutet.hdr')

img = open_image(imagefile)                                         ### opening image
img_open = img.open_memmap(writeable = True)

wavelengths= img.bands.centers                                      ### wavelengths from metadata

wavelength_list = []
for i in wavelengths:
    if i<=4.34:                                                     ### creating wav list from 1-4.34
        wavelength_list.append(i)

artefact_corr_data = np.zeros((len(img_open), len(img_open[0]), len(wavelength_list)))


#### VIR Corr Data
VIR_corr_txt = 'VIR_correction_factor_IR.txt'
VIR_corr_data = read_csv(VIR_corr_txt, location=location, delimiter=',', poprows=1)


for row in range(len(img_open)):
    if (row+1)%10==0:
        print("currently on row", row+1, "of", str(len(img_open)))                  ### to keep track of time

    for col in range(len(img_open[0])):
        img_band = img_open[row, col, :]

        VIR_corr_value = []
        for ib in range(len(wavelength_list)):
            VIR_corr_value.append(img_band[ib] * VIR_corr_data[1][ib])              ### VIR correction

        uncorr_reflectance = pixel_validation(wavelength_list, VIR_corr_value)      ### Negative Error

        if uncorr_reflectance!=0:
            uncorr_reflectance2 = oddeven(wavelength_list, uncorr_reflectance)          ### Odd-Even effect

            for wav in range(len(wavelength_list)):
                artefact_corr_data[row, col, wav] = uncorr_reflectance2[wav]              # Appending value



save_location = 'D:/Study Table/soumik_ceres/Report 6th sem/Results/ernutet_survey/'

envi.save_image(save_location+str("artifact_corrected.hdr"), artefact_corr_data, metadata = img.metadata, force = True)


print('Total time taken:', int((time.time()-start)/60), 'min', int((time.time()-start)%60), 'sec')


