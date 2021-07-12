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



save_location= 'D:/Study Table/soumik_ceres/Report 6th sem/Results/ernutet_survey/'

Temperature_file = save_location + str("Temperature Profile.csv")

data = np.genfromtxt(Temperature_file, delimiter=',')

heat_map = sb.heatmap(data, square=True, xticklabels=False, yticklabels=False, vmin = 170, vmax=230)
                                                                                    ### plotting heatmap
plt.title('Temperature Map, Ernutet Crater, Ceres (Survey Data)')
plt.show()