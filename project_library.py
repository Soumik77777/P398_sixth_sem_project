'''
Soumik Bhattacharyya
3rd year Int. M.Sc.
SPS, NISER

Functions used for:
Removing Thermal Component from Spectroscopy Data of Asteroid Ceres
'''


import os
import csv
import numpy as np
import statistics
import scipy
import scipy.constants as const
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as mpl_patches
from spectral import *


###############################################################

def read_csv(filename, location=None, poprows=None, delimiter=None):
    if location==0:
        filepath = filename
    elif location!= None and location!=0:
        filepath = str(location) + str(filename)
    else:
        filepath = str(os.getcwd()) + str('\\') + str(filename)

    if delimiter=='\t':
        delim = '\t'
    else:
        delim = ','

    with open(filepath, 'r') as infile:
        data = csv.reader(infile, delimiter=delim)
        datalist, rows = [], []
        for row in data:
            datalist.append(row)
        if poprows!=None:
            for i in range(poprows):
                datalist.pop(0)
        for j in range(len(datalist[0])):
            globals()['string%s' % j] = []
            for k in datalist:
                globals()['string%s' % j].append(float(k[j]))
            rows.append(globals()['string%s' % j])
        infile.close()

    return rows


def create_csv(filename, xlist, ylist, location=None, xheader=None, yheader=None):
    if location==0:
        filepath = filename
    elif location!= None and location!=0:
        filepath = str(location) + str(filename) + str('.csv')
    else:
        filepath = str(os.getcwd()) + str('\\') + str(filename) + str('.csv')

    if len(xlist) == len(ylist):
        with open(filepath, 'w', newline='') as file:
            writer = csv.writer(file)
            if xheader != None and yheader != None:
                writer.writerow([xheader, yheader])
            for k in range(len(xlist)):
                writer.writerow([xlist[k], ylist[k]])
    else:
        print("Length of xlist and ylist are not equal.")

    return 0


def blackbody_rad(wa,T):
    c1 = 3.74177*(10**-16)          #2*const.pi*const.h*(const.c**2)            #3.7417718521927573e-16
    c2 = 0.014387774                #const.c*const.h/const.Boltzmann            #0.014387768775039337

    radiance= ((c1/((wa/(10**6))**5))*(1/(np.exp(c2/((wa*T)/(10**6)))-1)))/((10**6)*np.pi)
    return radiance


def wavelength_range(xlist, ylist, l, u):
    x_in_range, y_in_range = [], []
    for i in range(len(xlist)):
        if l<xlist[i]<u:
            x_in_range.append(xlist[i])
            y_in_range.append(ylist[i])

    return x_in_range, y_in_range



def pixel_validation(wa, valuelist):
    neg_index_list = []
    for i in range(len(valuelist)):
        if valuelist[i]<=0:
            neg_index_list.append(i)                # makes an index list of all the negative values
    if len(neg_index_list)>=11:                     # rejects if more than 10 such values exist
        return 0
    elif any(2.5<wa[index]<=3.5 for index in range(len(neg_index_list))):
        return 0                                    # rejects if negative values are presnt in main absorption range
    else:
        for j in range(len(neg_index_list)):
            if int(neg_index_list[j]+1) in neg_index_list:
                if int(neg_index_list[j]+2) in neg_index_list:
                    if int(neg_index_list[j]+3) in neg_index_list:  # rejects if more than three consecutive values are neg
                        valuelist = 0
                        break
                    else:                           # correction if three consecutive entries are neg
                        valuelist = point_corr3(wa, valuelist, neg_index_list[j])
                else:                               # correction if two consecutive entries are neg
                    valuelist = point_corr2(wa, valuelist, neg_index_list[j])
            else:                                   # correction for single point
                valuelist = point_corr1(wa, valuelist, neg_index_list[j])

        return valuelist


def remove_neg(wav_list, value_list):  # if there are negative values even in the values present in
    # preceding and succseeding cells, it removes them

    neg_index_list = []
    if len(wav_list) == len(value_list):
        for i in range(len(wav_list)):
            if value_list[i] <= 0:
                neg_index_list.append(i)
    for ind in range(len(neg_index_list)):
        wav_list.pop(ind)
        value_list.pop(ind)
    return wav_list, value_list


def point_corr3(wa, valuelist, index):
    if index==0:
        wa_list = [wa[int(index)+3], wa[int(index)+4], wa[int(index)+5], wa[int(index)+6]]
        values = [valuelist[int(index)+3], valuelist[int(index)+4], valuelist[int(index)+5], valuelist[int(index)+6]]
    elif index==1:
        wa_list = [wa[int(index)-1], wa[int(index)+3], wa[int(index)+4], wa[int(index)+5]]
        values = [valuelist[int(index)-1], valuelist[int(index)+3], valuelist[int(index)+4], valuelist[int(index)+5]]
    elif index==len(valuelist)-3:
        wa_list = [wa[int(index)-4], wa[int(index)-3], wa[int(index)-2], wa[int(index)-1]]
        values = [valuelist[int(index)-4], valuelist[int(index)-3], valuelist[int(index)-2], valuelist[int(index)-1]]
    elif index==len(valuelist)-4:
        wa_list = [wa[int(index)-3], wa[int(index)-2], wa[int(index)-1], wa[int(index)+3]]
        values = [valuelist[int(index)-3], valuelist[int(index)-2], valuelist[int(index)-1], valuelist[int(index)+3]]
    else:
        wa_list = [wa[int(index)-2], wa[int(index)-1], wa[int(index)+3], wa[int(index)+4]]
        values = [valuelist[int(index)-2], valuelist[int(index)-1], valuelist[int(index)+3], valuelist[int(index)+4]]

    wa_list, values = remove_neg(wa_list, values)

    m, c = st_line_fit(wa_list, values, fig=0)

    valuelist[index] = m[0] * wa[index] + c[0]
    valuelist[index + 1] = m[0] * wa[index + 1] + c[0]
    valuelist[index + 2] = m[0] * wa[index + 2] + c[0]

    return valuelist


def point_corr2(wa, valuelist, index):
    if valuelist[index]>0 and valuelist[index + 1]>0:
        return valuelist
    else:
        if index == 0:
            wa_list = [wa[int(index) + 2], wa[int(index) + 3], wa[int(index) + 4], wa[int(index) + 5]]
            values = [valuelist[int(index) + 2], valuelist[int(index) + 3], valuelist[int(index) + 4],
                      valuelist[int(index) + 5]]
        elif index == 1:
            wa_list = [wa[int(index) - 1], wa[int(index) + 2], wa[int(index) + 3], wa[int(index) + 4]]
            values = [valuelist[int(index) - 1], valuelist[int(index) + 2], valuelist[int(index) + 3],
                      valuelist[int(index) + 4]]
        elif index == len(valuelist) - 2:
            wa_list = [wa[int(index) - 4], wa[int(index) - 3], wa[int(index) - 2], wa[int(index) - 1]]
            values = [valuelist[int(index) - 4], valuelist[int(index) - 3], valuelist[int(index) - 2],
                      valuelist[int(index) - 1]]
        elif index == len(valuelist) - 3:
            wa_list = [wa[int(index) - 3], wa[int(index) - 2], wa[int(index) - 1], wa[int(index) + 2]]
            values = [valuelist[int(index) - 3], valuelist[int(index) - 2], valuelist[int(index) - 1],
                      valuelist[int(index) + 2]]
        else:
            wa_list = [wa[int(index) - 2], wa[int(index) - 1], wa[int(index) + 2], wa[int(index) + 3]]
            values = [valuelist[int(index) - 2], valuelist[int(index) - 1], valuelist[int(index) + 2],
                      valuelist[int(index) + 3]]

        wa_list, values = remove_neg(wa_list, values)

        m, c = st_line_fit(wa_list, values, fig=0)

        valuelist[index] = m[0] * wa[index] + c[0]
        valuelist[index+1] = m[0] * wa[index+1] + c[0]

        return valuelist



def point_corr1(wa, valuelist, index):
    if valuelist[index]>0:
        return valuelist
    else:
        if index == 0:
            wa_list = [wa[int(index) + 2], wa[int(index) + 3], wa[int(index) + 4], wa[int(index) + 5]]
            values = [valuelist[int(index) + 2], valuelist[int(index) + 3], valuelist[int(index) + 4],
                      valuelist[int(index) + 5]]
        elif index == 1:
            wa_list = [wa[int(index) - 1], wa[int(index) + 1], wa[int(index) + 2], wa[int(index) + 3]]
            values = [valuelist[int(index) - 1], valuelist[int(index) + 1], valuelist[int(index) + 2],
                      valuelist[int(index) + 3]]
        elif index == len(valuelist) - 1:
            wa_list = [wa[int(index) - 4], wa[int(index) - 3], wa[int(index) - 2], wa[int(index) - 1]]
            values = [valuelist[int(index) - 4], valuelist[int(index) - 3], valuelist[int(index) - 2],
                      valuelist[int(index) - 1]]
        elif index == len(valuelist) - 2:
            wa_list = [wa[int(index) - 3], wa[int(index) - 2], wa[int(index) - 1], wa[int(index) + 1]]
            values = [valuelist[int(index) - 3], valuelist[int(index) - 2], valuelist[int(index) - 1],
                      valuelist[int(index) + 1]]
        else:
            wa_list = [wa[int(index) - 2], wa[int(index) - 1], wa[int(index) + 1], wa[int(index) + 2]]
            values = [valuelist[int(index) - 2], valuelist[int(index) - 1], valuelist[int(index) + 1],
                      valuelist[int(index) + 2]]

        wa_list, values = remove_neg(wa_list, values)

        m, c = st_line_fit(wa_list, values, fig=0)
        valuelist[index] = m[0] * wa[index] + c[0]

        return valuelist



def weighted_avg(x1, x2, y1, y2, x3):
    m = (y2 - y1) / (x2 - x1)
    c = y2 - (m * x2)
    y3 = (m * x3) + c
    return y3


def oddeven(wav_list, valuelist):
    wa1, un1 = wavelength_range(wav_list, valuelist, 1, 2.5)

    un11 = un1.copy()
    un12 = un1.copy()

    for i in range(len(un1)):
        if i % 2 == 0:
            if i != 0 and i != len(un1) - 1:
                un11[i] = weighted_avg(wa1[i - 1], wa1[i + 1], un1[i - 1], un1[i + 1], wa1[i])
        else:
            if i != len(un1) - 1:
                un12[i] = weighted_avg(wa1[i - 1], wa1[i + 1], un1[i - 1], un1[i + 1], wa1[i])

    un2 = []
    for k in range(len(un1)):
        un2.append((un11[k] + un12[k]) / 2)

    for j in range(len(un1)):
        valuelist[j] = un2[j]

    return valuelist



def temp_fitting(wav_list, uncorr_list, st_line, thermal_list, emissivity, ss_data):

    def ceres_thermal_part(wa, T):
        d = 2.8
        ceres_rad = (d ** 2) * blackbody_rad(wa, T) * emissivity
        return ceres_rad

    def sd(wav_list, corr_list, st_line, min=wav_list[0], max=wav_list[-1]):
        wa_range, corr_range = wavelength_range(wav_list, corr_list, min, max)
        wa_range, st_range = wavelength_range(wav_list, st_line, min, max)
        value=0
        for wa in range(len(wa_range)):
            value+= abs(corr_range[wa] - st_range[wa])
        return value

    def least_element(list):
        least_index, least_item = 0, list[0]
        for i in range(len(list)-1):
            if list[i+1]<least_item:
                least_index = i+1
                least_item = list[i+1]
        return least_index

    popt, error = func_fit(ceres_thermal_part, wav_list, thermal_list, p0=210, fig=0)

    init_temp = round(popt[0])
    init_temp_list = np.linspace(init_temp - 15, init_temp + 35, 11)
    temp_list = init_temp_list

    for loop in range(4):
        final_temp_list = temp_list.copy()
        daviation = []
        for temp in temp_list:
            ceres_emission, corrected_list = [], []
            for th in range(len(wav_list)):
                ceres_emission.append(ceres_thermal_part(wav_list[th], temp) * np.pi / ss_data[1][th])
                corrected_list.append(uncorr_list[th] - ceres_emission[-1])
            daviation.append(sd(wav_list, corrected_list, st_line, min=4, max=4.34))

        index = least_element(daviation)
        if index==0:
            lower_limit, upper_limit = temp_list[index], temp_list[index + 1]
        elif index==len(daviation)-1:
            lower_limit, upper_limit = temp_list[index - 1], temp_list[index]
        else:
            lower_limit, upper_limit = temp_list[index - 1], temp_list[index + 1]

        temp_list = np.linspace(lower_limit, upper_limit, 11)

    final_temperature = final_temp_list[index]

    ceres_emission, corrected_list = [], []
    for wav in range(len(wav_list)):
        ceres_emission.append(ceres_thermal_part(wav_list[wav], final_temperature) * np.pi / ss_data[1][wav])
        corrected_list.append(uncorr_list[wav] - ceres_emission[-1])

    return final_temperature, ceres_emission, corrected_list




def st_line_fit(xdata, ydata, fig=0, xlabel=None, ylabel=None, l=None, u=None,
                data_color='b', data_marker='o', data_ms='3', data_label='datapoints',
                fit_color='r', fit_ls='-', fit_lw=1, fit_label='Fitted Curve'
                ):

    def st_line(x,m,c):
        return m*x +c

    popt, pcov = curve_fit(st_line, xdata, ydata)
    error = np.sqrt(np.diag(pcov))
    m=[popt[0],error[0]]
    c=[popt[1],error[1]]

    if l!=None:
        l=l
    else:
        l=xdata[0]

    if u!=None:
        u=u
    else:
        u=xdata[-1]

    def fitted_plot(xdata, ydata, xlabel=xlabel, ylabel=ylabel,
                    data_color=data_color, data_marker=data_marker, data_ms=data_ms, data_label=data_label,
                    popt=popt, error=error, l=l, u=u,
                    fit_color=fit_color, fit_ls=fit_ls, fit_lw=fit_lw, fit_label=fit_label
                    ):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata, ydata, color=data_color, marker=data_marker, markersize=data_ms, label=data_label)

        slope = 'Slope(m)= ' + str(round(popt[0], 4)) + ' +/- ' + str(round(error[0], 4))
        intercept = str('Intercept(c)= ') + str(round(popt[1], 4)) + str(' +/- ') + str(round(error[1], 4))

        xfit = np.linspace(l, u, 100)
        ax.plot(xfit, st_line(xfit, *popt), color=fit_color, ls=fit_ls, lw=fit_lw, label=fit_label)

        # plt.axvline(x=0, color='k', linewidth=1)                    ###########
        # plt.axhline(y=0, color='k', linewidth=1)                    ###########
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(b=None, which='major', axis='both')
        leg1 = ax.legend(loc=4, bbox_to_anchor=(1, 0))  # bbox_to_anchor

        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 2
        labels = []
        labels.append(slope)
        labels.append(intercept)
        leg2 = ax.legend(handles, labels, loc=2, bbox_to_anchor=(0, 1), fontsize='small', fancybox=True, framealpha=0.7,
                         handlelength=0, handletextpad=0)
        ax.add_artist(leg1)

        plt.show()

        return 0

    if fig=='yes' or fig=='y' or fig==1:
        fitted_plot(xdata, ydata, xlabel=xlabel, ylabel=ylabel,
                    data_color=data_color, data_marker=data_marker, data_ms=data_ms, data_label=data_label,
                    popt=popt, error=error, l=l, u=u,
                    fit_color=fit_color, fit_ls=fit_ls, fit_lw=fit_lw, fit_label=fit_label
                    )
        return m,c
    elif fig==None or fig==0:
        return m,c
    else:
        print("Do you want to plot the fitted curve?")
        fig = str(input("Type Y for yes and N for no:"))

        if fig=='N':
            return m,c
        elif fig=='Y':
            fitted_plot(xdata, ydata, xlabel=xlabel, ylabel=ylabel,
                        data_color=data_color, data_marker=data_marker, data_ms=data_ms, data_label=data_label,
                        popt=popt, error=error, l=l, u=u,
                        fit_color=fit_color, fit_ls=fit_ls, fit_lw=fit_lw, fit_label=fit_label
                        )
            return m, c
        else:
            print("Please type 'Y' or define fig='yes' for Figure.")
            return m, c



def func_fit(func, xdata, ydata, params=None, p0=None, fig='ask', xlabel=None, ylabel=None, l=None, u=None,
                data_color='b', data_marker='o', data_ms='3', data_label='datapoints',
                fit_color='r', fit_ls='-', fit_lw=1, fit_label='Fitted Curve'
                ):
                    #def func(x, m, c):
                    #   return m*x + c

    popt, pcov = curve_fit(func, xdata, ydata, p0=p0)
    error = np.sqrt(np.diag(pcov))

    if l!=None:
        l=l
    else:
        l=xdata[0]

    if u!=None:
        u=u
    else:
        u=xdata[-1]

    if params==None:
        params=[]
        for i in range(len(popt)):
            params.append('Parametre'+str(i+1))

    def fitted_plot(func, xdata, ydata, params=params, xlabel=xlabel, ylabel=ylabel,
                    data_color=data_color, data_marker=data_marker, data_ms=data_ms, data_label=data_label,
                    popt=popt, error=error, l=l, u=u,
                    fit_color=fit_color, fit_ls=fit_ls, fit_lw=fit_lw, fit_label=fit_label
                    ):

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata, ydata, color=data_color, marker=data_marker, markersize=data_ms, label=data_label)

        xfit = np.linspace(l, u, 100)
        ax.plot(xfit, func(xfit, *popt), color=fit_color, ls=fit_ls, lw=fit_lw, label=fit_label)

        # plt.axvline(x=0, color='k', linewidth=1)                    ###########
        # plt.axhline(y=0, color='k', linewidth=1)                    ###########
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(b=None, which='major', axis='both')
        leg1 = ax.legend(loc=4, bbox_to_anchor=(1, 0))  # bbox_to_anchor

        handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)] * 2
        labels = []
        for i in range(len(popt)):
            labels.append(params[i]+'= '+str(round(popt[i], 4)) + ' +/- ' + str(round(error[i], 4)))
        leg2 = ax.legend(handles, labels, loc=2, bbox_to_anchor=(0, 1), fontsize='small', fancybox=True, framealpha=0.7,
                         handlelength=0, handletextpad=0)
        ax.add_artist(leg1)

        plt.show()

        return 0

    if fig=='yes' or fig=='y' or fig==1:
        fitted_plot(func, xdata, ydata, params=params, xlabel=xlabel, ylabel=ylabel,
                    data_color=data_color, data_marker=data_marker, data_ms=data_ms, data_label=data_label,
                    popt=popt, error=error, l=l, u=u,
                    fit_color=fit_color, fit_ls=fit_ls, fit_lw=fit_lw, fit_label=fit_label
                    )
        return popt, error
    elif fig==None or fig==0:
        return popt, error
    else:
        print("Do you want to plot the fitted curve?")
        fig = str(input("Type Y for yes and N for no:"))

        if fig=='N':
            return popt, error
        elif fig=='Y':
            fitted_plot(func, xdata, ydata, params=params, xlabel=xlabel, ylabel=ylabel,
                        data_color=data_color, data_marker=data_marker, data_ms=data_ms, data_label=data_label,
                        popt=popt, error=error, l=l, u=u,
                        fit_color=fit_color, fit_ls=fit_ls, fit_lw=fit_lw, fit_label=fit_label
                        )
            return popt, error
        else:
            print("Please type 'Y' or define fig='yes' for Figure.")
            return popt, error





