#!/bin/python3
# -*- coding: utf-8 -*-

#Created on Wed Apr 8 06:12:36 2020
#The purpose of this script is to take a tif of a crater, make a curve of the
#frequency of the intensity, and break that curve down into its component curves.
#Each component curve is then integrated to get the area.

#To use the code, start by initially running the first portion to get the curve.
#Then, select and run the portion of code that will create the amount of curves
#you would like the frequency curve to be broken into.

#In progress: finding a way to put the number of curves done into a loop.

#import all the proper packages
import matplotlib.pyplot as plt
import numpy as np
import os, scipy, seaborn
from osgeo import gdal
from numdifftools import Jacobian
from statistics import median

# set working directory
os.chdir("D:/CTX Images")

#os.list('D:/CTX Images') #might work might not

images = ['E-160N-44.tif', 'E084N-08.tif', 'E-068N-16.tif', 'E-072N-12.tif']

for im in images:
    # read the tiff and make it on object
    dunetest = gdal.Open(im)

    # convert the tiff to an array and clean out any n/a values
    ar = np.array(dunetest.GetRasterBand(1).ReadAsArray())
    ar = ar[~np.isnan(ar)]

    #median of the list
    print('MEDIAN: '+str(median(ar))+'\n')

    # make a density curve and plot it
    pdfplot = seaborn.distplot(ar)
    plt.show()

    #get the data points from the density curve
    pdfx, pdfy = seaborn.distplot(ar).get_lines()[0].get_data()
    de = list(range(102, 128))
    intensity = np.delete (pdfx, de)
    frequency = np.delete (pdfy, de)
    frequency *= ar.size

    area = []
    peak = []

    def fs(k, X):
        flist = []
        total = 0
        for i in range(k):
            flist.append(np.abs(X[i]) * np.exp(-(intensity - np.abs(X[k+i]))**2 / (2 * np.abs(X[(k*2)+i])**2)))
            total = flist[i] + total
        minimum = frequency - total
        peak.append(max(flist[k]))
        return minimum, flist, total

    def plot(k, bestX, total):
        plt.plot(intensity, frequency, color = "red")
        plt.plot(intensity,total(bestX), color="blue")
        for i in range(k):
            plt.plot(intensity, fs(k, bestX))
        plt.show()

    functions = [2,3,4,5]
    for k in functions:
        i = 0
        mins = 10000000000000000000000
        bestX = []
        flist = []
        total = 0
        while i < 500:
            X = []
            X[0:(k-1)] = [np.median(frequency)] * k
            X[k:((k*2)-1)] = np.random.randint(int(np.amin(intensity)), int(np.amax(intensity)), k)
            X[(k*2):((k*3)-1)] = [10] * k
            minimum, flist, total = fs(k, bestX)
            hi = scipy.optimize.least_squares(minimum, X, method = "trf")
            X = hi.x
            rss = hi.cost
            if rss < mins:
                mins = rss
                bestX = X
            i += 1
            #print('\n MINS: ' + str(mins) + '\n bestX: ' + str(bestX))
            
        plot(k, flist, total)
        for i in range(k):
            area.append(scipy.integrate.simps(flist[i](bestX), intensity))

    for i in peak:
        print(i)
    for i in area:
        print(i)
