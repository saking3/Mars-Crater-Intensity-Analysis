# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 06:12:36 2020
The purpose of this script is to take a tif of a crater, make a curve of the 
frequency of the intensity, and break that curve down into its component curves.
Each component curve is then integrated to get the area. 

To use the code, start by initially running the first portion to get the curve.
Then, select and run the portion of code that will create the amount of curves 
you would like the frequency curve to be broken into. 

In progress: finding a way to put the number of curves done into a loop. Calculating Jacobian to plug into minimizer fxn.
"""

#import all the proper packages 
import matplotlib.pyplot as plt
import numpy as np
import os 
from osgeo import gdal 
import scipy 
import seaborn
from numdifftools import Jacobian

# set working directory 
os.chdir("C:/Users/Sarah/OneDrive/Reccuring Slope Lineae on Mars/TIFS for Slope Analysis")

# read the tiff and make it on object 
dunetest = gdal.Open("cratertest.tif") #reads tiff 

# convert the tiff to an array and clean out any n/a values
ar = np.array(dunetest.GetRasterBand(1).ReadAsArray()) 
ar = ar[~np.isnan(ar)]

#calculate the length of the list
size = ar.size 

# make a density curve and plot it 
pdfplot = seaborn.distplot(ar)
plt.show()

#get the data points from the density curve
pdfx, pdfy = seaborn.distplot(ar).get_lines()[0].get_data()
de = list(range(102, 128)) 
intensity = np.delete (pdfx, de)
frequency = np.delete (pdfy, de)
frequency *= size #plot the data as a frequency and not a probability

################################ 2 CURVES! ################################
#define Gaussian curves 
def f1(X):
    return np.abs(X[0]) * np.exp(-(intensity - np.abs(X[2]))**2 / (2 * np.abs(X[4]**2)))
def f2(X): 
    return np.abs(X[1]) * np.exp(-(intensity - np.abs(X[3]))**2 / (2 * np.abs(X[5]**2)))  

#define total curve function
def f2_total(X):
    return f1(X) + f2(X)

#define minimum function
def minimum(X): 
    return frequency - f2_total(X) 

#define initial parameters
i = 0 #for the while loop
mins = 10000000000000000000000 #arbitrary large number for while loop to minimize
bestX = [] #will save the best variables to plug into the Gaussian fxns 
while i < 500: #setting i to a smaller number will make this faster
    X = []
    #set the Gaussian values
    X[0:1] = [np.median(frequency)] * 2 #1/(standard deviation * sqrt(2pi))
    X[2:3] = np.random.randint(int(np.amin(intensity)), int(np.amax(intensity)), 2) #peak values
    X[4:5] = [10] * 2 #standard deviation
    #minimize the minimum function
    hi = scipy.optimize.least_squares(minimum, X, method = "trf")
    #optimized parameters from the minimization
    X = hi.x
    #if residual sum of squares is better than the min, the X values get saved
    rss = hi.cost
    if rss < mins:
        mins = rss
        bestX = X
    i += 1 
    print(mins)
    print(bestX)
   
#plotting the results onto a graph 
plt.plot(intensity, frequency, color = "red")
plt.plot(intensity,f2_total(bestX), color="blue")
plt.plot(intensity,f1(bestX))
plt.plot(intensity,f2(bestX))
plt.show()

#integrate each curve
area1 = scipy.integrate.quad(f1(bestX), 0, 250)
area2 = scipy.integrate.quad(f2(bestX), 0, 250)

area1 = area1 * 25 

################################ 3 CURVES! ################################
#define Gaussian curves
def f1(X):
    return np.abs(X[0]) * np.exp(-(intensity - np.abs(X[3]))**2 / (2 * np.abs(X[6])**2))
def f2(X): 
    return np.abs(X[1]) * np.exp(-(intensity - np.abs(X[4]))**2 / (2 * np.abs(X[7]**2)))  
def f3(X): 
    return np.abs(X[2]) * np.exp(-(intensity - np.abs(X[5]))**2 / (2 * np.abs(X[8]**2)))

#define total curve function
def f3_total(X):
    return f1(X) + f2(X) + f3(X) 

#define minimum function
def minimum(X): 
    return frequency - f3_total(X)

#define initial parameters
i = 0 #for the while loop
mins = 10000000000000000000000 #arbitrary large number for while loop to minimize
bestX = [] #will save the best variables to plug into the Gaussian fxns 
while i < 1000: #setting i to a smaller number will make this faster
    X = []
    X[0:2] = [max(frequency)] * 3 #1/(standard deviation * sqrt(2pi))
    X[3:5] = np.random.randint(0, int(np.amax(intensity)), 3) #mean
    X[6:8] = [10] * 3 #standard deviation    
    #minimize the minimum function
    hi = scipy.optimize.least_squares(minimum, X, method = "lm") #least squares
    #optimized parameters from the minimization
    X = hi.x
    #if residual sum of squares is better than the min, the X values get saved
    rss = hi.cost
    if rss < mins:
        mins = rss
        bestX = X
    i += 1 
    print(mins)
    print(bestX)
   
#plotting the results onto a graph 
plt.plot(intensity, frequency, color = "red")
plt.plot(intensity,f3_total(bestX), color="blue")
plt.plot(intensity,f1(bestX))
plt.plot(intensity,f2(bestX))
plt.plot(intensity,f3(bestX))
plt.show()

#integrate each curve
area1 = scipy.integrate.quad(f1(bestX), 0, 250)
area2 = scipy.integrate.quad(f2(bestX), 0, 250)
area3 = scipy.integrate.quad(f3(bestX), 0, 250)
################################ 4 CURVES! ################################
#define Gaussian curves
def f1(X):
    return np.abs(X[0]) * np.exp(-(intensity - np.abs(X[4]))**2 / (2 * np.abs(X[8])**2))
def f2(X): 
    return np.abs(X[1]) * np.exp(-(intensity - np.abs(X[5]))**2 / (2 * np.abs(X[9])**2))  
def f3(X): 
    return np.abs(X[2]) * np.exp(-(intensity - np.abs(X[6]))**2 / (2 * np.abs(X[10])**2))
def f4(X):
    return np.abs(X[3]) * np.exp(-(intensity - np.abs(X[7]))**2 / (2 * np.abs(X[11])**2))

#define total curve function
def f4_total(X):
    return f1(X) + f2(X) + f3(X) + f4(X)

#define minimum function
def minimum(X): 
    return frequency - f4_total(X)

#define initial parameters
i = 0 #for the while loop
mins = 10000000000000000000000 #arbitrary large number for while loop to minimize
bestX = [] #will save the best variables to plug into the Gaussian fxns 
while i < 800:  #setting i to a smaller number will make this faster
    X = []
    X[0:3] = [max(frequency)] * 4 #1/(standard deviation * sqrt(2pi))
    X[4:7] = np.random.randint(0, int(np.amax(intensity)), 4) #mean
    X[8:11] = [10] * 4 #standard deviation    
    #minimize the minimum function
    hi = scipy.optimize.least_squares(minimum, X, jac = minjac, method = "trf") #least squares
    #optimized parameters from the minimization
    X = hi.x
    #if residual sum of squares is better than the min, the X values get saved
    rss = hi.cost
    if rss < mins:
        mins = rss
        bestX = X
    i += 1 
    print(mins)
    print(bestX)
   
#plotting the results onto a graph 
plt.plot(intensity, frequency, color = "red")
plt.plot(intensity,f4_total(bestX), color="blue")
plt.plot(intensity,f1(bestX))
plt.plot(intensity,f2(bestX))
plt.plot(intensity,f3(bestX))
plt.plot(intensity,f4(bestX))
plt.show()

#integrate each curve
area1 = scipy.integrate.quad(f1(bestX), 0, 250)
area2 = scipy.integrate.quad(f2(bestX), 0, 250)
area3 = scipy.integrate.quad(f3(bestX), 0, 250)
area4 = scipy.integrate.quad(f4(bestX), 0, 250)

################################ 5 CURVES! ################################
def f1(X):
    return np.abs(X[0]) * np.exp(-(intensity - np.abs(X[5]))**2 / (2 * np.abs(X[10])**2))
def f2(X): 
    return np.abs(X[1]) * np.exp(-(intensity - np.abs(X[6]))**2 / (2 * np.abs(X[11])**2))  
def f3(X): 
    return np.abs(X[2]) * np.exp(-(intensity - np.abs(X[7]))**2 / (2 * np.abs(X[12])**2))
def f4(X):
    return np.abs(X[3]) * np.exp(-(intensity - np.abs(X[8]))**2 / (2 * np.abs(X[13])**2))
def f5(X):
    return np.abs(X[4]) * np.exp(-(intensity - np.abs(X[9]))**2 / (2 * np.abs(X[14])**2))

#define total curve function
def f5_total(X):
    return f1(X) + f2(X) + f3(X) + f4(X) + f5(X)

#define minimum function
def minimum(X): 
    return frequency - f5_total(X)

#define initial parameters
i = 0 #for the while loop
mins = 10000000000000000000000 #arbitrary large number for while loop to minimize
bestX = [] #will save the best variables to plug into the Gaussian fxns 
while i < 1000:  #setting i to a smaller number will make this faster
    X = []
    X[0:4] = [max(frequency)] * 5 #1/(standard deviation * sqrt(2pi))
    X[5:9] = np.random.randint(0, int(np.amax(intensity)), 5) #mean
    X[10:14] = [10] * 5 #standard deviation    
    #minimize the minimum function
    hi = scipy.optimize.least_squares(minimum, X, method = "lm") #least squares
    #optimized parameters from the minimization
    X = hi.x
    #if residual sum of squares is better than the min, the X values get saved
    rss = hi.cost
    if rss < mins:
        mins = rss
        bestX = X
    i += 1 
    print(mins)
    print(bestX)
   
#plotting the results onto a graph 
plt.plot(intensity, frequency, color = "red")
plt.plot(intensity,f4_total(bestX), color="blue")
plt.plot(intensity,f1(bestX))
plt.plot(intensity,f2(bestX))
plt.plot(intensity,f3(bestX))
plt.plot(intensity,f4(bestX))
plt.plot(intensity, f5(bestX))
plt.show()

#integrate each curve
area1 = scipy.integrate.quad(f1(bestX), 0, 250)
area2 = scipy.integrate.quad(f2(bestX), 0, 250)
area3 = scipy.integrate.quad(f3(bestX), 0, 250)
area4 = scipy.integrate.quad(f4(bestX), 0, 250)
area5 = scipy.integrate.quad(f5(bestX), 0, 250)
