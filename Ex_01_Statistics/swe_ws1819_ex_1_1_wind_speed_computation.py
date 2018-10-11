#===============================================================================
'''
Project: Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Computation of wind speed

Author: mate.pentek@tum.de, kodakkal.anoop@tum.de 
         
Description: Gust wind speed estimation example 
        based upon 2.3 from Holmes, J.D., Wind loading 
        of structures“, 2nd ed.

Proposed exercises: See Jupyter Notebook file

Note:   ...
        Module dependencies:
            python
            matplotlib
            numpy
            custom_utilities.py

Created on:  30.11.2015
Last update: 10.10.2018
'''
#===============================================================================


#===================== BLOCK 1 =======================
# import
import matplotlib.pyplot as plt
import numpy as np

#===================== BLOCK 2 =======================
winddata = np.loadtxt('east_sale.dat')
year = winddata[:,0]
maxgust = winddata[:,1]

#===================== BLOCK 3 =======================
gust_m = np.mean(maxgust)
gust_std = np.std(maxgust)

#===================== BLOCK 4 =======================
plt.figure(num=1, figsize=(15, 4))
plt.plot(year, maxgust)
plt.scatter(year, maxgust)
plt.hlines(gust_m, min(year), max(year), color='r', label = 'mean')
plt.hlines(gust_m + gust_std, min(year), max(year), color='g', label = 'mean+std')
plt.hlines(gust_m - gust_std, min(year), max(year), color='g', label = 'mean-std')
plt.ylabel('Maximum annual gust wind speed m/s')
plt.xlabel('year')
plt.title('Gust Wind Speed - East Sale, Australia')
plt.legend()
plt.grid(True)

#===================== BLOCK 5 =======================
plt.figure(num=2, figsize=(8, 6))
plt.hist(maxgust, bins= np.arange(min(maxgust),max(maxgust)+1), rwidth=0.8)
plt.ylabel('Frequency')
plt.xlabel('Maximum annual gust wind speed m/s')
plt.title('Gust Wind Speed - East Sale, Australia histogram')
plt.grid(True)

#===================== BLOCK 6 =======================
gust_sorted = np.sort(maxgust)
max_rank = len(gust_sorted)
rank = np.arange(1, max_rank + 1)

#===================== BLOCK 7 =======================
gumbel_prob_nonexc = rank / (max_rank + 1)
gumbel_red_var = -np.log(-np.log(gumbel_prob_nonexc))

#===================== BLOCK 8 =======================
[gumbel_slope, gumbel_mode] = np.polyfit(gumbel_red_var, gust_sorted, 1)

#===================== BLOCK 9 =======================
plt.figure(num=3, figsize=(8, 6))
plt.scatter(gumbel_red_var, gust_sorted)
x = np.linspace(min(gumbel_red_var), max(gumbel_red_var),50)
plt.plot(x, gumbel_mode + gumbel_slope * x)
plt.text(0.1, 0.8,'Mode = ' + str(round(gumbel_mode,2)) +
         '\nSlope = ' + str(round(gumbel_mode,2)),
         transform=plt.gca().transAxes)
plt.ylabel('Maximum annual gust wind speed m/s')
plt.xlabel('Reduced variate')
plt.title('Gumbel Method')
plt.grid(True)

#===================== BLOCK 10 ======================
return_period = np.arange(10,1000,10)
gumbel_predicted_gustwind = gumbel_mode + gumbel_slope * (-np.log(-np.log(1-1/return_period)))

#===================== BLOCK 11 ======================
plt.figure(num=4, figsize=(15, 6))
plt.subplot(1,2,1)
plt.plot(return_period, gumbel_predicted_gustwind)
plt.ylabel('Predicted gust wind speed m/s')
plt.xlabel('Return period (Years)')
plt.title('Gumbel Method')
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(return_period, gumbel_predicted_gustwind)
plt.xscale('log')
plt.ylabel('Predicted gust wind speed m/s')
plt.xlabel('Return period (Years)')
plt.title('Gumbel Method')
plt.grid(True)

#===================== BLOCK 12 ======================
gringorten_prob_nonexc = (rank-0.44) /(max_rank+0.12)
gringorten_red_var = -np.log(-np.log(gringorten_prob_nonexc))
[gringorten_slope, gringorten_mode] = np.polyfit(gringorten_red_var, gust_sorted,1)
gringorten_predicted_gustwind = gringorten_mode + gringorten_slope * (-np.log(-np.log(1-1/return_period)))

#===================== BLOCK 13 ======================
plt.figure(num=5, figsize=(8, 6))
plt.scatter(gringorten_red_var, gust_sorted)
x = np.linspace(min(gringorten_red_var), max(gringorten_red_var),50)
plt.plot(x, gringorten_mode + gringorten_slope * x)
plt.text(0.1, 0.8,'Mode = ' + str(round(gringorten_mode,2)) +
         '\nSlope = ' + str(round(gringorten_slope,2)),
         transform=plt.gca().transAxes)
plt.ylabel('Maximum annual gust wind speed m/s')
plt.xlabel('Reduced variate')
plt.title('Gringorten Method')
plt.grid(True)

#===================== BLOCK 14 ======================
moments_slope = np.sqrt(6)/np.pi * gust_std
moments_mode  = gust_m - 0.5772 * moments_slope
moments_predicted_gustwind = moments_mode + moments_slope * (-np.log(-np.log(1-1/return_period)))

#===================== BLOCK 15 ======================
plt.figure(num=6, figsize=(15, 6))
plt.subplot(1,2,1)
plt.plot(return_period, gumbel_predicted_gustwind, label ='Gumbel Method')
plt.plot(return_period, gringorten_predicted_gustwind, label ='Gringgorten Method')
plt.plot(return_period, moments_predicted_gustwind, label ='Method of moments')
plt.ylabel('Predicted gust wind speed m/s')
plt.xlabel('Return period (Years)')
plt.title('Gust wind speed prediction')
plt.legend()
plt.grid(True)

plt.subplot(1,2,2)
plt.plot(return_period, gumbel_predicted_gustwind, label ='Gumbel Method')
plt.plot(return_period, gringorten_predicted_gustwind, label ='Gringgorten Method')
plt.plot(return_period, moments_predicted_gustwind, label ='Method of moments')
plt.xscale('log')
plt.ylabel('Predicted gust wind speed m/s')
plt.xlabel('Return period (Years)')
plt.title('Gust wind speed prediction')
plt.legend()
plt.grid(True)

plt.show()

#===================== BLOCK 16 ======================
# winddata = np.loadtxt('jeddah_airport.dat')
# year = winddata[:,0]
# maxgust = winddata[:,2]