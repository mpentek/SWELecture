#===============================================================================
'''
Project: Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Computation of wind speed

Author: mate.pentek@tum.de, kodakkal.anoop@tum.de 
         
Description: Wind profile generation based upon EC wind provisions
        and example 1.1

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
gust_windspeed = 40.12
# 1.4 is the code specified factor to convert from gust to mean wind speed 
mean_windspeed = gust_windspeed/1.4 

#===================== BLOCK 3 =======================
height_total = 200
storey_height = 4 # height of each floor 
air_density = 1.2 # airdensity in kg/m3
height = np.arange(0, height_total, storey_height)

#===================== BLOCK 4 =======================
a_mean  = 1.18
alpha_mean = 0.18
a_gust = 1.61
alpha_gust = 0.095

umean_1 = a_mean * mean_windspeed * (height/10)**alpha_mean
ugust_1 = a_gust * mean_windspeed * (height/10)**alpha_gust

#===================== BLOCK 5 =======================
plt.figure(num=1, figsize=(8, 6))
plt.plot(umean_1, height, label = 'Mean velocity profile')
plt.plot(ugust_1, height, label = 'Gust velocity profile')
plt.ylabel('Height')
plt.xlabel('wind speed m/s')
plt.title('Wind Speed along height - Terrain category I')
plt.legend()
plt.grid(True)

#===================== BLOCK 6 =======================
a_mean  = 0.56
alpha_mean = 0.3
a_gust = 1.05
alpha_gust = 0.2

umean_1 = a_mean * mean_windspeed * (height/10)**alpha_mean
ugust_1 = a_gust * mean_windspeed * (height/10)**alpha_gust

#===================== BLOCK 7 =======================
plt.figure(num=2, figsize=(8, 6))
plt.plot(umean_1, height, label = 'Mean velocity profile')
plt.plot(ugust_1, height, label = 'Gust velocity profile')
plt.ylabel('Height')
plt.xlabel('Wind speed m/s')
plt.title('Wind Speed along height - Terrain category IV')
plt.legend()
plt.grid(True)

plt.show()