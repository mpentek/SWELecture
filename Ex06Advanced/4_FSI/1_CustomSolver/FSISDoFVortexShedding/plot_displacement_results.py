#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: philipp.bucher@tum.de, mate.pentek@tum.de 

Description: Script for plotting data over time

Created on:  05.12.2015
Last update: 13.11.2087
'''
#===============================================================================


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from pylab import *

# ======================================================
import json
    
with open("ProjectParametersCSD.json") as parameter_file:
    project_parameters = json.load(parameter_file)

file2read = "sdof_" + project_parameters["problem_data"]["problem_name"].lower() + "_"+ project_parameters["model_data"]["dof_type"].lower() + ".dat"  
# ======================================================

# read the file
simul_time = loadtxt(file2read, skiprows=1, usecols = (0,))
time_data = loadtxt(file2read, skiprows=1, usecols = (1,))

# set up the plot
fig = plt.figure()
ax = plt.axes(xlim=(min(simul_time), max(simul_time)), ylim=(min(time_data), max(time_data)))
line, = ax.plot([], [], lw=2)
ax.set_xlabel("Time [s]")
ax.set_ylabel("Displacement [m]")
ax.set_title('Displacement Results for SDoF')
ax.plot(simul_time, time_data, "-b", lw=1.0)
plt.ticklabel_format(axis='y',style='sci',scilimits=(3,3)) # change y-axis numbers to scientific format
plt.grid(True)

plt.show()
