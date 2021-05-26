# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: philipp.bucher@tum.de, mate.pentek@tum.de 

Description: Script for plotting aerodynamic forces and moments

Created on:  05.12.2015
Last update: 06.12.2019
'''
# ===============================================================================

import json

import sys
from matplotlib.pylab import *
import numpy as np

with open("ProjectParametersCFD.json") as parameter_file:
    project_parameters = json.load(parameter_file)

file2read = project_parameters["auxiliar_process_list"][0]["Parameters"]["model_part_name"] + "_drag.dat"

# ======================================================

# *** read in the aerodynamic force/moment results ****************
simul_time = loadtxt(file2read, skiprows=3, usecols=(0,))
force_x = loadtxt(file2read, skiprows=3, usecols=(1,))
force_y = loadtxt(file2read, skiprows=3, usecols=(2,))
force_z = loadtxt(file2read, skiprows=3, usecols=(3,))

try:  # check if results for moments exist and read them if they exist
    moment_x = loadtxt(file2read, skiprows=3, usecols=(4,))
    moment_y = loadtxt(file2read, skiprows=3, usecols=(5,))
    moment_z = loadtxt(file2read, skiprows=3, usecols=(6,))
    moments_read = True
    number_of_subplots = 2
except:
    moments_read = False
    number_of_subplots = 1

# *** plot the aerodynamic forces/moments *************************
lw = 3  # line width in plot

plt.subplot(number_of_subplots, 1, 1)  # plot the forces
plt.plot(simul_time, force_x, 'b-', linewidth=lw,
         label='Aerodyn. Force X without Ramp')
plt.plot(simul_time, force_y, 'r-', linewidth=lw,
         label='Aerodyn. Force Y without Ramp')
plt.plot(simul_time, force_z, 'g-', linewidth=lw,
         label='Aerodyn. Force Z without Ramp')
plt.title('Aerodynamic Forces')
plt.ylabel('Force [N]')
# change y-axis numbers to scientific format
plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
plt.grid(True)
# apply the maximum computation time as upper x-axis limit
plt.xlim(xmax=simul_time[-1])
plt.legend()
if not moments_read:  # plot time label if second plot doesn't exist
    plt.xlabel('Time [s]')

if moments_read:  # plot the moments if they have been read
    plt.subplot(number_of_subplots, 1, 2)
    plt.plot(simul_time, moment_x, 'b-',
             linewidth=lw, label='Aerodyn. Moment X')
    plt.plot(simul_time, moment_y, 'r-',
             linewidth=lw, label='Aerodyn. Moment Y')
    plt.plot(simul_time, moment_z, 'k-',
             linewidth=lw, label='Aerodyn. Moment Z')
    plt.title('Aerodynamic Moments')
    plt.xlabel('Time [s]')
    plt.ylabel('Moment [Nm]')
    # change y-axis numbers to scientific format
    plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
    plt.grid(True)
    # apply the maximum computation time as upper x-axis limit
    plt.xlim(xmax=simul_time[-1])
    plt.legend()

subplots_adjust(hspace=.35)  # increase the space between the subplots
plt.show()
