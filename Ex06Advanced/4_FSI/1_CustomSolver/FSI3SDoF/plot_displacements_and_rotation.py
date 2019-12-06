# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS17-18
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek
        
Author: philipp.bucher@tum.de, mate.pentek@tum.de 

Description: Script for plotting data over time

Created on:  05.12.2015
Last update: 06.12.2019
'''
# ===============================================================================

from pylab import *  # import necessary module for plotting
import json

with open("ProjectParametersCSD.json") as parameter_file:
    project_parameters = json.load(parameter_file)

file_disp_x_vec_res = "sdof_" + project_parameters["problem_data"]["problem_name"].lower(
) + "_" + project_parameters["model_data"]["dof_type"][0].lower() + ".dat"
file_disp_y_vec_res = "sdof_" + project_parameters["problem_data"]["problem_name"].lower(
) + "_" + project_parameters["model_data"]["dof_type"][1].lower() + ".dat"
file_rot_vec_res = "sdof_" + project_parameters["problem_data"]["problem_name"].lower(
) + "_" + project_parameters["model_data"]["dof_type"][2].lower() + ".dat"

# file_disp_x_sep_res    = "displacementX_sepResidual.dat"
# file_disp_y_sep_res    = "displacementY_sepResidual.dat"
# file_rot_sep_res      = "rotation_sepResidual.dat"
# ======================================================


# *** read in the aerodynamic force/moment results ****************
time_dx_vr = loadtxt(file_disp_x_vec_res, skiprows=1, usecols=(0,))
time_dy_vr = loadtxt(file_disp_y_vec_res, skiprows=1, usecols=(0,))
time_r_vr = loadtxt(file_rot_vec_res, skiprows=1, usecols=(0,))
disp_x_vr = loadtxt(file_disp_x_vec_res, skiprows=1, usecols=(1,))
disp_y_vr = loadtxt(file_disp_y_vec_res, skiprows=1, usecols=(1,))
rot_vr = loadtxt(file_rot_vec_res, skiprows=1, usecols=(1,))

# time_dx_sr  = loadtxt(file_disp_x_sep_res, skiprows=1, usecols = (0,))
# time_dy_sr  = loadtxt(file_disp_y_sep_res, skiprows=1, usecols = (0,))
# time_r_sr   = loadtxt(file_rot_sep_res, skiprows=1, usecols = (0,))
# disp_x_sr   = loadtxt(file_disp_x_sep_res, skiprows=1, usecols = (1,))
# disp_y_sr   = loadtxt(file_disp_y_sep_res, skiprows=1, usecols = (1,))
# rot_sr     = loadtxt(file_rot_sep_res, skiprows=1, usecols = (1,))

lw = 1.5  # line width in plot

plt.subplot(3, 1, 1)  # plot displacement X
plt.plot(time_dx_vr, disp_x_vr, 'b-',
         linewidth=lw, label='Displacement vecRes')
# plt.plot(time_dx_sr, disp_x_sr, 'r-', linewidth=lw, label='Displacement sepRes')
plt.title('DISPLACEMENT X')
plt.grid(True)
# plt.xlim(xmax=max(time_dx_vr[-1], time_dx_sr[-1])) # apply the maximum computation time as upper x-axis limit
# change y-axis numbers to scientific format
plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
plt.ylabel('Displacement [m]')
plt.legend()

plt.subplot(3, 1, 2)  # plot displacement Y
plt.plot(time_dy_vr, disp_y_vr, 'b-',
         linewidth=lw, label='Displacement vecRes')
# plt.plot(time_dy_sr, disp_y_sr, 'r-', linewidth=lw, label='Displacement sepRes')
plt.title('DISPLACEMENT Y')
plt.grid(True)
# plt.xlim(xmax=max(time_dy_vr[-1], time_dy_sr[-1])) # apply the maximum computation time as upper x-axis limit
# change y-axis numbers to scientific format
plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
plt.ylabel('Displacement [m]')
plt.legend(loc='upper left')

plt.subplot(3, 1, 3)  # plot rotation
plt.plot(time_r_vr, rot_vr, 'b-', linewidth=lw, label='Rotation vecRes')
# plt.plot(time_r_sr, rot_sr, 'r-', linewidth=lw, label='Rotation sepRes')
plt.title('Rotation')
plt.grid(True)
# plt.xlim(xmax=max(time_r_vr[-1], time_r_sr[-1])) # apply the maximum computation time as upper x-axis limit
# change y-axis numbers to scientific format
plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
plt.ylabel('Rotation [rad]')
plt.xlabel('Time [s]')
plt.legend()

subplots_adjust(hspace=.35)  # increase the space between the subplots
plt.show()
