#===============================================================================
'''
Project: Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Example of vector calculus using symbolic expression
        Example of displaying a vector field, streamlines with in-built
        function, isocontours for certain values
        Example of pathline calculation with several time integration
        methods

Author: mate.pentek@tum.de, kodakkal.anoop@tum.de

Description: The first part deals with defining, setting up and visualizing
        a vector field. Here it represents a steady flow. Isocontours are 
        used to show the velocity magnitude.

        Numerical integration methods are implemented to show the
        differences for pathline calculation. This follows the Lagrangian
        view, i.e. particle tracking. While the flow domain is defined by a
        vector field independent of time, i.e. steady, the particles in the
        flow have a changing velocity as this depends by definition on the
        position. It is suggested to try different settings for deltaTime
        to see how it affects the outcome. 

Proposed exercises: See Jupyter Notebook file
        
Note:   ...
        Module dependencies: 
            python
            numpy
            sympy
            matplotlib.pyplot

Created on:  08.11.2015
Last update: 11.10.2018
'''
#===============================================================================


#===================== BLOCK 1 =======================
# import
import time
import matplotlib.pyplot as plt
import numpy as np
import sympy 
from matplotlib import pyplot as plt

#===================== BLOCK 2 =======================
# declare symbolic variables
x = sympy.Symbol('x')
y = sympy.Symbol('y')
z = sympy.Symbol('z')
t = sympy.Symbol('t')

#===================== BLOCK 3 =======================
##
## modify here for a different velocity field
##
symbolic_velocity_x = -2 #* y
symbolic_velocity_y = 2 #* x
symbolic_velocity_z = 0
print('u = [', symbolic_velocity_x,' ', symbolic_velocity_y,' ', symbolic_velocity_z,']')

#===================== BLOCK 4 =======================
# symbolic_velocity_x = 2
# symbolic_velocity_y = 2
# symbolic_velocity_z = 0
# print('u = [', symbolic_velocity_x,' ', symbolic_velocity_y,' ', symbolic_velocity_z,']')

#===================== BLOCK 5 =======================
symbolic_div = sympy.diff(symbolic_velocity_x, x) + sympy.diff(symbolic_velocity_y, y) + sympy.diff(symbolic_velocity_z, z);
print('div(u)= ', symbolic_div)

#===================== BLOCK 6 =======================
symbolic_curl_x = sympy.diff(symbolic_velocity_z, y) - sympy.diff(symbolic_velocity_y, z)
symbolic_curl_y = sympy.diff(symbolic_velocity_x, z) - sympy.diff(symbolic_velocity_z, x)
symbolic_curl_z = sympy.diff(symbolic_velocity_y, x) - sympy.diff(symbolic_velocity_x, y)
print('curl(u) = [', symbolic_curl_x,' ', symbolic_curl_y,' ', symbolic_curl_z,']')

#===================== BLOCK 7 =======================
# give grid size parameters
grid_spacing_start = -2
grid_spacing_end = 2
grid_spacing_size = 0.4

# generate grid size vectors 1D
grid_spacing_x = np.arange(grid_spacing_start, grid_spacing_end + grid_spacing_size, grid_spacing_size) 
grid_spacing_y = np.arange(grid_spacing_start, grid_spacing_end + grid_spacing_size, grid_spacing_size) 

meshgrid_x, meshgrid_y = np.meshgrid(grid_spacing_x, grid_spacing_y)

# define components and generate vector field components for the 2D
# meshgrid field
##
## modify here for a different velocity field
##
# Ux = -2 * y
# Uy = 2 * x
velocity_x_on_meshgrid = -2 * meshgrid_y #meshgrid_x - meshgrid_y
velocity_y_on_meshgrid  = 2 * meshgrid_x #meshgrid_x + meshgrid_y

#===================== BLOCK 8 =======================
plt.figure(num=1, figsize=(8, 8))
plt.title('Vector field for Vx = ' + str(symbolic_velocity_x) + ' and Vy = ' + str(symbolic_velocity_y))
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(grid_spacing_start-1, grid_spacing_end+1)
plt.ylim(grid_spacing_start-1, grid_spacing_end+1)

plt.quiver(meshgrid_x, meshgrid_y,
           velocity_x_on_meshgrid, velocity_y_on_meshgrid,
           angles='xy',scale_units='xy',scale=10)

#plt.axis('equal')
plt.grid()

#===================== BLOCK 9 =======================
plt.figure(num=2, figsize=(8, 8))
plt.streamplot(meshgrid_x, meshgrid_y,
               velocity_x_on_meshgrid, velocity_y_on_meshgrid,
               density= 0.5)

plt.title('Streamlines for vector field')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim([grid_spacing_start-1, grid_spacing_end+1])
plt.ylim([grid_spacing_start-1, grid_spacing_end+1])
plt.grid(True)

#===================== BLOCK 10 ======================
meshgrid_velocity_magnitude = np.sqrt((np.square(velocity_x_on_meshgrid) + np.square(velocity_y_on_meshgrid)))

# plot separately a contour line for the magnitude of velocity for an initial point
# x0 = 1, y0 = 1 -> looking for magnitude Ux(x0=1)= -2 and Uy(y0=1) = 2
# magnitude is 2*sqrt(2) = 2.8284
value_looked_for = np.sqrt((-2*1)**2 + (2*1)**2)

# plot first plot line of interest 
# after plot contours lines dashed 
plt.figure(num=3, figsize=(8, 8))
# plot a selected countour line
CS = plt.contour(meshgrid_x, meshgrid_y, meshgrid_velocity_magnitude, levels = [value_looked_for], linewidths= 2.5)
plt.clabel(CS, inline=1, fontsize=12)
# plot multiple countour lines
plt.contour(meshgrid_x, meshgrid_y, meshgrid_velocity_magnitude, linestyles='dashed',linewidths=1.5)

plt.title('Vector magnitude contour plot')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim([grid_spacing_start-1, grid_spacing_end+1])
plt.ylim([grid_spacing_start-1, grid_spacing_end+1])
plt.grid(True)

#===================== BLOCK 11 ======================
# start time
start_time = 0.0
# end time
end_time = np.pi
#delta time
delta_time = 0.1   
# time step
n_steps = round((end_time - start_time)/delta_time)
# time series
time_series = np.linspace(start_time, end_time, n_steps)

#===================== BLOCK 12 ======================
# defining the velocity as a function
# one can define with the lambda placeholder for variables
# velocities defined generically as functions depending on x, y, t 
# so positions x and y and time t
##
## modify here for a different velocity field
##
velocity_x = lambda t, x, y: - 2 * y
velocity_y = lambda t, x, y: 2 * x

# initialize with zero values
coord_x = np.zeros(len(time_series))
coord_y = np.zeros(len(time_series))
# initial position
coord_x[0] = 1   
coord_y[0] = 1

# length set to zero
pathline_length = 0

# calculation loop for forward Euler
# x(i) = x(i-1) + Vx(i-1) * delta_t
# y(i) = y(i-1) + Vy(i-1) * delta_t
for i in range(0, len(time_series)-1):
    # calculate new position
    coord_x[i+1] = coord_x[i] + velocity_x(time_series[i] + delta_time, 
                                           coord_x[i], 
                                           coord_y[i]) * delta_time
    coord_y[i+1] = coord_y[i] + velocity_y(time_series[i] + delta_time, 
                                           coord_x[i], 
                                           coord_y[i]) * delta_time 
    
    # increment pathline length
    pathline_length += np.sqrt((coord_x[i+1] - coord_x[i])**2 + (coord_y[i+1] - coord_y[i])**2) 

# store results under new name to be more specific
euler_coord_x = coord_x
euler_coord_y = coord_y
euler_pathline_length = pathline_length 

#===================== BLOCK 13 ======================
# initialize with zero values
coord_x = np.zeros(len(time_series))
coord_y = np.zeros(len(time_series))
# initial position
coord_x[0] = 1   
coord_y[0] = 1

# length set to zero
pathline_length = 0

# calculation loop for RK4
for i in range(0, len(time_series)-1):
    # calculate new position
    k_1 = velocity_x(time_series[i],
                    coord_x[i],
                    coord_y[i])
    k_2 = velocity_x(time_series[i] + 0.5 * delta_time,
                    coord_x[i] + 0.5 * delta_time * k_1,
                    coord_y[i] + 0.5 * delta_time * k_1)
    k_3 = velocity_x(time_series[i] + 0.5 * delta_time,
                    coord_x[i] + 0.5 * delta_time * k_2,
                    coord_y[i] + 0.5 * delta_time * k_2)
    k_4 = velocity_x(time_series[i] + delta_time,
                    coord_x[i] + k_3 * delta_time,
                    coord_y[i] + k_3 * delta_time)
     
    coord_x[i+1] = coord_x[i] + (1.0/6.0) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * delta_time  
     
    k_1 = velocity_y(time_series[i], 
                     coord_x[i], 
                     coord_y[i])
    k_2 = velocity_y(time_series[i] + 0.5 * delta_time, 
                     coord_x[i] + 0.5 * delta_time * k_1, 
                     coord_y[i] + 0.5 * delta_time * k_1)
    k_3 = velocity_y(time_series[i] + 0.5 * delta_time,
                     coord_x[i] + 0.5 * delta_time * k_2, 
                     coord_y[i] + 0.5 * delta_time * k_2)
    k_4 = velocity_y(time_series[i] + delta_time,
                     coord_x[i] + k_3 * delta_time,
                     coord_y[i] + k_3 * delta_time)
    
    coord_y[i+1] = coord_y[i] + (1.0/6.0) * (k_1 + 2 * k_2 + 2 * k_3 + k_4) * delta_time  
    
    # increment pathline length
    pathline_length += np.sqrt((coord_x[i+1] - coord_x[i])**2 + (coord_y[i+1] - coord_y[i])**2) 

# store results under new name to be more specific
rk4_coord_x = coord_x
rk4_coord_y = coord_y
rk4_pathline_length = pathline_length

#===================== BLOCK 14 ======================
# initialize with zero values
coord_x = np.zeros(len(time_series))
coord_y = np.zeros(len(time_series))
# initial position
coord_x[0] = 1   
coord_y[0] = 1

# length set to zero
pathline_length = 0

# predict first step with Euler as Adams-Bashforth needs more information to
# start because the method depends on 2 previous values (positions) of x
# and y
# note that it introduces an inital error, one could substitute it with RK4
# to predict instead and reduce the error
coord_x[1] = coord_x[0] + velocity_x(time_series[0] + delta_time,
                                    coord_x[0],
                                    coord_y[0]) * delta_time
coord_y[1] = coord_y[0] + velocity_y(time_series[0] + delta_time,
                                    coord_x[0],
                                    coord_y[0]) * delta_time

pathline_length += (np.sqrt((coord_x[1] - coord_x[0])**2 + (coord_y[1] - coord_y[0])**2)) 

# calculation loop for Adams-Bashforth
for i in range(0 + 1, len(time_series)-1):   
    # calculate new position
    coord_x[i+1] = coord_x[i] + 3.0/2.0 * velocity_x(time_series[i] + delta_time,
                                                    coord_x[i],
                                                    coord_y[i]) * delta_time \
                              - 1.0/2.0 * velocity_x(time_series[i-1] + delta_time,
                                                    coord_x[i-1],
                                                    coord_y[i-1]) * delta_time
    coord_y[i+1] = coord_y[i] + 3.0/2.0 * velocity_y(time_series[i] + delta_time,
                                                    coord_x[i],
                                                    coord_y[i]) * delta_time \
                              - 1.0/2.0 * velocity_y(time_series[i-1] + delta_time,
                                                    coord_x[i-1],
                                                    coord_y[i-1]) * delta_time

    # increment pathline length
    pathline_length += np.sqrt((coord_x[i+1] - coord_x[i])**2 + (coord_y[i+1] - coord_y[i])**2)

# store results under new name to be more specific
ab_coord_x = coord_x
ab_coord_y = coord_y
ab_pathline_length = pathline_length

#===================== BLOCK 15 ======================
# radius
r = np.sqrt((1)**2 + (1)**2)
# parametric time steps
# taking a center angle theta
theta = np.linspace(0, 2*np.pi, 1000)
# center of circle a,b
a = 0
b = 0
# point on the circle
coord_x = r * np.cos(theta) + a
coord_y = r * np.sin(theta) + b

# store results under new name to be more specific
exact_coord_x = coord_x
exact_coord_y = coord_y
exact_pathline_length = 2 * np.pi * r

#===================== BLOCK 16 ======================
plt.figure(num=4, figsize=(8, 8))
plt.title('Pathline for point x(t(0)) = y(t(0)) = 1, delta_time = ' + str(round(delta_time, 3)))
# Euler
plt.plot(euler_coord_x, euler_coord_y,'-r', label='Euler')
# RK4
plt.plot(rk4_coord_x, rk4_coord_y,'-.b', label='RK4')
# AB
plt.plot(ab_coord_x, ab_coord_y,'--k', label='AB')
# exact analytic
plt.plot(exact_coord_x, exact_coord_y,'-b', label='Exact')

plt.xlabel('x')
plt.ylabel('y')
plt.xlim([grid_spacing_start-1, grid_spacing_end+1])
plt.ylim([grid_spacing_start-1, grid_spacing_end+1])

plt.legend(loc='upper right', shadow=True)
plt.grid(True)

#===================== BLOCK 17 ======================
plt.figure(num=5, figsize=(15, 6))

plt.subplot(1,2,1)
plt.title('Pathline for point x(t(0)) = y(t(0)) = 1, delta_time = ' + str(round(delta_time, 3)))
# Euler
plt.plot(euler_coord_x, euler_coord_y,'-r', label='Euler')
# RK4
plt.plot(rk4_coord_x, rk4_coord_y,'-.b', label='RK4')
# AB
plt.plot(ab_coord_x, ab_coord_y,'--k', label='AB')
# exact analytic
plt.plot(exact_coord_x, exact_coord_y,'-b', label='Exact')

plt.xlabel('x')
plt.ylabel('y')
plt.xlim([0, 2])
plt.ylim([0, 2])
plt.legend(loc='upper right', shadow=True)
plt.grid(True)

plt.subplot(1,2,2)
plt.title('Pathline for point x(t(0)) = y(t(0)) = 1, delta_time = ' + str(round(delta_time, 3)))
# Euler
plt.plot(euler_coord_x, euler_coord_y,'-r', label='Euler')
# RK4
plt.plot(rk4_coord_x, rk4_coord_y,'-.b', label='RK4')
# AB
plt.plot(ab_coord_x, ab_coord_y,'--k', label='AB')
# exact analytic
plt.plot(exact_coord_x, exact_coord_y,'-b', label='Exact')

plt.xlabel('x')
plt.ylabel('y')
plt.xlim([0.95, 1.05])
plt.ylim([0.95, 1.05])
plt.legend(loc='upper right', shadow=True)
plt.grid(True)

plt.show()

#===================== BLOCK 18 ======================
print('Pathline length comparison for dt = ', str(round(delta_time,3)))
print('Euler: ', str(round(euler_pathline_length, 3)))
print('RK4: ', str(round(rk4_pathline_length, 3)))
print('AB: ', str(round(ab_pathline_length, 3)))
print('Exact: ', str(round(exact_pathline_length, 3)))
