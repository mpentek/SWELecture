#===============================================================================
'''
Project: Structural Wind Engineering WS19-20
        Chair of Structural Analysis @ TUM -R. Wuchner, M. Pentek
        
        MDoF system - structural parameters for a sample highrise

Author: mate.pentek@tum.de

        Based upon a collaborative project with A. Michalski, str.ucture GmbH 
        from Stuttgart
                
Description: This is a script containing a proposed build for the stiffness and 
        mass matrix for an MDoF model of a highrise. Data is derived from a full 
        FEM model in Sofistik.
 
Note:   ...

Created on:  24.11.2016
Last update: 27.09.2019
'''
#===============================================================================
import numpy as np


# function definition for external call
def get_mass():
    return mass

def get_stiffness():
    return stiffness

def get_height_coordinates():
    return z

# number of floors and building height
number_of_floors = 60

level_height = 3.5
z = np.zeros(number_of_floors+1)
for i in range(number_of_floors+1):
    z[i] = level_height * i

# setting up the mass matrix
mass = np.zeros((number_of_floors, number_of_floors))

# precalculated masses - for floors and columns in [kg]
mass_columns = 1500000
mass_floor   = 750000

[rows, columns] = mass.shape

# mass-lumping
for i in range(rows):
    if i == 0: #first element
        mass[i,i] = mass_floor + 1.0 * mass_columns
    elif ((i > 0) & (i <= (rows-2))): #other bottom level
        mass[i,i] = mass_floor + 1.0 * mass_columns 
    elif i == (rows-1): #top level
        mass[i,i] = mass_floor + 0.5 * mass_columns
            
print("Total mass check = " + str(np.sum(mass)) + " [kg]")

# setting up the stifness matrix
stiffness = np.zeros((number_of_floors, number_of_floors))

k_x = 1.5e10 # in [N/m]
# k_x = 5.0e9 # for lower stiffness
        
[rows,columns] = mass.shape 

for i in range(rows):
    if (i == 0): # first (bottom) row
        stiffness[i,i] = k_x * 2
        stiffness[i,i+1] = k_x * (-1) 
    elif ((i > 0) & (i <= (rows-2))): # intermediate rows
        stiffness[i,i-1] = k_x * (-1) 
        stiffness[i,i] = k_x * 2
        stiffness[i,i+1] = k_x * (-1) 
    elif i == (rows-1): # top row
        stiffness[i,i-1] = k_x * (-1) 
        stiffness[i,i] = k_x * 1