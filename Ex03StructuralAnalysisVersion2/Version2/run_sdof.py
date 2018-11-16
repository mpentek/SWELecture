#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Analysis type base class and derived classes specific types

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de, catharina.czech@tum.de, peter.kupas@tum.de

      
Note:   UPDATE: The script has been written using publicly available information and 
        data, use accordingly. It has been written and tested with Python 2.7.9.
        Tested and works also with Python 3.6 (already see differences in print).
        Module dependencies (-> line 61-74): 
            python
            numpy
            sympy
            matplotlib.pyplot

Created on:  22.11.2017
Last update: 16.11.2018
'''
#===============================================================================
import numpy as np
#import matplotlib.pyplot as plt


from source.analysis_type import*
from source.load_type import*
from source.custom_files import *
from source.mdof_solver import *


# =========================sdof system properties==========================  
parameter_file = open('cosim_sdof_parameters.json','r')
Parameters = json.loads(parameter_file.read())
solver_settings = Parameters['solver_settings']
sdof_solver = MDoFSolver(solver_settings,1)

#========================= static analysis==========================  
force_static = np.array([1.5]) # external static force acting on the system

static_analysis = StaticAnalysis(sdof_solver.model)
static_analysis.solve(force_static)
static_analysis.plot_solve_result()

#========================= eigenvalue analysis ==========================  
eigenvalue_analysis = EigenvalueAnalysis(sdof_solver.model)
eigenvalue_analysis.solve()
eigenvalue_analysis.plot_selected_eigenmode(1)
eigenvalue_analysis.animate_selected_eigenmode(1)

#========================= dynamic analysis ==========================  

# time parameters 
time_parameters = Parameters["problem_data"]
start_time = time_parameters["start_time"]
end_time = time_parameters["end_time"]
dt = time_parameters["time_step"]
array_time = np.arange (start_time,end_time + dt, dt)

# dynamic forces
"""
Choose from "signalSin", "signalRand", "signalConst", "signalSuperposed" or 
for free vibration choose "signalNone" 
"""
# external dynamic force acting on the system
freq = 10
force_dynamic = load_type("signalSin", array_time, 1, freq, force_static) 

dynamic_analysis = DynamicAnalysis(sdof_solver, force_dynamic, time_parameters) 
dynamic_analysis.solve()
dynamic_analysis.plot_selected_time_step(0.75)
dynamic_analysis.animate_time_history()