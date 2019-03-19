# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Sample MDoF shear run

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de, catharina.czech@tum.de, peter.kupas@tum.de
    
Note: ...

Created on:  22.11.2017
Last update: 16.11.2018
'''
# ===============================================================================
import numpy as np
import json

from source.analysis_type import*
from source.load_type import*
from source.custom_files import *
from source.mdof_solver import *

parameter_file = open('cosim_mdof_cantilever_shear_2d_parameters.json', 'r')
Parameters = json.loads(parameter_file.read())
solver_settings = Parameters['solver_settings']
mdof_solver = MDoFSolver(solver_settings, 1)

# =========================static analysis==========================
# static force definition
z = mdof_solver.model.nodal_coordinates["y0"]
# Generating wind profile
wind_velocity = 28.7  # m/s
velocity_vector = 1.05 * wind_velocity * pow(z/10, 0.2)
wind_force = 0.5 * 1.2 * velocity_vector**2 * 30
force_static = wind_force

static_analysis = StaticAnalysis(mdof_solver.model)
static_analysis.solve(force_static)
static_analysis.plot_solve_result()

# =========================eigenvalue analysis ==========================
eigenvalue_analysis = EigenvalueAnalysis(mdof_solver.model)
eigenvalue_analysis.solve()
eigenvalue_analysis.plot_selected_eigenmode(1)
eigenvalue_analysis.plot_selected_first_n_eigenmodes(3)
eigenvalue_analysis.animate_selected_eigenmode(1)

# =========================dynamic analysis ==========================
# time parameters
time_parameters = Parameters["problem_data"]
start_time = time_parameters["start_time"]
end_time = time_parameters["end_time"]
dt = time_parameters["time_step"]
array_time = np.arange(start_time, end_time + dt, dt)

# dynamic forces

# choose from "signalSin", "signalRand", "signalConst", "signalSuperposed" or
# for free vibration choose "signalNone"

# external dynamic force acting on the system
freq = 5
force_dynamic = load_type("signalSin", array_time, len(z), freq, force_static)

dynamic_analysis = DynamicAnalysis(mdof_solver, force_dynamic, time_parameters)
dynamic_analysis.solve()
dynamic_analysis.plot_selected_time_step(0.75)
dynamic_analysis.plot_result_dof(3, 'acceleration')
dynamic_analysis.animate_time_history()
