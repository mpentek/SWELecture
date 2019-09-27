# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Sample MDoF shear run

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de, catharina.czech@tum.de, peter.kupas@tum.de
    
Note: ...

Created on:  22.11.2017
Last update: 27.09.2019
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
# # Generating wind profile
# wind_velocity = 28.7  # m/s
# velocity_vector = 1.05 * wind_velocity * pow(z/10, 0.2)
# wind_force = 0.5 * 1.2 * velocity_vector**2 * 30
# a constant force on all dof sum of which is 1 
force_static_const = 1000 * np.ones(np.shape(z))
# a unit force at the top dof 
force_static_top = np.zeros(np.shape(z))
force_static_top[-1] = 1

static_analysis = StaticAnalysis(mdof_solver.model)
static_analysis.solve(force_static_top)
#static_analysis.plot_solve_result()

# =========================eigenvalue analysis ==========================
eigenvalue_analysis = EigenvalueAnalysis(mdof_solver.model)
eigenvalue_analysis.solve()
#eigenvalue_analysis.plot_selected_eigenmode(1)
#eigenvalue_analysis.plot_selected_first_n_eigenmodes(3)
#eigenvalue_analysis.animate_selected_eigenmode(1)

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
force_dynamic = load_type("signalSin", array_time, len(z), freq, force_static_const)

dynamic_analysis = DynamicAnalysis(mdof_solver, force_dynamic, time_parameters)
dynamic_analysis.solve()
displacement_ref = dynamic_analysis.displacement
print(displacement_ref)
#dynamic_analysis.plot_selected_time_step(0.75)
#dynamic_analysis.plot_result_at_dof(3, 'displacement')
#dynamic_analysis.animate_time_history()


#AK: Notes
# this is the reference solution with time history analysis 
# the load is the uniform sinusoidal with magnitude 1000 at each of the fiibe dof anf frequency is 5 
