# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        Sample MDoF generic run

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

parameter_file = open('cosim_mdof_generic_parameters.json', 'r')
Parameters = json.loads(parameter_file.read())
solver_settings = Parameters['solver_settings']
mdof_solver = MDoFSolver(solver_settings, 1)

# ============overwrite read-in values for a custom model===========
# # import or read-in -> sample data for a generic highrise
# import example_models.mdof_model_shear_highrise as m_highrise
# # stiffness matrix
# mdof_solver.model.k = m_highrise.get_stiffness()
# # # mass matrix
# mdof_solver.model.m = m_highrise.get_mass()
# # # height coordinates
# z = m_highrise.get_height_coordinates()
# z= z[:-1]
# # damping
# mdof_solver.model.b = np.zeros(np.shape(m_highrise.get_mass()))

# mdof_solver.model.nodal_coordinates["x0"] = np.zeros(len(z))
# mdof_solver.model.nodal_coordinates["y0"] = z
# mdof_solver.model.u0 = np.zeros(len(z))
# mdof_solver.model.v0 = np.zeros(len(z))
# mdof_solver.model.a0 = np.zeros(len(z))
# mdof_solver.model.f0 = np.zeros(len(z))
# mdof_solver.scheme.buffer = np.zeros((4,
#                                 mdof_solver.scheme.buffer_size,
#                                 len(z)))
# mdof_solver.nr_of_dofs = len(z)
# mdof_solver.buffer = np.zeros((3, mdof_solver.buffer_size, mdof_solver.nr_of_dofs))

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
dynamic_analysis.animate_time_history()
