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
import time

from source.analysis_type import*
from source.load_type import*
from source.custom_files import *
from source.mdof_solver import *

parameter_file = open('cosim_mdof_cantilever_shear_2d_parameters.json', 'r')
Parameters = json.loads(parameter_file.read())
solver_settings = Parameters['solver_settings']
mdof_solver = MDoFSolver(solver_settings, 1)
mdof_solver_modal = MDoFSolver(solver_settings, 1)
mdof_solver_modal_sdof = MDoFSolver(solver_settings, 1)

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
eigen_vector_norm = eigenvalue_analysis.eigenform
eigen_values = eigenvalue_analysis.frequency
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
freq = 10
force_dynamic = load_type("signalSin", array_time, len(z), freq, force_static_const)

# reference solution 
time_start = time.time()
dynamic_analysis = DynamicAnalysis(mdof_solver, force_dynamic, time_parameters)
dynamic_analysis.solve()
displacement_ref = dynamic_analysis.displacement
time_ref = time.time() - time_start


# dynamic analysis in modal coordinates
# solves the MDOF system itself : redundant : needs to be removed probably 
# modes_considered = 2
# # define here the modes to be considered for modal analysis
# eigen_vector_norm = eigen_vector_norm[:, : modes_considered] 
# eigen_values = eigen_values[:modes_considered]
# It might be difficult to consider less modes in here 

time_start = time.time()
dynamic_analysis_modal = DynamicAnalysisModal(mdof_solver_modal, force_dynamic, time_parameters, eigen_vector_norm, eigen_values)
dynamic_analysis_modal.solve()
displacement_modal = dynamic_analysis_modal.displacement
time_modal = time.time() - time_start


# dynamic analysis in modal coordinates woth SDOF 
# conoverts to the modal coordinates and solves sdof system and then transform back to the actual coordinates 
# TODO : check for a better way to construct SDOF solver 
modes_considered = 2
# define here the modes to be considered for modal analysis
eigen_vector_norm = eigen_vector_norm[:, : modes_considered] 
eigen_values = eigen_values[:modes_considered]

time_start = time.time()
dynamic_analysis_modal_sdof = DynamicAnalysisModalSDOF(mdof_solver_modal_sdof, force_dynamic, time_parameters, eigen_vector_norm, eigen_values)
dynamic_analysis_modal_sdof.solve()
displacement_modal_sdof = dynamic_analysis_modal_sdof.displacement
time_modal_sdof = time.time() - time_start


print('time for the reference solution', time_ref)

rms_error_modal = np.sqrt(((displacement_modal - displacement_ref) ** 2).mean())
print('rms_error in displacement for modal coordinates', rms_error_modal)
print('time for the modal coordinates', time_modal)

rms_error_modal_sdof = np.sqrt(((displacement_modal_sdof - displacement_ref) ** 2).mean())
print('rms_error in displacement for modal coordinates with sdof', rms_error_modal_sdof)
print('time for the modal coordinates', time_modal_sdof)
 

#dynamic_analysis.plot_selected_time_step(0.75)
dynamic_analysis.plot_result_at_dof(10, 'displacement')
dynamic_analysis_modal_sdof.plot_result_at_dof(10, 'displacement')
#dynamic_analysis.animate_time_history()


# this is the reference solution with time history analysis 
