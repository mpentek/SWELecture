# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        A python script to run a FSI simulation

Author: philipp.bucher@tum.de, a.winterstein@tum.de, mate.pentek@tum.de, anoop.kodakkal@tum.de
      
Note:   This script is tuned to work with the precompiled Kratos 6.0.0, so only basic python
        functionalities (like python math syntax YES, but numpy NOT) are permitted.

        Main Script for FSI with Kratos Mutliphysics

        ##

        This script is intended to be modified. Each solver can be imported and used as "BlackBox"

        Chair of Structural Analysis, Technical University of Munich
        All rights reserved

        This example is based on the dissertation of Daniel Mok
        "Partitionierte Lösungsansätze in der Strukturdynamik und der Fluid-Struktur-Interaktion"
        Chapter 7.3 "Flexible Klappe in Kanalströmung mit Einschnürung"

Created on:  16.01.2018
Last update: 13.11.2018
'''
# ===============================================================================


# ----- Importing the modules -----
import KratosMultiphysics
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
#from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import KratosMultiphysics.MeshMovingApplication.mesh_moving_analysis as MeshMovingAnalysis

# Import the "BlackBox" Solvers
import structure_3dof_solver as structure_solver

# here auxiliary functions e.g. for relaxation are declared
import fsi_utilities
import other_utilities

TOLERANCE = 1e-12
DIMENSION = 2

fluid_model = KratosMultiphysics.Model()
#structural_model = KratosMultiphysics.Model()

fluid_project_params_file_name = "ProjectParametersCFD.json"
with open(fluid_project_params_file_name, 'r') as parameter_file:
    parameters_fluid = KratosMultiphysics.Parameters(parameter_file.read())
structural_project_params_file_name = "ProjectParametersCSD.json"
with open(structural_project_params_file_name, 'r') as parameter_file:
    parameters_structure = KratosMultiphysics.Parameters(parameter_file.read())

'''
# --------------------------------------------------------
# ----- Setting up and initializing the Fluid Solver -----
# --------------------------------------------------------
'''
fluid_solver = FluidDynamicsAnalysis(fluid_model, parameters_fluid)

fluid_solver.Initialize()

fluid_model_part = fluid_model["FluidModelPart"]

print("======================================================================")
print("||||||||||||||||||||||| SETTING UP FLUID DONE ||||||||||||||||||||||||")
print("======================================================================")

'''
# -------------------------------------------------------------
# ----- Setting up and initializing the Structural Solver -----
# -------------------------------------------------------------
'''
structural_solver = structure_solver.CreateSolver(parameters_structure)

structural_solver.Initialize()

#structural_model_part = structural_model["Structure"]

print("======================================================================")
print("||||||||||||||||| SETTING UP STRUCTURAL DYNAMICS DONE ||||||||||||||||")
print("======================================================================")

'''
# ------------------------------------------------------
# ----- Setting up the FSI-related functionalities -----
# ------------------------------------------------------
# '''

fsi_project_params_file_name = "ProjectParametersFSI.json"
with open(fsi_project_params_file_name, 'r') as parameter_file:
    parameters_fsi = KratosMultiphysics.Parameters(parameter_file.read())

# ----- Setting up the time parameters -----
start_time = parameters_fsi["problem_data"]["start_time"].GetDouble()
end_time = parameters_fsi["problem_data"]["end_time"].GetDouble()
delta_time = parameters_fsi["problem_data"]["time_step"].GetDouble()

num_steps = int((end_time - start_time) / delta_time)

round_val = other_utilities.TimeRoundValue(delta_time)

time = start_time
step = 0

# ----- Setting up the FSI Parameters -----
# ---------------
# ----- ALE -----
# ---------------
fluid_solver._GetSolver().GetMeshMotionSolver().SetEchoLevel(
    0)  # Fix until the source of the prints is found

# -------------------
# ----- Mapping -----
# -------------------
# mapper = fsi_utilities.CreateMapper(fluid_model_part, structural_model_part, parameters_fsi["coupling_settings"]["mapper_settings"])
mapper = fsi_utilities.CreateMapper(
    fluid_model_part, parameters_fsi["coupling_settings"]["mapper_settings"])


# -------------------
# --- Convergence ---
# -------------------
convergence_accelerator = fsi_utilities.CreateConvergenceAccelerator(
    parameters_fsi["coupling_settings"]["convergence_accelerator_settings"])
relaxation_coefficient = convergence_accelerator.rel_coef_initial


##
##
##


# -------- OUTPUTS --------------------------------

# # Result output file names
# iterationFileName = "Iterations_Aitken.dat"
# cornerDispFileName = "cornerNode_XDisp_Aitken.dat"

# # clean up the Iterations.dat
# open(iterationFileName, "w").close()

# iterationFile = open(iterationFileName, "a")
# iterationFile.write("# time, iterations, relaxation coefficients" + "\n")
# iterationFile.close()

# # clean up the cornerNodeXDisp.dat
# open(cornerDispFileName, "w").close()

# cornerDispFile = open(cornerDispFileName, "a")
# cornerDispFile.write("# time, corner node displacement X" + "\n")
# cornerDispFile.close()


##
##
##


print("======================================================================")
print("|||||||||||||||||||||||| SETTING UP FSI DONE |||||||||||||||||||||||||")
print("======================================================================")

#file_writer = other_utilities.FileWriter("generic_fsi.dat", ["Time", structural_solver.dof_type, "Coupling_Iterations"])
#tip_node = structural_model_part.GetNode(1)

# ----- Solving the problem (time integration) -----
while(time <= end_time):
    new_time_fluid = fluid_solver._GetSolver().AdvanceInTime(time)
    #new_time_structure = structural_solver._GetSolver().AdvanceInTime(time)
    new_time_structure = structural_solver.AdvanceInTime(time)

    fluid_solver._GetSolver().Predict()
    # structural_solver._GetSolver().Predict()
    structural_solver.Predict()

    fluid_solver.InitializeSolutionStep()
    structural_solver.InitializeSolutionStep()

    time = time + delta_time
    if abs(time - new_time_fluid) > TOLERANCE:
        raise Exception("Fluid has wrong time!")
    if abs(time - new_time_structure) > TOLERANCE:
        raise Exception("Structure has wrong time!")
    step += 1

    print("\n--- Step =", step, "/", num_steps, "---")
    print("--- Time =", round(time, round_val), "/", end_time, "---")

    residual = 1
    # from the nodes in the structure (origin) on the interface
    # old_displacements = fsi_utilities.GetDisplacements(mapper.origin_interface.Nodes, DIMENSION)
    old_displacements = fsi_utilities.GetDisplacements(structural_solver)

    num_inner_iter = 1
    # Inner FSI Loop (executed once in case of explicit coupling)
    for k in range(convergence_accelerator.max_iter):

        # Apply Dirichlet B.C.'s from structural solver to mesh solver
        # map from the nodes in the structure (origin) on the interface
        # KratosStructuralMechanics.POINT_LOAD
        # to the interface nodes on the fluid (destination)
        # KratosMultiphysics.MESH_DISPLACEMENT
        # fsi_utilities.DisplacementToMesh(mapper)
        fsi_utilities.DisplacementToMesh(
            mapper, old_displacements, structural_solver)

        # Mesh and Fluid are currently solved independently, since the ALE solver does not copy the mesh velocity
        # Solve Mesh
        fluid_solver._GetSolver().SolveSolutionStep()

        # PMT: in the previous function called contained
        # solver_fluid.GetMeshMotionSolver().Solve()

        # PMT: done implicitly by the ALE solver for the model part set
        # in the project parameters: "ale_boundary_parts" : ["NoSlip2D_structure"],
        # Apply Mesh Velocity from mesh solver to fluid solver
        # SetMeshVelocityToFluid(structure_interface.Nodes)

        # Apply Neumann B.C.'s from fluid solver to structural solver
        # inverse map from the nodes in the fluid (destination) on the interface
        # KratosMultiphysics.REACTION
        # to the interface nodes on the structure (origin) with sign swap
        # KratosStructuralMechanics.POINT_LOAD
        fsi_utilities.NeumannToStructure(mapper, structural_solver, True)

        # # Solver Structure
        # structural_solver._GetSolver().SolveSolutionStep()
        structural_solver.SolveSolutionStep()

        # Convergence Checking (only for implicit coupling)
        if convergence_accelerator.max_iter > 1:
            #displacements = fsi_utilities.GetDisplacements(mapper.origin_interface.Nodes, DIMENSION)
            displacements = fsi_utilities.GetDisplacements(structural_solver)

            # Compute Residual
            old_residual = residual
            residual = convergence_accelerator.CalculateResidual(displacements, old_displacements)

            if (fsi_utilities.Norm(residual) <= convergence_accelerator.res_rel_tol):
                #fsi_utilities.SetDisplacements(displacements, mapper.origin_interface.Nodes, DIMENSION)
                fsi_utilities.SetDisplacements(displacements, structural_solver)

                print("******************************************************")
                print("************ CONVERGENCE AT INTERFACE ACHIEVED *******")
                print("******************************************************")
                break  # TODO check if this works bcs it is nested
            else:
                relaxation_coefficient = convergence_accelerator.ComputeRelaxationCoefficient(relaxation_coefficient, residual, old_residual, k)
                relaxed_displacements = convergence_accelerator.CalculateRelaxedSolution(relaxation_coefficient, old_displacements, residual)
                old_displacements = relaxed_displacements
                #fsi_utilities.SetDisplacements(relaxed_displacements,  mapper.origin_interface.Nodes, DIMENSION)
                fsi_utilities.SetDisplacements(relaxed_displacements, structural_solver)

                # iterations.append(round(relaxation_coefficient, 3))
                num_inner_iter += 1

            if (k+1 >= convergence_accelerator.max_iter):
                print("######################################################")
                print("##### CONVERGENCE AT INTERFACE WAS NOT ACHIEVED ######")
                print("######################################################")

            print("==========================================================")
            print("COUPLING RESIDUAL = ", fsi_utilities.Norm(residual))
            print("COUPLING ITERATION = ", k+1, "/", convergence_accelerator.max_iter)
            print("RELAXATION COEFFICIENT = ", relaxation_coefficient)
            print("==========================================================")

    # iterationFile = open(iterationFileName, "a")
    # iterationFile.write(str(time) + "\t" +str(k) + "\t" + str(iterations)[1:-1] + "\n")
    # iterationFile.close()

    # cornerDispFile = open(cornerDispFileName, "a")
    # for node in fluid_interface.Nodes:
    #     if (abs(node.X0+15) < 1e-4 and abs(node.Y0-190) < 1e-4):
    #         cornerDispFile.write(str(time) + "\t" + str(node.GetSolutionStepValue(DISPLACEMENT_X,0)) + "\n")
    # cornerDispFile.close()

    fluid_solver.FinalizeSolutionStep()
    structural_solver.FinalizeSolutionStep()

    fluid_solver.OutputSolutionStep()
    structural_solver.OutputSolutionStep()

    #disp = tip_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
    disp = structural_solver.GetDisplacement()
    #file_writer.WriteToFile([time, disp[0], disp[1], disp[2], num_inner_iter])
    #file_writer.WriteToFile([time, disp, num_inner_iter])

# TIME LOOP END
fluid_solver.Finalize()
structural_solver.Finalize()

# file_writer.CloseFile()
