#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS17-18 
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek
        
        A python script to run a FSI simulation

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de
      
Note:   This script is tuned to work with the precompiled Kratos 5.2.0, so only basic python
        functionalities (like python math syntax YES, but numpy NOT) are permitted.

Created on:  16.01.2018
Last update: 16.01.2018
'''
#===============================================================================
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# for FSI - mesh moving
from KratosMultiphysics.ALEApplication import *

# for structural sdof solver 
from structure_sdof_class_implementation import *

######################################################################################
######################################################################################
######################################################################################

## Parse the ProjectParameters
parameter_file = open("ProjectParametersFSIVortexShedding.json",'r')
ProjectParameters = Parameters( parameter_file.read())


# importing structural and fsi parameters
parameter_file_fsi = open("ProjectParametersFSICSD.json",'r')
ProjectParameterFsi = Parameters( parameter_file_fsi.read())

## Get echo level and parallel type
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

## Import KratosMPI if needed
if (parallel_type == "MPI"):
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.TrilinosApplication import *

## Fluid model part definition
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

## Solver construction
import python_solvers_wrapper_ale
solver = python_solvers_wrapper_ale.CreateSolver(main_model_part, ProjectParameters)

solver.AddVariables()

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()

## Initialize GiD  I/O
output_post  = ProjectParameters.Has("output_configuration")
if (output_post == True):
    if (parallel_type == "OpenMP"):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                      ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                      ProjectParameters["output_configuration"])
    elif (parallel_type == "MPI"):
        from gid_output_process_mpi import GiDOutputProcessMPI
        gid_output = GiDOutputProcessMPI(solver.GetComputingModelPart(),
                                         ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                         ProjectParameters["output_configuration"])

    gid_output.ExecuteInitialize()

## Creation of Kratos model (build sub_model_parts or submeshes)
FluidModel = Model()
FluidModel.AddModelPart(main_model_part)

## Get the list of the skin submodel parts in the object Model
for i in range(ProjectParameters["solver_settings"]["skin_parts"].size()):
    skin_part_name = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
    FluidModel.AddModelPart(main_model_part.GetSubModelPart(skin_part_name))

## Get the list of the no-skin submodel parts in the object Model (results processes and no-skin conditions)
for i in range(ProjectParameters["solver_settings"]["no_skin_parts"].size()):
    no_skin_part_name = ProjectParameters["solver_settings"]["no_skin_parts"][i].GetString()
    FluidModel.AddModelPart(main_model_part.GetSubModelPart(no_skin_part_name))

## Get the list of the initial conditions submodel parts in the object Model
for i in range(ProjectParameters["initial_conditions_process_list"].size()):
    initial_cond_part_name = ProjectParameters["initial_conditions_process_list"][i]["Parameters"]["model_part_name"].GetString()
    FluidModel.AddModelPart(main_model_part.GetSubModelPart(initial_cond_part_name))

## Get the gravity submodel part in the object Model
for i in range(ProjectParameters["gravity"].size()):
    gravity_part_name = ProjectParameters["gravity"][i]["Parameters"]["model_part_name"].GetString()
    FluidModel.AddModelPart(main_model_part.GetSubModelPart(gravity_part_name))

## Print model_part and properties
if (echo_level > 1) and ((parallel_type == "OpenMP") or (mpi.rank == 0)):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

## Processes construction
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
# Note 1: gravity is firstly constructed. Outlet process might need its information.
# Note 2: conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
list_of_processes =  process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["gravity"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(FluidModel).ConstructListOfProcesses( ProjectParameters["auxiliar_process_list"] )

if (echo_level > 1) and ((parallel_type == "OpenMP") or (mpi.rank == 0)):
    for process in list_of_processes:
        print(process)

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

## Solver initialization
solver.Initialize()

# ALE boundary conditions
# PMT note: we need model part from skin parts and sometimes from no slin parts as well
# for now temporary workaround
print("ALE BC for Fluid Model Part")
ale_fixed_boundaries = ["Outlet2D_outlet",
                        "Slip2D_wall",
                        "NoSlip2D_structure", 
                        "VelocityConstraints2D_inlet"]
#for i in range(ProjectParameters["solver_settings"]["skin_parts"].size())
#    skin_part = ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
for skin_part in ale_fixed_boundaries:
    print("skin_part: " + skin_part)
    for node in main_model_part.GetSubModelPart(skin_part).Nodes:
        node.Fix(MESH_DISPLACEMENT_X)
        node.Fix(MESH_DISPLACEMENT_Y)
        node.Fix(MESH_DISPLACEMENT_Z)
        
## Stepping and time settings
# AK modified 
# delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
Dt = ProjectParameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()
Nsteps = (end_time-ProjectParameters["problem_data"]["start_time"].GetDouble())/Dt

# AK added 

def RotateMirrorDisplace(theta, disp, nodes, dof_type):
    if dof_type == "DISPLACEMENT_X":
        for node in nodes:
            node.SetSolutionStepValue(MESH_DISPLACEMENT_X,disp)
            
    if dof_type == "DISPLACEMENT_Y":
        for node in nodes:
            node.SetSolutionStepValue(MESH_DISPLACEMENT_Y,disp)

def ExtractForce(nodes,dof_type):
    fx = 0.0
    fy = 0.0
    if dof_type == "DISPLACEMENT_X":
        for node in nodes:
            reaction = node.GetSolutionStepValue(REACTION, 0)
            fx += -reaction[0]
        return fx
    if dof_type == "DISPLACEMENT_Y":
        for node in nodes:
            reaction = node.GetSolutionStepValue(REACTION, 0)
            fy += -reaction[1]
        return fy
#AK added
# assigns mesh velocity to fluid - different then one got from Ivan - that was for monolithic        
def SetMeshVelocityToFluid(nodes):
    for node in nodes:
        u_mesh = node.GetSolutionStepValue(MESH_VELOCITY, 0) #get mesh velocity at current step
        node.SetSolutionStepValue(VELOCITY, 0, u_mesh) # assign it to fluid velocity at current step


#Structure properties
mass = ProjectParameterFsi["structure_data"]["mass"].GetDouble()
eigen_freq = ProjectParameterFsi["structure_data"]["eigen_freq"].GetDouble()
zeta = ProjectParameterFsi["structure_data"]["zeta"].GetDouble()
dof_type = ProjectParameterFsi["structure_data"]["dof_type"].GetString()

# Parameters:							   (Dt, M,    zeta, f,         B,     pInf, u0,  v0,  a0,  filename):
rho_inf = ProjectParameterFsi["structure_data"]["rho_inf"].GetDouble()
initial_disp = ProjectParameterFsi["structure_data"]["initial_disp"].GetDouble()
initial_vel = ProjectParameterFsi["structure_data"]["initial_vel"].GetDouble()
initial_acc = ProjectParameterFsi["structure_data"]["initial_acc"].GetDouble()
output_filename = ProjectParameterFsi["structure_data"]["output_filename"].GetString() # AK check if required 

structure_dispX_solver = StructureSDOF_Freq(Dt, mass, zeta, eigen_freq, rho_inf, initial_disp, initial_vel, initial_acc, output_filename)


# Getting the submodel structure for FSI iteration
struc_submodel_string = ProjectParameterFsi["structure_data"]["model_part_name"].GetString()
Structure = main_model_part.GetSubModelPart(struc_submodel_string)


# FSI iteration parameters
abs_residual = ProjectParameterFsi["FSI_parameters"]["abs_residual"].GetDouble()
rel_residual = ProjectParameterFsi["FSI_parameters"]["rel_residual"].GetDouble()
relax_coef =  ProjectParameterFsi["FSI_parameters"]["relax_coef"].GetDouble()
max_FSI_iteration = ProjectParameterFsi["FSI_parameters"]["max_FSI_iteration"].GetDouble()


if (output_post == True):
    gid_output.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Writing the full ProjectParameters file before solving
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 0):
    f = open("ProjectParametersOutput.json", 'w')
    f.write(ProjectParameters.PrettyPrintJsonString())
    f.close()

# --------CALCULATION LOOP-----------------------------------------------------

time = 0.0
out = 0
step = 0

while(time <= end_time):

    step += 1
    time = time + Dt
    main_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)
    sys.stdout.flush()
        
    if(step >= 3):
        for process in list_of_processes:
            process.ExecuteInitializeSolutionStep()

        gid_output.ExecuteInitializeSolutionStep()

        k = 0
        # predict structure rotation
        pitching_angle = 0
        displacement = structure_dispX_solver.predictRotation()

        # solve the fluid problem
        RotateMirrorDisplace(pitching_angle, displacement, Structure.Nodes, dof_type)
        # Solve Mesh
        print("## Solving the mesh")
        sys.stdout.flush()
        solver.GetMeshMotionSolver().Solve()
        # Apply Mesh Velocity from mesh solver to fluid solver
        print("## Setting mesh velocity to fluid")
        sys.stdout.flush()
        SetMeshVelocityToFluid(Structure.Nodes)
        # Solve Fluid
        print("## Solving the fluid")
        sys.stdout.flush()
        solver.GetFluidSolver().Solve()
        # solve the structure problem
        reaction = 0.0
        reaction += ExtractForce(Structure.Nodes,dof_type)
        structure_dispX_solver.solveStructure(reaction)
        
        # compute the residual
        initial_residual2 = structure_dispX_solver.getRotation() - displacement
        old_residual2 = initial_residual2
        residual2 = old_residual2
        while not (abs(residual2) < abs_residual or abs(residual2) < abs(initial_residual2) * rel_residual):
            #relax_coef = 1.0
            # update structure rotation
            pitching_angle = 0
            displacement += relax_coef * residual2
            k = k + 1

            # solve the fluid problem
            RotateMirrorDisplace(pitching_angle, displacement, Structure.Nodes, dof_type)
             # Solve Mesh
            print("## Solving the mesh")
            sys.stdout.flush()
            solver.GetMeshMotionSolver().Solve()
            # Apply Mesh Velocity from mesh solver to fluid solver
            print("## Setting mesh velocity to fluid")
            sys.stdout.flush()
            SetMeshVelocityToFluid(Structure.Nodes)
            # Solve Fluid
            print("## Solving the fluid")
            sys.stdout.flush()
            solver.GetFluidSolver().Solve()

            # solve the structure problem
            reaction = 0.0
            reaction += ExtractForce(Structure.Nodes,dof_type)
            structure_dispX_solver.solveStructure(reaction)

            # compute the residual
            old_residual2 = residual2
            residual2 = structure_dispX_solver.getRotation() - displacement
            print('Iteration[', k, ']: relax_coef = ', relax_coef, 'displacement rel. residual2 = ',abs(residual2 / initial_residual2), ' abs. residual2 = ',abs(residual2))
            if k >= max_FSI_iteration: 
                break
                
        structure_dispX_solver.printSupportOutput(time)
        structure_dispX_solver.incrementTimeStep()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    if (output_post == True):
        gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    if (gid_output.IsOutputStep()) and (output_post == True):
        gid_output.PrintOutput()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

for process in list_of_processes:
    process.ExecuteFinalize()

if (output_post == True):
    gid_output.ExecuteFinalize()
