# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        3DoF system solver using direct time integration - Generalized-Alpha Scheme

Author: mate.pentek@tum.de 
        
         
Description: This is a solver for direct numerical time integration for 3DoF systems
        as a combination of SDoFs.

        It assumes a linear SDoFs with a Generalized alpha scheme with fixed dt.
        
Note:   Mimics the structure of the code for analysis in KratosMultiphysics.

Created on:  15.11.2015
Last update: 13.11.2018
'''
# ===============================================================================

import KratosMultiphysics

from math import pi
import structure_sdof_solver


def CreateSolver(structure_settings):
    return Structure3DoF(structure_settings)


class Structure3DoF(object):
    def __init__(self, structure_settings):

        self.solvers = []
        # size() for Kratos.Prarameter array
        for i in range(structure_settings["model_data"]["dof_type"].size()):

            solver_settings = structure_settings.Clone()

            solver_settings["model_data"].RemoveValue("absolute_position")
            solver_settings["model_data"].AddEmptyValue("absolute_position").SetDouble(
                structure_settings["model_data"]["absolute_position"][i].GetDouble())

            solver_settings["model_data"].RemoveValue("mass")
            solver_settings["model_data"].AddEmptyValue("mass").SetDouble(
                structure_settings["model_data"]["mass"][i].GetDouble())

            solver_settings["model_data"].RemoveValue("eigen_freq")
            solver_settings["model_data"].AddEmptyValue("eigen_freq").SetDouble(
                structure_settings["model_data"]["eigen_freq"][i].GetDouble())

            solver_settings["model_data"].RemoveValue("zeta")
            solver_settings["model_data"].AddEmptyValue("zeta").SetDouble(
                structure_settings["model_data"]["zeta"][i].GetDouble())

            solver_settings["model_data"].RemoveValue("rho_inf")
            solver_settings["model_data"].AddEmptyValue("rho_inf").SetDouble(
                structure_settings["model_data"]["rho_inf"][i].GetDouble())

            solver_settings["model_data"].RemoveValue("initial")
            solver_settings["model_data"].AddEmptyValue("initial")

            solver_settings["model_data"]["initial"].AddEmptyValue("displacement").SetDouble(
                structure_settings["model_data"]["initial"]["displacement"][i].GetDouble())
            solver_settings["model_data"]["initial"].AddEmptyValue("velocity").SetDouble(
                structure_settings["model_data"]["initial"]["velocity"][i].GetDouble())
            solver_settings["model_data"]["initial"].AddEmptyValue("acceleration").SetDouble(
                structure_settings["model_data"]["initial"]["acceleration"][i].GetDouble())

            solver_settings["model_data"].RemoveValue("dof_type")
            solver_settings["model_data"].AddEmptyValue("dof_type").SetString(
                structure_settings["model_data"]["dof_type"][i].GetString())

            self.solvers.append(
                structure_sdof_solver.CreateSolver(solver_settings))

    def GetPosition(self):
        return [solver.GetPosition() for solver in self.solvers]

    def GetDisplacement(self):
        return [solver.GetDisplacement() for solver in self.solvers]

    def SetDisplacement(self, displacement):
        for i in range(len(self.solvers)):
            self.solvers[i].SetDisplacement(displacement[i])

    def Predict(self):
        for solver in self.solvers:
            solver.Predict()

    def PrintSupportOutput(self):
        for solver in self.solvers:
            solver.PrintSupportOutput()

    def SetExternalForce(self, ext_force):
        for i in range(len(self.solvers)):
            self.solvers[i].SetExternalForce(ext_force[i])

    def FinalizeSolutionStep(self):
        for solver in self.solvers:
            solver.FinalizeSolutionStep()

    def SolveSolutionStep(self):
        for solver in self.solvers:
            solver.SolveSolutionStep()

    def Initialize(self):
        print("3DoF: Initialize() called, needs to be implemented")
        for solver in self.solvers:
            solver.Initialize()
        pass

    def Finalize(self):
        print("3DoF: Finalize() called, needs to be implemented")
        for solver in self.solvers:
            solver.Finalize()
        pass

    def InitializeSolutionStep(self):
        print("3DoF: InitializeSolutionStep() called, needs to be implemented")
        for solver in self.solvers:
            solver.InitializeSolutionStep()
        pass

    def AdvanceInTime(self, time):
        times = []
        for solver in self.solvers:
            times.append(solver.AdvanceInTime(time))

        # maybe check first it all the same
        return times[0]

    def OutputSolutionStep(self):
        for solver in self.solvers:
            solver.OutputSolutionStep()
