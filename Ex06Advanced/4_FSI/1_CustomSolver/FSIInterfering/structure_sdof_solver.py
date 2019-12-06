# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        SDoF system solver using direct time integration - Generalized-Alpha Scheme

Author: mate.pentek@tum.de 
        Implementation adapted from I. Hanzlicek (2014). 
        
        Original implementation by M. Andre described in: Formulation of the Generalized-Alpha
        method for LAGRANGE. Technical Report, Chair of Structural Analysis @TUM, 2012.
        
        See J. Chung, G.M. Hulbert: A time integration algorithm for structural dynamics
        with improved numerical dissipation: the generalized-aplha mehod. ASME J. Appl. 
        Mech., 60:371-375,1993. 
         
Description: This is a solver for direct numerical time integration for SDoF systems.
        It assumes a linear SDOF with a Generalized alpha scheme with fixed dt.
        
Note:   Mimics the structure of the code for analysis in KratosMultiphysics.

Created on:  15.11.2015
Last update: 06.12.2019
'''
# ===============================================================================

from math import pi


def CreateSolver(structure_settings):
    return StructureSDoF(structure_settings)


class StructureSDoF(object):
    # """ Direct time integration of linear SDOF (Generalized-alpha method). """
    # This class takes the undamped eigenfrequency (f, [Hz]) as input
    def __init__(self, structure_settings):

        self.time = structure_settings["problem_data"]["start_time"].GetDouble(
        )
        self.dt = structure_settings["problem_data"]["time_step"].GetDouble()
        self.step = 0

        self.absolute_position = structure_settings["model_data"]["absolute_position"].GetDouble(
        )

        # structure moment of inertia, damping and spring stiffness
        self.m = structure_settings["model_data"]["mass"].GetDouble()
        f = structure_settings["model_data"]["eigen_freq"].GetDouble()
        omega = 2 * pi * f
        self.k = omega**2 * self.m
        zeta = structure_settings["model_data"]["zeta"].GetDouble()
        # zeta is the damping ratio - avoid this line, input explicitly
        self.b = 2.0 * (self.m * self.k)**0.5 * zeta

        p_inf = structure_settings["model_data"]["rho_inf"].GetDouble()
        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        self.alpha_m = (2.0 * p_inf - 1.0) / (p_inf + 1.0)
        self.alpha_f = p_inf / (p_inf + 1.0)
        self.beta = 0.25 * (1 - self.alpha_m + self.alpha_f)**2
        self.gamma = 0.5 - self.alpha_m + self.alpha_f
        # coefficients for LHS
        self.a1h = (1.0 - self.alpha_m) / (self.beta * self.dt**2)
        self.a2h = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alpha_f
        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alpha_m - 2.0 * self.beta) / (2.0 * self.beta)
        # coefficients for damping
        self.a1b = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.dt)
        self.a2b = (1.0 - self.alpha_f) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alpha_f) * \
            (0.5 * self.gamma / self.beta - 1.0) * self.dt
        # coefficient for stiffness
        self.a1k = -1.0 * self.alpha_f
        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt
        # coefficients for acceleration update
        self.a1a = self.a1v / (self.dt * self.gamma)
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)

        # initial displacement, velocity and acceleration
        self.u0 = structure_settings["model_data"]["initial"]["displacement"].GetDouble(
        )
        self.v0 = structure_settings["model_data"]["initial"]["velocity"].GetDouble(
        )
        self.a0 = structure_settings["model_data"]["initial"]["acceleration"].GetDouble(
        )

        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0

        # moment from a previous time step (initial moment)
        self.f0 = self.m * self.a0 + self.b * self.v0 + self.k * self.u0
        self.f1 = self.m * self.a0 + self.b * self.v0 + self.k * self.u0

        self.dof_type = structure_settings["model_data"]["dof_type"].GetString(
        )

        # external force
        self.ext_force = None

        # filename
        self.filename = "sdof_" + structure_settings["problem_data"]["problem_name"].GetString(
        ).lower() + "_" + self.dof_type.lower() + ".dat"

        # output
        self.support_output = open(self.filename, 'w')
        self.support_output.write(
            "# (1): time [s] (2): displacement/rotation [m/rad] (3): support force/ moment [N]/[N m]\n")
        self.support_output.flush()

    def GetDisplacement(self):
        return self.u1

    def SetDisplacement(self, displacement):
        self.u1 = displacement

    def Predict(self):
        return 2.0 * self.u1 - self.u0

    def PrintSupportOutput(self):
        self.support_output.write(
            str(self.time) + " " + str(self.u1) + " " + str(self.k * self.u1) + "\n")
        self.support_output.flush()

    def SetExternalForce(self, ext_force):
        # update external force
        self.ext_force = ext_force

    def _AssembleLHS(self):
        return self.a1h * self.m + self.a2h * self.b + self.a3h * self.k

    def _AssembleRHS(self):
        RHS = self.m * (self.a1m * self.u0 + self.a2m *
                        self.v0 + self.a3m * self.a0)
        RHS += self.b * (self.a1b * self.u0 + self.a2b *
                         self.v0 + self.a3b * self.a0)

        f = (1.0 - self.alpha_f) * self.ext_force + self.alpha_f * self.f0
        RHS += self.a1k * self.k * self.u0 + f

        return RHS

    def FinalizeSolutionStep(self):
        # update v1, a1
        self.v1 = self.a1v * (self.u1 - self.u0) + \
            self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + \
            self.a2a * self.v0 + self.a3a * self.a0

    def SolveSolutionStep(self):
        # assemble LHS
        LHS = self._AssembleLHS()

        # assemble RHS
        RHS = self._AssembleRHS()

        # solve
        # sys of eq reads: LHS * u1 = RHS
        self.u1 = RHS / LHS

    def _IncrementTimeStep(self):
        # update angular displacement, velocity and acceleration
        self.u0 = self.u1
        self.v0 = self.v1
        self.a0 = self.a1
        # update the moment
        self.f0 = self.f1

    def Initialize(self):
        print("SDoF: Initialize() called, needs to be implemented")
        pass

    def Finalize(self):
        print("SDoF: Finalize() called, needs to be implemented")
        pass

    def InitializeSolutionStep(self):
        print("SDoF: InitializeSolutionStep() called, needs to be implemented")
        pass

    def AdvanceInTime(self, time):
        self._IncrementTimeStep()
        self.step += 1
        self.time = self.time + self.dt

        return self.time

    def OutputSolutionStep(self):
        self.PrintSupportOutput()

    def GetPosition(self):
        # we are using absolute displacements
        # so the current position is the last displacement
        return self.absolute_position + self.u0
