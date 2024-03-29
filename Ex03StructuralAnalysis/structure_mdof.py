#===============================================================================
'''
Project: Structural Wind Engineering WS20-21
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek

        MDoF system solver using direct time integration - Generalized-Alpha Schemem,
		a monolithic formulation

Author: mate.pentek@tum.de

        Based upon the original implementation for a SDoF system by M. Andre described in:
	    Formulation of the Generalized-Alpha method for LAGRANGE. Technical Report, Chair
         of Structural Analysis @TUM, 2012.

        See J. Chung, G.M. Hulbert: A time integration algorithm for structural dynamics
        with improved numerical dissipation: the generalized-aplha mehod. ASME J. Appl.
        Mech., 60:371-375,1993.

Description: This is a solver for direct numerical time integration for a 2DoF system.
        It assumes linear DOFs with a Generalized alpha scheme with fixed dt.

Created on:  15.03.2016
Last update: 19.11.2020
'''
#===============================================================================

import numpy as np


class StructureMDoF:
    # constructor of the class
    def __init__(self, dt, m, b, k, p_inf, vu0, vv0, va0):
        # introducing and initializing properties and coefficients
        # construct an object self with the input arguments dt, M, B, K,
        # p_inf, u0, v0, a0

        # time step
        self.dt = dt

        # structure mass, damping and spring stiffness
        self.m = m
        self.b = b
        self.k = k

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

        #coefficients for damping
        self.a1b = (1.0 - self.alpha_f) * self.gamma / (self.beta * self.dt)
        self.a2b = (1.0 - self.alpha_f) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alpha_f) * (0.5 * self.gamma / self.beta - 1.0) * self.dt

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

		# structure
        # initial displacement, velocity and acceleration
        self.u0 = np.array(vu0)
        self.v0 = np.array(vv0)
        self.a0 = np.array(va0)

        # initialize for time integration
        self.u1 = np.array(vu0)
        self.v1 = np.array(vv0)
        self.a1 = np.array(va0)

		# force from a previous time step (initial force)
        self.f0 = np.dot(self.m,self.a0) + np.dot(self.b,self.v0) + np.dot(self.k,self.u0)
        self.f1 = np.dot(self.m,self.a1) + np.dot(self.b,self.v1) + np.dot(self.k,self.u1)

    def print_setup(self):
        print("Gen-Alpha time integration setup:")
        print("dt: ", self.dt)
        print("alphaM: ", self.alpha_f)
        print("alphaF: ", self.alpha_m)
        print("gamma: ", self.gamma)
        print("beta: ", self.beta)
        print()
        print("Structural setup:")
        print("mass: ", self.m)
        print("damping: ", self.b)
        print("stiffness: ", self.k)
        print()

    def print_values_at_current_step(self, n):
        print("Printing values at step no: ", n, " (+1)")
        print("u0: ", self.u1)
        print("v0: ", self.v1)
        print("a0: ", self.a1)
        print("f0: ", self.f1)
        print()

    def get_displacement(self):
        return self.u1

    def get_velocity(self):
        return self.v1

    def get_acceleration(self):
        return self.a1

    def get_old_displacement(self):
        return self.u0

    def get_old_velocity(self):
        return self.v0

    def get_old_acceleration(self):
        return self.a0

    def solve_structure(self, f1):
        # sys of eq reads: LHS * u1 = RHS

        F = (1.0 - self.alpha_f) * f1 + self.alpha_f * self.f0

        LHS = self.a1h * self.m + self.a2h * self.b + self.a3h * self.k
        RHS = np.dot(self.m,(self.a1m * self.u0 + self.a2m * self.v0 + self.a3m * self.a0))
        RHS += np.dot(self.b,(self.a1b * self.u0 + self.a2b * self.v0 + self.a3b * self.a0))
        RHS += np.dot(self.a1k * self.k, self.u0) + F

        # update self.f1
        self.f1 = np.array(f1)

        # updates self.u1,v1,a1
        self.u1 = np.linalg.solve(LHS, RHS)
        self.v1 = self.a1v * (self.u1 - self.u0) + self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + self.a2a * self.v0 + self.a3a * self.a0

    def update_structure_timestep(self):
        # update displacement, velocity and acceleration
        self.u0 = np.copy(self.u1)
        self.v0 = np.copy(self.v1)
        self.a0 = np.copy(self.a1)

        # update the force
        self.f0 = np.copy(self.f1)

    def predict_displacement(self):
        return 2.0 * self.u1 - self.u0