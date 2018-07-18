#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS17-18
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek
        
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
Last update: 02.11.2017
'''
#===============================================================================
## StructureMDOF class for a MultiDegreeOfFreedom dynamic system

import numpy as np

class StructureMDoF:
    # constructor of the class
    def __init__(self, dt, mM, mB, mK, pInf, vu0, vv0, va0):     
        # introducing and initializing properties and coefficients
        # construct an object self with the input arguments dt, M, B, K,
        # pInf, u0, v0, a0
        
        # time step
        self.dt = dt
        
        # mass, damping and spring stiffness
        self.M = mM
        self.B = mB
        self.K = mK 
			             
        # generalized alpha parameters (to ensure unconditional stability, 2nd order accuracy)
        self.alphaM = (2.0 * pInf - 1.0) / (pInf + 1.0)
        self.alphaF = pInf / (pInf + 1.0)
        self.beta = 0.25 * (1 - self.alphaM + self.alphaF)**2
        self.gamma = 0.5 - self.alphaM + self.alphaF

        # coefficients for LHS
        self.a1h = (1.0 - self.alphaM) / (self.beta * self.dt**2)
        self.a2h = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a3h = 1.0 - self.alphaF

        # coefficients for mass
        self.a1m = self.a1h
        self.a2m = self.a1h * self.dt
        self.a3m = (1.0 - self.alphaM - 2.0 * self.beta) / (2.0 * self.beta)
        
        #coefficients for damping
        self.a1b = (1.0 - self.alphaF) * self.gamma / (self.beta * self.dt)
        self.a2b = (1.0 - self.alphaF) * self.gamma / self.beta - 1.0
        self.a3b = (1.0 - self.alphaF) * (0.5 * self.gamma / self.beta - 1.0) * self.dt
        
        # coefficient for stiffness
        self.a1k = -1.0 * self.alphaF
        
        # coefficients for velocity update
        self.a1v = self.gamma / (self.beta * self.dt)
        self.a2v = 1.0 - self.gamma / self.beta
        self.a3v = (1.0 - self.gamma / (2 * self.beta)) * self.dt
        
        # coefficients for acceleration update
        self.a1a = self.a1v / (self.dt * self.gamma)
        self.a2a = -1.0 / (self.beta * self.dt)
        self.a3a = 1.0 - 1.0 / (2.0 * self.beta)   
        
		#structure
        # initial displacement, velocity and acceleration
        self.u0 = vu0
        self.v0 = vv0
        self.a0 = va0
		        
        # initial displacement, velocity and acceleration
        self.u1 = self.u0
        self.v1 = self.v0
        self.a1 = self.a0
			  
		# force from a previous time step (initial force)
        self.f0 = np.dot(self.M,self.a0) + np.dot(self.B,self.v0) + np.dot(self.K,self.u0)
        self.f1 = np.dot(self.M,self.a1) + np.dot(self.B,self.v1) + np.dot(self.K,self.u1)
		
    def printSetup(self):
        print("Printing Generalized Alpha Method integration scheme setup for 2DOFMonolithic:")
        print("dt: ", self.dt)
        print("alphaM: ", self.alphaF)
        print("alphaF: ", self.alphaM)
        print("gamma: ", self.gamma)
        print("beta: ", self.beta)
        print(" ")
        print("Printing structural setup:")
        print("mass: ", self.M)
        print("damping: ", self.B)
        print("stiffness: ", self.K)
        print(" ") 
 		
        
    def printValuesAtCurrentStep(self, n):
        print("Printing values at step no: ", n, " (+1)")
        print("For structure: ")
        print("u0: ", self.su1)
        print("v0: ", self.sv1)
        print("a0: ", self.sa1)
        print("f0: ", self.f1)
        print(" ")
        print("For damper: ")
        print("u0: ", self.du1)
        print("v0: ", self.dv1)
        print("a0: ", self.da1)
        print(" ")
		
    def getDisplacement(self):
        return self.u1

    def getVelocity(self):
        return self.v1
    
    def getAcceleration(self):
        return self.a1

    def getOldDisplacement(self):
        return self.u0

    def getOldVelocity(self):
        return self.v0
    
    def getOldAcceleration(self):
        return self.a0
		
    def solveStructure(self, f1):
        # sys of eq reads: LHS * u1 = RHS
		
        F = (1.0 - self.alphaF) * f1 + self.alphaF * self.f0
		
        LHS = self.a1h * self.M + self.a2h * self.B + self.a3h * self.K 
        RHS = np.dot(self.M,(self.a1m * self.u0 + self.a2m * self.v0 + self.a3m * self.a0))
        RHS += np.dot(self.B,(self.a1b * self.u0 + self.a2b * self.v0 + self.a3b * self.a0)) 
        RHS += np.dot(self.a1k * self.K, self.u0) + F

        # update self.f1
        self.f1 = f1
        
        # updates self.u1,v1,a1
        self.u1 = np.linalg.solve(LHS, RHS)
        self.v1 = self.a1v * (self.u1 - self.u0) + self.a2v * self.v0 + self.a3v * self.a0
        self.a1 = self.a1a * (self.u1 - self.u0) + self.a2a * self.v0 + self.a3a * self.a0
		
    def updateStructureTimeStep(self):
        # update displacement, velocity and acceleration 
        self.u0 = self.u1
        self.v0 = self.v1
        self.a0 = self.a1
        
        # update the force   
        self.f0 = self.f1

    def predictDisplacement(self):
        return 2.0 * self.u1 - self.u0		
# ===============================================================================

