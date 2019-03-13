# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek

        Kratos Multiphysics - current developments of the EmpireApplication
        for multiphysics coupling
        
        https://github.com/KratosMultiphysics/Kratos

Last update: 16.11.2018
'''
# ===============================================================================
# makes these scripts backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Importing the base class
from time_integration_base_scheme import TimeIntegrationBaseScheme

# Importing tools
from source.custom_files import ValidateAndAssignDefaults

# Other imports
import numpy as np
import json
import os


def CreateScheme(scheme_settings):
    return TimeIntegrationBdf2Scheme(scheme_settings)



class TimeIntegrationBdf2Scheme(TimeIntegrationBaseScheme):
    """
    A single-degree-of-freedom SDoF model

    Using for testing of the MDoF solver
    """

    def __init__(self, scheme_settings):

        default_settings = {
            "type": "bdf2",
            "settings": {}
        }

        # add buffer size - this is not user-specified
        # each derived scheme specifies it
        scheme_settings["settings"].update({"buffer_size": 5})

        ValidateAndAssignDefaults(default_settings, scheme_settings)

        # base scheme settings
        super(TimeIntegrationBdf2Scheme, self).__init__(
            scheme_settings["settings"])

    def _AssembleLHS(self, model):
        """
        """
        return (9 * model.m + 6 * model.b * self.dt + 4 * model.k * self.dt**2)

    def _AssembleRHS(self, model):
        """
        """
        self.buffer[3, 0, :] = self.force
        
        RHS = 4 * self.buffer[3, 0, :] * self.dt**2
        RHS -= np.dot(24 * model.m + 8 * model.b * self.dt, self.buffer[0, 1, :])
        RHS -= np.dot(-22 * model.m - 2 * model.b * self.dt, self.buffer[0, 2, :])
        RHS -= np.dot(8 * model.m, self.buffer[0, 3, :])
        RHS -= np.dot(model.m, self.buffer[0, 4, :])
        return RHS

        return RHS
    
    def Initialize(self, model):
        """
        """
        # call function from base
        super(TimeIntegrationBdf2Scheme, self).Initialize(model)
        # overwrite with scheme-specific values

        self.buffer[0,4,:] = model.u0  
        self.buffer[0,3,:] = model.u0
        # Euler backward integration scheme is used for the first time steps 

        LHS = (model.m + model.b * self.dt + model.k * self.dt**2)

        RHS2 = self.buffer[3, 2, :] * self.dt**2
        RHS2 -= np.dot(- 2 * model.m - model.b * self.dt, self.buffer[0, 3, :])
        RHS2 -= np.dot(model.m, self.buffer[0, 4, :])
        self.buffer[0, 2, :] = np.linalg.solve(LHS, RHS2)
        
        RHS1 = self.buffer[3, 1, :] * self.dt**2
        RHS1 -= np.dot(- 2 * model.m - model.b * self.dt, self.buffer[0, 2, :])
        RHS1 -= np.dot(model.m, self.buffer[0, 3, :])
        self.buffer[0, 1, :] = np.linalg.solve(LHS, RHS1)

    def UpdateDerivedValues(self):
        """
        """
        c0 =  3 * 0.5 / self.dt
        c1 =  -4 * 0.5 / self.dt
        c2 =  1 * 0.5 / self.dt
        
        self.buffer[1, 0, :] = (c0 * self.buffer[0, 0, :] + 
                                c1 * self.buffer[0, 1, :] +
                                c2 * self.buffer[0, 2, :]) 
        self.buffer[2, 0, :] = (c0 * self.buffer[1, 0, :] + 
                                c1 * self.buffer[1, 1, :] +
                                c2 * self.buffer[1, 2, :]) 

