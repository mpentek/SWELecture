# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de, catharina.czech@tum.de, peter.kupas@tum.de
    
Note: ...

Created on:  22.11.2017
Last update: 16.11.2018
'''
# ===============================================================================
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import sys

from source.custom_files import RecursivelyValidateAndAssignDefaults

# import visualize_results_utilities to use it as a python object
from source import visualize_result_utilities


class AnalysisType(object):
    """
    Base class for the different analysis types
    """

    def __init__(self, structure_model, name="DefaultAnalysisType"):
        self.name = name

        # the structure model - geometry and physics - has the Dirichlet BC
        # for the bottom node included
        self.structure_model = structure_model

        self.displacement = None
        self.rotation = None

        self.force = None
        self.reaction = None
        self.moment = None

    def solve(self):
        """
        Solve for something
        """
        print("Solving for something in AnalysisType base class \n")
        pass


class StaticAnalysis(AnalysisType):
    """
    Dervied class for the static analysis of a given structure model        
    """

    def __init__(self, structure_model, name="StaticAnalysis"):

        super().__init__(structure_model, name)

    def solve(self, ext_force):
        print("Solving for ext_force in StaticAnalysis derived class \n")
        self.force = ext_force

        # print(self.structure_model.category)

        self.result = np.linalg.solve(self.structure_model.k, self.force)

        if self.structure_model.category in ['SDoF', 'MDoFShear']:
            self.displacement = self.result
            self.reaction = -1/2 * \
                self.structure_model.k[0, 0] * self.displacement[0]
        elif self.structure_model.category in ['MDoFBeam', 'SDoFBridge']:
            self.displacement = self.result[::2]
            self.rotation = self.result[1::2]
            self.force = self.force[::2]
            # TODO this is a temporary fix the force is updated to only joint forces not moments for plotting
            self.reaction = -1/2 * \
                self.structure_model.k[0, 0] * self.displacement[0]
            self.moment = -1/2 * \
                self.structure_model.k[0, 1] * self.rotation[0]
        else:
            sys.exit('Selected structure not implemented!')

        self.structure_model.nodal_coordinates["x"] = np.add(
            self.structure_model.nodal_coordinates["x0"], self.displacement)
        self.structure_model.nodal_coordinates["y"] = self.structure_model.nodal_coordinates["y0"]

    def plot_solve_result(self):
        """
        Pass to plot function:
            from structure model undeformed geometry
            self.displacmenet
            self.force
            self.reaction_force
        """

        print("Plotting result in StaticAnalysis \n")

        origin_point = np.zeros(1)

        geometry = {"undeformed": [np.append(origin_point, self.structure_model.nodal_coordinates["x0"]),
                                   np.append(origin_point, self.structure_model.nodal_coordinates["y0"])],
                    "deformation": [np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["x"]),
                                                np.append(origin_point, self.structure_model.nodal_coordinates["x0"])),
                                    np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["y"]),
                                                np.append(origin_point, self.structure_model.nodal_coordinates["y0"]))],
                    "deformed": None}

        force = {"external": [np.append(origin_point, self.force), np.zeros(len(self.force) + 1)],
                 "base_reaction": [np.append(self.reaction, origin_point), np.zeros(len(self.force) + 1)]}

        scaling = {"deformation": 1,
                   "force": 1}

        plot_title = "Static Analysis : "

        visualize_result_utilities.plot_result(plot_title,
                                               geometry,
                                               force,
                                               scaling,
                                               1)


class EigenvalueAnalysis(AnalysisType):
    """
    Derived class for the (dynamic) eigenvalue analysis of a given structure model        
    """

    def __init__(self, structure_model, name="EigenvalueAnalysis"):
        super().__init__(structure_model, name)

        # adding additional attributes to the derived class
        self.eigenform = None
        self.frequency = None
        self.period = None

    def solve(self):

        eig_values_raw, eig_modes_raw = linalg.eigh(
            self.structure_model.k, self.structure_model.m)
        # rad/s
        eig_values = np.sqrt(np.real(eig_values_raw))
        # 1/s = Hz
        eig_freqs = eig_values / 2. / np.pi
        # s
        eig_pers = 1. / eig_freqs
        # sort eigenfrequencies
        eig_freqs_sorted_indices = np.argsort(eig_freqs)
        ##

        # normalize results
        # http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect19.d/IAST.Lect19.Slides.pdf
        # normalize - unit generalized mass - slide 23

        [rows, columns] = eig_modes_raw.shape

        eig_modes_norm = np.zeros((rows, columns))

        gen_mass_raw = np.zeros(columns)
        gen_mass_norm = np.zeros(columns)

        print("Generalized mass should be identity")
        for i in range(len(eig_values_raw)):
            gen_mass_raw[i] = (np.transpose(eig_modes_raw[:, i])).dot(
                self.structure_model.m).dot(eig_modes_raw[:, i])

            unit_gen_mass_norm_fact = np.sqrt(gen_mass_raw[i])

            eig_modes_norm[:, i] = eig_modes_raw[:, i]/unit_gen_mass_norm_fact

            gen_mass_norm[i] = (np.transpose(eig_modes_norm[:, i])).dot(
                self.structure_model.m).dot(eig_modes_norm[:, i])
            #print("norm ", i, ": ",gen_mass_norm[i])

        #print("Multiplication check: thethaT dot M dot theta: ",(np.transpose(eig_modes_norm)).dot(self.structure_model.m).dot(eig_modes_norm)," numerically 0 for off-diagonal terms")
        print()

        self.eigenform = np.zeros(eig_modes_norm.shape)
        self.frequency = np.zeros(eig_freqs.shape)
        self.period = np.zeros(eig_pers.shape)

        for index in range(len(eig_freqs)):
            self.eigenform[:, index] = eig_modes_norm[:,
                                                      eig_freqs_sorted_indices[index]]
            self.frequency[index] = eig_freqs[eig_freqs_sorted_indices[index]]
            self.period[index] = eig_pers[eig_freqs_sorted_indices[index]]

    def plot_selected_eigenmode(self, selected_mode):
        """
        Pass to plot function:
            from structure model undeformed geometry
            self.eigenform -> as displacement  
            self.frequency -> in legend
            self.period -> in legend

        """
        selected_mode = selected_mode - 1

        print("Plotting result for a selected eigenmode in EigenvalueAnalysis \n")

        if self.structure_model.category in ['SDoF', 'MDoFShear']:
            self.structure_model.nodal_coordinates["x"] = np.add(
                self.structure_model.nodal_coordinates["x0"], self.eigenform[:, selected_mode])

        elif self.structure_model.category in ['MDoFBeam', 'BDoFBridge']:
            self.structure_model.nodal_coordinates["x"] = np.add(
                self.structure_model.nodal_coordinates["x0"], self.eigenform[::2, selected_mode])

        else:
            sys.exit('Wrong structural type selected for eigenmodes')

        self.structure_model.nodal_coordinates["y"] = self.structure_model.nodal_coordinates["y0"]
        origin_point = np.zeros(1)

        geometry = {"undeformed": [np.append(origin_point, self.structure_model.nodal_coordinates["x0"]),
                                   np.append(origin_point, self.structure_model.nodal_coordinates["y0"])],
                    "deformation": [np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["x"]),
                                                np.append(origin_point, self.structure_model.nodal_coordinates["x0"])),
                                    np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["y"]),
                                                np.append(origin_point, self.structure_model.nodal_coordinates["y0"]))],
                    "deformed": None}

        force = {"external": None,
                 "base_reaction": None}

        scaling = {"deformation": 1,
                   "force": 1}

        plot_title = " Eigenmode: " + str(selected_mode+1)
        plot_title += "  Frequency: " + \
            '{0:.2f}'.format(self.frequency[selected_mode])
        plot_title += "  Period: " + \
            '{0:.2f}'.format(self.period[selected_mode])

        visualize_result_utilities.plot_result(plot_title,
                                               geometry,
                                               force,
                                               scaling,
                                               1)

    def plot_selected_first_n_eigenmodes(self, number_of_modes):
        """
        Pass to plot function:
            from structure model undeformed geometry
            self.eigenform -> as displacement  
            self.frequency -> in legend
            self.period -> in legend
        """

        print("Plotting result for selected first n eigenmodes in EigenvalueAnalysis \n")

        if self.structure_model.category in ['SDoF', 'MDoFShear']:
            self.structure_model.nodal_coordinates["x"] = self.eigenform

        elif self.structure_model.category in ['MDoFBeam', 'BDoFBridge']:
            self.structure_model.nodal_coordinates["x"] = self.eigenform[::2]

        else:
            sys.exit("Selected structure not implemented!")

        self.structure_model.nodal_coordinates["y"] = self.structure_model.nodal_coordinates["y0"]

        origin_point = np.zeros(1)
        origin_vector = np.zeros(len(self.eigenform))

        geometry = {"undeformed": [np.append(origin_point, self.structure_model.nodal_coordinates["x0"]),
                                   np.append(origin_point, self.structure_model.nodal_coordinates["y0"])],
                    "deformation": [np.vstack((origin_vector, self.structure_model.nodal_coordinates["x"])),
                                    np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["y"]),
                                                np.append(origin_point, self.structure_model.nodal_coordinates["y0"]))],
                    "deformed": None}

        force = {"external": None,
                 "base_reaction": None}

        scaling = {"deformation": 1,
                   "force": 1}
        #print("Geometry: ", geometry)
        #print("Self.Nodal coordinates: ", self.structure_model.nodal_coordinates["x"])

        plot_title = " "
        for selected_mode in range(number_of_modes):
            plot_title += "Eigenmode " + str(selected_mode + 1) + "  Frequency: " + str(np.round(
                self.frequency[selected_mode], 3)) + "  Period: " + str(np.round(self.period[selected_mode], 3)) + "\n"
        visualize_result_utilities.plot_result(plot_title,
                                               geometry,
                                               force,
                                               scaling,
                                               number_of_modes)

    def animate_selected_eigenmode(self, selected_mode):
        """
        Pass to plot function:
            from structure model undeformed geometry
            self.eigenform -> as displacement  
            self.frequency -> in legend
            self.period -> in legend  
        """
        skip = max(1, int(5-0.3*selected_mode))
        print(skip)
        selected_mode = selected_mode - 1

        print("Animating eigenmode in EigenvalueAnalysis \n")

        if self.structure_model.category in ['SDoF', 'MDoFShear']:
            displacement = self.eigenform[:, selected_mode]

        elif self.structure_model.category in ['MDoFBeam', 'MdoFBridge']:
            displacement = self.eigenform[::2, selected_mode]

        else:
            sys.exit("Selected structure not implemented!")

        time_steps = 100
        array_time = np.linspace(0, self.period[selected_mode], time_steps)

        displacement_time_history = [[] for i in range(len(array_time))]
        for i in range(len(array_time)):
            displacement_time_history[i] = [value * np.sin(
                2 * np.pi * self.frequency[selected_mode] * array_time[i]) for value in displacement]

        self.structure_model.nodal_coordinates["x"] = [
            [] for i in range(len(array_time))]
        self.structure_model.nodal_coordinates["y"] = [
            [] for i in range(len(array_time))]
        for i in range(len(array_time)):
            self.structure_model.nodal_coordinates["x"][i] = displacement_time_history[i]
            self.structure_model.nodal_coordinates["y"][i] = self.structure_model.nodal_coordinates["y0"]

        origin_point = np.zeros(1)

        geometry = {"undeformed": [np.append(origin_point, self.structure_model.nodal_coordinates["x0"]),
                                   np.append(origin_point, self.structure_model.nodal_coordinates["y0"])],
                    "deformation": [[[] for i in range(len(array_time))], [[] for i in range(len(array_time))]],
                    "deformed": [[[] for i in range(len(array_time))], [[] for i in range(len(array_time))]]}

        for i in range(len(array_time)):

            geometry["deformation"][0][i] = np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["x"][i]),
                                                        np.append(origin_point, self.structure_model.nodal_coordinates["x0"]))
            geometry["deformation"][1][i] = np.subtract(np.append(origin_point, self.structure_model.nodal_coordinates["y"][i]),
                                                        np.append(origin_point, self.structure_model.nodal_coordinates["y0"]))

        force = {"external": None,
                 "base_reaction": None}

        scaling = {"deformation": 1,
                   "force": 1}

        plot_title = "Eigenmode: " + str(selected_mode + 1)
        plot_title += "  Frequency: " + \
            '{0:.2f}'.format(self.frequency[selected_mode])
        plot_title += "  Period: " + \
            '{0:.2f}'.format(self.period[selected_mode])

        visualize_result_utilities.animate_result(plot_title,
                                                  array_time,
                                                  geometry,
                                                  force,
                                                  scaling)


class DynamicAnalysis(AnalysisType):
    """
    Dervied class for the dynamic analysis of a given structure model        

    """

    def __init__(self, structure_solver, force, time_settings, name="DynamicAnalysis"):

        super().__init__(structure_solver, name)

        default_settings = {
            "start_time": 0.0,
            "end_time": 5.0,
            "echo_level": 0,
            "time_step": 0.1,
            # "print_colors" : true
        }

        RecursivelyValidateAndAssignDefaults(
            default_settings, time_settings)   # check how to do this
        print("Force: ", len(force))
        # overwriting attribute from base constructors
        self.force = force
        self.dt = time_settings["time_step"]
        self.start_time = time_settings["start_time"]
        self.end_time = time_settings["end_time"]
        self.array_time = np.arange(
            self.start_time, self.end_time + self.dt, self.dt)
        rows = len(self.force)
        cols = len(self.array_time)
        self.structure_solver = structure_solver

        # adding additional attributes to the derived class
        self.displacement = np.zeros((rows, cols))
        self.velocity = np.zeros((rows, cols))
        self.acceleration = np.zeros((rows, cols))

        if self.force.shape[1] != len(self.array_time):
            err_msg = "The time step for forces does not match the time step defined"
            raise Exception(err_msg)

    def solve(self):
        print("Solving the structure for dynamic loads \n")
        # time loop
        self.structure_solver.Initialize()
        for i in range(1, len(self.array_time)):
            current_time = self.array_time[i]
            self.structure_solver.ApplyExternalLoad(self.force[:, i])
            self.structure_solver.SolveSolutionStep()

            # appending results to the list
            self.displacement[:, i] = self.structure_solver.GetSolutionStepValue(
                "DISPLACEMENT")
            self.velocity[:, i] = self.structure_solver.GetSolutionStepValue(
                "VELOCITY")
            self.acceleration[:, i] = self.structure_solver.GetSolutionStepValue(
                "ACCELERATION")
            # update results
            self.structure_solver.AdvanceInTime(current_time)

    def plot_result_at_dof(self, selected_dof, selected_result):
        """
        Pass to plot function:
            Plots the time series of required quantitiy 
        """
        print('Plotting result for selected dof in dynamic analysis \n')
        plot_title = selected_result.capitalize() + ' at DoF ' + str(selected_dof)

        # the angular results are not plotted for MDoFBeam model 
        if self.structure_solver.model.category in ['SDoF', 'MDoFShear']:
            dof = selected_dof - 1
        elif self.structure_solver.model.category in ['MDoFBeam', 'MDoFBridge']:
            dof = 2 * (selected_dof - 1)
        else:
            sys.exit()
    
        if selected_result == 'displacement':
            result_data = self.displacement[dof, :]
        elif selected_result == 'velocity': 
            result_data = self.velocity[dof, :]
        elif selected_result == 'acceleration':
            result_data = self.acceleration[dof, :]
        else: 
            sys.exit()

        visualize_result_utilities.plot_dynamic_result(plot_title, result_data, self.array_time)

    def plot_selected_time_step(self, selected_time_step):
        """
        Pass to plot function:
            from structure model undeformed geometry
            self.displacement -> here as time series -> select closes results to a requested time_step [s]  

        """
        print("Plotting result for a selected time step in DynamicAnalysis \n")

        # find closesed time step
        idx = (np.abs(self.array_time-selected_time_step)).argmin()

        # TODO for the dynamic analysis create self.structure_model.nodal_coordinates after solve not here
        if self.structure_solver.model.category in ['SDoF', 'MDoFShear']:
            self.structure_solver.model.nodal_coordinates["x"] = self.displacement[:, idx]
        elif self.structure_solver.model.category in ['MDoFBeam', 'MDoFBridge']:
            self.structure_solver.model.nodal_coordinates["x"] = self.displacement[::2, idx]
        else:
            sys.exit()

        self.structure_solver.model.nodal_coordinates["y"] = self.structure_solver.model.nodal_coordinates["y0"]

        origin_point = np.zeros(1)

        geometry = {"undeformed": [np.append(origin_point, self.structure_solver.model.nodal_coordinates["x0"]),
                                   np.append(origin_point, self.structure_solver.model.nodal_coordinates["y0"])],
                    "deformation": [np.subtract(np.append(origin_point, self.structure_solver.model.nodal_coordinates["x"]),
                                                np.append(origin_point, self.structure_solver.model.nodal_coordinates["x0"])),
                                    np.subtract(np.append(origin_point, self.structure_solver.model.nodal_coordinates["y"]),
                                                np.append(origin_point, self.structure_solver.model.nodal_coordinates["y0"]))],
                    "deformed": None}

        force = {"external": [np.append(origin_point, self.force), np.zeros(len(self.force) + 1)],
                 "base_reaction": [np.append(self.reaction, origin_point), np.zeros(len(self.force) + 1)]}

        scaling = {"deformation": 1,
                   "force": 1}

        plot_title = "Dyanimc Analyis: Deformation at t = " + \
            str(self.array_time[idx]) + " [s]"

        visualize_result_utilities.plot_result(plot_title,
                                               geometry,
                                               force,
                                               scaling,
                                               1)

    def animate_time_history(self):
        """
        Pass to plot function:
            from structure model undeformed geometry
            self.displacement -> here as time series  
        """

        print("Animating time history in DynamicAnalysis \n")

        self.structure_solver.model.nodal_coordinates["x"] = [
            [] for i in range(len(self.array_time))]
        self.structure_solver.model.nodal_coordinates["y"] = [
            [] for i in range(len(self.array_time))]
        for i in range(len(self.array_time)):
            if self.structure_solver.model.category in ['SDoF', 'MDoFShear']:
                self.structure_solver.model.nodal_coordinates["x"][i] = self.displacement[:, i]
            elif self.structure_solver.model.category in ['MDoFBeam', 'MDoFBridge']:
                self.structure_solver.model.nodal_coordinates["x"][i] = self.displacement[::2, i]
            else:
                sys.exit()

            self.structure_solver.model.nodal_coordinates["y"][
                i] = self.structure_solver.model.nodal_coordinates["y0"]

        origin_point = np.zeros(1)

        geometry = {"undeformed": [np.append(origin_point, self.structure_solver.model.nodal_coordinates["x0"]),
                                   np.append(origin_point, self.structure_solver.model.nodal_coordinates["y0"])],
                    "deformation": [[[] for i in range(len(self.array_time))], [[] for i in range(len(self.array_time))]],
                    "deformed": [[[] for i in range(len(self.array_time))], [[] for i in range(len(self.array_time))]]}

        for i in range(len(self.array_time)):
            geometry["deformation"][0][i] = np.subtract(np.append(origin_point, self.structure_solver.model.nodal_coordinates["x"][i]),
                                                        np.append(origin_point, self.structure_solver.model.nodal_coordinates["x0"]))
            geometry["deformation"][1][i] = np.subtract(np.append(origin_point, self.structure_solver.model.nodal_coordinates["y"][i]),
                                                        np.append(origin_point, self.structure_solver.model.nodal_coordinates["y0"]))

        force = {"external": None,
                 "base_reaction": None}

        scaling = {"deformation": 1,
                   "force": 1}

        plot_title = "Dyanimc Analyis: Deformation over time"

        visualize_result_utilities.animate_result(plot_title,
                                                  self.array_time,
                                                  geometry,
                                                  force,
                                                  scaling)
