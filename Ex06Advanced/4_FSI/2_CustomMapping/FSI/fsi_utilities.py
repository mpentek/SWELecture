# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS18-19 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
        FSI utilities

Author: philipp.bucher@tum.de, a.winterstein@tum.de, mate.pentek@tum.de, anoop.kodakkal@tum.de
      
Note: ...

Created on:  16.01.2018
Last update: 06.12.2019
'''
# ===============================================================================


import KratosMultiphysics
#import KratosMultiphysics.StructuralMechanicsApplication as KratosStructuralMechanics

import sys
import time
import math


'''
Additionals functionalities to make the FSI possible
This would otherwise be done by the EmpireApplication
or the FSIApplication
'''


def Norm(array):
    '''
    Return the L2 norm
    '''
    norm = 0
    for row in array:
        try:    # array is a matrix
            for entry in row:
                norm += entry**2
        except:  # array is a vector
            norm += row**2
    return pow(norm, 0.5)


def GetNodeid(structural_model_part, point):
    '''
    Gets the node id of the given coordinate 
    '''
    tol = 1e-6
    found_id = KratosMultiphysics.BruteForcePointLocator(
        structural_model_part).FindNode(point, tol)
    if found_id > -1:
        return found_id
    else:
        raise Exception("No Node found on the requested point")


def GetDisplacements(nodes_of_structure, dimension=3):
    '''
    Gets the (nodal) displacements of the nodes on the 
    submodel part of the structure on the interface with the fluid
    '''
    displacements = [0.0] * dimension * len(nodes_of_structure)
    index = 0
    for node in nodes_of_structure:
        displacements[dimension*index] = node.GetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0)
        displacements[dimension*index +
                      1] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0)
        if dimension == 3:
            displacements[dimension*index+2] = node.GetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Z, 0)
        index += 1

    return displacements


def SetDisplacements(displacements, nodes_of_structure, dimension=3):
    '''
    Sets the (nodal) displacements of the nodes on the 
    submodel part of the structure on the interface with the fluid
    '''
    index = 0
    for node in nodes_of_structure:
        node.SetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_X, 0, displacements[dimension*index])
        node.SetSolutionStepValue(
            KratosMultiphysics.DISPLACEMENT_Y, 0, displacements[dimension*index+1])
        if dimension == 3:
            node.SetSolutionStepValue(
                KratosMultiphysics.DISPLACEMENT_Z, 0, displacements[dimension*index+2])

        index += 1


def NeumannToStructure(mapper, flag):
    #mapper.InverseMap(KratosStructuralMechanics.POINT_LOAD, KratosMultiphysics.REACTION, flag)
    mapper.InverseMap(mapper.inverse_map["origin_variable"].GetString(
    ), mapper.inverse_map["destination_variable"].GetString(), flag)


def DisplacementToMesh(mapper):
    #mapper.Map(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.MESH_DISPLACEMENT)
    mapper.Map(mapper.map["origin_variable"].GetString(),
               mapper.map["destination_variable"].GetString())


'''
Additional functionalities to create here a custom convergence accelerator
This would otherwise be done by the EmpireApplication
'''


def CreateConvergenceAccelerator(convergence_accelerator_settings):
    if convergence_accelerator_settings["type"].GetString() == "aitken":
        return AitkenConvergenceAccelerator(convergence_accelerator_settings)
    else:
        raise Exception("the requested convergence accelerator type " +
                        convergence_accelerator_settings["type"].GetString() + "is not implemented.")


class ConvergenceAcceleratorBase(object):
    def __init__(self, convergence_accelerator_settings):
        self.max_iter = convergence_accelerator_settings["max_iterations"].GetInt(
        )
        self.res_rel_tol = convergence_accelerator_settings["residual_relative_tolerance"].GetDouble(
        )
        self.res_abs_tol = convergence_accelerator_settings["residual_absolute_tolerance"].GetDouble(
        )
        self.rel_coef_initial = convergence_accelerator_settings["relaxation_coefficient_initial_value"].GetDouble(
        )

    def CalculateResidual(self, solution, old_solution):
        '''
        Calculates the residual based upon the (current) solution
        and the old one
        '''
        residual = [0.0] * len(solution)
        for index in range(len(solution)):
            residual[index] = solution[index] - old_solution[index]
        return residual

    def CalculateRelaxedSolution(self, relaxation_coefficient, old_solution, residual):
        '''
        Calculates the relaxed (i.e. new) solution
        '''
        for index in range(len(old_solution)):
            old_solution[index] = old_solution[index] + \
                relaxation_coefficient * residual[index]
        return old_solution  # this is the relaxed new solution

    def ComputeRelaxationCoefficient(self):
        print("Function needs to be implemented in the derived class")


class AitkenConvergenceAccelerator(ConvergenceAcceleratorBase):
    def ComputeRelaxationCoefficient(self, old_coefficient, residual, old_residual, iteration, max_initial_coefficient=0.2, upper_limit=2.5, lower_limit=-0.5):
        if iteration < 1:
            new_coefficient = min(old_coefficient, max_initial_coefficient)
        else:
            numerator = 0
            denominator = 0
            for i in range(len(residual)):
                numerator += old_residual[i] * (residual[i] - old_residual[i])
                denominator += pow(residual[i]-old_residual[i], 2)
            new_coefficient = - old_coefficient * (numerator/denominator)

            # force some limit on the calculated coefficient
            if new_coefficient > upper_limit:
                new_coefficient = upper_limit
                print("WARNING: upper limit of " + str(upper_limit) +
                      "reached in Aitken: ComputeCustomRelaxation()")
            elif new_coefficient < lower_limit:
                new_coefficient = lower_limit
                print("WARNING: lower limit of " + str(lower_limit) +
                      "reached in Aitken: ComputeCustomRelaxation()")
        return new_coefficient


'''
Additional functionalities to create here a custom mapper
This would otherwise be done by the MappingApplication
'''


def CreateMapper(destination_model_part, origin_model_part, mapper_settings):
    return CustomMapper(destination_model_part, origin_model_part, mapper_settings)


class CustomMapper(object):
    def __init__(self, destination_model_part, origin_model_part, mapper_settings):
        # here DESTINATION = interface nodes on FLUID
        destination_interface_name = mapper_settings["interface_submodel_part_destination"].GetString(
        )
        self.destination_interface = destination_model_part.GetSubModelPart(
            destination_interface_name)

        # here ORIGIN = interface nodes on STRUCTURE
        origin_interface_name = mapper_settings["interface_submodel_part_origin"].GetString(
        )
        self.origin_interface = origin_model_part.GetSubModelPart(
            origin_interface_name)

        self.coupling_matrix = self._GetCouplingMatrix()

        self.map = mapper_settings["map"]
        self.inverse_map = mapper_settings["inverse_map"]

    def _GetCouplingMatrix(self):

        coupling_matrix = {}  # this is a map (key-value pair), not a list!

        for destination_node in self.destination_interface.Nodes:
            d_x = destination_node.X
            d_y = destination_node.Y
            id_1 = 0
            id_2 = 0

            # introduce a very/unlikely high initial value
            dist_min = 1e12
            dist_min2 = 1e12

            for origin_node in self.origin_interface.Nodes:
                o_x = origin_node.X
                o_y = origin_node.Y

                distance = ((d_x - o_x)**2 + (d_y - o_y)**2)**0.5

                # check if origin_node is closer to destination_node than previous origin_nodes
                if distance < dist_min:  # origin_node is closer than all other origin_nodes
                    id_2 = id_1
                    dist_min2 = dist_min

                    id_1 = origin_node.Id
                    dist_min = distance

                elif distance < dist_min2 and distance > dist_min:  # origin_node is closer than all second origin_node
                    id_2 = origin_node.Id
                    dist_min2 = distance

            # compute weights to struct nodes
            dist_a = dist_min
            dist_b = dist_min2
            dist_abs = dist_a + dist_b

            w_a = dist_b/dist_abs
            w_b = dist_a/dist_abs

            # node ID is the key, struct node IDs and weights are value
            coupling_matrix[destination_node.Id] = [id_1, id_2, w_a, w_b]

            print("NONMATCHING " + "(fluid) destination node", destination_node.Id,
                  "paired with (structure) origin node", id_1, "& ", id_2)

            print(
                "======================================================================")
            print(
                "||||||||||||||||||||||||INTERFACE PAIRING DONE||||||||||||||||||||||||")
            print(
                "======================================================================")

        return coupling_matrix

    # def InverseMap(self, origin_variable, destination_variable, swap_sign):
    def InverseMap(self, origin_var_name, destination_var_name, swap_sign):
        '''
        Used as a conservative (inverse)map - adequate for forces, reactions, etc. -
        to conserve the values in an integral sense.
        '''

        # for this case this will be: KratosStructuralMechanics.POINT_LOAD
        #origin_var_name = origin_variable.Name()
        origin_var_comp_x = KratosMultiphysics.KratosGlobals.GetVariable(
            origin_var_name + "_X")
        origin_var_comp_y = KratosMultiphysics.KratosGlobals.GetVariable(
            origin_var_name + "_Y")

        # for this case this will be: KratosMultiphysics.REACTION
        # here one has to change the reaction and take it with a negative sign!!!
        #destination_var_name = destination_variable.Name()
        destination_var_comp_x = KratosMultiphysics.KratosGlobals.GetVariable(
            destination_var_name + "_X")
        destination_var_comp_y = KratosMultiphysics.KratosGlobals.GetVariable(
            destination_var_name + "_Y")

        multiplicator = 1.0
        if swap_sign:
            multiplicator = -1.0

        for origin_node in self.origin_interface.Nodes:
            # set load to 0 to accumulate it in next step
            origin_node.SetSolutionStepValue(origin_var_comp_x, 0, 0)
            origin_node.SetSolutionStepValue(origin_var_comp_y, 0, 0)

        for destination_node in self.destination_interface.Nodes:
            # get coupling information of destination_node
            coupling_for_destination_node = self.coupling_matrix[destination_node.Id]
            node_1_id = coupling_for_destination_node[0]
            node_2_id = coupling_for_destination_node[1]
            node_1_w = coupling_for_destination_node[2]
            node_2_w = coupling_for_destination_node[3]

            # get component from destination and swap sign if necessary
            val_comp_x = multiplicator * \
                destination_node.GetSolutionStepValue(
                    destination_var_comp_x, 0)
            val_comp_y = multiplicator * \
                destination_node.GetSolutionStepValue(
                    destination_var_comp_y, 0)

            # calculate respective part component part for origin
            origin_node_1_val_comp_x = val_comp_x * node_1_w
            origin_node_1_val_comp_y = val_comp_y * node_1_w

            origin_node_2_val_comp_x = val_comp_x * node_2_w
            origin_node_2_val_comp_y = val_comp_y * node_2_w

            # accumulate component for origin from existing on origin
            origin_node_1_val_comp_x += self.origin_interface.Nodes[node_1_id].GetSolutionStepValue(
                origin_var_comp_x, 0)
            origin_node_1_val_comp_y += self.origin_interface.Nodes[node_1_id].GetSolutionStepValue(
                origin_var_comp_y, 0)

            origin_node_2_val_comp_x += self.origin_interface.Nodes[node_2_id].GetSolutionStepValue(
                origin_var_comp_x, 0)
            origin_node_2_val_comp_y += self.origin_interface.Nodes[node_2_id].GetSolutionStepValue(
                origin_var_comp_y, 0)

            # set value to origin
            self.origin_interface.Nodes[node_1_id].SetSolutionStepValue(
                origin_var_comp_x, 0, origin_node_1_val_comp_x)
            self.origin_interface.Nodes[node_1_id].SetSolutionStepValue(
                origin_var_comp_y, 0, origin_node_1_val_comp_y)

            self.origin_interface.Nodes[node_2_id].SetSolutionStepValue(
                origin_var_comp_x, 0, origin_node_2_val_comp_x)
            self.origin_interface.Nodes[node_2_id].SetSolutionStepValue(
                origin_var_comp_y, 0, origin_node_2_val_comp_y)

    # def Map(self, origin_variable, destination_variable):
    def Map(self, origin_var_name, destination_var_name):
        '''
        Used as a consistent  map - adequate for example for displacements (kinematics).
        '''

        # for this case this will be: KratosMultiphysics.DISPLACEMENT
        #origin_var_name = origin_variable.Name()
        origin_var_comp_x = KratosMultiphysics.KratosGlobals.GetVariable(
            origin_var_name + "_X")
        origin_var_comp_y = KratosMultiphysics.KratosGlobals.GetVariable(
            origin_var_name + "_Y")

        # for this case this will be: KratosMultiphysics.MESH_DISPLACEMENT
        #destination_var_name = destination_variable.Name()
        destination_var_comp_x = KratosMultiphysics.KratosGlobals.GetVariable(
            destination_var_name + "_X")
        destination_var_comp_y = KratosMultiphysics.KratosGlobals.GetVariable(
            destination_var_name + "_Y")

        for destination_node in self.destination_interface.Nodes:
            # get coupling information of destination_node
            coupling_for_destination_node = self.coupling_matrix[destination_node.Id]
            node_1_id = coupling_for_destination_node[0]
            node_2_id = coupling_for_destination_node[1]
            node_1_w = coupling_for_destination_node[2]
            node_2_w = coupling_for_destination_node[3]

            # get components from origin
            origin_node_1_val_comp_x = self.origin_interface.Nodes[node_1_id].GetSolutionStepValue(
                origin_var_comp_x, 0)
            origin_node_1_val_comp_y = self.origin_interface.Nodes[node_1_id].GetSolutionStepValue(
                origin_var_comp_y, 0)

            origin_node_2_val_comp_x = self.origin_interface.Nodes[node_2_id].GetSolutionStepValue(
                origin_var_comp_x, 0)
            origin_node_2_val_comp_y = self.origin_interface.Nodes[node_2_id].GetSolutionStepValue(
                origin_var_comp_y, 0)

            # linear interpolation of value for components
            destination_node_val_comp_x = origin_node_1_val_comp_x * \
                node_1_w + origin_node_2_val_comp_x * node_2_w
            destination_node_val_comp_y = origin_node_1_val_comp_y * \
                node_1_w + origin_node_2_val_comp_y * node_2_w

            # set value to destination
            destination_node.SetSolutionStepValue(
                destination_var_comp_x, 0, destination_node_val_comp_x)
            destination_node.SetSolutionStepValue(
                destination_var_comp_y, 0, destination_node_val_comp_y)
