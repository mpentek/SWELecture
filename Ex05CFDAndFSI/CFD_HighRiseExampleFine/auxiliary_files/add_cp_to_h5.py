#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: mate.pentek@tum.de 

Description: Helper script to postprocess H5 pressure output using additional ref pressure file

Created on:  27.01.2020
Last update: 27.01.2020
'''
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

################################################################
# user input

ref_point_file = os.path.join('ascii_output','reference_point_output.dat')
ref_p_first_t = 0.02
ref_p_delta_t = 0.02
rho = 1.2

deck_results_folder = os.path.join('hdf5_output','deck')
file_prefix = 'NoSlip3D_structure_T-'
h5_first_t = 0.00
h5_last_t = 0.26
h5_delta_t = 0.02
file_idxs = np.linspace(h5_first_t, h5_last_t, int(
    (h5_last_t+h5_delta_t-h5_first_t)/h5_delta_t), endpoint=True)

ramp_up_time = 20
ramp_up_fctr = 1.5
discarded_time = 0.1 #ramp_up_fctr * ramp_up_time

variable_labels = ['time', 'pres']

################################################################
# custom function definition


def get_discard_index(times_series, discarded_time):
    return np.where(times_series >= discarded_time)[0][0]


def get_id_and_position_from_header(file, precision=4):

    with open(file, 'r') as f:
        f.seek(0)
        first_line = f.readline()

    first_line = first_line.split()

    p_id = int(first_line[7])

    p_position = np.asarray([
        np.around(float(first_line[-5][:-1]), precision),
        np.around(float(first_line[-3][:-1]), precision),
        np.around(float(first_line[-1]), precision)])

    return p_id, p_position


def get_point_data_from_file(file_name, variable_labels):
    p_data = {}
    p_data['p_id'], p_data['p_position'] = get_id_and_position_from_header(
        file_name)
    for counter, vl in enumerate(variable_labels):
        p_data[vl] = np.loadtxt(file_name, usecols=(counter, ))

    return p_data


def get_mean(series):
    return np.mean(series)

################################################################
# read in raw date of time series


# from reference point
variable_labels.append('velx')
ref_point = get_point_data_from_file(
    ref_point_file, variable_labels)

################################################################
# postprocess raw data

discard_idx = get_discard_index(ref_point['time'], discarded_time)
print('Ref point raw data from time series discarded until ' +
      str(ref_point['time'][discard_idx]) + ' [s]')

ref_velx = get_mean(ref_point['velx'][discard_idx:])
print('Reference velocity of ' + str(ref_velx) +
      ' [m/s] at height ' + str(ref_point['p_position'][2]) + ' [m]')
fctr = 1 / (0.5 * rho * ref_velx**2)

step_ratio = int(h5_delta_t / ref_p_delta_t)

# from pressure on surface
for counter, idx in enumerate(file_idxs):
    file_name = os.path.join(deck_results_folder, file_prefix + "{0:.2f}".format(np.round(idx, 2)) + '.h5')
    print(file_name)

    with h5py.File(file_name, 'a') as h5f:
        pres_vals = np.array(
            h5f['ResultsData']['NodalSolutionStepData']['PRESSURE'])

        # setting correct time step and reading ref pres val
        if idx == 0.0:
            pres_ref = 0.0
        else:
            sel_step = counter * step_ratio - 1
            pres_ref = ref_point['pres'][sel_step]

            if (abs(idx - ref_point['time'][sel_step]) > 1e-5):
                msg = 'Time stampt from reference point '
                msg += str(ref_point['time'][sel_step]) + ' [s]\n'
                msg += 'does not match the one from H5 file '
                msg += '{0:.2f}'.format(np.round(idx, 2)) + ' [s]!\n'
                msg += 'Check syntax or input.'

                raise Exception(msg)
            else:
                print('H5 file time stamp ' +
                      '{0:.2f}'.format(np.round(idx, 2)) + ' [s]')
                print('Ref pressure time stamp ' +
                      str(ref_point['time'][sel_step]) + ' [s]')

        # calculate q
        q_vals = pres_vals - pres_ref      

        # add q to h5 
        if "ResultsData/NodalSolutionStepData/Q" in h5f:
            print("Q is available, overwriting existing data")

            data = h5f['ResultsData/NodalSolutionStepData/Q']
            data[...] = q_vals     
            np.allclose(h5f['ResultsData/NodalSolutionStepData/Q'].value, q_vals)

        else:
            print("Q is not available, creating dataset")

            group = h5f.get('ResultsData/NodalSolutionStepData')
            group.create_dataset('Q', data=q_vals)

        # calculate cp
        cp_vals = q_vals * fctr
        
        # add cp to h5
        if "ResultsData/NodalSolutionStepData/Cp" in h5f:
            print("Cp is available, overwriting existing data")

            data = h5f['ResultsData/NodalSolutionStepData/Cp']
            data[...] = cp_vals     
            np.allclose(h5f['ResultsData/NodalSolutionStepData/Cp'].value, cp_vals)
        else:
            print("Cp is not available, creating dataset")

            group = h5f.get('ResultsData/NodalSolutionStepData')
            group.create_dataset('Cp', data=cp_vals)

print("Done adding Q and Cp to H5 data")
