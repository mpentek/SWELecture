#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: mate.pentek@tum.de 

Description: Sample line plot for pressure results

Created on:  24.01.2020
Last update: 24.01.2020
'''
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt

################################################################
# user input

ref_point_file = 'ref_p.dat'
rho = 1.2

file_prefix = 'front_'
n_files = 7
file_idxs = list(range(1, n_files + 1))

ramp_up_time = 30
ramp_up_fctr = 1.5
discarded_time = ramp_up_fctr * ramp_up_time

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


def get_std(series):
    return np.std(series)


def get_rms(series):
    return np.sqrt(np.mean(series**2))


def get_cp(pres_i_series, pres_ref_series, fctr):
    return np.multiply(np.subtract(pres_i_series, pres_ref_series), fctr)


def get_euclidean_distance(coord_new, coord_old):
    return np.linalg.norm(coord_new - coord_old)

################################################################
# read in raw date of time series


# from pressure on surface
raw_data = []
for idx in file_idxs:
    file_name = file_prefix + str(idx) + '.dat'
    raw_data.append(get_point_data_from_file(
        file_name, variable_labels))


# from reference point
variable_labels.append('velx')
ref_point = get_point_data_from_file(
    ref_point_file, variable_labels)

################################################################
# postprocess raw data

discard_idx = get_discard_index(raw_data[0]['time'], discarded_time)
print('Raw data from time series discarded until ' +
      str(raw_data[0]['time'][discard_idx]) + ' [s]')

ref_velx = get_mean(ref_point['velx'][discard_idx:])
print('Reference velocity of ' + str(ref_velx) +
      ' [m/s] at height ' + str(ref_point['p_position'][2]) + ' [m]')
fctr = 1 / (0.5 * rho * ref_velx**2)

postprocessed_data = {
    'p_id': [],
    'running_coord': [],
    'cp_mean': [],
    'cp_std': [],
    'cp_rms': []
}

for counter, data in enumerate(raw_data):
    postprocessed_data['p_id'].append(data['p_id'])

    # using here a running coordinate
    # can be seen as a parametric coordinate in 1D for the line
    # set to 0.0 for the first point on the line
    if counter == 0:
        running_coord = 0.0
    else:
        running_coord += get_euclidean_distance(
            data['p_position'], last_position)
    last_position = data['p_position']
    postprocessed_data['running_coord'].append(running_coord)

    cp_series = get_cp(data['pres'], ref_point['pres'], fctr)

    postprocessed_data['cp_mean'].append(get_mean(cp_series[discard_idx:]))
    postprocessed_data['cp_std'].append(get_std(cp_series[discard_idx:]))
    postprocessed_data['cp_rms'].append(get_rms(cp_series[discard_idx:]))

################################################################
# plot raw data

fig, (ax1, ax2, ax3) = plt.subplots(3)
fig.suptitle('Pressure coefficient')

ax1.set_title('Mean')
ax1.plot(postprocessed_data['running_coord'],
         postprocessed_data['cp_mean'], 'ko-')
ax1.set_xlabel('Coord. [m]')
ax1.set_ylabel('Cp [-]')

ax2.set_title('Standard deviation')
ax2.plot(postprocessed_data['running_coord'],
         postprocessed_data['cp_std'], 'ko-')
ax2.set_xlabel('Coord. [m]')
ax2.set_ylabel('Cp [-]')

ax3.set_title('Root mean square')
ax3.plot(postprocessed_data['running_coord'],
         postprocessed_data['cp_rms'], 'ko-')
ax3.set_xlabel('Coord. [m]')
ax3.set_ylabel('Cp [-]')

plt.tight_layout()
plt.show()
