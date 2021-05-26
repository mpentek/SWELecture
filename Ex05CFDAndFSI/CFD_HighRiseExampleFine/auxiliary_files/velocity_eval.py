#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: mate.pentek@tum.de 

Description: Sample line plot for velocity results

Created on:  24.01.2020
Last update: 24.01.2020
'''
#===============================================================================

import numpy as np
import matplotlib.pyplot as plt

################################################################
# user input

file_prefix = 'line2Hup_'
n_files = 7
file_idxs = list(range(1, n_files + 1))

ramp_up_time = 30
ramp_up_fctr = 1.5
discarded_time = ramp_up_fctr * ramp_up_time

variable_labels = ['time', 'velx']

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


def get_turb_intens(vel_mean, vel_std, tol=1e-6):
    return 100 * vel_std/vel_mean


def get_euclidean_distance(coord_new, coord_old):
    return np.linalg.norm(coord_new - coord_old)

################################################################
# read in raw date of time series


raw_data = []
for idx in file_idxs:
    file_name = file_prefix + str(idx) + '.dat'
    raw_data.append(get_point_data_from_file(
        file_name, variable_labels))

################################################################
# postprocess raw data


discard_idx = get_discard_index(raw_data[0]['time'], discarded_time)
print('Raw data from time series discarded until ' +
      str(raw_data[0]['time'][discard_idx]) + ' [s]')

postprocessed_data = {
    'p_id': [],
    'running_coord': [],
    'vel_mean': [],
    'vel_ti': []
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

    postprocessed_data['vel_mean'].append(
        get_mean(data['velx'][discard_idx:]))
    postprocessed_data['vel_ti'].append(get_turb_intens(
        postprocessed_data['vel_mean'][-1], get_std(data['velx'][discard_idx:])))

################################################################
# plot raw data

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Velocity')

ax1.set_title('Mean velocity')
ax1.plot(postprocessed_data['vel_mean'],
         postprocessed_data['running_coord'], 'ko-')
ax1.set_ylabel('Height [m]')
ax1.set_xlabel('Velocity [m/s]')

ax2.set_title('Turbulence intensity')
ax2.plot(postprocessed_data['vel_ti'],
         postprocessed_data['running_coord'], 'ko-')
ax2.set_ylabel('Height [m]')
ax2.set_xlabel('Ti [%]')

plt.tight_layout()
plt.show()
