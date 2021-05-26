#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: mate.pentek@tum.de 

Description: Convert Kratos level forces into ParOptBeam format

Created on:  06.01.2020
Last update: 06.01.2020
'''
#===============================================================================

import os
import numpy as np
from matplotlib.pylab import plt


# ===============================================================
# USER INPUT

# TODO: user: change global file name and location
file_name_global = 'FluidModelPart.Drag_structure_global_force.dat'

# TODO: user: change number of levels (case) and file name prefix
# level force
case = 11
folder_name_level = 'level_force'
file_prefix_level = 'FluidModelPart.Drag_structure_level_force_'

# TODO: for ParOptBeam parametrization
number_of_sampling_interval_cases = [10+1]
start_point_position = [0.0, 0.0, 0.0]
end_point_position = [0.0, 0.0, 600.8]

# ===============================================================
# FUNCTION DEFINITION

# NOTE: assumed to be ordered from lowest to highers level
def get_cumulated_results(multiple_level_results):
    cumul_results = {
        # 't' : None,
        # 'fx' : None,
        # # 'fy' : None,
        # # 'fz' : None,
        # # 'mx' : None,
        # # 'my' : None,
        # # 'mz' : None
    }

    for i in range(len(multiple_level_results)):      
        if i == 0:
            # initialize
            cumul_results['t'] = multiple_level_results[i]['t']

            # flow-attached
            cumul_results['fx'] = multiple_level_results[i]['fx']
            cumul_results['fy'] = multiple_level_results[i]['fy']
            cumul_results['fz'] = multiple_level_results[i]['fz']
            cumul_results['mx'] = multiple_level_results[i]['mx'] - multiple_level_results[i]['fy'] * multiple_level_results[i]['z']
            cumul_results['my'] = multiple_level_results[i]['my'] + multiple_level_results[i]['fx'] * multiple_level_results[i]['z']
            cumul_results['mz'] = multiple_level_results[i]['mz']
        
            # body-attached
            cumul_results['fx\''] = multiple_level_results[i]['fx\'']
            cumul_results['fy\''] = multiple_level_results[i]['fy\'']
            cumul_results['fz\''] = multiple_level_results[i]['fz\'']
            cumul_results['mx\''] = multiple_level_results[i]['mx\''] - multiple_level_results[i]['fy\''] * multiple_level_results[i]['z']
            cumul_results['my\''] = multiple_level_results[i]['my\''] + multiple_level_results[i]['fx\''] * multiple_level_results[i]['z']
            cumul_results['mz\''] = multiple_level_results[i]['mz\'']

        else:
            # accumulate

            # flow-attached
            cumul_results['fx'] = np.add(cumul_results['fx'], multiple_level_results[i]['fx'])
            cumul_results['fy'] = np.add(cumul_results['fy'], multiple_level_results[i]['fy'])
            cumul_results['fz'] = np.add(cumul_results['fz'], multiple_level_results[i]['fz'])
            cumul_results['mx'] = np.add(cumul_results['mx'], multiple_level_results[i]['mx'] - multiple_level_results[i]['fy'] * multiple_level_results[i]['z'])
            cumul_results['my'] = np.add(cumul_results['my'], multiple_level_results[i]['my'] + multiple_level_results[i]['fx'] * multiple_level_results[i]['z'])
            cumul_results['mz'] = np.add(cumul_results['mz'], multiple_level_results[i]['mz'])
    
            # body-attached
            cumul_results['fx\''] = np.add(cumul_results['fx\''], multiple_level_results[i]['fx\''])
            cumul_results['fy\''] = np.add(cumul_results['fy\''], multiple_level_results[i]['fy\''])
            cumul_results['fz\''] = np.add(cumul_results['fz\''], multiple_level_results[i]['fz\''])
            cumul_results['mx\''] = np.add(cumul_results['mx\''], multiple_level_results[i]['mx\''] - multiple_level_results[i]['fy\''] * multiple_level_results[i]['z'])
            cumul_results['my\''] = np.add(cumul_results['my\''], multiple_level_results[i]['my\''] + multiple_level_results[i]['fx\''] * multiple_level_results[i]['z'])
            cumul_results['mz\''] = np.add(cumul_results['mz\''], multiple_level_results[i]['mz\''])
    

    return cumul_results


def calculate_and_print_error(cumul_res, ref_res,plot_results=False):
    labels = ['fx\'', 'fy\'', 'fz\'', 'mx\'', 'my\'', 'mz\'']
    calc_error = {}
    for label in labels:
        calc_error[label] = max([abs(x-y) for x,y in zip(ref_res[label],cumul_res[label])]) / max(abs(ref_res[label]))

    for key, value in calc_error.items():
        print('Error for ' + key + ': ' + str(value))

    print()

    if plot_results:
        # force
        plt.figure(1)
        plt.title('Fx\'')
        plt.plot(ref_res['t'], ref_res['fx\''], 'k.', alpha=0.5, label='ref')
        plt.plot(cumul_res['t'], cumul_res['fx\''], 'r--', label='cumul')
        plt.legend()

        plt.figure(2)
        plt.title('Fy\'')
        plt.plot(ref_res['t'], ref_res['fy\''], 'k.', alpha=0.5, label='ref')
        plt.plot(cumul_res['t'], cumul_res['fy\''], 'r--', label='cumul')
        plt.legend()

        plt.figure(3)
        plt.title('Fz\'')
        plt.plot(ref_res['t'], ref_res['fz\''], 'k.', alpha=0.5, label='ref')
        plt.plot(cumul_res['t'], cumul_res['fz\''], 'r--', label='cumul')
        plt.legend()

        # moments
        plt.figure(4)
        plt.title('Mx\'')
        plt.plot(ref_res['t'], ref_res['mx\''], 'k.', alpha=0.5, label='ref')
        plt.plot(cumul_res['t'], cumul_res['mx\''], 'r--', label='cumul')
        plt.legend()

        plt.figure(5)
        plt.title('My\'')
        plt.plot(ref_res['t'], ref_res['my\''], 'k.', alpha=0.5, label='ref')
        plt.plot(cumul_res['t'], cumul_res['my\''], 'r--', label='cumul')
        plt.legend()

        plt.figure(6)
        plt.title('Mz\'')
        plt.plot(ref_res['t'], ref_res['mz\''], 'k.', alpha=0.5, label='ref')
        plt.plot(cumul_res['t'], cumul_res['mz\''], 'r--', label='cumul')
        plt.legend()

        plt.show()
    

# ===============================================================
# REFERENCE RESULTS at base - global force file

global_results = {
    # time
    't':  np.array(np.loadtxt(file_name_global, usecols=(0,))),
    # flow-attached
    'fx': np.array(np.loadtxt(file_name_global, usecols=(1,))),
    'fy': np.array(np.loadtxt(file_name_global, usecols=(2,))),
    'fz': np.array(np.loadtxt(file_name_global, usecols=(3,))),
    'mx': np.array(np.loadtxt(file_name_global, usecols=(4,))),
    'my': np.array(np.loadtxt(file_name_global, usecols=(5,))),
    'mz': np.array(np.loadtxt(file_name_global, usecols=(6,))),
    # body-attached
    'fx\'': np.array(np.loadtxt(file_name_global, usecols=(7,))),
    'fy\'': np.array(np.loadtxt(file_name_global, usecols=(8,))),
    'fz\'': np.array(np.loadtxt(file_name_global, usecols=(9,))),
    'mx\'': np.array(np.loadtxt(file_name_global, usecols=(10,))),
    'my\'': np.array(np.loadtxt(file_name_global, usecols=(11,))),
    'mz\'': np.array(np.loadtxt(file_name_global, usecols=(12,)))
}

# ===============================================================
# READ-IN EXISTING - level force files

all_level_results = {}

# NOTE: for the actual run
for i in range(case):
    file_name_level = os.path.join(folder_name_level, file_prefix_level + str(i) + '.dat')

    with open(file_name_level) as f:
        for j, line in enumerate(f):
            if j == 1:
                third_line = f.readline()
                third_line = third_line.split()
                break

    # NOTE: not explicitly needed to store intermediately
    all_level_results[i] = {
        # level parameters
        'x': float(third_line[2][:-1]),
        'y': float(third_line[3][:-1]),
        'z': float(third_line[4]),
        # time
        't':  np.array(np.loadtxt(file_name_level, usecols=(0,))),
        # flow-attached
        'fx': np.array(np.loadtxt(file_name_level, usecols=(1,))),
        'fy': np.array(np.loadtxt(file_name_level, usecols=(2,))),
        'fz': np.array(np.loadtxt(file_name_level, usecols=(3,))),
        'mx': np.array(np.loadtxt(file_name_level, usecols=(4,))),
        'my': np.array(np.loadtxt(file_name_level, usecols=(5,))),
        'mz': np.array(np.loadtxt(file_name_level, usecols=(6,))),
        # body-attached
        'fx\'': np.array(np.loadtxt(file_name_level, usecols=(7,))),
        'fy\'': np.array(np.loadtxt(file_name_level, usecols=(8,))),
        'fz\'': np.array(np.loadtxt(file_name_level, usecols=(9,))),
        'mx\'': np.array(np.loadtxt(file_name_level, usecols=(10,))),
        'my\'': np.array(np.loadtxt(file_name_level, usecols=(11,))),
        'mz\'': np.array(np.loadtxt(file_name_level, usecols=(12,)))
    }

print('EVALUATED INITIAL RESULTS')
calculate_and_print_error(get_cumulated_results(all_level_results), global_results, True)
wait = input("check...")

# ===============================================================
# PARAMETRIZATION

# TODO: user: change number of sampling intervals


for number_of_sampling_intervals in number_of_sampling_interval_cases:
    print("RESULTS FOR CASE ", str(number_of_sampling_intervals))

    # setup the parametric space for the internal points on the line
    lower_bound = 0.0
    upper_bound = 1.0

    # initialize equidistant intervals
    my_param = [1.0] * (number_of_sampling_intervals)
    # overwrite first and last as endpoints are to be included
    my_param[0] = 0.5
    my_param[-1] = 0.5

    parametrized_internal_points = [lower_bound + x*(upper_bound-lower_bound)/(
        number_of_sampling_intervals-1) for x in range(number_of_sampling_intervals + 1)]

    parametrized_internal_points = [0.0] * (len(my_param)+1)

    for i in range(len(my_param)):
        for j in range(i+1):
            parametrized_internal_points[i+1] += my_param[j]

    parametrized_internal_points = [
        x / parametrized_internal_points[-1] for x in parametrized_internal_points]

    # print("Interval parameter: ", my_param)
    # print("Running parameter: ", parametrized_internal_points)

    # determining the positions of the output points
    direction_vector = [
        x - y for x, y in zip(end_point_position, start_point_position)]

    current_positions = []
    for k in range(len(parametrized_internal_points)):
        current_positions.append(
            [x + parametrized_internal_points[k]*y for x, y in zip(start_point_position, direction_vector)])

    levels = {}
    running_x_coords = []
    for idx in range(len(current_positions)-1):
        levels[idx] = {}
        levels[idx]['start'] = current_positions[idx]
        levels[idx]['end'] = current_positions[idx+1]

        # first interval
        if idx == 0:
            levels[idx]['center'] = [
                x1 for x1 in levels[idx]['start']]
        elif idx == (len(current_positions)-1)-1:
            levels[idx]['center'] = [
                x2 for x2 in levels[idx]['end']]
        else:
            levels[idx]['center'] = [
                (x1+x2)/2 for x1, x2 in zip(levels[idx]['start'], levels[idx]['end'])]

        levels[idx]['output_file'] = None

        # the center is where the forces are computed
        # this has to correspond to the x-coords of the beam
        running_x_coords.append(levels[idx]['center'][2])


    # ===============================================================
    # CHECK PARAMETRIZATION


    print('X coords: ', len(running_x_coords))
    print(running_x_coords)

    ints = [y-x for x, y in zip(running_x_coords[:-1], running_x_coords[1:])]
    print('Intervals: ', len(ints))
    print(ints)


    # ===============================================================
    # INITIALIZATION OF NEW DATA


    # level force
    case = number_of_sampling_intervals
    folder_name = str(case)
    file_prefix = 'level_'
    all_new_level_results = {}

    # NOTE: just for testing
    # for i in range(15):

    # NOTE: for the actual run
    for i in range(case):
        all_new_level_results[i] = {
            'x': 0.0,
            'y': 0.0,
            'z': running_x_coords[i],
            't':  all_level_results[0]['t'],

            # flow-attached forces
            'fx': np.zeros(len(all_level_results[0]['t'])),
            'fy': np.zeros(len(all_level_results[0]['t'])),
            'fz': np.zeros(len(all_level_results[0]['t'])),
            'mx': np.zeros(len(all_level_results[0]['t'])),
            'my': np.zeros(len(all_level_results[0]['t'])),
            'mz': np.zeros(len(all_level_results[0]['t'])),

            # only interested in body-attached forces
            'fx\'': np.zeros(len(all_level_results[0]['t'])),
            'fy\'': np.zeros(len(all_level_results[0]['t'])),
            'fz\'': np.zeros(len(all_level_results[0]['t'])),
            'mx\'': np.zeros(len(all_level_results[0]['t'])),
            'my\'': np.zeros(len(all_level_results[0]['t'])),
            'mz\'': np.zeros(len(all_level_results[0]['t']))
        }


    # ===============================================================
    # CALCULATION OF NEW LOADS


    # NOTE: assuming that
    # all_new_level_results -> has to have fewer points (levels)
    # all_level_results -> has to have more points (levels)
    for i in range(len(all_new_level_results)-1):
        
        for j in range(len(all_level_results)):
            if (j < len(all_level_results)-1 and all_new_level_results[i]['z']<= all_level_results[j]['z'] and all_level_results[j]['z'] < all_new_level_results[i+1]['z']) or (j == len(all_level_results)-1 and all_new_level_results[i]['z']<= all_level_results[j]['z'] and all_level_results[j]['z'] <= all_new_level_results[i+1]['z']):
                # print('Level ' + str(j) + ' is in between levels ' + str(i) + ' and ' + str(i+1) +'\n')

                # creating weighting factors for linear weighing of the results
                # dist to lower neighbour -> will be positive
                a = all_level_results[j]['z'] - all_new_level_results[i]['z']
                # dist to upper neighbour -> will be negative
                b = all_level_results[j]['z'] - all_new_level_results[i+1]['z']
                # above need to take abs() of value
                l = abs(a) + abs(b)
                fctr_lower = abs(b/l)
                fctr_upper = abs(a/l)

                # cumulate results

                # flow-attached forces
                # lower neighbour
                all_new_level_results[i]['fx'] = np.add(all_new_level_results[i]['fx'], fctr_lower * all_level_results[j]['fx'])
                all_new_level_results[i]['fy'] = np.add(all_new_level_results[i]['fy'], fctr_lower * all_level_results[j]['fy'])
                all_new_level_results[i]['fz'] = np.add(all_new_level_results[i]['fz'], fctr_lower * all_level_results[j]['fz'])
                all_new_level_results[i]['mx'] = np.add(all_new_level_results[i]['mx'], fctr_lower * (all_level_results[j]['mx'] - all_level_results[j]['fy'] * a))
                all_new_level_results[i]['my'] = np.add(all_new_level_results[i]['my'], fctr_lower * (all_level_results[j]['my'] + all_level_results[j]['fx'] * a))
                all_new_level_results[i]['mz'] = np.add(all_new_level_results[i]['mz'], fctr_lower * all_level_results[j]['mz'])
                                
                # upper neighbour
                all_new_level_results[i+1]['fx'] = np.add(all_new_level_results[i+1]['fx'], fctr_upper * all_level_results[j]['fx'])
                all_new_level_results[i+1]['fy'] = np.add(all_new_level_results[i+1]['fy'], fctr_upper * all_level_results[j]['fy'])
                all_new_level_results[i+1]['fz'] = np.add(all_new_level_results[i+1]['fz'], fctr_upper * all_level_results[j]['fz'])
                all_new_level_results[i+1]['mx'] = np.add(all_new_level_results[i+1]['mx'], fctr_upper * (all_level_results[j]['mx'] - all_level_results[j]['fy'] * b))
                all_new_level_results[i+1]['my'] = np.add(all_new_level_results[i+1]['my'], fctr_upper * (all_level_results[j]['my'] + all_level_results[j]['fx'] * b))
                all_new_level_results[i+1]['mz'] = np.add(all_new_level_results[i+1]['mz'], fctr_upper * all_level_results[j]['mz'])

                # only interested in body-attached forces
                # lower neighbour
                all_new_level_results[i]['fx\''] = np.add(all_new_level_results[i]['fx\''], fctr_lower * all_level_results[j]['fx\''])
                all_new_level_results[i]['fy\''] = np.add(all_new_level_results[i]['fy\''], fctr_lower * all_level_results[j]['fy\''])
                all_new_level_results[i]['fz\''] = np.add(all_new_level_results[i]['fz\''], fctr_lower * all_level_results[j]['fz\''])
                all_new_level_results[i]['mx\''] = np.add(all_new_level_results[i]['mx\''], fctr_lower * (all_level_results[j]['mx\''] - all_level_results[j]['fy\''] * a))
                all_new_level_results[i]['my\''] = np.add(all_new_level_results[i]['my\''], fctr_lower * (all_level_results[j]['my\''] + all_level_results[j]['fx\''] * a))
                all_new_level_results[i]['mz\''] = np.add(all_new_level_results[i]['mz\''], fctr_lower * all_level_results[j]['mz\''])
                                
                # upper neighbour
                all_new_level_results[i+1]['fx\''] = np.add(all_new_level_results[i+1]['fx\''], fctr_upper * all_level_results[j]['fx\''])
                all_new_level_results[i+1]['fy\''] = np.add(all_new_level_results[i+1]['fy\''], fctr_upper * all_level_results[j]['fy\''])
                all_new_level_results[i+1]['fz\''] = np.add(all_new_level_results[i+1]['fz\''], fctr_upper * all_level_results[j]['fz\''])
                all_new_level_results[i+1]['mx\''] = np.add(all_new_level_results[i+1]['mx\''], fctr_upper * (all_level_results[j]['mx\''] - all_level_results[j]['fy\''] * b))
                all_new_level_results[i+1]['my\''] = np.add(all_new_level_results[i+1]['my\''], fctr_upper * (all_level_results[j]['my\''] + all_level_results[j]['fx\''] * b))
                all_new_level_results[i+1]['mz\''] = np.add(all_new_level_results[i+1]['mz\''], fctr_upper * all_level_results[j]['mz\''])

            else:
                print('Level ' + str(j) + ' IS NOT in between levels ' + str(i) + ' and ' + str(i+1) +'\n')
                pass
    
    print('EVALUATED NEW RESULTS FOR CASE: ', str(number_of_sampling_intervals))
    calculate_and_print_error(get_cumulated_results(all_new_level_results), global_results, True)
    wait = input("check...")


    # ===============================================================
    # SAVE LOADS IN *.NPY format

    # here mapping from body-attached axis from CFD with Kratos
    # to body-attached axis of CSD with ParOptBeam hardcoded 

    dofs_per_node = 6
    dynamic_force = np.zeros([number_of_sampling_intervals * dofs_per_node, len(all_new_level_results[0]['t'])])
    counter = 0
    for i in range(len(all_new_level_results)):
        # 0 -> fz
        dynamic_force[i * dofs_per_node + 0, :] = all_new_level_results[i]['fz\'']
        # 1 -> fx
        dynamic_force[i * dofs_per_node + 1, :] = all_new_level_results[i]['fx\'']
        # 2 -> fy
        dynamic_force[i * dofs_per_node + 2, :] = all_new_level_results[i]['fy\'']
        # 3 -> mz
        dynamic_force[i * dofs_per_node + 3, :] = all_new_level_results[i]['mz\'']
        # 4 -> mx
        dynamic_force[i * dofs_per_node + 4, :] = all_new_level_results[i]['mx\'']
        # 5 -> my
        dynamic_force[i * dofs_per_node + 5, :] = all_new_level_results[i]['my\'']
    
    np.save('force_dynamic_' + str(number_of_sampling_intervals), dynamic_force) 

np.save('array_time', all_new_level_results[0]['t'])


