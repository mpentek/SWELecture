import os
import numpy as np
from matplotlib.pylab import plt

# TODO: user: change global file name and location
#result_folder = os.path.join(*['results','ascii_output','forces'])

#file_name_f = os.path.join(result_folder,'drag_structure.dat')
file_name_f = 'FluidModelPart.Drag_structure_drag.dat'

# TODO: user: change global file name and location
#file_name_fb  = os.path.join(result_folder,'force_structure.dat')
file_name_fb  = 'FluidModelPart.Drag_structure_global_force.dat'

# TODO: user: change number of levels (case) and file name prefix
# level force
case = 11
#folder_name_level  = os.path.join(result_folder,'level_force')
folder_name_level  = 'level_force'
file_prefix_level_fb = 'FluidModelPart.Drag_structure_level_force_'

# shift index because own processes start a step earlier to write
shift_idx = 1

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

def calculate_and_print_error(cumul_res, ref_res,consider_fb=False):
    labels = ['fx', 'fy', 'fz']
    if consider_fb:
        labels += ['mx', 'my', 'mz']
        labels += ['fx\'', 'fy\'', 'fz\'', 'mx\'', 'my\'', 'mz\'']

    calc_error = {}
    for label in labels:
        calc_error[label] = max([abs(x-y) for x,y in zip(ref_res[label],cumul_res[label])]) / max(abs(ref_res[label]))

    for key, value in calc_error.items():
        print('Error for ' + key + ': ' + str(value))

    print()    

# ===============================================================
global_results_f = {
    # time
    't':  np.array(np.loadtxt(file_name_f, usecols=(0,))),
    # flow-attached
    'fx': np.array(np.loadtxt(file_name_f, usecols=(1,))),
    'fy': np.array(np.loadtxt(file_name_f, usecols=(2,))),
    'fz': np.array(np.loadtxt(file_name_f, usecols=(3,)))
}

# REFERENCE RESULTS at base - global force file

global_results_fb = {
    # time
    't':  np.array(np.loadtxt(file_name_fb, usecols=(0,)))[shift_idx:],
    # flow-attached
    'fx': np.array(np.loadtxt(file_name_fb, usecols=(1,)))[shift_idx:],
    'fy': np.array(np.loadtxt(file_name_fb, usecols=(2,)))[shift_idx:],
    'fz': np.array(np.loadtxt(file_name_fb, usecols=(3,)))[shift_idx:],
    'mx': np.array(np.loadtxt(file_name_fb, usecols=(4,)))[shift_idx:],
    'my': np.array(np.loadtxt(file_name_fb, usecols=(5,)))[shift_idx:],
    'mz': np.array(np.loadtxt(file_name_fb, usecols=(6,)))[shift_idx:],
    # body-attached
    'fx\'': np.array(np.loadtxt(file_name_fb, usecols=(7,)))[shift_idx:],
    'fy\'': np.array(np.loadtxt(file_name_fb, usecols=(8,)))[shift_idx:],
    'fz\'': np.array(np.loadtxt(file_name_fb, usecols=(9,)))[shift_idx:],
    'mx\'': np.array(np.loadtxt(file_name_fb, usecols=(10,)))[shift_idx:],
    'my\'': np.array(np.loadtxt(file_name_fb, usecols=(11,)))[shift_idx:],
    'mz\'': np.array(np.loadtxt(file_name_fb, usecols=(12,)))[shift_idx:]
}

# ===============================================================
# READ-IN EXISTING - level force files
all_level_results = {}

# NOTE: for the actual run
for i in range(case):
    file_name = os.path.join(folder_name_level, file_prefix_level_fb + str(i) + '.dat')

    with open(file_name) as f:
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
        't':  np.array(np.loadtxt(file_name, usecols=(0,)))[shift_idx:],
        # flow-attached
        'fx': np.array(np.loadtxt(file_name, usecols=(1,)))[shift_idx:],
        'fy': np.array(np.loadtxt(file_name, usecols=(2,)))[shift_idx:],
        'fz': np.array(np.loadtxt(file_name, usecols=(3,)))[shift_idx:],
        'mx': np.array(np.loadtxt(file_name, usecols=(4,)))[shift_idx:],
        'my': np.array(np.loadtxt(file_name, usecols=(5,)))[shift_idx:],
        'mz': np.array(np.loadtxt(file_name, usecols=(6,)))[shift_idx:],
        # body-attached
        'fx\'': np.array(np.loadtxt(file_name, usecols=(7,)))[shift_idx:],
        'fy\'': np.array(np.loadtxt(file_name, usecols=(8,)))[shift_idx:],
        'fz\'': np.array(np.loadtxt(file_name, usecols=(9,)))[shift_idx:],
        'mx\'': np.array(np.loadtxt(file_name, usecols=(10,)))[shift_idx:],
        'my\'': np.array(np.loadtxt(file_name, usecols=(11,)))[shift_idx:],
        'mz\'': np.array(np.loadtxt(file_name, usecols=(12,)))[shift_idx:]
    }

# ===============================================================
print('EVALUATED ERROR: global_f vs global_fb')
calculate_and_print_error(global_results_f, global_results_fb, False)

print('EVALUATED ERROR: all level_fb vs global_fb')
calculate_and_print_error(get_cumulated_results(all_level_results), global_results_fb, True)




