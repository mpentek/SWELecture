# ===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek

Author: mate.pentek@tum.de, anoop.kodakkal@tum.de, catharina.czech@tum.de, peter.kupas@tum.de
    
Note: ...

Created on:  22.11.2017
Last update: 27.09.2019
'''
# ===============================================================================
import numpy as np
import math
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from matplotlib import animation

'''
Everything should boil down to 2 main visualizaiton types:
(1) "Static" (so not animated) image
(2) "Dynamic" (so animated) image

-> for both of these use a gray dashed line for the undeformed configurations
-> use with continuous black lines the deformed configuration
-> use dot marker for nodes
-> have an adjustable scaling factor for deformation
-> have an adjustable scaling factor for the given external forces
-> add the base reaction with the same factor as for forces
-> for forces use quiver plot


use version (1) for plotting 
    last static solve results 
    chosen eigenform and frequency
    a chosen time step for the dynamic simulation

use version (2) for animating 
    chosen eigenform and frequency
    the total duration of a dynamic simulation

'''

# only a maximum number of line plots is available
# these are formatted for : undefomed, (primary) deformed 0, (other) deformed 1, 2, 3
LINE_TYPE_SETUP = {"color":          ["grey", "black", "red", "green", "blue"],
                   "linestyle":      ["--",    "-",  "-",    "-",   "-"],
                   "marker":         ["o",    "s",  "^",    "p",   "x"],
                   "markeredgecolor": ["grey", "black", "red", "green", "blue"],
                   "markerfacecolor": ["grey", "black", "red", "green", "blue"],
                   "markersize":     [4,      4,    4,      4,    4, ]}

'''
geometry = {"undeformed":...,
            "deformation":...,
            "deformed": None}

where deformed = undeformed + deformation is done with these utilities 
taking scaling into consideration

geometry needs to contain the additional nodal information(s) when passed 
to the plot function

force = {"external":...,
         "base_reaction":...}
scaling_factor = {"deformation":...
                  "force":...}

plot_limits = {"x":[... , ...]
               "y":[... , ...]}
# defined based upon geometry["deformed"]
'''
def plot_dynamic_result(plot_title, result_data , array_time) : 

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(plot_title)
    plt.grid()
    plt.title(plot_title + ' Vs Time')    # set title
    # plot undeformed
    plt.plot(array_time, 
            result_data, 
            color=LINE_TYPE_SETUP["color"][1],
            linestyle=LINE_TYPE_SETUP["linestyle"][1],
            marker=LINE_TYPE_SETUP["marker"][1],
            markeredgecolor=LINE_TYPE_SETUP["markeredgecolor"][1],
            markerfacecolor=LINE_TYPE_SETUP["markerfacecolor"][1],
            markersize=LINE_TYPE_SETUP["markersize"][1])
    ax.legend()
    plt.show()


def plot_result(plot_title, geometry, force, scaling, n_data):

    # default parameter
    # if reaction_force is None:
    #     reaction_force = np.zeros(1)
    # if force is None:
    #     force = np.zeros((len(displacement), n_data))
    #     force_scaling_factor=0.0

    # TEST GEOMETRY
    #print('TESTING GEOMETRY:')
    # print(geometry)
    # Set up figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Handover data
    # x_undef = np.zeros(len(undeformed_geometry))
    # y_undef = undeformed_geometry
    # x_def = displacement*displacement_scaling_factor
    # x_force = force*force_scaling_factor

    # TODO: avoid try and except
    # use matrices straightforward
    try:
        # single / static case
        geometry["deformed"] = [np.add(geometry["undeformed"][0], geometry["deformation"][0] * scaling["deformation"]),
                                np.add(geometry["undeformed"][1], geometry["deformation"][1] * scaling["deformation"])]
    except:
        geometry["deformed"] = [np.add(geometry["undeformed"][0][:, np.newaxis], geometry["deformation"][0] * scaling["deformation"]),
                                np.add(geometry["undeformed"][1], geometry["deformation"][1] * scaling["deformation"])]

        pass
    # Axis, Grid and Label

    plot_limits = get_plot_limits(geometry["deformed"])
    ##

    ax.set_xlabel('Displacement(m)')
    ax.set_ylabel('Height(m)')

    # set axes, grid
    ax.set_xlim(plot_limits["x"][0], plot_limits["x"][1])
    ax.set_ylim(plot_limits["y"][0], plot_limits["y"][1])
    # ax.set_xticks(np.arange(xmin, xmax, 1))
    # ax.set_yticks(np.arange(ymin, ymax, 1))

    # Plot figure
    plt.grid()
    plt.title(plot_title)    # set title
    # plot undeformed
    plt.plot(geometry["undeformed"][0],
             geometry["undeformed"][1],
             color=LINE_TYPE_SETUP["color"][0],
             linestyle=LINE_TYPE_SETUP["linestyle"][0],
             marker=LINE_TYPE_SETUP["marker"][0],
             markeredgecolor=LINE_TYPE_SETUP["markeredgecolor"][0],
             markerfacecolor=LINE_TYPE_SETUP["markerfacecolor"][0],
             markersize=LINE_TYPE_SETUP["markersize"][0])

    if(n_data == 1):
        plt.plot(geometry["deformed"][0],
                 geometry["deformed"][1],
                 color=LINE_TYPE_SETUP["color"][1],
                 linestyle=LINE_TYPE_SETUP["linestyle"][1],
                 marker=LINE_TYPE_SETUP["marker"][1],
                 markeredgecolor=LINE_TYPE_SETUP["markeredgecolor"][1],
                 markerfacecolor=LINE_TYPE_SETUP["markerfacecolor"][1],
                 markersize=LINE_TYPE_SETUP["markersize"][1])

        try:
            plt.quiver(geometry["undeformed"][0],
                       geometry["undeformed"][1],
                       force["external"][0],
                       force["external"][1],
                       color="red")

            plt.quiver(geometry["undeformed"][0][0],
                       geometry["undeformed"][0][1],
                       force["base_reaction"][0][0],
                       force["base_reaction"][0][1],
                       color="green")

        except:
            # forces are None
            pass

    # multiple func in one plot
    elif (n_data < 4):
        for i in range(n_data):
            # TODO not using current formatting yet, needs update
            plt.plot(geometry["deformed"][0][:, i],
                     geometry["deformed"][1], label="mode " + str(i+1),
                     color=LINE_TYPE_SETUP["color"][i+1],
                     linestyle=LINE_TYPE_SETUP["linestyle"][i+1],
                     marker=LINE_TYPE_SETUP["marker"][i+1],
                     markeredgecolor=LINE_TYPE_SETUP["markeredgecolor"][i+1],
                     markerfacecolor=LINE_TYPE_SETUP["markerfacecolor"][i+1],
                     markersize=LINE_TYPE_SETUP["markersize"][i+1])

    else:
        raise Exception(" Plot format not supported for the request " +
                        str(n_data) + ", maximum 4 allowed.")

    # we do not have legend -> uncomnneted line ax.legend() to avoid waring: No labelleb objects found
    ax.legend()
    geometry = {"deformed": None}
    plt.show()


def animate_result(title, array_time, geometry, force, scaling, skip=1):

    # just for testing, where to get the undeformed structure?
    #
    # First set up the figure, the axis, and the plot element we want to
    # animate

    fig = plt.figure()

    # TODO: animate needs scaling as well
    # set min and max values
    #xmin = displacement_time_history.min()
    #xmax = displacement_time_history.max()
    #xmin = xmin - math.ceil((xmax-xmin)/10)
    #xmax = xmax + math.ceil((xmax-xmin)/10)

    # TODO extend and use plot limits

    for i in range(len(array_time)):
        geometry["deformed"][0][i] = np.add(
            geometry["undeformed"][0], geometry["deformation"][0][i] * scaling["deformation"])
        geometry["deformed"][1][i] = np.add(
            geometry["undeformed"][1], geometry["deformation"][1][i] * scaling["deformation"])

    xmin = np.min(geometry["deformed"][0])
    xmax = np.max(geometry["deformed"][0])
    #xmin = xmin - math.ceil((xmax-xmin)/30)
    #xmax = xmax + math.ceil((xmax-xmin)/30)

    ymin = np.min(geometry["deformed"][1])
    ymax = np.max(geometry["deformed"][1])
    ymin = ymin - math.ceil((ymax-ymin)/30)
    ymax = ymax + math.ceil((ymax-ymin)/30)

    ax = plt.axes(xlim=(xmin, xmax), ylim=(ymin, ymax))
    text = ax.text(0.02, 0.90, '', transform=ax.transAxes)
    text_mode = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    ax.set_xlabel("Displacement(m)")
    ax.set_ylabel("Height(m)")
    ax.set_title(title)
    ax.grid(True)

    undeformed_line, = ax.plot(geometry["undeformed"][0],
                               geometry["undeformed"][1],
                               color=LINE_TYPE_SETUP["color"][0],
                               linestyle=LINE_TYPE_SETUP["linestyle"][0],
                               marker=LINE_TYPE_SETUP["marker"][0],
                               markeredgecolor=LINE_TYPE_SETUP["markeredgecolor"][0],
                               markerfacecolor=LINE_TYPE_SETUP["markerfacecolor"][0],
                               markersize=LINE_TYPE_SETUP["markersize"][0])

    deformed_line, = ax.plot([], [],
                             color=LINE_TYPE_SETUP["color"][1],
                             linestyle=LINE_TYPE_SETUP["linestyle"][1],
                             marker=LINE_TYPE_SETUP["marker"][1],
                             markeredgecolor=LINE_TYPE_SETUP["markeredgecolor"][1],
                             markerfacecolor=LINE_TYPE_SETUP["markerfacecolor"][1],
                             markersize=LINE_TYPE_SETUP["markersize"][1])

    # initialization function: plot the background of each frame
    def init():
        deformed_line.set_data([], [])
        text.set_text('')
        text.set_text('')
        return deformed_line,

    # animation function.  This is called sequentially

    def animate(i):
        x = geometry["deformed"][0][skip*i]
        y = geometry["deformed"][1][skip*i]

        deformed_line.set_data(x, y)
        text.set_text('{0:.2f}'.format(array_time[skip*i]) + "[s]")
        return deformed_line, text
    # call the animator.  blit=True means only re-draw the parts that have
    # changed.
    # frames = number of columns in result
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=len(geometry["deformed"][0][::skip]), interval=50, blit=True)
    # interval is 50 for 20 frames per second
    plt.show()


def get_plot_limits(deformed_geometry, offset_factor=10.):
    try:
        # case of dynamic / multi plot
        x_min = np.matrix.min(deformed_geometry[0])
        x_max = np.matrix.max(deformed_geometry[0])
        y_min = np.matrix.min(deformed_geometry[1])
        y_max = np.matrix.max(deformed_geometry[1])
    except:
        # case of static / single plot
        x_min = np.amin(deformed_geometry[0])
        x_max = np.amax(deformed_geometry[0])
        y_min = np.amin(deformed_geometry[1])
        y_max = np.amax(deformed_geometry[1])

    plot_limits = {"x": [x_min - (x_max - x_min) / offset_factor,
                         x_max + (x_max - x_min) / offset_factor],
                   "y": [y_min - (y_max - y_min) / offset_factor,
                         y_max + (y_max - y_min) / offset_factor]}

    return plot_limits
