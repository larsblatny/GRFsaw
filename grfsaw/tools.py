# Author: Lars Blatny
import math
import numpy as np

def lcm(a, b):
    '''
    Finds the lowest common multiple (lcm) of a and b
    '''
    return abs(a*b) // math.gcd(a, b)

def num2txt(dir, name, num):
    '''
    Saves a number num to name.txt in the directory dir
    dir must end with "/"
    '''
    file = open(dir + name + ".txt", "w") # "w" for write to file
    file.write(str(num))
    file.close()

def set_axis_equal(ax):
    '''
    Make axes have equal scale for 3D plots.
    Inspired by: https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to/
    '''

    x_lim = ax.get_xlim3d()
    y_lim = ax.get_ylim3d()
    z_lim = ax.get_zlim3d()

    x_range = abs(x_lim[1] - x_lim[0])
    x_mid = np.mean(x_lim)
    y_range = abs(y_lim[1] - y_lim[0])
    y_mid = np.mean(y_lim)
    z_range = abs(z_lim[1] - z_lim[0])
    z_mid = np.mean(z_lim)

    rad = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mid-rad, x_mid+rad])
    ax.set_ylim3d([y_mid-rad, y_mid+rad])
    ax.set_zlim3d([z_mid-rad, z_mid+rad])
