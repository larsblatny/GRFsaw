# Author: Lars Blatny
import numpy as np
import multiprocessing as mp

def autocorrelation_xyz(binary_array, r_max, procs):
    """
    autocorrelation_xyz computes the (non-normalized) two-point correlation (autocorrelation) functions in the d=2 or 3 axis-aligned directions.

    Arguments:
    binary_array: binary two- or three-dimensional numpy array denoting solid (1) and void (0) points.
    r_max: maximum distance as an integer number of mesh points.
    procs: integer number of processors available for the computation.

    Returns:
    numpy arrays of the two-point correlation functions in x and y (and z if 3D) direction.
    """

    if binary_array.ndim == 2:
        threedim = False
        print("    r \t x-dir \t y-dir")
    else:
        threedim = True
        print("    r \t x-dir \t y-dir \t z-dir")

    num_loops = int(r_max/procs)
    pool = mp.Pool(processes=procs)

    if threedim:
        two_point_corr_x = []
        two_point_corr_y = []
        two_point_corr_z = []
        for i in range(0, num_loops):
            outputs = [pool.apply_async( loop_3d, args=(r, binary_array) ) for r in range(i*procs, (i+1)*procs)]
            two_point_corr_x.extend( [out.get()[0] for out in outputs] )
            two_point_corr_y.extend( [out.get()[1] for out in outputs] )
            two_point_corr_z.extend( [out.get()[2] for out in outputs] )
        pool.close()
        return two_point_corr_x, two_point_corr_y, two_point_corr_z
    else:
        two_point_corr_x = []
        two_point_corr_y = []
        for i in range(0, num_loops):
            outputs = [pool.apply_async( loop_2d, args=(r, binary_array) ) for r in range(i*procs, (i+1)*procs)]
            two_point_corr_x.extend( [out.get()[0] for out in outputs] )
            two_point_corr_y.extend( [out.get()[1] for out in outputs] )
        pool.close()
        return two_point_corr_x, two_point_corr_y

def loop_3d(r, binary_array):
    shape = binary_array.shape
    tmp_x   = 0
    tmp_y   = 0
    tmp_z   = 0
    count_x = 0
    count_y = 0
    count_z = 0
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):
            for k in range(0, shape[2]):

                if((i+r) < shape[0]):
                    tmp_x   += binary_array[i,j,k] * binary_array[i+r,j,k]
                    count_x += 1
                if((i-r) >= 0):
                    tmp_x   += binary_array[i,j,k] * binary_array[i-r,j,k]
                    count_x += 1

                if((j+r) < shape[1]):
                    tmp_y   += binary_array[i,j,k] * binary_array[i,j+r,k]
                    count_y += 1
                if((j-r) >= 0):
                    tmp_y   += binary_array[i,j,k] * binary_array[i,j-r,k]
                    count_y += 1

                if((k+r) < shape[2]):
                    tmp_z   += binary_array[i,j,k] * binary_array[i,j,k+r]
                    count_z += 1
                if((k-r) >= 0):
                    tmp_z   += binary_array[i,j,k] * binary_array[i,j,k-r]
                    count_z += 1

    tpc_x = tmp_x / count_x
    tpc_y = tmp_y / count_y
    tpc_z = tmp_z / count_z

    print("    ", str(r) + "\t" + str(round(tpc_x,3)) + "\t" + str(round(tpc_y,3)) + "\t" + str(round(tpc_z,3)))

    return (tpc_x, tpc_y, tpc_z)


def loop_2d(r, binary_array):
    shape = binary_array.shape
    tmp_x   = 0
    tmp_y   = 0
    count_x = 0
    count_y = 0
    for i in range(0, shape[0]):
        for j in range(0, shape[1]):

            if((i+r) < shape[0]):
                tmp_x   += binary_array[i,j] * binary_array[i+r,j]
                count_x += 1
            if((i-r) >= 0):
                tmp_x   += binary_array[i,j] * binary_array[i-r,j]
                count_x += 1

            if((j+r) < shape[1]):
                tmp_y   += binary_array[i,j] * binary_array[i,j+r]
                count_y += 1
            if((j-r) >= 0):
                tmp_y   += binary_array[i,j] * binary_array[i,j-r]
                count_y += 1

    tpc_x = tmp_x / count_x
    tpc_y = tmp_y / count_y

    print("    ", str(r) + "\t" + str(round(tpc_x,3)) + "\t" + str(round(tpc_y,3)))

    return (tpc_x, tpc_y)
