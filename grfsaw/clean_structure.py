# Author: Lars Blatny
import numpy as np

def clean_structure(binary_array):
    """
    clean_structure cleans a microstructure by removing standalone clusters.
    More precisely, it deletes any solid point (change label from 1 to 0) that do not belong to any spanning cluster.
    Accomplished by applying the Burning Method from all 6 sides without the "not-at-bottom"-check.
    Note that there can be more than one spanning cluster, and that if no spanning cluster exists all points are deleted.

    Arguments:
    binary_array: binary two- or three-dimensional numpy array denoting solid (1) and void (0) points.

    Returns:
    binary two- or three-dimensional numpy array denoting solid (1) and void (0) points. 
    """

    dim_increase = False
    if binary_array.ndim == 2:
        #print("    2D case detected: ", binary_array.shape)
        binary_array = binary_array[:, :, np.newaxis]
        dim_increase = True

    # Input: Solid points are labeled 1
    shape = binary_array.shape
    Nx = shape[0]
    Ny = shape[1]
    Nz = shape[2]

    binary_array = binary_array.astype(int)

    # Six different cases of labeling Top Surface 2
    for c in range(0,6):
        print("    Side: ", c+1, " / 6")
        if c==0:
            for i in range(0, Nx):
                for k in range(0, Nz):
                    if binary_array[i, 0, k] == 1:
                        binary_array[i, 0, k] = 2
        elif c==1:
            for i in range(0, Nx):
                for k in range(0, Nz):
                    if binary_array[i, Ny-1, k] == 1:
                        binary_array[i, Ny-1, k] = 2
        elif c==2:
            for i in range(0, Nx):
                for j in range(0, Ny):
                    if binary_array[i, j, 0] == 1:
                        binary_array[i, j, 0] = 2
        elif c==3:
            for i in range(0, Nx):
                for j in range(0, Ny):
                    if binary_array[i, j, Nz-1] == 1:
                        binary_array[i, j, Nz-1] = 2
        elif c==4:
            for j in range(0, Ny):
                for k in range(0, Nz):
                    if binary_array[0, j, k] == 1:
                        binary_array[0, j, k] = 2
        else: # c==5
            for j in range(0, Ny):
                for k in range(0, Nz):
                    if binary_array[Nx-1, j, k] == 1:
                        binary_array[Nx-1, j, k] = 2

        # Helper function for easy handling of boundaries
        def super_array(i,j,k):
            if ((i < 0) or (i >= Nx) or (j < 0) or (j >= Ny) or (k < 0) or (k >= Nz)):
                return 0
            else:
                return binary_array[i,j,k]

        label = 2
        number_of_unoccupied_neighbors = 1 # just a random positive value so we enter the while loop
        while (number_of_unoccupied_neighbors > 0):
            number_of_unoccupied_neighbors = 0
            for i in range(0, Nx):
                for j in range(0, Ny):
                    for k in range(0, Nz):
                        if binary_array[i,j,k] == label:

                            if super_array(i-1,j,k) == 1:
                                binary_array[i-1,j,k] = label+1
                                number_of_unoccupied_neighbors += 1
                            if super_array(i+1,j,k) == 1:
                                binary_array[i+1,j,k] = label+1
                                number_of_unoccupied_neighbors += 1

                            if super_array(i,j-1,k) == 1:
                                binary_array[i,j-1,k] = label+1
                                number_of_unoccupied_neighbors += 1
                            if super_array(i,j+1,k) == 1:
                                binary_array[i,j+1,k] = label+1
                                number_of_unoccupied_neighbors += 1

                            if super_array(i,j,k-1) == 1:
                                binary_array[i,j,k-1] = label+1
                                number_of_unoccupied_neighbors += 1
                            if super_array(i,j,k+1) == 1:
                                binary_array[i,j,k+1] = label+1
                                number_of_unoccupied_neighbors += 1

            #print("Label: ", label)
            #print("    Num of avail neighbors:  ", number_of_unoccupied_neighbors)
            label += 1

        # Set all points with label = 1 to label = 0.
        binary_array[binary_array == 1] = 0
        # Set all points with label > 1 to label = 1.
        binary_array[binary_array > 0] = 1

    if dim_increase:
        binary_array = np.squeeze(binary_array, axis=2)

    return binary_array
