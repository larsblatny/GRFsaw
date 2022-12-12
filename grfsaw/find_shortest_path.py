# Author: Lars Blatny
import numpy as np

###############################################################################
# Find shortest path of the spanning cluster
###############################################################################

def find_shortest_path(binary, void_phase=False):
    """
    find_shortest_path find shortest path of the spanning cluster by applying the Burning Method from all 6 sides.

    Arguments:
    binary: binary two- or three-dimensional numpy array denoting solid (1) and void (0) points.
    void_phase: boolean indicating wether calculating the path through the void phase (True) or solid phase (False, default).

    Returns:
    list of 2 (if 2D) or 3 (if 3D) integers representing the shortest path through the structure in terms of number of mesh points. 
    """

    binary = binary.astype(int)

    if void_phase:
        binary = 1 - binary

    if len(binary.shape) == 3:
        return find_shortest_path_3d(binary)
    else:
        return find_shortest_path_2d(binary)


def find_shortest_path_3d(binary_in):

    Nx = binary_in.shape[0]
    Ny = binary_in.shape[1]
    Nz = binary_in.shape[2]

    dir = ['x', 'y', 'z']
    shortest_path = np.zeros(3).astype(int)
    for case in range(0,3):

        binary = np.copy(binary_in)

        ####### Label Top Surface = 2 #######
        if case == 0:
            for j in range(0, Ny):
                for k in range(0, Nz):
                    binary[0, j, k] = 2
        elif case == 1:
            for i in range(0, Nx):
                for k in range(0, Nz):
                    binary[i, 0, k] = 2
        else: # case == 2
            for i in range(0, Nx):
                for j in range(0, Ny):
                    binary[i, j, 0] = 2
        #####################################

        # Helper function for easy handling of boundaries
        def super_array(i,j,k):
            if ((i < 0) or (i >= Nx) or (j < 0) or (j >= Ny) or (k < 0) or (k >= Nz)):
                return 0
            else:
                return binary[i,j,k]

        label = 2
        number_of_unoccupied_neighbors = 1 # a random positive value so we enter the while loop
        not_at_bottom = True
        while (not_at_bottom and number_of_unoccupied_neighbors > 0):
            number_of_unoccupied_neighbors = 0
            for i in range(0, Nx):
                for j in range(0, Ny):
                    for k in range(0, Nz):
                        if binary[i,j,k] == label:
                            if case == 0:
                                if i == Nx-1:
                                    not_at_bottom = False
                            elif case == 1:
                                if j == Ny-1:
                                    not_at_bottom = False
                            else: # case == 2
                                if k == Nz-1:
                                    not_at_bottom = False

                            if super_array(i-1,j,k) == 1:
                                binary[i-1,j,k] = label+1
                                number_of_unoccupied_neighbors += 1
                            if super_array(i+1,j,k) == 1:
                                binary[i+1,j,k] = label+1
                                number_of_unoccupied_neighbors += 1

                            if super_array(i,j-1,k) == 1:
                                binary[i,j-1,k] = label+1
                                number_of_unoccupied_neighbors += 1
                            if super_array(i,j+1,k) == 1:
                                binary[i,j+1,k] = label+1
                                number_of_unoccupied_neighbors += 1

                            if super_array(i,j,k-1) == 1:
                                binary[i,j,k-1] = label+1
                                number_of_unoccupied_neighbors += 1
                            if super_array(i,j,k+1) == 1:
                                binary[i,j,k+1] = label+1
                                number_of_unoccupied_neighbors += 1
            #print("    Label: ", label, ". Bottom reached: ", not not_at_bottom, " (available neighbors: ", number_of_unoccupied_neighbors, ")")
            label += 1

        if not not_at_bottom:
            shortest_path[case] = label-2
            print("    Shortest path in " + dir[case] + "-direction: " + str(shortest_path[case]) + ". Side length is " + str(binary.shape[case]) )
        else:
            shortest_path[case] = -1
            print("    No spanning cluster in " + dir[case] + "-direction.")

    return shortest_path




def find_shortest_path_2d(binary_in):

    Nx = binary_in.shape[0]
    Ny = binary_in.shape[1]

    dir = ['x', 'y']
    shortest_path = np.zeros(2).astype(int)
    for case in range(0,2):

        binary = np.copy(binary_in)

        ####### Label Top Surface = 2 #######
        if case == 0:
            for j in range(0, Ny):
                binary[0, j] = 2
        else: # case == 1
            for i in range(0, Nx):
                binary[i, 0] = 2
        #####################################

        # Helper function for easy handling of boundaries
        def super_array(i,j):
            if ( (i < 0) or (i >= Nx) or (j < 0) or (j >= Ny) ):
                return 0
            else:
                return binary[i,j]

        label = 2
        number_of_unoccupied_neighbors = 1 # a random positive value so we enter the while loop
        not_at_bottom = True
        while (not_at_bottom and number_of_unoccupied_neighbors > 0):
            number_of_unoccupied_neighbors = 0
            for i in range(0, Nx):
                for j in range(0, Ny):
                    if binary[i,j] == label:
                        if case == 0:
                            if i == Nx-1:
                                not_at_bottom = False
                        else: # case == 1
                            if j == Ny-1:
                                not_at_bottom = False

                        if super_array(i-1,j) == 1:
                            binary[i-1,j] = label+1
                            number_of_unoccupied_neighbors += 1
                        if super_array(i+1,j) == 1:
                            binary[i+1,j] = label+1
                            number_of_unoccupied_neighbors += 1

                        if super_array(i,j-1) == 1:
                            binary[i,j-1] = label+1
                            number_of_unoccupied_neighbors += 1
                        if super_array(i,j+1) == 1:
                            binary[i,j+1] = label+1
                            number_of_unoccupied_neighbors += 1

            #print("    Label: ", label, ". Bottom reached: ", not not_at_bottom, " (available neighbors: ", number_of_unoccupied_neighbors, ")")
            label += 1

        if not not_at_bottom:
            shortest_path[case] = label-2
            print("    Shortest path in " + dir[case] + "-direction: " + str(shortest_path[case]) + ". Side length is " + str(binary.shape[case]) )
        else:
            shortest_path[case] = -1
            print("    No spanning cluster in " + dir[case] + "-direction.")

    return shortest_path
