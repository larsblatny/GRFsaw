# Author: Lars Blatny
import numpy as np

def autocorrelation_fft(binary_array):
    """
    autocorrelation_fft computes the (non-normalized) angular-averaged two-point correlation (autocorrelation) function.
    It is based on utilizing Fast Fourier Transform (FFT).
    According to the Wiener-Khinchin Theorem, the autocorrelation is given by
    the Fourier transform of the absolute square of the fourier coefficients of binary_array.
    Assumption: the microstructure is wrapped around its boundaries (i.e., PBC).

    Arguments:
    binary_array: binary two- or three-dimensional numpy array denoting solid (1) and void (0) points.

    Returns:
    numpy array of the two-point correlation function. 
    """

    shape = binary_array.shape
    Nx = shape[0]
    Ny = shape[1]

    if binary_array.ndim == 2:
        # print("    2D case detected: ", binary_array.shape)
        threedim = False
    else:
        threedim = True
        Nz = shape[2]

    fft_coeff_sq = np.abs(np.fft.fftn(binary_array))**2
    inv_fft = np.fft.ifftn(fft_coeff_sq).real.round() # this is the n-D autocorrelaton

    ### the max value of the autocorrelaton is found at r=(0,0,0)
    if threedim:
        inv_fft /= inv_fft[0,0,0]
    else:
        inv_fft /= inv_fft[0,0]

    ### the distances from the origin (0,0,0) are given by dist_3d, assuming dx = 1
    dist_x_sq  = ( np.fft.fftfreq(Nx) * Nx )**2
    dist_y_sq  = ( np.fft.fftfreq(Ny) * Ny )**2
    if threedim:
        dist_z_sq  = ( np.fft.fftfreq(Nz) * Nz )**2
        dist_3d_sq = np.zeros( (Nx,Ny,Nz) )
        for i in range(0, Nx):
            for j in range(0, Ny):
                for k in range(0, Nz):
                     dist_3d_sq[i, j, k] = dist_x_sq[i] + dist_y_sq[j] + dist_z_sq[k]
    else:
        dist_3d_sq = np.zeros( (Nx,Ny) )
        for i in range(0, Nx):
            for j in range(0, Ny):
                dist_3d_sq[i, j] = dist_x_sq[i] + dist_y_sq[j]

    dist_3d = np.sqrt(dist_3d_sq)

    ### the unique distances (points with the same magnitude r from r=0)
    corr_distances, index = np.unique(dist_3d, return_inverse=True)

    ### the angular-averaged autocorrelation for each distance r
    corr_fft_values = np.bincount(index, weights=inv_fft.flatten(order='C')) / np.bincount(index)
    # Here the
    #       Denominator = number of occurances of same distance r (i.e., number of points on the same spherical shell)
    #       Nominator = sum that have the same unique distance r

    return corr_fft_values, corr_distances
