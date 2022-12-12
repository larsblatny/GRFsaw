# Author: Lars Blatny
import numpy as np

def local_wave_sum_2d(p, N_local, q0, xgrid, ygrid, q_hat_1, q_hat_2, phase):
    '''
    Helper function for the (parallellized) summation of standing 2D sinusoidal waves.
    '''
    S_local = np.zeros_like(xgrid)
    for n in range(p*N_local, (p+1)*N_local):
    	S_local = S_local + np.cos( q0[n]*(q_hat_1[n] * xgrid + q_hat_2[n] * ygrid) + phase[n] )
    return S_local

def local_wave_sum_3d(p, N_local, q0, xgrid, ygrid, zgrid, q_hat_1, q_hat_2, q_hat_3, phase):
    '''
    Helper function for the (parallellized) summation of standing 3D sinusoidal waves. 
    '''
    S_local = np.zeros_like(xgrid)
    for n in range(p*N_local, (p+1)*N_local):
    	S_local = S_local + np.cos( q0[n]*(q_hat_1[n] * xgrid + q_hat_2[n] * ygrid + q_hat_3[n] * zgrid) + phase[n] )
    return S_local
