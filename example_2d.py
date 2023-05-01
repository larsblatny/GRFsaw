# Author: Lars Blatny
from grfsaw.microstructure import Microstructure
import numpy as np
########################### Parameters ########################################

# Name of output folder
folder = "output/example_2d"

# Number of processors
procs = 4

# Choose if single-cut or double-cut structure
single_cut = True

# Number of mesh points in x,y,z-direction
shape = [300, 300]

# Solid volume fraction
phi = 0.6
phi_tol = 0.01 # set error tolerance on phi, see below while-loop

# Number of waves. Must be a multiple of procs. Must be divisible by 2 if anisotropy is used.
N = 10000

# Type of sampling distrbution, 'n' for normal or 'g' for gamma distribution
distr = 'n'

# Number of microstructural elements per side length in x-direction
m_mean = 13

# Standard deviation
m_std = 0.01

# Anisotropy parameter between 0 (anisotropic) and 1 (isotropic)
a = 1.0

# Anisotropy preferred direction, either vertical "v" or horizontal "h". Any other input is treated as isotropy.
anitype = 'v'

###############################################################################

phi_obs = -1
iter = 1
while (abs(phi_obs-phi) > phi_tol or iter == 1):
    print("Iteration = ", iter)
    iter += 1
    seed = np.random.randint(1, 1e5)
    m = Microstructure(shape, phi, folder, single_cut, N, distr, m_mean, m_std, anitype, a, seed, procs)
    m.create()
    phi_obs = np.sum(m.binary_array) / m.num_mesh_points

m.save_parameters()
m.create()
m.plot_distributions()

m.save(vtk=False) # VTK is only for 3D structures atm
m.plot()
m.autocorrelation(cleaned=False, one_dim=True)
m.find_shortest_path()

m.cleaning()
m.save(cleaned=True, vtk=False)
m.plot(cleaned=True)
