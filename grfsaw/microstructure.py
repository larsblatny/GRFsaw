# Author: Lars Blatny

import os, sys, math
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from scipy.special import erfinv

from grfsaw import *

class Microstructure:
    '''
    Class for the microstructure to be generated.

    Attributes given by user:
    -------------------------
    shape : list
        list of 2 (if 2D) or 3 (if 3D) integers representing the (mesh) size of the microstructure as number of mesh points in each direction
    phi : float
        desired solid volume fraction = 1 - porosity
    folder : str
        name of directory to save data and plots (default: "output/")
    single_cut : boolean
        generate single-cut or double-cut structure (default True)
    N : int
        number of waves to use, the more the better but slower (default 10 000)
    distr : char
        type of distribution to use for wavevector magnitude sampling, 'n' for normal or 'g' for gamma distribution (default 'n')
    m_mean : float
        average number of microstructural elements per length in the x-direction ~ a measure of the grain size of the structure (default 3)
    m_std : float
        the standard deviation of m (default 0.001)
    anitype : char
        type of preferred direction, 'h' for horizontal or 'v' for vertical, anything else is interpreted as isotropic (default isotropic)
    a : float
        measure of the anisotropy level between 0 and 1 (default 1, meaning isotropic)
    seed : int
        seed for RNG (default np.random.randint(1, 1e5))
    procs : int
        number of cpus to use for the parallellized computations (default max number of cpus)

    Internal attributes:
    -----------------------
    q0 : 1d float numpy array
        list of wave vector magnitudes
    theta : 1d float numpy array
        list of wave vector directons
    binary_array : 2d or 3d boolean numpy array
        the level-cut GRF
    binary_array_cleaned : 2d or 3d boolean numpy array
        the level-cut GRF containing only the parts belonging to the spanning cluster

    Methods:
    -----------------------
    '''

    def __init__(self, shape, phi, folder="output/", single_cut=True, Nw=10000, distr='n', m_mean=3, m_std=0.001, anitype='i', a=1, seed=np.random.randint(1, 1e5), procs=mp.cpu_count()):
        self.shape = shape
        self.phi = phi
        self.folder = folder
        self.single_cut = single_cut
        self.Nw = int(max(procs*2, round(Nw / lcm(procs,2)) * lcm(procs,2) ))
        self.distr = distr
        self.m_mean = m_mean
        self.m_std = m_std
        self.a = a
        self.anitype = anitype
        self.seed = seed
        self.procs = procs

        self.q0 = None
        self.theta = None
        self.binary_array = None
        self.binary_array_cleaned = None

        print("Using ", self.procs, " processors and ", self.Nw, " waves")

        ### checks on params
        if (self.Nw % self.procs != 0):
            print("ERROR: Nw not divisble by procs")
            return

        if ((self.anitype == 'h' or self.anitype == 'v') and self.Nw % 2 != 0):
            print("ERROR: Nw must be divisble by 2 when using anisotropy")
            return

        if len(self.shape) == 2:
            self.num_mesh_points = self.shape[0] * self.shape[1]
        elif len(self.shape) == 3:
            self.num_mesh_points = self.shape[0] * self.shape[1] * self.shape[2]
        else:
            print("ERROR: Dimensions of structure must be 2 or 3")
            return

        if (self.folder[-1] != '/'):
            print("Adding trailing backslash to directory name")
            self.folder = self.folder + "/"

        m_std_max = self.m_mean / np.sqrt(51)
        if (self.distr == 'g' and self.m_std < m_std_max):
            self.m_std = m_std_max
            print("NB: m_std increased to", round(m_std_max, 3), "as too small STD not possible with gamma distribution")

        ### create directories
        os.makedirs(os.path.dirname(self.folder), exist_ok=True)
        os.makedirs(os.path.dirname(self.folder + "plots/"), exist_ok=True)

    def save_parameters(self):
        '''
        Saves setup parameters to folder for future reference.
        '''
        print("Saving setup parameters...")

        dir = self.folder + "params/"
        os.makedirs(os.path.dirname(dir), exist_ok=True)

        num2txt(dir, "phi",      self.phi)
        num2txt(dir, "Nw",       self.Nw)
        num2txt(dir, "distr",    self.distr)
        num2txt(dir, "m_mean",   self.m_mean)
        num2txt(dir, "m_std",    self.m_std)
        num2txt(dir, "a",        self.a)
        num2txt(dir, "anitype",  self.anitype)
        num2txt(dir, "seed",     self.seed)

        file = open(dir + "mesh.txt", "w") # "w" for write to file
        for i in range(0, len(self.shape)):
            file.write(str(self.shape[i]) + "\n")
        file.close()


    def create(self):
        '''
        Creates the microstructure with the given user parameters.
        '''
        print("Generating structure...")

        print("    Sampling wave vector magnitudes...")
        self.__sample_wave_magnitude()

        if len(self.shape) == 2:
            self.__create_2d()
        else:
            self.__create_3d()


    def __create_2d(self):
        '''
        Helper for create() for 2D structures.
        '''

        ### sample wave vector directions theta between 0 and pi (180deg)
        if self.anitype == 'h':
            theta1 = np.random.uniform((1-self.a)*np.pi/2, (1+self.a)*np.pi/2, int(self.Nw/2))
            theta2 = np.random.uniform((3-self.a)*np.pi/2, (3+self.a)*np.pi/2, int(self.Nw/2))
            self.theta = np.concatenate((theta1, theta2))
        elif self.anitype == 'v':
            theta1 = np.random.uniform(   -self.a*np.pi/2,     self.a*np.pi/2, int(self.Nw/2))
            theta2 = np.random.uniform((2-self.a)*np.pi/2, (2+self.a)*np.pi/2, int(self.Nw/2))
            self.theta = np.concatenate((theta1, theta2))
        else:
            print("    Assuming isotropy")
            self.theta = np.random.uniform(0, 2*np.pi, self.Nw) # if uniform sampling - isotropic sample

        ### sampling directions
        q_hat_1 = np.cos(self.theta)
        q_hat_2 = np.sin(self.theta)

        ### mesh
        x = np.arange(self.shape[0]).astype(int)
        y = np.arange(self.shape[1]).astype(int)
        xm, ym = np.meshgrid(x, y, indexing='ij')

        ### sampling phases
        phase = np.random.uniform(0, 2*np.pi, self.Nw)

        ### adding waves
        print("    Creating GRF...")
        N_local = int(self.Nw/self.procs)
        pool = mp.Pool(processes=self.procs)
        outputs = [pool.apply_async( local_wave_sum_2d, args=(p, N_local, self.q0, xm, ym, q_hat_1, q_hat_2, phase) ) for p in range(0, self.procs)]
        S_local = [out.get() for out in outputs]
        S = sum(S_local) * np.sqrt(1.0/self.Nw)
        pool.close()

        ### perform level-cut of GRF
        print("    Performing the level-cut of the GRF...")
        self.__levelcut(S)


    def __create_3d(self):
        '''
        Helper for create() for 3D structures.
        '''

        ### sample wave vector directions theta between 0 and pi (180deg)
        if self.anitype == 'h':
            self.theta = np.arccos(np.random.uniform(-self.a, self.a, self.Nw))
        elif self.anitype == 'v':
            theta1 = np.arccos(np.random.uniform(-1,       -(1-self.a), int(self.Nw/2)))
            theta2 = np.arccos(np.random.uniform(1-self.a,           1, int(self.Nw/2)))
            self.theta = np.concatenate((theta1, theta2))
        else:
            print(    "Assuming isotropy")
            self.theta = np.arccos(np.random.uniform(-1, 1, self.Nw)) # cos^{-1}( 2 * U(0,1) - 1)

        ### sampling directions
        azimuth_angle = np.random.uniform(0, 2*np.pi, self.Nw)
        q_hat_1 = np.sin(self.theta) * np.cos(azimuth_angle)
        q_hat_2 = np.sin(self.theta) * np.sin(azimuth_angle)
        q_hat_3 = np.cos(self.theta)

        ### mesh
        x = np.arange(self.shape[0]).astype(int)
        y = np.arange(self.shape[1]).astype(int)
        z = np.arange(self.shape[2]).astype(int)
        xm, ym, zm = np.meshgrid(x, y, z, indexing='ij')

        ### sampling phases
        phase = np.random.uniform(0, 2*np.pi, self.Nw)

        ### adding waves
        print("    Creating GRF...")
        N_local = int(self.Nw/self.procs)
        pool = mp.Pool(processes=self.procs)
        outputs = [pool.apply_async( local_wave_sum_3d, args=(p, N_local, self.q0, xm, ym, zm, q_hat_1, q_hat_2, q_hat_3, phase) ) for p in range(0, self.procs)]
        S_local = [out.get() for out in outputs]
        S = sum(S_local) * np.sqrt(1.0/self.Nw)
        pool.close()

        ### perform level-cut of GRF
        print("    Performing the level-cut of the GRF...")
        self.__levelcut(S)

    def __sample_wave_magnitude(self):
        '''
        Helper for create() - samples wave vector magnitudes from the given probability distribution.
        '''

        q_mean  = 2*np.pi * self.m_mean / (self.shape[0]-1)
        delta_q = 2*np.pi * self.m_std / (self.shape[0]-1)

        np.random.seed(self.seed)
        if self.distr == 'n':
            # self.q0 = np.random.normal(loc=q_mean, scale=delta_q, size=self.Nw)
            self.q0 = normal_distribution(mean=q_mean, std=delta_q).rvs(self.Nw)
        elif self.distr == 'g':
            b = (self.m_mean/self.m_std )**2 - 1
            self.q0 = gamma_distribution(a=0.0).rvs(size=self.Nw, b=b, mean=q_mean)
        else:
            print("    ERROR: Invalid sampling distrubution, distr must be 'n' or 'g'")
            return

    def __levelcut(self, S):
        '''
        Helper for create() - performs the level-cut of the Gaussian random field S.
        '''

        ### ALT 1: level cut -> solid phase ABOVE xi
        xi = erfinv(1 - 2 * self.phi)

        ### ALT 2: level cut -> solid phase BELOW xi
        # xi = erfinv(2 * phi - 1)

        if self.single_cut:

            self.binary_array = (S > xi)   # ALT 1
            # self.binary_array = (S < xi) # ALT 2
        else: # two-cut
            upper = erfinv(self.phi)
            lower = erfinv(-self.phi)

            self.binary_array = np.logical_and((S > lower), (S < upper))

        phi_observed = np.sum(self.binary_array) / self.num_mesh_points
        num2txt(self.folder, "phi_observed", phi_observed)
        print( "    Observed solid volume fraction: " + str(round(phi_observed,4)) )


    def plot_distributions(self, figsize=(13, 8), dpi=100, fontsize=22):
        '''
        Saves histograms of the sampling of wave lengths and directions.
        The optional parameters figsize, dpi and fontsize relates to the size and resolution of the plots.
        '''
        print("Plotting histograms of wave vector magnitudes and directions...")

        plt.rcParams.update({'font.size': fontsize})
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')

        plt.figure(figsize=figsize, dpi=dpi)
        plt.hist(self.q0*(self.shape[0]-1)/(2*np.pi), bins=50, density=False, color='k')
        plt.xlabel(r'$L_x/\lambda$')
        plt.ylabel(r'Number of samples')
        plt.savefig(self.folder + "plots/histogram_wavelengths.png", bbox_inches = 'tight')
        plt.close()

        plt.figure(figsize=figsize, dpi=dpi)
        plt.hist(np.degrees(self.theta), bins=50, density=False, color='k')
        plt.xlabel(r'$\theta$ (deg)')
        plt.ylabel(r'Number of samples')
        if (len(self.shape) == 3):
            plt.xlim(left=0, right=180)
        else:
            plt.xlim(left=-90, right=270)
        plt.savefig(self.folder + "plots/histogram_orientations.png", bbox_inches = 'tight')
        plt.close()


    def save(self, cleaned=False, vtk=False):
        '''
        Saves the structure to file.
        If cleaned=True it uses the cleaned version of the structure (i.e., the spanning cluster only).
        If vtk=True the structure is saved as a VTK-file (this requires PyVista and works only for 3D structures),
        otherwise it is saved as a CSV-file with columns representing the x, y and z coordinates where
        x = 0, ..., shape[0]-1, y = 0, ..., shape[1]-1 and z = 0, ..., shape[2]-1.
        '''
        print("Saving structure to file...")

        if cleaned:
            binary = self.binary_array_cleaned
            name = "_cleaned"
        else:
            binary = self.binary_array
            name = ""

        # np.save(self.folder + "binary" + name + ".npy", binary)
        pos = np.where(binary==1)
        if len(self.shape) == 2:
            point_cloud = np.column_stack((pos[0], pos[1]))
        elif len(self.shape) == 3:
            point_cloud = np.column_stack((pos[0], pos[1], pos[2]))

        if len(self.shape)==3 and vtk:
            import pyvista
            pdata = pyvista.PolyData(point_cloud, force_float=False)
            pc = pdata.glyph(scale=False, geom=pyvista.Cube(), orient=False)
            pc.save(self.folder + "xyz.vtk")
        else:
            np.savetxt(self.folder + "xyz" + name + ".csv", point_cloud, fmt='%d', delimiter=',')


    def plot(self, slicelist=[0], cleaned=False, figsize=(13, 8), dpi=100, fontsize=22):
        '''
        Saves plots of the structure.
        If 3D structure, user-defined slices are saved.
        The argument slicelist is only used for 3D structures and is a list of x-, y- and z-indices of slices to be plotted.
        If cleaned=True it uses the cleaned version of the structure (i.e., the spanning cluster only).
        The other optional parameters figsize, dpi and fontsize relates to the size and resolution of the plot.
        '''
        print("Plotting structure...")

        plt.rcParams.update({'font.size': fontsize})
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')

        if cleaned:
            binary = self.binary_array_cleaned
            name = "_cleaned"
        else:
            binary = self.binary_array
            name = ""

        if (len(self.shape) == 3):

            if True in (slicelist >= np.min(self.shape)):
                print("    ERROR: slicelist elements must be smaller than min(shape)")
                return

            for i in slicelist:
                plt.figure(figsize=figsize, dpi=dpi)
                plt.imshow(binary[:,:,i].T, cmap='Greys', interpolation='nearest', origin='lower')
                plt.xlabel(r'$x$')
                plt.ylabel(r'$y$')
                plt.savefig(self.folder + "plots/plot_slice_z_n" + str(i) + name + ".png", bbox_inches = 'tight')
                plt.close()

                plt.figure(figsize=figsize, dpi=dpi)
                plt.imshow(binary[:,i,:].T, cmap='Greys', interpolation='nearest', origin='lower')
                plt.xlabel(r'$x$')
                plt.ylabel(r'$z$')
                plt.savefig(self.folder + "plots/plot_slice_y_n" + str(i) + name + ".png", bbox_inches = 'tight')
                plt.close()

                plt.figure(figsize=figsize, dpi=dpi)
                plt.imshow(binary[i,:,:].T, cmap='Greys', interpolation='nearest', origin='lower')
                plt.xlabel(r'$y$')
                plt.ylabel(r'$z$')
                plt.savefig(self.folder + "plots/plot_slice_x_n" + str(i) + name + ".png", bbox_inches = 'tight')
                plt.close()

        else:
            plt.figure(figsize=figsize, dpi=dpi)
            plt.imshow(binary[:,:].T, cmap='Greys', interpolation='nearest', origin='lower')
            plt.xlabel(r'$x$')
            plt.ylabel(r'$y$')
            plt.savefig(self.folder + "plots/plot" + name + ".png", bbox_inches = 'tight')
            plt.close()



    def visualize_3d(self, cleaned=False, pyvista=True, from_file=False):
        '''
        Visualizes the full 3D structure.
        If cleaned=True it uses the cleaned version of the structure (i.e., the spanning cluster only).
        If pyvista=True (recommended) PyVista is used for the visualization.
        If from_file=True the structure is read from the CSV-file in the directory.
        The other optional parameters figsize, dpi and fontsize relates to the size and resolution of the plot.
        '''
        print("Visualizing 3D structure...")

        if (len(self.shape) != 3):
            print("    ERROR: visualize3D is only for 3D structures")
            return

        if from_file:
            if cleaned:
                name = "_cleaned"
            else:
                name = ""
            point_cloud = np.loadtxt(self.folder + "xyz" + name + ".csv", delimiter=',')
        else:
            if cleaned:
                binary = self.binary_array_cleaned
            else:
                binary = self.binary_array
            pos = np.where(binary==1)
            point_cloud = np.column_stack((pos[0], pos[1], pos[2]))

        if pyvista:
            import pyvista
            pdata = pyvista.PolyData(point_cloud, force_float=False)
            pc = pdata.glyph(scale=False, geom=pyvista.Cube(), orient=False)
            pc.plot()
        else:
            from mpl_toolkits.mplot3d import Axes3D

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.set_aspect('equal')
            ax.scatter(point_cloud[:,0], point_cloud[:,1], point_cloud[:,2], c='black')
            # ax.voxels(m, facecolors='k', edgecolor='k') # this is very slow!
            set_axis_equal(ax)
            plt.show()


    def cleaning(self):
        '''
        Cleans the microstructure by removing clusters that do not belong to any spanning clusters.
        Uses the so-called Burning Method, see clean_structure.py for more details.
        '''
        print("Cleaning structure...")
        self.binary_array_cleaned = clean_structure(self.binary_array)

        phi_observed_cleaned = self.binary_array_cleaned.sum() / self.num_mesh_points
        print("    Observed solid volume fraction (cleaned): " + str(phi_observed_cleaned) )
        num2txt(self.folder, "phi_observed_cleaned", round(phi_observed_cleaned,4))


    def find_shortest_path(self, cleaned=False):
        '''
        Calculates the shortest path of the spanning cluster by applying the Burning Method from all sides.
        See find_shortest_path.py for more details.
        If cleaned=True it uses the cleaned version of the structure (i.e., the spanning cluster only).
        '''
        print("Calculating shortest path through structure...")

        if cleaned:
            binary = self.binary_array_cleaned
            name = "_cleaned"
        else:
            binary = self.binary_array
            name = ""

        shortest_path = find_shortest_path(binary)
        np.savetxt(self.folder + "shortest_path" + name + ".txt", shortest_path, fmt='%d')


    def autocorrelation(self, cleaned=False, one_dim=False, figsize=(13, 8), dpi=100, fontsize=22):
        '''
        Calculates and plots the (normalized) autocorrelaton function, aka two-point correlation function.
        The specific surface area per unit solid volume (SSA) is also estimated from this function.
        If cleaned=True it uses the cleaned version of the structure (i.e., the spanning cluster only).
        If one_dim=True the one-dimensional autocorrelaton in all (2 or 3) directons is also calculated.
        The other optional parameters figsize, dpi and fontsize relates to the size and resolution of the plot.
        '''
        print("Calculating autocorrelation...")

        if cleaned:
            binary = self.binary_array_cleaned
            name = "_cleaned"
        else:
            binary = self.binary_array
            name = ""

        phi_observed = np.sum(binary) / self.num_mesh_points

        autocorr_r, dist_autocorr_r = autocorrelation_fft(binary)
        autocorr_r *= phi_observed

        ### surface area per unit total (solid+void) volume
        SSA = -4 * ( autocorr_r[1] - autocorr_r[0] ) # Note that dist_autocorr_r[1] - dist_autocorr_r[0] = 1 - 0 = 1

        ### convert from "surface area per unit total volume" to "surface area per unit solid volume"
        SSA /= phi_observed

        print("    SSA: ", round(SSA,4))
        num2txt(self.folder, "SSA" + name, SSA)

        ### normalize autocorrelation
        autocorr_r   = (autocorr_r  - phi_observed**2) / ( phi_observed - phi_observed**2 )

        if (one_dim):
            r_max = int(min(self.shape)/3)
            autocorr_xyz = autocorrelation_xyz(binary, r_max, self.procs)
            dist_autocorr_x = np.arange(0, len(autocorr_xyz[0]))

            ### normalize autocorrelation
            autocorr_xyz  = (autocorr_xyz - phi_observed**2) / ( phi_observed - phi_observed**2 )

        lambda_mean = (self.shape[0]-1) / self.m_mean

        plt.rcParams.update({'font.size': fontsize})
        plt.rc('text', usetex=False)
        plt.rc('font', family='serif')
        plt.figure(figsize=figsize, dpi=dpi)
        plt.plot( dist_autocorr_r / lambda_mean, autocorr_r,   'k-',  label=r'FFT'   )
        if (one_dim):
            label = [r'$x$-dir', r'$y$-dir', r'$z$-dir']
            for i in range(0, len(self.shape)):
                plt.plot( dist_autocorr_x / lambda_mean, autocorr_xyz[i], '*-', label=label[i] )
        plt.legend()
        plt.xlabel(r'$r/\langle \lambda \rangle$')
        plt.ylabel(r'$g_{norm}$')
        plt.xlim(0, min(2, np.max(dist_autocorr_r / lambda_mean)))
        plt.savefig(self.folder + "plots/autocorrelation" + name + ".png", bbox_inches = 'tight')


    def analytic_ssa(self):
        '''
        Calculates the analytic prediction of the specific surface area per unit solid volume (SSA).
        Such analytic expression is only available for isotropic single-cut structures based on the gamma-distribution.
        '''
        print("Calculating SSA analytically...")
        q_mean  = 2*np.pi * self.m_mean / (self.shape[0]-1)
        delta_q = 2*np.pi * self.m_std / (self.shape[0]-1)
        b = (q_mean/delta_q)**2 - 1
        if len(self.shape) == 3 and self.single_cut and ((self.anitype != 'v' and self.anitype !='h') or self.a == 1) and (self.distr == 'g'):
            xi = erfinv(1 - 2 * self.phi)
            SSA = 2 * q_mean / ( np.pi * np.sqrt(3.0)) * np.sqrt( (b+2.0) / (b+1.0) ) * np.exp(-2*xi**2)
            SSA /= self.phi
            print("    Analytic SSA: ", round(SSA,4))
            num2txt(self.folder, "SSA_analytic", SSA)
        else:
            print("    Analytic expression for SSA is not available with these parameters.")
