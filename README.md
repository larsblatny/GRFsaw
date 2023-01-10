# GRFsaw

GRFsaw is a simple software for artificially constructing binary microstructures (in two or three spatial dimensions) with user-defined
* porosity (solid volume fraction)
* heterogeinity ("grain" size distribution)
* anisotropy (preferred "grain" elongation)

It is based on thresholding Gaussian random fields (GRFs).

Several post-processing schemes are available for
* computing the angular-averaged two-point correlation function using Fast Fourier Transforms
* computing the one-dimensional two-point correlation function in different directions
* finding the spanning (percolating) cluster(s) using the Burning Method
* finding the shortest path through the structure
* plotting and visualizing in 2D and 3D
* saving structures as CSV- or VTK-files

## Dependencies
This code relies only on the standard libraries [numpy](https://numpy.org/), [scipy](https://scipy.org/), [matplotlib](https://matplotlib.org/) and [multiprocessing](https://docs.python.org/3/library/multiprocessing.html)

Use of [PyVista](https://docs.pyvista.org/) is optional, but needed for visualizing 3D microstructures and saving structures as VTK-files.

Tested on Python >= 3.6.

## Usage

```python
from grfsaw.microstructure import Microstructure
```
### Example minimal usage
```python
shape = [200, 200]    
solid_fraction = 0.65
m = Microstructure(shape, solid_fraction)
m.create()
m.plot()
```

### List of user options
These are all the options the user may specify for the desired microstructure.   
Only two parameters are required, ```shape``` and ```phi```, the rest are given their default values below:
```python
shape = <REQUIRED> # list of 2 (if 2D structure) or 3 (if 3D) integers
#   the size of the structure as number of mesh points in each direction

phi = <REQUIRED> # float
#   desired solid volume fraction = 1 - porosity

folder = "output/" # str
#   name of directory to save data and plots

single_cut = True # boolean
#     decides if generate single-cut or double-cut structure

N = 10000 # int
#   number of waves to use, the more the better but slower

distr = 'n' # char
#     type of distribution to use for wavevector magnitude sampling,
#     'n' for normal or 'g' for gamma distribution

m_mean = 3 # float
#   average number of microstructural elements per length in the x-direction
#   = a measure of the grain size of the structure

m_std = 0.001 # float
#   the standard deviation of m

anitype = 'i' # char
#   type of preferred orientation (type of anisotropy)
#   'h' for horizontal or 'v' for vertical, anything else means fully isotropic

a = 1 # float
#   measure of the anisotropy level between 0 and 1, meaning isotropic (default 1)

seed = numpy.random.randint(1, 1e5) # int
#   seed for RNG

procs = multiprocessing.cpu_count() # int
#   number of cpus to use for the computations (default is all cpus on computer)


m = Microstructure(shape,
                   phi,
                   folder,
                   single_cut,
                   N,
                   distr,
                   m_mean,
                   m_std,
                   anitype,
                   a,
                   seed,
                   procs)
```

### Example: 2D microstructures

![Microstructures2D](figures/pic_2d.png?raw=true "micro2d")

From left to right (all using `N = 10000` and  `distr = 'g'`):   
1) `phi = 0.63`,  `single_cut = True`,  `m_mean = 13`, `m_std = 1.8`,  `anitype = 'i'`   
2) `phi = 0.63`,  `single_cut = True`,  `m_mean = 13`, `m_std = 7.5`,  `anitype = 'i'`     
3) `phi = 0.63`,  `single_cut = True`,  `m_mean = 13`, `m_std = 1.8`,  `anitype = 'h'`, `a = 0.6`    
4) `phi = 0.35`,  `single_cut = False`, `m_mean = 9`,  `m_std = 1.3`,  `anitype = 'i'`   
5) `phi = 0.35`,  `single_cut = False`, `m_mean = 9`,  `m_std = 5.2`,  `anitype = 'i'`

### Example: 3D microstructures

![Microstructures3D](figures/pic_3d.png?raw=true "micro3d")
Both structures: `phi = 0.42`,  `single_cut = True`,  `N = 10000`,  `distr = 'g'`,  `anitype = 'i'` and  `m_mean = 9`  
Left: `m_std = 1.3`, Right: `m_std = 5.2`

### List of methods
In all the below methods, the argument ```cleaning``` (default ```False```) refers to if you are using the cleaned version of the structure or not. In the cleaned version, all parts of the structure that do not belong to a spanning cluster have been removed.
In addition, the arguments ```figsize```, ```dpi``` and ```fontsize``` refer to size, resolution and font size, respectively, of the plots to be generated.

```python
m.create()
#   creates the structures

m.save_parameters()
#   saves setup parameters to file for future reference

m.save(cleaned=False, vtk=False)
#   saves the structure as CSV-file or as VTK-file (if vtk=True)

m.plot(slicelist=[0], cleaned=False, figsize=(13, 8), dpi=100, fontsize=22)
#   plots the structure, or, if 3D structure, a slices of the structuregiven by the arg slicelist

m.visualize_3d(cleaned=False, pyvista=True, from_file=False)
#   visualizes the 3D structure using PyVista (if pyvista=True, recommended) or matplotlib

m.plot_distributions(figsize=(13, 8), dpi=100, fontsize=22)
#   plots histograms of the samples of the wave lengths and direction

m.autocorrelation(cleaned=False, one_dim=False, figsize=(13, 8), dpi=100, fontsize=22)
#   computes and plots the autocorrelation (two-point correlation) of the structure using FFT
#   if one-dim=True, it also computes the 1D autocorrelation in the x-, y- and z-directions

m.find_shortest_path(cleaned=False)
#   calculates the shortest path of the spanning cluster

m.cleaning()
#   removes parts of the structure that do not belong to any side-to-side spanning clusters

m.analytic_ssa()
#   calculates the analytic specific surface area per unit solid volume (SSA)
```

The functions are further documented here: https://larsblatny.github.io/GRFsaw/

## Example: Cleaning
A microstructure before (left) and after (right) applying `cleaning()`. We see that stand-alone clusters are removed.
![MicrostructuresCleaned](figures/cleaned.png?raw=true "Cleaned")

## Example: Autocorrelation
A plot of the two-point correlation function (aka spatial autocorrelation function) produced by `autocorrelation(one_dim=True)`. The angular-averaged two-point correlation as obtained from FFT is plotted in black. This was a highly anisotropic structure, noticed by the differing 1D correlations in the x- and y directions.
![Autocorrelation](figures/autocorr.png?raw=true "Correlation")


## Contributing
Pull requests are welcome. Or feel free to open an issue to discuss what can be changed.

## License
This code is licensed under the [MIT License](https://choosealicense.com/licenses/mit/)
