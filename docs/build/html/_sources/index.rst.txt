.. GRFsaw documentation master file, created by
   sphinx-quickstart on Wed Dec  7 15:04:48 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation of GRFsaw!
==========================================

GRFsaw is a Python code for artificially constructing bicontinuous microstructures (in two or three spatial dimensions) with user-defined

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

For it's usage and dependencies, please see the README.md file.

For the mathematical background, please see the paper.pdf file in the paper directory.

GRFsaw is licensed under the MIT License.

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
