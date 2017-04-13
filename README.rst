disk_data_analysis - Analysis and Plotting Tools for HDF5 simulation data of meshless/moving-mesh hydrodynamical disks 
==================================================
.. sectnum::

.. class:: no-web
           
   .. image:: example_figures/disk_images.png
      :height: 100px
      :width: 400 px

Overview
--------

A Python package!

Installation
--------

You need to have git installed. In addition, you need the NumPy and SciPy Python packages.

.. code::
   
   git clone https://github.com/djmunoz/disk_data_analysis.git

   cd disk_data_analysis
   
   sudo python setup.py install

That is all!


Examples
--------


Reading-in and plotting data
~~~~~~~

First of all, load the package and check the functions/attributes in it:

.. code:: python
	  
   import disk_data_analysis.circumbinary as dda

   dir(dda)
   
The most basic function in that list is :code:`get_snapshot_data`,
which can be used to give

.. code:: python
	  
   snap = dda.get_snapshot_data('./data/snap_',0,['POS','VEL'])

which will only contain the requested quantities :code:`POS`, (positions)
and :code:`VEL` (velocities). You can double-check the attributes of
the  :code:`snap` data structure by doing

.. code:: python
	  
   dir(snap.gas)

   
You can see the distribution of particles/mesh generating points by plotting
the positions

.. code:: python

   import matplotlib.pyplot as plt

   x = snap.gas.POS[:,0]
   y = snap.gas.POS[:,1]
   box = snap.header.boxsize
   plt.plot(x,y,'b.',ms=2.4)
   plt.xlim(0.5 * box - 2.5, 0.5 * box + 2.5)
   plt.ylim(0.5 * box - 2.5, 0.5 * box + 2.5)
   plt.xlabel(r'$x$',size=18)
   plt.ylabel(r'$y$',size=18)
   plt.show()
   
.. image:: example_figures/figure_points.png
   :align: right
   :scale: 50
   :height: 40px
   :width: 40 px
   
Computing radial profiles
~~~~~~~


Mapping onto polar grids
~~~~~~~

Often, SPH and moving-mesh simulations of disks will be compared to simulation
results obtained with polar-grid codes. In particular, comparison of
azimuthally-averaged quantities is common practice. While azimuthal averaging
is trivial in polar grids, it requires some tinkering when computational cells
(or particles) are not placed in regular intervals in radius. One way around
this is to remap the fluid's primitive variables into a structured grid of
points.


.. code:: python


   import numpy as np

   # We need the density and positions
   snap = dda.get_snapshot_data('./data/snap_',0,['POS','RHO'])
   x, y, dens = snap.gas.pos[:,0], snap.gas.pos[:,1], snap.gas.rho

   # get a sense of the dynamical range in radius in the simulation
   Rmin, Rmax = 1.0, 80.0
   NR, Nphi = 200, 300
   R, phi = np.meshgrid(np.logspace(np.log10(Rmin),np.log10(Rmax),NR),\
	                np.linspace(0,2*np.pi,Nphi))
   X, Y = R * np.cos(phi), R * np.sin(phi)
   
   # interpolate Z values on defined grid
   Z = griddata(np.vstack((x.flatten(),y.flatten())).T, \
	        np.vstack(z.flatten()),(X,Y),method='linear').reshape(X.shape)
			
   
   density_gridded = da.

Perhaps you would rather use an unevenly sampled grid loosely based
on the actual positioning of the cells/particles

.. code:: python
	  
   polar_grid =
   


Displaying 2-D fields without pre-computed image data
~~~~~~~
