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


Reading-in data
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
and :code:`VEL` (velocities)

You can see the distribution of particles/mesh generating points by plotting
the positions

.. code:: python

   import matplotlib.pyplot as plt

   x = snap.gas.pos[:,0]
   y = snap.gas.pos[:,1]
   box = snap.header.boxsize
   plt.plot(x,y,'b.',ms=1.0)
   plt.xlim(0.5 * box - 2.5, 0.5 * box + 2.5)
   plt.ylim(0.5 * box - 2.5, 0.5 * box + 2.5)
   plt.xlabel(r'$x$',size=18)
   plt.ylabel(r'$y$',size=18)
   plt.show()
   

   
Computing radial profiles
~~~~~~~


Mapping onto polar grids
~~~~~~~


.. code:: python
	  
   polar_grid =
   
   density_gridded = da.

Perhaps you would rather use an unevenly sampled grid loosely based
on the actual positioning of the cells/particles

.. code:: python
	  
   polar_grid =
   


Displaying 2-D fields without pre-computed image data
~~~~~~~
