======
Computing maps of the torque distribution/density
======
.. sectnum::

   
Torque distribution map for a single snapshot
-----

First: Compute gravitational accelerations
~~~~~

Sometimes, your code output might have the gravitational
acceleration of cells already stored in. However, if you need the
back-reaction torque acting on, for example, a central binary, unless
the binary elements are a live part of the output, you might have to
recompute the accelerations due to the members of the binary
*separately*.


Accelerations
....

.. code:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import disk_data_analysis.circumbinary as dda

   snap = dda.get_snapshot_data('./data/snap_',0,['POS','VELX','VELY','RHO','ACCE'])

	  
In order to compute the accelerations due to two objects orbiting each other in a Keplerian fashion, we need to compute their positions as a function of time (this only if your binary elements ARE NOT present as particles in the simulation snapshot).

.. code:: python

   from disk_data_analysis.orbital import orbit_in_time
   
   # Compute forces

   # First, 
   time = snap.header.time
   eb = 0.0
   x, y, _, _ = orbit_in_time(time + np.pi, eb)
   accel = dda.compute_external_gravforce(snap,XYZ=[0]



.. image:: ./doc_images/accretion_regions.png


Separation of torques in intensity maps
....
	   
Recap: intensity images
''''

Recall that we can visualize previously-generated image data.


.. code:: python

   import matplotlib.pyplot as plt
   from disk_data_analysis.plotting import plot_slice, ImageData

   image = ImageData('../data/density_field_000')
   fig = plt.figure(figsize = (6.0,6.0))
   ax = fig.add_subplot(111)
   ax = plot_slice(ax,image)
   ax.set_xlabel(r'$x$',size=18)
   ax.set_ylabel(r'$y$',size=18)
   ax.set_aspect('equal')


.. image:: ./doc_images/density_field.png


Three color maps in one figure
''''

To do this, we need to do the spatial separation as above but not in cell-center coordinates, but in
pixel coordinates

.. code:: python


   pixelxcoord, pixelycoord = np.meshgrid(np.arange(image.pixelx)*1.0/image.pixelx,\
	                                  np.arange(image.pixely)*1.0/image.pixely)

   Lx, Ly = extent[1] - extent[0], extent[3] - extent[2]
   pixelxcoord, pixelycoord = (pixelxcoord.T - 0.5) * Lx + 0.5 * snap.header.boxsize,\
                               (pixelycoord.T - 0.5) * Ly + 0.5 * snap.header.boxsize
    
    
   dxpixel1, dypixel1 = pixelxcoord - pos1[0], pixelycoord - pos1[1]
   dxpixel2, dypixel2 = pixelxcoord - pos2[0], pixelycoord - pos2[1]
   dxpixel, dypixel = pixelxcoord - 0.5 * snap.header.boxsize, pixelycoord - 0.5 * snap.header.boxsize
   
   drpixel1 = np.sqrt(dxpixel1**2 + dypixel1**2)
   drpixel2 = np.sqrt(dxpixel2**2 + dypixel2**2)
   drpixel = np.sqrt(dxpixel**2 + dypixel**2)
   
   #we do the same region separation as before
   csd_region = (drpixel1 < csd_trunc) | (drpixel2 < csd_trunc)
   cbd_region = drpixel > cbd_trunc
   str_region = ((drpixel <= cbd_trunc) ) &\
   ((drpixel1 >= csd_trunc) & (drpixel2 >= csd_trunc))
   
   image_csd = ImageData()
   image_cbd = ImageData()
   image_str = ImageData()
   
   image_csd.read('../data/density_field_000')
   image_cbd.read('../data/density_field_000')
   image_str.read('../data/density_field_000')
   
   image_cbd.data[np.invert(cbd_region)] = np.nan    
   image_csd.data[np.invert(csd_region)] = np.nan
   image_str.data[np.invert(str_region)] = np.nan
   
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax = plot_slice(ax,image_cbd,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,
   cmap=cm.get_cmap('Blues'))
   ax = plot_slice(ax,image_csd,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,
   cmap=cm.get_cmap('Purples'))
   ax = plot_slice(ax,image_str,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,
   cmap=cm.get_cmap('Reds'))
   ax.set_xlabel(r'$x$',size=18)
   ax.set_ylabel(r'$y$',size=18)
   ax.set_aspect('equal')
   
   plt.show()


.. image:: ./doc_images/density_field_colormaps.png
