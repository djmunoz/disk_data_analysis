Computing angular momentum balance
======



.. code:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import disk_data_analysis.circumbinary as dda

   snap = dda.get_snapshot_data('./data/snap_',0,['POS','VELX','VELY','RHO','ACCE','R'])

   # if we do not have gradient data, we can compute it
   gradientvx = dda.compute_snapshot_gradient(snap,'VELX')
   gradientvy = dda.compute_snapshot_gradient(snap,'VELY')

   snap.add_data(gradientvx,'GRVX')
   snap.add_data(gradientvy,'GRVY')

As in other examples, in order to get radial profiles, it is useful to create a regularly spaced
polar grid

.. code:: python
	  
   # as in other cases, map onto a regular grid
   Rmin, Rmax = 1.0, 80.0
   NR, Nphi = 200, 400
   grid = dda.grid_polar(NR = NR, Nphi = Nphi,Rmin=1.0,Rmax= 80.0,scale='log')
   grid.X, grid.Y = grid.X + snap.header.boxsize * 0.5, grid.Y  +  snap.header.boxsize * 0.5


Onto this grid, we want to map |inlineq1| to compute the advective angular momentum transfer rate:

.. |inlineq1| image:: ./doc_images/inline_eq1.png
   :align: middle
		      
.. image:: ./doc_images/equation1.png


	   
Similarly, we want to map |inlineq2| onto the grid to compute

.. |inlineq2| image:: ./doc_images/inline_eq2.png
   :align: middle

.. image:: ./doc_images/equation2.png	      


	      
.. code:: python
	  
   #interpolating ...
   rho_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],\
	                                                  quantities=['RHO'],method = 'linear')[0]

   vx_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],\
	                                                 quantities=['VELX'],method = 'linear')[0]

   vy_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],\
	                                                 quantities=['VELY'],method = 'linear')[0]
   

Now, let us pay specially attention on how we carry out the mapping for higher-order primitive
quantities such as the velocity gradients.


.. code:: python
	  
   #interpolating ...
   dvxdx_interp = dda.disk_interpolate_gradient_quantities(snap,[grid.X,grid.Y],\
	                                                   quantities=['GRVX'],method = 'linear')[0]





   quantities=['VELX'],method = 'linear')[0]

   vy_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],\
	                                                 quantities=['VELY'],method = 'linear')[0]
   
   
   grid = dda.grid_cartesian(Xmin=-80.0,Xmax=80.0,Ymin=-80.0,Ymax=80.0,NX=512,NY=512,mask= '(R < 1.0) | (R > 80.0)')
   grid.X, grid.Y =  grid.X + snap.header.boxsize * 0.5, grid.Y + snap.header.boxsize * 0.5
   
   gradvx = dda.compute_gradient_on_grid(snap.gas.POS[:,0], snap.gas.POS[:,1],\
                                         snap.gas.VEL[:,0],grid)
   
   snap.add_snapshot_gradient('VEL','GRAV')
