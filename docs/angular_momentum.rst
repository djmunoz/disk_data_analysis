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
   gradvx_interp = dda.disk_interpolate_gradient_quantities(snap,[grid.X,grid.Y],\
	                                                   quantities=['GRVX'],method = 'nearest')[0]

							   
   gradvy_interp = dda.disk_interpolate_gradient_quantities(snap,[grid.X,grid.Y],\
	                                                   quantities=['GRVX'],method = 'nearest')[0]

With this quantities re-mapped, we can use Equation(1) to compute the advective angular momentum
transfer rate:

.. code:: python
	  
   grid.X, grid.Y = grid.X - snap.header.boxsize * 0.5, grid.Y  -  snap.header.boxsize * 0.5
   # compute the advective angular momentum transfer rate
   jdot_adv = -2 * np.pi * rho_interp * (grid.X * grid.Y * (vy_interp**2 - vx_interp**2) +\
                                         vx_interp * vy_interp * (grid.X**2 - grid.Y**2))

   # average out the azimuthal axis
   jdot_adv = jdot_adv.mean(axis=0)

To compute the viscous transfer rate, we need one more element: the kinematic viscosity coefficient
*nu* as a function of radius on the grid.

.. image:: ./doc_images/equation3.png

.. code:: python

   alpha = 0.1
   h0 = 0.1
   GM = 1.0
   
   def nu(R):
	  return alpha * h0**2 * np.sqrt(GM) * R**(0.5)

   nu_grid = nu(grid.R)

.. code:: python
   
   # Similarly, compute the viscous angular momentum transfer rate
   jdot_visc = (-2 * np.pi * nu_grid * rho_interp * \
	        (2 * grid.X * grid.Y * (gradvy_interp[1] - gradvx_interp[0]) + \
		 (grid.X**2 - grid.Y**2) * (gradvx_interp[1] + gradvy_interp[0]))).mean(axis=0)


It is useful to normalize the angular momentum flux in units of:

.. image:: ./doc_images/equation4.png

.. code:: python


   mdot = -2 * np.pi * (rho_interp * (grid.X * vx_interp + grid.Y * vy_interp)).mean(axis=0)
   # if you do not know mdot0 from your simulation setup, it can be re-computed as 
   mdot0 = mdot[(grid.R.mean(axis = 0) < 62) & (grid.R.mean(axis = 0) > 50)].mean()

   jdotnorm = mdot0
   
   # and plot it
   plt.plot(grid.R.mean(axis=0),jdot_adv)
   plt.xlim(0,15)
   plt.xlabel(r'$R$')
   plt.ylabel(r'$\dot{J}_{\rm adv}$')
   plt.show()
					 

Of course, there is still one more term in the angular momentum balance equation, and that is
the external gravitational torque:


.. image:: ./doc_images/equation5.png


We can treat the gravitational acceleration as a gradient evaluated at the center of a cell:


.. code:: python

   snap.add_data(snap.gas.ACCE[:,0:2],'GRPHI')

   # and interpolate using the 'nearest' method
   gradphi_interp = dda.disk_interpolate_gradient_quantities(snap,[grid.X,grid.Y],\
	                                                     quantities=['GRPHI'],method = 'nearest')[0]

Then the gravitational torque density and the integrated gravitational torque are:
   
.. code:: python


   dTgravdR = -2 * np.pi * (grid.R * rho_interp * (gradphi_interp[1] * grid.X -\
                                                   gradphi_interp[0] * grid.Y)).mean(axis = 0)

   # Before integrating, make sure anomalous values are not taken into account
   Rmax = 70
   dTgravdR[grid.R.mean(axis = 0) > Rmax] = 0.0
						   
   # now we integrate
   from scipy.integrate import cumtrapz, quad

   # First, the slow way
   Tgrav_slow = [quad(dTgravdR,R,Rmax) for i in grid.R.mean(axis = 0)]

   plt.plot(grid.R.mean(axis=0),Tgrav_slow)
   plt.show()
   
   Tgrav = cumtrapz(dTgravdR[::-1],x=-grid.R.mean(axis=0)[::-1],initial=0)


   
	  
   vy_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],\
	                                                 quantities=['VELY'],method = 'linear')[0]
   
   
   grid = dda.grid_cartesian(Xmin=-80.0,Xmax=80.0,Ymin=-80.0,Ymax=80.0,NX=512,NY=512,mask= '(R < 1.0) | (R > 80.0)')
   grid.X, grid.Y =  grid.X + snap.header.boxsize * 0.5, grid.Y + snap.header.boxsize * 0.5
   
   gradvx = dda.compute_gradient_on_grid(snap.gas.POS[:,0], snap.gas.POS[:,1],\
                                         snap.gas.VEL[:,0],grid)
   
   snap.add_snapshot_gradient('VEL','GRAV')
