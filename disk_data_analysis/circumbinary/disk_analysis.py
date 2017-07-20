from disk_hdf5 import readsnapHDF5 as rs
import numpy as np
import sys
from scipy.integrate import cumtrapz, trapz
from disk_interpolate_primitive import *
from disk_simulation_data import *

def disk_compute_radial_balance(snapshot, griddata, nu_gridded, Rmin = None, Rmax = None):

    if (Rmin is None): Rmin  = griddata.R.min()
    if (Rmax is None): Rmax  = griddata.R.max()
    
    # Make sure we have all the necessary gradient information
    gradientvx = compute_snapshot_gradient(snapshot,'VELX')
    gradientvy = compute_snapshot_gradient(snapshot,'VELY')

    # Add the gradients
    snapshot.add_data(gradientvx,'GRVX')
    snapshot.add_data(gradientvy,'GRVY')
    snapshot.add_data(snapshot.gas.ACCE[:,0:2],'GRPHI')


    #interpolate primitive quantities
    rho_interp = disk_interpolate_primitive_quantities(snapshot,[griddata.X,griddata.Y],\
                                                           quantities=['RHO'],method = 'linear')[0]
    
    vx_interp = disk_interpolate_primitive_quantities(snapshot,[griddata.X,griddata.Y],\
                                                          quantities=['VELX'],method = 'linear')[0]
    
    vy_interp = disk_interpolate_primitive_quantities(snapshot,[griddata.X,griddata.Y],\
                                                          quantities=['VELY'],method = 'linear')[0]

    # interpolate derivatives
    gradvx_interp = disk_interpolate_gradient_quantities(snapshot,[griddata.X,griddata.Y],\
                                                        quantities=['GRVX'],method = 'nearest')[0]

    gradvy_interp = disk_interpolate_gradient_quantities(snapshot,[griddata.X,griddata.Y],\
                                                             quantities=['GRVX'],method = 'nearest')[0]
    
    gradphi_interp = disk_interpolate_gradient_quantities(snapshot,[griddata.X,griddata.Y],\
                                                          quantities=['GRPHI'],method = 'nearest')[0]


    
    gridR = griddata.R.mean(axis=0)
    gridX, gridY = griddata.X - snapshot.header.boxsize * 0.5, griddata.Y  -  snapshot.header.boxsize * 0.5

    # Mass transfer rate
    mdot = -2 * np.pi * (rho_interp * (gridX * vx_interp + gridY * vy_interp)).mean(axis=0)

    # Angular momentum transfer rate
    jdot_adv = -2 * np.pi * (rho_interp * (gridX * gridY * (vy_interp**2 - vx_interp**2) +\
                                          vx_interp * vy_interp * (gridX**2 - gridY**2))).mean(axis=0)
    
    jdot_visc = (-2 * np.pi * nu_gridded * rho_interp * \
                 (2 * gridX * gridY * (gradvy_interp[1] - gradvx_interp[0]) + \
                  (gridX**2 - gridY**2) * (gradvx_interp[1] + gradvy_interp[0]))).mean(axis=0)


    dTgravdR = -2 * np.pi * (griddata.R * rho_interp * (gradphi_interp[1] * gridX -\
                                                    gradphi_interp[0] * gridY)).mean(axis = 0)
    dTgravdR[(gridR > Rmax) | (gridR < Rmin)] = 0.0

    Tgrav = np.asarray([trapz(dTgravdR[gridR > R],x=gridR[gridR > R]) for R in gridR])


    return mdot, jdot_adv, jdot_visc, Tgrav
