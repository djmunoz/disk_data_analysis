import numpy as np
from scipy.integrate import cumtrapz

"""
Set of routines to compute angular momentum transfer rates as defined
in Miranda, Munoz & Lai (2017).

"""

def compute_angular_momentum_transfer(snapnum):

    snap = dda.get_snapshot_data('./data/snap_',snapnum,['POS','VEL','RHO'])
    snap.add_snapshot_data('VEL')
    
    torque_adv,torque_visc,torque_grav = angular_momentum_transfer(snap,grid)


def angular_momentum_transfer(snapshot,grid):
    """
    Compute the angular momentum transfer as a function of
    radius from a simulation snapshot and a prescribed grid
    data structure containing the coordinates of grid points
    where the torque balance should be evaluated.

    """
    
    NR, Nphi = grid.NR, grid.Nphi
    X, Y = grid.X, grid.Y
    radii = grid.R

    # Compute the re-gridded quantities first
    rho_interp,velx_interp,vely_interp  = dda.disk_interpolate_primitive_quantities(snapshot,[X,Y],quantities=['RHO','VELX','VELY'])

    # Make sure we have the velocity gradients
    gradvx,gradvy =  dda.disk_compute_derivative_quantities(snapshot,[X,Y],quantities=['VELX','VELY'])

    # And also make sure we have the gravitational acceleration
    
    
    # Compute the advection transfer rate
    vrvphi = x * y * (vy * vy - vx * vx) +  vx * vx
    Tadv = - 2 * np.pi * ((rho * x * y * vx).rshape(NR,Nphi)).mean(axis=0)
    
    # Compute the radial torque density...
    dTdR = ((radii * rho * dPhi_dphi).rshape(NR,Nphi)).mean(axis=0)
    #... and integrate it into a torque
    Tgrav = cumptrapz(dTdr,x=R)






    return Tadv, Tvisc, Tgrav
