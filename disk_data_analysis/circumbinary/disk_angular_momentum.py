import numpy as np
from scipy.integrate import cumtrapz
from disk_interpolate_primitive import *
from disk_simulation_data import *

"""
Set of routines to compute angular momentum transfer rates as defined
in Miranda, Munoz & Lai (2017).

"""

def compute_angular_momentum_transfer(snapshot, Rmin, Rmax, NR = None, Nphi = None, alpha=0.1, h0=0.1):
    """
    Compute the angular momentum transfer as a function of
    radius from a simulation snapshot and a prescribed grid
    data structure containing the coordinates of grid points
    where the torque balance should be evaluated.

    """
    

    if (NR is None):
        NR = 1024
    if (Nphi is None):
        Nphi = int(2 * np.pi / (10**((np.log10(Rmax) - np.log10(Rmin))/NR) - 1))

    # Create a polar grid
    grid = grid_polar(NR = NR, Nphi = Nphi, Rmin= Rmin,Rmax = Rmax,scale='log')

    mdot = mass_advection(snapshot,grid)
    
    torque_adv = angular_momentum_advection(snapshot,grid)

    torque_visc = angular_momentum_viscosity(snapshot,grid, alpha = alpha, h0 = h0)

    torque_grav = angular_momentum_gravity(snapshot,grid)
    
    return grid.R.mean(axis=0),mdot, torque_adv,torque_visc,torque_grav


def mass_advection(snapshot,grid):
    '''
    Compute the mass flux due to advection
    '''

    X0, Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    gridX, gridY = grid.X + X0, grid.Y + Y0


    # Compute the cell-centered quantities
    mdot_per_cell = -snapshot.gas.RHO * ((snapshot.gas.POS[:,0] - X0) * snapshot.gas.VELX + \
                                         (snapshot.gas.POS[:,1] - Y0) * snapshot.gas.VELY) / snapshot.gas.R
      
    snapshot.add_data(mdot_per_cell,'MASSFLUX')
    # interpolate onto the grid
    mdot_interp = disk_interpolate_primitive_quantities(snapshot,[gridX,gridY],\
                                                        quantities=['MASSFLUX'],method = 'nearest')[0]

    # Take the azimuthal integral to get the profile

    mdot = (mdot_interp * grid.R).mean(axis=0) * 2 * np.pi
    
    return mdot

def angular_momentum_advection(snapshot,grid):
    '''
    Compute the angular momentum flux due to advection
    '''

    X0, Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    gridX, gridY = grid.X + X0, grid.Y + Y0


    # Compute the cell-centered quantities
    jdot_per_cell = -snapshot.gas.RHO * ((snapshot.gas.POS[:,0] - X0) * (snapshot.gas.POS[:,1] - Y0) * \
                                         (snapshot.gas.VELY**2 - snapshot.gas.VELX**2) + \
                                         snapshot.gas.VELX * snapshot.gas.VELY * \
                                         ((snapshot.gas.POS[:,0] - X0)**2 - (snapshot.gas.POS[:,1] - Y0)**2)) / \
                                         snap.gas.R
    snapshot.add_data(jdot_per_cell,'TORQUEDENS')
    # interpolate onto the grid
    jdot_interp = disk_interpolate_primitive_quantities(snapshot,[gridX,gridY],\
                                                        quantities=['TORQUEDENS'],method = 'nearest')[0]

    # Take the azimuthal integral to get the profile

    jdot = (jdot_interp * grid.R).mean(axis=0) * 2 * np.pi
    
    return jdot


def angular_momentum_viscosity(snapshot,grid, alpha=0.1, h0=0.1):
    '''
    Compute the angular momentum flux due to viscous diffusion of momentum
    '''

    X0, Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    gridX, gridY = grid.X + X0, grid.Y + Y0


    
    # First, we need cell-centered velocity gradients
    gradientvx = compute_snapshot_gradient(snapshot,'VELX')
    gradientvy = compute_snapshot_gradient(snapshot,'VELY')
    # Add the gradients
    snapshot.add_data(gradientvx,'GRVX')
    snapshot.add_data(gradientvy,'GRVY')
    
    
    # Compute the cell-centered quantities
    GM = 1.0
    def nu(R): return alpha * h0**2 * np.sqrt(GM) * R**(0.5)
    nu_cell = nu(snapshot.gas.R)
    jdot_per_cell = -snapshot.gas.RHO * (2 * (snapshot.gas.POS[:,0] - X0) * (snapshot.gas.POS[:,1] - Y0) * \
                                         (snapshot.gas.GRVY[:,1] - snapshot.gas.GRVX[:,0]) + \
                                         ((snapshot.gas.POS[:,0] - X0)**2 - (snapshot.gas.POS[:,1] - Y0)**2) * \
                                         (snapshot.gas.GRVX[:,1] + snapshot.gas.GRVY[:,0])) * nu_cell / snapshot.gas.R 
    
    snapshot.add_data(jdot_per_cell,'TORQUEDENS')
    # interpolate onto the grid
    jdot_interp = disk_interpolate_primitive_quantities(snapshot,[gridX,gridY],\
                                                        quantities=['TORQUEDENS'],method = 'nearest')[0]


    # Take the azimuthal integral to get the profile

    jdot = (jdot_interp * grid.R).mean(axis=0) * 2 * np.pi
    
    return jdot


    
def angular_momentum_gravity(snapshot,grid):
    '''
    Compute the angular momentum flux due to an external (non-axisymmetric) gravitational source

    '''
    
    X0, Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    gridX, gridY = grid.X + X0, grid.Y + Y0

    
    # Compute the cell-centered quantities
    jdot_per_cell = -snapshot.gas.RHO * ((snapshot.gas.POS[:,0] - X0) * snapshot.gas.ACCE[:,1] - \
                                         (snapshot.gas.POS[:,1] - Y0) * snapshot.gas.ACCE[:,0])
    snapshot.add_data(jdot_per_cell,'TORQUEDENS')
    # interpolate onto the grid
    jdot_interp = disk_interpolate_primitive_quantities(snapshot,[gridX,gridY],\
                                                        quantities=['TORQUEDENS'],method = 'nearest')[0]


    # Take the azimuthal integral to get the profile
    djdotdR = (jdot_interp * grid.R).mean(axis=0) * 2 * np.pi
    
    
    # In the case of gravity, we need to carry out an additional integration step
    gridR = grid.R.mean(axis=0)
    jdot = cumtrapz(djdotdR[::-1],x = -gridR[::-1],initial=0)[::-1]
    
    
    return jdot
