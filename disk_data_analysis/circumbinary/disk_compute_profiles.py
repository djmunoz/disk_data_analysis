from __future__ import print_function
import numpy as np
from scipy.integrate import cumtrapz
from .disk_orbital_elements import compute_disk_eccentriciy_vector,\
    compute_disk_angular_momentum_vector, compute_disk_semimajor
from .disk_simulation_data import *

def compute_cell_size_profile(radii,cell_vol,reference_radius):
# In unstructured meshes, we need to estimate a representative cell size
# as a function of radius

    if (np.array(reference_radius).shape[0] == 1):
        cell_ind = (radii < 1.25*reference_radius[0]) & (radii > 0.85*reference_radius[0])
        #for 2-D simulations
        cell_size = np.sqrt(cell_vol[cell_ind]/np.pi).mean()
    else:
        cell_size = []
        for ref_rad in  reference_radius:
            cell_ind = (radii < 1.25*ref_rad) & (radii > 0.85*ref_rad)
            #for 2-D simulations
            cell_size.append(np.sqrt(cell_vol[cell_ind]/np.pi).mean())
            
    return np.array(cell_size)


def compute_smoothed_radial_placements(radii,volumes,Rmin,Rmax,NRmin=128):

    '''
    Function to compute a collection of effective radial zones from unstructured-mesh data.

    '''

    # First attempt: a logarithmically-space radial grid
    nzones = (np.log10(Rmax)-np.log10(Rmin))/np.sqrt(volumes.min()/np.pi) * Rmin/2
    radius_arr = np.logspace(np.log10(Rmin),np.log10(Rmax),nzones)

    while True:
        # Bin the original radial placements
        
        digitized = np.digitize(radii, radius_arr)
        bin_pop = np.array([radii[digitized == i].shape[0] for i in range(1, len(radius_arr))])

        # Remove bins that underpopulated, but not more than 5 at a time
        ind = (bin_pop <= 0.03 * bin_pop.max()) & (bin_pop <= sorted(bin_pop)[4])
        if (bin_pop[ind].shape[0] == 0): break
        radius_arr = np.append(np.append(Rmin,radius_arr[1:][np.invert(ind)]),Rmax)
        # Smooth the sharp transitions
        w = np.hanning(5)
        radius_arr = np.append(Rmin,np.convolve(w/w.sum(),radius_arr[1:],mode='valid'))

        if (radius_arr.shape[0] <= NRmin): break
    
    return radius_arr
        

def compute_weighted_azimuthal_average(quantity,radii,cell_vol,reference_radius,bin_width):
    # The azimuthal average is weighted by cells volume

    cell_ind = (radii < (reference_radius+3.0*bin_width)) & (radii > (reference_radius-3.0*bin_width)) 
    #gaussian kernel function strongly peaked at r=reference_radius
    kernel=np.exp(-(radii[cell_ind]-(reference_radius))**2/2/(0.5*bin_width)**2)
    weight = kernel * cell_vol[cell_ind] / radii[cell_ind]
    norm = weight.sum()

    #compute averaged quantity
    mu = (quantity[cell_ind]*weight).sum()/norm
    #s = np.sqrt((rho[cell_ind]**2*weight).sum()/norm - mu**2)
    
    return mu


def compute_profiles(quantity,radii,volumes,rad_list,semimajor=False):
    #Loop over values of radius in rad_list
    
    def average_quantity(r):
                #now, bin cells using mean cell size as a guess for bin width
                bin_width = 5* compute_cell_size_profile(radii,volumes,[r])
                quantity_av=compute_weighted_azimuthal_average(quantity,radii,volumes,r,bin_width)
                return quantity_av
    
    
    quantity_av = [average_quantity(rad) for rad in rad_list]
    #quantity_av = average_quantity(rad_list[2])

    return quantity_av
    


def compute_mdot_profile(snapshot,rad_list,code="AREPO"):


    vol = snapshot.gas.MASS/snapshot.gas.RHO
    mdot = -snapshot.gas.RHO * snapshot.gas.VELR 
    ind = snapshot.gas.ID > -2

    mdot_av=compute_profiles(mdot[ind],snapshot.gas.R[ind],vol[ind],rad_list)

    mdot_av[:] *= 2 * np.pi *  rad_list[:]
    return np.array(mdot_av)


def compute_jdot_adv_profile(snapshot,rad_list,code="AREPO"):


    vol = snapshot.gas.MASS/snapshot.gas.RHO
    jdot = -snapshot.gas.RHO * snapshot.gas.VELR * snapshot.gas.VELPHI * snapshot.gas.R
    ind = snapshot.gas.ID > -2

    jdot_av=compute_profiles(jdot[ind],snapshot.gas.R[ind],vol[ind],rad_list)

    jdot_av[:] *= 2 * np.pi *  rad_list[:]
    return np.array(jdot_av)


def compute_jdot_visc_profile(snapshot,rad_list,code="AREPO", alpha = 0.1, h0 = 0.1):

    X0, Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    vol = snapshot.gas.MASS/snapshot.gas.RHO
    # add gradients
    gradientvx = compute_snapshot_gradient(snapshot,'VELX')
    gradientvy = compute_snapshot_gradient(snapshot,'VELY')
    snapshot.add_data(gradientvx,'GRVX')
    snapshot.add_data(gradientvy,'GRVY')
    GM = 1.0
    def nu(R): return alpha * h0**2 * np.sqrt(GM) * R**(0.5)
    nu_cell = nu(snapshot.gas.R)
    jdot = -snapshot.gas.RHO * (2 * (snapshot.gas.POS[:,0] - X0) * (snapshot.gas.POS[:,1] - Y0) * \
                                (snapshot.gas.GRVY[:,1] - snapshot.gas.GRVX[:,0]) + \
                                ((snapshot.gas.POS[:,0] - X0)**2 - (snapshot.gas.POS[:,1] - Y0)**2) * \
                                (snapshot.gas.GRVX[:,1] + snapshot.gas.GRVY[:,0])) * nu_cell / snapshot.gas.R

    ind = snapshot.gas.ID > -2
    jdot_av=compute_profiles(jdot[ind],snapshot.gas.R[ind],vol[ind],rad_list)

    jdot_av[:] *= 2 * np.pi *  rad_list[:]
    return np.array(jdot_av)


def compute_jdot_grav_profile(snapshot,rad_list,code="AREPO", alpha = 0.1, h0 = 0.1):
    
    X0, Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    vol = snapshot.gas.MASS/snapshot.gas.RHO
    # add gradients
    jdotdens = -snapshot.gas.RHO * ((snapshot.gas.POS[:,0] - X0) * (-snapshot.gas.ACCE[:,1]) - \
                                    (snapshot.gas.POS[:,1] - Y0) * (-snapshot.gas.ACCE[:,0])) 
    ind = snapshot.gas.ID > -2
    jdotdens_av=compute_profiles(jdotdens[ind],snapshot.gas.R[ind],vol[ind],rad_list)
    jdot_av = 2 * np.pi *  cumtrapz((jdotdens_av[:] * rad_list[:])[::-1],x = -rad_list[::-1],initial=0)[::-1]
    #jdot_av = 2 * np.pi *  cumtrapz((jdotdens_av[:] * rad_list[:]),x = rad_list,initial=0)
    
    return np.array(jdot_av)



    
    # Compute the cell-centered quantities
    jdotdens_per_cell = -snapshot.gas.RHO * ((snapshot.gas.POS[:,0] - X0) * (-snapshot.gas.ACCE[:,1]) - \
                                             (snapshot.gas.POS[:,1] - Y0) * (-snapshot.gas.ACCE[:,0])) 
    snapshot.add_data(jdotdens_per_cell,'TORQUEDENS')
    # interpolate onto the grid
    jdotdens_interp = disk_interpolate_primitive_quantities(snapshot,[gridX,gridY],\
                                                            quantities=['TORQUEDENS'],method = 'nearest')[0]

    del snapshot.gas.TORQUEDENS
    
    # In the case of gravity, we need to carry out an additional integration step
    gridR = grid.R.mean(axis=0)
    jdot_interp = cumtrapz((jdotdens_interp * grid.R)[:,::-1],x = -gridR[::-1],initial=0,axis=1)[:,::-1] / grid.R

    
    return jdot_interp


def compute_eccentricity_profile(snapshot,rad_list,semimajor=False,code="AREPO"):

    ind = snapshot.gas.ID > -2
    ecc = compute_disk_eccentriciy_vector(snapshot.gas.POS,snapshot.gas.VEL)
    if (semimajor):
        radii = compute_disk_semimajor(snapshot.gas.POS-0.5 *snapshot.header.boxsize,snapshot.gas.VEL)[ind]
    else:
        radii = snapshot.gas.R[ind]

    ex_av = np.array(compute_profiles(ecc[ind,0],radii,snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],
                                      rad_list,semimajor=semimajor))
    ey_av = np.array(compute_profiles(ecc[ind,1],radii,snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],
                                      rad_list,semimajor=semimajor))

    
    return np.sqrt(ex_av**2 + ey_av**2)

def compute_density_profile(snapshot,rad_list,semimajor=False,code="AREPO"):

    ind = snapshot.gas.ID > -2
    if (semimajor):
        radii = compute_disk_semimajor(snapshot.gas.POS,snapshot.gas.VEL)[ind]
    else:
        radii = snapshot.gas.R[ind]

    rho = np.array(compute_profiles(snapshot.gas.RHO[ind],radii,snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],
                                    rad_list,semimajor=semimajor))
    
    return rho

def compute_orbital_elements_profile(snapshot,rad_list,code="AREPO"):

    X0,Y0 = 0.5 * snapshot.header.boxsize, 0.5 * snapshot.header.boxsize
    snapshot.gas.POS[:,0]-=X0
    snapshot.gas.POS[:,1]-=Y0
    ecc = compute_disk_eccentriciy_vector(snapshot.gas.POS,snapshot.gas.VEL)
    j = compute_disk_angular_momentum_vector(snapshot.gas.POS,snapshot.gas.VEL)
    
    ind = snapshot.gas.ID > -2
    ex_av = np.array(compute_profiles(ecc[ind,0]*snapshot.gas.RHO[ind],
                                      snapshot.gas.R[ind],snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],rad_list))
    ey_av = np.array(compute_profiles(ecc[ind,1]*snapshot.gas.RHO[ind],
                                      snapshot.gas.R[ind],snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],rad_list))
    j_av = np.array(compute_profiles(j[ind,2]*snapshot.gas.RHO[ind],
                                     snapshot.gas.R[ind],snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],rad_list))
    rho_av = np.array(compute_profiles(snapshot.gas.RHO[ind],
                                       snapshot.gas.R[ind],snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],rad_list))
    
    return ex_av/rho_av,ey_av/rho_av,j_av/rho_av
