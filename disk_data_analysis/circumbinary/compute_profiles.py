import numpy as np
import matplotlib.pyplot as plt
#from joblib import Parallel, delayed
#import multiprocessing
from string import split
#import disk_simulation_data 
from disk_hdf5 import readsnapHDF5 as rs

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
    nzones = (np.log10(Rmax)-np.log10(Rmin))/np.sqrt(volumes.min()) * Rmin
    radius_arr = np.logspace(np.log10(Rmin),np.log10(Rmax),nzones)

    while True:
        # Bin the original radial placements
        
        digitized = np.digitize(radii, radius_arr)
        bin_pop = np.array([radii[digitized == i].shape[0] for i in range(1, len(radius_arr))])

        # Remove bins that underpopulated, but not more than 5 at a time
        ind = (bin_pop <= 0.1 * bin_pop.max()) & (bin_pop < sorted(bin_pop)[5])
        if (bin_pop[ind].shape[0] == 0): break
        radius_arr = np.append(np.append(radius_arr[0],radius_arr[1:][np.invert(ind)]),radius_arr[-1])
        # Smooth the sharp transitions
        w = np.hanning(7)
        radius_arr = np.convolve(w/w.sum(),radius_arr,mode='valid')

        if (radius_arr.shape[0] <= NRmin): break
        
    return radius_arr
        

def compute_weighted_azimuthal_average(quantity,radii,cell_vol,reference_radius,bin_width):
    # The azimuthal average is weighted by cells volume

    cell_ind = (radii < (reference_radius+1.5*bin_width)) & (radii > (reference_radius-1.5*bin_width)) 
    #gaussian kernel function strongly peaked at r=reference_radius
    kernel=np.exp(-(radii[cell_ind]-(reference_radius))**2/2/(0.4*bin_width)**2)
    weight = kernel * cell_vol[cell_ind] / radii[cell_ind]
    norm = weight.sum()

    #compute averaged quantity
    mu = (quantity[cell_ind]*weight).sum()/norm
    #s = np.sqrt((rho[cell_ind]**2*weight).sum()/norm - mu**2)
    
    return mu


def compute_profiles(quantity,radii,volumes,rad_list,num_cores=1):
    #Loop over values of radius in rad_list
    
    def average_quantity(r):
                #now, bin cells using mean cell size as a guess for bin width
                bin_width = 5* compute_cell_size_profile(radii,volumes,[r])
                quantity_av=compute_weighted_azimuthal_average(quantity,radii,volumes,r,bin_width)
                return quantity_av
    
    
    quantity_av = [average_quantity(rad) for rad in rad_list]
    #quantity_av = average_quantity(rad_list[2])

    return quantity_av
    

def compute_mdot_profile(snapshot,rad_list,code="AREPO",delta_t=0):


    time = snapshot.header.time + delta_t

    mdot = 2 * np.pi * snapshot.gas.RHO * snapshot.gas.VELR * snapshot.gas.R
    ind = snapshot.gas.ID > -2
        
    mdot_av=compute_profiles(mdot[ind],snapshot.gas.R[ind],snapshot.gas.MASS[ind]/snapshot.gas.RHO[ind],rad_list)
    
    return np.array(mdot_av)

