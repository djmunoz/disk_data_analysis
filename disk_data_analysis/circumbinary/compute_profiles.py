import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
from string import split
import disk_simulation_data 
import readsnapHDF5 as rs

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
    print len(quantity_av)
    #quantity_av = average_quantity(rad_list[2])

    return quantity_av
    

def compute_mdot_profile(rad_list,num,path_to_files,snapshot_base,code="AREPO",delta_t=0):


    filename_prefix=path_to_files+snapshot_base
    time = disk_simulation_data.get_snapshot_time(filename_prefix,num,code=code,snapinterval=delta_t)
    print "SNAPSHOT #",num," TIME=",time
    dens,radius,vr,vol,ids = disk_simulation_data.get_snapshot_data(filename_prefix,num,["RHO","R","VELR","VOL","ID"],code=code)
    mdot = 2*np.pi * dens * vr * radius
    ind = ids > -2
        
    mdot_av=compute_profiles(mdot[ind],radius[ind],vol[ind],rad_list)
    
    print len(mdot_av),len([time]+ mdot_av)
    return np.array([time]+mdot_av)

