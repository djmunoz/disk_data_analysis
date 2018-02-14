import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import h5py

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"

'''
Script to compute the torque density map from a snapshot when
a reconstructed image for the density field is not available

'''



NPIXELS  = 1024

ROTATING_FRAME = True

# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out



if __name__ == "__main__":

 #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h=simset.h
    alpha=simset.alpha
    if (len(sys.argv) > 3):
        eta = float(sys.argv[3])
    else:
        eta = 1.0

    

    if (len(sys.argv) < 5): init_snap = 0
    else:   init_snap = int(sys.argv[4])
    if (len(sys.argv) < 6): final_snap = 1001
    else:  final_snap = int(sys.argv[5])
    if (len(sys.argv) < 7): snap_step = 1
    else:  snap_step = int(sys.argv[6])


    snap_list = np.arange(init_snap,final_snap+1,snap_step)    

    
    run_name= simset.run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h,alpha,eta))
    
    print "Reading from simulation run:", run_name
    directory = simset.run_path+run_name
    snap_path = directory + simset.base + simset.snap_base
    
    mu = qb / (1.0 + qb)


    # Create empty images
    density_average = np.zeros([NPIXELS,NPIXELS])
    torquedensity_average = np.zeros([NPIXELS,NPIXELS])

    # For a cartesian grid...
    NX, NY = NPIXELS, NPIXELS
    Xmin, Xmax = -6.0, 6.0
    Ymin, Ymax = -6.0, 6.0
    
    
    for snapnum in snap_list:

        # Read snapshot information
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS','VELX','VELY','RHO','ACCE','R','ID','MASS'],parttype=0)

        # Create grid
        grid = dda.grid_cartesian(NX = NX, NY = NY,Xmin = Xmin, Xmax = Xmax, Ymin = Ymin, Ymax = Ymax)
        grid.X, grid.Y = grid.X + snap.header.boxsize * 0.5, grid.Y  +  snap.header.boxsize * 0.5

        
        # Compute location of binary components
        time = snap.header.time
        print "SNAPSHOT #%i; time=%.6f " % (snapnum,time)
        
        X0, Y0 = 0.5 * snap.header.boxsize, 0.5 * snap.header.boxsize
        xb, yb, _, _ = orbit_in_time(time + np.pi, eb)
        if (ROTATING_FRAME):
            angle = np.arctan2(yb,xb)
            xb, yb = 1.0, 0.0
            X, Y = snap.gas.POS[:,0] - X0, snap.gas.POS[:,1] - Y0
            snap.gas.POS[:,0] = X * np.cos(angle) + Y * np.sin(angle) + X0
            snap.gas.POS[:,1] = -X * np.sin(angle) + Y * np.cos(angle) + Y0
            del X, Y
            
        pos1 = [mu * xb + X0, mu * yb + Y0, 0] # primary
        pos2 = [-(1 - mu) * xb + X0, -(1 - mu) * yb + Y0, 0] # secondary

        
        # Compute forces
        accel1 = (1 - mu) * dda.compute_external_gravforce_from_snapshot(snap,XYZ = pos1,softening=0.026 * 2.8 * (1- mu))
        accel2 = mu * dda.compute_external_gravforce_from_snapshot(snap,XYZ=pos2,softening=0.026 * 2.8 * mu)

        # 1. Reconstruct density by interpolation...
        
        rho_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],quantities=['RHO'])[0]

        # 2. Compute the torque associated to each cell
        torquedens_per_cell = - snap.gas.RHO[:]/(1- mu) * (xb * accel1[:,1] - yb * accel1[:,0])\
                              + snap.gas.RHO[:]/mu * (xb * accel2[:,1] - yb * accel2[:,0])
        torquedens_per_cell /= Mdot0
        snap.add_data(torquedens_per_cell,'TORQUEDENS')
        # and interpolate using the 'nearest' method onto the same pixel grid
        torque_interp = dda.disk_interpolate_primitive_quantities(snap,[grid.X,grid.Y],\
                                                                  quantities=['TORQUEDENS'],method = 'nearest')[0]
    
        density_average += rho_interp
        torquedensity_average += torque_interp


    # Divide by the total number of snapshots added
    density_average /= len(snap_list)
    torquedensity_average /= len(snap_list)
        
    from disk_data_analysis.plotting import plot_slice, ImageData
    import matplotlib.cm as cm

    # Save image data into HDF5 file
    with h5py.File('averaged_density.hdf5', 'w') as hf:
        hf.create_dataset("Density",  data = density_average)
    
    with h5py.File('averaged_torquedensity.hdf5', 'w') as hf:
        hf.create_dataset("TorqueDensity",  data = torquedensity_average)
    
    
    # Create image data structure for the density field...
    densityimage = ImageData()
    densityimage.data = density_average.T
    densityimage.pixelx = grid.X.shape[0]
    densityimage.pixely = grid.Y.shape[0]

    # ... and plot it   
    fig = plt.figure()
    fig.subplots_adjust(top=0.97,right=0.88,bottom=0.12,left=0.10)
    ax = fig.add_subplot(111)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.08)

    
    minval,maxval = 10**(-4.9),10**(0.25)
    extent = [74.0,86.0,74.0,86.0]
    ax, im = plot_slice(ax,densityimage,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,scale='log')
    ax.set_xlabel(r'$x/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$y/a_{\rm b}$',size=22)
    ax.set_aspect('equal')

    
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_xticks()]
    ax.set_xticklabels(ticklabels)
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_yticks()]
    ax.set_yticklabels(ticklabels)

    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label(r'$\Sigma/\Sigma_{\rm out}$',size=20)

    figfilename = './density_field_averaged.png'
    print "Saving figure",figfilename
    fig.savefig(figfilename)
    fig.clf()
    
    # Create image data structure for the torque density field...
    torqueimage = ImageData()
    torqueimage.data = torquedensity_average.T
    torqueimage.pixelx = grid.X.shape[0]
    torqueimage.pixely = grid.Y.shape[0]
    

    # ... and plot it
    fig = plt.figure()
    fig.subplots_adjust(top=0.97,right=0.88,bottom=0.12,left=0.10)
    ax = fig.add_subplot(111)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    vmax,vmin,absvmin = 5e3,-5e3,0.01
    ax,im = plot_slice(ax,torqueimage,extent=extent,scale='semilog',cmap=cm.get_cmap('bwr'),vmax=vmax,vmin=vmin,absvmin=absvmin)
    cbar = plt.colorbar(im,cax=cax)
    maxlog=int(np.ceil(np.log10(vmax)))
    minlog=int(np.ceil(np.log10(-vmin)))
    logthresh = int(np.round(np.log10(absvmin)))
    
    #generate logarithmic ticks 
    ctick_locs=([-(10**x) for x in xrange(minlog,logthresh-1,-2)]+[0.0] +[(10**x) for x in xrange(logthresh,maxlog+1,2)] )
    cbar.set_ticks(ctick_locs[1:-1])
    ctick_labels=[]
    for tick in ctick_locs:
        if (tick != 0):
            ctick_labels += [r"$%i^{%i}$" % (int(np.sign(tick))*10,int(np.log10(np.abs(tick))))]
        else:
             ctick_labels += [r"$%i$" % tick]

    cbar.set_ticklabels(ctick_labels[1:-1])
    cbar.set_label(r'$(dt_{\rm grav}/dA) /\,[(\dot{M}_0/M_{\rm b})\Omega_{\rm b}]}$',size=20)
    ax.set_xlabel(r'$x/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$y/a_{\rm b}$',size=22)

    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_xticks()]
    ax.set_xticklabels(ticklabels)
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_yticks()]
    ax.set_yticklabels(ticklabels)
    
    
    figfilename = './torque_density_field_averaged.png'
    print "Saving figure",figfilename
    fig.savefig(figfilename)
    fig.clf()
    
    print ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    

