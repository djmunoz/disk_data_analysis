import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from disk_data_analysis.plotting import plot_slice, ImageData
import matplotlib.cm as cm
import sys

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"



# Diego J. Munoz
# 2017

'''
Script to generate 2-D polar maps of the angular momentum flux
distribution in hydro simulations of disk

'''

# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out


NPIXELS=1024

ROTATING_FRAME = False

if __name__ == "__main__":

    
    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h0=simset.h
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


    
    run_name= simset.run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h0,alpha,eta))
    print "Reading from simulation run:", run_name
    directory = simset.run_path+run_name
    snap_path = directory + simset.base + simset.snap_base
    
    mu = qb / (1.0 + qb)


    # create a grid
    NR, Nphi = NPIXELS, NPIXELS
    Rmin, Rmax = 0.75, 16.0
    grid = dda.grid_polar(NR = NR, Nphi = Nphi,Rmin= Rmin,Rmax= Rmax,scale='linear')

    for snapnum in snap_list:
        print "SNAPSHOT #",snapnum
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS','VELX','VELY','RHO','ACCE','R'],parttype=0)


        # compute angular momentum flux maps
        fluxj_grav = dda.compute_angular_momentum_flux_gravity(snap,grid) / Mdot0
        fluxj_visc = dda.compute_angular_momentum_flux_viscosity(snap,grid) / Mdot0
        fluxj_adv = dda.compute_angular_momentum_flux_advection(snap,grid) / Mdot0


        #######################################################################
        # Plotting...
        fig = plt.figure(1,figsize=(7.0,12.0))
        
        # First map
        ax1 = fig.add_axes([0.10,0.77,0.77,0.22])
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        
        fluximage = ImageData()
        fluximage.data = fluxj_adv.T
        fluximage.pixelx = grid.X.shape[0]
        fluximage.pixely = grid.Y.shape[0]
        
        minval, maxval = None, None
        extent = [Rmin,Rmax,0,2*np.pi]
        cmap = cm.get_cmap('jet')
        ax, im = plot_slice(ax1,fluximage,normalize=6.0e-4,cmap = cmap, 
                            vmin=minval,vmax=maxval,extent=extent,scale='linear')
        cbar = plt.colorbar(im,cax=cax)
        
        ax1.set_xticklabels([])
        ticks = [0,np.pi/2,np.pi,1.5 * np.pi, 2 * np.pi ]
        ax1.set_yticks(ticks)
        ticklabels = [r'$0$',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'$\pi$']
        ax1.set_yticklabels(ticklabels)
        [tick.label.set_fontsize(16) for tick in ax1.yaxis.get_major_ticks()]
        ax1.set_aspect('auto')
        
        # Second map
        ax2 = fig.add_axes([0.10,0.53,0.77,0.22])
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        
        del fluximage
        fluximage = ImageData()
        fluximage.data = fluxj_visc.T
        fluximage.pixelx = grid.X.shape[0]
        fluximage.pixely = grid.Y.shape[0]
        
        minval, maxval = None, None
        extent = [Rmin,Rmax,0,2*np.pi]
        cmap = cm.get_cmap('jet')
        ax, im = plot_slice(ax2,fluximage,normalize=6.0e-4,cmap = cmap, 
                            vmin=minval,vmax=maxval,extent=extent,scale='linear')
        cbar = plt.colorbar(im,cax=cax)
        
        ax2.set_xticklabels([])
        ticks = [0,np.pi/2,np.pi,1.5 * np.pi, 2 * np.pi ]
        ax2.set_yticks(ticks)
        ticklabels = [r'$0$',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'$\pi$']
        ax2.set_yticklabels(ticklabels)
        [tick.label.set_fontsize(16) for tick in ax2.yaxis.get_major_ticks()]
        ax2.set_aspect('auto')
        
        #Third map
        ax3 = fig.add_axes([0.10,0.29,0.77,0.22])
        divider = make_axes_locatable(ax3)
        cax = divider.append_axes("right", size="3%", pad=0.08)
        
        del fluximage
        fluximage = ImageData()
        fluximage.data = fluxj_grav.T
        fluximage.pixelx = grid.X.shape[0]
        fluximage.pixely = grid.Y.shape[0]
        
        minval, maxval = None, None
        extent = [Rmin,Rmax,0,2*np.pi]
        cmap = cm.get_cmap('jet')
        ax, im = plot_slice(ax3,fluximage,normalize=6.0e-4,cmap = cmap, 
                            vmin=minval,vmax=maxval,extent=extent,scale='linear')
        cbar = plt.colorbar(im,cax=cax)
        
        #ax3.set_xticklabels([])
        ticks = [0,np.pi/2,np.pi,1.5 * np.pi, 2 * np.pi ]
        ax3.set_yticks(ticks)
        ticklabels = [r'$0$',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'$\pi$']
        ax3.set_yticklabels(ticklabels)
        [tick.label.set_fontsize(16) for tick in ax3.yaxis.get_major_ticks()]
        ax3.set_aspect('auto')

        # Radial profiles panel
        ax = fig.add_axes([0.10,0.07,0.735,0.18])
        
        ax.plot(grid.R.mean(axis=0),(fluxj_adv * grid.R).mean(axis=0)*2*np.pi)
        ax.plot(grid.R.mean(axis=0),(fluxj_visc * grid.R).mean(axis=0)*2*np.pi)
        ax.plot(grid.R.mean(axis=0),(fluxj_grav * grid.R).mean(axis=0)*2*np.pi)
        
        #ax.plot(radii,torque_grav)
        
        ax.set_xlim(extent[0],extent[1])
        #ax.set_ylim(-2,3.9)
        ax.set_xlabel(r'$R/a_{\rm b}$',size=22)
        ax.set_ylabel(r'$\mu_{\rm b}(dt_{\rm grav}/dA) /\,(\dot{M}_0\Omega_{\rm b})}$',size=16,labelpad=0)

    

    
        ##############################################
        # save figure
        figfilename = './angular_momentum_flux_maps_%i.png' % snapnum
        print "Saving figure",figfilename
        fig.savefig(figfilename)
        fig.clf()

        print "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    
