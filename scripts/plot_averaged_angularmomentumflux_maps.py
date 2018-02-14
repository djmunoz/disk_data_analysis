import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.plotting import plot_slice, ImageData
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import h5py

import simulation_settings as simset


plt.rcParams['mathtext.fontset'] = "stix"


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

    if (len(sys.argv) < 5): orbit_init = 0
    else:   orbit_init = int(sys.argv[4])
    if (len(sys.argv) < 6): orbit_final = 100
    else:  orbit_final = int(sys.argv[5])

    filename = 'averaged_angmomflux_norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f.hdf5' % (orbit_init,orbit_final,qb,eb,h0,alpha,eta)
    print "Reading HDF5 file:",filename

    with h5py.File(filename, 'r') as hf:
        fluxj_adv_average =  hf['FluxAdvection'][:]
        fluxj_visc_average = hf['FluxViscosity'][:]
        fluxj_grav_average = hf['FluxGravity'][:]

    # create a grid
    NR, Nphi = fluxj_adv_average.shape[0],fluxj_adv_average.shape[1]
    Rmin, Rmax = 0.75, 16.0
    grid = dda.grid_polar(NR = NR, Nphi = Nphi,Rmin= Rmin,Rmax= Rmax,scale='linear')

    # Plotting...
    fig = plt.figure(1,figsize=(7.0,12.0))
    
    # First map
    ax1 = fig.add_axes([0.10,0.77,0.77,0.22])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    
    fluximage = ImageData()
    fluximage.data = fluxj_adv_average.T
    fluximage.pixelx = grid.R.shape[0]
    fluximage.pixely = grid.phi.shape[0]
    
    minval, maxval = -7, 7
    extent = [Rmin,Rmax,0,2*np.pi]
    cmap = cm.get_cmap('jet')
    ax, im = plot_slice(ax1,fluximage,cmap = cmap, 
                        vmin=minval,vmax=maxval,extent=extent,scale='linear')
    cbar = plt.colorbar(im,cax=cax)

    #c = ax1.contour(grid.R[0,:],grid.phi[:,0],fluxj_adv_average,
    #                levels=np.arange(1,10),colors='k',
    #                inline=1,fmt='%1.1f')
    c = ax1.contour(fluxj_adv_average,np.logspace(-2,0.5,10),colors='k',
                    extent=extent)

    ax1.set_xticklabels([])
    ticks = [0,np.pi/2,np.pi,1.5 * np.pi, 2 * np.pi ]
    ax1.set_yticks(ticks)
    ticklabels = [r'$0$',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'$\pi$']
    ax1.set_yticklabels(ticklabels)
    [tick.label.set_fontsize(16) for tick in ax1.yaxis.get_major_ticks()]
    ax1.set_aspect('auto')

    ax1.text(0.45,0.88,"inward angular momentum flux due to",transform = ax1.transAxes)
    ax1.text(0.45,0.78,"                     advection",transform = ax1.transAxes,size=16)
    
    # Second map
    ax2 = fig.add_axes([0.10,0.53,0.77,0.22])
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    
    del fluximage
    fluximage = ImageData()
    fluximage.data = fluxj_visc_average.T
    fluximage.pixelx = grid.R.shape[0]
    fluximage.pixely = grid.phi.shape[0]
    
    minval, maxval = None, None
    extent = [Rmin,Rmax,0,2*np.pi]
    cmap = cm.get_cmap('jet')
    ax2, im = plot_slice(ax2,fluximage,normalize=1.0,cmap = cmap, 
                        vmin=minval,vmax=maxval,extent=extent,scale='linear')
    cbar = plt.colorbar(im,cax=cax)

    c = ax2.contour(fluxj_visc_average,np.linspace(0.02,0.05,10),colors='k',
                    extent=extent)

    ax2.set_xticklabels([])
    ticks = [0,np.pi/2,np.pi,1.5 * np.pi, 2 * np.pi ]
    ax2.set_yticks(ticks)
    ticklabels = [r'$0$',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'$\pi$']
    ax2.set_yticklabels(ticklabels)
    [tick.label.set_fontsize(16) for tick in ax2.yaxis.get_major_ticks()]
    ax2.set_aspect('auto')

    ax2.text(0.45,0.88,"outward angular momentum flux due to",transform = ax2.transAxes)
    ax2.text(0.45,0.78,"                      viscosity",transform = ax2.transAxes,size=16)
    
    #Third map
    ax3 = fig.add_axes([0.10,0.29,0.77,0.22])
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    
    del fluximage
    fluximage = ImageData()
    fluximage.data = fluxj_grav_average.T
    fluximage.pixelx = grid.R.shape[0]
    fluximage.pixely = grid.phi.shape[0]
    
    minval, maxval = None, None
    extent = [Rmin,Rmax,0,2*np.pi]
    cmap = cm.get_cmap('jet')
    ax3, im = plot_slice(ax3,fluximage,normalize=1.0,cmap = cmap, 
                        vmin=minval,vmax=maxval,extent=extent,scale='linear')
    cbar = plt.colorbar(im,cax=cax)

    c = ax3.contour(fluxj_grav_average,np.logspace(-6,-2,10),colors='k',
                    extent=extent)

    #ax3.set_xticklabels([])
    ticks = [0,np.pi/2,np.pi,1.5 * np.pi, 2 * np.pi ]
    ax3.set_yticks(ticks)
    ticklabels = [r'$0$',r'$\frac{1}{2}\pi$',r'$\pi$',r'$\frac{3}{2}\pi$',r'$\pi$']
    ax3.set_yticklabels(ticklabels)
    [tick.label.set_fontsize(16) for tick in ax3.yaxis.get_major_ticks()]
    ax3.set_aspect('auto')
    
    ax3.text(0.45,0.88,"outward angular momentum flux due to",transform = ax3.transAxes)
    ax3.text(0.45,0.78,"      gravitational torques",transform = ax3.transAxes,size=16)
    

    # Radial profiles panel
    ax = fig.add_axes([0.10,0.07,0.735,0.18])
    ax.plot(grid.R.mean(axis=0),(fluxj_adv_average * grid.R).mean(axis=0)*2*np.pi,
            color = 'cornflowerblue')
    ax.plot(grid.R.mean(axis=0),(fluxj_visc_average * grid.R).mean(axis=0)*2*np.pi,
            color = 'forestgreen')
    ax.plot(grid.R.mean(axis=0),(fluxj_grav_average * grid.R).mean(axis=0)*2*np.pi,
            color = 'orange')

    fluxj_net = (fluxj_adv_average - fluxj_visc_average - fluxj_grav_average)
    ax.plot(grid.R.mean(axis=0),(fluxj_net * grid.R).mean(axis=0)*2*np.pi,
            color = 'firebrick')
    jdot_mean = ((fluxj_net * grid.R).mean(axis=0)*2*np.pi)[grid.R.mean(axis=0) < 10].mean()
    jdot_mean_other = 0.54
    ax.plot([1,50],[jdot_mean,jdot_mean],ls=':',color = 'firebrick')
    ax.plot([1,50],[jdot_mean_other,jdot_mean_other],ls='--',color = 'k',lw=0.9)
  
    ax.plot([1,50],[0.0,0.0],ls=':',color = 'gray',lw=0.4,zorder=0)

    ax.text(6,3.5,r'$\dot{J}_{\rm adv}\equiv 2\pi R\langle F_{J,\rm{adv}}\rangle_\phi$',
            transform = ax.transData, color = 'cornflowerblue',size=16)
    ax.text(10.5,1.8,r'$\dot{J}_{\rm visc}\equiv 2\pi R\langle F_{J,\rm{visc}}\rangle_\phi$',
            transform = ax.transData, color = 'forestgreen',size=16)
    ax.text(8,0.16,r'$\dot{J}_{\rm grav}\equiv 2\pi R\langle F_{J,\rm{grav}}\rangle_\phi$',
            transform = ax.transData, color = 'orange',size=16)
    ax.text(6,1.0,r'$\dot{J}_{\rm adv}-\dot{J}_{\rm visc}-\dot{J}_{\rm grav}$',
            transform = ax.transData, color = 'firebrick',size=18)

    ax.text(extent[1],0.95*jdot_mean,'%.2f' % jdot_mean,color = 'firebrick',
            horizontalalignment='right',verticalalignment='top',
            transform=ax.transData)
    ax.text(extent[1],0.95*jdot_mean_other,'%.2f' % jdot_mean_other,
            color = 'k',
            horizontalalignment='right',verticalalignment='top',
            transform=ax.transData)

    #ax.set_xlim(extent[0],extent[1])
    ax.set_xlim(1.0,extent[1])
    #ax.set_ylim(-2,3.9)
    ax.set_xlabel(r'$R/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$\langle\dot{J}\rangle_T\,/(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$',size=16,labelpad=0)

    ##############################################
    # save figure
    figfilename = './angular_momentum_flux_averaged_maps_norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f.png' % (orbit_init,orbit_final,qb,eb,h0,alpha,eta)
    print "Saving figure",figfilename
    fig.savefig(figfilename)
    fig.clf()

    ##############################################
    # Plot just the angular momentum transfer profile
    # Radial profiles panel

    fig = plt.figure(figsize=(7.0,4.5))
    fig.subplots_adjust(bottom=0.15,left=0.1,right=0.97,top=0.98)
    ax = fig.add_subplot(111)
    
    ax.plot(grid.R.mean(axis=0),(fluxj_adv_average * grid.R).mean(axis=0)*2*np.pi,
            color = 'cornflowerblue')
    ax.plot(grid.R.mean(axis=0),(fluxj_visc_average * grid.R).mean(axis=0)*2*np.pi,
            color = 'forestgreen')
    ax.plot(grid.R.mean(axis=0),(fluxj_grav_average * grid.R).mean(axis=0)*2*np.pi,
            color = 'orange')

    fluxj_net = (fluxj_adv_average - fluxj_visc_average - fluxj_grav_average)
    ax.plot(grid.R.mean(axis=0),(fluxj_net * grid.R).mean(axis=0)*2*np.pi,
            color = 'firebrick')
    jdot_mean = ((fluxj_net * grid.R).mean(axis=0)*2*np.pi)[grid.R.mean(axis=0) < 10].mean()
    jdot_mean_other = 0.54
    jdot_mean, jdot_mean_other = 1.05, 0.75
    ax.plot([1,50],[jdot_mean,jdot_mean],ls=':',color = 'firebrick')
    ax.plot([1,50],[jdot_mean_other,jdot_mean_other],ls='--',color = 'k',lw=0.9)
  
    ax.plot([1,50],[0.0,0.0],ls=':',color = 'gray',lw=0.4,zorder=0)

    ax.text(6,3.5,r'$\dot{J}_{\rm adv}\equiv 2\pi R\langle F_{J,\rm{adv}}\rangle_\phi$',
            transform = ax.transData, color = 'cornflowerblue',size=16)
    ax.text(10.5,1.8,r'$\dot{J}_{\rm visc}\equiv 2\pi R\langle F_{J,\rm{visc}}\rangle_\phi$',
            transform = ax.transData, color = 'forestgreen',size=16)
    ax.text(8,0.16,r'$\dot{J}_{\rm grav}\equiv 2\pi R\langle F_{J,\rm{grav}}\rangle_\phi$',
            transform = ax.transData, color = 'orange',size=16)
    ax.text(6,jdot_mean,r'$\dot{J}_{\rm adv}-\dot{J}_{\rm visc}-\dot{J}_{\rm grav}$',
            transform = ax.transData, color = 'firebrick',
            size=18,verticalalignment='bottom')

    ax.text(extent[1],0.95*jdot_mean,'%.2f' % jdot_mean,color = 'firebrick',
            horizontalalignment='right',verticalalignment='top',
            transform=ax.transData)
    ax.text(extent[1],0.95*jdot_mean_other,'%.2f' % jdot_mean_other,
            color = 'k',
            horizontalalignment='right',verticalalignment='top',
            transform=ax.transData)

    #ax.set_xlim(extent[0],extent[1])
    ax.set_xlim(1.0,extent[1])
    ax.set_ylim(-0.5,4.2)
    ax.set_xlabel(r'$R/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$\langle\dot{J}\rangle_T\,/(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$',size=16,labelpad=0)

    # save figure
    figfilename = './angular_momentum_transfer_average_norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f.pdf' % (orbit_init,orbit_final,qb,eb,h0,alpha,eta)
    print "Saving figure",figfilename
    fig.savefig(figfilename)
    fig.clf()

    
