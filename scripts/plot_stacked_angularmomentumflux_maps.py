import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.plotting import plot_slice, ImageData
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import h5py

plt.rcParams['mathtext.fontset'] = "stix"


filename_list = ['averaged_angmomflux_norbits3500-3550.hdf5',
                 'averaged_angmomflux_norbits3550-3600.hdf5',
                 'averaged_angmomflux_norbits3600-3650.hdf5']
                 
NPIXELS=1024


if __name__ == "__main__":

    

    # create a grid
    NR, Nphi = NPIXELS, NPIXELS
    Rmin, Rmax = 0.75, 16.0
    grid = dda.grid_polar(NR = NR, Nphi = Nphi,Rmin= Rmin,Rmax= Rmax,scale='linear')

    nmin,nmax = 1e30,0

    for kk,filename in enumerate(filename_list):

        if ('norbits' in filename):
            n0 = int(filename.split('.hdf5')[0].split('norbits')[1].split('-')[0])
            n1 = int(filename.split('.hdf5')[0].split('norbits')[1].split('-')[1])
            if (n0 < nmin): nmin = n0
            if (n1 > nmax): nmax = n1


        with h5py.File(filename, 'r') as hf:
            fluxj_adv_average =  hf['FluxAdvection'][:]
            fluxj_visc_average = hf['FluxViscosity'][:]
            fluxj_grav_average = hf['FluxGravity'][:]
        
        if (kk == 0):
            fluxj_adv_stacked =  fluxj_adv_average
            fluxj_visc_stacked = fluxj_visc_average
            fluxj_grav_stacked = fluxj_grav_average
        else:
            fluxj_adv_stacked +=  fluxj_adv_average
            fluxj_visc_stacked += fluxj_visc_average
            fluxj_grav_stacked += fluxj_grav_average


    fluxj_adv_stacked /= len(filename_list)
    fluxj_visc_stacked /= len(filename_list)
    fluxj_grav_stacked /= len(filename_list)

    # Plotting...
    fig = plt.figure(1,figsize=(7.0,12.0))
    
    # First map
    ax1 = fig.add_axes([0.10,0.77,0.77,0.22])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.08)
    
    fluximage = ImageData()
    fluximage.data = fluxj_adv_stacked.T
    fluximage.pixelx = grid.R.shape[0]
    fluximage.pixely = grid.phi.shape[0]
    
    minval, maxval = -4, 4
    extent = [Rmin,Rmax,0,2*np.pi]
    cmap = cm.get_cmap('jet')
    ax, im = plot_slice(ax1,fluximage,cmap = cmap, 
                        vmin=minval,vmax=maxval,extent=extent,scale='linear')
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label(r'$F_{J,{\rm adv}}\,/(\dot{M}_0a_{\rm b}\Omega_{\rm b})$',size=20)

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
    fluximage.data = fluxj_visc_stacked.T
    fluximage.pixelx = grid.R.shape[0]
    fluximage.pixely = grid.phi.shape[0]
    
    minval, maxval = 0.0, 0.06
    extent = [Rmin,Rmax,0,2*np.pi]
    cmap = cm.get_cmap('jet')
    ax2, im = plot_slice(ax2,fluximage,normalize=1.0,cmap = cmap, 
                        vmin=minval,vmax=maxval,extent=extent,scale='linear')
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label(r'$F_{J,{\rm visc}}\,/(\dot{M}_0a_{\rm b}\Omega_{\rm b})$',size=20)

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
    fluximage.data = fluxj_grav_stacked.T
    fluximage.pixelx = grid.R.shape[0]
    fluximage.pixely = grid.phi.shape[0]
    
    minval, maxval = None, None
    extent = [Rmin,Rmax,0,2*np.pi]
    cmap = cm.get_cmap('jet')
    ax3, im = plot_slice(ax3,fluximage,normalize=1.0,cmap = cmap, 
                        vmin=minval,vmax=maxval,extent=extent,scale='linear')
    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label(r'$F_{J,{\rm grav}}\,/(\dot{M}_0a_{\rm b}\Omega_{\rm b})$',size=20)

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
    
    ax.plot(grid.R.mean(axis=0),(fluxj_adv_stacked * grid.R).mean(axis=0)*2*np.pi,
            color = 'cornflowerblue')
    ax.plot(grid.R.mean(axis=0),(fluxj_visc_stacked * grid.R).mean(axis=0)*2*np.pi,
            color = 'forestgreen')
    ax.plot(grid.R.mean(axis=0),(fluxj_grav_stacked * grid.R).mean(axis=0)*2*np.pi,
            color = 'orange')
    
    ax.text(6,3.5,r'$\dot{J}_{\rm adv}\equiv 2\pi R\langle F_{J,\rm{adv}}\rangle_\phi$',
            transform = ax.transData, color = 'cornflowerblue',size=16)
    ax.text(10.5,1.8,r'$\dot{J}_{\rm visc}\equiv 2\pi R\langle F_{J,\rm{visc}}\rangle_\phi$',
            transform = ax.transData, color = 'forestgreen',size=16)
    ax.text(8,0.2,r'$\dot{J}_{\rm grav}\equiv 2\pi R\langle F_{J,\rm{grav}}\rangle_\phi$',
            transform = ax.transData, color = 'orange',size=16)

    
    ax.set_xlim(extent[0],extent[1])
    #ax.set_ylim(-2,3.9)
    ax.set_xlabel(r'$R/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$\langle\dot{J}\rangle_T\,/(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$',size=16,labelpad=0)
    
    

    
    ##############################################
    tag = '_norbits%i-%i' % (nmin,nmax)
    # save figure
    figfilename = './angular_momentum_flux_averaged_maps'+tag+'.png'
    print "Saving figure",figfilename
    fig.savefig(figfilename)
    fig.clf()

    
    ##############################################
    # Plot just the angular momentum transfer profile
    # Radial profiles panel

    fig = plt.figure(figsize=(7.0,4.5))
    fig.subplots_adjust(bottom=0.15,left=0.1,right=0.97,top=0.98)
    ax = fig.add_subplot(111)
    
    ax.plot(grid.R.mean(axis=0),(fluxj_adv_stacked * grid.R).mean(axis=0)*2*np.pi,
            color = 'cornflowerblue')
    ax.plot(grid.R.mean(axis=0),(fluxj_visc_stacked * grid.R).mean(axis=0)*2*np.pi,
            color = 'forestgreen')
    ax.plot(grid.R.mean(axis=0),(fluxj_grav_stacked * grid.R).mean(axis=0)*2*np.pi,
            color = 'orange')

    fluxj_net = (fluxj_adv_stacked - fluxj_visc_stacked - fluxj_grav_stacked)
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

    # save figure
    figfilename = './angular_momentum_transfer_average'+tag+'.pdf'
    print "Saving figure",figfilename
    fig.savefig(figfilename)
    fig.clf()
