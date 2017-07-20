import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable

if __name__ == "__main__":

    snap = dda.get_snapshot_data('../data/snap_',0,['POS','VELX','VELY','RHO','ACCE','R','ID','MASS'])

    # Binary properties
    eb, qb = 0.0 , 1.0 # eccentricity and mass ratio
    mu = qb / (1.0 + qb)

    # Normalization of the accretion rate
    Mdot0 = 8.63857932377e-06
    # Normalization of the outer density
    Sigma_out = 1.0238548710e-4
    
    # Compute forces
    time = snap.header.time
    X0, Y0 = 0.5 * snap.header.boxsize, 0.5 * snap.header.boxsize
    xb, yb, _, _ = orbit_in_time(time + np.pi, eb)
    pos1 = [mu * xb + X0, mu * yb + Y0, 0] # primary
    accel1 = (1 - mu) * dda.compute_external_gravforce(snap,XYZ = pos1,softening=0.026 * 2.8 * (1- mu))
                                           
    pos2 = [-(1 - mu) * xb + X0, -(1 - mu) * yb + Y0, 0] # secondary
    accel2 = mu * dda.compute_external_gravforce(snap,XYZ=pos2,softening=0.026 * 2.8 * mu)

    ind = snap.gas.ID > 0

    print (np.abs((accel1 + accel2)[ind,:2] - snap.gas.ACCE[ind,:2])/snap.gas.ACCE[ind,:2]).mean(axis=0)
    
    # Separate regions
    dx, dy = snap.gas.POS[:,0] - X0, snap.gas.POS[:,1] - Y0
    dx1, dy1 = snap.gas.POS[:,0] - pos1[0], snap.gas.POS[:,1] - pos1[1]
    dx2, dy2 = snap.gas.POS[:,0] - pos2[0], snap.gas.POS[:,1] - pos2[1] 
    dr = np.sqrt(dx * dx  + dy * dy)
    dr1 = np.sqrt(dx1 * dx1  + dy1 * dy1)
    dr2 = np.sqrt(dx2 * dx2  + dy2 * dy2)
    # circum-single disks
    csd_trunc = 0.23
    cbd_trunc = 2.35
    csd_region = (dr1 < csd_trunc) | (dr2 < csd_trunc)
    cbd_region = (dr1 > cbd_trunc) & (dr2 > cbd_trunc)
    str_region = ((dr1 <= cbd_trunc) | (dr2 <= cbd_trunc)) & ((dr1 >= csd_trunc) & (dr2 >= csd_trunc))

    fig = plt.figure(1,figsize=(6.0,6.0))
    plt.plot(snap.gas.POS[csd_region,0],snap.gas.POS[csd_region,1],'b.',ms=1.6)
    plt.plot(snap.gas.POS[cbd_region,0],snap.gas.POS[cbd_region,1],'r.',ms=1.6)
    plt.plot(snap.gas.POS[str_region,0],snap.gas.POS[str_region,1],'.',color='orange',ms=1.6)
    plt.xlim(0.5 * snap.header.boxsize - 5, 0.5 * snap.header.boxsize + 5)
    plt.ylim(0.5 * snap.header.boxsize - 5, 0.5 * snap.header.boxsize + 5)
    plt.xlabel(r'$x$',size=18)
    plt.ylabel(r'$y$',size=18)
    plt.axes().set_aspect('equal')
    plt.savefig('./accretion_regions.png')
    fig.clf()

    from disk_data_analysis.plotting import plot_slice, ImageData
    import matplotlib.cm as cm
    
    fig = plt.figure(2,figsize = (12.0,5.5))
    fig.subplots_adjust(top=0.98,bottom=0.12,left=0.07,right=0.92)
    
    image = ImageData()
    image.read('../data/density_field_close_000')
    image.interpolateneg()
    image.replaceneg()
    
    #ax = fig.add_subplot(121)
    ax = fig.add_axes([-0.17,0.13,0.85,0.85])
    minval,maxval = 10**(-4.12),10**(1.03)
    extent = [77.5,82.5,77.5,82.5]
    ax,im = plot_slice(ax,image,normalize=Sigma_out,
                       vmin=minval,vmax=maxval,extent=extent,scale='log')
    ax.set_xlabel(r'$x/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$y/a_{\rm b}$',size=22)
    ax.set_aspect('equal')


    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_xticks()]
    ax.set_xticklabels(ticklabels)
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_yticks()]
    ax.set_yticklabels(ticklabels)
    

    image = ImageData()
    image.read('../data/density_field_000')
    image.interpolateneg()
    image.replaceneg()

    
    #ax = fig.add_subplot(122)
    ax = fig.add_axes([0.3,0.13,0.85,0.85])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.08)

    
    minval,maxval = 10**(-4.9),10**(0.25)
    extent = [74.0,86.0,74.0,86.0]
    ax, im = plot_slice(ax,image,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,scale='log')
    ax.set_xlabel(r'$x/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$y/a_{\rm b}$',size=22)
    ax.set_aspect('equal')

    
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_xticks()]
    ax.set_xticklabels(ticklabels)
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_yticks()]
    ax.set_yticklabels(ticklabels)

    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label(r'$\Sigma/\Sigma_{\rm out}$',size=20)

    
    fig.savefig('./density_field.png')
    
    #plt.show()


    pixelxcoord, pixelycoord = np.meshgrid(np.arange(image.pixelx)*1.0/image.pixelx,\
                                       np.arange(image.pixely)*1.0/image.pixely)

    Lx, Ly = extent[1] - extent[0], extent[3] - extent[2]
    pixelxcoord, pixelycoord = (pixelxcoord.T - 0.5) * Lx + 0.5 * snap.header.boxsize,\
                               (pixelycoord.T - 0.5) * Ly + 0.5 * snap.header.boxsize
    
    
    dxpixel1, dypixel1 = pixelxcoord - pos1[0], pixelycoord - pos1[1]
    dxpixel2, dypixel2 = pixelxcoord - pos2[0], pixelycoord - pos2[1]
    dxpixel, dypixel = pixelxcoord - 0.5 * snap.header.boxsize, pixelycoord - 0.5 * snap.header.boxsize
    
    drpixel1 = np.sqrt(dxpixel1**2 + dypixel1**2)
    drpixel2 = np.sqrt(dxpixel2**2 + dypixel2**2)
    drpixel = np.sqrt(dxpixel**2 + dypixel**2)
    
    #we do the same region separation as before
    csd_region = (drpixel1 < csd_trunc) | (drpixel2 < csd_trunc)
    cbd_region = drpixel > cbd_trunc
    str_region = ((drpixel <= cbd_trunc) ) &\
                 ((drpixel1 >= csd_trunc) & (drpixel2 >= csd_trunc))
    
    image_csd = ImageData()
    image_cbd = ImageData()
    image_str = ImageData()

    image_csd.read('../data/density_field_000')
    image_cbd.read('../data/density_field_000')
    image_str.read('../data/density_field_000')

    image_cbd.data[np.invert(cbd_region)] = np.nan    
    image_csd.data[np.invert(csd_region)] = np.nan
    image_str.data[np.invert(str_region)] = np.nan

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax,im = plot_slice(ax,image_cbd,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,
                       cmap=cm.get_cmap('Blues'),scale='log')
    ax,im = plot_slice(ax,image_csd,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,
                       cmap=cm.get_cmap('Greens'),scale='log')
    ax,im = plot_slice(ax,image_str,normalize=6.0e-4,vmin=minval,vmax=maxval,extent=extent,
                       cmap=cm.get_cmap('Reds'),scale='log')
    ax.set_xlabel(r'$x$',size=18)
    ax.set_ylabel(r'$y$',size=18)
    ax.set_aspect('equal')
    
    fig.savefig('./density_field_colormaps.png')
    #plt.show()
    fig.clf()



    ############
    # 3. Compute spatial distribution of torque

    # a. First, create a torque density map

    # We assume we already have the mass density map, so we only need the acceleration map

    # treat the cell-centered acceleration as a reconstructible quantity
    print "hehe"
    dens_data = image.data
    snap.add_data(accel1,'GRPHI')

    # and interpolate using the 'nearest' method onto the same pixel grid
    gradphi_interp1 = dda.disk_interpolate_gradient_quantities(snap,[pixelxcoord,pixelycoord],\
                                                               quantities=['GRPHI'],method = 'nearest')[0]

    snap.add_data(accel2,'GRPHI')
    gradphi_interp2 = dda.disk_interpolate_gradient_quantities(snap,[pixelxcoord,pixelycoord],\
                                                               quantities=['GRPHI'],method = 'nearest')[0]
    
    # we need to take the cross product between the binary's (relative) position vector and
    # the acceleration vector for each pixel.
    torque_dens_data = np.zeros([image.pixelx,image.pixely])
    torque_dens_data = dens_data[:,:] * (-(xb * gradphi_interp1[1][:,:] - yb * gradphi_interp1[0][:,:])/(1 - mu)
                                         +(xb * gradphi_interp2[1][:,:] - yb * gradphi_interp2[0][:,:])/mu)
    

    
    print torque_dens_data.max(),torque_dens_data.min()
    torqueimage = ImageData()
    torqueimage.data = torque_dens_data/Mdot0
    torqueimage.pixelx = torque_dens_data.shape[0]
    torqueimage.pixely = torque_dens_data.shape[1]
    
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
            print  np.log10(np.abs(tick)),np.sign(tick)
            ctick_labels += [r"$%i^{%i}$" % (int(np.sign(tick))*10,int(np.log10(np.abs(tick))))]
        else:
             ctick_labels += [r"$%i$" % tick]
        print tick,ctick_labels[-1]

    cbar.set_ticklabels(ctick_labels[1:-1])
    cbar.set_label(r'$(dt_{\rm grav}/dA) /\,[(\dot{M}_0/M_{\rm b})\Omega_{\rm b}]}$',size=20)
    ax.set_xlabel(r'$x/a_{\rm b}$',size=22)
    ax.set_ylabel(r'$y/a_{\rm b}$',size=22)

    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_xticks()]
    ax.set_xticklabels(ticklabels)
    ticklabels = ["%i" % (t - 0.5 * snap.header.boxsize)  for t in ax.get_yticks()]
    ax.set_yticklabels(ticklabels)

    fig.savefig('./torque_density_field.png')
    fig.clf()

    # pixel area
    pixelarea = (Lx/torqueimage.pixelx) * (Ly/torqueimage.pixely)
    print torqueimage.pixelx,torqueimage.pixely
    
    # total mass
    print image.data.sum() * pixelarea
    print image.data[cbd_region].sum() * pixelarea
    print image.data[csd_region].sum() * pixelarea
    print image.data[str_region].sum() * pixelarea

    # total torque
    print "\n", torqueimage.data[csd_region].sum() * pixelarea
    print torqueimage.data[cbd_region].sum() * pixelarea
    print torqueimage.data[str_region].sum() * pixelarea 
    print torqueimage.data.sum() * pixelarea   
    
    torque_per_cell = - snap.gas.MASS[:]/(1- mu) * (xb * accel1[:,1] - yb * accel1[:,0])\
                      + snap.gas.MASS[:]/mu * (xb * accel2[:,1] - yb * accel2[:,0])
    
    torque_per_cell /= Mdot0
    
    csd_region = (dr1 < csd_trunc) | (dr2 < csd_trunc)
    cbd_region = (dr > cbd_trunc) 
    str_region = ((dr <= cbd_trunc)) & ((dr1 >= csd_trunc) & (dr2 >= csd_trunc))

    print "\n", torque_per_cell[csd_region].sum()
    print torque_per_cell[cbd_region].sum()
    print torque_per_cell[str_region].sum()
    print torque_per_cell.sum()
