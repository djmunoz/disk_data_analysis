from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from ..disk_hdf5 import snapHDF5 as rs
from string import *
import os
from matplotlib import cm
from scipy.ndimage import geometric_transform
from scipy.interpolate import griddata
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid import make_axes_locatable
import matplotlib.axes as maxes
from scipy.ndimage.filters import gaussian_filter

# path to output files
directory_base = '/n/ghernquist/dmunoz/PLANET_WAKES_2D/SIMULATIONS/'

dir_array = [directory_base+'NEPTUNE_Nrings.128/',
             directory_base+'NEPTUNE_Nrings.256/',
             directory_base+'NEPTUNE_Nrings.512/',
             ]

dir_array = [directory_base+'NEPTUNE-VISCOUS_Nrings.128_NU.1E-5/',
             directory_base+'NEPTUNE-VISCOUS_Nrings.256_NU.1E-5/',
             directory_base+'NEPTUNE-VISCOUS_Nrings.512_NU.1E-5/'
             ]


dir_array = [directory_base+'JUPITER_Nrings.128/',
             directory_base+'JUPITER_Nrings.256/',
             directory_base+'JUPITER_Nrings.512/',
             ]

dir_array = [directory_base+'JUPITER-VISCOUS_Nrings.128_NU.1E-5/',
             directory_base+'JUPITER-VISCOUS_Nrings.256_NU.1E-5/',
             directory_base+'JUPITER-VISCOUS_Nrings.512_NU.1E-5/'
             ]

current_dir = split(os.getcwd(), "/")[len(split(os.getcwd(), "/"))-1]
# directory='./'

base = "output/"

rangeX = [1.5, 6.5]
rangeY = [1.5, 6.5]
range_R = [0.4, 2.5]
range_theta = [0, 2*np.pi]

Lx = rangeX[1] - rangeX[0]
Ly = rangeY[1] - rangeY[0]


snap_list = [0, 500, 1000]
snap_list = [0, 100, 200]
snap_list = [0, 100]

inner_radius = 0.25
outer_radius = 2.5


skip_cartesian = 1

frame_labels = True


def cartesian2polar(coords, inputshape, origin):
    """Coordinate transform for converting a polar array to Cartesian coordinates. 
    inputshape is a tuple containing the shape of the polar array. origin is a
    tuple containing the x and y indices of where the origin should be in the
    output array."""

    r_index, theta_index = coords

    r = r_index * (rangeX[1] - rangeX[0])/2.0/inputshape[0]
    theta = theta_index * 2.0*np.pi/inputshape[1] + np.pi

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    i = np.round(x/Lx*inputshape[0]) + origin[0]
    j = np.round(y/Ly*inputshape[0]) + origin[1]

    return (i, j)


for num in snap_list:
    # prepare figures first
    if not (skip_cartesian):
        fig = plt.figure(1, figsize=(7.5*len(dir_array)+1.0, 6.5))
        fig.subplots_adjust(top=0.96, bottom=0.08,
                            right=0.90, left=0.05, hspace=0.0)
    if (frame_labels):
        fig_polar = plt.figure(1, figsize=(7.5*len(dir_array)+2.0, 6.0))
        fig_polar.subplots_adjust(
            top=0.99, bottom=0.03, right=0.98, left=0.02, hspace=0.0)
    else:
        fig_polar = plt.figure(1, figsize=(7.0*len(dir_array), 4.5))
        fig_polar.subplots_adjust(
            top=0.99, bottom=0.03, right=0.98, left=0.02, hspace=0.0)

    count_dir = 0
    for directory in dir_array:
        count_dir += 1

        if (('NEPTUNE' in directory) | ('NEPTUNE' in current_dir)):
            min_scale = -0.1
            max_scale = 1.0
            mass_ratio = 0.0001
            tag = 'nep'
        elif (('JUPITER' in directory) | ('JUPITER' in current_dir)):
            min_scale = -0.5
            max_scale = 2.0
            mass_ratio = 0.001
            tag = 'jup'

        if (('128' in directory) | ('128' in current_dir)):
            radial_zones = 128
        elif (('256' in directory) | ('256' in current_dir)):
            radial_zones = 256
        elif (('512' in directory) | ('512' in current_dir)):
            radial_zones = 512
        elif (('1024' in directory) | ('1024' in current_dir)):
            radial_zones = 1024

        print("SNAPSHOT #", num)
        # open the snapshot header
        filename = directory+base+"snap_"+str(num).zfill(3)
        header = rs.snapshot_header(filename)
        time = header.time/2.0/math.pi
        BoxX, BoxY = header.boxsize, header.boxsize
        pos = rs.read_block(filename, "POS ", parttype=0)
        dens = rs.read_block(filename, "RHO ", parttype=0)
        vel = rs.read_block(filename, "VEL ", parttype=0)

        x, y = pos[:, 0], pos[:, 1]
        vx, vy = vel[:, 0], vel[:, 1]

        # define regular grid spatially covering input data
        n = 1024
        xg = np.linspace(rangeX[0], rangeX[1], n)
        yg = np.linspace(rangeY[0], rangeY[1], n)
        delta_x = np.diff(xg).mean()
        delta_y = np.diff(yg).mean()
        X, Y = np.meshgrid(xg, yg)

        # data to interpolate
        z = vx

        # interpolate Z values on defined grid
        Z = griddata(np.vstack((x.flatten(), y.flatten())).T,
                     np.vstack(z.flatten()), (X, Y), method='linear').reshape(X.shape)
        # mask nan values, so they will not appear on plot
        Zm = np.ma.masked_where(np.isnan(Z), Z)

        #dZdY = np.gradient(Zm.flatten('F'),2*delta_y).reshape(Y.shape).T
        dZdY = np.gradient(Zm, delta_y)[0].reshape(Zm.shape)

        # data to interpolate
        z = vy

        # interpolate Z values on defined grid
        Z = griddata(np.vstack((x.flatten(), y.flatten())).T,
                     np.vstack(z.flatten()), (X, Y), method='linear').reshape(X.shape)
        # mask nan values, so they will not appear on plot
        Zm = np.ma.masked_where(np.isnan(Z), Z)

        #dZdX = np.gradient(Zm.flatten('C'),2*delta_x).reshape(X.shape)
        dZdX = np.gradient(Zm, delta_x)[1].reshape(Zm.shape)

        # data to interpolate
        z = dens

        # interpolate Z values on defined grid
        Z = griddata(np.vstack((x.flatten(), y.flatten())).T,
                     np.vstack(z.flatten()), (X, Y), method='linear').reshape(X.shape)
        # mask nan values, so they will not appear on plot
        Zm = np.ma.masked_where(np.isnan(Z), Z)

        Omega = 1.0
        Vortensity = (dZdX - dZdY + 2*Omega)/Zm

        smooth_scale = (outer_radius - inner_radius) / \
            radial_zones / delta_x / 1.5
        Vortensity = gaussian_filter(Vortensity, smooth_scale, mode='nearest')

        if (num == 0):
            Vortensity0 = Vortensity

        print("Doing geometric transformation to R-theta plane...")
        # create polar-coordinate array
        Vortensity_polar = geometric_transform(Vortensity.T, cartesian2polar, output_shape=(Vortensity.T.shape[0], Vortensity.T.shape[0]),
                                               extra_keywords={'inputshape': Vortensity.T.shape, 'origin': (Vortensity.T.shape[0]/2, Vortensity.T.shape[1]/2)})
        if (num == 0):
            Vortensity0_polar = Vortensity_polar

        ####################################################################
        # plot
        if not (skip_cartesian):

            ax = fig.add_subplot(1, len(dir_array), count_dir)
            divider = make_axes_locatable(ax)

            cmap = cm.get_cmap('jet')
            im = ax.imshow(np.log10(Vortensity/Vortensity0), origin='lower',
                           vmin=min_scale, vmax=max_scale,
                           extent=[rangeX[0], rangeX[1], rangeY[0], rangeY[1]], cmap=cmap)

            xlabel("$x$", fontsize=16)
            ylabel("$y$", fontsize=16)
            ax.set_xlim(rangeX[0], rangeX[1])
            ax.set_ylim(rangeY[0], rangeY[1])

            xticks, yticks = ax.xaxis.get_majorticklocs(), ax.yaxis.get_majorticklocs()
            ax.xaxis.set_ticklabels(['%d' % (xticks[n] - 0.5*BoxX)
                                     for n in range(len(xticks))])
            ax.yaxis.set_ticklabels(['%d' % (yticks[n] - 0.5 * BoxY)
                                     for n in range(len(yticks))])

            R = inner_radius
            phi = arange(0, 2*pi, 0.01)
            ax.fill(0.5*(rangeX[0] + rangeX[1]) + cos(phi)*R,
                    0.5*(rangeY[0] + rangeY[1]) + sin(phi)*R, fc="white", fill=True, ec='none')
            R = outer_radius
            phi = arange(0, pi, 0.01)
            ax.fill(np.append(0.5*(rangeX[0] + rangeX[1]) + cos(phi)*R, [rangeX[0], rangeX[0], rangeX[1], rangeX[1]]),
                    np.append(0.5*(rangeY[0] + rangeY[1]) + sin(phi)*R, [0.5*(rangeY[0]+rangeY[1]), rangeY[1], rangeY[1],
                                                                         0.5*(rangeY[0]+rangeY[1])]), fc="w", fill=True, ec='none')
            phi = arange(pi, 2*pi, 0.01)
            ax.fill(np.append(0.5*(rangeX[0] + rangeX[1]) + cos(phi)*R, [rangeX[1], rangeX[1], rangeX[0], rangeX[0]]),
                    np.append(0.5*(rangeY[0] + rangeY[1]) + sin(phi)*R, [0.5*(rangeY[0]+rangeY[1]), rangeY[0], rangeY[0],
                                                                         0.5*(rangeY[0]+rangeY[1])]), fc="w", fill=True, ec='none')

            cax = divider.new_horizontal("5%", pad=0.05, axes_class=maxes.Axes)
            fig.add_axes(cax)
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label(r"log$_{10}$($\omega/\omega_0$)", fontsize=16)

            ax.text(0.5, 0.9, "Time= %4.1f orbits " %
                    time, fontsize=20, transform=ax.transAxes)
            ########################
            # write image (change output format by changing filename extension)
            fig.savefig(directory+"vortensity-interpolated_" +
                        str(num).zfill(3)+".png")
            fig.clf()
     ####################################################################
        # plot polar image
        print("Creating image in polar coordinates...")

        if (frame_labels):
            margin_x = 0.04
            axes = [margin_x+(count_dir-1)*(1.0 - margin_x) / len(dir_array),
                    0.135, (1.0 - margin_x)/len(dir_array), 0.78]
            ax = fig_polar.add_axes(axes)

        else:
            ax = fig_polar.add_subplot(1, len(dir_array), count_dir)

        # plot field (change here what you want to plot)
        cmap = cm.get_cmap('jet')

        im = ax.imshow(np.log10(Vortensity_polar.T/Vortensity0_polar.T),
                       vmin=min_scale, vmax=max_scale, origin='lower', interpolation='nearest', cmap=cmap,
                       extent=[0, range_R[1], range_theta[0], range_theta[1]])

        ###############
        ax.xaxis.tick_top()
        ax.set_xlim(range_R[0], range_R[1])
        ax.set_ylim(range_theta[0], range_theta[1])

        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0)/1.6)

        if ((frame_labels) & (count_dir == 1)):
            xlabel("$R$", fontsize=20)
            ylabel(r"$\phi$", fontsize=22)
            ax.tick_params(axis='x', labelsize=14)
            ax.tick_params(axis='y', labelsize=16)

            ax.xaxis.set_label_position('top')
            ax.yaxis.set_ticks([0, 0.5*np.pi, np.pi, 1.5*np.pi, 2*np.pi])
            ax.yaxis.set_ticklabels(
                [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])

            cax = fig_polar.add_axes([0.3, 0.08, 0.95, 0.03])
            cax.hold(True)
            cbar = plt.colorbar(im, cax=cax, orientation='horizontal')
            tick_interval = 0.2
            cbar.set_ticks(np.arange(
                tick_interval*int(min_scale/tick_interval), max_scale*1.1, tick_interval))
            cbar.set_label(r"$\log_{10}(\zeta'/\zeta'_0)$", fontsize=26)
            cbar.ax.tick_params(labelsize=12)
        else:
            ax.axes.get_xaxis().set_visible(False)
            ax.axes.get_yaxis().set_visible(False)

        ax.text(0.55, 0.9, "Time= %4.1f orbits " %
                time, fontsize=20, transform=ax.transAxes, color='w')

    if (frame_labels):
        append = ''
    else:
        append = 'nolabels_'
    ###############
    # write image (change output format by changing filename extension)
    fig_polar.savefig("vortensity-array_polar_"+tag +
                      "_"+append+str(num).zfill(3)+".png")
    fig_polar.clf()
