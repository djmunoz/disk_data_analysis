import numpy as np
import matplotlib.pyplot as plt
import sys
import disk_data_analysis.circumbinary as dda

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"

from disk_data_analysis.plotting import plot_slice,ImageData

if __name__ == '__main__':

    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h0=simset.h
    alpha=simset.alpha
    eta = float(sys.argv[3])
    snap_num = int(sys.argv[4])
    
    snap_path = simset.find_simulation(qb,eb,h0,alpha,eta)

    fig = plt.figure(2,figsize = (12.0,5.5))
    fig.subplots_adjust(top=0.98,bottom=0.12,left=0.07,right=0.92)
    
    image = ImageData()
    image.read(snap_path+'density_field_%03i' % snap_num)
    image.interpolateneg()
    image.replaceneg()

    ax = fig.add_subplot(111)
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

    cbar = plt.colorbar(im,cax=cax)
    cbar.set_label(r'$\Sigma/\Sigma_{\rm out}$',size=20)

    
    fig.savefig('./density_field.png')
