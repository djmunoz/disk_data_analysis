# Script to plot a collection of different 2-D slices (or projections)
# as different panels into the same figure file.

# Diego J. Munoz
#
##############################
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid import make_axes_locatable
import  matplotlib.axes as maxes 
import matplotlib.cm as cm
import my_colortable as m_ct
import numpy as np
from pylab import * 
import readsnapHDF5 as rs
from mpl_toolkits.axes_grid import make_axes_locatable 
import  matplotlib.axes as maxes 
from matplotlib.colors import LogNorm
import os
from string import split
import sys
import glob
###################################################################


min_scale= -5.6
max_scale= 0.8
sigma0=6.0e-4
racc=0.02


if __name__=="__main__":
	#first, identify simulation and directory with the data
        #binary and disk properties
	snap_init=int(sys.argv[1])
	snap_end=int(sys.argv[2])
        if (len(sys.argv) > 3):
            Lx = int(sys.argv[3])
            Ly = int(sys.argv[4])
        else:
	    Lx, Ly = None, None 
        # paths to files
	directory = "./"
	base = "output/"
	

	#prepare figure
	fig_sizex= 5.0
	fig_sizey= 4.0
	fig = plt.figure(1, figsize=(fig_sizex,fig_sizey))
	fig.subplots_adjust(top=0.91,left=0.17,right=0.99,bottom=0.06,hspace=0.02,wspace=-0.003)
	ax = fig.add_subplot(1,1,1)
	cax=fig.add_axes([0.175,0.92,0.81,0.02])

	print "Reading from output directory:",base
	snap_base="snap_"
	
	print snap_init,snap_end
	snap_list = range(snap_init,snap_end)
	time_list = []
        for num in snap_list:
                filename=directory+base+snap_base+str(num).zfill(3)
                #open the snapshot header
                header = rs.snapshot_header(filename)
                time_list.append(header.time)
                if (Lx is None) & (Ly is None):
                        Lx, Ly = header.boxsize, header.boxsize
                X0, X1 = 0.5 * header.boxsize - 0.5 * Lx,0.5 * header.boxsize + 0.5 * Lx
                Y0, Y1 = 0.5 * header.boxsize - 0.5 * Ly,0.5 * header.boxsize + 0.5 * Ly
                
	        #look for image file
	        projfilename=directory+base+"density_field_"+str(num).zfill(3)
                #check if the file exists
                if not os.path.isfile(projfilename):
		        print "\nERROR: projection not found"
                        continue

		f = open(projfilename,"rb")
		pixelx = np.fromfile(f, dtype=np.int32, count=1)[0]
		pixely = np.fromfile(f, dtype=np.int32, count=1)[0]
		dens = np.fromfile(f, dtype=np.float32, count=pixelx*pixely)
		f.close()
		dens=dens.reshape(pixelx,pixely)
		dens[(dens < 0.0)] = 1.0e-20
		dens[(dens < 0.0)] = 1.0e-20
		s, gamma, hue, r  = 1.7, 0.90, 1.1, 1.6 #old
		#s, gamma, hue, r  = 0.5, 1.0, 1.0, -1.5 
		s, gamma, hue, r  = 0.4, 0.7, 1.0, -1.6 
		cmap,nlow,nhigh = m_ct.my_colortable(start=s,rots=r,hue=hue,gamma=gamma)
		im=ax.imshow(np.log10(dens/sigma0).T,vmin=min_scale, vmax=max_scale,
			     origin='lower',
			     interpolation='nearest',
		     extent=[X0,X1,Y0,Y1], cmap=cmap)
		cmap.set_under(color='k')

		cbar=plt.colorbar(im, cax=cax,orientation='vertical',ticklocation='right')
		cbar.set_clim(-5,1)
		cbar.set_ticks([-5,-4,-3,-2,-1,0,1])
		cbar.set_label(r"log$_{10}$($\Sigma/\Sigma_0$)",fontsize=16,labelpad=10)
		for label in cbar.ax.get_xticklabels(): label.set_fontsize(11)
					
			
		ax.set_xlim(X0,X1)
		ax.set_ylim(Y0,Y1)
		ax.locator_params(nbins = 6)
		ax.tick_params('both',length=8,width=1.2,which='major')
		xticks=ax.get_xticks()
		yticks=ax.get_yticks()
		
		ax.set_xticks(xticks[1:-1])
		ax.set_yticks(yticks[1:-1])

		figname="./figure-disk_slice_%i" % num
		print "Saving figure",figname+".pdf"
		#plt.savefig(figname+".pdf")
		plt.savefig(figname+".png",dpi=200)
		plt.clf()

	exit()
	
