# Script to plot a collection of different 2-D slices (or projections)
# as different panels into the same figure file.

# Diego J. Munoz
#
##############################
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
import orbital
###################################################################


min_scale= -5.6
max_scale= 0.8
sigma0=6.0e-4
racc=0.02


if __name__=="__main__":
	#first, identify simulation and directory with the data
        #binary and disk properties
	qb=float(sys.argv[1])
	eb=float(sys.argv[2])
	h=0.1
	alpha=0.1
	
	snap = int(sys.argv[3]) #initial snapshot number for which we want to make density plots

	if (len(sys.argv) > 4):
		HORIZONTAL=int(sys.argv[4])
	else:
		HORIZONTAL = 0

        # paths to files
	run_path = "/data2/djmunoz/CIRCUMBINARY_DISKS_2D/PRODUCTION_RUNS/"
	
	#for comparison, we need two different runs of similar configuration
	run_base ="restart-woboundary-standardres-binary_norbits1950"
	run_name= run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f/" % (qb,eb,h,alpha))
	print "Reading from simulation run:", run_name
	directory = run_path+run_name

	output_base_label1 = "sink_short_2200"
	output_base_label2 = "sink_short_vhires_2200"
	output_base_label3 = "sink_short_vvhires_2200"
	output_base_label2 = "sink-s0.025_short_vhires_2200"
	output_base_label3 = "sink-s0.025_short_vvhires_2200"
	
	output_base_label_list = [output_base_label1,output_base_label2,output_base_label3]
	resolution_factor = [1,4,16]
	accretion_radius = [0.02,0.01,0.005]
	accretion_radius = [0.02,0.02,0.02]

	#prepare figure
	if not HORIZONTAL:
		fig_sizex= 4.0
		fig_sizey= 11.5
		fig = plt.figure(1, figsize=(fig_sizex,fig_sizey))
		fig.subplots_adjust(top=0.91,left=0.17,right=0.99,bottom=0.06,hspace=0.02,wspace=-0.003)
		ax1 = fig.add_subplot(3,1,1)
		ax2 = fig.add_subplot(3,1,2)
		ax3 = fig.add_subplot(3,1,3)
		cax=fig.add_axes([0.175,0.92,0.81,0.02])
	else:
		fig_sizex= 11.5
		fig_sizey= 4.0
		fig = plt.figure(1, figsize=(fig_sizex,fig_sizey))
		fig.subplots_adjust(top=0.86,left=0.08,right=0.91,bottom=0.0,hspace=0.0,wspace=-0.003)
		ax1 = fig.add_subplot(1,3,1)
		ax2 = fig.add_subplot(1,3,2)
		ax3 = fig.add_subplot(1,3,3)
		cax=fig.add_axes([0.92,0.03,0.02,0.80])


	ax_list = [ax1,ax2,ax3]


	for jj,base_label in enumerate(output_base_label_list):
		base = "output_restart_"+base_label+"/"
		print "Reading from output directory:",base
		snap_base="snap_"

                #check the number of orbits at the zeroth snapshot
		if (jj == 0): #use first output directory as reference
			filename=directory+base+snap_base+str(snap).zfill(3)
			header = rs.snapshot_header(filename)
			BoxX,BoxY = header.boxsize,header.boxsize
			time0 = header.time
			norbits=int(time0/(2*np.pi))
			ecc_anom = orbital.KeplerEquation(time0 + np.pi,eb)
			posx1, posy1 = -qb/(1+qb)*(np.cos(ecc_anom) - eb),-qb/(1+qb)*(np.sqrt(1 - eb * eb) * np.sin(ecc_anom))
			lims = -0.7+posx1,0.7+posx1,-0.7+posy1,0.7+posy1
				

		time = 1.0e30
		for snapfile in sorted(glob.glob(directory+base+snap_base+"*")):
			if (np.abs(time0-rs.snapshot_header(snapfile).time) < np.abs(time0-time)):
				time = rs.snapshot_header(snapfile).time
				num = int(split(split(snapfile,snap_base)[1],".hdf5")[0])
		print num
		print time
	
		filename=directory+base+snap_base+str(num).zfill(3)
		header = rs.snapshot_header(filename)
		time = header.time
		
		ax = ax_list[jj]
		ax.set_xlim(lims[0],lims[1])
		ax.set_ylim(lims[2],lims[3])
		ax.set_aspect(1)
		X0,X1 = ax.get_xlim()
		Y0,Y1 = ax.get_ylim()

		#look for image file
		projfilename=directory+base+"density_field_"+str(num).zfill(3)
                #check if the file exists
		if not os.path.isfile(projfilename):
			print "\nERROR: projection not found"
			print "execute: "
			print "./Arepo param.txt 4 %i 1024 1024 0 1 2  %.2f  %.2f  %.2f %.2f 0.0 > OUTPUT_IM.%i" \
			    % (num,(X0+0.5*BoxX),(X1+0.5*BoxX),(Y0+0.5*BoxY),(Y1+0.5*BoxY),num)
		else:
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

			if (jj==0):
				if not HORIZONTAL:
					cbar=plt.colorbar(im, cax=cax,orientation='horizontal',ticklocation='top')
				else:
					cbar=plt.colorbar(im, cax=cax,orientation='vertical',ticklocation='right')
				cbar.set_clim(-5,1)
				cbar.set_ticks([-5,-4,-3,-2,-1,0,1])
				cbar.set_label(r"log$_{10}$($\Sigma/\Sigma_0$)",fontsize=16,labelpad=10)
				for label in cbar.ax.get_xticklabels(): label.set_fontsize(11)
					
		#use time stamp to identify the location of the primary and secondary of the binary
		if ("innerboundary" in run_base):
			ecc_anom = orbital.KeplerEquation(time + np.pi,eb)
			posx1, posy1 = -qb/(1+qb)*(np.cos(ecc_anom) - eb),-qb/(1+qb)*(np.sqrt(1 - eb * eb) * np.sin(ecc_anom))
			posx2, posy2 = 1.0/(1+qb)*(np.cos(ecc_anom) - eb),1.0/(1+qb)*(np.sqrt(1 - eb * eb) * np.sin(ecc_anom))
			ax.plot([posx1,posx2],[posy1,posy2],'b.',ms=8.0,zorder=100)
			
		if (resolution_factor[jj] > 1):
			ax.text(0.05,0.07,r"$m_\mathrm{gas}/%i\;\;r_\mathrm{acc}=%.3fa_b$" % (resolution_factor[jj],accretion_radius[jj]),
				size=18,color='white',transform=ax.transAxes)
		
		if not (HORIZONTAL):
			if (jj== len(output_base_label_list)-1):
				ax.set_xlabel(r"$x/a_b$",size=22)
			else:
				ax.set_xticklabels([])
			ax.set_ylabel(r"$y/a_b$",size=22,labelpad=-5)
		else:
			if (jj== 0):
				ax.set_ylabel(r"$y/a_b$",size=22,labelpad=-2)
			else:
				ax.set_yticklabels([])
			#ax.xaxis.tick_top()
			ax.xaxis.set_tick_params(labeltop="on",labelbottom="off")
			ax.xaxis.set_label_position("top")
			ax.set_xlabel(r"$x/a_b$",size=22,labelpad=7)	
			
		ax.set_xlim(X0,X1)
		ax.set_ylim(Y0,Y1)
		ax.locator_params(nbins = 6)
		ax.tick_params('both',length=8,width=1.2,which='major')
		xticks=ax.get_xticks()
		yticks=ax.get_yticks()
		
		ax.set_xticks(xticks[1:-1])
		ax.set_yticks(yticks[1:-1])

	if HORIZONTAL:
		figname="./figure-resolution-comparison-panels_q%.1f_e%.1f_horizontal" % (qb,eb)
	else:
		figname="./figure-resolution-comparison-panels_q%.1f_e%.1f" % (qb,eb)
	print "Saving figure",figname+".pdf"
	plt.savefig(figname+".pdf")
	plt.savefig(figname+".png",dpi=200)
	plt.clf()

	exit()
	
