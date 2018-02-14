# Script to compute the angular momentum balance as a function of radius
# in 2-D disk simulations

# Diego J. Munoz
# 2017
##############################
import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
import readsnapHDF5 as rs
import glob
from string import split
import sys

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"


'''
This script computes the angular momentum transfer rates (i.e., torques)
from 3 different sources: gravitational, viscous, advective stored
in a precomputed file. The file itself is produced by a companion
script called compute_angular_momentum_balance.py
'''


# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out



if __name__ == '__main__':

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

    # open a file
    #infilename = 'jdot_balance_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f.txt' % (qb,eb,h0,alpha,eta)
    infilename = 'jdot_balance.txt'

    print "Opening file...",infilename
    f = open(infilename,'r')
    # read the first line
    f.readline()
    # read the second line
    radii = np.asarray(f.readline().split()[2:]).astype(float)
    f.close()
    # read the rest of the file
    data = np.loadtxt(infilename,skiprows=2)
    ind_mdot = data[:,1]== 0
    ind_adv = data[:,1]== 1
    ind_visc = data[:,1]== 2
    ind_grav = data[:,1]== 3

    time = data[ind_mdot,0]
    mdot = data[ind_mdot,2:]
    jdot_adv = data[ind_adv,2:]
    jdot_visc = data[ind_visc,2:]
    jdot_grav = data[ind_grav,2:]

    print mdot[:,(radii < 20) & (radii > 8)].mean(axis=1).mean(),Mdot0
    Mdot0 = mdot[:,(radii < 20) & (radii > 8)].mean(axis=1).mean()
    print Mdot0
    # Prepare figure
    fig = plt.figure(1,figsize=(12,7))
    fig.subplots_adjust(hspace=0.7,wspace=0.3,top=0.98,left=0.07,right=0.98)

    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
    ax =  fig.add_axes([0.15,0.1,0.7,0.45])
    
    for k,t in enumerate(time):
        ax1.plot(radii,jdot_adv[k,:] / Mdot0,color='cornflowerblue',alpha=0.4,lw=0.7)
        ax2.plot(radii,jdot_visc[k,:] / Mdot0,color='forestgreen',alpha=0.4,lw=0.7)
        ax3.plot(radii,jdot_grav[k,:]/ Mdot0,color='orange',alpha=0.4,lw=0.7)
    
    ax1.set_xlim(0,16)
    ax2.set_xlim(0,16)
    ax3.set_xlim(0,16)
    
    ax1.set_ylim(-100,50)
    ax2.set_ylim(-0.1,4)
    ax3.set_ylim(-2.5,2.5)

    ax1.set_xlabel(r'$R/a_{\rm b}$',size=18)
    ax2.set_xlabel(r'$R/a_{\rm b}$',size=18)
    ax3.set_xlabel(r'$R/a_{\rm b}$',size=18)

    ax1.set_ylabel(r'$\dot{J}_{\rm adv}/(\dot{M}_0a^2_{\rm b}\Omega_{\rm b})$',size=18,labelpad=-2)
    ax2.set_ylabel(r'$\dot{J}_{\rm visc}/(\dot{M}_0a^2_{\rm b}\Omega_{\rm b})$',size=18,labelpad=2)
    ax3.set_ylabel(r'$\dot{T}_{\rm grav}^{>R}/(\dot{M}_0a^2_{\rm b}\Omega_{\rm b})$',size=18,labelpad=0)
    
    ax1.text(0.4,0.9,r'%i snapshots' % (jdot_adv.shape[0]),transform = ax1.transAxes)
    ax1.text(0.4,0.8,r'$t\in [%iP_{\rm b},%iP_{\rm b}]$' % (time[0]/2/np.pi,time[-1]/2/np.pi),\
             transform = ax1.transAxes,size=14)

    ticks = [0,4,8,12,16]
    ax1.set_xticks(ticks)
    ax2.set_xticks(ticks)
    ax3.set_xticks(ticks)

    label_adv = r'$\langle\dot{J}_{\rm adv}\rangle_T$'
    label_visc = r'$\langle\dot{J}_{\rm visc}\rangle_T$'
    label_grav = r'$\langle{T}_{\rm grav}^{>R}\rangle_T$'
    ax.plot(radii,jdot_adv.mean(axis=0) / Mdot0,color='cornflowerblue',label=label_adv)
    ax.plot(radii,jdot_visc.mean(axis=0) / Mdot0,color='forestgreen',label=label_visc)
    ax.plot(radii,jdot_grav.mean(axis=0) / Mdot0,color='orange',label=label_grav)
    ax.plot(radii,(jdot_adv-jdot_visc-jdot_grav).mean(axis=0)/ Mdot0,\
            color='firebrick',label=None)

    ax.plot([0,80],[0,0],'k:',lw=0.5,zorder=0)

    ax.legend(loc='upper left',frameon=False,prop={'size':16})

    ax.set_xlim(0,10)
    ax.set_xticks([0,1,2,4,6,8,10])
    ax.set_ylim(-0.5,3)
    ax.set_xlabel(r'$R/a_{\rm b}$',size=24)
    ax.set_ylabel(r'$\langle\dot{J}\rangle_{T}\,/\,(\dot{M}_0a^2_{\rm b}\Omega_{\rm b})$',size=24)



    figfilename = './jdot_balance.pdf'
    print 'Saving figure...',figfilename
    fig.savefig(figfilename)
    fig.clf()
