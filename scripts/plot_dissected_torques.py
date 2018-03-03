import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from disk_data_analysis.plotting import plot_slice, ImageData
import matplotlib.cm as cm
import sys
import os

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"



# Diego J. Munoz
# 2017

'''
Script to plot the spatially dissected torques from
circumbinary disk simulations

'''

READ_SNAP = False

# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out


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
    print "Reading output from simulation run:", run_name
    directory = simset.run_path+run_name
    snap_path = directory + simset.base + simset.snap_base

    #check the number of orbits at the zeroth and last snapshots
    orbit_range = []
    for snapnum in snap_list:
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS'],parttype=0)
        time = snap.header.time
        orbit_range.append(int(time/(2*np.pi)))
    orbit_init = orbit_range[0]
    orbit_final = orbit_range[-1]

    file_base = 'binary_dissected_torques_'
    if (eb == 0):
        run_tag = 'norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f' % (orbit_init,orbit_final,qb,eb,h0,alpha)
    elif (np.floor(np.log10(eb)) == -1):
        run_tag = 'norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f' % (orbit_init,orbit_final,qb,eb,h0,alpha)
    elif (np.floor(np.log10(eb)) == -2):
        run_tag = 'norbits%i-%i_q%.1f_e%.2f_h%.1f_alpha%.1f' % (orbit_init,orbit_final,qb,eb,h0,alpha)

    torque_filename = file_base+run_tag+'.txt'

    print "Looking for file:",torque_filename

    if (os.path.isfile(torque_filename)):
        print "Reading file:",torque_filename
    else:
        print "File does not exist."
        exit()

    torque_data = np.loadtxt(torque_filename)
    t = torque_data[:,0]
    torque0 = torque_data[:,1]
    torque1 = torque_data[:,2]
    torque2 = torque_data[:,3]
    torquetotal = torque_data[:,4]


    fig  = plt.figure(figsize=(8,10))

    ax = fig.add_axes([0.13,0.8,0.86,0.195])
    ax.plot(t/2/np.pi,torque0,color='salmon',marker='o',ms=2.0)
    ax.text(0.82,0.85,r'$\langle T\rangle=%.3f$' % torque0.mean(),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylabel(r'$T_{\rm b,grav}^{\rm (in)}\,/\,(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$',size=19)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)
        
    ax = fig.add_axes([0.13,0.57,0.86,0.195])
    ax.plot(t/2/np.pi,torque1,color='cornflowerblue') 
    ax.text(0.82,0.85,r'$\langle T\rangle=%.3f$' % torque1.mean(),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylabel(r'$T_{\rm b,grav}^{\rm (cav)}\,/\,(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$',size=19)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    ax = fig.add_axes([0.13,0.34,0.86,0.195])
    ax.plot(t/2/np.pi,torque2,color='gold') 
    ax.text(0.82,0.85,r'$\langle T\rangle=%.3f$' % torque2.mean(),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=22,labelpad=0)
    ax.set_ylabel(r'$T_{\rm b,grav}^{\rm (out)}\,/\,(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$',size=19)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    ax = fig.add_axes([0.13,0.06,0.86,0.195])
    [i.set_linewidth(2.0) for i in ax.spines.itervalues()]
    ax.plot(t/2/np.pi,torquetotal,color='navy',lw=1.5) 
    ax.text(0.82,0.85,r'$\langle T\rangle=%.3f$' % torquetotal.mean(),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=22)
    ax.set_ylabel(r'$T_{\rm b,grav}\,/\,(\dot{M}_0a_{\rm b}^2\Omega_{\rm b})$'
                  ,size=19)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    fig.savefig('./dissected_torques_'+run_tag+'.png',dpi=600)
    fig.clf()
