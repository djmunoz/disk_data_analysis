import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from disk_data_analysis.plotting import plot_slice, ImageData
import matplotlib.cm as cm
import sys
from plot_binary_torque_contributions import*

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"



# Diego J. Munoz
# 2017

'''
Script to compute the external torque acting directly on
an accreting binary from simulation output in the form of
snapshot files and ASCII files

'''

READ_SNAP = False

# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out

# Local functions
def read_binary_accretion_file(filename):
    """
    Read from disk a precomputed file with the accretion rates
    onto each component of the binary

    """
    accretion_data = np.loadtxt(filename)
    time = accretion_data[:,0]
    mdot1 = accretion_data[:,1]
    mdot2 = accretion_data[:,2]

    return time, mdot1, mdot2
    
def read_binary_forcing_file(filename,gravity_force = True, accretion_force = True):
    
    """ 
    Read from disk a precomputed file with the external forces
    acting on the binary

    """

    force_data = np.loadtxt(filename)
    time = force_data[:,0]
    x2 = force_data[:,1]
    y2 = force_data[:,2]
    x1 = force_data[:,3]
    y1 = force_data[:,4]
    vx2 = force_data[:,5]
    vy2 = force_data[:,6]
    vx1 = force_data[:,7]
    vy1 = force_data[:,8]
    fx2_acc = force_data[:,9]
    fy2_acc = force_data[:,10]
    fx1_acc = force_data[:,11]
    fy1_acc = force_data[:,12]
    fx2_grav = force_data[:,13]
    fy2_grav = force_data[:,14]
    fx1_grav = force_data[:,15]
    fy1_grav = force_data[:,16]

    return time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_grav,fy1_grav,fx2_grav,fy2_grav,fx1_acc,fy1_acc,fx2_acc,fy2_acc

def compute_external_torques(x,y,fx,fy):
    """
    Compute torque due to external forces on a binary orbit given
    the position of both elements of the binary and the forces acting
    on each of those two elements

    x,y: numpy arrays -- (x,y) position of the binary relative coordinate

    fx,fy: numpy arrays -- x- and y-components of the net external force 
    **per unit mass** acting on the binary

    """
    
    t =  x * fy - y * fx
    
    return t


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
    
    mu = qb / (1.0 + qb) #mass ratio secondary-to-total
    reduced_mass = qb / (1 + qb)**2 # actual reduced mass
    print reduced_mass
    
    #check the number of orbits at the zeroth and last snapshots
    orbit_range = []
    for snapnum in snap_list:
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS'],parttype=0)
        time = snap.header.time
        orbit_range.append(int(time/(2*np.pi)))
    orbit_init = orbit_range[0]
    orbit_final = orbit_range[-1]

    # Go through the snapshots:
    total_torque_disk = []
    total_torque_bin = []
    time_array = []
    if (READ_SNAP):
        #for snapnum in snap_list:
        for snapnum in range(0,3001):

            # Read snapshot information
            snap = dda.get_snapshot_data(snap_path,snapnum,['POS','VELX','VELY','RHO','ACCE','R','ID','MASS'],parttype=0)
            
            # Make sure we have accelerations due to each binary component
            # Compute forces
            time = snap.header.time
            print "SNAPSHOT #%i; time=%.6f " % (snapnum,time)
            time_array.append(time)
            X0, Y0 = 0.5 * snap.header.boxsize, 0.5 * snap.header.boxsize
            
            xb, yb, _, _ = orbit_in_time(time + np.pi, eb)
            
            pos1 = [mu * xb + X0, mu * yb + Y0, 0] # primary
            accel1 = (1 - mu) * dda.compute_external_gravforce_from_snapshot(snap,XYZ = pos1,softening=0.026 * 2.8 * (1- mu))
            
            pos2 = [-(1 - mu) * xb + X0, -(1 - mu) * yb + Y0, 0] # secondary
            accel2 = mu * dda.compute_external_gravforce_from_snapshot(snap,XYZ=pos2,softening=0.026 * 2.8 * mu)
            
            
            # Compute the torque associated to each cell
            torque_bin_per_cell = mu * (1.0 -mu) *(- snap.gas.MASS[:]/(1- mu) * (xb * accel1[:,1] - yb * accel1[:,0])\
                                                   + snap.gas.MASS[:]/mu * (xb * accel2[:,1] - yb * accel2[:,0]))
            
            torque_disk_per_cell = snap.gas.MASS[:] * (+(snap.gas.POS[:,0] - X0) * snap.gas.ACCE[:,1]  \
                                                       
                                                       -(snap.gas.POS[:,1] - Y0) * snap.gas.ACCE[:,0])
            
            ind_select = (snap.gas.R >= 1.0) & (snap.gas.R < 70)
            total_torque_bin.append(torque_bin_per_cell[ind_select].sum())
            total_torque_disk.append(torque_disk_per_cell[ind_select].sum())

        total_torque_disk = np.asarray(total_torque_disk)
        total_torque_bin = np.asarray(total_torque_bin)

    acc_file_base = 'binary-accretion-rate_'
    tor_file_base = 'binary-forcing-rate_'
    if (eb == 0):
        run_tag = 'norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f.txt' % (orbit_init,orbit_final,qb,eb,h0,alpha)
    elif (np.floor(np.log10(eb)) == -1):
        run_tag = 'norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f.txt' % (orbit_init,orbit_final,qb,eb,h0,alpha)
    elif (np.floor(np.log10(eb)) == -2):
        run_tag = 'norbits%i-%i_q%.1f_e%.2f_h%.1f_alpha%.1f.txt' % (orbit_init,orbit_final,qb,eb,h0,alpha)

    accretion_filename = acc_file_base+run_tag
    torque_filename = tor_file_base+run_tag

    print "Reading files..."
    print "                ",accretion_filename
    print "                ",torque_filename

    
    time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_g,fy1_g,fx2_g,fy2_g,fx1_a,fy1_a,fx2_a,fy2_a = read_binary_forcing_file(torque_filename)
    _, mdot1, mdot2 =  read_binary_accretion_file(accretion_filename)    
    figfilename = './torque_evolution_q%.1f_e%.2f.png' % (qb,eb)
    
    plot_binary_torque_contributions(np.array([time,x1,y1,x2,y2]).T,
                                     np.array([mdot1,mdot2]).T,\
                                     np.array([fx1_g,fy1_g,fx2_g,fy2_g,fx1_a,fy1_a,fx2_a,fy2_a]).T,\
                                     figfilename, qb = qb, eb = eb)

    exit()

    mdot = mdot1 + mdot2
    qdot = (1 + qb) * (mdot2 - qb * mdot1)
    fx1, fy1 = fx1_a + fx1_g,fy1_a + fy1_g
    fx2, fy2 = fx2_a + fx2_g,fy2_a + fy2_g
    
    torque_g = compute_external_torques(x2 - x1,y2 - y1,fx2_g - fx1_g, fy2_g - fy1_g)
    torque_a = compute_external_torques(x2 - x1,y2 - y1,fx2_a - fx1_a, fy2_a - fy1_a)
    
    Jdot = reduced_mass * ((1.0 - qb)/(1.0 + qb) * qdot /qb + mdot + (torque_a+torque_g))

    fig  = plt.figure(figsize=(12,12))
    if (READ_SNAP): npanels = 6
    else : npanels = 5

    ax = fig.add_subplot(npanels,1,1)
    ax.plot(time/2/np.pi,mdot/np.abs(Mdot0),color='k')
    ax.set_ylabel(r'$\dot{M}_{\rm b}/\dot{M}_0$')
    ax.text(0.85,0.05,'mean=%.3f' % (mdot.mean()/np.abs(Mdot0)),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(0,2.0)
    #
    ax = fig.add_subplot(npanels,1,2)
    ax.plot(time/2/np.pi,qdot/np.abs(Mdot0),color='k')
    ax.set_ylabel(r'$\dot{q}_{\rm b} / (\dot{M}_0/ M_{\rm b})$')
    ax.text(0.85,0.05,'mean=%.3f' % (qdot.mean()/np.abs(Mdot0)),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(-1.5,1.5)
    #
    ax = fig.add_subplot(npanels,1,3)
    ax.plot(time/2/np.pi,torque_g/np.abs(Mdot0),color='steelblue')
    ax.set_ylabel(r'$\dot{l}_{\rm b,grav}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$')
    ax.text(0.85,0.05,'mean=%.3f' % (torque_g.mean()/np.abs(Mdot0)),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(-30,30)
    #
    ax = fig.add_subplot(npanels,1,4)
    if READ_SNAP:

        ax.plot(np.asarray(time_array)/2/np.pi,total_torque_bin/np.abs(Mdot0),color='steelblue')
        ax.plot(np.asarray(time_array)/2/np.pi,-total_torque_disk/np.abs(Mdot0),color='pink')
        ax.set_ylabel(r'$\dot{l}_{\rm b,grav}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$')
        ax.text(0.85,0.05,'mean=%.3f' % (total_torque_bin/np.abs(Mdot0)).mean(),transform=ax.transAxes,color='steelblue')
        ax.text(0.65,0.05,'mean=%.3f' % (-total_torque_disk/np.abs(Mdot0)).mean(),transform=ax.transAxes,color='pink')
        ax.set_xlim(orbit_init,orbit_final)
        ax.set_ylim(-6,6)
        ax = fig.add_subplot(npanels,1,5)
    #
    ax.plot(time/2/np.pi,torque_a/np.abs(Mdot0),color='royalblue')
    ax.set_ylabel(r'$\dot{l}_{\rm b,acc-ani}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$')
    ax.text(0.85,0.05,'mean=%.3f' % (torque_a.mean()/np.abs(Mdot0)),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(-3,3)
    #
    ax = fig.add_subplot(npanels,1,npanels)
    ax.plot(time/2/np.pi,Jdot/np.abs(Mdot0),color='darkblue')
    #ax.plot(time/2/np.pi,compute_external_torques(x2 - x1,y2 - y1,fx2 - fx1, fy2 - fy1))
    ax.set_xlabel(r'$t/P_{\rm b}$',size=18)
    ax.set_ylabel(r'$\dot{J}_{\rm b}/\left[(\dot{M}_0 a_{\rm b}^2\Omega_{\rm b}\right]$')
    ax.text(0.85,0.05,'mean=%.3f' % (Jdot.mean()/np.abs(Mdot0)),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)

    
    # Saving figure
    figfilename = './torque_evolution_q%.1f_e%.2f.png' % (qb,eb)
    print "Saving figure ",figfilename
    fig.savefig(figfilename)

    ##################################################
    
    fig  = plt.figure(figsize=(12,5))
    ax = fig.add_subplot(211)
    ax.plot(time/2/np.pi,mdot/np.abs(Mdot0),color='k')
    ax.set_ylabel(r'$\dot{M}_{\rm b}/\dot{M}_0$')
    ax.text(0.85,0.05,'mean=%.3f' % (mdot/np.abs(Mdot0)).mean(),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(0,2.0)

    #
    ax = fig.add_subplot(212)
    ax.plot(time/2/np.pi,Jdot/np.abs(Mdot0),color='k')
    ax.set_ylabel(r'$\dot{J}_{\rm b}/(\dot{M}_0 a_{\rm b}\Omega_{\rm b})$')
    ax.text(0.85,0.05,'mean=%.3f' % (Jdot/np.abs(Mdot0)).mean(),transform=ax.transAxes)
    ax.set_xlim(orbit_init,orbit_final)
    #ax.set_ylim(0,2.0)
    # Saving figure

    figfilename = './torque_eigenvalue_q%.1f_e%.2f.png' % (qb,eb)
    print "Saving figure ",figfilename
    fig.savefig(figfilename)

    print "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"

    
