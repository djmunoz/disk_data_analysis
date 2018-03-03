import numpy as np
import matplotlib.pyplot as plt
import sys

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"


# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out

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


def plot_binary_torque_contributions(binary_data,accretion_data,force_data,figfilename,
                                     qb=1.0,eb=0.0,orbit_init = None, orbit_final = None):
    '''
    Function to plot the different contributions to binary torque
    from a file of precomputed quantities.
    
    '''

    # Read quantities
    time = binary_data[:,0]
    x1 = binary_data[:,1]
    y1 = binary_data[:,2]
    x2 = binary_data[:,3]
    y2 = binary_data[:,4]

    mdot1 = accretion_data[:,0]
    mdot2 = accretion_data[:,1]

    fx1_grav = force_data[:,0]
    fy1_grav = force_data[:,1]
    fx2_grav = force_data[:,2]
    fy2_grav = force_data[:,3]
    
    fx1_acc = force_data[:,4]
    fy1_acc = force_data[:,5]
    fx2_acc = force_data[:,6]
    fy2_acc = force_data[:,7]

    # Combine data
    mdot = mdot1 + mdot2
    qdot = (1 + qb) * (mdot2 - qb * mdot1)
    fx1, fy1 = fx1_acc + fx1_grav,fy1_acc + fy1_grav
    fx2, fy2 = fx2_acc + fx2_grav,fy2_acc + fy2_grav
    
    torque_grav = compute_external_torques(x2 - x1,y2 - y1,fx2_grav - fx1_grav, fy2_grav - fy1_grav)
    torque_acc = compute_external_torques(x2 - x1,y2 - y1,fx2_acc - fx1_acc, fy2_acc - fy1_acc)
    
    reduced_mass = qb * 1.0 / (1 + qb) / (1 + qb)
    lb = np.sqrt(1 - eb**2)
    Jdot = reduced_mass * ((1.0 - qb)/(1.0 + qb) * qdot /qb * lb + mdot * lb + (torque_acc+torque_grav))
    
    if (orbit_init is None):
        orbit_init = time[0]/2 / np.pi
    if (orbit_final is None):
        orbit_final = time[-1]/2 / np.pi

    ind_interval = (time/2/np.pi >= orbit_init) & (time/2/np.pi <= orbit_final)

    fig  = plt.figure(figsize=(12,16))

    ax = fig.add_axes([0.08,0.83,0.86,0.16])
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    ax.plot(time/2/np.pi,mdot/np.abs(Mdot0),color='lightslategray',alpha=0.9,lw=1.5)
    ax.set_ylabel(r'$\dot{M}_{\rm b}/\dot{M}_0$',size=19)
    ax.text(0.85,0.85,r'$\langle \dot{M}\rangle=%.3f$' %  (mdot[ind_interval].mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=16,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    if ((mdot[ind_interval]/np.abs(Mdot0)).max() > 2.0):
        ax.set_ylim(0,(mdot[ind_interval]/np.abs(Mdot0)).max())  
    else:
        ax.set_ylim(0,2.0)  
    ticks = ax.get_yticks()
    ax.set_yticks(ticks[1:])
    [tick.label.set_fontsize(16) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #
    ax = fig.add_axes([0.08,0.64,0.87,0.16])
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    ax.plot(time/2/np.pi,qdot/np.abs(Mdot0),color='royalblue',alpha=0.9)
    ax.set_ylabel(r'$\dot{q}_{\rm b} / (\dot{M}_0/ M_{\rm b})$',size=19)
    ax.text(0.85,0.85,r'$\langle \dot{q}\rangle=%.3f$' %  (qdot[ind_interval].mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=16,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    if (np.abs(qdot[ind_interval]/np.abs(Mdot0)).max() > 3):
        ax.set_ylim(-5,2)
    else:
        ax.set_ylim(-2.5,2.5)
    ticks = ax.get_yticks()
    ax.set_yticks(ticks[1:])
    [tick.label.set_fontsize(16) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #
    ax = fig.add_axes([0.08,0.45,0.87,0.16])
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    ax.plot(time/2/np.pi,torque_grav/np.abs(Mdot0),color='darkblue')
    ax.set_ylabel(r'$\dot{l}_{\rm b,grav}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$',size=18)
    ax.text(0.85,0.85,r'$\langle \dot{l}\rangle=%.3f$' %  (torque_grav[ind_interval].mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=16,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(-30,35)
    ticks = ax.get_yticks()
    ax.set_yticks(ticks[1:])
    [tick.label.set_fontsize(16) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #
    ax = fig.add_axes([0.08,0.26,0.87,0.16])
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    ax.plot(time/2/np.pi,torque_acc/np.abs(Mdot0),color='steelblue')
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=26,labelpad=0)
    ax.set_ylabel(r'$\dot{l}_{\rm b,acc-ani}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$',size=18)
    ax.text(0.85,0.85,r'$\langle \dot{l}\rangle=%.3f$' %  (torque_acc[ind_interval].mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=16,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    if (np.abs(torque_acc[ind_interval]/np.abs(Mdot0)).max() > 4):
        ax.set_ylim(-25,20)
    else:
        ax.set_ylim(-5,4)
    ticks = ax.get_yticks()
    ax.set_yticks(ticks[1:])
    [tick.label.set_fontsize(16) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #
    ax = fig.add_axes([0.08,0.04,0.87,0.16])
    [i.set_linewidth(2.0) for i in ax.spines.itervalues()]
    ax.plot(time/2/np.pi,Jdot/np.abs(Mdot0),color='k')
    #ax.plot(time/2/np.pi,compute_external_torques(x2 - x1,y2 - y1,fx2 - fx1, fy2 - fy1))
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=26)
    ax.set_ylabel(r'$\dot{J}_{\rm b}/\left[(\dot{M}_0 a_{\rm b}^2\Omega_{\rm b}\right]$',size=19)
    ax.text(0.85,0.85,r'$\langle \dot{J}\rangle=%.3f$' %  (Jdot[ind_interval].mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=16,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_ylim(-8,8)
    ticks = ax.get_yticks()
    ax.set_yticks(ticks[1:])
    [tick.label.set_fontsize(16) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    
    # Saving figure
    print "Saving figure ",figfilename
    fig.savefig(figfilename,dpi=150)

    ##################################################
    print "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"

    


    return


if __name__ == "__main__":
    
    
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    orbit_init = int(sys.argv[3])
    orbit_final = int(sys.argv[4])

    qb = float(filename1.split('_q')[1].split('_')[0])
    eb = float(filename1.split('_e')[1].split('_')[0])
    figfilename = './torque_evolution_norbits%i-%i_q%.1f_e%.2f.png' % (orbit_init,orbit_final,qb,eb)


    time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_g,fy1_g,fx2_g,fy2_g,fx1_a,fy1_a,fx2_a,fy2_a = read_binary_forcing_file(filename2)
    _, mdot1, mdot2 =  read_binary_accretion_file(filename1)    


    
    plot_binary_torque_contributions(np.array([time,x1,y1,x2,y2]).T,
                                     np.array([mdot1,mdot2]).T,\
                                     np.array([fx1_g,fy1_g,fx2_g,fy2_g,fx1_a,fy1_a,fx2_a,fy2_a]).T,\
                                     figfilename, qb = qb, eb = eb, orbit_init = orbit_init,orbit_final=orbit_final)
