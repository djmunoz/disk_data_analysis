import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from disk_data_analysis.plotting import plot_slice, ImageData
import matplotlib.cm as cm
import sys

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


def compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,mdot1,mdot2,
                                  fx1_ext,fy1_ext,fx2_ext,fy2_ext,
                                  G = 1.0 ,Mb = 1.0, ab = 1.0, qb = 1.0, eb = 0,etype=0):
    """
    Compute the change in both orbital energy and specific angular
    momentum of a binary subject to mass changes and external 
    torques

    """
    lb = np.sqrt(G * Mb * ab * (1 - eb * eb))
    Eb = -G * Mb / 2 / ab

    x,y = x2 - x1, y2 - y1
    vx,vy = vx2 - vx1, vy2 - vy1
    fxext, fyext = fx2_ext - fx1_ext, fy2_ext - fy1_ext
    mdot = mdot1 + mdot2
    qdot = qb * (1 + qb) / Mb * (mdot2 / qb - mdot1)
    lb = x * vy - y * vx
    r = np.sqrt(x * x + y * y)
    Eb = 0.5 * (vx * vx + vy * vy) - G * Mb / r
    Jb = qb * Mb / (1 + qb)**2 * lb

    Edot = -G * mdot / r + (vx * fxext + vy * fyext)

    ldot = compute_external_torques(x,y,fxext,fyext)

    Jdot = Jb * ( (1.0 - qb)/(1.0 + qb) * qdot / qb + mdot/Mb + ldot / lb)
    


    sinphi, cosphi = y/r, x/r
    fr_ext = fxext * cosphi + fyext * sinphi
    fphi_ext = -fxext * sinphi + fyext * cosphi



    if (eb == 0):
        edot = np.sqrt(ab * (1.0 - eb * eb)/ G / Mb) * ((eb + 2 * cosphi + eb * cosphi * cosphi) / (1 + eb * cosphi) * fphi_ext +\
                                                        sinphi * fr_ext) - \
                                                        mdot / Mb * (eb + cosphi)
        adotovera = 2 * np.sqrt(ab / (1.0 - eb * eb)/ G / Mb) * (eb * sinphi * fr_ext + (1 + eb * cosphi) * fphi_ext) -\
                  mdot / Mb * (1 + 2 * eb * cosphi + eb * eb) / (1 - eb * eb)
    else:
        edot= 1.0/(G**2 * Mb**2 / lb**2 / Eb + 2) * np.sqrt(1 + 2 * lb**2 * Eb / G**2 / Mb**2) *  (2 * ldot / lb +  Edot / Eb - 2 * mdot / Mb)
        adotovera = ab * (-Edot / Eb + mdot / Mb)

    return Edot, ldot, Jdot, adotovera, edot

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

    if (len(sys.argv) < 5): 
        orbit_init = None
    else:   
        orbit_init = int(sys.argv[4])
    if (len(sys.argv) < 6): 
        orbit_final = None
    else:  
        orbit_final = int(sys.argv[5])



    
    run_name= simset.run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h0,alpha,eta))
    print "Reading from simulation run:", run_name
    
    mu = qb / (1.0 + qb) #mass ratio secondary-to-total
    reduced_mass = qb / (1 + qb)**2 # actual reduced mass


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
    mdot = mdot1 + mdot2
    qdot = (1 + qb) * (mdot2 - qb * mdot1)
    fx1, fy1 = fx1_a + fx1_g,fy1_a + fy1_g
    fx2, fy2 = fx2_a + fx2_g,fy2_a + fy2_g
    
    

    torque_g = compute_external_torques(x2 - x1,y2 - y1,fx2_g - fx1_g, fy2_g - fy1_g)
    torque_a = compute_external_torques(x2 - x1,y2 - y1,fx2_a - fx1_a, fy2_a - fy1_a)
    
    Jdot = reduced_mass * ((1.0 - qb)/(1.0 + qb) * qdot /qb + mdot + (torque_a+torque_g))

    

    Edot, ldot, Jdot0, adot, edot = compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,mdot1,mdot2,fx1,fy1,fx2,fy2,eb=eb)
    
    #separate the contributions to adot
    _, _, _, adot_ani, edot_ani = compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,0.0,0.0,fx1_a,fy1_a,fx2_a,fy2_a,eb=eb)
    _, _, _, adot_iso, edot_iso = compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,mdot1,mdot2,0.0,0.0,0.0,0.0,eb=eb)
    _, _, _, adot_grav,edot_grav = compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,0.0,0.0,fx1_g,fy1_g,fx2_g,fy2_g,eb=eb)



    ########################################################
    # Create figures
    fig  = plt.figure(figsize=(10,5))

    ax = fig.add_axes([0.09,0.58,0.85,0.4])
    ax.plot(time/2/np.pi,adot/np.abs(Mdot0),color='firebrick')
    ax.set_ylabel(r'$\dot{a}_{\rm b}/(a_{\rm b}\dot{M}_0/ M_{\rm b})$',size=19)
    ax.text(0.82,0.85,r'$\langle \dot{a}\rangle=%.3f$' % (adot.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)
    #ax.set_ylim(0,2.0)
    #
    ax = fig.add_axes([0.09,0.11,0.85,0.4])
    ax.plot(time/2/np.pi,edot/np.abs(Mdot0),color='navy')
    ax.set_ylabel(r'$\dot{e}_{\rm b} / (\dot{M}_0/ M_{\rm b})$',size=19)
    ax.text(0.82,0.85,r'$\langle \dot{e}\rangle=%.3f$' % (edot.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.set_xlim(orbit_init,orbit_final)
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=22,labelpad=0)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)
    #ax.set_ylim(-1.5,1.5)
    #

    
    # Saving figure
    figfilename = './orbital_evolution_q%.1f_e%.2f.pdf' % (qb,eb)
    print "Saving figure ",figfilename
    fig.savefig(figfilename)

    ##################################################
    fig  = plt.figure(figsize=(10,8))

    ax = fig.add_axes([0.11,0.7,0.85,0.25])
    ax.plot(time/2/np.pi,adot_iso/np.abs(Mdot0),color='red')
    ax.text(0.82,0.85,r'$\langle \dot{a}\rangle=%.3f$' % (adot_iso.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.text(0.02,0.9,'accretion only',size=14,transform=ax.transAxes)
    ax.set_ylabel(r'$\dot{a}_{\rm b}\,_{_{\rm acc}}\, /\,(a_{\rm b}\dot{M}_0/ M_{\rm b})$',size=19)
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #ax.set_ylim(0,2.0)
    #
    ax = fig.add_axes([0.11,0.4,0.85,0.25])
    ax.plot(time/2/np.pi,adot_ani/np.abs(Mdot0),color='salmon')
    ax.text(0.82,0.85,r'$\langle \dot{a}\rangle=%.3f$' % (adot_ani.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.text(0.02,0.9,'anisotropic accretion torque only',size=14,transform=ax.transAxes)
    ax.set_ylabel(r'$\dot{a}_{\rm b}\,_{_{\rm acc-ani}}\, /\, (a_{\rm b}\dot{M}_0/ M_{\rm b})$',size=19)
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #ax.set_ylim(-1.5,1.5)
    #
    ax = fig.add_axes([0.11,0.1,0.85,0.25])
    ax.plot(time/2/np.pi,adot_grav/np.abs(Mdot0),color='indianred')
    ax.set_ylabel(r'$\dot{a}_{\rm b}\,_{_{\rm grav}}\, /\,(a_{\rm b}\dot{M}_0/ M_{\rm b})$',size=19)
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=22)
    ax.text(0.82,0.85,r'$\langle \dot{a}\rangle=%.3f$' % (adot_grav.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.text(0.02,0.9,'gravitational torque only',size=14,transform=ax.transAxes)
    ax.set_ylim(-99,99)
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #ax.set_ylim(-1.5,1.5)
    #
    
    # Saving figure
    figfilename = './semimajor_axis_evolution_q%.1f_e%.2f.pdf' % (qb,eb)
    print "Saving figure ",figfilename
    fig.savefig(figfilename)

    ##################################################
    
    ##################################################
    fig  = plt.figure(figsize=(10,8))

    ax = fig.add_axes([0.11,0.7,0.85,0.25])
    ax.plot(time/2/np.pi,edot_iso/np.abs(Mdot0),color='blue')
    ax.text(0.82,0.85,r'$\langle \dot{e}\rangle=%.3f$' % (edot_iso.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.text(0.02,0.9,'accretion only',size=14,transform=ax.transAxes)
    ax.set_ylabel(r'$\dot{e}_{\rm b}\,_{_{\rm acc}}\, /\,(\dot{M}_0/ M_{\rm b})$',size=19)
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #ax.set_ylim(0,2.0)
    #
    ax = fig.add_axes([0.11,0.4,0.85,0.25])
    ax.plot(time/2/np.pi,edot_ani/np.abs(Mdot0),color='cornflowerblue')
    ax.text(0.82,0.85,r'$\langle \dot{e}\rangle=%.3f$' % (edot_ani.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.text(0.02,0.9,'anisotropic accretion torque only',size=14,transform=ax.transAxes)
    ax.set_ylabel(r'$\dot{e}_{\rm b}\,_{_{\rm acc-ani}}\, /\, (\dot{M}_0/ M_{\rm b})$',size=19)
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #ax.set_ylim(-1.5,1.5)
    #
    ax = fig.add_axes([0.11,0.1,0.85,0.25])
    ax.plot(time/2/np.pi,edot_grav/np.abs(Mdot0),color='slateblue')
    ax.set_ylabel(r'$\dot{e}_{\rm b}\,_{_{\rm grav}}\, /\,(\dot{M}_0/ M_{\rm b})$',size=19)
    ax.set_xlabel(r'$t\,/\,P_{\rm b}$',size=22)
    ax.text(0.82,0.85,r'$\langle \dot{e}\rangle=%.3f$' % (edot_grav.mean()/np.abs(Mdot0)),
            transform=ax.transAxes,size=14,
            bbox=dict(facecolor='w',edgecolor='k'))
    ax.text(0.02,0.9,'gravitational torque only',size=14,transform=ax.transAxes)
    ax.set_ylim(-99,99)
    ax.set_xlim(orbit_init,orbit_final)
    [tick.label.set_fontsize(12) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(12) for tick in ax.yaxis.get_major_ticks()]
    ax.tick_params(axis='both',which='major',direction='in',length=6)

    #ax.set_ylim(-1.5,1.5)
    #
    
    # Saving figure
    figfilename = './eccentricity_evolution_q%.1f_e%.2f.pdf' % (qb,eb)
    print "Saving figure ",figfilename
    fig.savefig(figfilename)

    ##################################################
