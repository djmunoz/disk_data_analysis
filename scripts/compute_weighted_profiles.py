import numpy as np
import sys
import disk_data_analysis.circumbinary as dda
import matplotlib.pyplot as plt

import simulation_settings as simset




# Diego J. Munoz
# 2017

'''
Script to compute radial profiles from disks simulations

'''

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
    snap_path = simset.find_simulation(qb,eb,h0,alpha,eta)
    mu = qb / (1.0 + qb)

    accretion_file = './mdot_profile.txt'
    radii_list = np.logspace(0,np.log10(20),100)
    delta_radii_list = np.diff(radii_list)
    radii_list = 0.5 * (radii_list[1:] + radii_list[:-1])
    f = open(accretion_file,'w+')
    f.write('-1111111   ')
    f.write('  '.join(radii_list.astype('|S').tolist())+"\n")

    for snapnum in snap_list:
        print snapnum
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS','MASS','VELR','RHO','R','ACCE','VELPHI','VOL','VELX','VELY'],parttype=0)
        time = snap.header.time
        X0, Y0 = 0.5 * snap.header.boxsize, 0.5 * snap.header.boxsize


        r1,r2 = np.meshgrid(snap.gas.R,radii_list)
        kernel_fac = np.exp(-(r1-r2)**2/delta_radii_list[:,None]**2/2)/np.sqrt(2*np.pi)/delta_radii_list[:,None]
	normalize = 1.0 / (kernel_fac * snap.gas.VOL / snap.gas.R).mean(axis=1)
        # Mass transfer
        massflux = snap.gas.R * snap.gas.RHO * snap.gas.VELR * 2 * np.pi
        mdot = (massflux[None,:] * kernel_fac * snap.gas.VOL / snap.gas.R).mean(axis=1) * normalize

        # Angular momentum transfer
        # ...gravity
        angmomflux_grav =  snap.gas.MASS[:] * ((snap.gas.POS[:,0] - X0) * snap.gas.ACCE[:,1]  \
                                               
                                              -(snap.gas.POS[:,1] - Y0) * snap.gas.ACCE[:,0])
        ii = np.argsort(snap.gas.R)
        aflux,r2 = np.meshgrid(angmomflux_grav,radii_list)
        aflux[r1 < r2] = 0
        jdot_grav = aflux.sum(axis=1)
        # ...advection
        angmomflux_adv =  -2 * np.pi * snap.gas.R[:] * (snap.gas.RHO[:] * snap.gas.R[:] * snap.gas.VELR[:] * snap.gas.VELPHI[:]) 
        jdot_adv = (angmomflux_adv[None,:] * kernel_fac * snap.gas.VOL / snap.gas.R).mean(axis=1) * normalize
        # ...viscosity
        gradientvx = dda.compute_snapshot_gradient(snap,'VELX')
        gradientvy = dda.compute_snapshot_gradient(snap,'VELY')
        snap.add_data(gradientvx,'GRVX')
        snap.add_data(gradientvy,'GRVY')
        GM = 1.0
        def nu(R): return alpha * h0**2 * np.sqrt(GM) * np.sqrt(R)
        angmomflux_visc =  -2 * np.pi * snap.gas.RHO[:] * nu(snap.gas.R[:]) * (2 * (snap.gas.POS[:,0] - X0) * (snap.gas.POS[:,1] - Y0) * \
                                                                               (snap.gas.GRVY[:,1] - snap.gas.GRVX[:,0]) + \
                                                                               ((snap.gas.POS[:,0] - X0)**2 - (snap.gas.POS[:,1] - Y0)**2) * \
                                                                               (snap.gas.GRVX[:,1] + snap.gas.GRVY[:,0]))
                                                                              
        jdot_visc = (angmomflux_visc[None,:] * kernel_fac * snap.gas.VOL / snap.gas.R).mean(axis=1) * normalize


        #indx = np.argsort(snap.gas['R'])
        f.write('%.8f   '% time)
        f.write(' '.join(mdot.astype('|S').tolist())+"\n")
        f.write('%.8f   '% time)
        f.write(' '.join(jdot_grav.astype('|S').tolist())+"\n")
        f.write('%.8f   '% time)
        f.write(' '.join(jdot_adv.astype('|S').tolist())+"\n")
        f.write('%.8f   '% time)
        f.write(' '.join(jdot_visc.astype('|S').tolist())+"\n")

    f.close()
