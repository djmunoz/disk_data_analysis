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

###################################################################

time_offset=0

# Normalization of the accretion rate
Mdot0 = simset.Mdot0
# Normalization of the outer density
Sigma_out = simset.Sigma_out

GM  = 1.0
Rmax = simset.Rmax

if __name__ == '__main__':
######################################################################
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
    
    mu = qb / (1.0 + qb)



    #Check the number of orbits at the zeroth and last snapshots
    orbit_range = []
    for snapnum in [init_snap,final_snap]:
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS'],parttype=0)
        time = snap.header.time + time_offset
        orbit_range.append(int(time/(2*np.pi)))
    norbits = orbit_range[0]



    # Prepare grid for re-mapping of primitive variables

    NR, Nphi = 500, 600
    grid = dda.grid_polar(NR = NR, Nphi = Nphi,Rmin=1.0,Rmax= 80.0,scale='log')
    grid.X, grid.Y = grid.X + 80.0, grid.Y  +  80.0
    
    def nu(R): return (alpha * h0**2 * np.sqrt(GM) * R**(0.5))
    
    nu_grid = nu(grid.R)

    
    # open a file
    outfilename = 'jdot_balance_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f.txt' % (qb,eb,h0,alpha,eta)

    f = open(outfilename,'w')
    
    Rmin,Rmax = 0.75, 50.0
    for snapnum in snap_list:
        print "SNAPSHOT #",snapnum
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS','VELX','VELY','RHO','ACCE','R'],parttype=0)
        radii, mdot, jdot_adv, jdot_visc, torque_grav = dda.compute_angular_momentum_transfer(snap,Rmin,Rmax,NR=300)
        
        if(snapnum == snap_list[0]):
            f.write("time\t type\t\t radii\n")
            f.write("-\t-\t"+' '.join(radii.astype('|S7'))+"\n")


        f.write("%12.6f 0\t" % (snap.header.time))
        f.write(' '.join(mdot.astype('|S').tolist())+"\n")

        f.write("%12.6f 1\t" % (snap.header.time))
        f.write(' '.join(jdot_adv.astype('|S').tolist())+"\n")

        f.write("%12.6f 2\t" % (snap.header.time))
        f.write(' '.join(jdot_visc.astype('|S').tolist())+"\n")
        
        f.write("%12.6f 3\t" % (snap.header.time))
        f.write(' '.join(torque_grav.astype('|S').tolist())+"\n")

    
    print "Saved angular momentum transfer rate data to file:",outfilename
    f.close()
    

    
