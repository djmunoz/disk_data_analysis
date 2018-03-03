import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys


import simulation_settings as simset

'''
Script to compute the total torque and different torque contributions
directly from simulation snapshots

'''

    

if __name__ == "__main__":

    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h=simset.h
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

    
    run_name= simset.run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h,alpha,eta))
    
    print "Reading from simulation run:", run_name
    directory = simset.run_path+run_name
    snap_path = directory + simset.base + simset.snap_base
    
    mu = qb / (1.0 + qb)

    #check the number of orbits at the zeroth and last snapshots
    time_list = []
    for snapnum in snap_list:
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS'],parttype=0)
        time = snap.header.time
        time_list.append(time)
    orbit_range = (np.asarray(time_list)/2/np.pi).astype(int)
    orbit_init = orbit_range[0]
    orbit_final = orbit_range[-1]

    

    total_torque = []
    torque0 = []
    torque1 = []
    torque2 = []
    time_array = []

    outfilename1="binary_dissected_torques_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f.txt"\
        % (orbit_init,orbit_final,qb,eb,h,alpha)
    print "Will save dissected torque data to file:"
    print "                                  ",outfilename1


    for snapnum in snap_list:

        # Read snapshot information
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS','VELX','VELY','RHO','ACCE','R','ID','MASS'],parttype=0)
        # Position of the binary as a function of time
        time = snap.header.time
        X0, Y0 = 0.5 * snap.header.boxsize, 0.5 * snap.header.boxsize
        xb, yb, _, _ = orbit_in_time(time + np.pi, eb)
        time_array.append(time)
        pos1 = [mu * xb + X0, mu * yb + Y0, 0] # primary
        pos2 = [-(1 - mu) * xb + X0, -(1 - mu) * yb + Y0, 0] # secondary

        # Make sure we have accelerations due to each binary component
        # Compute forces

        print "SNAPSHOT #%i; time=%.6f " % (snapnum,time)
        '''

        accel1 = (1 - mu) * dda.compute_external_gravforce_from_snapshot(snap,XYZ = pos1,softening=0.026 * 2.8 * (1- mu))
        

        
        accel2 = mu * dda.compute_external_gravforce_from_snapshot(snap,XYZ=pos2,softening=0.026 * 2.8 * mu)


        # Compute the torque associated to each cell
        torque_per_cell = - snap.gas.MASS[:]/(1- mu) * (xb * accel1[:,1] - yb * accel1[:,0])\
                          + snap.gas.MASS[:]/mu * (xb * accel2[:,1] - yb * accel2[:,0])
        '''
    
        torque_per_cell = snap.gas.MASS[:] * (-(snap.gas.POS[:,0] - X0) * snap.gas.ACCE[:,1] +\
                                              (snap.gas.POS[:,1] - Y0) * snap.gas.ACCE[:,0])
        
        torque_per_cell /= simset.Mdot0

        # Separate regions
        dx, dy = snap.gas.POS[:,0] - X0, snap.gas.POS[:,1] - Y0
        #dx1, dy1 = snap.gas.POS[:,0] - pos1[0], snap.gas.POS[:,1] - pos1[1]
        #dx2, dy2 = snap.gas.POS[:,0] - pos2[0], snap.gas.POS[:,1] - pos2[1] 
        dr = np.sqrt(dx * dx  + dy * dy)
        #dr1 = np.sqrt(dx1 * dx1  + dy1 * dy1)
        #dr2 = np.sqrt(dx2 * dx2  + dy2 * dy2)

        # circum-single disks
        #csd_trunc = 0.23
        #cbd_trunc = 2.35
        #region0 = (dr1 < csd_trunc) | (dr2 < csd_trunc)
        #region2 = (dr > cbd_trunc) 
        #region1 = ((dr <= cbd_trunc)) & ((dr1 >= csd_trunc) & (dr2 >= csd_trunc))
        
        boundary1, boundary2 = 1.0, 2.0
        region0 = dr < boundary1
        region2 = dr > boundary2
        region1 = (dr <= boundary2) & (dr >= boundary1)
        

        torque0.append(torque_per_cell[region0].sum())
        torque1.append(torque_per_cell[region1].sum())
        torque2.append(torque_per_cell[region2].sum())

        total_torque.append(torque_per_cell.sum())

    time_array = np.asarray(time_array)
    total_torque = np.asarray(total_torque)
    torque0 = np.asarray(torque0)
    torque1 = np.asarray(torque1)
    torque2 = np.asarray(torque2)


    print "Saving dissected torque data to file:",outfilename1

    np.savetxt(outfilename1,np.array([time_array,torque0,torque1,torque2,total_torque]).T,
               fmt='%12f %.16e %.16e %.16e %.16e')




    print "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"

    
