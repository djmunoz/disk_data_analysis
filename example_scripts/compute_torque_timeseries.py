import numpy as np
import matplotlib.pyplot as plt
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

'''
Script to compute the total torque and different torque contributions
directly from simulation snapshots

'''

# Binary properties
eb, qb = 0.0 , 1.0 # eccentricity and mass ratio
mu = qb / (1.0 + qb)

# Normalization of the accretion rate
Mdot0 = 8.63857932377e-06
# Normalization of the outer density
Sigma_out = 1.0238548710e-4
    

if __name__ == "__main__":

    snapinit = int(sys.argv[1])
    snapfinal = int(sys.argv[2])


    for snapnum in range(snapinit,snapfinal):

        # Read snapshot information
        snap = dda.get_snapshot_data('../data/snap_',snapnum,['POS','VELX','VELY','RHO','ACCE','R','ID','MASS'],parttype=0)

        # Make sure we have accelerations due to each binary component
        # Compute forces
        time = snap.header.time
        X0, Y0 = 0.5 * snap.header.boxsize, 0.5 * snap.header.boxsize
        xb, yb, _, _ = orbit_in_time(time + np.pi, eb)
        pos1 = [mu * xb + X0, mu * yb + Y0, 0] # primary
        accel1 = (1 - mu) * dda.compute_external_gravforce_from_snapshot(snap,XYZ = pos1,softening=0.026 * 2.8 * (1- mu))
        
        pos2 = [-(1 - mu) * xb + X0, -(1 - mu) * yb + Y0, 0] # secondary
        accel2 = mu * dda.compute_external_gravforce_from_snapshot(snap,XYZ=pos2,softening=0.026 * 2.8 * mu)


        # Compute the torque associated to each cell
        torque_per_cell = - snap.gas.MASS[:]/(1- mu) * (xb * accel1[:,1] - yb * accel1[:,0])\
                          + snap.gas.MASS[:]/mu * (xb * accel2[:,1] - yb * accel2[:,0])
        
        torque_per_cell /= Mdot0

        # Separate regions
        dx, dy = snap.gas.POS[:,0] - X0, snap.gas.POS[:,1] - Y0
        dx1, dy1 = snap.gas.POS[:,0] - pos1[0], snap.gas.POS[:,1] - pos1[1]
        dx2, dy2 = snap.gas.POS[:,0] - pos2[0], snap.gas.POS[:,1] - pos2[1] 
        dr = np.sqrt(dx * dx  + dy * dy)
        dr1 = np.sqrt(dx1 * dx1  + dy1 * dy1)
        dr2 = np.sqrt(dx2 * dx2  + dy2 * dy2)
        # circum-single disks
        csd_trunc = 0.23
        cbd_trunc = 2.35
        csd_region = (dr1 < csd_trunc) | (dr2 < csd_trunc)
        cbd_region = (dr > cbd_trunc) 
        str_region = ((dr <= cbd_trunc)) & ((dr1 >= csd_trunc) & (dr2 >= csd_trunc))
        
        print "\n", torque_per_cell[csd_region].sum()
        print torque_per_cell[cbd_region].sum()
        print torque_per_cell[str_region].sum()
        print torque_per_cell.sum()
