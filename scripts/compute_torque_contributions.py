import matplotlib.pyplot as plt
import numpy as np
import disk_data_analysis.circumbinary as dda
import disk_data_analysis.orbital as orbit_in_time
import sys

#Global variables
#  Normalization of the accretion rate
Mdot0 = 8.63857932377e-06
#  Binary properties
eb, qb = 0.0, 1.0
mu = qb / (1.0 + qb)
#  Simulation parameters
eta,h,alpha = 1.0, 0.1, 0.1

run_directory = '../restart-3000-pluto-woboundary-standardres-binary_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/' % (qb,eb,h,alpha,eta)

output_directory = 'output_restart_short/'

snap_base = 'snap_'

if __name__ == "__main__":

    snap_init = int(sys.argv[1])
    snap_final = int(sys.argv[2])

    for snap in range(snap_init,snap_final):
        
        snap_files = run_directory+output_directory + snap_base
        snap = ddf.get_snapshot_data(snap_files,snap_num,['POS','MASS'])

        # Compute forces
