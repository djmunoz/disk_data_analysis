# paths to files
run_path = "/1TB/dmunoz/CIRCUMBINARY_DISKS_2D/RESTART_PLUTO_RUNS/"
run_base ="restart-3000-pluto-woboundary-standardres-binary"
base = "output_restart_3500/"
#base = "output_restart_3500_hires/"
snap_base="snap_"

# Normalization of the accretion rate
Mdot0 = 8.63857932377e-06
# Normalization of the outer density
Sigma_out = 1.0238548710e-4

#Unchanged parameters
alpha = 0.1
h = 0.1
Rmax = 70

def find_simulation(qb,eb,h0,alpha,eta):
    run_name= run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h0,alpha,eta))
    print "Reading from simulation run:", run_name
    directory = run_path+run_name
    snap_path = directory + base + snap_base
    return snap_path
