# Script to plot the accretion rate of central stars as a function of
# time in 2-D disk simulations

# Diego J. Munoz
# 2017

import numpy as np
import sys
import disk_data_analysis.circumbinary as dda



if __name__ == '__main__':
     
    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h0=simset.h
    alpha=simset.alpha
    eta = float(sys.argv[3])
    s = float(sys.argv[4])
    
    if (len(sys.argv) < 6): orbit_init = None
    else:   orbit_init = int(sys.argv[5])
    if (len(sys.argv) < 7): orbit_final = None
    else:  orbit_final = int(sys.argv[6])

    snap_path,snap_base = simset.find_simulation(qb,eb,h0,alpha,eta,s)

    outfilename1="binary-forcing-rate_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f_eta%.2f_s%.2f.txt"\
        % (orbit_init,orbit_final,qb,eb,h0,alpha,eta,s)
    outfilename2="binary-accretion-rate_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f_eta%.2f_s%.2f.txt"\
        % (orbit_init,orbit_final,qb,eb,h0,alpha,eta,s)

    accretionfile = snap_path+'circumstellarsink.txt'
    dda.write_binary_externalforces_file(accretionfile,outfilename1,outfilename2,qb=qb,eb=eb,
                                         orbit_init = orbit_init,orbit_final = orbit_final)





