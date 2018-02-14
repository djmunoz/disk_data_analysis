import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
import os

#rarely changed parameters
outer_rad = 70.0
outer_density = 1.0238548701E-4 # *fixed* density at outer radius
l = 1.0  #  power-law index of temperature profile
outer_accretion_coeff=1.07
Mdot0 = 1.0


def read_accretion_rate_file(filename):
    
    if (not "binary-accretion-rate" in filename):
        print "ERROR: data file of the wrong type"
        exit()
    if ("_q" in filename): qb= float(split(split(filename,"_q")[1],"_")[0])
    if ("_e" in filename): eb= float(split(split(filename,"_e")[1],"_")[0])
    if ("_h" in filename):h= float(split(split(filename,"_h")[1],"_")[0])
    if ("_alpha" in filename):alpha= float(split(split(filename,"_alpha")[1],"_")[0])
    if ("_eta" in filename):eta= float(split(split(filename,"_eta")[1],".txt")[0])
    
    Mdot0 = -3 * np.pi * outer_density * alpha * h**2 * outer_rad**(1.5-l)
    Mdot0 *= outer_accretion_coeff

    times =  np.loadtxt(filename)[:,0]/(2*np.pi)
    rates = np.loadtxt(filename)[:,1:]/np.abs(Mdot0)

    return times, rates[:,0], rates[:,1]

if __name__ == '__main__':


    file_a = 'binary-accretion-rate_norbits3499-3550_q1.0_e0.0_h0.1_alpha0.1_eta1.00.txt'
    file_b = 'binary-accretion-rate_norbits3499-3550_q1.0_e0.0_h0.1_alpha0.1_eta0.10.txt'
  
    
    time_a, mdot1_a, mdot2_a = read_accretion_rate_file(file_a)
    time_b, mdot1_b, mdot2_b = read_accretion_rate_file(file_b)
