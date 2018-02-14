# Script to plot the accretion rate of central stars as a function of
# time in 2-D disk simulations

# Diego J. Munoz
# 2016


##############################
import numpy as np
import readsnapHDF5 as rs
import glob
from smooth import smooth
from string import split
import sys
###################################################################

time_offset=0

MAX_NLINES=20000 #cap the size of the accretion data file
SMOOTH=0


if __name__ == '__main__':
######################################################################
    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h=0.1
    alpha=0.1
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

    # paths to files
    run_path = "/home/dmunoz/Documents/RESEARCH_PROJECTS/CIRCUMBINARY_DISKS_2D/RESTART_PLUTO_RUNS/"
    
    run_base ="restart-3000-pluto-woboundary-standardres-binary"

    run_name= run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h,alpha,eta))

    print "Reading from simulation run:", run_name

    # OTHER PARAMETERS
    if ("norbits" in run_name): norbits= float(split(split(run_name,"_norbits")[1],"_")[0])
    if ("_q" in run_name): qb= float(split(split(run_name,"_q")[1],"_")[0])
    if ("_e" in run_name): eb= float(split(split(run_name,"_e")[1],"_")[0])
    if ("_h" in run_name):h= float(split(split(run_name,"_h")[1],"_")[0])
    if ("_alpha" in run_name):alpha= float(split(split(run_name,"_alpha")[1],"_")[0])
    if ("_eta" in run_name):eta= float(split(split(run_name,"_eta")[1],"/")[0])

    directory = run_path+run_name
    base = "output_restart_3500/"

    snap_base="snap_"

    if ("sink" in base): sink_label='-sink' 
    else: sink_label=''
    if ("hires" in base): 
        if ("_hires" in base):res_label='-hires' 
        elif ("_vhires" in base): res_label='-vhires' 
        elif ("_vvhires" in base): res_label='-vvhires' 
    else: res_label=''

    labels=sink_label+res_label

    #check the number of orbits at the zeroth and last snapshots
    orbit_range = []
    for snap in [init_snap,final_snap]:
        filename=directory+base+snap_base+str(snap).zfill(3)
        header = rs.snapshot_header(filename)
        time = header.time + time_offset
        orbit_range.append(int(time/(2*np.pi)))
    
    norbits = orbit_range[0]

    outfilename="binary-accretion-rate"+labels+"_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f_eta%.2f.txt"\
	% (orbit_range[0],orbit_range[1],qb,eb,h,alpha,eta)
    print "Saving accretion rate data to file:"
    print "                                  ",outfilename

    
    # Now read accretion data file
    print "Reading in data..."
    accretionfile=directory+base+'circumstellarsink.txt'
    #accretionfile=directory+base+'circumstellarsink_justmasses.txt'
    accretion_data=np.loadtxt(accretionfile)
    accretion_time=accretion_data[:,0]
    accretion_m1=accretion_data[:,1]
    accretion_m2=accretion_data[:,2]
    
    print accretion_time.shape,accretion_m1.shape,accretion_m2.shape

    if (accretion_time.shape[0] > 1.5*MAX_NLINES):
	    skip = int(np.ceil(accretion_time.shape[0]/(1.5*MAX_NLINES)))
	    accretion_time = accretion_time[::skip]
	    accretion_m1 = accretion_m1[::skip]
	    accretion_m2 = accretion_m2[::skip]

    
    ############################################
    m1 = accretion_m1[-1]
    m2 = accretion_m2[-1]
    m1_final=accretion_m1[-1]
    m2_final=accretion_m2[-1]

    #check that accretion times do not have discontinuities
    while (np.any(np.diff(accretion_time) <= 0)):
        print "Warning"
        #find where there was a discontinuity (due to a restart)
        t0= accretion_time[np.diff(accretion_time) <= 0][0]
        ind0 = np.where(np.diff(accretion_time) <= 0)[0][0]-1
        #find the masses at that point
        m1_0 = accretion_m1[ind0]
        m2_0 = accretion_m2[ind0]
        #skip the discontinuity and reinstate the mass right before it
        ind1 = np.where(accretion_time > t0)[0][0]
        accretion_time = np.append(accretion_time[:ind0],accretion_time[ind1:])
        accretion_m1= np.append(accretion_m1[:ind0],accretion_m1[ind1:]+m1_0)
        accretion_m2= np.append(accretion_m2[:ind0],accretion_m2[ind1:]+m2_0)
    #similarly
    m1_shift=0
    m2_shift=0
    #while (np.any(np.diff(accretion_m1) < 0)):
    #    print "Warning"
    #    #find where there was a discontinuity (due to a restart)
    #    t0= accretion_time[np.diff(accretion_m1) <= 0][0]
    #    ind0 = np.where(np.diff(accretion_m1) <= 0)[0][0]-1
    #    #find the masses at that point
    #    m1_0 = accretion_m1[ind0]
    #    m2_0 = accretion_m2[ind0]
    #    #skip the discontinuity and reinstate the mass right before it
    #    ind1 = np.where(accretion_time > t0)[0][0]
    #    t1 = accretion_time[ind1]
    #    m1_1 = accretion_m1[ind1]
    #    m2_1 = accretion_m2[ind1]
    #    print t0/2/np.pi,m1_0,m2_0,t1/2/np.pi,m1_1,m2_1
    #    accretion_time = np.append(accretion_time[:ind0-1],accretion_time[ind0+2:])
    #    accretion_m1= np.append(accretion_m1[:ind0-1],accretion_m1[ind0+2:]+(m1_0-m1_shift))
    #    accretion_m2= np.append(accretion_m2[:ind0-1],accretion_m2[ind0+2:]+(m2_0-m2_shift))        
    #    t1 = accretion_time[ind1]
    #    m1_1 = accretion_m1[ind1]
    #    m2_1 = accretion_m2[ind1]
    #    m1_shift=+m1_0
    #    m2_shift=+m2_0
    #    print t0,m1_0,m2_0,t1,m1_1,m2_1
   

    if (SMOOTH):
            #Smooth data over to filter out fluctuations over 
	    #irrelevantly short time-scales
	    print "Smoothing time-series..."
	    softlen=0.01
	    wlen1= (2*np.pi/np.sqrt(m1/softlen**3)/np.diff(accretion_time).mean())
	    accretion_m1=smooth(accretion_m1,window_len=wlen1,window="hamming")
	    wlen2= (2*np.pi/np.sqrt(m2/softlen)/np.diff(accretion_time).mean())
	    accretion_m2=smooth(accretion_m2,window_len=wlen2,window="hamming")
    
    for kk in range(accretion_time.shape[0]):
        print accretion_time[kk],accretion_m1[kk],accretion_m2[kk]


    #############
        
    P_orbit = 2*np.pi
    time_list=[]
    for num in snap_list:
	    filename=directory+base+snap_base+str(num).zfill(3)
	        #open the snapshot header
	    header = rs.snapshot_header(filename)
	    time_list.append(header.time)
	    print "SNAPSHOT #",num," TIME=",header.time
    
    #Take the time derivative to obtain the accretion rate

    dt=np.gradient(accretion_time)
    mdot_m1 = np.gradient(accretion_m1)/dt
    mdot_m2 = np.gradient(accretion_m2)/dt
    
    ind = (accretion_time > time_list[0]) & (accretion_time < time_list[-1]) & (mdot_m1 > 0) & (mdot_m2 > 0)
    ind = (accretion_time > time_list[0])  & (mdot_m1 > 0) & (mdot_m2 > 0)

    np.savetxt(outfilename,np.array([accretion_time[ind],mdot_m1[ind],mdot_m2[ind]]).T,fmt='%12f %.8g %.8g')


    print "Saved accretion rate data to file:",outfilename

    

    
