# Script to plot the accretion rate of central stars as a function of
# time in 2-D disk simulations

# Diego J. Munoz
# 2016


##############################
import numpy as np
import matplotlib.pyplot as plt
import glob
from smooth import smooth
from string import split
import sys
import disk_data_analysis.circumbinary as dda
from disk_data_analysis.orbital import orbit_in_time

import simulation_settings as simset

plt.rcParams['mathtext.fontset'] = "stix"



time_offset=0

MAX_NLINES=20000 #cap the size of the accretion data file
SMOOTH=0




if __name__ == '__main__':
######################################################################
    #binary and disk properties

     
    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h0=simset.h
    alpha=simset.alpha
    if (len(sys.argv) > 3):
        eta = float(sys.argv[3])
    else:
        eta = 1.0

    if (len(sys.argv) < 5): orbit_init = None
    else:   orbit_init = int(sys.argv[4])
    if (len(sys.argv) < 6): orbit_final = None
    else:  orbit_final = int(sys.argv[5])

    
    run_name= simset.run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f_eta%.2f/" % (qb,eb,h0,alpha,eta))
    print "Reading from simulation run:", run_name
    directory = simset.run_path+run_name
    snap_path = directory + simset.base + simset.snap_base
    print "Reading from simulation run:", run_name


    if ("sink" in simset.base): sink_label='-sink' 
    else: sink_label=''
    if ("hires" in simset.base): 
        if ("_hires" in simset.base):res_label='-hires' 
        elif ("_vhires" in simset.base): res_label='-vhires' 
        elif ("_vvhires" in simset.base): res_label='-vvhires' 
    else: res_label=''


    labels=sink_label+res_label

    #check the number of orbits at the zeroth and last snapshots
    time_list = []
    snap_list = []
    for snapfile in glob.iglob(snap_path+'*.hdf5'):
        snapnum = int(snapfile.split(snap_path)[1].split('.hdf5')[0])
        snap_list.append(snapnum)

    for snapnum in [min(snap_list),max(snap_list)]:
        snap = dda.get_snapshot_data(snap_path,snapnum,['POS'],parttype=0)
        time = snap.header.time
        time_list.append(time)
    if (orbit_init is None):
        orbit_init = int(min(time_list)/2/np.pi)
    if (orbit_final is None):
        orbit_final = int(max(time_list)/2/np.pi)
    time_min = orbit_init * 2 * np.pi
    time_max = orbit_final * 2 * np.pi




    # Now read accretion data file
    print "Reading in data..."
    #accretionfile=directory+base+'circumstellarsink.txt'
    accretionfile=directory+simset.base+'circumstellarsink_justmasses.txt'
    accretionfile=directory+simset.base+'circumstellarsink.txt'
    accretion_data=np.loadtxt(accretionfile)
    accretion_time=accretion_data[:,0]
    accretion_m1=accretion_data[:,1]
    accretion_m2=accretion_data[:,2]
    accretion_v1x=accretion_data[:,3]
    accretion_v1y=accretion_data[:,4]
    accretion_v2x=accretion_data[:,5]
    accretion_v2y=accretion_data[:,6]
    accretion_fext1x=accretion_data[:,7]
    accretion_fext1y=accretion_data[:,8]
    accretion_fext2x=accretion_data[:,9]
    accretion_fext2y=accretion_data[:,10]
    if (accretion_data.shape[1] == 15):
        accretion_v1x_g=accretion_data[:,11]
        accretion_v1y_g=accretion_data[:,12]
        accretion_v2x_g=accretion_data[:,13]
        accretion_v2y_g=accretion_data[:,14]

    if (accretion_time.shape[0] > 1.5*MAX_NLINES):
	    skip = int(np.ceil(accretion_time.shape[0]/(1.5*MAX_NLINES)))
	    accretion_time = accretion_time[::skip]
	    accretion_m1 = accretion_m1[::skip]
	    accretion_m2 = accretion_m2[::skip]
            accretion_v1x=accretion_v1x[::skip]
            accretion_v1y=accretion_v1y[::skip]
            accretion_v2x=accretion_v2x[::skip]
            accretion_v2y=accretion_v2y[::skip]
            accretion_fext1x=accretion_fext1x[::skip]
            accretion_fext1y=accretion_fext1y[::skip]
            accretion_fext2x=accretion_fext2x[::skip]
            accretion_fext2y=accretion_fext2y[::skip]
            if (accretion_data.shape[1] == 15):
                accretion_v1x_g = accretion_v1x_g[::skip]
                accretion_v1y_g = accretion_v1y_g[::skip]
                accretion_v2x_g = accretion_v2x_g[::skip]
                accretion_v2y_g = accretion_v2y_g[::skip]
    ############################################
    if (snap_list is None):
        orbit_init = int(accretion_time[0]/2/np.pi)
        orbit_final = int(accretion_time[-1]/2/np.pi)
        time_min = accretion_time[0]
        time_max = accretion_time[-1]
        
    outfilename1="binary-forcing-rate"+labels+"_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f.txt"\
        % (orbit_init,orbit_final,qb,eb,h0,alpha)
    outfilename2="binary-accretion-rate"+labels+"_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f.txt"\
        % (orbit_init,orbit_final,qb,eb,h0,alpha)

    print "Saving accretion force data to file:"
    print "                                  ",outfilename1
    
    print "Saving accretion data to file:"
    print "                                  ",outfilename2



    #######################################################
    m1 = accretion_m1[-1]
    m2 = accretion_m2[-1]
    m1_final=accretion_m1[-1]
    m2_final=accretion_m2[-1]

    #check that accretion times do not have discontinuities
    while (np.any(np.diff(accretion_time) <= 0)):
        print "Warning"
        #find where there was a discontinuity (due to a restart)
        '''
        t0= accretion_time[1:][np.diff(accretion_time) <= 0][0]
        ind0 = np.where(np.diff(accretion_time) <= 0)[0][0]-1
        #find the masses at that point
        m1_0 = accretion_m1[ind0]
        m2_0 = accretion_m2[ind0]
        #skip the discontinuity and reinstate the mass right before it
        ind1 = np.where(accretion_time > t0)[0][0]
        accretion_time = np.append(accretion_time[:ind0],accretion_time[ind1:])
        accretion_m1= np.append(accretion_m1[:ind0],accretion_m1[ind1:]+m1_0)
        accretion_m2= np.append(accretion_m2[:ind0],accretion_m2[ind1:]+m2_0)
        '''
        t0= accretion_time[1:][np.diff(accretion_time) <= 0][0]
        ind0 = np.argwhere(np.diff(accretion_time) <= 0)[0][0]
        # Find the masses and velocities at the break
        m1_0 = accretion_m1[ind0]
        m2_0 = accretion_m2[ind0]
        v1x_0,v1y_0 = accretion_v1x[ind0],accretion_v1y[ind0]
        v2x_0,v2y_0 = accretion_v2x[ind0],accretion_v2y[ind0]
        # skip the discontinuity but keep track of the mass
        accretion_time = np.append(accretion_time[:ind0],accretion_time[ind0+1:])
        accretion_m1= np.append(accretion_m1[:ind0],accretion_m1[ind0+1:]+m1_0)
        accretion_m2= np.append(accretion_m2[:ind0],accretion_m2[ind0+1:]+m2_0)
        accretion_v1x= np.append(accretion_v1x[:ind0],accretion_v1x[ind0+1:]+v1x_0)
        accretion_v2x= np.append(accretion_v2x[:ind0],accretion_v2x[ind0+1:]+v2x_0)
        accretion_v1y= np.append(accretion_v1y[:ind0],accretion_v1y[ind0+1:]+v1y_0)
        accretion_v2y= np.append(accretion_v2y[:ind0],accretion_v2y[ind0+1:]+v2y_0)
        accretion_fext1x = np.append(accretion_fext1x[:ind0],accretion_fext1x[ind0+1:])
        accretion_fext1y = np.append(accretion_fext1y[:ind0],accretion_fext1y[ind0+1:])
        accretion_fext2x = np.append(accretion_fext2x[:ind0],accretion_fext2x[ind0+1:])
        accretion_fext2y = np.append(accretion_fext2y[:ind0],accretion_fext2y[ind0+1:])
        if (accretion_data.shape[1] == 15):
            v1x_g_0,v1y_g_0 = accretion_v1x_g[ind0],accretion_v1y_g[ind0]
            v2x_g_0,v2y_g_0 = accretion_v2x_g[ind0],accretion_v2y_g[ind0]
            accretion_v1x_g= np.append(accretion_v1x_g[:ind0],accretion_v1x_g[ind0+1:]+v1x_g_0)
            accretion_v2x_g= np.append(accretion_v2x_g[:ind0],accretion_v2x_g[ind0+1:]+v2x_g_0)
            accretion_v1y_g= np.append(accretion_v1y_g[:ind0],accretion_v1y_g[ind0+1:]+v1y_g_0)
            accretion_v2y_g= np.append(accretion_v2y_g[:ind0],accretion_v2y_g[ind0+1:]+v2y_g_0)
        if (SMOOTH):
            #Smooth data over to filter out fluctuations over 
	    #irrelevantly short time-scales
	    print "Smoothing time-series..."
	    softlen=0.01
	    wlen1= (2*np.pi/np.sqrt(m1/softlen**3)/np.diff(accretion_time).mean())
	    accretion_m1=smooth(accretion_m1,window_len=wlen1,window="hamming")
	    wlen2= (2*np.pi/np.sqrt(m2/softlen)/np.diff(accretion_time).mean())
	    accretion_m2=smooth(accretion_m2,window_len=wlen2,window="hamming")
    
    #for kk in range(accretion_time.shape[0]):
    #    print accretion_time[kk],accretion_m1[kk],accretion_m2[kk]


    #############
        
    #Take the time derivative to obtain the accretion rate
    dt=np.gradient(accretion_time)
    mdot_m1 = np.gradient(accretion_m1)/dt
    mdot_m2 = np.gradient(accretion_m2)/dt


    dv1xdt_acc = np.gradient(accretion_v1x)/dt
    dv1ydt_acc = np.gradient(accretion_v1y)/dt
    dv2xdt_acc = np.gradient(accretion_v2x)/dt
    dv2ydt_acc = np.gradient(accretion_v2y)/dt

    if (accretion_data.shape[1] == 15):
       dv1xdt_grav = np.gradient(accretion_v1x_g)/dt
       dv1ydt_grav = np.gradient(accretion_v1y_g)/dt
       dv2xdt_grav = np.gradient(accretion_v2x_g)/dt
       dv2ydt_grav = np.gradient(accretion_v2y_g)/dt 


    ind = (accretion_time >= time_min) & (accretion_time <= time_max) & (mdot_m1 > 0) & (mdot_m2 > 0)
    accretion_time = accretion_time[ind]
    mdot_m1,mdot_m2 = mdot_m1[ind], mdot_m2[ind]
    dv1xdt_acc, dv1ydt_acc = dv1xdt_acc[ind],dv1ydt_acc[ind]
    dv2xdt_acc, dv2ydt_acc = dv2xdt_acc[ind],dv2ydt_acc[ind]
    accretion_fext1x, accretion_fext1y = accretion_fext1x[ind],accretion_fext1y[ind]
    accretion_fext2x, accretion_fext2y = accretion_fext2x[ind],accretion_fext2y[ind] 
    if (accretion_data.shape[1] == 15):
        dv1xdt_grav, dv1ydt_grav = dv1xdt_grav[ind],dv1ydt_grav[ind]
        dv2xdt_grav, dv2ydt_grav = dv2xdt_grav[ind],dv2ydt_grav[ind]

    mu =  qb / ( 1.0 + qb)
    xb, yb, vxb, vyb = np.vectorize(orbit_in_time)(accretion_time + np.pi, eb)
    x1, y1 = mu * xb, mu * yb
    x2, y2 = -(1.0 - mu) * xb, -(1.0 - mu) * yb
    vx1, vy1 = mu * vxb, mu * vyb
    vx2, vy2 = -(1.0 - mu) * vxb, -(1.0 - mu) * vyb
    '''
    x,y,vx,vy = [],[],[],[]
    for t in accretion_time:
        ecc_anom = orbital.KeplerEquation(t + np.pi,eb)
        x.append(np.cos(ecc_anom) - eb)
        y.append(np.sqrt(1 - eb * eb) * np.sin(ecc_anom))
        vx.append(-np.sin(ecc_anom) / (1 - eb * np.cos(ecc_anom)))
        vy.append(np.sqrt(1 - eb * eb) * np.cos(ecc_anom)/ (1 - eb * np.cos(ecc_anom)))
    x,y,vx,vy = np.array(x),np.array(y),np.array(vx),np.array(vy)
    mu =  qb / ( 1.0 + qb)
    x1, y1 = mu * x, mu * y
    x2, y2 = -(1.0 - mu) * x, -(1.0 - mu) * y
    vx1, vy1 = mu * vx, mu * vy
    vx2, vy2 = -(1.0 - mu) * vx, -(1.0 - mu) * vy
    '''


    if (accretion_data.shape[1] == 15):
        np.savetxt(outfilename1,np.array([accretion_time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,
                                          dv1xdt_acc,dv1ydt_acc,dv2xdt_acc,dv2ydt_acc,
                                          dv1xdt_grav,dv1ydt_grav,dv2xdt_grav,dv2ydt_grav,
                                      ]).T,
                   fmt='%12f %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g')
    else:
        np.savetxt(outfilename1,np.array([accretion_time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,
                                     dv1xdt_acc,dv1ydt_acc,dv2xdt_acc,dv2ydt_acc,
                                     accretion_fext1x,accretion_fext1y,accretion_fext2x,accretion_fext2y]).T,
                   fmt='%12f %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g')
  

    print "Saved accretion rate data to file:",outfilename1

    np.savetxt(outfilename2,np.array([accretion_time,mdot_m1,mdot_m2]).T,
               fmt='%12f %.8g %.8g')


    print "Saved accretion rate data to file:",outfilename2


    print "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"

    
