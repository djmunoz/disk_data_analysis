# Script to plot the gravitational force acting on twp
# central stars as a function of  time in 2-D disk simulations

# Diego J. Munoz
# 2017


##############################
import numpy as np
import readsnapHDF5 as rs
import glob
from smooth import smooth
from string import split
import sys
###################################################################

outer_rad = 70.0
outer_density = 1.0238548701E-4
l=1.0
h= 0.1
alpha= 0.1
Mdot0 = -3 * np.pi * outer_density * alpha * h**2 * outer_rad**(1.5-l)
outer_accretion_coeff = 1.07
Mdot0 *= outer_accretion_coeff


def KeplerEquation(mean_anom,e):
  
  while (np.abs(mean_anom) > 2.0*np.pi):
    mean_anom-=2.0*np.pi*np.sign(mean_anom)
  if (mean_anom < 0.0): mean_anom = 2*np.pi + mean_anom

  k = 0.85
  ecc_anom = mean_anom + np.sign(np.sin(mean_anom))* k * e
  #ecc_anom = mean_anom
  if (e > 0.8):
    ecc_anom = np.pi
    
  abstol,reltol = 1.0e-8, 1.0e-8
  iter = 0
  while(True):
    f = ecc_anom -e * np.sin(ecc_anom) - mean_anom
    fprime = 1.0 - e * np.cos(ecc_anom)
    fprime2 = e * np.sin(ecc_anom)
    fprime3 = e * np.cos(ecc_anom)
    delta1 = - f / fprime
    delta2 = - f /(fprime + 0.5 * delta1 * fprime2)
    delta3 = - f /(fprime + 0.5 * delta2 * fprime2 + 0.16666666666 * delta2**2 * fprime3) 

    if (delta3 == 0): break
    
    if (np.abs(ecc_anom) > 0.0):
      abserr,relerr = np.abs(delta3),np.abs(delta3)/np.abs(ecc_anom)
    else:
      abserr,relerr = np.abs(delta3),1.0e40

    ecc_anom+=delta3
    #print iter,ecc_anom,e,delta3
    
    if (np.abs(ecc_anom) > abstol/reltol):
      if (abserr < abstol): break
    else:
      if (relerr < reltol): break
    iter+=1      

  return ecc_anom % (2*np.pi)


if __name__ == '__main__':
######################################################################
    #binary and disk properties
    qb=float(sys.argv[1])
    eb=float(sys.argv[2])
    h=0.1
    alpha=0.1

    if (len(sys.argv) < 4): init_snap = 0
    else:   init_snap = int(sys.argv[3]) 
    if (len(sys.argv) < 5): final_snap = 1001
    else:  final_snap = int(sys.argv[4]) 
    if (len(sys.argv) < 6): snap_step = 1
    else:  snap_step = int(sys.argv[5]) 
    
    snap_list = np.arange(init_snap,final_snap+1,snap_step)

    # paths to files
    run_path = "/home/dmunoz/Documents/RESEARCH_PROJECTS/CIRCUMBINARY_DISKS_2D/RESTART_PLUTO_RUNS/"
    
    run_base ="restart-3000-pluto-woboundary-standardres-binary"

    run_name= run_base+("_q%.1f_e%.1f_h%.1f_alpha%.1f/" % (qb,eb,h,alpha))

    print "Reading from simulation run:", run_name

    # OTHER PARAMETERS
    if ("norbits" in run_name): norbits= float(split(split(run_name,"_norbits")[1],"_")[0])
    if ("_q" in run_name): qb= float(split(split(run_name,"_q")[1],"_")[0])
    if ("_e" in run_name): eb= float(split(split(run_name,"_e")[1],"_")[0])
    if ("_h" in run_name):h= float(split(split(run_name,"_h")[1],"_")[0])
    if ("_alpha" in run_name):alpha= float(split(split(run_name,"_alpha")[1],"/")[0])

    directory = run_path+run_name
    base = "output_restart_3500/"

    snap_base="snap_"


    #check the number of orbits at the zeroth and last snapshots
    orbit_range = []
    time = []
    for snap in range(init_snap,final_snap+1):
        filename=directory+base+snap_base+str(snap).zfill(3)
        header = rs.snapshot_header(filename)
        time.append(header.time)
        orbit_range.append(int(header.time/(2*np.pi)))
    
    norbits = orbit_range[0]

    outfilename="gravity-force_norbits%i-%i_q%3.1f_e%3.1f_h%3.1f_alpha%3.1f.txt"\
	% (orbit_range[0],orbit_range[-1],qb,eb,h,alpha)
    print "Saving gravitational forces to file:"
    print "                                  ",outfilename

    
    ecc_anom = np.vectorize(KeplerEquation)(np.array(time)  + np.pi,eb)
    x = np.cos(ecc_anom) - eb
    y = np.sqrt(1 - eb * eb) * np.sin(ecc_anom)
    vx = -np.sin(ecc_anom) / (1 - eb * np.cos(ecc_anom))
    vy = np.sqrt(1 - eb * eb) * np.cos(ecc_anom)/ (1 - eb * np.cos(ecc_anom))
    x2, x1 = -qb/(1 + qb) * x, 1.0/(1 + qb) * x
    y2, y1 = -qb/(1 + qb) * y, 1.0/(1 + qb) * y

    print "Looping over snapshots..."
    force1x_disk,force1y_disk, force2x_disk, force2y_disk = np.zeros(len(time)),np.zeros(len(time)),np.zeros(len(time)),np.zeros(len(time))
    force1x_cav,force1y_cav, force2x_cav, force2y_cav = np.zeros(len(time)),np.zeros(len(time)),np.zeros(len(time)),np.zeros(len(time))
    for i,snap in enumerate(range(init_snap,final_snap+1)):
        print "SNAPSHOT #",snap
        filename=directory+base+snap_base+str(snap).zfill(3)
        header = rs.snapshot_header(filename)
        pos = rs.read_block(filename,"POS ", parttype=0)
        mass = rs.read_block(filename,"MASS", parttype=0)
        acc = rs.read_block(filename,"ACCE", parttype=0)
        ids = rs.read_block(filename,"ID  ", parttype=0)

        r1 = np.sqrt((pos[:,0] - (x1[i] + 0.5  * header.boxsize))**2 +\
                     (pos[:,1] - (y1[i] + 0.5  * header.boxsize))**2)
        r2 = np.sqrt((pos[:,0] - (x2[i] + 0.5  * header.boxsize))**2 +\
                     (pos[:,1] - (y2[i] + 0.5  * header.boxsize))**2)
        r =  np.sqrt((pos[:,0] - 0.5  * header.boxsize)**2 +\
                     (pos[:,1] - 0.5  * header.boxsize)**2)

        ind = (ids >= -2) & (r > (1 + eb))
        force1x_disk[i] = (-mass[ind] / r1[ind] / r1[ind] / r1[ind] * ((x1[i] + 0.5  * header.boxsize) - pos[ind,0])).sum()
        force1y_disk[i] = (-mass[ind] / r1[ind] / r1[ind] / r1[ind] * ((y1[i] + 0.5  * header.boxsize) - pos[ind,1])).sum()
        force2x_disk[i] = (-mass[ind] / r2[ind] / r2[ind] / r2[ind] * ((x2[i] + 0.5  * header.boxsize) - pos[ind,0])).sum()
        force2y_disk[i] = (-mass[ind] / r2[ind] / r2[ind] / r2[ind] * ((y2[i] + 0.5  * header.boxsize) - pos[ind,1])).sum()

        ind = (ids >= -2) & (r <= (1 + eb))
        force1x_cav[i] = (-mass[ind] / r1[ind] / r1[ind] / r1[ind] * ((x1[i] + 0.5  * header.boxsize) - pos[ind,0])).sum()
        force1y_cav[i] = (-mass[ind] / r1[ind] / r1[ind] / r1[ind] * ((y1[i] + 0.5  * header.boxsize) - pos[ind,1])).sum()
        force2x_cav[i] = (-mass[ind] / r2[ind] / r2[ind] / r2[ind] * ((x2[i] + 0.5  * header.boxsize) - pos[ind,0])).sum()
        force2y_cav[i] = (-mass[ind] / r2[ind] / r2[ind] / r2[ind] * ((y2[i] + 0.5  * header.boxsize) - pos[ind,1])).sum()
        
    force1x = force1x_cav + force1x_disk
    force1y = force1y_cav + force1y_disk
    force2x = force2x_cav + force2x_disk
    force2y = force2y_cav + force2y_disk
    
    torque_disk = (x1 - x2) * (force1y_disk - force2y_disk) - (y1 - y2) * (force1x_disk - force2x_disk)
    torque_cav = (x1 - x2) * (force1y_cav - force2y_cav) - (y1 - y2) * (force1x_cav - force2x_cav)
    torque = (x1 - x2) * (force1y - force2y) - (y1 - y2) * (force1x - force2x)
    import matplotlib.pyplot as plt
    plt.plot(time,torque)
    plt.show()
    print torque.mean()/np.abs(Mdot0)
    print torque_disk.mean()/np.abs(Mdot0)
    print torque_cav.mean()/np.abs(Mdot0)

    np.savetxt(outfilename,np.array([time,x1,y1,x2,y2,force1x,force1y,force2x,force2y]).T,fmt='%12f %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g')

    print "Saved data to file:",outfilename

    

    
