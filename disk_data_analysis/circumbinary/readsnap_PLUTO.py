# Python routines to read snapshot output files produced by
# the PLUTO code (Mignone 2007) in ACII, HDF5 or binary formats
#
# by Diego J. Munoz
# 2015
#
#

import numpy as np
import os
import sys
from string import split

#define a header-like class for the PLUTO files
class pluto_header:
	def __init__(self, *args, **kwargs):
		if (len(args) == 1):
                    filename = args[0]
                    f = open(filename,'r')
                    self.filename=filename

                    self.x1in, self.x1out,self.x1points = 0,1,1
                    self.x2in, self.x2out,self.x2points = 0,1,1
                    self.x3in, self.x3out,self.x3points = 0,1,1


                    for line in f:
                        if (line[0] != '#'): continue
                        if ('time' in line):
                            time=0
                            self.time = time
                        if ('DIMENSIONS' in line):
                            self.ndim = int(split(line,"DIMENSIONS:")[1])
                        if ('X1' in line):
                            self.x1in = float(split(split(line,"[")[1],",")[0])
                            self.x1out = float(split(split(line,"]")[0],",")[1])
                            self.x1points = int(split(split(line,"point")[0],",")[-1])
                        if ('X2' in line):
                            self.x2in = float(split(split(line,"[")[1],",")[0])
                            self.x2out = float(split(split(line,"]")[0],",")[1])
                            self.x2points = int(split(split(line,"point")[0],",")[-1])
                        if ('X3' in line):
                            self.x3in = float(split(split(line,"[")[1],",")[0])
                            self.x3out = float(split(split(line,"]")[0],",")[1])
                            self.x3points = int(split(split(line,"point")[0],",")[-1])

                    f.close()
        

        def get_pluto_grid(self):
            f = open(self.filename,'r')
            dim,k = 0,0
            self.x1zones=[]
            self.x2zones=[]
            self.x3zones=[]
            for line in f:
                if (line[0] == '#'): continue
                if (str(self.x1points) == line.rstrip()):
                    dim=1
                    continue
                if (str(self.x2points) == line.rstrip()):
                    dim=2
                    continue
                if (str(self.x3points) == line.rstrip()):
                    if (self.ndim < 3):break
                    dim=3
                    continue

                if (dim == 1):
                    self.x1zones.append([float(split(line)[1]),float(split(line)[2])])
                if (dim == 2):
                    self.x2zones.append([float(split(line)[1]),float(split(line)[2])])    
                
                k+=1

            self.x1zones=np.array(self.x1zones)
            self.x2zones=np.array(self.x2zones)
            self.x3zones=np.array(self.x3zones)
            f.close()
                        
def pluto_read_block(path,snapshot_number,block,header_file,coordinates='cartesian'):
        
        head=pluto_header(path+header_file)
        head.get_pluto_grid()
        
        rcenters=0.5*(head.x1zones[:,0]+head.x1zones[:,1])
        phicenters=0.5*(head.x2zones[:,0]+head.x2zones[:,1])
        Deltar = (head.x1zones[:,1]-head.x1zones[:,0])
        Deltaphi = (head.x2zones[:,1]-head.x2zones[:,0])
        
        r,phi = np.meshgrid(rcenters,phicenters)
        dr,dphi = np.meshgrid(Deltar,Deltaphi)
        r=r.flatten()
        phi=phi.flatten()
        dr=dr.flatten()
        dphi=dphi.flatten()

        if (block == "R   "):
                ret_val = r
                
        if (block == "PHI "):
                ret_val = phi
                
        if (block == "POS "):
                if (coordinates=='cartesian'):
                        x=r*np.cos(phi)
                        y=r*np.sin(phi)
                        z=np.zeros(x.shape[0])
                        ret_val=np.array([x,y,z]).T
                elif (coordinates=='polar'):
                        z=np.zeros(r.shape[0])
                        ret_val = np.array([r,phi,z]).T

        if (block == "VEL ") | (block == "VR  ") | (block == "VPHI"):
                if (block == "VR  "):
                        # Velocity in radial direction
                        filename=path+"vx1."+str(snapshot_number).zfill(4)+".dbl.-001.ascii"
                        velr=np.loadtxt(filename)
                        ret_val = velr
                if (block == "VPHI"):
                        # Velocity in azimuthal direction
                        filename=path+"vx2."+str(snapshot_number).zfill(4)+".dbl.-001.ascii"
                        velphi=np.loadtxt(filename)
                        ret_val = velphi        
                if (block == "VEL "):
                        filename=path+"vx1."+str(snapshot_number).zfill(4)+".dbl.-001.ascii"
                        velr=np.loadtxt(filename)
                        filename=path+"vx2."+str(snapshot_number).zfill(4)+".dbl.-001.ascii"
                        velphi=np.loadtxt(filename)
                        
                        velx=np.cos(phi)*velr - np.sin(phi) * velphi
                        vely=np.sin(phi)*velr + np.cos(phi) * velphi
                        velz=np.zeros(velx.shape[0])
                        ret_val=np.array([velx,vely,velz]).T


        if (block == "RHO "):
                filename=path+"rho."+str(snapshot_number).zfill(4)+".dbl.-001.ascii"
                ret_val=np.loadtxt(filename)

        if (block == "VOL "):
                ret_val = np.array(r * dr * dphi)

        if (block == "MASS"):
                filename=path+"rho."+str(snapshot_number).zfill(4)+".dbl.-001.ascii"
                ret_val= np.loadtxt(filename) * np.array(r * dr * dphi)
                
        if (block == "ID  "):
                ret_val = np.arange(1,r.shape[0]+1).astype(int)

        
        return ret_val
