from disk_hdf5 import readsnapHDF5 as rs
from disk_voronoi import voronoi_simulation_data 
#import readsnap_PLUTO as rspluto
#import pluto_data_utils as rspluto
import numpy as np
import sys


disk_data_fields = ["POS","VEL","MASS","U","RHO","R","PHI","VELR","VELPHI"]


datablocks = {"POS ":["Coordinates",3], 
	      "VEL ":["Velocities",3],
	      "ID  ":["ParticleIDs",1],
	      "MASS":["Masses",1],
              "U   ":["InternalEnergy",1],
              "RHO ":["Density",1],
              "VOL ":["Volume",1],
	      "PRES":["Pressure",1],
              "CMCE":["CenterOfMass",3],
              "AREA":["SurfaceArea",1],
              "ACCE":["Acceleration",3],
              "GRAR":["DensityGradient", 3],
              "GRAV":["VelocityGradient", 9]
              }


STAR_PARTTYPE = 4


class gas_data():
    def __init__(self,*args,**kwargs):
        #self.pos=kwargs.get("POS")
        #self.vel=kwargs.get("VEL")
        #self.dens=kwargs.get("RHO")
        #self.utherm=kwargs.get("U")
        #self.ids=kwargs.get("ID")

        #self.R = kwargs.get("R")
        #self.phi = kwargs.get("PHI")
        #self.velR = kwargs.get("VELR")
        #self.velphi = kwargs.get("VELPHI")

        for key in kwargs:
            setattr(self, key, kwargs[key])
        
            
class particle_data():
    def __init__(self,*args,**kwargs):
        self.pos=kwargs.get("pos")
        self.vel=kwargs.get("vel")
        self.mass=kwargs.get("mass")
        self.ids=kwargs.get("ids")

        if (self.pos is None):
            self.pos = np.empty([0,3])
        if (self.vel is None):
            self.vel = np.empty([0,3])
        if (self.mass is None):
            self.mass = np.empty([0])
        if (self.ids is None):
            self.ids = np.empty([0])


class snapshot():
    """
    Snapshot class, containing -- in a single data structure -- the simulation snapshot information

    """
    
    def __init__(self,*args,**kwargs):

        self.header = kwargs.get("header")
        
        self.parttype=kwargs.get("parttype")
        
        if (self.parttype == 0):
            self.gas = gas_data(**kwargs)
        if (self.parttype == 1):
            self.particle = particle_data(**kwargs)
        

def get_snapshot_data(filename_prefix,snap_num,quantities,parttype= 0 ,code="AREPO"):

    nquant=len(quantities)
    outquant = []
    if (parttype is None): types = [0,1,2,3,4]
    else: types = [parttype]

    for quant in quantities:
        if (code == "AREPO"):
            header = rs.snapshot_header(filename_prefix+str(snap_num).zfill(3))
            for parttype in types:
                if (quant.ljust(4) in datablocks):
                    outquant.append(rs.read_block(filename_prefix+str(snap_num).zfill(3),quant.ljust(4),parttype=parttype))
                elif (quant in disk_data_fields):
                    BoxX,BoxY = header.boxsize,header.boxsize
                    pos = rs.read_block(filename_prefix+str(snap_num).zfill(3),"POS ", parttype=parttype)
                    pos[:,0],pos[:,1] = (pos[:,0]- BoxX/2), (pos[:,1]- BoxY/2)
                    radius = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
                    phi = np.arctan2(pos[:,1],pos[:,0])
                    if (quant == "R"):
                        outquant.append(radius)
                    if (quant == "PHI"):
                        outquant.append(phi)
                    if ((quant == "VELR") | (quant == "VELPHI")):
                        vel = rs.read_block(filename_prefix+str(snap_num).zfill(3),"VEL ", parttype=parttype)
                        vphi = -np.sin(phi) * vel[:,0] + np.cos(phi) * vel[:,1]
                        vr   =  np.cos(phi) * vel[:,0] + np.sin(phi) * vel[:,1]
                        if (quant == "VELR"):
                            outquant.append(vr)
                        if (quant == "VELPHI"):
                            outquant.append(vphi)

                else: 
                    print "[error] Quantity type ", quant, "not known!"
                    sys.exit()  

        elif (code == "PLUTO"):
            if (quant.ljust(4) in datablocks):
                outquant.append(rspluto.read_data_block(filename_prefix,snap_num,quant.ljust(4)))
            elif (quant in disk_data_fields):
                if ((quant == "R") | (quant == "PHI")):
                    pos = rspluto.read_data_block(filename_prefix,snap_num,"POS ")
                    radius = pos[:,0]
                    phi = pos[:,1]
                    if (quant == "R"):
                        outquant.append(radius)
                    if (quant == "PHI"):
                        outquant.append(phi)
                if ((quant == "VELR") | (quant == "VELPHI")):
                    vel = rspluto.read_data_block(filename_prefix,snap_num,"VEL ")
                    vphi = vel[:,1]
                    vr   = vel[:,0]
                    if (quant == "VELR"):
                        outquant.append(vr)
                    if (quant == "VELPHI"):
                        outquant.append(vphi)
            else: 
                print "[error] Quantity type ", quant, "not known!"
                sys.exit()  
                
        else:
            print "[error] Simulations code ", code, "not known or supported!"
            sys.exit()  


    attributes = dict(zip(quantities,outquant))
    snap = snapshot(parttype=parttype,header=header,**attributes)
    
    return snap


def get_snapshot_time(filename,snapnum,code="AREPO",snapinterval=0):
    if (code == "AREPO"):
        header = rs.snapshot_header(filename)
        time = rs.snapshot_header(filename).time
    elif (code == "PLUTO"):
        time = snapnum*snapinterval

    return time


def compute_snapshot_gradient(snapshot,field):

    # Read-in mesh-generating points
    points = np.array([snapshot.gas.POS[:,0],snapshot.gas.POS[:,1]]).T
    # Create Voronoi mesh
    mesh = voronoi_simulation_data.VoronoiMesh(points)
    # Compute the Voronoi-Green-Gauss gradients
    print getattr(snapshot.gas,field)
    gradientx,gradienty = voronoi_simulation_data.compute_voronoi_gradients(mesh,getattr(snapshot.gas,field))
    # Compute the slope limiters
    limiter = voronoi_simulation_data.compute_voronoi_limiter(mesh,getattr(snapshot.gas,field),np.array([gradientx,gradienty]).T)
    
    # Return limited gradients

    return np.array([gradientx*limiter, gradienty * limiter]).T
