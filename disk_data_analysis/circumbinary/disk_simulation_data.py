from __future__ import print_function
from disk_hdf5 import snapHDF5 as rs
try:
    from disk_voronoi import voronoi_simulation_data 
except ImportError:
    None
import readsnap_PLUTO as rspluto
#import pluto_data_utils as rspluto
import numpy as np
import sys


disk_data_fields = ["POS","VEL","MASS","U","RHO","R","PHI","VELR","VELPHI","VELX","VELY",
                    "ACCE","GRAP","VISC"]


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
              "GRAV":["VelocityGradient", 9],
              "AREA":["SurfaceArea",1],
              "GRAP":["PressureGradient", 3],
              "VISC":["Viscosity",1],
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
        
    def add_gas_data(self,data,fieldname):
        setattr(self, fieldname, data)
        

            
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

    def add_particle_data(self,data,fieldname):
        setattr(self, fieldname, data)
        
class snapshot():
    """
    Snapshot class, containing -- in a single data structure -- the simulation snapshot information

    """
    
    def __init__(self,*args,**kwargs):

        self.header = kwargs.get("header")
        
        self.parttype=kwargs.get("parttype")
        
        if (self.parttype == 0):
            self.gas = gas_data(**kwargs)
        else:
            self.gas = None
        if (self.parttype > 1) & (self.parttype <= 5):
            self.particle = particle_data(**kwargs)
        else:
            self.particle = None
            
    def add_data(self,data,fieldname,parttype=0):
        if (parttype == 0):
            self.gas.add_gas_data(data,fieldname)
        if (parttype > 1) & (parttype <= 5):
            self.gas.add_particle_data(data,fieldname)
            
def get_snapshot_data(filename_prefix='./',snap_num = None,quantities=None,parttype= None ,
                      code="AREPO",filename = None):

    nquant=len(quantities)
    outquant = []
    if (parttype is None): types = [0,1,2,3,4]
    else: types = [parttype]

    if (filename is None):
        path = filename_prefix+str(snap_num).zfill(3)
    else:
        path  = filename
        
    if (code == "AREPO"):
        header = rs.snapshot_header(path)

    if (code == "PLUTO"):
        header = rs.snapshot_header()

        
    if (nquant == 0):
        snap = snapshot(parttype=parttype,header=header)
        return snap

    
    for quant in quantities:
        if (code == "AREPO"):
            for parttype in types:
                if (quant.ljust(4) in datablocks):
                    outquant.append(rs.read_block(path,quant.ljust(4),parttype=parttype))
                elif (quant in disk_data_fields):
                    BoxX,BoxY = header.boxsize,header.boxsize
                    pos = rs.read_block(path,"POS ", parttype=parttype)
                    pos[:,0],pos[:,1] = (pos[:,0]- BoxX/2), (pos[:,1]- BoxY/2)
                    radius = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
                    phi = np.arctan2(pos[:,1],pos[:,0])
                    if (quant == "R"):
                        outquant.append(radius)
                    if (quant == "PHI"):
                        outquant.append(phi)
                    if ((quant == "VELR") | (quant == "VELPHI") | (quant == "VELX") | (quant == "VELY" )):
                        vel = rs.read_block(path,"VEL ", parttype=parttype)
                        vphi = -np.sin(phi) * vel[:,0] + np.cos(phi) * vel[:,1]
                        vr   =  np.cos(phi) * vel[:,0] + np.sin(phi) * vel[:,1]
                        if (quant == "VELR"):
                            outquant.append(vr)
                        if (quant == "VELPHI"):
                            outquant.append(vphi)
                    if (quant == "VELX"):
                        outquant.append(vel[:,0])
                    if (quant == "VELY"):
                        outquant.append(vel[:,1])
                else: 
                    print("[error] Quantity type ", quant, "not known!")
                    sys.exit()  

        elif (code == "PLUTO"):
            header_file = 'grid.out'
            if (quant.ljust(4) in datablocks):
                outquant.append(rspluto.pluto_read_block(filename_prefix,snap_num,quant.ljust(4),
                                                         header_file=header_file))
            elif (quant in disk_data_fields):
                if (quant == "R"):
                    radius = rspluto.pluto_read_block(filename_prefix,snap_num,"R   ",header_file=header_file)
                    outquant.append(radius)
                if (quant == "PHI"):
                    phi = rspluto.pluto_read_block(filename_prefix,snap_num,"PHI ",header_file=header_file)
                    outquant.append(phi)
                if (quant == "VELR"):
                    vr = rspluto.pluto_read_block(filename_prefix,snap_num,"VR  ",header_file=header_file)
                    outquant.append(vr)
                if (quant == "VELPHI"):
                    vphi =rspluto.pluto_read_block(filename_prefix,snap_num,"VPHI",header_file=header_file)
                    outquant.append(vphi)
            else: 
                print("[error] Quantity type ", quant, "not known!")
                sys.exit()  
        else:
            print("[error] Simulations code ", code, "not known or supported!")
            sys.exit()  


    attributes = dict(zip(quantities,outquant))
    snap = snapshot(parttype=parttype,header=header,**attributes)
    
    return snap


def write_snapshot(snapshot,filename="./disk.dat.hdf5",time=0, \
                   relax_density_in_input = False):

    f=rs.openfile(filename)
    rs.writeheader(f, snapshot.header)
    rs.write_block(f, "POS ", 0, snapshot.gas.POS)
    rs.write_block(f, "VEL ", 0, snapshot.gas.VEL)
    if (relax_density_in_input):
        rs.write_block(f, "MASS", 0, snapshot.gas.RHO)
    else:
        if (snapshot.gas.MASS is not None):
            rs.write_block(f, "MASS", 0, snapshot.gas.MASS)
        rs.write_block(f, "RHO ", 0, snapshot.gas.RHO)
    #rs.write_block(f, "U   ", 0, snapshot.gas.UTHERM)
    rs.write_block(f, "ID  ", 0, snapshot.gas.ID)
    
    if (snapshot.particle is not None):
        rs.write_block(f, "POS ", STAR_PARTTYPE, snapshot.particle.POS)
        rs.write_block(f, "VEL ", STAR_PARTTYPE, snapshot.particle.VEL)
        rs.write_block(f, "MASS", STAR_PARTTYPE, snapshot.particle.MASS)
        rs.write_block(f, "ID  ", STAR_PARTTYPE, snapshot.particle.ID)

    rs.closefile(f)


def get_snapshot_time(filename,snapnum,code="AREPO",snapinterval=0):
    if (code == "AREPO"):
        header = rs.snapshot_header(filename)
        time = rs.snapshot_header(filename).time
    elif (code == "PLUTO"):
        time = snapnum*snapinterval

    return time


def compute_snapshot_gradient(snapshot,fieldname = None, fielddata = None):

    # Read-in mesh-generating points
    points = np.array([snapshot.gas.POS[:,0],snapshot.gas.POS[:,1]]).T
    # Create Voronoi mesh
    mesh = voronoi_simulation_data.VoronoiMesh(points)
    # Compute the Voronoi-Green-Gauss gradients
    if (fieldname is not None):
        data = getattr(snapshot.gas,fieldname)
    elif (fielddata is not None):
        data = fielddata
    else:
        raise ValueError('You need to provide either "fieldname" or a data vector "fielddata"')
    gradientx,gradienty = voronoi_simulation_data.compute_voronoi_gradients(mesh,data)
    # Compute the slope limiters
    limiter = voronoi_simulation_data.compute_voronoi_limiter(mesh,data,np.array([gradientx,gradienty]).T)
    
    # Return limited gradients

    return np.array([gradientx*limiter, gradienty * limiter]).T

def compute_external_gravforce_from_snapshot(snapshot, XYZ = [0.0,0.0,0.0], softening=0.0):

    # Read-in mesh-generating points
    x,y,z = snapshot.gas.POS[:,0], snapshot.gas.POS[:,1], snapshot.gas.POS[:,2]

    forces = compute_external_gravforce(np.asarray([x,y,z]).T, source = XYZ, softening = softening)

    return forces 
    
    
def compute_external_gravforce(positions, source = [0.0,0.0,0.0], softening=0.0):

    dx = positions[:,0] - source[0]
    dy = positions[:,1] - source[1]
    dz = positions[:,2] - source[2]
    dR = np.sqrt(dx * dx + dy * dy + dz * dz)
    
    def softenedDistance(R,h):
        r2 = R * R
        fac = R * 0
        h_inv = 1.0 / h
        h3_inv = h_inv * h_inv * h_inv
        u = R * h_inv
        
        fac = 1.0 / (r2 * R)
        ind = u < 0.5
        fac[ind] = h3_inv * (10.666666666667 + u[ind] * u[ind] * (32.0 * u[ind] - 38.4))
        ind = (u >= 0.5) & (u < 1)
        fac[ind] = h3_inv * (21.333333333333 - 48.0 * u[ind] + 38.4 * u[ind] * u[ind] - 10.666666666667 * u[ind] * u[ind] * u[ind] - 0.066666666667 / \
                  (u[ind] * u[ind] * u[ind]))
        
        return fac

    oneoverR3 =  softenedDistance(dR,softening)

    force = -np.asarray([dx, dy, dz]).T *  oneoverR3[:,None]

    return force
