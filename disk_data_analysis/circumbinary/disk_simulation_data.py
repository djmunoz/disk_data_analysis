from disk_hdf5 import readsnapHDF5 as rs
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
              }

def get_snapshot_data(filename_prefix,snap_num,quantities,code="AREPO"):
    print code,"hehe"
    print quantities
    nquant=len(quantities)
    outquant = []
    for quant in quantities:
        if (code == "AREPO"):
            if (quant.ljust(4) in datablocks):
                outquant.append(rs.read_block(filename_prefix+str(snap_num).zfill(3),quant.ljust(4),parttype=0))
            elif (quant in disk_data_fields):
                header = rs.snapshot_header(filename_prefix+snap_num)
                BoxX,BoxY = header.boxsize,header.boxsize
                if ((quant == "R") | (quant == "PHI")):
                    pos = rs.read_block(filename_prefix+snap_num,"POS ", parttype=0)
                    pos[:,0],pos[:,1] = (pos[:,0]- BoxX/2), (pos[:,1]- BoxY/2)
                    radius = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
                    phi = np.arctan2(pos[:,1],pos[:,0])
                    if (quant == "R"):
                        outquant.append(radius)
                    if (quant == "PHI"):
                        outquant.append(phi)
                if ((quant == "VELR") | (quant == "VELPHI")):
                    pos = rs.read_block(filename_prefix+snap_num,"POS ", parttype=0)
                    vel = rs.read_block(filename_prefix+snap_num,"VEL ", parttype=0)
                    pos[:,0],pos[:,1] = (pos[:,0]- BoxX/2), (pos[:,1]- BoxY/2)
                    radius = np.sqrt(pos[:,0]**2 + pos[:,1]**2)
                    phi = np.arctan2(pos[:,1],pos[:,0])
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
                    print pos.shape
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

    return outquant


def get_snapshot_time(filename,snapnum,code="AREPO",snapinterval=0):
    if (code == "AREPO"):
        header = rs.snapshot_header(filename)
        time = rs.snapshot_header(filename).time
    elif (code == "PLUTO"):
        time = snapnum*snapinterval

    return time
