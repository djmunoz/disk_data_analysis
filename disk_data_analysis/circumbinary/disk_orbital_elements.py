import numpy as np





def compute_disk_eccentriciy_vector(pos,vel,GM=1,coordinates='cartesian'):

    if (coordinates=='cartesian'):
        x , y , z = pos[:,0], pos[:,1], pos[:,2]
        r = np.sqrt(x**2+y**2+z**2)
        vx , vy ,vz = vel[:,0], vel[:,1], vel[:,2]
        ex = (vy * (x * vy - y * vx) - vz * (z * vx - x * vz)) / GM - x/r
        ey = (vz * (y * vz - z * vy) - vx * (x * vy - y * vx)) / GM - y/r
        ez = (vx * (z * vx - x * vz) - vy * (y * vz - z * vy)) / GM - z/r


    elif (coordinates=='polar'):
        r , phi , z = pos[:,0], pos[:,1], pos[:,2]
        vr , vphi ,vz = vel[:,0], vel[:,1], vel[:,2]        
        er = r * (vphi**2 + vz**2) - 1
        ephi = -r * vr * vphi
        ez = -r * vr * vz
        ex = er * np.cos(phi) - ephi * np.sin(phi) 
        ey = er * np.sin(phi) + ephi * np.cos(phi) 
           
    e = np.array([ex,ey,ez]).T

    return e


def compute_disk_angular_momentum_vector(pos,vel,GM=1,coordinates='cartesian'):

    if (coordinates=='cartesian'):
        x , y , z = pos[:,0], pos[:,1], pos[:,2]
        vx , vy ,vz = vel[:,0], vel[:,1], vel[:,2]
        jx = y * vz - z * vy
        jy = z * vx - x * vz
        jz = x * vy - y * vx

    j = np.array([jx,jy,jz]).T

    return j

def compute_disk_semimaj(pos,vel,GM=1,coordinates='cartesian'):
    if (coordinates=='cartesian'):
        x , y , z = pos[:,0], pos[:,1], pos[:,2]
        r = np.sqrt(x**2+y**2+z**2)
        vx , vy ,vz = vel[:,0], vel[:,1], vel[:,2]
        v2 = vx**2+vy**2+vz**2

    elif (coordinates=='polar'):
        r , phi , z = pos[:,0], pos[:,1], pos[:,2]
        vr , vphi ,vz = vel[:,0], vel[:,1], vel[:,2]   
        v2 = vr**2+vphi**2+vz**2
           
    a = np.array(1.0/(2.0/r - v2/GM))

    return a



