import numpy as np
from disk_data_analysis.orbital import orbit_in_time

accretion_variables = ['m1','m2','v1x_a','v1y_a','v2x_a','v2y_a',
                       'fext1x','fext1y','fext2x','fext2y',
                       'v1x_g','v1y_g','v2x_g','v2y_g']


def compute_time_derivative(x,time):

    '''
    Central difference approximation to the time derivative of an time-series array

    '''
    
    dt=np.gradient(time)
    dxdt = np.gradient(x)/dt

    return dxdt


def remove_discontinuity(x, time, skiptype = 0):

    '''
    Remove/skip jumpts or non-monotonic changes in the time variable due to errors in simulation output


    '''
    ind0 = np.argwhere(np.diff(time) <= 0)[0][0]

    if (skiptype == 0):
        xprime = np.append(x[:ind0],x[ind0+1:])
    elif (skiptype == 1):
        x0 = x[ind0]
        xprime = np.append(x[:ind0],x[ind0+1:]+x0)
    else:
        xprime = x
    
    return xprime

def compute_binary_force_timeseries_from_accretionfile(filename,variables=accretion_variables):

    '''
    Function to compute the accretion rate and forcing rate of a binary from mass growth and
    velocity changes

    '''

    read_variable = np.ones(len(accretion_variables)).astype(bool)
    for kk in range(len(accretion_variables)):
        if not (accretion_variables[kk] in variables):
            read_variable[kk] = False
    
    data=np.loadtxt(filename)

    time = data[:,0]
    count_variable = 1
    data_variables = []
    for kk in range(len(accretion_variables)):
        if (read_variable[kk] == True):
            data_variables.append(data[:,count_variable])
            count_variable+=1

    m1 = data[:,1]
    m2 = data[:,2]
    v1x_a = data[:,3]
    v1y_a = data[:,4]
    v2x_a = data[:,5]
    v2y_a = data[:,6]
    fext1x = data[:,7]
    fext1y = data[:,8]
    fext2x = data[:,9]
    fext2y = data[:,10]
    if (data.shape[1] == 15):
        v1x_g = data[:,11]
        v1y_g = data[:,12]
        v2x_g = data[:,13]
        v2y_g = data[:,14]


    #check that accretion times do not have discontinuities
    while (np.any(np.diff(time) <= 0)):
        m1 = remove_discontinuity(m1, time, skiptype = 1)
        m2 = remove_discontinuity(m2, time, skiptype = 1)
        v1x_a = remove_discontinuity(v1x_a, time, skiptype = 1)
        v1y_a = remove_discontinuity(v1y_a, time, skiptype = 1)
        v2x_a = remove_discontinuity(v2x_a, time, skiptype = 1)
        v2y_a = remove_discontinuity(v2y_a, time, skiptype = 1)
        fext1x = remove_discontinuity(fext1x, time, skiptype = 0)
        fext1y = remove_discontinuity(fext1y, time, skiptype = 0)
        fext2x = remove_discontinuity(fext2x, time, skiptype = 0)
        fext2y = remove_discontinuity(fext2y, time, skiptype = 0)
        if (data.shape[1] == 15):
            v1x_g= remove_discontinuity(v1x_g, time, skiptype = 1)
            v1y_g= remove_discontinuity(v1y_g, time, skiptype = 1)
            v2x_g= remove_discontinuity(v2x_g, time, skiptype = 1)
            v2y_g= remove_discontinuity(v2y_g, time, skiptype = 1)
        time = remove_discontinuity(time, time, skiptype = 0)
        

    # Take the time derivatives
    dm1dt = compute_time_derivative(m1,time)
    dm2dt = compute_time_derivative(m2,time)

    dv1xdt_a = compute_time_derivative(v1x_a,time)
    dv1ydt_a = compute_time_derivative(v1y_a,time)
    dv2xdt_a = compute_time_derivative(v2x_a,time)
    dv2ydt_a = compute_time_derivative(v2y_a,time)

    if (data.shape[1] == 15):
        dv1xdt_g = compute_time_derivative(v1x_g,time)
        dv1ydt_g = compute_time_derivative(v1y_g,time)
        dv2xdt_g = compute_time_derivative(v2x_g,time)
        dv2ydt_g = compute_time_derivative(v2y_g,time)
    else:
        dv1xdt_g = fext1x
        dv1ydt_g = fext1y
        dv2xdt_g = fext2x
        dv2ydt_g = fext2y


    force_data = np.asarray([time,dm1dt,dm2dt,
                             dv1xdt_a,dv1ydt_a,dv2xdt_a,dv2ydt_a,
                             dv1xdt_g,dv1ydt_g,dv2xdt_g,dv2ydt_g]).T

    return force_data
    
def write_binary_externalforces_file(accretionfile,outfilename1,
                                     outfilename2 = None,
                                     qb = 1.0, eb = 0.0,
                                     orbit_init = None,orbit_final = None, maxlines = None):

    '''
    Compute force data and save to a file on disk.

    '''

    force_data = compute_binary_force_timeseries_from_accretionfile(accretionfile)
    if (maxlines is not None):
        if (force_data.shape[0] > 1.5 * maxlines):
            skip = int(np.ceil(force_data.shape[0]/(1.5 * maxlines)))
            force_data = force_data[::skip,:]
            
    time = force_data[:,0]
    if (orbit_init is None):
        orbit_init = int(time[0]/2/np.pi)
    if (orbit_final is None):
        orbit_final = int(time[-1]/2/np.pi)
        
    time_min, time_max = orbit_init * 2 * np.pi, orbit_final * 2 * np.pi
    ind = (time >= time_min) & (time <= time_max) 
    time = time[ind]
    
    dm1dt = force_data[:,1][ind]
    dm2dt = force_data[:,2][ind]
    dv1xdt_a = force_data[:,3][ind]
    dv1ydt_a = force_data[:,4][ind]
    dv2xdt_a = force_data[:,5][ind]
    dv2ydt_a = force_data[:,6][ind]
    dv1xdt_g = force_data[:,7][ind]
    dv1ydt_g = force_data[:,8][ind]
    dv2xdt_g = force_data[:,9][ind]
    dv2ydt_g = force_data[:,10][ind]

    
    mu =  qb / ( 1.0 + qb)
    xb, yb, vxb, vyb = np.vectorize(orbit_in_time)(time + np.pi, eb)
    x1, y1 = mu * xb, mu * yb
    x2, y2 = -(1.0 - mu) * xb, -(1.0 - mu) * yb
    vx1, vy1 = mu * vxb, mu * vyb
    vx2, vy2 = -(1.0 - mu) * vxb, -(1.0 - mu) * vyb


    np.savetxt(outfilename1,np.array([time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,
                                      dv1xdt_a,dv1ydt_a,dv2xdt_a,dv2ydt_a,
                                      dv1xdt_g,dv1ydt_g,dv2xdt_g,dv2ydt_g]).T,
               fmt='%12f %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g')
  

    print "Saved accretion rate data to file:",outfilename1

    if (outfilename2 is not None):
        np.savetxt(outfilename2,np.array([time,dm1dt,dm2dt]).T,
                   fmt='%12f %.8g %.8g')
        print "Saved accretion rate data to file:",outfilename2


    return



def read_binary_externalforces_file(forcefilename,
                                    qb = 1.0, eb = 0.0,
                                    orbit_init = None,orbit_final = None, maxlines = None):
    '''
    Read from disk a precomputed file with the external forces
    acting on the binary

    '''

    force_data = np.loadtxt(forcefilename)
    time = force_data[:,0]
    x1 = force_data[:,1]
    y1 = force_data[:,2]
    x2 = force_data[:,3]
    y2 = force_data[:,4]
    vx1 = force_data[:,5]
    vy1 = force_data[:,6]
    vx2 = force_data[:,7]
    vy2 = force_data[:,8]
    fx1_a = force_data[:,9]
    fy1_a = force_data[:,10]
    fx2_a = force_data[:,11]
    fy2_a = force_data[:,12]
    fx1_g = force_data[:,13]
    fy1_g = force_data[:,14]
    fx2_g = force_data[:,15]
    fy2_g = force_data[:,16]
    

    return time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_a,fy1_a,fx2_a,fy2_a,fx1_g,fy1_g,fx2_g,fy2_g


def read_binary_accretion_file(accretionfilename,
                               qb = 1.0, eb = 0.0,
                               orbit_init = None,orbit_final = None, maxlines = None):
    '''
    Read from disk a precomputed file with the accretion rates
    onto each component of the binary

    '''

    
    accretion_data = np.loadtxt(accretionfilename)
    time = accretion_data[:,0]
    mdot1 = accretion_data[:,1]
    mdot2 = accretion_data[:,2]

    return time, mdot1, mdot2
    
def compute_external_torques(x,y,fx,fy):
    """
    Compute torque due to external forces on a binary orbit given
    the position of both elements of the binary and the forces acting
    on each of those two elements

    x,y: numpy arrays -- (x,y) position of the binary relative coordinate

    fx,fy: numpy arrays -- x- and y-components of the net external force 
    **per unit mass** acting on the binary

    """
    
    t =  x * fy - y * fx
    
    return t


def compute_binary_torque_contributions_from_files(accretionfilename,forcefilename,
                                                   qb=1.0,eb=0.0,orbit_init = None, orbit_final = None):
    '''
    Function to compute the different contributions to binary torque
    from a file of precomputed quantities.
    
    '''

    # Read quantities
    time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_a,fy1_a,fx2_a,fy2_a,fx1_g,fy1_g,fx2_g,fy2_g = read_binary_externalforces_file(forcefilename,
                                                                                                                       qb=qb,eb=eb,
                                                                                                                       orbit_init = orbit_init,
                                                                                                                       orbit_final = orbit_final)

    _, mdot1, mdot2 =  read_binary_accretion_file(accretionfilename,
                                                  qb=qb,eb=eb,
                                                  orbit_init = orbit_init,
                                                  orbit_final = orbit_final)
    
    # Combine data
    mdot = mdot1 + mdot2
    qdot = (1 + qb) * (mdot2 - qb * mdot1)
    fx1, fy1 = fx1_a + fx1_g, fy1_a + fy1_g
    fx2, fy2 = fx2_a + fx2_g, fy2_a + fy2_g
    
    torque_grav = compute_external_torques(x1 - x2, y1 - y2, fx1_g - fx2_g, fy1_g - fy2_g)
    torque_acc = compute_external_torques(x1 - x2, y1 - y2, fx1_a - fx2_a, fy1_a - fy2_a)
    
    reduced_mass = qb * 1.0 / (1 + qb) / (1 + qb)
    lb = np.sqrt(1 - eb**2)
    Jb = reduced_mass * lb
    Jdot = Jb * ((1.0 - qb)/(1.0 + qb) * qdot /qb  + mdot  + (torque_acc+torque_grav) / lb)

    return mdot,qdot,torque_grav,torque_acc,Jdot


def compute_binary_angular_momentum_change(x1,y1,x2,y2,mdot1,mdot2,
                                           fx1_ext,fy1_ext,fx2_ext,fy2_ext,
                                           qb = 1.0, eb = 0,
                                           G = 1.0 ,Mb = 1.0, ab = 1.0):
    '''
    Compute the change in both orbital energy and specific angular
    momentum of a binary subject to mass changes and external 
    torques
    '''


    # Combine data
    mdot = mdot1 + mdot2
    qdot = (1 + qb) * (mdot2 - qb * mdot1)
    fx1, fy1 = fx1_a + fx1_g, fy1_a + fy1_g
    fx2, fy2 = fx2_a + fx2_g, fy2_a + fy2_g
    
    ldot = compute_external_torques(x1 - x2, y1 - y2, fx1_ext - fx2_ext, fy1_ext - fy2_ext)
    
    reduced_mass = qb * 1.0 / (1 + qb) / (1 + qb) * Mb
    lb = np.sqrt((1 - eb**2) * G * Mb *  ab)
    Jb = reduced_mass * lb
    Jdot = Jb * ((1.0 - qb)/(1.0 + qb) * qdot /qb  + mdot/Mb  + ldot / lb)

    return mdot, qdot, ldot, Jdot



def compute_binary_orbital_change(x1,y1,x2,y2,mdot1,mdot2,
                                  fx1_ext,fy1_ext,fx2_ext,fy2_ext,
                                  G = 1.0 ,Mb = 1.0, ab = 1.0, qb = 1.0, eb = 0):
    '''
    Compute the change in both orbital energy and specific angular
    momentum of a binary subject to mass changes and external 
    torques

    '''
    
    lb = np.sqrt(G * Mb * ab * (1 - eb * eb))
    Eb = -G * Mb / 2 / ab

    x,y = x1 - x2, y1 - y2
    fxext, fyext = fx1_ext - fx2_ext, fy1_ext - fy2_ext
    mdot = mdot1 + mdot2
    qdot = qb * (1 + qb) / Mb * (mdot2 / qb - mdot1)
    r = np.sqrt(x * x + y * y)

    Jb = qb * Mb / (1 + qb)**2 * lb

    Edot = -G * mdot / r + (vx * fxext + vy * fyext)

    ldot = compute_external_torques(x,y,fxext,fyext)

    Jdot = Jb * ( (1.0 - qb)/(1.0 + qb) * qdot / qb + mdot/Mb + ldot / lb)
    


    sinphi, cosphi = y/r, x/r
    fr_ext = fxext * cosphi + fyext * sinphi
    fphi_ext = -fxext * sinphi + fyext * cosphi



    if (eb == 0):
        edot = np.sqrt(ab * (1.0 - eb * eb)/ G / Mb) * ((eb + 2 * cosphi + eb * cosphi * cosphi) / (1 + eb * cosphi) * fphi_ext +\
                                                        sinphi * fr_ext) - \
                                                        mdot / Mb * (eb + cosphi)
        adotovera = 2 * np.sqrt(ab / (1.0 - eb * eb)/ G / Mb) * (eb * sinphi * fr_ext + (1 + eb * cosphi) * fphi_ext) -\
                  mdot / Mb * (1 + 2 * eb * cosphi + eb * eb) / (1 - eb * eb)
    else:
        edot= 1.0/(G**2 * Mb**2 / lb**2 / Eb + 2) * np.sqrt(1 + 2 * lb**2 * Eb / G**2 / Mb**2) *  (2 * ldot / lb +  Edot / Eb - 2 * mdot / Mb)
        adotovera = ab * (-Edot / Eb + mdot / Mb)

    return Edot, ldot, Jdot, adotovera, edot
