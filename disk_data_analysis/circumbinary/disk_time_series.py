import numpy as np
from orbital import orbit_in_time

accretion_variables = ['m1','m2','v1x_a','v1y_a','v2x_a','v2y_a',
                       'fext1x','fext1y','fext2x','fext2y',
                       'v1x_g','v1y_g','v2x_g','v2y_g']


def compute_time_derivative(x,time):

    '''
    Central difference approximation to the time derivative of an time-series array

    '''
    
    dt=np.gradient(accretion_time)
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
    
    data=np.loadtxt(accretionfile)

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
    
def create_binary_externalforces_files(accretionfile,outfilename1,
                                       outfilename2 = None,
                                       qb = 1.0, eb = 0.0,
                                       orbit_init = None,orbit_final = None):

    '''
    Compute force data and save to a file on disk.

    '''

    force_data = compute_binary_force_timeseries_from_accretionfile(accretionfile)
    
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

    if (outfilename is not None):
        np.savetxt(outfilename2,np.array([time,dm1dt,dm2dt]).T,
                   fmt='%12f %.8g %.8g')
        print "Saved accretion rate data to file:",outfilename2

