import numpy as np
import sys
import os
import matplotlib.pyplot as plt

'''
outer_rad = 70.0
outer_density = 1.0238548701E-4
l=1.0
h= 0.1
alpha= 0.1
Mdot0 = -3 * np.pi * outer_density * alpha * h**2 * outer_rad**(1.5-l)
outer_accretion_coeff = 1.07
Mdot0 *= outer_accretion_coeff
#Mdot0 = 7.302e-06
'''

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


def read_binary_accretion_file(filename):
    """
    Read from disk a precomputed file with the accretion rates
    onto each component of the binary

    """
    accretion_data = np.loadtxt(filename)
    time = accretion_data[:,0]
    mdot1 = accretion_data[:,1]
    mdot2 = accretion_data[:,2]

    return time, mdot1, mdot2
    
def read_binary_forcing_file(filename,gravity_force = True, accretion_force = True):
    
    """ 
    Read from disk a precomputed file with the external forces
    acting on the binary

    """

    force_data = np.loadtxt(filename)
    time = force_data[:,0]
    x2 = force_data[:,1]
    y2 = force_data[:,2]
    x1 = force_data[:,3]
    y1 = force_data[:,4]
    vx2 = force_data[:,5]
    vy2 = force_data[:,6]
    vx1 = force_data[:,7]
    vy1 = force_data[:,8]
    fx2_acc = force_data[:,9]
    fy2_acc = force_data[:,10]
    fx1_acc = force_data[:,11]
    fy1_acc = force_data[:,12]
    fx2_grav = force_data[:,13]
    fy2_grav = force_data[:,14]
    fx1_grav = force_data[:,15]
    fy1_grav = force_data[:,16]

    return time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_grav,fy1_grav,fx2_grav,fy2_grav,fx1_acc,fy1_acc,fx2_acc,fy2_acc

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

def compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,mdot1,mdot2,
                                  fx1_ext,fy1_ext,fx2_ext,fy2_ext,
                                  G = 1.0 ,Mb = 1.0, ab = 1.0, qb = 1.0, eb = 0,etype=0):
    """
    Compute the change in both orbital energy and specific angular
    momentum of a binary subject to mass changes and external 
    torques

    """
    lb = np.sqrt(G * Mb * ab * (1 - eb * eb))
    Eb = -G * Mb / 2 / ab

    x,y = x2 - x1, y2 - y1
    vx,vy = vx2 - vx1, vy2 - vy1
    fxext, fyext = fx2_ext - fx1_ext, fy2_ext - fy1_ext
    mdot = mdot1 + mdot2
    qdot = qb * (1 + qb) / Mb * (mdot2 / qb - mdot1)
    lb = x * vy - y * vx
    r = np.sqrt(x * x + y * y)
    Eb = 0.5 * (vx * vx + vy * vy) - G * Mb / r
    Jb = qb * Mb / (1 + qb)**2 * lb

    Edot = -G * mdot / r + (vx * fxext + vy * fyext)

    ldot = compute_external_torques(x,y,fxext,fyext)

    Jdot = Jb * ( (1.0 - qb)/(1.0 + qb) * qdot / qb + mdot/Mb + ldot / lb)
    


    sinphi, cosphi = y/r, x/r
    fr_ext = fxext * cosphi + fyext * sinphi
    fphi_ext = -fxext * sinphi + fyext * cosphi



    if (eb == 0):
        #edot= np.sqrt(ab / G / Mb) * ((2 * cosphi * fphi_ext + sinphi * fr_ext + mdot / Mb * cosphi)
        edot = np.sqrt(ab * (1.0 - eb * eb)/ G / Mb) * ((eb + 2 * cosphi + eb * cosphi * cosphi) / (1 + eb * cosphi) * fphi_ext +\
                                                        sinphi * fr_ext) - \
                                                        mdot / Mb * (eb + cosphi)
        #adot = 2 * np.sqrt(ab**3 / G / Mb) * fphi_ext - ab * mdot / Mb
        adotovera = 2 * np.sqrt(ab / (1.0 - eb * eb)/ G / Mb) * (eb * sinphi * fr_ext + (1 + eb * cosphi) * fphi_ext) -\
                  mdot / Mb * (1 + 2 * eb * cosphi + eb * eb) / (1 - eb * eb)
    else:
        edot= 1.0/(G**2 * Mb**2 / lb**2 / Eb + 2) * np.sqrt(1 + 2 * lb**2 * Eb / G**2 / Mb**2) *  (2 * ldot / lb +  Edot / Eb - 2 * mdot / Mb)
        adotovera = ab * (-Edot / Eb + mdot / Mb)

    return Edot, ldot, Jdot, adotovera/(np.abs(Mdot0/Mb)), edot/(np.abs(Mdot0/Mb))



if __name__ == "__main__":

    qb = float(sys.argv[1])
    eb = float(sys.argv[2])

    h0 = 0.1
    alpha = 0.1
    Mb = 1.0

    orbit_init = 3499
    orbit_final = 3550
    

    acc_file_base = 'binary-accretion-rate_'
    tor_file_base = 'binary-forcing-rate_'
    if (eb == 0):
        run_tag = 'norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f.txt' % (orbit_init,orbit_final,qb,eb,h0,alpha)
    elif (np.floor(np.log10(eb)) == -1):
        run_tag = 'norbits%i-%i_q%.1f_e%.1f_h%.1f_alpha%.1f.txt' % (orbit_init,orbit_final,qb,eb,h0,alpha)
    elif (np.floor(np.log10(eb)) == -2):
        run_tag = 'norbits%i-%i_q%.1f_e%.2f_h%.1f_alpha%.1f.txt' % (orbit_init,orbit_final,qb,eb,h0,alpha)

    accretion_filename = acc_file_base+run_tag
    torque_filename = tor_file_base+run_tag

    print "Reading files..."
    print "                ",accretion_filename
    print "                ",torque_filename



    time,x1,y1,x2,y2,vx1,vy1,vx2,vy2,fx1_g,fy1_g,fx2_g,fy2_g,fx1_a,fy1_a,fx2_a,fy2_a = read_binary_forcing_file(torque_filename)
    _, mdot1, mdot2 =  read_binary_accretion_file(accretion_filename)
    mdot = mdot1 + mdot2
    fx1, fy1 = fx1_a + fx1_g,fy1_a + fy1_g
    fx2, fy2 = fx2_a + fx2_g,fy2_a + fy2_g
    
    torque_g = compute_external_torques(x2 - x1,y2 - y1,fx2_g - fx1_g, fy2_g - fy1_g)
    torque_a = compute_external_torques(x2 - x1,y2 - y1,fx2_a - fx1_a, fy2_a - fy1_a)
    
    fig  = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(411)
    ax.plot(time/2/np.pi,mdot/np.abs(Mdot0),color='k')
    ax.set_ylabel(r'$\dot{M}_{\rm b}/\dot{M}_0$')
    ax.text(0.85,0.05,'mean=%.3f' % (mdot/np.abs(Mdot0)).mean(),transform=ax.transAxes)
    ax.set_xlim(3500,3550)
    ax.set_ylim(0,2.0)
    #
    ax = fig.add_subplot(412)
    ax.plot(time/2/np.pi,torque_g/np.abs(Mdot0),color='steelblue')
    ax.set_ylabel(r'$\dot{l}_{\rm b,grav}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$')
    ax.text(0.85,0.05,'mean=%.3f' % (torque_g/np.abs(Mdot0)).mean(),transform=ax.transAxes)
    ax.set_xlim(3500,3550)
    #
    ax = fig.add_subplot(413)
    ax.plot(time/2/np.pi,torque_a/np.abs(Mdot0),color='royalblue')
    ax.set_ylabel(r'$\dot{l}_{\rm b,acc}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$')
    ax.text(0.85,0.05,'mean=%.3f' % (torque_a/np.abs(Mdot0)).mean(),transform=ax.transAxes)
    ax.set_xlim(3500,3550)
    #
    ax = fig.add_subplot(414)
    ax.plot(time/2/np.pi,(torque_a+torque_g)/np.abs(Mdot0),color='darkblue')
    #ax.plot(time/2/np.pi,compute_external_torques(x2 - x1,y2 - y1,fx2 - fx1, fy2 - fy1))
    ax.set_xlabel(r'$t/P_{\rm b}$',size=18)
    ax.set_ylabel(r'$\dot{l}_{\rm b}/\left[(\dot{M}_0/M_{\rm b})a_{\rm b}^2\Omega_{\rm b}\right]$')
    ax.text(0.85,0.05,'mean=%.3f' % ((torque_g+torque_a)/np.abs(Mdot0)).mean(),transform=ax.transAxes)
    ax.set_xlim(3500,3550)
    
    fig.savefig('./torque_evolution_q%.1f_e%.2f.pdf' % (qb,eb))
    plt.show()

    #plt.plot(time/2/np.pi,x1)
    #plt.plot(time /2/np.pi,x2)
    #x,y = x2 - x1, y2 - y1
    #vx,vy = vx2 - vx1, vy2 - vy1
    #r = np.sqrt(x * x + y * y)
    #sinphi, cosphi = y/r, x/r

    #Mb = 1.0
    #eb = 0.1
    #qb = 1.0
    #ecc_anom = np.vectorize(KeplerEquation)(time  + np.pi,eb)
    #x = np.cos(ecc_anom) - eb
    #y = np.sqrt(1 - eb * eb) * np.sin(ecc_anom)
    #vx = -np.sin(ecc_anom) / (1 - eb * np.cos(ecc_anom))
    #vy = np.sqrt(1 - eb * eb) * np.cos(ecc_anom)/ (1 - eb * np.cos(ecc_anom))
    #x1, x2 = -qb/(1 + qb) * x, 1.0/(1 + qb) * x
    #y1, y2 = -qb/(1 + qb) * y, 1.0/(1 + qb) * y
    #vx1, vx2 = -qb/(1 + qb) * vx, 1.0/(1 + qb) * vx
    #vy1, vy2 = -qb/(1 + qb) * vy, 1.0/(1 + qb) * vy
    #plt.plot(time /2/np.pi,x1)
    #plt.plot(time /2/np.pi,x2)
    #plt.show()

    ind = (time > (3501 * 2 * np.pi) ) & (time < (3547 * 2 * np.pi) )
    print mdot.mean(),mdot[ind].mean(),Mdot0
    Edot, ldot, Jdot, adot, edot = compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,mdot1,mdot2,fx1,fy1,fx2,fy2,eb=eb,etype=0)
    Edot, ldot, Jdot, adot0, edot0 = compute_binary_orbital_change(x1,y1,x2,y2,vx1,vy1,vx2,vy2,mdot1,mdot2,fx1,fy1,fx2,fy2,eb=eb,etype=1)

    #plt.plot(x1,y1,'bo')
    #plt.plot(x2,y2,'ro')

    print edot.mean(),adot.mean()
    print Jdot.mean()/mdot.mean()
    #plt.plot(time /2/np.pi,ldot)
    plt.plot(time /2/np.pi,edot,'b')
    plt.plot(time /2/np.pi,edot0,'r')
    #plt.plot(time /2/np.pi,mdot)
    plt.show()

    plt.plot(time /2/np.pi,adot)
    plt.plot(time /2/np.pi,adot0)
    #plt.plot(time /2/np.pi,mdot)
    plt.show()
    
    plt.plot(time[ind] /2/np.pi,mdot[ind]/np.abs(Mdot0))
    


    m1_init,m2_init =  2.1365034469e-07, 2.9898899901e-07
    m1_final,m2_final = 1.3752485758e-03, 1.3797994548e-03
    m_init, m_final = m1_init + m2_init, m1_final + m2_final

    print (m_final - m_init)/(time[-1]-time[0])
    plt.plot(time/2/np.pi,  np.repeat((m_final - m_init)/(time[-1]-time[0]),time.shape[0])/np.abs(Mdot0))
    plt.show()
    

    m1_init,m2_init = 6.5010669008e-03, 6.5016160482e-03
    m1_final,m2_final = 2.1133099476e-03, 2.1061396094e-03
    m_init, m_final = m1_init + m2_init, m1_final + m2_final
    t_init,t_final = 18846.6062682, 21991.1485751
    print (m_final - m_init)/(t_final-t_init)
