from numpy import sqrt, cos, sin, tan, arccos, arctan2, pi, cosh, sinh, arctanh, abs, sign, e
import matplotlib.pyplot as plt
import numpy.random as rd


G = 1

def elements(x,y,z,vx,vy,vz,mu):

  R2 = x**2 + y**2 + z**2
  R = sqrt(R2)
  V2 = vx**2 + vy**2 + vz**2
  RtimesRdot = x * vx + y * vy + z * vz
  hx = y * vz - z * vy
  hy = z * vx - x * vz
  hz = x * vy - y * vx
  h2 = hx**2 + hy**2 + hz**2
  if (RtimesRdot > 0 ):
    Rdot = sqrt(V2 - h2/R2) 
  else:
    Rdot = -sqrt(V2 - h2/R2) 

  #eccentricity and pericenter distance
  mu_1 = 1.0/mu
  temp = 1.0  +  h2 * mu_1 * (V2 *mu_1  -  2.0 / R)
  if (temp <= 0): ecc = 0.0
  else: ecc = sqrt(temp)

  if (ecc < 1e-8): ecc = 1.e-8

  peridist = h2 * mu_1 / (1.0 + ecc)
  semimaj = 1.0 / (2.0/R - V2 * mu_1)

  #inclination
  incl = arccos(hz/sqrt(h2))
  if (incl != 0.0):
    if (hz > 0): node = arctan2(hx /sqrt(h2)/sin(incl),
                                   -hy/sqrt(h2)/sin(incl))
    else: node = arctan2(-hx /sqrt(h2)/sin(incl),
                            hy/sqrt(h2)/sin(incl))
  else:
    node = 0.0

  #true longitude (argument of pericenter plus true anomaly)   
  if ((incl > 1.e-3) & (incl < pi-1.0e-3)):
    sinomegaplusf = z/R/sin(incl)
    cosomegaplusf = 1.0/cos(node)*(x/R + sin(node) * sinomegaplusf *
                                      cos(incl))
    
    periargplustrue_anom = arctan2(sinomegaplusf,cosomegaplusf)
  else:
    periargplustrue_anom = arctan2(y,x)*cos(incl)


  #true anomaly and argument of pericenter
  true_anom = arctan2(semimaj*(1.0 -ecc**2)/sqrt(h2)/ecc * Rdot,
                         1.0/ecc*(semimaj/R*(1.0 - ecc**2) - 1.0))
  periarg = periargplustrue_anom - true_anom
    

  periarg = periarg % (2*pi)
  true_anom = true_anom % (2*pi)
  if (true_anom < 0): true_anom = 2.0*pi + true_anom
  if (periarg < 0): periarg = 2.0*pi + periarg


  if (ecc < 1.0):
    tanecc_anom_2 = tan(0.5 * true_anom)* sqrt((1.0 - ecc)/(1.0 + ecc))
    tanecc_anom = 2.0 * tanecc_anom_2 / (1.0 - tanecc_anom_2**2)
    cosecc_anom = (1.0 - R / semimaj )/ ecc
    ecc_anom = arctan2(tanecc_anom * cosecc_anom,cosecc_anom)
    if (ecc_anom < 0): ecc_anom = 2.0*pi + ecc_anom
    
    mean_anom = ecc_anom - ecc * sin(ecc_anom)
    
    return peridist,ecc,incl,periarg,node,true_anom,mean_anom,ecc_anom
      
  else:
    tanhhyp_anom_2 = tan(0.5 * true_anom)* sqrt((ecc - 1.0)/(ecc + 1.0))
    tanhhyp_anom = 2.0 * tanhhyp_anom_2 / (tanhhyp_anom_2**2 + 1.0)
    hyp_anom = arctanh(tanhhyp_anom)
    
    mean_anom = ecc * sinh(hyp_anom) - hyp_anom
    
    return peridist,ecc,incl,periarg,node,true_anom,mean_anom,hyp_anom


###################################################################################
def orbit(peri,e,I,omega,Omega,M,mu):

  
  
  if (e < 1):
    a = peri/(1.0 - e)
    ecc_anom = keplerEquation(M,e)
    x = a * (cos(ecc_anom) - e)
    y = a * sqrt(1 - e * e) * sin(ecc_anom)
    xdot = -sqrt(mu/a) / (1.0 - e*cos(ecc_anom))* sin(ecc_anom)
    ydot = sqrt(mu/a) / (1.0 - e*cos(ecc_anom))* cos(ecc_anom) * sqrt(1.0 - e * e)
    f = arctan2(sqrt(1.0 - e * e) * sin(ecc_anom),cos(ecc_anom) - e)
    radius = a*(1.0 - e*cos(ecc_anom))
    
  elif (e > 1):
    a = peri/(1.0 - e)
    hyp_anom = hyperbolicKeplerEquation(M,e)
    x = a * (cosh(hyp_anom) - e)
    y = -a * sqrt(e * e - 1.0) * sinh(hyp_anom)
    xdot = -sqrt(mu/-a) / (e*cosh(hyp_anom)- 1.0)* sinh(hyp_anom)
    ydot = sqrt(mu/-a) / (e*cosh(hyp_anom) - 1.0)* cosh(hyp_anom) * sqrt(e * e - 1.0)
    f = arctan2(sqrt(e * e - 1.0) * sinh(hyp_anom),e - cosh(hyp_anom))
    radius = a*(1.0 - e * cosh(hyp_anom))
    
  elif (e == 1):
    tan_trueanom_2 = parabolicKeplerEquation(M)
    x = peri * (1.0 - tan_trueanom_2**2)
    y = 2.0* peri * tan_trueanom_2
    xdot = -sqrt(2.0 * mu / peri) / (1.0 + tan_trueanom_2**2) * tan_trueanom_2
    ydot =  sqrt(2.0 * mu / peri) / (1.0 + tan_trueanom_2**2)
    f = arctan2(2.0 * tan_trueanom_2/(1.0 + tan_trueanom_2**2),(1.0 - tan_trueanom_2**2)/(1.0 + tan_trueanom_2**2))
    radius = peri * (1.0 + tan_trueanom_2**2)
        

  #rotation matrix
  d11 =  cos(omega) * cos(Omega) - sin(omega) * sin(Omega) * cos(I)
  d12 =  cos(omega) * sin(Omega) + sin(omega) * cos(Omega) * cos(I)
  d13 =  sin(omega) * sin(I)
  d21 = -sin(omega) * cos(Omega) - cos(omega) * sin(Omega) * cos(I)
  d22 = -sin(omega) * sin(Omega) + cos(omega) * cos(Omega) * cos(I)
  d23 =  cos(omega) * sin(I)

  X = d11 * x + d21 * y
  Y = d12 * x + d22 * y
  Z = d13 * x + d23 * y
  VX = d11 * xdot + d21 * ydot
  VY = d12 * xdot + d22 * ydot
  VZ = d13 * xdot + d23 * ydot
      
  return X, Y, Z, VX, VY, VZ

###################################################################################3
def keplerEquation(mean_anom,e):
 
  k = 0.85
  iter = 0
  abstol = 1.0e-8
  reltol = 1.0e-8
  


  while (abs(mean_anom) > 2.0*pi):
    mean_anom-=2.0*pi*sign(mean_anom)
  if (mean_anom < 0.0): mean_anom = 2*pi + mean_anom

  k = 0.85
  ecc_anom = mean_anom + sign(sin(mean_anom))* k * e
  #ecc_anom = mean_anom
  if (e > 0.8):
    ecc_anom = pi
    
  while(True):
    f = ecc_anom -e * sin(ecc_anom) - mean_anom
    fprime = 1.0 - e * cos(ecc_anom)
    fprime2 = e * sin(ecc_anom)
    fprime3 = e * cos(ecc_anom)
    delta1 = - f / fprime
    delta2 = - f /(fprime + 0.5 * delta1 * fprime2)
    delta3 = - f /(fprime + 0.5 * delta2 * fprime2 + 0.16666666666 * delta2**2 * fprime3) 

    if (delta3 == 0): break
    
    if (abs(ecc_anom) > 0.0):
      abserr,relerr = abs(delta3),abs(delta3)/abs(ecc_anom)
    else:
      abserr,relerr = abs(delta3),1.0e40

    ecc_anom+=delta3
    #print iter,ecc_anom,e,delta3
    
    if (abs(ecc_anom) > abstol/reltol):
      if (abserr < abstol): break
    else:
      if (relerr < reltol): break
    iter+=1      

  return ecc_anom % (2*pi)

###################################################################################3
def hyperbolicKeplerEquation(mean_anom,e):


  if (mean_anom > 0):
    hyp_anom = log(mean_anom/e + 1.8)
  else:
    hyp_anom = -log(abs(mean_anom)/e + 1.8)
  abstol,reltol = 1.0e-10, 1.0e-10
  iter = 0
  while(True):
    f = e * sinh(hyp_anom) - hyp_anom - mean_anom
    fprime =  e * cosh(hyp_anom) - 1.0
    fprime2 = e * sinh(hyp_anom) 
    fprime3 = e * cosh(hyp_anom)
    delta1 = - f / fprime
    delta2 = - f /(fprime + 0.5 * delta1 * fprime2)
    delta3 = - f /(fprime + 0.5 * delta2 * fprime2 + 0.16666666666 * delta2**2 * fprime3) 

    if (delta3 == 0): break
    
    if (abs(hyp_anom) > 0.0):
      abserr,relerr = abs(delta3),abs(delta3)/abs(hyp_anom)
    else:
      abserr,relerr = abs(delta3),1.0e40
      
    hyp_anom+=delta3
    #print iter,hyp_anom,e,delta3
    
    if (abs(hyp_anom) > abstol/reltol):
      if (abserr < abstol): break
    else:
      if (relerr < reltol): break
    iter+=1

      
  return hyp_anom


##################################################################
def parabolicKeplerEquation(mean_anom): #jerome cardan's method

  B = 1.5 * mean_anom
  tan_trueanom_2 = (B + sqrt(1 + B *  B))**(1.0/3) - 1.0/(B + sqrt(1 + B *  B))**(1.0/3)
  return tan_trueanom_2


##################################################################
def semimajorAxis(x,y,z,vx,vy,vz,mu):

    mu_1 = 1.0/mu
    R2 = x**2 + y**2 + z**2
    R = sqrt(R2)
    V2 = vx**2 + vy**2 + vz**2
    semimaj = 1.0 / (2.0/R - V2 * mu_1)
    
    return semimaj

def eccentricity(x,y,z,vx,vy,vz,mu):

    R2 = x**2 + y**2 + z**2
    R = sqrt(R2)
    V2 = vx**2 + vy**2 + vz**2
    RtimesRdot = x * vx + y * vy + z * vz
    hx = y * vz - z * vy
    hy = z * vx - x * vz
    hz = x * vy - y * vx
    h2 = hx**2 + hy**2 + hz**2
    if (RtimesRdot > 0 ):
        Rdot = sqrt(V2 - h2/R2) 
    else:
        Rdot = -sqrt(V2 - h2/R2) 
        
    mu_1 = 1.0/mu
    temp = 1.0  +  h2 * mu_1 * (V2 *mu_1  -  2.0 / R)
    if (temp <= 0): eccentricity = 0.0
    else: eccentricity = sqrt(temp)

    return eccentricity


##################################################################


def orbit_in_time(time,eccentricity):

    ecc_anom = keplerEquation(time,eccentricity)
    x = cos(ecc_anom) - eccentricity
    y = sqrt(1 - eccentricity * eccentricity) * sin(ecc_anom)
    vx = -sin(ecc_anom) / (1 - eccentricity * cos(ecc_anom))
    vy = sqrt(1 - eccentricity * eccentricity) * cos(ecc_anom) / \
         (1 - eccentricity * cos(ecc_anom))

    return x, y, vx, vy
