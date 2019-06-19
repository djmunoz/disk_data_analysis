import numpy as np
from scipy.integrate import quad
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t
cdef DTYPE_t pi = 3.1415926535897932384626433832795

cdef extern from "math.h":
    double sin(double)
    double cos(double)
    double tanh(double)
    double sqrt(double)
    double atan2(double,double)
    double acos(double)
    double abs(double)
    double log(double)
    double ceil(double)
    double cosh(double)
    double sinh(double)
    
def cosi_integrand(DTYPE_t y, DTYPE_t k, DTYPE_t z):
    """
    Integrand to Eq. (11) of Morton & Winn (2014)
    """
    return cosh(k*sqrt(1-y*y)) / sqrt(1-y*y) / sqrt(1-(z/y)*(z/y))

cpdef _cosi_pdf(DTYPE_t z, DTYPE_t kappa):
    """
    Equation (11) of Morton & Winn (2014)
    """
   
    return 2*kappa/(pi * sinh(kappa)) * quad(cosi_integrand,z,1,args=(kappa,z,))[0] 



