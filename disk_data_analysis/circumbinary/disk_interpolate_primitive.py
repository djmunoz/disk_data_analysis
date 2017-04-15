import numpy as np
from scipy.interpolate import griddata




def disk_interpolate_primitive_quantities(snapshot, RPHI, quantities = None):

    """
    Interpolate primitive quantity data (stored in a snapshot data structure) usually 
    evaluated in spatial locations placed in an unstructured fashion into a structured grid

    """


    
    if (quantities is None):
        quantities = ["RHO"]

        
    R,PHI = RPHI[0],RPHI[1]
    x,y = snapshot.gas.POS[:,0] - snapshot.header.boxsize * 0.5, snapshot.gas.POS[:,1] - snapshot.header.boxsize * 0.5
    rvals, phivals = np.sqrt(x**2 + y**2), np.arctan2(y,x)
    
    
    interp_quant = []
    for quant in quantities:

        z = getattr(snapshot.gas,quant)
        # interpolate Z values on defined grid
        Z = griddata(np.vstack((rvals.flatten(),phivals.flatten())).T, \
                          np.vstack(z.flatten()),(R,PHI),method='linear').reshape(R.shape)
        interp_quant.append(Z)

    return interp_quant
        
        
