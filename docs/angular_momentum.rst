Computing angular momentum balance
======



.. code:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import disk_data_analysis.circumbinary as dda

   snap = dda.get_snapshot_data('./data/snap_',0,['POS','VEL','RHO','ACCE'])

   grid = dda.grid_cartesian(Xmin=-80.0,Xmax=80.0,Ymin=-80.0,Ymax=80.0,NX=1024,NY=1024,mask= '(R < 1.0) | (R > 80.0)')
   grid.X, grid.Y =  grid.X + snap.header.boxsize * 0.5, grid.Y + snap.header.boxsize * 0.5
   
   gradvx = dda.compute_gradient_on_grid(snap.gas.POS[:,0], snap.gas.POS[:,1],\
                                         snap.gas.VEL[:,0],grid)
   
   snap.add_snapshot_gradient('VEL','GRAV')
