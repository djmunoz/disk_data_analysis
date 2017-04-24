Computing angular momentum balance
======



.. code:: python

   import numpy as np
   import disk_data_analysis.circumbinary as dda

   snap = dda.get_snapshot_data('./data/snap_',0,['POS','VEL','RHO','ACCE'])

   snap.add_snapshot_gradient('VEL','GRAV')
