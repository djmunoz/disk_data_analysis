======
Computing maps of the torque distribution/density
======
.. sectnum::

   
Torque distribution map for a single snapshot
-----

First: Compute gravitational accelerations
~~~~~

Sometimes, your code output might have the gravitational
acceleration of cells already stored in. However, if you need the
back-reaction torque acting on, for example, a central binary, unless
the binary elements are a live part of the output, you might have to
recompute the accelerations due to the members of the binary
*separately*.


bye
....

.. code:: python

   import numpy as np
   import matplotlib.pyplot as plt
   import disk_data_analysis.circumbinary as dda
