




.. code:: python

   import disk_data_analysis.circumbinary as dda
   from scipy.spatial import Voronoi, voronoi_plot_2d

   # Load the simulation data and, for the sake of simplicity, select a subset of cells
   snap = dda.get_snapshot_data('./data/snap_',0,['POS','RHO','R','VOL'])
   ind = snap.gas.R < 3.8
   points = np.array([snap.gas.POS[ind,0],snap.gas.POS[ind,1]]).T
   quant = snap.gas.RHO[ind]
   cellvol = snap.gas.VOL[ind]

   # Compute the Voronoi tessellation...
   vor = Voronoi(points)

   # ... and compute some quantities associated with the mesh
   rvectors = points[vor.ridge_points[:,1]]-points[vor.ridge_points[:,0]]
   rvecsize = np.sqrt(rvectors[:,0]**2 + rvectors[:,1]**2)

   cvectors = 0.5 * (vor.vertices[np.asarray(vor.ridge_vertices)[:,1]]+vor.vertices[np.asarray(vor.ridge_vertices)[:,0]]) - 0.5 *( points[vor.ridge_points[:,1]]+points[vor.ridge_points[:,0]])

   faceareas = vor.vertices[np.asarray(vor.ridge_vertices)[:,1]] - vor.vertices[np.asarray(vor.ridge_vertices)[:,0]]
   faceareas = np.sqrt(faceareas[:,0]**2 + faceareas[:,1]**2)

   delta_quant = quant[vor.ridge_points[:,1]]- quant[vor.ridge_points[:,0]]
   mean_quant = 0.5 * (quant[vor.ridge_points[:,1]] + quant[vor.ridge_points[:,0]])

   # We use these quantities to compute the Voronoi-Green-Gauss gradients
   gradient_sum_0 = faceareas[:] * (delta_quant[:] * cvectors[:,0] - mean_quant[:] * rvectors[:,0]) / rvecsize[:]
   gradient_sum_1 = faceareas[:] * (delta_quant[:] * cvectors[:,1] - mean_quant[:] * rvectors[:,1]) / rvecsize[:]

	  
   def sum_by_group(values, groups):
	  order = np.argsort(groups)
	  groups = groups[order]
	  values = values[order]
	  values.cumsum(out=values)
	  index = np.ones(len(groups), 'bool')
	  index[:-1] = groups[1:] != groups[:-1]
	  values = values[index]
	  groups = groups[index]
	  values[1:] = values[1:] - values[:-1]
	  return values, groups
   
   print  sum_by_group(gradient_sum_0,vor.ridge_points[:,0])
	  
