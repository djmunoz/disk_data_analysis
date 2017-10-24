from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np

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

def group_by_group(values, groups):
       order = np.argsort(groups)
       groups = groups[order]
       values = values[order]
       index = np.ones(len(groups), 'bool')
       index[:-1] = groups[1:] != groups[:-1]
       values = np.split(values,np.argwhere(index == True).flatten()[:-1]+1)
       groups = groups[index]
       values = np.asarray([[max(element),min(element)] for element in values])
       return values[:,0],values[:,1], groups

def combine_left_and_right(leftquant,rightquant,maskleft,maskright,kind='sum'):
    if (kind == 'sum'):
        combined = rightquant[maskright] + leftquant[maskleft]
    elif (kind == 'min'):
        combined = np.minimum(rightquant[maskright],leftquant[maskleft])
    elif (kind == 'max'):
        combined = np.maximum(rightquant[maskright],leftquant[maskleft])
           
    return np.append(np.append(combined,rightquant[~maskright]),leftquant[~maskleft])



class Cells():

    def __init__(self,*args,**kwargs):
        self.pos = kwargs.get("pos")
        self.centroids = kwargs.get("centroids")
        self.vols = kwargs.get("vols")
        self.areas = kwargs.get("areas")

class Faces():

    def __init__(self,*args,**kwargs):
        self.centroids = kwargs.get("centroids")
        self.areas = kwargs.get("areas")
        self.midpoints = kwargs.get("midpoints")
        self.rvectors = kwargs.get("rvectors")
        self.left_point_ind = kwargs.get("left_point_ind")
        self.right_point_ind = kwargs.get("right_point_ind") 
        
class VoronoiMesh():

    def __init__(self,points):
        vor = Voronoi(points)
        centroids, areas, midpoints, rvectors, left_point_ind, right_point_ind = self.compute_face_properties(vor)
        self.faces = Faces(centroids = centroids, areas = areas, \
                           midpoints = midpoints, rvectors = rvectors,\
                           left_point_ind = left_point_ind,right_point_ind = right_point_ind)
        centroids, vols, areas = self.compute_cell_properties(vor)
        centroids[vols <= 0,:] = points[vols <= 0,:] 
        self.cells = Cells(pos=points,centroids = centroids, vols=vols,areas=areas)


    def compute_face_properties(self,vor):

        rvectors = vor.points[vor.ridge_points[:,0]] - vor.points[vor.ridge_points[:,1]]
        faceareas = vor.vertices[np.asarray(vor.ridge_vertices)[:,1]] - vor.vertices[np.asarray(vor.ridge_vertices)[:,0]]
        faceareas = np.sqrt(faceareas[:,0]**2 + faceareas[:,1]**2)    
        facemidpoints =  0.5 * (vor.points[vor.ridge_points[:,1]]+vor.points[vor.ridge_points[:,0]])
        facecentroidsx =  0.5 * (vor.vertices[np.asarray(vor.ridge_vertices)[:,1]][:,0]+vor.vertices[np.asarray(vor.ridge_vertices)[:,0]][:,0])
        facecentroidsy =  0.5 * (vor.vertices[np.asarray(vor.ridge_vertices)[:,1]][:,1]+vor.vertices[np.asarray(vor.ridge_vertices)[:,0]][:,1])
        facecentroids = np.array([facecentroidsx,facecentroidsy]).T

        return facecentroids, faceareas, facemidpoints, rvectors, vor.ridge_points[:,0], vor.ridge_points[:,1]
        
    def compute_cell_properties(self,vor):

        # indices connecting faces to cell centers
        left_ind, right_ind = np.unique(np.sort(vor.ridge_points[:,0])), np.unique(np.sort(vor.ridge_points[:,1]))
        index_overlap  = np.intersect1d(left_ind,right_ind)
        mask_right, mask_left = np.in1d(right_ind,index_overlap),np.in1d(left_ind,index_overlap)
        global_ind = np.append(np.append(left_ind[mask_left],right_ind[~mask_right]),left_ind[~mask_left])

        # Some geometric properties
        rvecsize = np.sqrt(self.faces.rvectors[:,0]**2 + self.faces.rvectors[:,1]**2)
        triangvols_left = 0.5 * self.faces.areas * 0.5 * rvecsize # volume of a triangle
        triangvols_right = triangvols_left
        triangcentroid_left = 1.0 /3 * (vor.points[vor.ridge_points[:,0]] +  vor.vertices[np.asarray(vor.ridge_vertices)[:,0]] + vor.vertices[np.asarray(vor.ridge_vertices)[:,1]])
        triangcentroid_right = 1.0 /3 * (vor.points[vor.ridge_points[:,1]] +  vor.vertices[np.asarray(vor.ridge_vertices)[:,0]] + vor.vertices[np.asarray(vor.ridge_vertices)[:,1]])

        # Volumes
        cellvol_left, _ = sum_by_group(triangvols_left,vor.ridge_points[:,0])
        cellvol_right, _ = sum_by_group(triangvols_right,vor.ridge_points[:,1])
        cellvols = combine_left_and_right(cellvol_left,cellvol_right,mask_left,mask_right)

        # Centroids
        centroidx_left, _ = sum_by_group(triangcentroid_left[:,0] * triangvols_left,vor.ridge_points[:,0])
        centroidx_right,_ = sum_by_group(triangcentroid_right[:,0] * triangvols_right,vor.ridge_points[:,1])
        centroidy_left, _ = sum_by_group(triangcentroid_left[:,1] * triangvols_left,vor.ridge_points[:,0])
        centroidy_right, _ = sum_by_group(triangcentroid_right[:,1] * triangvols_right,vor.ridge_points[:,1])
        centroidsx = combine_left_and_right(centroidx_left,centroidx_right,mask_left,mask_right,kind='sum')
        centroidsy = combine_left_and_right(centroidy_left,centroidy_right,mask_left,mask_right,kind='sum')
        centroids = np.array([centroidsx,centroidsy]).T
        centroids[cellvols > 0,:] = centroids[cellvols > 0,:]/cellvols[cellvols > 0,None]
        centroids[cellvols <= 0,:] = np.nan
        
        # Areas
        cellarea_left, _ = sum_by_group(self.faces.areas,vor.ridge_points[:,0])
        cellarea_right, _ = sum_by_group(self.faces.areas,vor.ridge_points[:,1])
        cellareas = combine_left_and_right(cellarea_left,cellarea_right,mask_left,mask_right,kind='sum')

        return centroids[global_ind.argsort()], cellvols[global_ind.argsort()], cellareas[global_ind.argsort()]




def compute_voronoi_gradients(VoronoiMesh,quant):

       rvecsize = np.sqrt(VoronoiMesh.faces.rvectors[:,0]**2 + VoronoiMesh.faces.rvectors[:,1]**2)
       cvectors = VoronoiMesh.faces.centroids - VoronoiMesh.faces.midpoints
       
       delta_quant =  quant[VoronoiMesh.faces.right_point_ind] -  quant[VoronoiMesh.faces.left_point_ind]
       try:
              gradient_sum_x_left = VoronoiMesh.faces.areas * (delta_quant * cvectors[:,0] - 0.5 * quant[VoronoiMesh.faces.right_point_ind] * VoronoiMesh.faces.rvectors[:,0]) \
                                    / rvecsize /  VoronoiMesh.cells.vols[VoronoiMesh.faces.left_point_ind]
       except RuntimeWarning:
              gradient_sum_x_left = 0

       try:
              gradient_sum_y_left = VoronoiMesh.faces.areas * (delta_quant * cvectors[:,1] - 0.5 * quant[VoronoiMesh.faces.right_point_ind] * VoronoiMesh.faces.rvectors[:,1]) / rvecsize /  VoronoiMesh.cells.vols[VoronoiMesh.faces.left_point_ind]
       except RuntimeWarning:
              gradient_sum_y_left = 0
       try:
              gradient_sum_x_right = VoronoiMesh.faces.areas * (-delta_quant * cvectors[:,0] + 0.5 * quant[VoronoiMesh.faces.left_point_ind] * VoronoiMesh.faces.rvectors[:,0]) / rvecsize /  VoronoiMesh.cells.vols[VoronoiMesh.faces.right_point_ind]
       except RuntimeWarning:
              gradient_sum_x_right = 0
       try:
              gradient_sum_y_right = VoronoiMesh.faces.areas * (-delta_quant * cvectors[:,1] + 0.5 * quant[VoronoiMesh.faces.left_point_ind] * VoronoiMesh.faces.rvectors[:,1]) / rvecsize /  VoronoiMesh.cells.vols[VoronoiMesh.faces.right_point_ind]
       except RuntimeWarning:
              gradient_sum_y_right = 0  
              
       # Sum the gradient contributions by groups
       gradx_left, left_ind = sum_by_group(gradient_sum_x_left,VoronoiMesh.faces.left_point_ind)
       gradx_right, right_ind = sum_by_group(gradient_sum_x_right,VoronoiMesh.faces.right_point_ind)
       grady_left, _ = sum_by_group(gradient_sum_y_left,VoronoiMesh.faces.left_point_ind)
       grady_right, _ = sum_by_group(gradient_sum_y_right,VoronoiMesh.faces.right_point_ind)
       
       # indices
       index_overlap  = np.intersect1d(left_ind,right_ind)
       mask_right, mask_left = np.in1d(right_ind,index_overlap),np.in1d(left_ind,index_overlap)
       global_ind =  np.append(np.append(left_ind[mask_left],right_ind[~mask_right]),left_ind[~mask_left])
       
       
       # Combine left-and-right information
       gradx_global = combine_left_and_right(gradx_left,gradx_right,mask_left,mask_right,kind='sum')
       grady_global = combine_left_and_right(grady_left,grady_right,mask_left,mask_right,kind='sum')
       
       gradientx = gradx_global[global_ind.argsort()]
       gradienty = grady_global[global_ind.argsort()]
       
       return gradientx, gradienty


def compute_voronoi_limiter(VoronoiMesh,quant,gradquant):

    
    #
    max_left, min_left, left_ind = group_by_group(quant[VoronoiMesh.faces.right_point_ind],VoronoiMesh.faces.left_point_ind)
    max_right, min_right, right_ind = group_by_group(quant[VoronoiMesh.faces.left_point_ind],VoronoiMesh.faces.right_point_ind)

    # indices
    index_overlap  = np.intersect1d(left_ind,right_ind)
    mask_right, mask_left = np.in1d(right_ind,index_overlap),np.in1d(left_ind,index_overlap)
    global_ind =  np.append(np.append(left_ind[mask_left],right_ind[~mask_right]),left_ind[~mask_left])

    
    max_global = combine_left_and_right(max_left,max_right,mask_left,mask_right,kind='max')
    min_global = combine_left_and_right(min_left,min_right,mask_left,mask_right,kind='min')
    max_global = np.maximum(max_global[global_ind.argsort()],quant)
    min_global = np.minimum(min_global[global_ind.argsort()],quant)
    
    # Linear local change (per Voronoi face) given the current gradient estimates
    delta_left =  gradquant[:,0][VoronoiMesh.faces.left_point_ind] * (VoronoiMesh.faces.centroids[:,0] - VoronoiMesh.cells.centroids[VoronoiMesh.faces.left_point_ind][:,0]) + \
                          gradquant[:,1][VoronoiMesh.faces.left_point_ind] * (VoronoiMesh.faces.centroids[:,1] - VoronoiMesh.cells.centroids[VoronoiMesh.faces.left_point_ind][:,1])
    delta_right = gradquant[:,0][VoronoiMesh.faces.right_point_ind] * (VoronoiMesh.faces.centroids[:,0] - VoronoiMesh.cells.centroids[VoronoiMesh.faces.right_point_ind][:,0]) + \
                          gradquant[:,1][VoronoiMesh.faces.right_point_ind] * (VoronoiMesh.faces.centroids[:,1] - VoronoiMesh.cells.centroids[VoronoiMesh.faces.right_point_ind][:,1])

    limiter_left,limiter_right = np.ones(delta_left.shape[0]),np.ones(delta_right.shape[0])
    ind_positive,ind_negative =  (delta_left > 0) & (VoronoiMesh.faces.areas > 1.0e-5 * VoronoiMesh.cells.areas[VoronoiMesh.faces.left_point_ind]),\
                                 (delta_left < 0) & (VoronoiMesh.faces.areas > 1.0e-5 * VoronoiMesh.cells.areas[VoronoiMesh.faces.right_point_ind])

    quant_maximum_left, quant_maximum_right = max_global[VoronoiMesh.faces.left_point_ind],\
                                              max_global[VoronoiMesh.faces.right_point_ind]
    quant_minimum_left, quant_minimum_right = min_global[VoronoiMesh.faces.left_point_ind],\
                                              min_global[VoronoiMesh.faces.right_point_ind]
    
    limiter_left[ind_positive] = (quant_maximum_left[ind_positive] - quant[VoronoiMesh.faces.left_point_ind][ind_positive])\
                                 / delta_left[ind_positive]
    limiter_left[ind_negative] = (quant_minimum_left[ind_negative] - quant[VoronoiMesh.faces.left_point_ind][ind_negative])\
                                 / delta_left[ind_negative]

    ind_positive,ind_negative =  (delta_right > 0) & (VoronoiMesh.faces.areas > 1.0e-5 * VoronoiMesh.cells.areas[VoronoiMesh.faces.right_point_ind]),\
                                 (delta_right < 0) & (VoronoiMesh.faces.areas > 1.0e-5 * VoronoiMesh.cells.areas[VoronoiMesh.faces.left_point_ind])
    limiter_right[ind_positive] = (quant_maximum_right[ind_positive] - quant[VoronoiMesh.faces.right_point_ind][ind_positive])\
                                 / delta_right[ind_positive]
    limiter_right[ind_negative] = (quant_minimum_right[ind_negative] - quant[VoronoiMesh.faces.right_point_ind][ind_negative])\
                                 / delta_right[ind_negative]

    # Need to find just one limiter per cell
    _, min_limiter_left, _ = group_by_group(limiter_left,VoronoiMesh.faces.left_point_ind)
    _, min_limiter_right, _ = group_by_group(limiter_right,VoronoiMesh.faces.right_point_ind)

    min_limiter_global = combine_left_and_right(min_limiter_left,min_limiter_right,mask_left,mask_right,kind='min')
    limiter = np.minimum(np.ones(min_limiter_global.shape[0]), min_limiter_global[global_ind.argsort()])

    return limiter
