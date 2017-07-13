import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from my_colortable import my_colortable

s, gamma, hue, r  = 0.45, 0.8, 1.1, -1.46
mycmap,nlow,nhigh = my_colortable(start=s,rots=r,hue=hue,gamma=gamma)


class ImageData():
    def __init__(self,*args,**kwargs):
        self.data = kwargs.get("data")

        if (self.data is not None):
            self.pixelx = self.data.shape[0]
            self.pixely = self.data.shape[1]
        
    def read(self, imagefile):

        f = open(imagefile,"rb")
	self.pixelx = np.fromfile(f, dtype=np.int32, count=1)[0]
	self.pixely = np.fromfile(f, dtype=np.int32, count=1)[0]
	self.data = np.fromfile(f, dtype=np.float32, count=self.pixelx*self.pixely).reshape(self.pixelx,self.pixely)
	f.close()

    def replacenan(self,value):
        self.data[np.isnan(self.data)] = value
        return
    
    def replaceneg(self,value=None):
        if (value is None):
            value = self.data[self.data > 0].min()
        self.data[self.data < 0] = value
        return

    def interpolateneg(self):
        for index in np.argwhere(self.data < 0):
            self.data[index[0],index[1]] = 0.25 * (self.data[index[0]-1,index[1]-1] +
                                                   self.data[index[0]-1,index[1]+1] +
                                                   self.data[index[0]+1,index[1]-1] +
                                                   self.data[index[0]+1,index[1]+1])
        return

        
def plot_slice(ax,image,vmin=None,vmax=None, scale='log',normalize=1.0,
               cmap = mycmap,extent=None):

    imagedata = image.data/normalize

    
    if (vmin is None):
        vmin = np.sort(imagedata).min()
    if (vmax is None):
        vmax= np.sort(imagedata).max()


    if (scale == 'log'):
        imagedata = np.log10(imagedata)
        vmax = np.log10(vmax)
        vmin = np.log10(vmin)

        
    im=ax.imshow(imagedata.T,vmin=vmin, vmax=vmax,
		 origin='lower',
		 interpolation='nearest',
		 extent=extent, cmap=cmap)
    cmap.set_under(color='k')
    
    


    return ax
    
