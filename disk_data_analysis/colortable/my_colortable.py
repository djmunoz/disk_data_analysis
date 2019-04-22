import pylab as m
import numpy as np


def my_colortable(start=0.0,rots=-1.5,hue=1.0,gamma=1.0):

    cdict = {
        'red'  :  ( ),
        'green':  ( ),
        'blue' :  ( ),
        }

    nlow = 0
    nhigh = 0
    nlev = 1024
    for i in range(0,nlev):
        fract = float(i)/float(nlev-1)
        angle = 2.0*np.pi*(start/3.0 +1.0 +rots*fract)
        fract = fract**gamma
        amp = hue*fract*(1-fract)/2.0

        red =   fract+amp*(-0.14861*np.cos(angle)+1.78277*np.sin(angle))
        green = fract+amp*(-0.29227*np.cos(angle)-0.90649*np.sin(angle))
        blue =  fract+amp*(+1.97294*np.cos(angle))


        if red <= 0:
            red = 0
            nlow = nlow+1
        if green <= 0:
            green = 0
            nlow = nlow +1
        if blue <= 0:
            blue = 0
            nlow = nlow +1     
        if red > 1:
            red = 1
            nhigh = nhigh +1
        if green > 1:
            green = 1
            nhigh = nhigh +1
        if blue > 1:
            blue = 1
            nhigh = nhigh +1
            
        cdict['red'] = cdict['red'] + ((fract,red,red),)
        cdict['green'] = cdict['green'] + ((fract,green,green),)
        cdict['blue'] = cdict['blue']  + ((fract,blue,blue),)
    #generate the colormap with 1024 interpolated values
    my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

    return my_cmap,nlow,nhigh
