#Damian Rumble, UoE
#25/09/2013
#run11.py

#This is a copy of 'run7.py' with files changes to plot test temp. maps instead of dust maps. 

#####################################################
#import maths, ploting and astrophysical packages

import numpy 
import aplpy
import matplotlib.pyplot as plt
import atpy
from pylab import *

####################################################
#plot background DUST map - subtitute in as appropriate 

single_ratio='output/map/single/single450850.fits'

####################################################
#plot process map - subtitute in as appropriate 

convolve_450='output/450/single_s450/s450convolve.fits'
convolve_850='output/850/single_s850/s850convolve.fits'

####################################################
#plot background TEMP map - subtitute in as appropriate 

singleT='output/map/single/single450850temp1_8snr5%50.fits'

####################################################
#SET map

s850 = aplpy.FITSFigure(convolve_450)

#####################################################################
#show contours

#show_contour(data, hdu=0, layer=None, levels=5, filled=False, cmap=None, colors=None, returnlevels=False, convention=None, dimensions=[0, 1], slices=[], smooth=None, kernel='gauss', overlap=False, **kwargs)

#####################################################
#map items

#show_label(cset,inline=True,fmt='%1.1f',fontsize=10)

s850.show_colorscale(cmap='gist_yarg',vmax=0.1,vmin=0.000006,stretch='log')
#s850.show_grayscale(invert=True)
s850.tick_labels.set_font(size='medium')
s850.add_colorbar()

#####################################################################
#show plot

#led = s850.show_contour.legend(loc='lower right' )
#plt.legend(('Model length', 'Data length', 'Total message length'),'upper center', shadow=True, fancybox=True)

plt.show()
