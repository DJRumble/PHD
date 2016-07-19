#####################################################
#import maths, ploting and astrophysical packages

import numpy 
import aplpy
import matplotlib.pyplot as plt
import atpy
import c2dtable
from pylab import *
from astropy import units as u
import matplotlib.pyplot as mpl
import time
import astropy.io.fits as pyfits
import os

####################################################
#maps


ratio = '/data/damian/temp_MKR/output/beam/map/SerpensMain_20141223-auto/SerpensMain_20141223-auto_Sratio.fits'

#sub region position coords
pos1 = [277.47625,1.23358333333,0.13,0.13] #Main

Fig1 = aplpy.FITSFigure(ratio)
Fig1.show_colorscale(vmax=10,vmin=4, stretch='linear',cmap='bone')
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig1.add_label(0.2, 0.9, 'SCUBA-2 flux ratio', relative=True)
Fig1.ticks.show()
Fig1.add_colorbar()
Fig1.colorbar.set_axis_label_text('SCUBA-2 450/850 flux ratio')

plt.show()
