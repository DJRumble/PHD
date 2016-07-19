#20160103

#This plot plots SCUBA-2 data at each stage of the Daul-beam method for Serpens Main.


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

########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

####################################################
#maps

ratioK = '/data/damian/temp_MKR/output/kernel/map/SerpensMain_20141223-auto/SerpensMain_20141223-auto_Sratio.fits'
ratioDB = '/data/damian/temp_MKR_v10/output//map/SerpensMain_20141223_IR2extmask_s2_cal_JypixJH/SerpensMain_20141223_IR2extmask_s2_cal_JypixJH450850.fits'

tempK = '/data/damian/temp_MKR/output/kernel/map/SerpensMain_20141223-auto/SerpensMain_20141223-autotemperature.fits'
tempDB = '/data/damian/temp_MKR_v10/output//map/SerpensMain_20141223_IR2extmask_s2_cal_JypixJH/SerpensMain_20141223_IR2extmask_s2_cal_JypixJH450850temp%5.fits'

#Noise levels
sig850 = 0.0031
sig450 = 0.0218

#sub region position coords
pos1 = [277.47625,1.23358333333,0.04,0.04] #Main

#subregion figure sizes: ix,iy,dx,dy
#plot8 = [0.15,0.05,0.4,0.22]
#plot7 = [0.55,0.05,0.4,0.22]
#plot6 = [0.15,0.27,0.4,0.22]
#plot5 = [0.55,0.27,0.4,0.22]
plot4 = [0.15,0.05,0.4,0.45]
plot3 = [0.55,0.05,0.4,0.45]
plot2 = [0.15,0.50,0.4,0.45]
plot1 = [0.55,0.50,0.4,0.45]

### Deffine large map ###

#deffine 'big' fig
big = plt.figure(figsize=(8,16))

################
### 450 MAPs ###
################

#deffine 
Fig1 = aplpy.FITSFigure(ratioK,figure=big,subplot=plot2)
Fig2 = aplpy.FITSFigure(ratioDB,figure=big,subplot=plot4)

#Show colour-scale
Fig1.show_colorscale(vmax=10.,vmin=4.0)
#Fig1.add_colorbar()
#Fig1.colorbar.show(location='top',ticks=None,axis_label_text='450$\mu m$ flux density (Jy/pixel)')
Fig2.show_colorscale(vmax=10.,vmin=4.0)
#Fig2.add_colorbar()
#Fig2.colorbar.show(location='top',ticks=None,axis_label_text=None)

#PLOT contours
#SCUBA-2
#Fig1.show_contour(input450, levels=(sig850,1000), linewidth=3, colors=('k'))
#Fig2.show_contour(input450, levels=(sig850,1000), linewidth=3, colors=('k'))

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos1[0],pos1[1],pos1[2],pos1[3])

#Label
Fig1.add_label(0.3, 0.9, 'Kernel-convolution', relative=True)
Fig2.add_label(0.2, 0.9, 'Two-component beam', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')
Fig2.tick_labels.set_xformat('hh:mm:ss')

Fig1.add_colorbar()
Fig1.colorbar.set_axis_label_text('SCUBA-2 450/850 flux ratio')
Fig1.colorbar.set_location('top')
#Fig2.add_colorbar()
#Fig2.colorbar.hide()

#Fig2.axis_labels.hide_x()
Fig1.axis_labels.hide_x()
#Fig2.tick_labels.hide_x()
Fig1.tick_labels.hide_x()   
 
#################
#### 850 MAPS ####
#################

#deffine 
Fig1 = aplpy.FITSFigure(tempK,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(tempDB,figure=big,subplot=plot3)

#Show colour-scale
Fig1.show_colorscale(vmax=40.,vmin=10.0, stretch='linear',cmap='CMRmap')
#Fig1.add_colorbar()
#Fig1.colorbar.show(ticks=None,axis_label_text='850$\mu m$ flux density (Jy/pixel)')
Fig2.show_colorscale(vmax=40.,vmin=10.0, stretch='linear',cmap='CMRmap')
#Fig2.add_colorbar()
#Fig2.colorbar.show(location='top',ticks=None,axis_label_text=None)

#PLOT contours
#SCUBA-2
#Fig1.show_contour(input850, levels=(sig850,1000), linewidth=3, colors=('k'))
#Fig2.show_contour(input850, levels=(sig850,1000), linewidth=3, colors=('k'))

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos1[0],pos1[1],pos1[2],pos1[3])

#Label
Fig1.add_label(0.3, 0.9, 'Kernel-convolution', relative=True)
Fig2.add_label(0.2, 0.9, 'Two-component beam', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')
Fig2.tick_labels.set_xformat('hh:mm:ss')

Fig1.axis_labels.hide_y()
Fig2.axis_labels.hide_y()
Fig1.tick_labels.hide_y()
Fig2.tick_labels.hide_y()

Fig1.add_colorbar()
Fig1.colorbar.set_axis_label_text('Temp. (K)')
Fig1.colorbar.set_location('top')
#Fig2.add_colorbar()
#Fig2.colorbar.hide()

#Fig2.axis_labels.hide_x()
Fig1.axis_labels.hide_x()
#Fig2.tick_labels.hide_x()
Fig1.tick_labels.hide_x() 

big.canvas.draw()

plt.show()
