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

input450 = '/data/damian/maps/serpens/serpensmain/IR2/SerpensMain_20141223_s450_IR2extmask_s2_cal_JypixJH.fits'
convolve450 = '/data/damian/temp_MKR/output/beam/map/SerpensMain_20141223-auto/s450convolve.fits'
align450 = '/data/damian/temp_MKR/output/beam/process/SerpensMain_20141223-auto/gas450align.fits'
thresh450 = '/data/damian/temp_MKR/output/beam/process/SerpensMain_20141223-auto/s450alignTH.fits'

input850 = '/data/damian/maps/serpens/serpensmain/IR2/SerpensMain_20141223_s850_IR2extmask_s2_cal_JypixJH.fits'
convolve850 = '/data/damian/temp_MKR/output/beam/map/SerpensMain_20141223-auto/s850convolve.fits'
align850 = '/data/damian/temp_MKR/output/beam/process/SerpensMain_20141223-auto/s850collapse.fits'
thresh850 = '/data/damian/temp_MKR/output/beam/process/SerpensMain_20141223-auto/s850collapseTH.fits'

#Noise levels
sig850 = 0.0031
sig450 = 0.0218

#sub region position coords
pos1 = [277.47625,1.23358333333,0.13,0.13] #Main

#subregion figure sizes: ix,iy,dx,dy
plot8 = [0.15,0.05,0.4,0.22]
plot7 = [0.55,0.05,0.4,0.22]
plot6 = [0.15,0.27,0.4,0.22]
plot5 = [0.55,0.27,0.4,0.22]
plot4 = [0.15,0.49,0.4,0.22]
plot3 = [0.55,0.49,0.4,0.22]
plot2 = [0.15,0.71,0.4,0.22]
plot1 = [0.55,0.71,0.4,0.22]

### Deffine large map ###

#deffine 'big' fig
big = plt.figure(figsize=(8,16))

################
### 450 MAPs ###
################

#deffine 
Fig1 = aplpy.FITSFigure(input450,figure=big,subplot=plot2)
Fig2 = aplpy.FITSFigure(convolve450,figure=big,subplot=plot4)
Fig3 = aplpy.FITSFigure(align450,figure=big,subplot=plot6)
Fig4 = aplpy.FITSFigure(thresh450,figure=big,subplot=plot8)

###Show colour-scale
Fig1.show_grayscale(vmax=10.*sig450,vmin=0.0, stretch='sqrt',invert=True)
Fig1.add_colorbar()
Fig1.colorbar.show(ticks=None,axis_label_text='450$\mu m$ flux density (Jy/pixel)')
Fig2.show_grayscale(vmax=10.*sig450,vmin=0.0, stretch='sqrt',invert=True)
Fig2.add_colorbar()
Fig2.colorbar.hide()#show(location='top',ticks=None,axis_label_text=None)
Fig3.show_grayscale(vmax=10.*sig450,vmin=0.0, stretch='sqrt',invert=True)
Fig3.add_colorbar()
Fig3.colorbar.hide()#show(location='top',ticks=None,axis_label_text=None)
Fig4.show_grayscale(vmax=10.*sig450,vmin=0.0, stretch='sqrt',invert=True)
Fig4.add_colorbar()
Fig4.colorbar.hide()#show(location='top',ticks=None,axis_label_text=None)

#PLOT contours
#SCUBA-2
Fig1.show_contour(input450, levels=(sig450,1000), linewidth=3, colors=('k'))
Fig2.show_contour(input450, levels=(sig450,1000), linewidth=3, colors=('k'))
Fig3.show_contour(input450, levels=(sig450,1000), linewidth=3, colors=('k'))
Fig4.show_contour(input450, levels=(sig450,1000), linewidth=3, colors=('k'))

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig3.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig4.recenter(pos1[0],pos1[1],pos1[2],pos1[3])

#Label
Fig1.add_label(0.3, 0.9, r'450$\mathrm{\mu} m$ input', relative=True)
Fig2.add_label(0.3, 0.9, r'450$\mathrm{\mu} m$ convolve', relative=True)
Fig3.add_label(0.3, 0.9, r'450$\mathrm{\mu} m$ align', relative=True)
Fig4.add_label(0.3, 0.9, r'450$\mathrm{\mu} m$ mask', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')
Fig3.ticks.show()
Fig3.ticks.set_color('black')
Fig3.tick_labels.set_font(size='small')
Fig4.ticks.show()
Fig4.ticks.set_color('black')
Fig4.tick_labels.set_font(size='small')
Fig4.tick_labels.set_xformat('dd:mm:ss')

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()
Fig1.axis_labels.hide_x()
Fig2.tick_labels.hide_x()
Fig3.tick_labels.hide_x()
Fig1.tick_labels.hide_x()   
 
#################
#### 850 MAPS ###
#################

#deffine 
Fig1 = aplpy.FITSFigure(input850,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(convolve850,figure=big,subplot=plot3)
Fig3 = aplpy.FITSFigure(align850,figure=big,subplot=plot5)
Fig4 = aplpy.FITSFigure(thresh850,figure=big,subplot=plot7)

#Show colour-scale
Fig1.show_grayscale(vmax=10.*sig850,vmin=0.0, stretch='sqrt',invert=True)
Fig1.add_colorbar()
Fig1.colorbar.show(ticks=None,axis_label_text='850$\mu m$ flux density (Jy/pixel)')
Fig2.show_grayscale(vmax=10.*sig850,vmin=0.0, stretch='sqrt',invert=True)
Fig2.add_colorbar()
Fig2.colorbar.hide()#show(location='top',ticks=None,axis_label_text=None)
Fig3.show_grayscale(vmax=10.*sig850,vmin=0.0, stretch='sqrt',invert=True)
Fig3.add_colorbar()
Fig3.colorbar.hide()#show(location='top',ticks=None,axis_label_text=None)
Fig4.show_grayscale(vmax=10.*sig850,vmin=0.0, stretch='sqrt',invert=True)
Fig4.add_colorbar()
Fig4.colorbar.hide()#show(location='top',ticks=None,axis_label_text=None)

#PLOT contours
#SCUBA-2
Fig1.show_contour(input850, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig2.show_contour(input850, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig3.show_contour(input850, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig4.show_contour(input850, levels=(sig850,1000), linewidth=3, colors=('k'))

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig3.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig4.recenter(pos1[0],pos1[1],pos1[2],pos1[3])

#Label
Fig1.add_label(0.3, 0.9, r'850$\mathrm{\mu} m$ input', relative=True)
Fig2.add_label(0.3, 0.9, r'850$\mathrm{\mu} m$ convolve', relative=True)
Fig3.add_label(0.3, 0.9, r'850$\mathrm{\mu} m$ collapse', relative=True)
Fig4.add_label(0.3, 0.9, r'850$\mathrm{\mu} m$ mask', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')
Fig3.ticks.show()
Fig3.ticks.set_color('black')
Fig3.tick_labels.set_font(size='small')
Fig4.ticks.show()
Fig4.ticks.set_color('black')
Fig4.tick_labels.set_font(size='small')
Fig4.tick_labels.set_xformat('dd:mm:ss')

Fig1.axis_labels.hide_y()
Fig2.axis_labels.hide_y()
Fig3.axis_labels.hide_y()
Fig4.axis_labels.hide_y()
Fig1.tick_labels.hide_y()
Fig2.tick_labels.hide_y()
Fig3.tick_labels.hide_y()
Fig4.tick_labels.hide_y()

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()
Fig1.axis_labels.hide_x()
Fig2.tick_labels.hide_x()
Fig3.tick_labels.hide_x()
Fig1.tick_labels.hide_x() 



big.canvas.draw()

plt.show()
