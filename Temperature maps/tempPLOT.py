
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

import ndf2fits
import YSO

########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

####################################################
#maps

temp_Rumble = 'PerseusWest_20150317-autotemperature.fits'
temp_Rumble2 = 'PerseusWest_20150317-autotemperatureWCSALN_SMT.fits'
temp_Chen = 'tempMap_Her+850_clean.fits' #'tempMap_Her+850_betaErrLess20per_final.fits'
beta_Chen = 'betaMap_Her+850_betaErrLess20per_final.fits'

pos1 = [52.2841666667,31.3251666667,0.03,0.05] #NGC1333

plot1 = [0.1,0.1,0.4,0.9]
plot2 = [0.5,0.1,0.4,0.9]

big = plt.figure(figsize=(16, 7))

#deffine 
Fig1 = aplpy.FITSFigure(temp_Rumble,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(temp_Chen,figure=big,subplot=plot2)

#Show colour-scale
Fig1.show_colorscale(vmax=28.0,vmin=8.0, stretch='linear',cmap='CMRmap')
Fig2.show_colorscale(vmax=28.0,vmin=8.0, stretch='linear',cmap='CMRmap')

#Show contours
Fig1.show_contour(beta_Chen, levels=(1.6,2.0),colors=('b','k'))
Fig2.show_contour(beta_Chen, levels=(1.6,2.0),colors=('cyan','w'))

#Recenter
Fig1.recenter(pos1[0],pos1[1],0.125,0.3)
Fig2.recenter(pos1[0],pos1[1],0.125,0.3)

#Label
Fig1.add_label(0.15, 0.95, "SCUBA-2 ratio", relative=True,size='x-large')
Fig2.add_label(0.35, 0.95, "SPIRE, PACS, SCUBA-2 SED fitting", relative=True,size='x-large')

#Deffine OB star data
data1 = np.loadtxt('OBstars_list.txt')
RA_OB = data1[:,0]	
DEC_OB = data1[:,1]

#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)

#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)

Fig2.axis_labels.hide_y()
Fig2.tick_labels.hide_y()

Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='x-large')
Fig1.axis_labels.set_font(size='x-large')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='x-large')
Fig2.axis_labels.set_font(size='x-large')
#Colour bar - remove if necessary

Fig1.add_colorbar()
Fig1.colorbar.hide()
Fig2.add_colorbar()
Fig2.colorbar.set_axis_label_text('Dust temperature (K)')
Fig2.colorbar.set_font(size='x-large')
Fig2.colorbar.set_axis_label_font(size='x-large')

big.canvas.draw()

plt.show()
