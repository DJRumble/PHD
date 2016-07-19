#12/07/2016

#This script is designed to plot temp maps of Oph.

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

###########
#PLOT-flux#
###########

s850 = 's850/OphSco_Main_20150325_850_IR2_noco_HKJypix_col.fits'
temp = 'temperature_maps/OphSco_Main-autotemperature.fits'
proto = 'YSO/surface/OphSco_Main_proto_surface_r300_units.fits'

sig850 = 0.002

pos1 = [246.667083333,-24.44925,0.32,0.32] 

#Deffine OB star data
data1 = np.loadtxt('OBstars_list.txt')
RA_OB = data1[:,0]	
DEC_OB = data1[:,1]

data2 = np.loadtxt('OBstars.txt',dtype='string')
RA_OBx = map(float, data2[:,2])	
DEC_OBx = map(float, data2[:,3])

Fig1 = aplpy.FITSFigure(temp)
Fig1.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')

#PLOT contours
#SCUBA-2
Fig1.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))

#Plot YSOs
YSO.c2dGBS_v2(Fig1)

#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])

#Label
Fig1.add_label(0.65, 0.49, 'A', relative=True)
Fig1.add_label(0.22, 0.42, 'B', relative=True)
Fig1.add_label(0.45, 0.26, 'C', relative=True)
Fig1.add_label(0.4, 0.15, 'E', relative=True)
Fig1.add_label(0.25, 0.08, 'F', relative=True)

Fig1.add_label(0.5, 0.62, 'S1', relative=True)
Fig1.add_label(0.9, 0.5, 'HD 147889', relative=True)
Fig1.add_label(0.1, 0.95, 'Oph L1688', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='large')
Fig1.tick_labels.set_xformat('hh:mm:ss')

Fig1.add_colorbar()
Fig1.colorbar.set_axis_label_text('T$_{d}$ (K)') 
Fig1.add_scalebar(371 * u.arcsecond,label=('0.25 pc'),corner=('top right'))

plt.show()
