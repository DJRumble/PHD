#Damian Rumble, UoE
#17/01/2013
#plotting_script.py

#This script is designed to plot YSO found in various research papaers on to existing SCUBA-2 data maps.

#!!!!!Archived!!!!!

#####################################################
#import maths, ploting and astrophysical packages

import numpy
import aplpy
import matplotlib.pyplot as plt

#####################################################
#plot background maps

s850 = aplpy.FITSFigure('coadd_serpensmain_s850_auto_mask_cal537_thresh_0_10.fits')
s850.show_grayscale(invert='true')
s850.tick_labels.set_font(size='small')

####################################################
#load and plot YSO

data = numpy.loadtxt('w07_yso_0_1_wcs_deg_only.txt')
ra,dec = data[:,0], data[:,1]

s850.show_markers(ra,dec, layer='marker_set_1', edgecolor='red',facecolor='red',marker='x',s=50, alpha=0.7)

####################################################
#Plot & Save image

plt.show()

s850.save('plot.png')
