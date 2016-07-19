#Damian Rumble, UoE
#17/01/2013
#run3.py

#This script is designed to plot YSO found in various research papaers on to existing SCUBA-2 data maps.
#files mentioned in this script can be found in directory 'run1' where this script should be stored

#!!!!!!! ARCHIVED !!!!!!!!!

#####################################################
#import maths, ploting and astrophysical packages

import numpy
import aplpy
import matplotlib.pyplot as plt

#####################################################
#plot background maps

s850 = aplpy.FITSFigure('serpensmwc297_s450_ext_mask2_cal491_mos_th.fits')
s850.show_colorscale()
s850.tick_labels.set_font(size='small')

####################################################
#load and plot YSO

data = numpy.loadtxt('20130131Evans09_yso_c1_serpens_deg.txt')
ra,dec = data[:,1], data[:,2]
s850.show_markers(ra,dec, layer='c1', edgecolor='yellow',facecolor='yellow',marker='o',s=20, alpha=1)

data = numpy.loadtxt('20130131Evans09_yso_c2_serpens_deg.txt')
ra,dec = data[:,1], data[:,2]
s850.show_markers(ra,dec, layer='c2', edgecolor='orange',facecolor='orange',marker='o',s=30, alpha=0.8)

data = numpy.loadtxt('20130131Evans09_yso_c3_serpens_deg.txt')
ra,dec = data[:,1], data[:,2]
s850.show_markers(ra,dec, layer='c3', edgecolor='red',facecolor='red',marker='o',s=40, alpha=0.6)

data = numpy.loadtxt('20130131Evans09_yso_FS_serpens_deg.txt')
ra,dec = data[:,1], data[:,2]
s850.show_markers(ra,dec, layer='cFS', edgecolor='ora',facecolor='orange',marker='o',s=50, alpha=0.4)

#data = numpy.loadtxt('w07_yso_3_wcs_deg_only.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='YSO_3', edgecolor='yellow',facecolor='yellow',marker='x',s=10, alpha=0.2)

####################################################
#Plot contours

s850.show_contour('SER_90asec_Av.fits', colors='blue')


####################################################
#Plot & Save image

plt.show()

s850.save('plot.png')
