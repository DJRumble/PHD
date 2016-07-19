#Damian Rumble, UoE
#06/02/2013
#run4.py

#The script is designed to produce a map of each of the 4 regions of aquila, plot their masks in contours on top. Plot catalogues of YSO over them to check how accurate the masks are. 

#!!!!!!!! ARCHIVED !!!!!!!!!

#####################################################
#import maths, ploting and astrophysical packages

import numpy
import aplpy
import matplotlib.pyplot as plt

#####################################################
#plot background maps

s850 = aplpy.FITSFigure('serpensmwc297_s850_ext_mask2_cal537_mos_th.fits')
s850.show_grayscale(invert='false')
s850.tick_labels.set_font(size='small')
s850.add_colorbar()

####################################################
#load and plot YSO

#data = numpy.loadtxt('20130207_G09_C0_alpha_mwc297.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='c0', edgecolor='yellow',facecolor='yellow',marker='o',s=50, alpha=1)

data = numpy.loadtxt('20130207_G09_CI_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='c1', edgecolor='orange',facecolor='orange',marker='o',s=50, alpha=1)

data = numpy.loadtxt('20130207_G09_CII_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='c2', edgecolor='red',facecolor='red',marker='o',s=40, alpha=0.8)

#data = numpy.loadtxt('20130207_G09_CIII_alpha_mwc297.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='c3', edgecolor='magenta',facecolor='magenta',marker='o',s=30, alpha=0.6)

data = numpy.loadtxt('20130125_d06_mwc297_xraysources_deg.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='cx', edgecolor='black',facecolor='black',marker='x',s=50, alpha=1)

#data = numpy.loadtxt('20130206Gutermuth09_yso_II_serpens_deg.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='ctd', edgecolor='blue',facecolor='blue',marker='o',s=20, alpha=0.4)

#data = numpy.loadtxt('20130207_G09_CII_alpha_serpens.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='c3', edgecolor='red',facecolor='red',marker='o',s=10, alpha=0.4)

####################################################
#Plot contours

s850.show_contour('serpensmwc297_s2+h2+ard_mask_cont.fits', colors='blue')


####################################################
#Plot & Save image

plt.show()

s850.save('plot.png')
