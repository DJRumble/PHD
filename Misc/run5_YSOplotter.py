#Damian Rumble, UoE
#08/02/2013
#run5.py

#The script is designed to produce a map of each of the 4 regions of aquila, plot their masks in contours on top. Plot catalogues of YSO over them to check how accurate the masks are. 

#this version is for all regions but only YSO class 1 data. Colors are all Red, set is specified by icon.
#GBS - rhombus
#Connolly - plus
#Gutermuth - triangle
#Gutermuth alpha - x 
#Winston - circle
#Evans - square

#Additionally this script uses atpy to load data from the GBS catalogue.

#!!!!! ARCHIVED !!!!!!!

#####################################################
#import maths, ploting and astrophysical packages

import numpy
import aplpy
import matplotlib.pyplot as plt
import atpy
import c2dtable

#####################################################
#plot background maps

s850 = aplpy.FITSFigure('coadd_aquila_850_u20_cal556_th.fits')
s850.show_grayscale(invert='false')
s850.tick_labels.set_font(size='small')
s850.add_colorbar()

####################################################
#load GBS catalogue - DONE MANUALLY ON PYTHON TERMINAL

#tbl = atpy.Table()
#c2dtable.c2d_read('catalog-Serp-YSOc+head.tbl',tbl,3)
#ra_t=tbl.data['RA']
#dec_t=tbl.data['DEC']

####################################################
#Plot YSO from table - NOT WORKING ATM

#ra = 200.0+ra_t[:]
#dec = dec_t[:]
#s850.show_markers(ra,dec, layer='GBS', edgecolor='black',facecolor='lime',marker='^',s=50, alpha=1)

####################################################
#load and plot YSO from txt files

YSO = numpy.loadtxt('20130211c10_c_0_I_south_deg.txt')
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='1', edgecolor='red',facecolor='red',marker='+',s=25, alpha=0.8)

YSO = numpy.loadtxt('GBS_YSOcat_radec.txt')
ra,dec = YSO[:,0], YSO[:,1]
s850.show_markers(ra,dec, layer='GBS', edgecolor='black',facecolor='red',marker='d',s=25, alpha=0.8)

data = numpy.loadtxt('20130206Gutermuth09_yso_I_0_serpens_deg.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='2', edgecolor='black',facecolor='red',marker='^',s=25, alpha=0.8)

data = numpy.loadtxt('20130131Evans09_yso_c1_serpens_deg.txt')
ra,dec = data[:,1], data[:,2]
s850.show_markers(ra,dec, layer='3', edgecolor='black',facecolor='red',marker='s',s=25, alpha=0.8)

data = numpy.loadtxt('20130207_G09_CI_alpha_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='4', edgecolor='red',facecolor='red',marker='x',s=25, alpha=0.8)

data = numpy.loadtxt('20130207_G09_CI_alpha_serpens.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='5', edgecolor='red',facecolor='red',marker='x',s=25, alpha=0.8)

data = numpy.loadtxt('20130207_G09_CI_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='6', edgecolor='red',facecolor='red',marker='^',s=25, alpha=0.8)

data = numpy.loadtxt('w07_yso_0_1_wcs_deg_only.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='7', edgecolor='black',facecolor='red',marker='o',s=25, alpha=0.8)

####################################################
#Plot contours

s850.show_contour('aquila_s2_mask_cont.fits', colors='blue')


####################################################
#Plot & Save image

plt.show()


s850.save('plot.png')
