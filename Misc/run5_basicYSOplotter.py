#Damian Rumble, UoE
#08/02/2013
#run5.py

#The script is designed to produce a map of each of the 4 regions of aquila, plot their masks in contours on top. Plot catalogues of YSO over them to check how accurate the masks are. 

#this version is for serpens south and East region

#Additionally this script uses atpy to load data from the GBS catalogue - sorts YSO by alpha classification and plots. 

#!!!!! ARCHIVED !!!!!!

#####################################################
#import maths, ploting and astrophysical packages

import numpy 
import aplpy
import matplotlib.pyplot as plt
import atpy
import c2dtable
import alpha_func as AF

#####################################################
#plot background maps

s850 = aplpy.FITSFigure('coadd_aquila_850_u20_cal556_th.fits')
s850.show_grayscale(invert='false')
s850.tick_labels.set_font(size='small')
s850.add_colorbar()

####################################################
#load GBS catalogue 

tbl = atpy.Table()
c2dtable.c2d_read('catalog-Serp-YSOc+head.tbl',tbl,3)
ra_t=tbl.data['RA']
dec_t=tbl.data['DEC']
alpha = tbl.data['alpha']
print alpha

###################################################
#Classification by Colour

CYSO = []
i = 0
while (i < len(alpha)):
#iterate through length of the array
    a = alpha[i]
    if a >= 0.3:
#Class I as yellow
       colour = 'yellow'
    elif 0.3 > a >= -0.3:
#Class FS as Orange
       colour = 'orange'
    elif -0.3 > a >= -1.6:
#Class II as Red
        colour = 'none'
    elif -1.6 > a:
#Class III as Magenta
       colour = 'none'
    else:
#else(anything not a number)
        colour = 'cyan'
    CYSO.append(colour)
    i = i + 1

YSOcounts = []
YSOcounts.append(len(alpha))
YSOcounts.append(CYSO.count('yellow'))
YSOcounts.append(CYSO.count('orange'))
YSOcounts.append(CYSO.count('red'))
YSOcounts.append(CYSO.count('magenta'))
YSOcounts.append(CYSO.count('cyan'))
print YSOcounts

####################################################
#Plot YSO from table - APPARENTLY THIS WORKS NOW

ra = 200.0+ra_t
dec = dec_t

s850.show_markers(ra,dec, layer='GBS',edgecolor='none' ,facecolor=CYSO,marker='o',s=50, alpha=0.75,label='GBS YSO (c1 yellow, cFS orange)')

##################################################
#Load All YSO data

'''
key:
Class 0 (only) - lime
Class 0/1 - yellow
Class 1 - Yellow
Class FS - Orange
Class II - Red (not fetured)
Class III - magenta (not featured)
'''

##Connelly +
SO = numpy.loadtxt('20130211c10_c_0_I_south_deg.txt')
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='c01', edgecolor='black',facecolor='yellow',marker='+',s=80, alpha=1, label='Connelly10 NIR cI YSO')

#Gutermuth ^ & x
#MWC297 0
data = numpy.loadtxt('20130207_G09_C0_alpha_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='2', edgecolor='black',facecolor='lime',marker='^',s=25, alpha=0.8)
#MWC297 1
data = numpy.loadtxt('20130207_G09_CI_alpha_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='4', edgecolor='yellow',facecolor='yellow',marker='x',s=25, alpha=0.8)
data = numpy.loadtxt('20130207_G09_CI_mwc297.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='6', edgecolor='black',facecolor='yellow',marker='^',s=25, alpha=0.8)
#serpens 0
data = numpy.loadtxt('20130206Gutermuth09_yso_I_0_serpens_deg.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='2', edgecolor='black',facecolor='yellow',marker='^',s=25, alpha=0.8)
#Serpens 1
data = numpy.loadtxt('20130207_G09_CI_alpha_serpens.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='5', edgecolor='yellow',facecolor='yellow',marker='x',s=25, alpha=0.8)

##Evans Square
# 1
data = numpy.loadtxt('20130131Evans09_yso_c1_serpens_deg.txt')
ra,dec = data[:,1], data[:,2]
s850.show_markers(ra,dec, layer='3', edgecolor='black',facecolor='yellow',marker='s',s=25, alpha=0.8)
# FS
data = numpy.loadtxt('20130131Evans09_yso_FS_serpens_deg.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='c2', edgecolor='red',facecolor='red',marker='o',s=30, alpha=0.6)

##Wintson Circle
#0/1
data = numpy.loadtxt('w07_yso_0_1_wcs_deg_only.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='7', edgecolor='black',facecolor='yellow',marker='o',s=25, alpha=0.8)
#FS
data = numpy.loadtxt('w07_yso_FS_wcs_deg_only.txt')
ra,dec = data[:,0], data[:,1]
s850.show_markers(ra,dec, layer='7', edgecolor='black',facecolor='yellow',marker='o',s=25, alpha=0.8

####################################################
#load and plot YSO from txt files




#data = numpy.loadtxt('')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='ctd', edgecolor='magenta',facecolor='magenta',marker='o',s=20, alpha=0.4)

#data = numpy.loadtxt('')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='c3', edgecolor='blue',facecolor='blue',marker='o',s=20, alpha=0.4)

#data = numpy.loadtxt('20130211w10_cIII_south_deg.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='c3', edgecolor='red',facecolor='red',marker='o',s=10, alpha=0.4)

YSO = numpy.loadtxt('20130211k10_xc_w40_deg.txt')
ra,dec = YSO[:,1], YSO[:,2]
s850.show_markers(ra,dec, layer='cx', edgecolor='black',facecolor='black',marker='x',s=25, alpha=1, label='Kuhn10 X-Ray sources')

#data = numpy.loadtxt('20130211k10_non_xc_w40_deg.txt')
#ra,dec = data[:,0], data[:,1]
#s850.show_markers(ra,dec, layer='cnx', edgecolor='black',facecolor='lime',marker='o',s=50, alpha=0.5)

####################################################
#Plot contours

s850.show_contour('aquila_s2_mask_cont.fits', colors='blue',label='mask')


####################################################
#Plot image and Legend.  Save image manually

led = s850._ax1.legend(loc='lower right' )

plt.show()

