#Damian Rumble, UoE
#24/04/2013
#run7.py

#This is a copy of 'final_maps.py' with files changes to plot YSOs on temp. maps instead of dust maps. 

#The script is designed to produce a map of each of the 4 regions of aquila, Plot catalogues of YSO.

#this version should be able to handle all regions and all data set from one common directory - run6

#Additionally this script uses atpy to load data from the GBS catalogue - sorts YSO by alpha classification and plots. 

#YSO are classed as:
#letters indicate the authour who produced it
#numbers indicate the classification of the data (1 Protostars, 2 PMS, 3 other)

#Protostars - Lime
##Class 0  
##Class I 
##Class FS 

#Pre-Main-Sequence Stars - red
##Class II 
##Cass TD
##Class III 

#####################################################
#import maths, ploting and astrophysical packages

import numpy 
import aplpy
import matplotlib.pyplot as plt
import atpy
import c2dtable

####################################################
#plot background DUST map - subtitute in as appropriate 


Dsouth='maps/serpens_south/coadd_aquila_850_u20_cal556_th.fits'
Dmwc='maps/MWC297/serpensmwc297_s850_ext_mask2_cal537_mos_th.fits'
Dmain='maps/serpensmain/serpensmain_s850_ext_mask_cal537_mos_v2_th.fits'
DE1='maps/serpens_E1/serpensE1_s850_ext_mask1_cal537_mos_th.fits'

####################################################
#plot background TEMP map - subtitute in as appropriate 


Tsouth='output/FITS/aquila_emtest1000_cal491_mos450850temp1_8snr3.fits'
Tmwc='output/FITS/serpensmwc297_ext_mask2_cal491_mos450850temp1_8snr3.fits'
Tmain='output/FITS/serpensmain_ext_mask_cal491_mos_v2450850temp1_8snr3.fits'
TE1='output/FITS/serpensE1_ext_mask1_cal491_mos450850temp1_8snr3.fits'

####################################################
#SET map

s850 = aplpy.FITSFigure(Tsouth)

#####################################################################
#show contours

s850.show_contour(Dsouth, colors=('blue'),levels=(0.1,0.3),linestyles=('dotted','solid'))


#####################################################
#map items

s850.show_grayscale(invert='true')
s850.tick_labels.set_font(size='small')
s850.add_colorbar()

####################################################
#load GBS catalogue 

tbl = atpy.Table()
c2dtable.c2d_read('catalog-Serp-YSOc+head.tbl',tbl,3)
ra_t=tbl.data['RA']
dec_t=tbl.data['DEC']
alpha= tbl.data['alpha'] 

###################################################
#Classification by Colour

'''
key:
Protostars - yellow
Class 0  
Class I 
Class FS 

Pre-Main-Sequence Stars - red
Class II 
Cass TD
Class III 
'''

CYSO = []
i = 0
while (i < len(alpha)):
#iterate through length of the array
    a = alpha[i]
    if a >= 0.3:
#Class I as yellow
       colour = 'lime'
    elif 0.3 > a >= -0.3:
#Class FS as Orange
       colour = 'lime'
    elif -0.3 > a >= -1.6:
#Class II as Red
        colour = 'red'
    elif -1.6 > a:
#Class III as Magenta
       colour = 'red'
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

s850.show_markers(ra,dec, layer='GBS',edgecolor=CYSO ,facecolor='none',marker='d',s=25, alpha=1)


####################################################
#Load various data catalogs

'''
letters indicate the authour who produced it
numbers indicate the classification of the data (1 Protostars, 2 PMS, 3 other)
'''

c1 = 'data/connelly10/20130211c10_c_0_I_south_deg.txt'
g2 = 'data/gutermuth09/20130411_G09_PMS_all.txt'
g1 = 'data/gutermuth09/20130411_G09_protostars_all.txt'
e2 = 'data/evans09/20130411_E09_PMS.txt'  #note these are columns 1,2
e1 = 'data/evans09/20130411_E09_protostars.txt' #note these are columns 1,2
w2 = 'data/winston07/20130411_W07_PMS.txt'
w1 = 'data/winston07/20130411_W07_protostars.txt'
k1 = 'data/kuhn10/20130419_k10_protostars.txt'#note these are columns 1,2
k2 = 'data/kuhn10/20130419_k10_PMS.txt'#note these are columns 1,2
d2 = 'data/damiani06/20130125_d06_mwc297_xraysources_deg.txt'
i3 = 'data/damiani06/20130208_IRO_serpens_E1.txt'
s3 = 'data/damiani06/20130208_combined_serpens_E1.txt'
x3 = 'data/stars.txt'

##################################################
#Load All YSO data

'''
key:

Displays Protostars (classes 0,I,FS) as Lime 
Displays PMS stars (classes II,D and III) as Red

Markers - 
Connelly +
gutermuth v 
gutermuth (GBS as above) d 
evans s
winston o
Kuhn x
damiani x
IRO d
h2o h
'''
#PMS (2)

YSO = numpy.loadtxt(k2)
ra = YSO[:,1]
dec = YSO[:,2]
s850.show_markers(ra,dec, layer='k2', edgecolor='red',facecolor='none',marker='x',s=50, alpha=1)

YSO = numpy.loadtxt(d2)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='d2', edgecolor='red',facecolor='none',marker='x',s=50, alpha=1)

YSO = numpy.loadtxt(g2)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='g2', edgecolor='red',facecolor='none',marker='v',s=50, alpha=1)

YSO = numpy.loadtxt(e2)
ra = YSO[:,1]
dec = YSO[:,2]
s850.show_markers(ra,dec, layer='e2', edgecolor='red',facecolor='none',marker='s',s=50, alpha=1)

YSO = numpy.loadtxt(w2)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='w2', edgecolor='red',facecolor='none',marker='o',s=50, alpha=1)

#Protostars (1)

YSO = numpy.loadtxt(k1)
ra = YSO[:,1]
dec = YSO[:,2]
s850.show_markers(ra,dec, layer='k1', edgecolor='lime',facecolor='none',marker='x',s=50, alpha=1)

YSO = numpy.loadtxt(c1)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='g1', edgecolor='lime',facecolor='none',marker='v',s=50, alpha=1)

YSO = numpy.loadtxt(e1)
ra = YSO[:,1]
dec = YSO[:,2]
s850.show_markers(ra,dec, layer='e1', edgecolor='lime',facecolor='none',marker='s',s=50, alpha=1)

YSO = numpy.loadtxt(w1)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='w1', edgecolor='lime',facecolor='none',marker='o',s=50, alpha=1)

YSO = numpy.loadtxt(c1)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='c1', edgecolor='lime',facecolor='none',marker='+',s=50, alpha=1)



#other (3)

YSO = numpy.loadtxt(i3)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='i3', edgecolor='magenta',facecolor='none',marker='d',s=50, alpha=1)

YSO = numpy.loadtxt(s3)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='s3', edgecolor='cyan',facecolor='none',marker='h',s=50, alpha=1)

YSO = numpy.loadtxt(x3)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='xa', edgecolor='cyan',facecolor='cyan',marker='v',s=80, alpha=1)

YSO = numpy.loadtxt(x3)
ra = YSO[:,0]
dec = YSO[:,1]
s850.show_markers(ra,dec, layer='xb', edgecolor='cyan',facecolor='cyan',marker='^',s=80, alpha=1)



#####################################################################
#show plot


plt.show()
