#Damian Rumble, UoE
#24/04/2013
#run7.py

#This is a copy of 'final_maps.py' with files changes to plot YSOs on temp. maps instead of dust maps. 

#The script is designed to produce a map of each of the 4 regions of aquila, Plot catalogues of YSO.

#this version should be able to handle all regions and all data set from one common directory - run6

#Additionally this script uses atpy to load data from the GBS catalogue - sorts YSO by alpha classification and plots. 

#14/06/2013 - added contours for creating anaylsis maps.

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
from pylab import *

####################################################
#plot background DUST map - subtitute in as appropriate 


Dsouth='maps/serpens_south/coadd_aquila_850_u20_cal556_th.fits'
Dmwc_DJR='maps/MWC297/SerpensMWC297_20130918_232648_s850_DJR_cal.fits'
Dmain='maps/serpensmain/serpensmain_s850_ext_mask_cal537_mos_v2_th.fits'
DE1='maps/serpens_E1/serpensE1_s850_ext_mask1_cal537_mos_th.fits'

s450convolve = 'analysis/4b/s450convolve.fits'
s850convolve = 'analysis/4b/s850convolve.fits'

####################################################
#plot background TEMP map - subtitute in as appropriate 

Tsouth='output_6/map/aquila_skyloopfull1_jypix/aquila_skyloopfull1_jypix450850temp1_8snr5.fits'
Tmwc='output7/map/SerpensMWC297_20130425_IR1_JH/SerpensMWC297_20130425_IR1_JH450850temp1_8snr5%5.sdf'
Tmwc_DJR='output7/map/SerpensMWC297_20130918_DJR_33/SerpensMWC297_20130916_233448_DJR_cal450850temp1_8snr5%3.fits'
Tmwc_old='output/map/serpensmwc297_ext_mask2_cal491_mos/serpensmwc297_ext_mask2_cal491_mos450850temp1_8snr5.fits'
Tmain='output7/map/SerpensMain_20130425_IR1_JH/SerpensMain_20130425_IR1_JH450850temp1_8snr5%3.fits'
TE1='output_6/map/SerpensE1_20130425_IR1_JH/SerpensE1_20130425_IR1_JH450850temp1_8snr5.fits '

delTmain='analysis/4b/FsubPtemp.fits'

####################################################
#SET map

s850 = aplpy.FITSFigure(Dmwc_DJR)

#####################################################################
#show contours

mwc850='maps/MWC297/SerpensMWC297_20121221_s850_IR1_JH.fits'
mwc450='maps/MWC297/SerpensMWC297_20130425_s450_IR1_JH.fits'

Amain_pos='analysis/4b/FsubPtemp_posmask.fits'
Amain_ng='analysis/4b/FsubPtemp_negmask.fits'
Amain_rel='analysis/4b/SerpensMain_relative_100delT.fits'

pos_mask='analysis/4b/FsubPtemp_posmask.fits'
neg_mask='analysis/4b/FsubPtemp_negmask.fits'

s450convolve_msk='analysis/4b/s450convolve_mask.fits'
s850convolve_msk='analysis/4b/s850convolve_mask.fits'
s450convolve_PB_msk='analysis/4b/s450convolve_PB_mask.fits'
s850convolve_PB_msk='analysis/4b/s850convolve_PB_mask.fits'

temp_mask='analysis/4b/SerpensMain_20130617_FB_mask.fits'
temp_edge='analysis/4b/SerpensMain_relative_100delT_mask.fits'

Terr_mwc='output7/map/SerpensMWC297_20130425_IR1_JH/SerpensMWC297_20130425_IR1_JH450850temp1_8snr5%5_errorpercent.fits'
Terr_main='output7/map/SerpensMain_20130425_IR1_JH450850temp1_8snr5%3_errorpercent.fits'

ext_450='analysis/4b/serpensmain_s2+h_mask_450_cont.fits'
ext_850='analysis/4b/serpensmain_s2+h_mask_850_cont.fits'

#s850.show_contour(Amain_pos, colors=('blue'),levels=(10,15,20,25))
#s850.show_contour(Amain_neg, colors=('red'),levels=())
#s850.show_contour(Amain_rel, colors=('cyan','blue','purple','magenta'),levels=(2,4,8,16))
#,linestyles=('dotted')


#s850.show_contour(neg_mask, colors=('blue'),fill=('blue'), levels=(1))
#s850.show_contour(pos_mask, colors=('red'),fill=('red'), levels=(1))

#s850.show_contour(s450convolve_msk, colors=('blue'), levels=(1))
#s850.show_contour(s450convolve_PB_msk, colors=('blue'),levels=(1))

#s850.show_contour(s850convolve_msk, colors=('blue'),fill=('blue'), levels=(1))
#s850.show_contour(s850convolve_PB_msk, colors=('white'),fill=('red'), levels=(1))

#s850.show_contour(ext_450, colors=('red'), levels=(1))
#s850.show_contour(ext_850, colors=('red'), levels=(1))

#s850.show_contour(temp_mask, colors=('green'), levels=(1))
#s850.show_contour(temp_edge, colors=('yellow'), levels=(1))

#cset1 = s850.show_contour(mwc450, levels=(0.6, 0.9, 1.2, 1.5), linewidth=3, colors=('red'))
#cset2 = s850.show_contour(mwc850, levels=(0.1, 0.2, 0.3, 0.4), linewidth=3, colors=('blue'))
cset2 = s850.show_contour(Tmwc_old, levels=(5,10,15,35), linewidth=(1,2,3), colors=('red'),linestlye=(1,2,3),)


#show_contour(data, hdu=0, layer=None, levels=5, filled=False, cmap=None, colors=None, returnlevels=False, convention=None, dimensions=[0, 1], slices=[], smooth=None, kernel='gauss', overlap=False, **kwargs)

#####################################################
#map items

#show_label(cset,inline=True,fmt='%1.1f',fontsize=10)

s850.show_colorscale(cmap='gist_yarg',vmax=0.08,vmin=0.0,stretch='linear')
#s850.show_grayscale(invert=True)
s850.tick_labels.set_font(size='medium')
s850.add_colorbar()
plt.title('Original temp map')

####################################################
#load GBS catalogue 

tbl = atpy.Table()
c2dtable.c2d_read('catalog-Serp-YSOc+head.tbl',tbl,3)
ra_t=tbl.data['RA']
dec_t=tbl.data['DEC']
alpha= tbl.data['alpha'] 

###################################################
#Classification by Colour

#key:
#Protostars - yellow
#Class 0  
#Class I 
#Class FS 

#Pre-Main-Sequence Stars - red
#Class II 
#Cass TD
#Class III 

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


#letters indicate the authour who produced it
#numbers indicate the classification of the data (1 Protostars, 2 PMS, 3 other)


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

#
#key:

#Displays Protostars (classes 0,I,FS) as Lime 
#Displays PMS stars (classes II,D and III) as Red

#Markers - 

#Connelly +
#gutermuth v 
#gutermuth (GBS as above) d 
#evans s
#winston o
#Kuhn x
#damiani x
#IRO d
#h2o h

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

#led = s850.show_contour.legend(loc='lower right' )
#plt.legend(('Model length', 'Data length', 'Total message length'),'upper center', shadow=True, fancybox=True)


plt.show()
