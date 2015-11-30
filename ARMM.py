#Damian Rumble, UoE
#31/07/2013
#V2.0
#09/01/2014 
#ARMM.py

#This script is my primary plotting tool for maps. V2.0 incoorperates modules for Markers and a number of setup options in the form of functions (which may eventually end up in seperate modules them selves)

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
#import montage_wrapper as montage

#My modules
import YSO

#1 sigma levels
'''
STDV_south_450 = 0.00847368513332   #Jy/pixel skyloop
STDV_south_850 = 0.00151656703904 #Jy/pixel skyloop
STDV_main_450 = 0.00837840813063  #Jy/pixel
STDV_main_850 = 0.00167943638908  #Jy/pixel
STDV_main_850noco = 0.00168595226123  #Jy/pixel
STDV_mwc297_450 = 0.0209323405711 #0.0168632057269 #Jy/pixel
STDV_mwc297_850 = 0.00261840219553 #0.00218242691874 #Jy/pixel
STDV_E1_450 = 0.0185238342182  #Jy/pixel
STDV_E1_850 = 0.00202344639093  #Jy/pixel
STDV_ophA_450 = 0.0131259511768  #Jy/pixel
STDV_ophA_850 = 0.0187261 #Jy/pixel
'''

####################################################
#Plotting functions - wavelength specific 

def A(fig,l):
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)


    #FW1 
    ra = 277.03773
    dec = -3.80334
    Fig.show_markers(ra,dec, layer='1', edgecolor='red', marker='+',s=100, alpha=1)

    #FW 4
    ra = 277.2795
    dec = -3.72182
    Fig.show_markers(ra,dec, layer='4', edgecolor='red', marker='+',s=100, alpha=1)

    #FW 7
    ra = 277.2727
    dec = -3.71267
    Fig.show_markers(ra,dec, layer='7', edgecolor='red', marker='+',s=100, alpha=1)



    #MWC 297
    ra = 276.914708333
    dec = -3.83116666667
    Fig.show_markers(ra,dec, layer='mwc297', edgecolor='red', marker='+',s=100, alpha=1)
  
    #show SMM5 peak
    ra = 276.91551653269602
    dec = -3.8311826194934908
    Fig.show_markers(ra,dec, layer='SMM5', edgecolor='lime', marker='+',s=100, alpha=1)

    if l == '450':

    #Deffine contours
        cset = Fig.show_contour(fig, levels=(0.0494,0.247), linewidth=3, colors=('black','yellow'))

    #Figure properties
        Fig.show_grayscale(vmax=0.31,vmin=0.0,invert=True, stretch='sqrt')
        Fig.add_label(0.07, 0.95, r'b) $450\mu m$',size='large', relative=True)

    elif l == '850':

    #Deffine contours
        cset = Fig.show_contour(fig, levels=(0.00655,0.03328), linewidth=3, colors=('black','yellow'))

    #Figure properties
        Fig.show_grayscale(vmax=0.039,vmin=0.0,invert=True, stretch='sqrt')
        Fig.add_label(0.07, 0.95, r'a) $850\mu m$',size='large', relative=True)

    else:
        cset = 1

    #Deffine Markers
    #YSO.oth(Fig)
    #YSO.GBS(Fig)
    #YSO.mwc297(Fig)
    #YSO.mwc297_oth(Fig)

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig.ticks.show()
    Fig.ticks.set_color('white')

    return

def B(fig,contourA,contourB,contourC):
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    #Deffine contours
    cset1 = Fig.show_contour(contourA, levels=(0.00655,0.03328), linewidth=3, colors=('red'))
    cset2 = Fig.show_contour(contourB, levels=(92.4,462), linewidth=3, colors=('blue'))
    cset3 = Fig.show_contour(contourC, levels=(92.4,462), linewidth=3, colors=('blue'))

    #Deffine Markers
    #YSO.oth(Fig)
    YSO.GBS(Fig)
    YSO.mwc297(Fig)
    #YSO.mwc297_oth(Fig)

    #Figure properties
    Fig.show_grayscale(vmax=35.0,vmin=5.0,invert=True, stretch='sqrt')
    #Fig.add_label(0.1, 0.95, r'$450\mu m$', relative=True)
    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.ticks.show()
    Fig.ticks.set_color('black')

    return

def C(fig):
    #plot temp map
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    #Deffine contours
    cset = Fig.show_contour(fig, levels=(12,20,30), linewidth=3, colors=('black'))

    #Figure properties
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    Fig.show_colorscale(vmax=35.0,vmin=8.0, stretch='sqrt')
    Fig.add_label(0.05, 0.95, r'c) $T_{d}$',size='large', relative=True)

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.recenter(277.074583333,-3.74686111111, height=0.25 ,width=0.65)

    ra = 276.914708333
    dec = -3.83116666667

    return

def D(s450,s850,temp):
    #tri panel map - should tile all maps at once. 
    big = plt.figure(figsize=(12, 17))

    #deffine image in APLPY
    Fig_850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.1,0.64,0.8,0.26])
    Fig_450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.37,0.8,0.26])
    Fig_temp = aplpy.FITSFigure(temp,figure=big,subplot=[0.1,0.1,0.8,0.26])

    #Deffine contours
    cset1 = Fig_850.show_contour(s850, levels=(0.011,0.03328), linewidth=3, colors=('black','yellow'))
    cset2 = Fig_450.show_contour(s450, levels=(0.082,0.247), linewidth=3, colors=('black','yellow'))
    cset3 = Fig_temp.show_contour(temp, levels=(12,20,30), linewidth=3, colors=('black'))

    #MWC 297
    ra = 276.914708333
    dec = -3.83116666667
    Fig_450.show_markers(ra,dec, layer='mwc297', edgecolor='red', marker='+',s=100, alpha=1)
    Fig_850.show_markers(ra,dec, layer='mwc297', edgecolor='red', marker='+',s=100, alpha=1)

    #Figure properties
    Fig_450.show_grayscale(vmax=0.25,vmin=0.0,invert=True, stretch='sqrt')
    Fig_450.add_label(0.07, 0.95, r'b) $450\mu m$',size='large', relative=True)
    Fig_850.show_grayscale(vmax=0.034,vmin=0.0,invert=True, stretch='sqrt')
    Fig_850.add_label(0.07, 0.95, r'a) $850\mu m$',size='large', relative=True)
    Fig_temp.show_colorscale(vmax=45.0,vmin=8.0, stretch='sqrt')
    Fig_temp.add_label(0.05, 0.95, r'c) $T_{d}$',size='large', relative=True)

    Fig_450.tick_labels.set_font(size='medium')
    Fig_450.add_colorbar()
    Fig_450.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig_850.tick_labels.set_font(size='medium')
    Fig_850.add_colorbar()
    Fig_850.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig_temp.tick_labels.set_font(size='medium')
    Fig_temp.add_colorbar()
    Fig_temp.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)

    Fig_450.hide_xaxis_label()
    Fig_850.hide_xaxis_label()
    Fig_450.hide_xtick_labels()
    Fig_850.hide_xtick_labels()

    Fig_850.hide_yaxis_label()
    Fig_temp.hide_yaxis_label()

    Fig_450.ticks.show()
    Fig_450.ticks.set_color('black')
    Fig_850.ticks.show()
    Fig_850.ticks.set_color('black')
    Fig_temp.ticks.show()
    Fig_temp.ticks.set_color('black')

    big.canvas.draw()

    return


def E(fig,s850):
    #FELLWALKER CLUMPS

    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    #MWC 297
    ra = 276.914708333
    dec = -3.83116666667
    #Fig.show_markers(ra,dec, layer='mwc297', edgecolor='red', marker='+',s=100, alpha=1)

    #Deffine contours
    cset = Fig.show_contour(s850, levels=(0.011,1000), linewidth=3, colors=('black','yellow'))

    #Figure properties
    Fig.show_colorscale(vmax=25,vmin=0.0, stretch='linear')

    Fig.tick_labels.set_font(size='medium')
    #Fig.add_colorbar()
    Fig.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig.ticks.show()
    Fig.ticks.set_color('black')

    return

def F(fig,mask):
    #YSO MAP WITH CONTOURS FOR DISTRBUTIONS

    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    #Deffine Markers
    #YSO.GBS(Fig)
    YSO.mwc297(Fig)
    #YSO.damiani(Fig)

    #Deffine contours
    cset1 = Fig.show_contour(mask, levels=(0.5,10), linewidth=4, colors=('black'))
    cset2 = Fig.show_contour(fig, levels=(0.00655,10), linewidth=2, colors=('black'))

    #Figure properties
    Fig.show_grayscale(vmax=0.039,vmin=0.0,invert=True, stretch='sqrt')
    Fig.add_label(0.07, 0.95, r'$850\mu m$',size='large', relative=True)

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()

    Fig.ticks.show()
    Fig.ticks.set_color('black')
 
    return

def G(fig,scuba):
    #Spitzer maps - with SCUBA contours
    
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig,north=True)

    #Deffine contours
    cset = Fig.show_contour(scuba, levels=(0.011,0.022,0.033), linewidth=4, colors=('black'))

    #MWC 297
    ra = 276.914708333
    dec = -3.83116666667
    Fig.show_markers(ra,dec, layer='mwc297', edgecolor='red', marker='x',s=500, alpha=1)


    #Figure properties
    Fig.show_grayscale(vmax=120,vmin=20.0,invert=False, stretch='log')
    Fig.add_label(0.07, 0.95, r'$24\mu m$',size='large', relative=True)

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()

    Fig.ticks.show()
    Fig.ticks.set_color('black')
 
    return

def H(fig,fig_old):
    #plot temp map

    big = plt.figure(figsize=(18, 12))

    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig,figure=big,subplot=[0.1,0.1,0.4,0.8])
    Fig_old = aplpy.FITSFigure(fig_old,figure=big,subplot=[0.5,0.1,0.4,0.8])


    #Deffine contours
    cset = Fig.show_contour(fig, levels=(11,21,33), linewidth=(6,4,2), colors=('magenta'))
    cset_old = Fig_old.show_contour(fig_old, levels=(11,21,33), linewidth=(6,4,2), colors=('magenta'))

    #Figure properties
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    Fig.show_colorscale(vmax=55.0,vmin=10.0, stretch='sqrt')
    Fig.add_label(0.1, 0.95, r'a) Pre',size='large', relative=True)

    Fig.tick_labels.set_font(size='medium')
    #Fig.add_colorbar()
    Fig.recenter(276.914708333,-3.82116666667, height=0.04 ,width=0.04)

    Fig_old.axis_labels.hide_y()
    Fig_old.tick_labels.hide_y()
    Fig_old.ticks.show()
    Fig_old.ticks.set_color('black')
    Fig_old.show_colorscale(vmax=55.0,vmin=10.0, stretch='sqrt')
    Fig_old.add_label(0.1, 0.95, r'b) Post',size='large', relative=True)

    Fig_old.tick_labels.set_font(size='medium')
    Fig_old.add_colorbar()
    Fig_old.recenter(276.914708333,-3.82116666667, height=0.04 ,width=0.04)

    big.canvas.draw()

    return





####################################################
#plot background DUST map - subtitute in as appropriate 

#450 files
south_s450 = 'serpens_south/aquila_s450_skyloopfull1_jypix.fits'
main_s450 = 'serpensmain/SerpensMain_20130425_s450_IR1_J_threshed_snr.fits'
mwc297_s450 = 'MWC297/IR1/SerpensMWC297_20140331_IR1_s450_freefree_3d.fits'    
E1_s450 = 'serpens_E1/SerpensE1_20130425_s450_IR1_JH.fits'

#850 files
south_s850 = 'serpens_south/aquila_s850_skyloopfull1_jypix.fits'
main_s850 = 'serpensmain/fits_files/SerpensMain_20130425_s850_IR1_JH.fits'
mwc297_s850 = 'MWC297/IR1/SerpensMWC297_20140331_IR1_s850_freefree_3d.fits'
mwc297_850_section = 'MWC297/IR1/section_collapse.fits'
#mwc297_850_section = 'SerpensMWC297_20140501_s850_IR1_PSfreefree+jet103_section.fits'
     
E1_s850 = 'serpens_E1/SerpensE1_20121214_s850_IR1_JH.fits'

#24uM SPITZER
mwc297_s24a = 'spitzer/SER_E_2_A_MIPS1_mosaic.fits'
mwc297_s24b = 'spitzer/SER_E_3_A_MIPS1_mosaic.fits'

#temp. files
mwc297_temp='MWC297/IR1/SerpensMWC297_20140331_IR1_freefree_3d450850temp%5.fits'
mwc297_temp_old='MWC297/IR1/SerpensMWC297_20130425_IR1_JH450850temp%5.fits'

#cump files
clump = 'MWC297/SerpensMWC297_IR1_collapse.fits'

#Mask files 
mask = 'MWC297/serpensmwc297_s2+h_mask_850_magic.fits'

####################################################
############        MAKE MAPS HERE     #############
####################################################

#A(mwc297_s850,'850')
#A(mwc297_s450,'450')

#B(mwc297_temp,mwc297_s850,mwc297_s24a,mwc297_s24b)

#C(mwc297_temp)


#D(mwc297_s450, mwc297_s850, mwc297_temp)

#E(clump,mwc297_s850)

F(mwc297_850_section,mask)

#G(mwc297_s24a,mwc297_s850)

#H(mwc297_temp_old,mwc297_temp)

plt.show()
