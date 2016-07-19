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
import os
import matplotlib.pyplot as plt
import atpy
import c2dtable
from pylab import *
from astropy import units as u
import matplotlib.pyplot as mpl
import astropy.io.fits as pyfits
#import montage_wrapper as montage

#My modules
import YSO
import ndf2fits 

########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

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

def A(s450,s850,YSO,HII1,mask):
    #deffine image in APLPY
    big = plt.figure(figsize=(18, 16))

    S450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.1,0.49,0.8])
    S850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.5,0.1,0.49,0.8])

    S450.show_grayscale(vmax=0.0213,vmin=0., stretch='linear',invert='TRUE')
    S850.show_grayscale(vmax=0.0015,vmin=0., stretch='linear',invert='TRUE')

    S850.axis_labels.hide_y()
    S850.tick_labels.hide_y()
    S850.ticks.show()
    S450.ticks.show()
    S450.ticks.set_color('black')
    S850.ticks.set_color('black')

    #SS
    #S450.recenter(277.406666667,-1.90194444444, height=0.85 ,width=0.28)
    #S850.recenter(277.406666667,-1.90194444444, height=0.85 ,width=0.28)
    #W40
    #S450.recenter(277.8995,-2.09200, height=0.60 ,width=0.40)
    #S850.recenter(277.8995,-2.09200, height=0.60 ,width=0.40)
    #East
    S450.recenter(279.4275,-1.36108333333, height=0.75 ,width=0.75)
    S850.recenter(279.4275,-1.36108333333, height=0.75 ,width=0.75)

    #277.406666667,-1.90194444444 SS
    #277.8995,-2.09200 W40

    S450.add_colorbar()
    S450.colorbar.set_location('top')
    S450.colorbar.set_axis_label_text('450$\mathrm{\mu}$m (Jy/pixel)')
    S450.add_beam(9.8*u.arcsecond,9.8*u.arcsecond,0,color='black')
    S450.add_scalebar(825 * u.arcsecond,label=('1 pc'),corner=('bottom left'))
    #S450.scalebar.show()
    #S450.scalebar.set_label('1 pc')
    #S450.scalebar.set_corner('bottom left')

    S850.add_colorbar()
    S850.colorbar.set_location('top')
    S850.colorbar.set_axis_label_text('850$\mathrm{\mu}$m (Jy/pixel)')
    S850.add_beam(14.5*u.arcsecond,14.5*u.arcsecond,0,color='black')
    S850.add_scalebar(825 * u.arcsecond,label=('1 pc'),corner=('bottom left'))
    #S850.scalebar.show()
    #S850.scalebar.set_label('1 pc')
    #S850.scalebar.set_corner('bottom left')


    #W40
    cset450 = S450.show_contour(s450, levels=(0.0173,100), linewidth=3, colors=('black'))
    cset850 = S850.show_contour(s850, levels=(0.0025,100), linewidth=3, colors=('black'))
    csetCO = S850.show_contour(mask, levels=(0.5,2), linewidth=2, colors=('blue'))

    #cset5 = S450.show_contour(YSO,levels=(15,45,75),colors=('r'))
    #cset1 = S850.show_contour(HII1, levels=(0.01,0.05,0.075,0.1,0.2), linewidth=4, colors=('blue')) #Units: Jy/Beam, Resolution: 45''
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    S450.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
    S850.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
    
    return

def a(s450,s850,s450mask,s850mask,south_s450_var,south_s850_var):
    #Plot SCUBA-2 data - both wavelengths and both variance maps
    #deffine image in APLPY
    big = plt.figure(figsize=(10, 10))

    S450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.55,0.4,0.4])
    S850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.1,0.1,0.4,0.4])
    V450 = aplpy.FITSFigure(south_s450_var,figure=big,subplot=[0.55,0.55,0.4,0.4])
    V850 = aplpy.FITSFigure(south_s850_var,figure=big,subplot=[0.55,0.1,0.4,0.4])

    S450.show_grayscale(vmax=0.185,vmin=0., stretch='sqrt',invert='TRUE')
    S850.show_grayscale(vmax=0.027,vmin=0., stretch='sqrt',invert='TRUE')
    V450.show_grayscale(vmax=5E-05,vmin=0., stretch='sqrt',invert='TRUE')
    V850.show_grayscale(vmax=1E-6,vmin=0., stretch='sqrt',invert='TRUE')

    S850.ticks.show()
    S450.ticks.show()
    S450.ticks.set_color('black')
    S850.ticks.set_color('black')
    S450.tick_labels.set_font(size='medium')
    S450.tick_labels.set_style('colons')
    S850.tick_labels.set_font(size='medium')
    S850.tick_labels.set_style('colons')
    V850.ticks.show()
    V450.ticks.show()
    V450.ticks.set_color('black')
    V850.ticks.set_color('black')
    V450.tick_labels.set_font(size='medium')
    V450.tick_labels.set_style('colons')
    V850.tick_labels.set_font(size='medium')
    V850.tick_labels.set_style('colons')

    S450.add_colorbar()
    #S450.colorbar.set_location('top')
    S450.colorbar.set_axis_label_text('450$\mu m$ (Jy/pixel)')
    #S450.add_beam(9.8*u.arcsecond,9.8*u.arcsecond,0,color='black')

    S850.add_colorbar()
    #S850.colorbar.set_location('top')
    S850.colorbar.set_axis_label_text('850$\mu m$ (Jy/pixel)')
    #S850.add_beam(14.5*u.arcsecond,14.5*u.arcsecond,0,color='black')

    V450.add_colorbar()
    V450.colorbar.set_axis_label_text('Variance 450$\mu m$ (Jy/pixel)')
    V850.add_colorbar()
    V850.colorbar.set_axis_label_text('Variance 850$\mu m$ (Jy/pixel)')

    cset450 = S450.show_contour(s450mask, levels=(0.5,100), linewidth=3, colors=('blue'))
    cset850 = S850.show_contour(s850mask, levels=(0.5,100), linewidth=3, colors=('blue'))

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

def C(fig,H):
    #plot temp map
    bigfig = mpl.figure(figsize=(12, 16))

    #deffine image in APLPY
    Fig1 = aplpy.FITSFigure(fig,figure=bigfig)
    Fig2 = aplpy.FITSFigure(fig,figure=bigfig,subplot=[0.635,0.13,0.16,0.16])
    #Deffine contours
    cset = Fig1.show_contour(H, levels=(300,1200,4800,12000), linewidth=3, colors=('k'))
    cset = Fig2.show_contour(H, levels=(300,1200,4800,12000), linewidth=3, colors=('k'))

    #Figure properties
    Fig1.ticks.show()
    Fig1.ticks.set_color('black')
    Fig1.show_colorscale(vmax=45.0,vmin=10.0, stretch='linear')

    Fig2.set_tick_labels_font(size='x-small')
    Fig2.set_axis_labels_font(size='small')
    Fig2.recenter(277.8325,-2.11200, height=0.04 ,width=0.04)
    Fig2.show_colorscale(vmax=45.0,vmin=10.0, stretch='linear')
    Fig2.ticks.show()
    Fig2.ticks.set_color('black')
    Fig2.hide_yaxis_label()
    #Fig2.hide_ytick_labels()
    Fig2.hide_xaxis_label()
    #Fig2.hide_xtick_labels()
    
    Fig1.tick_labels.set_font(size='medium')
    Fig1.add_colorbar()
    Fig1.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    Fig1.colorbar.set_axis_label_text('Dust temperature (K)')

    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig1.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
    Fig2.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
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


def E(fig):
    #FELLWALKER CLUMPS
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)
    YSO.w40(Fig)
    #Figure properties
    Fig.show_colorscale(vmax=82,vmin=0.0, stretch='sqrt')
    Fig.tick_labels.set_font(size='medium')
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    return

def F(fig,mask):
    #YSO MAP WITH CONTOURS FOR DISTRBUTIONS
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    #Deffine Markers
    #YSO.GBS(Fig)
    #YSO.mwc297(Fig)
    #YSO.damiani(Fig)
    #YSO.w40(Fig)
    #YSO.oth(Fig)

    OBmarkers = np.loadtxt('/data/damian/run27/OBstars_list.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=250, alpha=1)
    #OBmarkers = np.loadtxt('/data/damian/data/aquila/VLA3objects.txt')
    #ra = OBmarkers[:,0]
    #dec = OBmarkers[:,1]
    #Fig.show_markers(ra,dec, edgecolor='k', facecolor='red', marker='*',s=100, alpha=1)
    #Deffine contours
    Fig.show_contour(mask,levels=(50,1000),colors=('k'))

    #Figure properties
    #Fig.recenter(277.9525,-2.12600, height=0.35 ,width=0.45)
    Fig.show_contour(fig, levels=(0.0018,100), linewidth=3, colors=('black'))
    Fig.show_grayscale(vmax=0.026,vmin=0., stretch='linear',invert='TRUE')
    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.colorbar.set_axis_label_text('850$\mathrm{\mu}$m (Jy/pixel)')
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

def I(fig,red,blue):
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)
    #Deffine contours
    cset_red = Fig.show_contour(red, levels=(5,10,15,20,25,30,35), linewidth=3, colors=('r'))
    cset_red2 = Fig.show_contour(red, levels=(30,35,100), linewidth=3, filled=('TRUE'), colors=('r'))
    cset_blue = Fig.show_contour(blue, levels=(5,10,15,20,25,30,35), linewidth=3, colors=('b'))
    cset_blue2 = Fig.show_contour(blue, levels=(30,35,100), linewidth=3, filled=('TRUE'), colors=('b'))
    cset = Fig.show_contour(fig, levels=(0.0027,0.0081,0.027,0.081), linewidth=3, colors=('black'))
    #Figure properties
    Fig.show_grayscale(vmax=0.027,vmin=0., stretch='linear',invert='TRUE')
    #Fig.add_label(0.07, 0.95, r'$850\mu m$',size='large', relative=True)
    #Deffine Markers
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=25, alpha=1)
    YSO.w40_radio(Fig)

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    #Fig.colorbar.set_location('top')
    Fig.colorbar.set_axis_label_text('850$\mu m$ (Jy/pixel)')
    Fig.recenter(277.86083,-2.0885,height=0.1 ,width=0.1)
    Fig.ticks.show()
    Fig.ticks.set_color('white')
    return

def J(fig,l):
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)
    
    if l == '70':
    #Figure properties
        Fig.show_grayscale(vmax=2000.0,vmin=0.0,invert=True, stretch='sqrt')
        Fig.add_label(0.07, 0.95, r'$70\mu m$',size='large', relative=True)
    elif l == '160':
    #Figure properties
        Fig.show_grayscale(vmax=1000.,vmin=0.0,invert=True)
        Fig.add_label(0.07, 0.95, r'$160\mu m$',size='large', relative=True)
    elif l == '250':
       #Figure properties
        Fig.show_grayscale(vmax=2500.,vmin=0.0,invert=True)
        Fig.add_label(0.07, 0.95, r'$250\mu m$',size='large', relative=True)
    elif l == '350':
    #Figure properties
        Fig.show_grayscale(vmax=750.,vmin=0.0,invert=True,stretch='sqrt')
        Fig.add_colorbar()
        Fig.colorbar.set_axis_label_text('350$\mu m$ (MJy/Sr)')
        #Fig.add_label(0.07, 0.95, r'$350\mu m$',size='large', relative=True)
    elif l == '500':
    #Figure properties
        Fig.show_grayscale(vmax=500.,vmin=0.0,invert=True)
        Fig.add_label(0.07, 0.95, r'$350\mu m$',size='large', relative=True)
    else:
        cset = 1
    #Deffine Markers
    #YSO.GBS(Fig)
    #YSO.w40(Fig)
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=200, alpha=1)
    #MWC 297
    ra = 276.914708333
    dec = -3.83116666667
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=200, alpha=1)

    Fig.axis_labels.set_font(size='x-large')
    Fig.tick_labels.set_font(size='large')
    Fig.colorbar.set_font(size='large')
    Fig.colorbar.set_axis_label_font(size='x-large')
    
    #Fig.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig.ticks.show()
    Fig.ticks.set_color('black')

    return

def K(before,after,freefree):
    big = plt.figure(figsize=(18, 12))
    #deffine image in APLPY
    Before = aplpy.FITSFigure(before,figure=big,subplot=[0.1,0.1,0.4,0.8])
    After = aplpy.FITSFigure(after,figure=big,subplot=[0.5,0.1,0.4,0.8])
    YSO.w40_radio(Before)
    YSO.w40_radio(After)
    YSO.w40(Before)
    YSO.w40(After)
    #Figure properties
    Before.show_grayscale(vmax=0.032,vmin=0.0,invert=True,stretch='sqrt')
    #After.show_grayscale(vmax=0.032,vmin=0.0,invert=True,stretch='sqrt')
    After.show_grayscale(vmax=0.022,vmin=0.001375, stretch='sqrt',invert='TRUE')
    #Deffine contours
    cset = Before.show_contour(freefree, levels=(0.001,0.005,0.01), linewidth=(6,4,2), colors=('r'))
    cset2 = After.show_contour(freefree, levels=(0.001,0.005,0.01), linewidth=(6,4,2), colors=('r'))
    #Before.recenter(277.8925,-2.08200, height=0.1 ,width=0.15)
    #After.recenter(277.8925,-2.08200, height=0.1 ,width=0.15)
    Before.recenter(277.8525,-2.10200, height=0.03 ,width=0.035)
    After.recenter(277.8525,-2.10200, height=0.03 ,width=0.035)

    Before.tick_labels.set_font(size='medium')
    After.tick_labels.set_font(size='medium')
    Before.add_label(0.1, 0.95, r'SCUBA-2 450um',size='large', relative=True)
    After.add_label(0.1, 0.95, r'VLA 3.6cm',size='large', relative=True)
    Before.ticks.show()
    Before.ticks.set_color('black')
    After.axis_labels.hide_y()
    After.tick_labels.hide_y()
    After.ticks.show()
    After.ticks.set_color('black')
    #After.add_colorbar()
    #After.colorbar.set_axis_label_text('850$\mu m$ (Jy/pixel)')
    big.canvas.draw()
    return

def L(fig,HII1):#,mask):
    #W40 with CO observstion areas
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)
    #Deffine Markers 
    #YSO.GBS(Fig)
    #YSO.oth(Fig)
    YSO.w40_radio(Fig)
    #YSO.w40(Fig)
    #Deffine contour
    cset1 = Fig.show_contour(HII1, levels=(0.01,0.05,0.1), linewidth=4, colors=('blue')) #Units: Jy/Beam, Resolution: 45''
    #cset2 = Fig.show_contour(HII2, levels=(0.0015,0.015,0.0175,0.02), linewidth=2, colors=('cyan')) #Units: Jy/Beam, Resolution: 45''
    cset3 = Fig.show_contour(fig, levels=(0.0027,0.0081,0.027), linewidth=2, colors=('black','black','magenta'))
    #Figure properties
    Fig.show_grayscale(vmax=0.054,vmin=0.0,invert=True)
    Fig.add_label(0.07, 0.95, r'$850\mu m$',size='large', relative=True)
    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    return

def M(fig,flux,HII1,HII2,H70):#,mask):
    #W40 21cm map with Herschel Contours
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(HII1)
    #Deffine Markers
    #YSO.GBS(Fig)
    #YSO.w40(Fig)
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='yellow', facecolor='yellow', marker='*',s=100, alpha=1)
    #Deffine contours
    #cset1 = Fig.show_contour(HII1, levels=(0.01,0.1), linewidth=4, colors=('blue')) #Units: Jy/Beam, Resolution: 45''
    #cset2 = Fig.show_contour(HII2, levels=(0.0015,0.015,0.0175,0.02), linewidth=2, colors=('black')) #Units: Jy/Beam, Resolution: 45''
    #cset3 = Fig.show_contour(COmask, levels=1, linewidth=2, colors=('black')) #Mask of VanDerWeil2012 CO
    cset4 = Fig.show_contour(H70,levels=(300,1200,4800,12000),colors=('r')) #Herschel 70um data.
    cset3 = Fig.show_contour(flux, levels=(0.0025,1), linewidth=2, colors=('b'))
    #Figure properties
    #Fig.show_colorscale(vmax=35.0,vmin=15, stretch='sqrt')
    #Fig.add_label(0.9, 0.95, r'$T_{d}\,(K)$',size='large', relative=True)
    Fig.show_grayscale(vmax=0.05,vmin=0., stretch='linear', invert='True')
    #Fig.add_label(0.9, 0.95, r'$S_{21\,cm}$ (Jy/pix)',size='large', relative=True)
    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.colorbar.set_axis_label_text('$S_{21\,cm}$ (Jy/beam)')
    Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    Fig.add_beam(45*u.arcsecond,45*u.arcsecond,0,color='black')
    return

def N(fig,SCUBA2):
    Fig = aplpy.FITSFigure(fig)
    YSO.w40(Fig)
    YSO.w40_radio(Fig)
    #3.6cm
    Fig.show_grayscale(vmax=0.022,vmin=0.001375, stretch='sqrt',invert='TRUE')
    #21cm
    #Fig.show_colorscale(vmax=0.025,vmin=0.0015, stretch='linear')
    cset3 = Fig.show_contour(SCUBA2, levels=(0.0027,0.0081,0.027), linewidth=2, colors=('black','black','magenta'))
    Fig.add_colorbar()
    Fig.colorbar.set_axis_label_text('3.6 cm (Jy/Beam)')
    Fig.ticks.show()
    return

def m(s850,s21,mask):
    big = plt.figure(figsize=(18, 12))
    #deffine image in APLPY
    S850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.1,0.1,0.4,0.8])
    S21 = aplpy.FITSFigure(s21,figure=big,subplot=[0.5,0.1,0.4,0.8])
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    S850.show_markers(ra,dec, layer='W40a', edgecolor='k', facecolor='y', marker='*',s=200, alpha=1)
    S21.show_markers(ra,dec, layer='W40b', edgecolor='k', facecolor='y', marker='*',s=200, alpha=1)

    S850.show_colorscale(vmax=0.5,vmin=0.001627, stretch='log')
    S21.show_colorscale(vmax=0.5,vmin=0.001627, stretch='log')
    #S850.show_colorscale(vmax=25,vmin=15, stretch='linear')
    #S21.show_colorscale(vmax=25,vmin=15, stretch='linear')
    #cset1 = S850.show_contour(mask, levels=(0.5,10), linewidth=2, colors=('white'))
    #cset2 = S21.show_contour(mask, levels=(0.5,10), linewidth=2, colors=('white'))
    S21.axis_labels.hide_y()
    S21.tick_labels.hide_y()
    S21.ticks.show()
    S850.ticks.show()
    S850.recenter(277.8605,-2.1100, height=0.2 ,width=0.2)
    S21.recenter(277.8605,-2.1100, height=0.2 ,width=0.2)
    S850.add_label(0.1, 0.05, r'SCUBA-2',size='large', relative=True, color='white')
    S21.add_label(0.1, 0.05, r'VLA 21cm',size='large', relative=True, color='white')
    #S21.add_colorbar()
    #S21.colorbar.set_axis_label_text('850$\mu m$ flux density (Jy/15as pixel)')
    return

def O(fig,YSO): #Plot map of spectral index

    Fig = aplpy.FITSFigure(fig)
    Fig.show_colorscale(vmax=4.0,vmin=1.5, stretch='linear',cmap='cubehelix') #alpha
    #Fig.show_colorscale(vmax=6.5E+22,vmin=0.0, stretch='linear') #CD
    Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    Fig.add_colorbar()
    Fig.ticks.show()
    #Markers
    markers = np.loadtxt('/data/damian/run26/SMM/S2-FF-CO/SMM-cores2.tab')
    ra = markers[:,1]
    dec = markers[:,2]
    #Fig.show_markers(ra,dec, edgecolor='k', facecolor='k', marker='*',s=200, alpha=0.25)
    #Fig.show_markers(ra,dec, edgecolor='k', marker='*',s=200, alpha=1)

    Fig.colorbar.set_axis_label_text('Spectal index') #alpha
    #Fig.colorbar.set_axis_label_text('Column Density (10$^{22}$ H$_{2}$ cm$^{-2}$)') #CD
    #cset5 = Fig.show_contour(YSO,levels=(20,60,110,160,210),colors=('k'))
    Fig.axis_labels.set_font(size='x-large')
    Fig.tick_labels.set_font(size='x-large')
    Fig.colorbar.set_font(size='x-large')
    Fig.colorbar.set_axis_label_font(size='x-large')
    return

def P(s450,s850,s450ff,s850ff):

    big = plt.figure(figsize=(16, 12))

    #deffine image in APLPY
    S450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.1,0.4,0.8])
    S850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.5,0.1,0.4,0.8])

    #Figure properties
    S450.show_grayscale(vmax=0.185,vmin=0., stretch='linear',invert='TRUE')
    S850.show_grayscale(vmax=0.027,vmin=0., stretch='linear',invert='TRUE')
    
    #Deffine contours
    cset = S450.show_contour(s450ff, levels=(0.0104,0.0173,0.0519,100), colors=('r','k'),filled='TRUE')
    cset2 = S850.show_contour(s850ff, levels=(0.0015,0.0025,0.0075,100), colors=('r','k'),filled='TRUE')

    cset450 = S450.show_contour(s450, levels=(0.0104,0.0173,0.0519,0.1110), colors=('k'))
    cset850 = S850.show_contour(s850, levels=(0.0015,0.0025,0.0075,0.0161), colors=('k'))

    #recenter
    S450.recenter(277.8525,-2.10200, height=0.025 ,width=0.05)
    S850.recenter(277.8525,-2.10200, height=0.025 ,width=0.05)

    YSO.w40_radio(S450)
    #YSO.w40_radio(S850)

    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    S850.show_markers(ra,dec, layer='W40a', edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)
    S450.show_markers(ra,dec, layer='W40b', edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)

    ra = 277.849937
    dec = -2.091536
    S450.show_markers(ra,dec, edgecolor='k', facecolor='blue', marker='*',s=100, alpha=1)
    S850.show_markers(ra,dec, edgecolor='k', facecolor='blue', marker='*',s=100, alpha=1)



    #850 set up
    S850.axis_labels.hide_y()
    S850.tick_labels.hide_y()
    S850.ticks.show()
    S850.ticks.set_color('black')
    S850.add_colorbar()
    S850.colorbar.set_font(size='small')
    S850.colorbar.set_location('top')
    S850.colorbar.set_axis_label_text('850$\mu m$ (Jy/pixel)')
    S850.add_beam(14.5*u.arcsecond,14.5*u.arcsecond,0,color='black')
    S850.tick_labels.set_font(size='small')

    #450 set up
    S450.add_colorbar()
    S450.colorbar.set_location('top')
    S450.colorbar.set_axis_label_text('450$\mu m$ (Jy/pixel)')
    S450.colorbar.set_font(size='small')
    S450.add_beam(9.8*u.arcsecond,9.8*u.arcsecond,0,color='black')
    S450.ticks.show()
    S450.ticks.set_color('black')
    S450.tick_labels.set_font(size='small')
 
    big.canvas.draw()



    return

def Q(fig,SCUBA2):
    #VLA 3.6cm map
    Fig = aplpy.FITSFigure(fig)

    YSO.w40_radio(Fig)
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=200, alpha=1)
 
    cset850 = Fig.show_contour(SCUBA2, levels=(0.0025,0.0075,0.025),colors='k')

    Fig.show_grayscale(vmax=0.022,vmin=0.001375, stretch='sqrt',invert='TRUE')
    Fig.add_colorbar()
    Fig.recenter(277.8525,-2.09200, height=0.06 ,width=0.07)
    Fig.colorbar.set_axis_label_text('$S_{3.6\,cm}$ (Jy/beam)')
    Fig.add_beam(9.97*u.arcsecond,9.97*u.arcsecond,24.42*u.degree,color='black')
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    return

def R(fig,surfaceA,surfaceB,H,mask):
    #produce surface map of YSO distributions
    Fig = aplpy.FITSFigure(fig)
    
    levels = (15,45,75,100,150)

    cset_H = Fig.show_contour(H, levels=(200,10000),colors='blue')

    #cset1 = Fig.show_contour(surfaceA, levels=levels,colors='green')
    #cset2 = Fig.show_contour(surfaceB, levels=levels,colors='k')

    cset_mask = Fig.show_contour(mask, levels=(0.5,10),colors='k')

    #YSO.w40(Fig)
    YSO.c2dGBS(Fig)

    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=250, alpha=1)

    ra = 277.299583333
    dec = -2.06388888889
    Fig.show_markers(ra,dec, layer='W40a', edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)

    Fig.show_grayscale(vmax=0.027,vmin=0.0015, stretch='sqrt',invert='TRUE')
    #Fig.recenter(277.8995,-2.09200, height=0.6 ,width=0.6)
    Fig.recenter(277.4406666667,-1.90194444444, height=0.85 ,width=0.5)
    Fig.add_colorbar()
    #Fig.colorbar.set_location('top')
    Fig.colorbar.set_axis_label_text('850$\mu m$ (Jy/pixel)')
    #Fig.add_beam(9.8*u.arcsecond,9.8*u.arcsecond,0,color='black')

def S(fig,H70,HII1,HII2):#,mask):
    #W40 temp map with Herschel Contours

    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    #Deffine Markers
    #YSO.GBS(Fig)
    #YSO.w40(Fig)
    YSO.w40_radio(Fig)

    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)	

    #Deffine contours
    cset1 = Fig.show_contour(HII1, levels=(0.05,0.075,0.1,0.2), linewidth=4, colors=('blue')) #Units: Jy/Beam, Resolution: 45''
    #cset2 = Fig.show_contour(HII2, levels=(0.001,0.005,0.0075,0.01), linewidth=2, colors=('magenta')) #Units: Jy/Beam, Resolution: 45''
    #cset3 = Fig.show_contour(COmask, levels=1, linewidth=2, colors=('black')) #Mask of VanDerWeil2012 CO
    #cset4 = Fig.show_contour(H70,levels=(300,1200,4800,12000),colors=('k')) #Herschel 70um data.
    #cset3 = Fig.show_contour(flux, levels=(0.0027,0.0081,0.027), linewidth=2, colors=('black','black','magenta'))
    #cset5 = Fig.show_contour(YSO,levels=(15,45,75),colors=('magenta'))

    #Figure properties
    Fig.show_colorscale(vmax=30.0,vmin=10, stretch='sqrt')

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    Fig.colorbar.set_axis_label_text('$T_{d}\,(K)$')
    Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    return

def T(fig,H70,s21):
    #poster image - 850 + H70/21cm contours + CO

    #deffine image in APLPY
    Fig = aplpy.FITSFigure(fig)

    YSO.w40_radio(Fig)

    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)

    #cset4 = Fig.show_contour(H70,levels=(300,1200,4800,12000),colors=('blue')) #Herschel 70um data.
    cset1 = Fig.show_contour(s21, levels=(0.01,0.03,0.05,0.1), linewidth=4, colors=('red')) #Units: Jy/Beam, Resolution: 45''
    
    #Figure properties
    Fig.show_grayscale(vmax=0.027,vmin=0., stretch='linear',invert='TRUE')

    Fig.tick_labels.set_font(size='medium')
    Fig.add_colorbar()
    #Fig.colorbar.set_axis_label_text('$850\,(K)$')
    Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    return
    

def U(small,large):
    #deffine image in APLPY
    Fig = aplpy.FITSFigure(large)
    cset1 = Fig.show_contour(small, levels=(0.001,0.005,0.0075,0.01), linewidth=4, colors=('red')) #Units: Jy/Beam, Resolution: 45''
    cset2 = Fig.show_contour(large, levels=(0.001,0.005,0.0075,0.01), linewidth=4, colors=('lime')) #Units: Jy/Beam, Resolution: 45''
    
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=30, alpha=1)
    YSO.w40_radio(Fig)
    
    Fig.ticks.set_color('black')
    return

def V(COblue,COred,S2,mask):
    #produce maps of large CO emission for W40 
    Fig = aplpy.FITSFigure(S2)
    #Fig.show_colorscale(vmax=100,vmin=0, stretch='linear')
    Fig.show_grayscale(vmax=0.027,vmin=0.00054, stretch='sqrt',invert=True)
    Fig.add_colorbar()
    #markers
    #YSO.w40(Fig)
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
    #contours
    cset1 = Fig.show_contour(COred, levels=(5,20,40,60,80), linewidth=2, colors=('red'))
    cset2 = Fig.show_contour(COblue, levels=(20,40,60,80), linewidth=2, colors=('blue'))
    #cset3 = Fig.show_contour(S2, levels=(0.0027,0.0081,0.027), linewidth=2, colors=('black','black','black'))
    #cset3 = Fig.show_contour(mask, levels=(0.5,2), linewidth=2, colors=('k'))
    #figure properties
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    #Fig.colorbar.set_axis_label_text('line intensity\,(K)')
    Fig.colorbar.set_axis_label_text('SCUBA-2 Flux density (Jy per pixel)')
    return

def W(CO,S2,mask):
   #produce maps of large CO emission for W40 
    Fig = aplpy.FITSFigure(CO)
    Fig.show_colorscale(vmax=0.5,vmin=0, stretch='linear') 
    Fig.add_colorbar()
    cset3 = Fig.show_contour(S2, levels=(0.0027,0.0081,0.027), linewidth=2, colors=('black','black','black'))
    cset3 = Fig.show_contour(mask, levels=(0.5,2), linewidth=2, colors=('k'))
    #figure properties
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    Fig.colorbar.set_axis_label_text('Intergrated line intensity\,(K kms$^{-1}$) per pixel')
    Fig.recenter(277.8615,-2.0700, height=0.30 ,width=0.375)
    return

def X(H,S2,VLA,blue,red):
    #Herschel maps with contours
    Fig = aplpy.FITSFigure(H)
    Fig.show_colorscale(vmax=10000,vmin=0, stretch='linear') 
    Fig.add_colorbar()
    #markers
    #YSO.w40(Fig)
    markers = np.loadtxt('/data/damian/run26/SMM/S2-FF-CO/SMM-cores2.tab')
    ra = markers[:,1]
    dec = markers[:,2]
    #Fig.show_markers(ra,dec, edgecolor='k', facecolor='k', marker='*',s=200, alpha=0.25)
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='k', marker='*',s=200, alpha=1)
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    YSO.w40(Fig)
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)
    #contours
    cset3 = Fig.show_contour(S2, levels=(0.0025,0.0075,0.025), linewidth=2, colors=('black','black','black'))
    cset4 = Fig.show_contour(VLA, levels=(0.05,0.25), linewidth=2, colors=('white'))
    #cset1 = Fig.show_contour(red, levels=(5,25,75), linewidth=2, colors=('red'))
    #cset2 = Fig.show_contour(blue, levels=(20,40,60,80), linewidth=2, colors=('blue'))
    #figure properties
    Fig.ticks.show()
    Fig.ticks.set_color('black')
    Fig.colorbar.set_axis_label_text('S$_{70\mu m}$ flux (MJy/Sr)')
    Fig.recenter(277.8615,-2.0700, height=0.30 ,width=0.375)
    return

def Y(s450,s850,m1,m2):
    #deffine image in APLPY
    big = plt.figure(figsize=(18, 10))

    S450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.1,0.4,0.8])
    S850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.5,0.1,0.4,0.8])

    S450.show_grayscale(vmax=0.0213,vmin=-10., stretch='sqrt')
    S850.show_grayscale(vmax=0.0213,vmin=-10., stretch='sqrt')

    S850.axis_labels.hide_y()
    S850.tick_labels.hide_y()
    S850.ticks.show()
    S450.ticks.show()
    S450.ticks.set_color('black')
    S850.ticks.set_color('black')

    S450.add_label(0.1, 0.95, r'IR1 850$\mathrm{\mu}m$', relative=True)
    S850.add_label(0.1, 0.95, r'IR2 850$\mathrm{\mu}m$', relative=True)

    IR1sig = 0.0150
    IR2sig = 0.0111

    S450.show_contour(s450, levels=(IR1sig,3*IR1sig), linewidth=2, colors=('red'))
    S850.show_contour(s850, levels=(IR2sig,3*IR2sig), linewidth=2, colors=('red'))
   
    S450.show_contour(m1, levels=(0.5,2), linewidth=2, colors=('blue'))
    S850.show_contour(m2, levels=(0.5,2), linewidth=2, colors=('blue'))

    #S450.add_colorbar()
    #S450.colorbar.set_location('top')
    #S450.colorbar.set_axis_label_text('IR1 850$\mathrm{\mu}$m (Jy/pixel)')

    #S850.add_colorbar()
    #S850.colorbar.set_location('top')
    #S850.colorbar.set_axis_label_text('IR2 850$\mathrm{\mu}$m (Jy/pixel)')
    return

def Z(CD): #OB star checking
    Fig = aplpy.FITSFigure(CD)
    Fig.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
    Fig.add_colorbar()
    Fig.colorbar.set_axis_label_text('Column Density (H$_{2}$ cm$^{-2}$)') 

    #Plot YSOs
    #YSO.c2dGBS_v2(Fig)

#PLOT FW contours
#OPEN up the region specific clump list

    clumpcat = '/data/damian/run31/FW/cat/OrionA_20150227_clump_cat.fit'
    clumps = '/data/damian/run31/FW/clumps/OrionA_20150227_clumps.sdf'

    data = pyfits.open(clumpcat)
    DATA = data[1].data

    for j in DATA:
        I = int(j[0])
        print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
        thresh = '/data/damian/run31/mask/thresh'+str(I)+'.sdf'
        mask = '/data/damian/run31/mask/mask_'+str(I)+'.sdf'
        maskFITS = '/data/damian/run31/mask/mask_'+str(I)+'.fits'
        magic = '/data/damian/run31/mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumps,thresh,I,I,'bad','bad','> /dev/null')
        os.system(cmd)
    ### process 2 - remove 0value data
        cmd = '%s/nomagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
        os.system(cmd)
    ### recreate mask of clump i - mask0
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
        os.system(cmd)
        
        ndf2fits.ndf2fits(mask)

        Fig.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

        cmd = 'rm %s %s %s'%(thresh,mask,magic)
        os.system(cmd)

    #Deffine OB star data
    data1 = np.loadtxt('/data/damian/run27/OBstars_list.txt')
    RA_OB = data1[:,0]	
    DEC_OB = data1[:,1]

    data2 = np.loadtxt('/data/damian/run27/OBstars.txt',dtype='string')
    RA_OBx = map(float, data2[:,2])	
    DEC_OBx = map(float, data2[:,3])

    data1 = np.loadtxt('/data/damian/data/OBstars/obcatalog_deg.txt')
    RA_OB_R = data1[:,0]	
    DEC_OB_R = data1[:,1]


    Fig.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
    Fig.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
    Fig.show_markers(RA_OB_R,DEC_OB_R, edgecolor='r', facecolor='none', marker='*',s=100, alpha=1)
    
    Fig.tick_labels.set_xformat('hh:mm:ss')

    return

def AA(IR1,IR2,IR3):
    big = plt.figure(figsize=(20, 10))

    SIR1 = aplpy.FITSFigure(IR1,figure=big,subplot=[0.1,0.1,0.28,0.8])
    SIR2 = aplpy.FITSFigure(IR2,figure=big,subplot=[0.38,0.1,0.28,0.8])
    SIR3 = aplpy.FITSFigure(IR3,figure=big,subplot=[0.66,0.1,0.28,0.8])

    SIR1.show_grayscale(vmax=0.7,vmin=-0.07, stretch='sqrt',invert='TRUE')
    SIR2.show_grayscale(vmax=0.03,vmin=-0.001, stretch='sqrt',invert='TRUE')
    SIR3.show_grayscale(vmax=0.03,vmin=-0.001, stretch='sqrt',invert='TRUE')

    center = [277.8995,-2.08200,0.2,0.20]

    SIR1.recenter(center[0],center[1], height=center[2], width=center[3])
    SIR2.recenter(center[0],center[1], height=center[2], width=center[3])
    SIR3.recenter(center[0],center[1], height=center[2], width=center[3])

    SIR1.show_contour(IR1, levels=((0.0093*5.),(0.0093*15.)), colors=('blue'))
    SIR2.show_contour(IR2, levels=((0.0005215*5.),(0.0005215*15.)), colors=('blue'))
    SIR3.show_contour(IR3, levels=((0.000584*5.),(0.000584*15.)), colors=('blue'))

    SIR2.axis_labels.hide_y()
    SIR2.tick_labels.hide_y()
    SIR3.axis_labels.hide_y()
    SIR3.tick_labels.hide_y()

    SIR1.ticks.show()
    SIR1.ticks.show()
    SIR1.ticks.set_color('black')
    SIR1.ticks.set_color('black')

    SIR2.ticks.show()
    SIR2.ticks.show()
    SIR2.ticks.set_color('black')
    SIR2.ticks.set_color('black')

    SIR3.ticks.show()
    SIR3.ticks.show()
    SIR3.ticks.set_color('black')
    SIR3.ticks.set_color('black')

    SIR1.add_colorbar()
    SIR1.colorbar.hide()
    SIR1.add_beam(14.5*u.arcsecond,14.5*u.arcsecond,0,color='black')
    SIR1.add_scalebar(412 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))

    SIR2.add_colorbar()
    SIR2.colorbar.hide()
    SIR2.add_beam(14.5*u.arcsecond,14.5*u.arcsecond,0,color='black')
    SIR2.add_scalebar(412 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))

    SIR3.add_colorbar()
    SIR3.colorbar.set_axis_label_text('850$\mathrm{\mu}$m (Jy/pixel)')
    SIR3.add_beam(14.5*u.arcsecond,14.5*u.arcsecond,0,color='black')
    SIR3.add_scalebar(412 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))

    SIR1.add_label(0.9, 0.95, 'IR1', relative=True)
    SIR2.add_label(0.9, 0.95, 'IR1', relative=True)
    SIR3.add_label(0.9, 0.95, 'IR1', relative=True)

    return 

####################################################
#plot background DUST map - subtitute in as appropriate 

#450 files
south_s450 = 'aquila/serpens_south/skyloop/aquila_s450_skyloopfull1_jypix.fits'
south_s450_IR2 = 'aquila/serpens_south/IR2/Aquila_extS2nosm_s450-4am.fits'
main_s450 = 'maps/serpens/serpensmain/IR2/SerpensMain_20141223_s450_IR2extmask_s2_cal_JypixJH.fits'
mwc297_s450 = 'serpens/MWC297/IR1/SerpensMWC297_20140331_IR1_s450_freefree_3d.fits'    
E1_s450 = 'serpens_E1/SerpensE1_20130425_s450_IR1_JH.fits'
south_s450_var ='aquila/serpens_south/IR2/Aquila_20141105_s450_VAR.fits'
south450_freefree = 'aquila/serpens_south/IR2/Aquila_extS2nosm_s450-4am-fullfreefree.fits'
E_450_IR2 = '/data/damian/run27/input_s450/SerpensE_20141219_s450_IR2extmask_s2_cal_JypixJH.fits'
N_450_IR2 = 'serpens/serpensN/IR2/SerpensN_20141223_s450_IR2extmask_s2h2_cal_JypixJH.fits'

#850 files
south_s850 = 'aquila/serpens_south/skyloop/aquila_s850_skyloopfull1_jypix.fits'
main_s850 = 'serpens/serpensmain/IR2/SerpensMain_20150326_850_IR2_noco_HKJypix_col.fits'
mwc297_s850 = 'serpens/MWC297/IR1/SerpensMWC297_20140414_IR1_s450_freefree+jet.fits'
mwc297_s850_IR1 = 'serpens/MWC297/IR1/SerpensMWC297_20121221_s850_IR1_JH.fits'
mwc297_nonegIR2 = 'serpens/MWC297/IR2/SerpensMWC297_20141219_s850_IR2extmask_s2_cal_JypixJH_negMSK.fits'
mwc297_nonegIR1 = 'serpens/MWC297/IR1/SerpensMWC297_20121221_s850_IR1_JH_negMSK.fits'
south_s850_IR2 = 'aquila/serpens_south/IR2/Aquila_noco_extS2nosm_s850-4am.fits'

south_s850_IR1 = 'aquila/serpens_south/IR1/Aquila_20130416_s850_IR1_JH.fits'
south_s850_IR2o = 'aquila/serpens_south/IR2/Aquila_20140811_s850_IR2_auto_JH.fits'
south_s850_IR3 = 'aquila/serpens_south/IR3/Aquila_s850_automask_IR3_test_cal_JypixJH.fits'

south_s850_IR2freefree = 'aquila/serpens_south/IR2/Aquila_noco_extS2nosm_s850-4am-fullfreefree.fits'
E1_s850 = 'serpens_E1/SerpensE1_20121214_s850_IR1_JH.fits'
E_850_IR2 = '/data/damian/run27/input_s850/SerpensE_20141219_s850_IR2extmask_s2_cal_JypixJH_col.fits'
N_850_IR2 = 'serpens/serpensN/IR2/SerpensN_20141223_s850_IR2extmask_s2h2_cal_JypixJH.fits'

south_s850_var ='aquila/serpens_south/IR2/Aquila_20141105_s850_VAR.fits'
south_s850_largpix = '/data/damian/run25/maps/s21/freefreefraction/Aquila-850ats21specs.fits'

#Herschel data
H70 = 'aquila/serpens_south/Herschel/aquilaM2-070.fits'
H160 = 'aquila/serpens_south/Herschel/aquilaM2-160.fits'
H250 = 'aquila/serpens_south/Herschel/aquilaM2-250.fits'
H350 = 'aquila/serpens_south/Herschel/aquilaM2-350.fits'
H500 = 'aquila/serpens_south/Herschel/aquilaM2-500.fits'

#24uM SPITZER
mwc297_s24a = 'spitzer/SER_E_2_A_MIPS1_mosaic.fits'
mwc297_s24b = 'spitzer/SER_E_3_A_MIPS1_mosaic.fits'

#CO files
W40_mask = '/data/damian/maps/aquila/serpens_south/CO/COmask.fits'
W40_red = '/data/damian/maps/aquila/serpens_south/CO/W40/jcmth20120228_north.fits'
W40_blue = '/data/damian/maps/aquila/serpens_south/CO/W40/jcmth20120228_south.fits'
CO_mask = '/data/damian/maps/aquila/serpens_south/CO/CO_mask_cont.fits'
W40_int = '/data/damian/maps/aquila/serpens_south/CO/W40/W4020150704_reduced001_CO3-2_intmap.fits'
W40_int_blue = '/data/damian/maps/aquila/serpens_south/CO/W40/outflows/W4020150704_reduced001_CO3-2_intmap_BLUE.fits'
W40_int_red = '/data/damian/maps/aquila/serpens_south/CO/W40/outflows/W4020150704_reduced001_CO3-2_intmap_RED.fits'
W40_CO = '/data/damian/reduction/Aquila_noco/mosaics/GOOD/20120705_00034_s850_extmask_s2nosm_cal_Jypix_mos.fits'

#freefree
HII_large = '/data/damian/maps/aquila/serpens_south/radio/aquila_VLAss_21cm.fits'
HII_small = '/data/hatchell/SCUBA2/GBS/work/VLA_8_7GHz_W40_JHedit.fits'
s21as850 = '/data/damian/run25/maps/s21/freefreefraction/Aquila_S21scaled850.fits'

large850 = '/data/damian/run25/maps/s21/FF/Aquila_S21scaled850.fits'
small850 = '/data/damian/run25/maps/s21/s850freefree_VLAconv.fits'

freefree450 = '/data/damian/analysis/Freefree_contamination/maps/tests/comb/s450freefree_comb_v2.fits'
freefree850 = '/data/damian/analysis/Freefree_contamination/maps/tests/comb/s850freefree_comb_v2.fits'

#temp. files
mwc297_temp='serpens/MWC297/IR1/SerpensMWC297_20140331_IR1_freefree_3d450850temp%5.fits'
mwc297_temp_old='/data/damian/maps/tempmaps/SerpensMWC297_20140414_IR1_freefree+jet450850temp%5_cshV.fits'
w40_temp = '/data/damian/run26/input/S2-FF-CO/Aquila_noco_s2nosm-4am-FFtemperature.fits'
w40_temp_VLA = '/data/damian/run25/maps/s21/temp_21.fit'
w40_temp_FF_VLA = '/data/damian/run25/maps/s21/temp-FF_21.fit'

#cump files
clump = 'serpens/MWC297/IR1/SerpensMWC297_IR1_collapse.fits'
clumpW40 = '/data/damian/run26/output/S2-FF-CO/Aquila_noco_s850-4am-FF_clumps5COL.fits'

#Mask files 
maskIR1 = 'serpens/MWC297/masks/serpensmwc297_s2+h_mask_850_magic.fits'
maskIR2 = 'serpens/MWC297/IR2/SerpensMWC297_20140820_850_IR2_mask_JH_cont.fits'
maskW40450 = 'aquila/serpens_south/IR2/Aquila_s2nosm_mask_450_cont.fits'
maskW40850 = 'aquila/serpens_south/IR2/Aquila_s2nosm_mask_850_cont.fits'
masklargepix = '/data/damian/run25/maps/s21/ratiomask.fits'

#alpha
alpha = '/data/damian/temp_MKR_kernel/output/kernel/map/Aquila_noco_s2nosm-4am-FF/Aquila_noco_s2nosm-4am-FF_Salpha.fits'

#CD
#CDW40 = '/data/damian/run26/input/S2-FF-CO/Aquila-noco4am-ColumnD.fits'
CDW40 = '/data/damian/run26/input/S2-FF-CO/combinedCDmap.fits'
CDOrionB2023 = '/data/damian/run27/CDmaps/OrionB_N2023_CDcomb.fits'
CDOrionB2068 = '/data/damian/run27/CDmaps/OrionB_N2068_CDcomb.fits'
CDOrionA = '/data/damian/run27/CDmaps/OrionA_CDcomb_MSK.fits'

#YSOs
W40_YSOs = '/data/damian/run26/output/surface/Aquila850_YSOsurface_r300_units.fits'
W40_proto = '/data/damian/run26/output/surface/Aquila850_protostarsurface_r300_units.fits'
W40_PMS = '/data/damian/run26/output/surface/Aquila850_PMSsurface_r300_units.fits'

SS_YSO = '/data/damian/run27/YSO/surface/Aquila_YSOsurface_r300_units.fits'

####################################################
############        MAKE MAPS HERE     #############
####################################################

#A Plot publication quality 450 and 850 data together
#A(south_s450_IR2,south_s850_IR2,W40_YSOs,HII_large,W40_mask)
#A(E_450_IR2,E_850_IR2,W40_YSOs,HII_large,W40_mask)
#A(main_s450,main_s850,W40_YSOs,HII_large,W40_mask)

#a Appendix SCUBA-2 data + Variance 
#a(south_s450_IR2,south_s850_IR2,maskW40450,maskW40850,south_s450_var,south_s850_var)

#B(mwc297_temp,mwc297_s850,mwc297_s24a,mwc297_s24b)

#plots temp map, with Herschel Contours
#C(w40_temp,H70)

#D(mwc297_s450, mwc297_s850, mwc297_temp)

#plot clump map
#E(clumpW40)

#Map with distrubution of YSO markers
#F(E_450_IR2,W40_YSOs)

#SCUBA-2 plus Spitzer24
#G(mwc297_s24a,mwc297_s850)

#H(mwc297_temp_old,mwc297_temp)

#Plot CO outflow lines over a map. Required specific maps of red and blue shifted CO emission
#I(south_s850_IR2,W40_red,W40_blue) 

#Herschel map plotter
#J(H350,'350')

#Free-free contamination
#K(south_s450_IR2,HII_small,south450_freefree)

#W40 with CO observation regions + HII and YSOs
#L(south_s850_IR2,HII_large)

#Hii large + Herschel
#M(w40_temp,south_s850_IR2,HII_large,HII_small,H70)

#Radio plot of W40
#N(HII_small,south_s850_IR2)

#plot W40 850um dust and HII next to each other
#m(south_s850_largpix,s21as850,masklargepix)
#m(w40_temp_VLA,w40_temp_FF_VLA,masklargepix)

#plot map of alpha/CD
#O(CDW40,W40_YSOs)
#O(alpha,W40_YSOs)

### Plot W40 free sources on 450/850 maps with common contours
#P(south_s450_IR2,south_s850_IR2,freefree450,freefree850)

#Plot of 3.6cm
#Q(HII_small,south_s850_IR2)

#Surface plot of YSOs
#R(south_s850_IR2,W40_proto,SS_YSO,H70,maskW40850)

#temp map + H70
#S(w40_temp,H70,HII_large,freefree850)

###Poster image
#T(south_s850_IR2,H70,HII_large)

#U(small850,large850)

#CO map with additional features
#V(W40_int_blue,W40_int_red,south_s850_IR2,W40_mask)
#W(W40_CO,south_s850_IR2,W40_mask)

#plot Hershcel data with contours
#X(H70,south_s850_IR2,HII_large,W40_int_blue,W40_int_red)

#PLOT comparaison of IR1IR2
#Y(mwc297_nonegIR1,mwc297_nonegIR2,maskIR1,maskIR2)

#PLOT CD maps with OB stars
Z(CDOrionA)

#PLOT comparision IR1,IR2,IR3
#AA(south_s850_IR1,south_s850_IR2o,south_s850_IR3)

plt.show()
