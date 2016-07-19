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

import YSO

####################################################
#Plotting functions - wavelength specific 


def S2map(s450,s850,mask):
    #deffine image in APLPY
    big = plt.figure(figsize=(16,12))
    #18,11.5 #main
    #18,17 #oOuth
    #16,12 East
    
    S450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.1,0.49,0.8])
    S850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.5,0.1,0.49,0.8])

    #W40
    S450.show_grayscale(vmax=0.213,vmin=0., stretch='linear',invert='TRUE')
    S850.show_grayscale(vmax=0.015,vmin=0., stretch='linear',invert='TRUE')
    #NOISE 450,850
    #AQUILA 0.185,0.026
    #MAIN 0.131,0.018
    #Serp E 0.0213, 0.0015
    #Serp N 0.0108, 0.0016

    S850.axis_labels.hide_y()
    S850.tick_labels.hide_y()
    S850.ticks.show()
    S450.ticks.show()
    S450.ticks.set_color('black')
    S850.ticks.set_color('black')
    
    #North
    #S450.recenter(279.73625,0.41488888888, height=0.45 ,width=0.25)
    #S850.recenter(279.73625,0.41488888888, height=0.45 ,width=0.25)
    #East
    S450.recenter(279.4275,-1.46108333333, height=0.8 ,width=0.8)
    S850.recenter(279.4275,-1.46108333333, height=0.8 ,width=0.8)
    #NH3
    #S450.recenter(277.27425,0.556111111111, height=0.5 ,width=0.55)
    #S850.recenter(277.27425,0.556111111111, height=0.5 ,width=0.55)
    #Main
    #S450.recenter(277.45625,1.25441666667, height=0.25 ,width=0.4)
    #S850.recenter(277.45625,1.25441666667, height=0.25 ,width=0.4)
    #SS
    #S450.recenter(277.406666667,-1.90194444444, height=0.85 ,width=0.28)
    #S850.recenter(277.406666667,-1.90194444444, height=0.85 ,width=0.28)
    #W40
    #S450.recenter(277.8995,-2.09200, height=0.60 ,width=0.40)
    #S850.recenter(277.8995,-2.09200, height=0.60 ,width=0.40)

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

    n450, n850 = 0.0213, 0.0015
    #W40 - 0.0173,0.0025
    #Main - 0.0131,0.0018
    
    #contours
    cset450 = S450.show_contour(s450, levels=(n450,100), linewidth=3, colors=('black'))
    cset850 = S850.show_contour(s850, levels=(n850,100), linewidth=3, colors=('black'))

    #CO
    csetCO = S850.show_contour(mask, levels=(0.5,2), linewidth=2, colors=('blue'))

    #cset5 = S450.show_contour(YSO,levels=(15,45,75),colors=('r'))
    #cset1 = S850.show_contour(HII1, levels=(0.01,0.05,0.075,0.1,0.2), linewidth=4, colors=('blue')) #Units: Jy/Beam, Resolution: 45''
    OBmarkers = np.loadtxt('/data/damian/run27/OBstars_list.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    S450.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)
    S850.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)

    #YSO.c2dGBS_v2(S850)
    
    return

def YSOmap(fig,surfaceA,H,mask):
    #produce surface map of YSO distributions
    Fig = aplpy.FITSFigure(fig)
    
    levels = (15,45,75,100,150)

    #cset_H = Fig.show_contour(H, levels=(200,10000),colors='blue')
    cset1 = Fig.show_contour(surfaceA, levels=levels,colors='k')
    #cset2 = Fig.show_contour(surfaceB, levels=levels,colors='k')
    #cset_mask = Fig.show_contour(mask, levels=(0.5,10),colors='k')

    #YSO.w40(Fig)
    YSO.c2dGBS(Fig)

    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, edgecolor='k', facecolor='yellow', marker='*',s=250, alpha=1)

    #SerpSouth
    #ra = 277.299583333
    #dec = -2.06388888889
    #Fig.show_markers(ra,dec, layer='W40a', edgecolor='k', facecolor='yellow', marker='*',s=150, alpha=1)

    Fig.show_grayscale(vmax=0.036,vmin=0.0015, stretch='sqrt',invert='TRUE')
    #Fig.recenter(277.8995,-2.09200, height=0.6 ,width=0.6) #W40
    #Fig.recenter(277.4406666667,-1.90194444444, height=0.85 ,width=0.5) #SS
    Fig.recenter(277.45625,1.25441666667, height=0.4 ,width=0.4) #Main
    Fig.add_colorbar()
    #Fig.colorbar.set_location('top')
    Fig.colorbar.set_axis_label_text('850$\mu m$ (Jy/pixel)')



####################################################
#plot background DUST map - subtitute in as appropriate 

#450 files
south_s450_IR2 = 'aquila/serpens_south/IR2/Aquila_extS2nosm_s450-4am.fits'
main_s450 = 'serpens/serpensmain/IR2/SerpensMain_20141223_s450_IR2extmask_s2_cal_JypixJH.fits'
mwc297_s450 = 'serpens/MWC297/IR1/SerpensMWC297_20140331_IR1_s450_freefree_3d.fits'    
south450_freefree = 'aquila/serpens_south/IR2/Aquila_extS2nosm_s450-4am-fullfreefree.fits'
E_450_IR2 = '/data/damian/run27/input_s450/SerpensE_20141219_s450_IR2extmask_s2_cal_JypixJH.fits'
N_450_IR2 = '/data/damian/run27/input_s450/SerpensN_20141219_s450_IR2extmask_s2_cal_JypixJH.fits'

#850 files
main_s850 = 'serpens/serpensmain/IR2/SerpensMain_20150326_850_IR2_noco_HKJypix_col.fits'
mwc297_s850 = 'serpens/MWC297/IR1/SerpensMWC297_20140414_IR1_s450_freefree+jet.fits'
south_s850_IR2freefree = 'aquila/serpens_south/IR2/Aquila_noco_extS2nosm_s850-4am-fullfreefree.fits'
E1_s850 = 'serpens_E1/SerpensE1_20121214_s850_IR1_JH.fits'
E_850_IR2 = '/data/damian/run27/input_s850/SerpensE_20141219_s850_IR2extmask_s2_cal_JypixJH_col.fits'
N_850_IR2 = '/data/damian/run27/input_s850/SerpensN_20141219_s850_IR2extmask_s2_cal_JypixJH_col.fits'

#Herschel data
H70 = 'aquila/serpens_south/Herschel/aquilaM2-070.fits'

#Mask files 
mask = 'serpens/MWC297/serpensmwc297_s2+h_mask_850_magic.fits'
maskW40450 = 'aquila/serpens_south/IR2/Aquila_s2nosm_mask_450_cont.fits'
maskW40850 = 'aquila/serpens_south/IR2/Aquila_s2nosm_mask_850_cont.fits'
Main_CO = 'serpens/serpensmain/CO/SerpMain20070705_CO3-2_mask_cont.fits'
W40_CO = '/data/damian/maps/aquila/serpens_south/CO/COmask.fits'

#YSOs
W40_YSOs = '/data/damian/run26/output/surface/Aquila850_YSOsurface_r300_units.fits'
W40_proto = '/data/damian/run26/output/surface/Aquila850_protostarsurface_r300_units.fits'
W40_PMS = '/data/damian/run26/output/surface/Aquila850_PMSsurface_r300_units.fits'

SS_YSO = '/data/damian/run27/YSO/surface/Aquila_YSOsurface_r300_units.fits'
Main_YSO = '/data/damian/run27/YSO/surface/SerpensMain_YSOsurface_r300_units.fits'

####################################################
############        MAKE MAPS HERE     #############
####################################################

#A Plot publication quality 450 and 850 data together
#S2map(main_s450,main_s850,W40_CO)
#S2map(south450_freefree,south_s850_IR2freefree,W40_CO)
S2map(E_450_IR2,E_850_IR2,W40_CO)
#S2map(N_450_IR2,N_850_IR2,W40_CO)
#S2map(main_s450,main_s850,Main_CO)

#Surface plot of YSOs
#YSOmap(main_s850,Main_YSO,H70,maskW40850)

plt.show()
