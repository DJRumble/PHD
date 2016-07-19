#20150930

#A mapping script desgined specifically for the paper2 analysis

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

#My modules
import YSO

####################################################
#Plotting functions - wavelength specific 


#TRI PANNEL PLOT FUNCTION
def A(s450,s850,temp,sig450,sig850,surface):
    #tri panel map - should tile all maps at once. 1

    big = plt.figure(figsize=(12, 17))
    #deffine image in APLPY
    Fig_850 = aplpy.FITSFigure(s850,figure=big,subplot=[0.1,0.64,0.8,0.26])
    Fig_450 = aplpy.FITSFigure(s450,figure=big,subplot=[0.1,0.37,0.8,0.26])
    Fig_temp = aplpy.FITSFigure(temp,figure=big,subplot=[0.1,0.1,0.8,0.26])
    #Deffine contours
    cset1 = Fig_850.show_contour(s850, levels=(sig850,300*sig850), linewidth=3, colors=('black','black'))
    cset2 = Fig_450.show_contour(s450, levels=(sig450,300*sig450), linewidth=3, colors=('black','black'))
    #cset3 = Fig_temp.show_contour(temp, levels=(12,20,30), linewidth=3, colors=('black'))

    #Figure properties
    Fig_450.show_grayscale(vmax=10*sig450,vmin=0.0,invert=True, stretch='sqrt')
    Fig_450.add_label(0.07, 0.95, r'b) $450\mu m$',size='large', relative=True)
    Fig_850.show_grayscale(vmax=10*sig850,vmin=0.0,invert=True, stretch='sqrt')
    Fig_850.add_label(0.07, 0.95, r'a) $850\mu m$',size='large', relative=True)
    Fig_temp.show_colorscale(vmax=35.0,vmin=10.0, stretch='sqrt',cmap='gist_heat')
    Fig_temp.add_label(0.05, 0.95, r'c) $T_{d}$',size='large', relative=True)

    #PLOT ysos
    YSO.c2dGBS(Fig_850)

    #PLOT contours
    cset850 = Fig_temp.show_contour(s850, levels=(sig850,100), linewidth=3, colors=('black'))
    csetYSO = Fig_temp.show_contour(surface, levels=(20,60,100), linewidth=3, colors=('red'))

    Fig_450.tick_labels.set_font(size='medium')
    Fig_450.add_colorbar()
    #Fig_450.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig_850.tick_labels.set_font(size='medium')
    Fig_850.add_colorbar()
    #Fig_850.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)
    Fig_temp.tick_labels.set_font(size='medium')
    Fig_temp.add_colorbar()
    #Fig_temp.recenter(277.074583333,-3.74686111111,height=0.25 ,width=0.65)

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

#TEMPERATURE plot
def B(s850,temp,sig850,pos):
    Fig_temp = aplpy.FITSFigure(temp)
    Fig_temp.show_colorscale(vmax=35.0,vmin=10.0, stretch='linear',cmap='gist_heat')

    #PLOT ysos
    YSO.c2dGBS(Fig_temp)

    data = np.loadtxt('OBstars_list.txt')
    RA_OB = data[:,0]	
    DEC_OB = data[:,1]
    Fig_temp.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)

    #PLOT contours
    cset850 = Fig_temp.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('b'))

    Fig_temp.ticks.show()
    Fig_temp.ticks.set_color('black')
    Fig_temp.tick_labels.set_font(size='medium')
    Fig_temp.add_colorbar()
    Fig_temp.colorbar.set_axis_label_text('T$_{d}$ (K)') 
    #Fig_temp.recenter(pos[0],pos[1], height=pos[2] ,width=pos[3])




####################################################
############        MAKE MAPS HERE     #############
####################################################

#ENTER MAP NAME

Map = 'CrA'

if Map == 'N':
    N_450 = 'input_s450/SerpensN_20141219_s450_IR2extmask_s2_cal_JypixJH.fits'
    N_850 = 'input_s850/SerpensN_20141219_s850_IR2extmask_s2_cal_JypixJH.fits'
    N = 'input_temp/SerpensN_20141219-autotemperature.fits'
    #N
    sig450, sig850 = 0.0108,0.00161
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [279.74,0.46,0.41,0.35]

    #Plot maps
    #TEMP PLOT
    B(N_850,N,sig850,pos)

elif Map == 'E':
    E_450 = 'input_s450/SerpensE_20141219_s450_IR2extmask_s2_cal_JypixJH.fits'
    E_850 = 'input_s850/SerpensE_20141219_s850_IR2extmask_s2_cal_JypixJH.fits'
    E = 'input_temp/SerpensE_20141219-autotemperature.fits'
    #E
    sig450 = 0.0212848097458
    sig850 = 0.00151780863501
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [279.4,-1.45,1,0.8]

    #Plot maps
    #TEMP PLOT
    B(E_850,E,sig850,pos)

elif Map == 'MWC297':
    #MWC297_IR2_450 = 'input_s450/SerpensMWC297_20141219_s450_IR2extmask_s2_cal_JypixJH.fits'
    MWC297_IR1_450 = 'input_s450/SerpensMWC297_20140414_IR1_s450_freefree+jet.fits'
    #MWC297_IR2_850 = 'input_s850/SerpensMWC297_20141219_s850_IR2extmask_s2_cal_JypixJH.fits'
    MWC297_IR1_850 = 'input_s850/SerpensMWC297_20140414_IR1_s850_freefree+jet.fits'
    #MWC297_IR2 = 'input_temp/SerpensMWC297_20141219-autotemperature.fits'
    MWC297_IR1 = 'input_temp/SerpensMWC297_IR1_freefree+jettemperature.fits'
    #MWC297
    sig450 = 0.082 #IR1 #0.0260826450613
    sig850 = 0.011 #IR1 #0.00177450362315
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [277.074583333,-3.74686111111,0.25,0.65]

    #Plot maps
    #TEMP PLOT
    B(MWC297_IR1_850,MWC297_IR1,sig850,pos)

elif Map == 'Aquila':
    Aquila_450 = 'input_s450/Aquila_extS2nosm_s450-4am-fullfreefree.fits'
    Aquila_850 = 'input_s850/Aquila_noco_extS2nosm_s850-4am-fullfreefree.fits'
    Aquila = 'input_temp/Aquila_noco-autotemperature.fits'
    #Aquila
    sig450 = 0.0104480643359
    sig850 = 0.00150916751884
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [277.9025,-2.13200,0.65,0.48] #W40
    #pos = [277.485,-1.92200,0.65,0.52] #SerpSouth

    #Plot maps
    #TEMP PLOT
    B(Aquila_850,Aquila,sig850,pos)

elif Map == 'Main':
    main_450 = 'input_s450/SerpensMain_20141223_s450_IR2extmask_s2_cal_JypixJH.fits'
    main_850 = 'input_s850/SerpensMain_20141223_s850_IR2extmask_s2_cal_JypixJH.fits'
    main_noco850 = 'input_s850/SerpensMain_20150326_850_IR2_noco_JypixHK.fits'
    mainnoco = 'input_temp/SerpensMain_20150326-autotemperature.fits' #???????
    main = 'input_temp/SerpensMain_20141223-autotemperature.fits'
    #Aquila
    sig450 = 0.0130978557306 #IR2
    sig850 = 0.00182226342022 #IR2
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    #pos = [277.4767916667,1.2355,0.25,0.25] #Main
    pos = [277.3,0.5,0.5,0.5] #NH3 ????

    #Plot maps
    #TEMP PLOT
    B(main_850,main,sig850,pos)

elif Map == 'Oph':
    Oph_450 = 'input_s450/OphSco_Main_20150318_450_IR2_ext_JypixHK.fits'
    Oph_850 = 'input_s850/OphSco_Main_20150318_850_IR2_ext_JypixHK.fits'
    Oph = 'input_temp/OphSco_Main-autotemperature.fits'
    #Aquila
    sig450 = 0.0163624005264
    sig850 = 0.00240057348201
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [277.9025,-2.13200,0.65,0.48] #
    
    #Plot maps
    #TEMP PLOT
    B(Oph_850,Oph,sig850,pos)
elif Map == 'IC5146':
    s450 = 'input_s450/IC5146_20150225_450_IR2_ext_HKJypix_col.fits'
    s850 = 'input_s850/IC5146_20150225_850_IR2_ext_HKJypix_col.fits'
    temp = 'input_temp/IC5146_20150225-autotemperature.fits'
    surface = 'YSO/surface/IC5146_YSOsurface_r300_units.fits'
    #Aquila
    sig450 = 0.0181
    sig850 = 0.0016
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [277.9025,-2.13200,0.65,0.48] #
    
    #Plot maps
    #TEMP PLOT
    A(s450,s850,temp,sig450,sig850,surface)
    #B(s850,temp,sig850,pos)
elif Map == 'Lupus':
    s450 = 'input_s450/Lupus_20150318_850_IR2_ext_HKJypix_col.fits'
    s850 = 'input_s850/Lupus_20150318_850_IR2_ext_HKJypix_col.fits'
    temp = 'input_temp/Lupus_20150318-autotemperature.fits'
    surface = 'YSO/surface/Lupus_YSOsurface_r300_units.fits'
    #Aquila
    sig450 = 0.0392
    sig850 =  0.0027 #0.0019
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [235.7242,-34.1525,0.1,0.12] #
    
    #Plot maps
    #TEMP PLOT
    #A(s450,s850,temp,sig450,sig850,surface)
    B(s850,temp,sig850,pos)
elif Map == 'CrA':
    s450 = 'input_s450/Lupus_20150318_850_IR2_ext_HKJypix_col.fits'
    s850 = 'input_s850/CrA_20150323_850_IR2_ext_HKJypix_col.fits'
    temp = 'input_temp/CrA_20150323-autotemperature.fits'
    surface = 'YSO/surface/Lupus_YSOsurface_r300_units.fits'
    #Aquila
    sig450, sig850 =  0.0247, 0.0018
    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [235.7242,-34.1525,0.1,0.12] #
    
    #Plot maps
    #TEMP PLOT
    #A(s450,s850,temp,sig450,sig850,surface)
    B(s850,temp,sig850,pos)
elif Map == 'Perseus':
    #s450 = 'input_s450/PerseusWest_20150227_450_IR2_ext_HKJypix_col.fits'
    #s850 = 'input_s850/PerseusWest_20150317_850_noco_IR2_HKJypix_col.fits'
    #temp = 'input_temp/PerseusWest_20150317-autotemperature.fits'
    #surface = 'YSO/surface/PerseusWest_YSOsurface_r300_units.fits'
    s450 = 'input_s450/PerseusIC348_20150227_450_IR2_ext_HKJypix_col.fits'
    s850 = 'input_s850/PerseusIC348_20150317_850_noco_IR2_HKJypix_col.fits'
    temp = 'input_temp/PerseusIC348_20150317-autotemperature.fits'
    #surface = 'YSO/surface/PerseusWest_YSOsurface_r300_units.fits'
    #PW
    sig450,sig850 = 0.0293, 0.0050 #IC 348
    #sig450,sig850 = 0.0148,0.0017 #NGC 1333

    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    #pos = [52.2933,31.36639,0.1,0.12] #NGC133
    pos = [56.1417,32.1628,0.1,0.12] #IC348
    #Plot maps
    #TEMP PLOT
    #A(s450,s850,temp,sig450,sig850,surface)
    B(s850,temp,sig850,pos)
elif Map == 'Auriga':
    #LKHA101
    #s450 = 'input_s450/PerseusWest_20150227_450_IR2_ext_HKJypix_col.fits'
    s850 = 'input_s850/Auriga_LKHa101_20150318_850_IR2_noco_HKJypix_col.fits'
    temp = 'input_temp/Auriga_LKHa101_20150318-autotemperature.fits'
    #Main
    #s450 = 'input_s450/PerseusWest_20150227_450_IR2_ext_HKJypix_col.fits'
    #s850 = 'input_s850/AurigaMain_20150318_850_IR2_noco_HKJypix_col.fits'
    #temp = 'input_temp/AurigaMain_20150318-autotemperature.fits'
    #CN
    #s450 = 'input_s450/PerseusIC348_20150227_450_IR2_ext_HKJypix_col.fits'
    #s850 = 'input_s850/Auriga_CN_20150224_850_IR2_ext_HKJypix_col.fits'
    #temp = 'input_temp/Auriga_CN-autotemperature.fits'

    #surface = 'YSO/surface/PerseusWest_YSOsurface_r300_units.fits'

    #sig450,sig850 =0.0345,0.0020 #AC
    sig450,sig850 =0.0115,0.0017 #Amain

    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    #pos = [67.5583,35.2733,0.1,0.12] #AC
    pos = [67.5583,35.2733,0.1,0.12] #Main
    #Plot maps
    #TEMP PLOT
    #A(s450,s850,temp,sig450,sig850,surface)
    B(s850,temp,sig850,pos)
elif Map == 'Ceph':
    #L1228
    #s450 = 'input_s450/PerseusWest_20150227_450_IR2_ext_HKJypix_col.fits'
    s850 = 'input_s850/Cepheus_L1228_20150224_850_IR2_ext_HKJypix_col.fits'
    temp = 'input_temp/Cepheus_L1228-autotemperature.fits'
    #L1251
    #s450 = 'input_s450/PerseusWest_20150227_450_IR2_ext_HKJypix_col.fits'
    #s850 = 'input_s850/Cepheus_L1251_20150224_850_IR2_ext_HKJypix_col.fits'
    #temp = 'input_temp/Cepheus_L1251-autotemperature.fits'
    #South
    #s450 = 'input_s450/PerseusIC348_20150227_450_IR2_ext_HKJypix_col.fits'
    #s850 = 'input_s850/Cepheus_South_20150225_850_IR2_noco_HKJypix_col.fits'
    #temp = 'input_temp/Cepheus_South-autotemperature.fits'

    #Surface map
    #surface = 'YSO/surface/PerseusWest_YSOsurface_r300_units.fits'
    
    #Noise
    #sig450,sig850 = 0.0256,0.0059 #South
    sig450,sig850 = 0.0124,0.0019 #L1228
    #sig450,sig850 = 0.0170,0.0018 #L1251

    #centre - fmt: xc(deg), yc(deg), h(frac), w(frac)
    pos = [315.4042,68.1633,0.1,0.12] #South
    
    #Plot maps
    #A(s450,s850,temp,sig450,sig850,surface)
    B(s850,temp,sig850,pos)

plt.show()
