#20/01/2016

#This script is designed to plot temp and CD maps for the LUPUS region for thesis chapter 8.

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
import time
import astropy.io.fits as pyfits
import os

import ndf2fits
import YSO

########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

####################################################
#maps

s850 = 's850/Lupus_20150318_850_IR2_ext_HKJypix_col.fits'
temp = 'temperature/Lupus_20150318-autotemperature.fits'
proto = 'YSO/surface/Lupus_proto_surface_r300_units.fits'
pms = 'YSO/surface/Lupus_PMS_surface_r300_units.fits'
clumps = 'FW/clumps/Lupus_20150318_clumps.fits'
clumpcat = 'FW/cat/Lupus_20150318_clump_cat.fit'
CD = 'CDmaps/Lupus_CDcomb.fits'

#Noise levels
sig850 = 0.0019

#sub region position coords
pos1 = [235.7242,-34.1525,0.05,0.05] #Central
pos2 = [236.321666667 ,-34.3401111111,0.1,0.1] #East

#sub region figure sizes: ix,iy,dx,dy
plot1 = [0.1,0.1,0.38,0.4]
plot2 = [0.51,0.1,0.38,0.4]
plot3 = [0.1,0.5,0.38,0.4]
plot4 = [0.51,0.5,0.38,0.4]

#ENTER MAP VERSION HERE
version = 'temp'

if version == 'temp':
    #deffine 'big' fig
    big = plt.figure(figsize=(18, 12))
    #deffine 
    Fig1 = aplpy.FITSFigure(temp,figure=big,subplot=plot1)
    Fig2 = aplpy.FITSFigure(temp,figure=big,subplot=plot2)

    #Show colour-scale
    Fig1.show_colorscale(vmax=35.0,vmin=10.0, stretch='linear',cmap='CMRmap')
    Fig2.show_colorscale(vmax=35.0,vmin=10.0, stretch='linear',cmap='CMRmap')

    #PLOT contours
    #SCUBA-2
    Fig1.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))
    Fig2.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))
    #YSO
    levels = (15,45,75,100,150)
    Fig1.show_contour(proto, levels=(levels),colors='g')
    Fig2.show_contour(proto, levels=(levels),colors='g')
    Fig1.show_contour(pms, levels=(levels),colors='r')
    Fig2.show_contour(pms, levels=(levels),colors='r')
    
    #PLOT OB stars - primary stars in cyan
    data = np.loadtxt('OBstars_list.txt')
    RA_OB = data[:,0]	
    DEC_OB = data[:,1]
    Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
    Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)

    data = np.loadtxt('OBstars.txt',dtype='string')
    RA_OB = map(float, data[:,2])	
    DEC_OB = map(float, data[:,3])
    Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
    Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)

    #Recenter
    Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
    Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])

    #Label
    Fig1.add_label(0.15, 0.95, 'Lupus I Central', relative=True)
    Fig2.add_label(0.15, 0.95, 'Lupus I East', relative=True)

    #format axes
    Fig2.hide_yaxis_label()

    Fig1.ticks.show()
    Fig1.ticks.set_color('black')
    Fig1.tick_labels.set_font(size='medium')
    Fig2.ticks.show()
    Fig2.ticks.set_color('black')
    Fig2.tick_labels.set_font(size='medium')
    
    #Colour bar - remove if necessary
    Fig2.add_colorbar()
    Fig2.colorbar.set_axis_label_text('T$_{d}$ (K)') 

    big.canvas.draw()

    #Save fig specific date and time
    date = str(time.strftime("%Y%m%d"))
    #savefig('plots/%s_Lupus_temp.pdf'%(date))


#if version == 'CD':
    #deffine 'big' fig
    #big = plt.figure(figsize=(17, 12))
    #deffine map
    Fig1 = aplpy.FITSFigure(CD,figure=big,subplot=plot3)
    Fig2 = aplpy.FITSFigure(CD,figure=big,subplot=plot4)

    #Show colour-scale
    Fig1.show_colorscale(vmax=3.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
    Fig2.show_colorscale(vmax=3.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')

    #PLOT FW contours
    #OPEN up the region specific clump list
    data = pyfits.open(clumpcat)
    DATA = data[1].data

    for j in DATA:
        I = int(j[0])
        #print I
        ### SETUP temp. file names
        thresh = 'mask/thresh'+str(I)+'.sdf'
        mask = 'mask/mask'+str(I)+'.sdf'
        magic = 'mask/magic'+str(I)+'.sdf'

        ### isolate clump i ###
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumps,thresh,I,I,'bad','bad','> /dev/null')
        #os.system(cmd)
         ### process 2 - remove 0value data
        cmd = '%s/nomagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
        #os.system(cmd)
        ### recreate mask of clump i - mask0
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
        #os.system(cmd)
        
        
        #ndf2fits.ndf2fits(mask)

        maskFITS = 'mask/mask'+str(I)+'.fits'

        Fig1.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
        Fig2.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

        cmd = 'rm %s %s %s'%(thresh,mask,magic)
        #os.system(cmd)

    #Plot YSOs
    YSO.c2dGBS_v2(Fig1)
    YSO.c2dGBS_v2(Fig2)
    
    #PLOT OB stars - primary stars in cyan
    data = np.loadtxt('OBstars_list.txt')
    RA_OB = data[:,0]	
    DEC_OB = data[:,1]
    Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)
    Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)

    data = np.loadtxt('OBstars.txt',dtype='string')
    RA_OB = map(float, data[:,2])	
    DEC_OB = map(float, data[:,3])
    Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
    Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)

    #Recenter
    Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
    Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])

    #Label
    Fig1.add_label(0.15, 0.95, 'Lupus I Central', relative=True)
    Fig2.add_label(0.15, 0.95, 'Lupus I East', relative=True)

    #format axes
    Fig2.hide_yaxis_label()

    Fig1.ticks.show()
    Fig1.ticks.set_color('black')
    Fig1.tick_labels.set_font(size='medium')
    Fig2.ticks.show()
    Fig2.ticks.set_color('black')
    Fig2.tick_labels.set_font(size='medium')
    Fig1.axis_labels.hide_x()
    Fig2.axis_labels.hide_x()
    Fig1.tick_labels.hide_x()
    Fig2.tick_labels.hide_x()
    
    #Colour bar - remove if necessary
    Fig2.add_colorbar()
    Fig2.colorbar.set_axis_label_text('Column Density (H$_{2}$ cm$^{-2}$)') 

    big.canvas.draw()

    #Save fig specific date and time
    date = str(time.strftime("%Y%m%d"))
    savefig('plots/%s_Lupus.pdf'%(date))



plt.show()



