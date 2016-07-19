#20/01/2016

#This script is designed to plot temp and CD maps for the LKHa-101 region of Auriga for thesis chapter 8.

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

###########
#PLOT-flux#
###########

s850 = 's850/CrA_20150323_850_IR2_ext_HKJypix_col.fits'
temp = 'temperature_maps/CrA_20150323-autotemperature.fits'
clumpcat = 'FW/cat/CrA_20150323_clump_cat.fit'
clumps = 'FW/clumps/CrA_20150323_clumps.sdf'
proto = 'YSO/surface/CrA_proto_surface_r300_units.fits'
pms = 'YSO/surface/CrA_PMS_surface_r300_units.fits'
CD = 'CDmaps/CrA_CDcomb.fits'

s850_b = 's850/PerseusNGC1333_20150317_850_noco_IR2_HKJypix_col.fits'
temp_b = 'temperature_maps/PerseusWest_20150317-autotemperature.fits'
clumpcat_b = 'FW/cat/PerseusNGC1333_20150317_clump_cat.fit'
clumps_b = 'FW/clumps/PerseusNGC1333_20150317_clumps.sdf'
CD_b = 'CDmaps/PerseusNGC1333_CDcomb.fits'

sig850_b = 0.0033 #NGC1333
sig850_a = 0.0030 #CrA

pos1 = [52.2683333333,31.3017222222,0.12,0.12] #Per NGC 1333
pos2 = [285.470833333,-36.9565277778,0.05,0.05] #CrA

#subregion figure sizes: ix,iy,dx,dy
plot1 = [0.1,0.1,0.35,0.4]
plot2 = [0.55,0.1,0.35,0.4]
plot3 = [0.1,0.55,0.35,0.4]
plot4 = [0.55,0.55,0.35,0.4]

#Deffine OB star data
data1 = np.loadtxt('OBstars_list.txt')
RA_OB = data1[:,0]	
DEC_OB = data1[:,1]

data2 = np.loadtxt('OBstars.txt',dtype='string')
RA_OBx = map(float, data2[:,2])	
DEC_OBx = map(float, data2[:,3])

### Deffine large map ###

#deffine 'big' fig
big = plt.figure(figsize=(10, 15))

################
### TEMP MAP ###
################
#deffine 
Fig1 = aplpy.FITSFigure(temp_b,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(temp,figure=big,subplot=plot3)

#Show colour-scale
Fig1.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig2.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')

#PLOT contours
#SCUBA-2
Fig1.show_contour(s850_b, levels=(sig850_b,1000), linewidth=3, colors=('k'))
Fig2.show_contour(s850, levels=(sig850_a,1000), linewidth=3, colors=('k'))

#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])

#Label
Fig1.add_label(0.22, 0.9, 'NGC 1333', relative=True)
Fig2.add_label(0.22, 0.9, 'CrA', relative=True)

#Fig3.add_label(0.75, 0.25, 'IRAS 18352-0148', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')

Fig2.axis_labels.hide_x()
#Fig2.tick_labels.hide_x()

Fig1.tick_labels.set_xformat('hh:mm:ss')
Fig2.tick_labels.set_xformat('hh:mm:ss')
    
#Colour bar - remove if necessary
Fig2.add_colorbar()
#Fig3.colorbar.set_location('top')
Fig2.colorbar.set_axis_label_text('T$_{d}$ (K)') 
#box = [0.1,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')
Fig1.add_colorbar()
Fig1.colorbar.hide()

#Add Scale Bar (1pc)
Fig1.add_scalebar(219 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
Fig2.add_scalebar(248 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
#################
#### CD MAPS ####
#################

Fig1 = aplpy.FITSFigure(CD_b,figure=big,subplot=plot2)
Fig2 = aplpy.FITSFigure(CD,figure=big,subplot=plot4)

#Show colour-scale
Fig1.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig2.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')

#PLOT FW contours
#OPEN up the region specific clump list
#Perseus NGC 1333
data = pyfits.open(clumpcat)
DATA = data[1].data
#'''
for j in DATA:
    I = int(j[0])
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_NGC1333_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_NGC1333_'+str(I)+'.fits'
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

    Fig1.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig2.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    #os.system(cmd)

#CrA
data = pyfits.open(clumpcat_b)
DATA = data[1].data
#'''
for j in DATA:
    I = int(j[0])
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_CrA_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_CrA_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumps_b,thresh,I,I,'bad','bad','> /dev/null')
    #os.system(cmd)
    ### process 2 - remove 0value data
    cmd = '%s/nomagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
    #os.system(cmd)
    ### recreate mask of clump i - mask0
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
    #os.system(cmd)
        
    #ndf2fits.ndf2fits(mask)

    Fig1.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig2.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    #os.system(cmd)
#'''
#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])

#Label
Fig1.add_label(0.22, 0.9, 'NGC 1333', relative=True)
Fig2.add_label(0.22, 0.9, 'CrA', relative=True)
#Fig3.add_label(0.75, 0.25, 'IRAS 18352-0148', relative=True)

#Format axis
#format axes
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')
#Fig1.axis_labels.hide_x()
Fig2.axis_labels.hide_x()
#Fig1.tick_labels.hide_x()
#Fig2.tick_labels.hide_x()

Fig1.axis_labels.hide_y()
Fig2.axis_labels.hide_y()
Fig1.tick_labels.hide_y()
Fig2.tick_labels.hide_y()

Fig1.tick_labels.set_xformat('hh:mm:ss')
Fig2.tick_labels.set_xformat('hh:mm:ss')
    
#Colour bar - remove if necessary
Fig2.add_colorbar()
Fig2.colorbar.set_axis_label_text('Column Density (H$_{2}$ cm$^{-2}$)') 
Fig1.add_colorbar()
Fig1.colorbar.hide()

big.canvas.draw()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
print 'saving map as: '+str('plots/%s_CrA+NGC1333.pdf'%(date))
savefig('plots/%s_CrA+NGC1333.pdf'%(date))


plt.show()

