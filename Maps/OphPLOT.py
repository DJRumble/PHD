#20/01/2016

#This script is designed to plot temp and CD maps for thesis chapter 8.

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

s850 = 's850/OphSco_Main_20150325_850_IR2_noco_HKJypix_col.fits'
temp = 'temperature_maps/OphSco_Main-autotemperature.fits'
clumpcat1 = 'FW/cat/OphSco_L1688_20150318_clump_cat.fit'
clumpcat2 = 'FW/cat/OphSco_L1689_20150318_clump_cat.fit'
clumps1 = 'FW/clumps/OphSco_L1688_20150318_clumps.sdf'
clumps2 = 'FW/clumps/OphSco_L1689_20150318_clumps.sdf'
proto = 'YSO/surface/OphSco_Main_proto_surface_r300_units.fits'
pms = 'YSO/surface/OphSco_Main_PMS_surface_r300_units.fits'
CD = 'CDmaps/OphSco_L1688_CDcomb.fits'

sig850 = 0.002

pos1 = [246.617083333,-24.34925,0.12,0.12] 
pos2 = [246.894583333,-24.5659166667,0.18,0.18] 
pos3 = [248.114166667,-24.4758055556,0.07,0.07] 

#subregion figure sizes: ix,iy,dx,dy
plot1 = [0.10,0.05,0.35,0.28]
plot2 = [0.55,0.05,0.35,0.28]
plot3 = [0.10,0.36,0.35,0.28]
plot4 = [0.55,0.36,0.35,0.28]
plot5 = [0.10,0.67,0.35,0.28]
plot6 = [0.55,0.67,0.35,0.28]

#Deffine OB star data
data1 = np.loadtxt('OBstars_list.txt')
RA_OB = data1[:,0]	
DEC_OB = data1[:,1]

data2 = np.loadtxt('OBstars.txt',dtype='string')
RA_OBx = map(float, data2[:,2])	
DEC_OBx = map(float, data2[:,3])

### Deffine large map ###

#deffine 'big' fig
big = plt.figure(figsize=(12, 15))

################
### TEMP MAP ###
################

#deffine 
Fig1 = aplpy.FITSFigure(temp,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(temp,figure=big,subplot=plot3)
Fig3 = aplpy.FITSFigure(temp,figure=big,subplot=plot5)

#Show colour-scale
Fig1.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig2.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig3.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')

#PLOT contours
#SCUBA-2
Fig1.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig2.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig3.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))

levels = (15,45,75,100,150)
#Fig1.show_contour(proto, levels=(levels),colors='g')
#Fig2.show_contour(proto, levels=(levels),colors='g')
#Fig3.show_contour(proto, levels=(levels),colors='g')
#Fig1.show_contour(pms, levels=(levels),colors='r')
#Fig2.show_contour(pms, levels=(levels),colors='r')
#Fig3.show_contour(pms, levels=(levels),colors='r')

#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
YSO.c2dGBS_v2(Fig3)
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig3.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig3.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])
Fig3.recenter(pos3[0],pos3[1],pos3[2],pos3[3])

#Label
Fig3.add_label(0.15, 0.9, 'L1689', relative=True)
Fig2.add_label(0.22, 0.9, 'L1688 B/C/E/F', relative=True)
Fig1.add_label(0.15, 0.9, 'L1688 A', relative=True)

#Fig3.add_label(0.75, 0.25, 'IRAS 18352-0148', relative=True)

#Format axis
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')

Fig3.ticks.show()
Fig3.ticks.set_color('black')
Fig3.tick_labels.set_font(size='small')

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()

Fig1.tick_labels.set_xformat('hh:mm:ss')
Fig2.tick_labels.set_xformat('hh:mm:ss')
Fig3.tick_labels.set_xformat('hh:mm:ss')
    
#Colour bar - remove if necessary
Fig3.add_colorbar()
#Fig3.colorbar.set_location('top')
Fig3.colorbar.set_axis_label_text('T$_{d}$ (K)') 
#box = [0.1,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')
Fig1.add_colorbar()
Fig1.colorbar.hide()
Fig2.add_colorbar()
Fig2.colorbar.hide()

#Add Scale Bar (1pc)
Fig1.add_scalebar(371 * u.arcsecond,label=('0.25 pc'),corner=('top right'))
Fig2.add_scalebar(371 * u.arcsecond,label=('0.25 pc'),corner=('bottom left'))
Fig3.add_scalebar(371 * u.arcsecond,label=('0.25 pc'),corner=('bottom left'))

#################
#### CD MAPS ####
#################

Fig1 = aplpy.FITSFigure(CD,figure=big,subplot=plot2)
Fig2 = aplpy.FITSFigure(CD,figure=big,subplot=plot4)
Fig3 = aplpy.FITSFigure(CD,figure=big,subplot=plot6)

#Show colour-scale
Fig1.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig2.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig3.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')

#PLOT FW contours
#OPEN up the region specific clump list
#Perseus Oph 1688
data = pyfits.open(clumpcat1)
DATA = data[1].data
#'''
for j in DATA:
    I = int(j[0])
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_1688_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_1688_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumps1,thresh,I,I,'bad','bad','> /dev/null')
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
    Fig3.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    #os.system(cmd)

#Oph 1689
data = pyfits.open(clumpcat2)
DATA = data[1].data
#'''
for j in DATA:
    I = int(j[0])
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_1689_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_1689_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumps2,thresh,I,I,'bad','bad','> /dev/null')
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
    Fig3.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    
    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    #os.system(cmd)
#'''
#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
YSO.c2dGBS_v2(Fig3)
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig3.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig3.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])
Fig3.recenter(pos3[0],pos3[1],pos3[2],pos3[3])

#Label
Fig3.add_label(0.15, 0.9, 'L1689', relative=True)
Fig2.add_label(0.22, 0.9, 'L1688 B/C/E/F', relative=True)
Fig1.add_label(0.15, 0.9, 'L1688 A', relative=True)

#format axes
Fig1.ticks.show()
Fig1.ticks.set_color('black')
Fig1.tick_labels.set_font(size='small')
Fig2.ticks.show()
Fig2.ticks.set_color('black')
Fig2.tick_labels.set_font(size='small')
Fig3.ticks.show()
Fig3.ticks.set_color('black')
Fig3.tick_labels.set_font(size='small')

Fig1.axis_labels.hide_y()
Fig2.axis_labels.hide_y()
Fig3.axis_labels.hide_y()
Fig1.tick_labels.hide_y()
Fig2.tick_labels.hide_y()
Fig3.tick_labels.hide_y()
Fig1.tick_labels.set_xformat('hh:mm:ss')
Fig2.tick_labels.set_xformat('hh:mm:ss')
Fig3.tick_labels.set_xformat('hh:mm:ss')

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()

#Colour bar - remove if necessary
Fig3.add_colorbar()
#Fig3.colorbar.set_location('top')
Fig3.colorbar.set_axis_label_text('Column Density (H$_{2}$ cm$^{-2}$)') 
#box = [0.5,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')
Fig1.add_colorbar()
Fig1.colorbar.hide()
Fig2.add_colorbar()
Fig2.colorbar.hide()

big.canvas.draw()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
print 'saving map as: '+str('plots/%s_Oph.pdf'%(date))
savefig('plots/%s_Oph.pdf'%(date))

plt.show()
