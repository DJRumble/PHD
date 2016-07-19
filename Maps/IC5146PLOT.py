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

s850 = 's850/IC5146_20150225_850_IR2_ext_HKJypix_col.fits'
temp = 'temperature/IC5146_20150225-autotemperature.fits'
clumpcat = 'FW/cat/IC5146_20150225_clump_cat.fit'
clumps = 'FW/clumps/IC5146_20150225_clumps.sdf'
proto = 'YSO/surface/IC5146_proto_surface_r300_units.fits'
pms = 'YSO/surface/IC5146_PMS_surface_r300_units.fits'
CD = 'CDmaps/IC5146_CDcomb.fits'

sig850 = 0.0026

pos1 = [328.368,47.2667,0.10,0.18] #EAST
pos2 = [326.797083333,47.621,0.2,0.3] #WEST

plot1 = [0.1,0.1,0.35,0.4]
plot2 = [0.5,0.1,0.3,0.4]
plot3 = [0.1,0.5,0.35,0.4]
plot4 = [0.5,0.5,0.3,0.4]
#Deffine OB star data
data1 = np.loadtxt('OBstars_list.txt')
RA_OB = data1[:,0]	
DEC_OB = data1[:,1]

data2 = np.loadtxt('OBstars.txt',dtype='string')
RA_OBx = map(float, data2[:,2])	
DEC_OBx = map(float, data2[:,3])

### Deffine large map ###

#deffine 'big' fig
big = plt.figure(figsize=(22, 6))


################
### TEMP MAP ###
################

#deffine 
Fig1 = aplpy.FITSFigure(temp,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(temp,figure=big,subplot=plot2)

Fig1.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig2.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')

Fig1.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig2.show_contour(s850, levels=(sig850,1000), linewidth=3, colors=('k'))

levels = (15,45,75,100,150)

Fig1.show_contour(proto, levels=(levels),colors='g')
Fig2.show_contour(proto, levels=(levels),colors='g')
Fig1.show_contour(pms, levels=(levels),colors='r')
Fig2.show_contour(pms, levels=(levels),colors='r')

#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])

#Label
Fig2.add_label(0.15, 0.9, 'IC5146 East', relative=True)
Fig1.add_label(0.15, 0.9, 'IC5146 West', relative=True)

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

#Add Scale Bar (1pc at 450 pcs)
Fig1.add_scalebar(217 * u.arcsecond,label=('1 pc'),corner=('bottom left'))
Fig2.add_scalebar(217 * u.arcsecond,label=('1 pc'),corner=('bottom left'))

################
### CD MAP ###
################

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
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_IC5146'+str(I)+'.sdf'
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

    maskFITS = 'mask/mask_IC5146'+str(I)+'.fits'

    Fig1.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig2.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    #os.system(cmd)

#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)


#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])

#Label
#Label
Fig2.add_label(0.15, 0.9, 'IC5146 East', relative=True)
Fig1.add_label(0.15, 0.9, 'IC5146 West', relative=True)

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
print 'saving map as: '+str('plots/%s_IC5146.pdf'%(date))
savefig('plots/%s_IC5146.pdf'%(date))

plt.show()
