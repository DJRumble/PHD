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

#L1228
s850a = 's850/Cepheus_L1228_20150224_850_IR2_ext_HKJypix_col.fits'
tempa = 'temperature/Cepheus_L1228-autotemperature.fits'
clumpcata = 'FW/cat/Cepheus_L1228_20150224_clump_cat.fit '
clumpsa = 'FW/clumps/Cepheus_L1228_20150224_clumps.sdf'
protoa = 'YSO/surface/Cepheus_L1228_proto_surface_r300_units.fits'
pmsa = 'YSO/surface/Cepheus_L1228_PMS_surface_r300_units.fits'
CDa = 'CDmaps/Cepheus_L1228_CDcomb.fits'
#L1251
s850b = 's850/Cepheus_L1251_20150224_850_IR2_ext_HKJypix_col.fits'
tempb = 'temperature/Cepheus_L1251-autotemperature.fits'
clumpcatb = 'FW/cat/Cepheus_L1251_20150224_clump_cat.fit'
clumpsb = 'FW/clumps/Cepheus_L1251_20150224_clumps.sdf'
protob = 'YSO/surface/Cepheus_L1251_proto_surface_r300_units.fits'
pmsb = 'YSO/surface/Cepheus_L1251_PMS_surface_r300_units.fits'
CDb = 'CDmaps/Cepheus_L1251_CDcomb.fits'


sig850 = 0.0041

pos1 = [314.300833333,77.5938333333,0.15,0.15] #L1228
pos2 = [339.324583333,75.2455416667,0.15,0.15] #L1251B
pos3 = [337.629166667,75.2355833333,0.15,0.15] #L1251A

#subregion figure sizes: ix,iy,dx,dy
plot1 = [0.15,0.05,0.4,0.28]
plot2 = [0.55,0.05,0.4,0.28]
plot3 = [0.15,0.36,0.4,0.28]
plot4 = [0.55,0.36,0.4,0.28]
plot5 = [0.15,0.67,0.4,0.28]
plot6 = [0.55,0.67,0.4,0.28]

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
Fig1 = aplpy.FITSFigure(tempa,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(tempb,figure=big,subplot=plot3)
Fig3 = aplpy.FITSFigure(tempb,figure=big,subplot=plot5)

#Show colour-scale
Fig1.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig2.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig3.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')

#PLOT contours
#SCUBA-2
Fig1.show_contour(s850a, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig2.show_contour(s850b, levels=(sig850,1000), linewidth=3, colors=('k'))
Fig3.show_contour(s850b, levels=(sig850,1000), linewidth=3, colors=('k'))

levels = (15,45,75,100,150)
Fig1.show_contour(protoa, levels=(levels),colors='g')
Fig2.show_contour(protob, levels=(levels),colors='g')
Fig3.show_contour(protob, levels=(levels),colors='g')
Fig1.show_contour(pmsa, levels=(levels),colors='r')
Fig2.show_contour(pmsb, levels=(levels),colors='r')
Fig3.show_contour(pmsb, levels=(levels),colors='r')
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig3.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig3.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])
Fig3.recenter(pos3[0],pos3[1],pos3[2],pos3[3])

#Label
Fig3.add_label(0.15, 0.9, 'L1251A', relative=True)
Fig2.add_label(0.15, 0.9, 'L1251B', relative=True)
Fig1.add_label(0.15, 0.9, 'L1228', relative=True)

#Fig1.add_label(0.45, 0.65, 'LkHa 101', relative=True)

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
    
#Colour bar - remove if necessary
#Fig3.add_colorbar()
#Fig3.colorbar.set_location('top')
#Fig3.colorbar.set_axis_label_text('T$_{d}$ (K)') 
#box = [0.1,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')

#Add Scale Bar (1pc)
Fig1.add_scalebar(516 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
Fig2.add_scalebar(344 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
Fig3.add_scalebar(344 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))

#################
#### CD MAPS ####
#################

Fig1 = aplpy.FITSFigure(CDa,figure=big,subplot=plot2)
Fig2 = aplpy.FITSFigure(CDb,figure=big,subplot=plot4)
Fig3 = aplpy.FITSFigure(CDb,figure=big,subplot=plot6)

#Show colour-scale
Fig1.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig2.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig3.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')

#PLOT FW contours
#OPEN up the region specific clump list
data = pyfits.open(clumpcata)
DATA = data[1].data

for j in DATA:
    I = int(j[0])
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_Cephb_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_Cephb_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumpsa,thresh,I,I,'bad','bad','> /dev/null')
    os.system(cmd)
    ### process 2 - remove 0value data
    cmd = '%s/nomagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
    os.system(cmd)
    ### recreate mask of clump i - mask0
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
    os.system(cmd)
        
    ndf2fits.ndf2fits(mask)

    Fig1.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig2.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig3.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    os.system(cmd)

#PLOT FW contours
#OPEN up the region specific clump list
data = pyfits.open(clumpcatb)
DATA = data[1].data

for j in DATA:
    I = int(j[0])
    print str(I)+'/'+str(len(DATA))
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_Cephc_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_Cephc_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumpsb,thresh,I,I,'bad','bad','> /dev/null')
    os.system(cmd)
    ### process 2 - remove 0value data
    cmd = '%s/nomagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
    os.system(cmd)
    ### recreate mask of clump i - mask0
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
    os.system(cmd)
        
    ndf2fits.ndf2fits(mask)

    Fig1.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig2.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))
    Fig3.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    os.system(cmd)


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
Fig3.add_label(0.15, 0.9, 'L1251A', relative=True)
Fig2.add_label(0.15, 0.9, 'L1251B', relative=True)
Fig1.add_label(0.15, 0.9, 'L1228', relative=True)

#Fig1.add_label(0.45, 0.65, 'LkHa 101', relative=True)

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

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()

#Colour bar - remove if necessary
#Fig3.add_colorbar()
#Fig3.colorbar.set_location('top')
#Fig3.colorbar.set_axis_label_text('Column Density (H$_{2}$ cm$^{-2}$)') 
#box = [0.5,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')

big.canvas.draw()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
print 'saving map as: '+str('plots/%s_CephB.pdf'%(date))
savefig('plots/%s_CephB.pdf'%(date))


plt.show()

