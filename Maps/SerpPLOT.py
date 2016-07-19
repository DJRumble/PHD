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

####################################################
#maps

#Main
s850a = 's850/SerpensMain_20150326_850_IR2_noco_HKJypix_col.fits'
tempa = 'temperature/SerpensMain_20150326-autotemperature.fits'
clumpcata = 'FW/cat/SerpensMain_20150326_clump_cat.fit'
clumpsa = 'FW/clumps/SerpensMain_20150326_clumps.sdf'
protoa = 'YSO/surface/SerpensMain_YSOsurface_r300_units.fits' ##??
pmsa = 'YSO/surface/SerpensMain_YSOsurface_r300_units.fits' ##??
CDa = 'CDmaps/Serpens_Main_CDcomb.fits'
#MWC 297
s850b = 's850/SerpensMWC297_20140414_IR1_s850_freefree+jet_col.fits'
tempb = 'temperature/SerpensMWC297_20140414-autotemperature.fits'
clumpcatb = 'FW/cat/SerpensMWC297_20140414_IR1_sclump_cat.fit'
clumpsb = 'FW/clumps/SerpensMWC297_20140414_IR1_sclumps.sdf'
protob = 'YSO/surface/SerpensMWC297_proto_surface_r300_units.fits'
pmsb = 'YSO/surface/SerpensMWC297_PMS_surface_r300_units.fits'
CDb = 'CDmaps/SerpensMWC297_CDcomb.fits'
#North
s850c = 's850/SerpensN_20141219_s850_IR2extmask_s2_cal_JypixJH_col.fits'
tempc = 'temperature/SerpensN_20141219-autotemperature.fits'
clumpcatc = 'FW/cat/SerpensN_20141219_sclump_cat_v2.fit'
clumpsc = 'FW/clumps/SerpensN_20141219_sclumps_v2.sdf'
protoc = 'YSO/surface/SerpensN_proto_surface_r300_units.fits'
pmsc = 'YSO/surface/SerpensN_PMS_surface_r300_units.fits'
CDc = 'CDmaps/SerpensN_CDcomb.fits'

#Noise levels
sig850a = 0.0031
sig850b = 0.0111
sig850c = 0.0027

#sub region position coords
pos4 = [277.47625,1.23358333333,0.13,0.13] #Main
pos3 = [277.322083333,0.543861111111,0.12,0.12] #NH3
pos2 = [277.079166667,-3.74194444444,0.15,0.15] #MWC 297
pos1 = [279.759583333,0.465833333333,0.15,0.15] #North

#subregion figure sizes: ix,iy,dx,dy
plot1 = [0.15,0.04,0.4,0.20]
plot2 = [0.55,0.04,0.4,0.20]
plot3 = [0.15,0.28,0.4,0.20]
plot4 = [0.55,0.28,0.4,0.20]
plot5 = [0.15,0.52,0.4,0.20]
plot6 = [0.55,0.52,0.4,0.20]
plot7 = [0.15,0.76,0.4,0.20]
plot8 = [0.55,0.76,0.4,0.20]

#Deffine OB star data
data1 = np.loadtxt('OBstars_list.txt')
RA_OB = data1[:,0]	
DEC_OB = data1[:,1]

data2 = np.loadtxt('OBstars.txt',dtype='string')
RA_OBx = map(float, data2[:,2])	
DEC_OBx = map(float, data2[:,3])

### Deffine large map ###

#deffine 'big' fig
big = plt.figure(figsize=(10,16))

################
### TEMP MAP ###
################

#deffine 
Fig1 = aplpy.FITSFigure(tempc,figure=big,subplot=plot1)
Fig2 = aplpy.FITSFigure(tempb,figure=big,subplot=plot3)
Fig3 = aplpy.FITSFigure(tempa,figure=big,subplot=plot5)
Fig4 = aplpy.FITSFigure(tempa,figure=big,subplot=plot7)

#Show colour-scale
Fig1.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig2.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig3.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')
Fig4.show_colorscale(vmax=40.0,vmin=10.0, stretch='linear',cmap='CMRmap')


#PLOT contours
#SCUBA-2
Fig1.show_contour(s850c, levels=(sig850c,1000), linewidth=3, colors=('k'))
Fig2.show_contour(s850b, levels=(sig850b,1000), linewidth=3, colors=('k'))
Fig3.show_contour(s850a, levels=(sig850a,1000), linewidth=3, colors=('k'))
Fig4.show_contour(s850a, levels=(sig850a,1000), linewidth=3, colors=('k'))


levels = (15,45,75,100,150)
#Fig1.show_contour(protoc, levels=(levels),colors='g')
#Fig2.show_contour(protob, levels=(levels),colors='g')
#Fig3.show_contour(protoa, levels=(levels),colors='g')
#Fig4.show_contour(protoa, levels=(levels),colors='g')
#Fig1.show_contour(pmsc, levels=(levels),colors='r')
#Fig2.show_contour(pmsb, levels=(levels),colors='r')
#Fig3.show_contour(pmsa, levels=(levels),colors='r')
#Fig4.show_contour(pmsa, levels=(levels),colors='r')
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig3.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)
Fig4.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=100, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig3.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)
Fig4.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=100, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])
Fig3.recenter(pos3[0],pos3[1],pos3[2],pos3[3])
Fig4.recenter(pos4[0],pos4[1],pos4[2],pos4[3])

#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
YSO.c2dGBS_v2(Fig3)
YSO.c2dGBS_v2(Fig4)

#Label
Fig4.add_label(0.2, 0.9, 'Main', relative=True)
Fig3.add_label(0.2, 0.9, 'NH$_{3}$', relative=True)
Fig2.add_label(0.2, 0.9, 'MWC 297', relative=True)
Fig1.add_label(0.2, 0.9, 'North', relative=True)

Fig2.add_label(0.85, 0.05, 'MWC 297', relative=True)

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
Fig4.ticks.show()
Fig4.ticks.set_color('black')
Fig4.tick_labels.set_font(size='small')

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()
Fig4.axis_labels.hide_x()
    
#Colour bar - remove if necessary
#Fig4.add_colorbar()
#Fig4.colorbar.set_location('top')
#Fig4.colorbar.set_axis_label_text('T$_{d}$ (K)') 
#box = [0.1,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')

#Add Scale Bar (1pc at 450 pcs)
Fig1.add_scalebar(412 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
Fig2.add_scalebar(412 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
Fig3.add_scalebar(240 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))
Fig4.add_scalebar(240 * u.arcsecond,label=('0.5 pc'),corner=('bottom left'))

#################
#### CD MAPS ####
#################

Fig1 = aplpy.FITSFigure(CDc,figure=big,subplot=plot2)
Fig2 = aplpy.FITSFigure(CDb,figure=big,subplot=plot4)
Fig3 = aplpy.FITSFigure(CDa,figure=big,subplot=plot6)
Fig4 = aplpy.FITSFigure(CDa,figure=big,subplot=plot8)

#Show colour-scale
Fig1.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig2.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig3.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')
Fig4.show_colorscale(vmax=8.5E22,vmin=2.5E21, stretch='linear',cmap='afmhot_r')


#PLOT FW contours
#OPEN up the region specific clump list
data = pyfits.open(clumpcatc)
DATA = data[1].data

for j in DATA:
    I = int(j[0])
    print I
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_SerpN_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_SerpN_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumpsc,thresh,I,I,'bad','bad','> /dev/null')
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

#PLOT FW contours
#OPEN up the region specific clump list
data = pyfits.open(clumpcatb)
DATA = data[1].data

for j in DATA:
    I = int(j[0])
    print I
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_MWC297_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_MWC297_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumpsb,thresh,I,I,'bad','bad','> /dev/null')
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

#PLOT FW contours
#OPEN up the region specific clump list
data = pyfits.open(clumpcata)
DATA = data[1].data

for j in DATA:
    I = int(j[0])
    print I
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask_Serp_'+str(I)+'.sdf'
    maskFITS = 'mask/mask_Serp_'+str(I)+'.fits'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumpsa,thresh,I,I,'bad','bad','> /dev/null')
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

#Plot YSOs
YSO.c2dGBS_v2(Fig1)
YSO.c2dGBS_v2(Fig2)
YSO.c2dGBS_v2(Fig3)
YSO.c2dGBS_v2(Fig4)
    
#PLOT OB stars - primary stars in cyan
Fig1.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig3.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)
Fig4.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='none', marker='*',s=250, alpha=1)

Fig1.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig2.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig3.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)
Fig4.show_markers(RA_OBx,DEC_OBx, edgecolor='k', facecolor='cyan', marker='*',s=250, alpha=1)

#Recenter
Fig1.recenter(pos1[0],pos1[1],pos1[2],pos1[3])
Fig2.recenter(pos2[0],pos2[1],pos2[2],pos2[3])
Fig3.recenter(pos3[0],pos3[1],pos3[2],pos3[3])
Fig4.recenter(pos4[0],pos4[1],pos4[2],pos4[3])

#Label
Fig4.add_label(0.2, 0.9, 'Main', relative=True)
Fig3.add_label(0.2, 0.9, 'NH$_{3}$', relative=True)
Fig2.add_label(0.2, 0.9, 'MWC 297', relative=True)
Fig1.add_label(0.2, 0.9, 'North', relative=True)

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
Fig4.ticks.show()
Fig4.ticks.set_color('black')
Fig4.tick_labels.set_font(size='small')

Fig1.axis_labels.hide_y()
Fig2.axis_labels.hide_y()
Fig3.axis_labels.hide_y()
Fig4.axis_labels.hide_y()
Fig1.tick_labels.hide_y()
Fig2.tick_labels.hide_y()
Fig3.tick_labels.hide_y()
Fig4.tick_labels.hide_y()

Fig2.axis_labels.hide_x()
Fig3.axis_labels.hide_x()
Fig4.axis_labels.hide_x()

#Colour bar - remove if necessary
#Fig4.add_colorbar()
#Fig4.colorbar.set_location('top')
#Fig4.colorbar.set_axis_label_text('Column Density (H$_{2}$ cm$^{-2}$)') 
#box = [0.5,0.95,0.4,0.05]
#Fig3.colorbar.set_box(box, box_orientation='horizontal')

big.canvas.draw()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
print 'saving map as: '+str('plots/%s_Serp.pdf'%(date))
savefig('plots/%s_Serp.pdf'%(date))

plt.show()


