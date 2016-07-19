#20151209

#A mapping script to plot FW contours over a flux map.

import aplpy
import matplotlib.pyplot as plt
import ndf2fits
from astropy.io import fits
import os
import YSO
import time
from pylab import *

kapdir = '/star/bin/kappa'

###########
#PLOT-flux#
###########

s850 = 'input_s850/IC5146_20150225_850_IR2_ext_HKJypix_col.fits'
s450 = 'input_s450/IC5146_20150225_450_IR2_ext_HKJypix_col.fits'
temp = 'input_temp/IC5146_20150225-autotemperature.fits'
clumpcat = 'FW/cat/IC5146_20150225_clump_cat.fit'
clumps = 'FW/clumps/IC5146_20150225_clumps.sdf'
surface = 'YSO/surface/IC5146_YSOsurface_r300_units.fits'

pos = [328.368,47.2667,0.25,0.35]

sig850 = 0.0016
sig450 = 0.0181

big = plt.figure(figsize=(25, 6))

Fig = aplpy.FITSFigure(s850,figure=big,subplot=[0.1,0.1,0.375,0.8])
Fig.show_grayscale(vmax=10*sig850,vmin=(1./5.)*sig850,invert=True, stretch='sqrt')

###############
#PLOT-contours#
###############

#OPEN up the region specific clump list
data = fits.open(clumpcat)
DATA = data[1].data

#PLOT ysos
YSO.c2dGBS(Fig)

#OB star position
RA_OB = 328.3708	
DEC_OB = 47.2667
Fig.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)


Fig.show_contour(s850, levels=(sig850,300*sig850), linewidth=3, colors=('black','black'))

###########
#PLOT-misc#
###########

Fig.ticks.show()
Fig.ticks.set_color('black')
Fig.tick_labels.set_font(size='medium')
Fig.add_colorbar()
Fig.colorbar.set_axis_label_text('S$_{850\mu m}}$ (Jy/pix)') 
Fig.recenter(pos[0],pos[1], height=pos[2] ,width=pos[3])


#temp zoom
Fig_temp = aplpy.FITSFigure(temp,figure=big,subplot=[0.525,0.1,0.375,0.8])
Fig_temp.show_colorscale(vmax=35.0,vmin=10.0, stretch='linear')

#OB star position
RA_OB = 328.3708	
DEC_OB = 47.2667
Fig_temp.show_markers(RA_OB,DEC_OB, edgecolor='k', facecolor='yellow', marker='*',s=100, alpha=1)

for j in DATA:
    I = int(j[0])
    print I
    ### SETUP temp. file names
    thresh = 'mask/thresh'+str(I)+'.sdf'
    mask = 'mask/mask'+str(I)+'.sdf'
    magic = 'mask/magic'+str(I)+'.sdf'

    ### isolate clump i ###
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%s %s'%(kapdir,clumps,thresh,I,I,'bad','bad','> /dev/null')
    #os.system(cmd)
    ###process 2 - remove 0value data
    cmd = '%s/nomagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
    #os.system(cmd)
    ### recreate mask of clump i - mask0
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
    #os.system(cmd)
    #if I <= 13:
        
    #ndf2fits.ndf2fits(mask)

    maskFITS = 'mask/mask'+str(I)+'.fits'

    Fig_temp.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s %s'%(thresh,mask,magic,maskFITS)
    #os.system(cmd)


Fig_temp.show_contour(surface, levels=(20,60,100), linewidth=3, colors=('green'),linestyle='--')

Fig_temp.ticks.show()
Fig_temp.ticks.set_color('black')
Fig_temp.tick_labels.set_font(size='medium')
Fig_temp.hide_ytick_labels()
Fig_temp.hide_yaxis_label()
Fig_temp.add_colorbar()
Fig_temp.colorbar.set_axis_label_text('T$_{d}$ (K)') 
Fig_temp.recenter(pos[0],pos[1], height=pos[2] ,width=pos[3])

big.canvas.draw()

plt.show()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_IC5146_HD46.pdf'%(date))
