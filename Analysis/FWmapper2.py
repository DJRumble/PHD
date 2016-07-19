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

s850 = 'input_s850/Serpens_Main_20150326_850_IR2_noco_HKJypix_col.fits'
#s450 = 'input_s450/IC5146_20150225_450_IR2_ext_HKJypix_col.fits'
#temp = 'input_temp/IC5146_20150225-autotemperature.fits'
clumpcat = 'FW/cat/Serpens_Main_20150326_clump_cat.fit '
clumps = 'FW/clumps/Serpens_Main_20150326_clumps.sdf'
#surface = 'YSO/surface/IC5146_YSOsurface_r300_units.fits'

pos = [277.47625,1.23358333333,0.25,0.25]

sig850 = 0.0031
sig450 = 0.0218


###############
#PLOT-contours#
###############

#OPEN up the region specific clump list
data = fits.open(clumpcat)
DATA = data[1].data


###########
#PLOT-misc#
###########

Fig = aplpy.FITSFigure(s850)
Fig.show_grayscale(vmax=15*sig850,vmin=(1./5.)*sig850,invert=True, stretch='sqrt')

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

    Fig.show_contour(maskFITS, levels=(0.5,1000), linewidth=3, colors=('black'))

    cmd = 'rm %s %s %s'%(thresh,mask,magic)
    #os.system(cmd)

Fig.ticks.show()
Fig.ticks.set_color('black')
Fig.tick_labels.set_font(size='medium')
Fig.add_colorbar()
Fig.colorbar.set_axis_label_text('S$_{850\mu m}}$ (Jy/pix)') 
Fig.recenter(pos[0],pos[1], height=pos[2] ,width=pos[3])

plt.show()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_Main.pdf'%(date))
