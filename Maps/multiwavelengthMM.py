#20150220 - DJR, UoE
#Multiwavelenght map maker

#####################################################
#import maths, ploting and astrophysical packages

import numpy 
import aplpy
import matplotlib.pyplot as plt
import atpy
import c2dtable
from pylab import *

####################################################
#Plotting functions 

def MAP(i,map,Xo,Yo,X,Y,max,min,label):
    #4 pannel image of W40 at different wavelengths

   #deffine images in APLPY
    Fig = aplpy.FITSFigure(map,figure=big,subplot=[Xo,Yo,X,Y])

    #OB markers
    OBmarkers = np.loadtxt('/data/damian/data/aquila/OBobjects.txt')
    ra = OBmarkers[:,0]
    dec = OBmarkers[:,1]
    Fig.show_markers(ra,dec, layer='mwc297', edgecolor='k', facecolor='y', marker='o',s=50, alpha=1)

    #Contours (of 850um)
    s850 = 'aquila/serpens_south/IR2/Aquila_extS2_850-4-30am.fits'
    Fig.show_contour(s850, levels=(0.01,10), colors=('black'))
    
    Fig.add_label(0.25, 0.95,label,size='medium', relative=True)
    Fig.tick_labels.set_font(size='small')
    Fig.axis_labels.set_font(size='small', weight='medium', \
                         stretch='normal', family='sans-serif', \
                         style='normal', variant='normal')
    Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)
    Fig.ticks.show()


    if (i == 1) or (i == 3):
        Fig.axis_labels.hide_y()

    if (i == 2) or (i == 3):
        Fig.axis_labels.hide_x()


    if (i == 1) or (i == 2):
        Fig.show_grayscale(vmax=max,vmin=min, stretch='sqrt',invert='TRUE')
    else:
        Fig.show_colorscale(vmax=max,vmin=min, stretch='linear')
    
    Fig.add_colorbar()
    Fig.colorbar.set_font(size='small', weight='medium', \
                      stretch='normal', family='sans-serif', \
                      style='normal', variant='normal')

    return




####################################################
#maps

#24uM SPITZER
mwc297_s24b = 'spitzer/SER_E_3_A_MIPS1_mosaic_lrg.fits'

#Herschel data
H70 = '/data/damian/data/aquila/Herschel/aquilaM2-070.fits'
H160 = '/data/damian/data/aquila/Herschel/aquilaM2-160.fits'
H250 = '/data/damian/data/aquila/Herschel/aquilaM2-250.fits'
H350 = '/data/damian/data/aquila/Herschel/aquilaM2-350.fits'
H500 = '/data/damian/data/aquila/Herschel/aquilaM2-500.fits'

#SCUBA-2
w40_s450_IR2 = 'aquila/serpens_south/IR2/Aquila_extS2_450-4-30am.fits'
w40_s850_IR2 = 'aquila/serpens_south/IR2/Aquila_extS2_850-4-30am.fits'

#VLA - 3.6cm
HII_small = '/data/damian/data/aquila/radio/VLA_8_7GHz_W40_lrg.fits'

#21cm
HII_large = '/data/damian/data/aquila/radio/aquila_VLAss_21cm.fits'

####################################################
############        MAKE MAPS HERE     #############
####################################################

figures = [H70,w40_s450_IR2,w40_s850_IR2,HII_large]
max = [0.0,0.0,0.00,0.0]
min = [6000.0,0.18,0.06,0.05]
label = ['$70\mu m$ (MJy/Sr)','$450\mu m$ (Jy/pix)','$850\mu m$ (Jy/pix)','$21cm$ (Jy/beam)']

x = 0.4
dx = 0.05
y = 0.4
dy = 0.05

i = 0
Xo = 0
Yo = 0

big = plt.figure(figsize=(10, 13))

for i in range (len(figures)):
    print i
    if i%2 == 0:
        Yo = Yo + dy
        Xo = dx
        MAP(i,figures[i],Xo,Yo,x,y,min[i],max[i],label[i])
    elif i%2 == 1:
        Xo = Xo+x+dx
        Yo = Yo
        MAP(i,figures[i],Xo,Yo,x,y,min[i],max[i],label[i])
        Yo = Yo + y
    i = i + 1

plt.show()
