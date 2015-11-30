#Damian Rumble, UoE
#20141120

#This code is for plotting minimaps of my cores. It needs to; Centre on the core location in the temp map, plot contours of 850Micron & Spitzer data. Show the clump (or Core) mask. And plot YSO on it. 

#####################################################
#import maths, ploting and astrophysical packages

import aplpy
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import matplotlib.lines as line
import numpy

#My modules
import YSO

def MAP(Ra,Dec,Xo,Yo,X,Y,i):
    #deffine image in APLPY
###################################################
#plot background DUST map - subtitute in as appropriate 

    print 'Mapping SMM'+str(i)

    temp = 'SerpensMWC297_20140414_IR1_freefree+jet450850temp%5.fits' 
    s850 = 'SerpensMWC297_20140414_IR1_s850_freefree+jet.fits'
    s450 = 'SerpensMWC297_20140414_IR1_s450_freefree+jet.fits'

    Fig = aplpy.FITSFigure(temp,figure=big,subplot=[Xo,Yo,X,Y])

#peak marks
    Fig.show_markers(Ra,Dec, layer='peak', edgecolor='black', marker='+',s=100, alpha=1)
    Fig.show_markers(Ra,Dec, layer='peak', edgecolor='black', marker='x',s=100, alpha=1)
    Fig.show_markers(Ra,Dec, layer='aperture', edgecolor='red', marker='o',s=2500, alpha=1)

#Deffine contours
    #cset1 = Fig.show_contour(s850, levels=(0.0079,0.013,0.026,0.053,0.079), linewidth = ('1','1','2','1','1'), colors=('black'))
    #cset1 = Fig.show_contour(s850, levels=(0.0066,0.011,0.022,0.044,0.066), linewidth = ('1','1','2','1','1'), colors=('black'))
    cset2 = Fig.show_contour(s450, levels=(0.08,0.16,0.32,0.48), linewidth = ('1','1','2','1','1'), colors=('black'))


#Deffine Markers
    YSO.GBS(Fig)
    #YSO.mwc297(Fig)

#Figure properties
    Fig.show_colorscale(vmax=45,vmin=8.0, stretch='linear')
    Fig.add_label(0.07, 0.85, '['+str(i)+']',size='medium', relative=True)
    Fig.tick_labels.set_font(size='x-small')
    Fig.recenter(Ra,Dec,height=0.02 ,width=0.02)
    Fig.ticks.show()
    #Fig.add_colorbar()

    Fig.ticks.set_color('black')
    Fig.hide_xaxis_label()
    Fig.hide_yaxis_label()

    return

####################################################
############        MAKE MAPS HERE     #############
####################################################

data = numpy.loadtxt('minimaps_clumps.txt')

x = 0.35
dx = 0.1
y = 0.09
dy = 0.019

i = 0
Xo = 0
Yo = 0

big = plt.figure(figsize=(10, 17))

for i in range (len(data)):
    print i
    if i%2 == 0:
        Xo = dx
        Yo = dy+y+Yo

        Ra = float(data[i][0])
        Dec = float(data[i][1])
        SMM = int(data[i][2])
        print Xo,Yo,SMM
        MAP(Ra,Dec,Xo,Yo,x,y,SMM)
    elif i%2 == 1:
        Xo = Xo+x+dx
        Yo = Yo
        Ra = float(data[i][0])
        Dec = float(data[i][1])
        SMM = int(data[i][2])
        print Xo,Yo, SMM
        MAP(Ra,Dec,Xo,Yo,x,y,SMM)
        
    i = i + 1

big.canvas.draw()
#big.add_colorbar()
plt.savefig('20141120_450b_minimap.pdf',format='pdf')
plt.show()
