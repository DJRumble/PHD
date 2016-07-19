#!/usr/local64/python2.6/bin/python

import aplpy
import numpy as np
import matplotlib.pyplot as plt
import cat

s850 = aplpy.FITSFigure('FITS/ngc1333_s2sro_s850_m2-co_cal653.fits')
s850.show_colorscale(cmap='gist_heat',vmax=1,vmin=-0.1,stretch='sqrt')
a = np.arange(0,10,1.0)
clevels = 0.02*5.0**a
print clevels
s850.show_contour('./FITS/ngc1333_s2sro_850_m50_cal500.fits',levels=clevels, colors="black")
# Overplot SCUBA cores
cat103 = cat.Cat('/h/hatchell/Perseus/jcmt10_work/Catalogues/cat103.apl',namep=1,ra=2,dec=3)
cat103.read()
ra = cat103.column('ra',type='float')
dec = cat103.column('dec',type='float')
s850.show_markers(ra,dec,edgecolor='white',facecolor='none',marker='o',s=20)
#... add_label works on Reduce version
#for source in cat103.rownames:
#    s850.add_label(ra,dec,source)

# Overplot Spitzer cores
cat103 = cat.Cat('/h/hatchell/Perseus/jcmt10_work/Catalogues/gutermuth08_ysos.apl',namep=1,ra=2,dec=3)
cat103.read()
ra = cat103.column('ra',type='float')
dec = cat103.column('dec',type='float')
s850.show_markers(ra,dec,edgecolor='yellow',facecolor='none',marker='x',s=20)
