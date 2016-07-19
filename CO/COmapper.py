map = 'mosaics/GOOD/Aquila_COcontamination_mskTH.fits'
S2 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_noco_extS2nosm_s850-4am.fits'
mask = 'COmask.fits'

import numpy 
import aplpy
import matplotlib.pyplot as plt

Fig = aplpy.FITSFigure(map)
Fig.show_colorscale(vmax=0.5,vmin=-0.01, stretch='sqrt') 
Fig.add_colorbar()
cset1 = Fig.show_contour(mask, levels=(0.0025,1111), linewidth=2, colors=('black'))
cset3 = Fig.show_contour(S2, levels=(0.0025,1111), linewidth=2, colors=('white'))
Fig.recenter(277.8925,-2.13200, height=0.51 ,width=0.42)

plt.show()
