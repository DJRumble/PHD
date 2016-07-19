import aplpy
import matplotlib.pyplot as plt
fig = 'COratio_msk.fits'
cnt = 'Aquila_extS2_450-4-30am.fits'
Fig = aplpy.FITSFigure(fig)
Fig.show_colorscale(vmax=0.2,vmin=-0.02)
Fig.add_colorbar()
Fig.colorbar.set_axis_label_text('CO contamination fraction')
cset1 = Fig.show_contour(cnt, levels=(0.0185,1000), linewidth=3, colors=('w'))

plt.show()
