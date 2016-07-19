#!/usr/astro64/python2.6/bin/python
#plotsed.py
# plot Spitzer SED

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import readc2dcat
import atpy,aplpy

fields = ['SCO 1','SCO 2','SCO 3','SCO 4','SCO 5','SCO 6','CB68','L234E', 'ALL', 'IRAC YSOc']
#catfiles=['/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-SCO_1-v1.tbl','/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-SCO_2-v1.tbl','/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-SCO_3-v2.tbl','/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-SCO_4-v1.tbl','/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-SCO_5-v1.tbl','/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-SCO_6-v1.tbl','C2DCORES/CB68/catalog-CB68-v2.tbl','C2DCORES/L234E/catalog-L234E-v2.tbl', 'Catalogues/master4_edited2.tbl', '/h/hatchell/SCRATCH/SPITZER_SCO/CATALOGS/catalog-Sco-YSOc.tbl']
# FITS versions
catfiles=['./Catalogues/catalog-SCO_1-v2.fits','./Catalogues/catalog-SCO_2-v1.fits','./Catalogues/catalog-SCO_3-v2.fits','./Catalogues/catalog-SCO_4-v1.fits','./Catalogues/catalog-SCO_5-v1.fits','./Catalogues/catalog-SCO_6-v1.fits','./Catalogues/catalog-CB68-v2.fits','./Catalogues/catalog-L234E-v2.fits', 'Catalogues/catalog-SCO_ALL.fits', './Catalogues/catalog-Sco-YSOc.fits']
catalogues=dict(zip(fields,catfiles))

# Choose field here.  'all' takes forever.
newfield='CB68'
field = ''
if (field != newfield):
    field = newfield
    print "Reading catalogue ",catalogues[field]
    tbl=atpy.Table(catalogues[field],type='FITS')
    cat=tbl.data
    # RA fix for fits catalogue reading
    ra = 200.0+cat['RA']
    dec = cat['DEC']

# print "Reading catalogue..."
# cat=readc2dcat.readc2dcat('C2DCORES/CB68/catalog-CB68-v2.tbl')
# ra = cat['RA']
# dec = cat['DEC']


#from sed_chart_pub.pro
bandname =   ['B',  'V',  'Rc',  'J',    'H',   'K',    'I1',  'I2',  'I3',  'I4',  'M1', 'M2'] 
wave =       [0.44, 0.55, 0.64,  1.235,  1.662, 2.159,  3.6,   4.5,   5.8,   8.0,   24.0, 70.]
zeroflux =   [4130.,3781.,3080., 1594.,  1024., 666.7,  277.5, 179.5, 116.6, 63.1,  7.3,  0.81]
zerocgs  = [6.39339e-05,  3.74600e-05,  2.25360e-05,  3.13214e-06,  1.11103e-06,  4.28659e-07,  6.41719e-08,  2.65660e-08,  1.03879e-08,  2.95485e-09,  3.79828e-11,  4.95422e-13]

#make dictionaries
waves = dict(zip(bandname, wave))
zerofluxes = dict(zip(bandname,zeroflux))
zerocgss = dict(zip(bandname,zerocgs))
clight = 2.99e8

#Wavelengths to plot
lams =np.array([1.235,  1.662, 2.159,  3.6,   4.5,   5.8,   8.0,   24.0, 70.])

def makemag(fluxes,band):
    return -2.5 * np.log10(fluxes)*1e-3 / zerofluxes[band]

def makeergs(fluxes,band):
    return fluxes/waves[band]**2 * 1.e-3 * (clight*100 * 1.e-19) 

# F(mJy) -> Mag   (Cohen et al. 2003)
#jm = np.where(((cat['J_flux_c'] != -99.) & (cat['J_D_flux_c'] != 0)),makemag(cat['J_flux_c'],'J'),0)

# F(mJy) -> F(erg/cm2/s/mu) 
jff =  np.where(((cat['J_flux_c'] != -99.) & (cat['J_D_flux_c'] != 0)),makeergs(cat['J_flux_c'],'J'),-99.)
hff =  np.where(((cat['H_flux_c'] != -99.) & (cat['H_D_flux_c'] != 0)),makeergs(cat['H_flux_c'],'H'),-99.)
kff =  np.where(((cat['Ks_flux_c'] != -99.) & (cat['Ks_D_flux_c'] != 0)),makeergs(cat['Ks_flux_c'],'K'),-99.)
i1ff =  np.where(((cat['IR1_flux_c'] != -99.) & (cat['IR1_D_flux_c'] != 0)),makeergs(cat['IR1_flux_c'],'I1'),-99.)
i2ff =  np.where(((cat['IR2_flux_c'] != -99.) & (cat['IR2_D_flux_c'] != 0)),makeergs(cat['IR2_flux_c'],'I2'),-99.)
i3ff =  np.where(((cat['IR3_flux_c'] != -99.) & (cat['IR3_D_flux_c'] != 0)),makeergs(cat['IR3_flux_c'],'I3'),-99.)
i4ff =  np.where(((cat['IR4_flux_c'] != -99.) & (cat['IR4_D_flux_c'] != 0)),makeergs(cat['IR4_flux_c'],'I4'),-99.)
m1ff =  np.where(((cat['MP1_flux_c'] != -99.) & (cat['MP1_D_flux_c'] != 0)),makeergs(cat['MP1_flux_c'],'M1'),-99.)
m2ff =  np.where(((cat['MP2_flux_c'] != -99.) & (cat['MP2_D_flux_c'] != 0)),makeergs(cat['MP2_flux_c'],'M2'),-99.)

stars=np.where((ra>254.328) & (ra < 254.336) & (dec > -16.162) & (dec < -16.150))
pstar=np.where(cat['object_type'][stars] =='YSOc_red',stars,0)
pstar= pstar[np.nonzero(pstar)]
nomips=np.where(cat['object_type'][stars] =='star_no_MP1',stars,0)    
nomips=nomips[np.nonzero(nomips)]
print pstar,nomips
print 'stars = ', cat['object_type'][stars]

fluxes_pstar=np.array([jff[pstar],hff[pstar],kff[pstar],i1ff[pstar],i2ff[pstar],i3ff[pstar],i4ff[pstar],m1ff[pstar],m2ff[pstar]])
fluxes_nomips=np.array([jff[nomips],hff[nomips],kff[nomips],i1ff[nomips],i2ff[nomips],i3ff[nomips],i4ff[nomips],m1ff[nomips],m2ff[nomips]])
fluxes_pstar=fluxes_pstar.flatten()
fluxes_nomips=fluxes_nomips.flatten()
print fluxes_nomips
flam_pstar=fluxes_pstar*lams
flam_nomips=fluxes_nomips*lams


#fluxes_nomips=np.array([cat['J_flux_c'][nomips],cat['H_flux_c'][nomips],cat['Ks_flux_c'][nomips],cat['IR1_flux_c'][nomips],cat['IR2_flux_c'][nomips],cat['IR3_flux_c'][nomips],cat['IR4_flux_c'][nomips],cat['MP1_flux_c'][nomips],cat['MP2_flux_c'][nomips]])

plt.figure(1)
plt.clf()
plt.subplot(221)
plt.loglog(lams,flam_pstar,\
           color='red',marker='o',linestyle='None')
#           label='IRAS F16544$-$1604')
#plt.figtext(0.25,0.75,'IRAS F16544$-$1604')
plt.loglog(lams,flam_nomips,\
           color='green',marker='+',linestyle='None')
#           label='SE star',font='serif')
#plt.figtext(0.15,0.68,'SE star')
plt.figtext(0.4,0.83,'CB68')
print plt.gca()
axes=plt.gca()
plt.setp(axes,'xlim',(1,100))
plt.setp(axes,'ylim',(1e-14,1e-6))
plt.xlabel('$\lambda$ ($\mu$m)')
plt.ylabel('$\lambda F_\lambda$ (erg/s/cm$^{2}$)')
#plt.show()

# Fig.9: CB68 image
img = aplpy.FITSFigure('MOSAICS/Rotated/CB68_COMB_IRAC2_mosaic_rot.fits')
print "Setting small scale x tick spacing"
xtickspacing = 1/16.
img.set_tick_labels_format('hh:mm:ss','dd:mm:ss')
img.set_tick_labels_font(size='x-large')
img.set_axis_labels_size('x-large')

img.show_colorscale(cmap='gist_gray',stretch='log',pmin=20,pmax=98)
img.show_colorbar()
img.show_contour('MOSAICS/Rotated/CB68_FBCD_MIPS2_rotated.fits',levels=[50,100,500],colors=['black','black','black'])
img.show_markers(ra[pstar],dec[pstar],marker='o',edgecolor='r',s=100)
img.show_markers(ra[nomips],dec[nomips],marker='+',edgecolor='g',s=280)
# img.show_markers(ra[np.where(stardusts)],dec[np.where(stardusts)],marker='x',edgecolor='r',s=180)



# if (__name__=='__main__'):
#     main()
