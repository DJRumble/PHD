#!/usr/astro64/python2.6/bin/python
#plotsed.py
# plot Spitzer SED

import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import readc2dcat
import atpy,aplpy

wcf='Catalogues/wise_prelim2.tbl'
scf = 'Catalogues/kmipstbl.fits'
print "Reading catalogue ",wcf
wtbl=atpy.Table(wcf,type='IPAC')
wcat=wtbl.data
stbl=atpy.Table(scf,type='FITS')
scat = stbl.data

#Spitzer and wise conversions
bandname =   ['B',  'V',  'Rc',  'J',    'H',   'K',    'I1',  'I2',  'I3',  'I4',  'M1', 'M2', 'W1', 'W2','W3','W4'] 
wave =       [0.44, 0.55, 0.64,  1.235,  1.662, 2.159,  3.6,   4.5,   5.8,   8.0,   24.0, 70.,3.4, 4.6, 12.0, 22.0]
zeroflux =   [4130.,3781.,3080., 1594.,  1024., 666.7,  277.5, 179.5, 116.6, 63.1,  7.3,  0.81, 309540.0, 171787.0, 31674.0, 8363.0]
#Wise zero fluxes from Jarrett et al. 2011 / wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html converted from Jy to mJy

#zerocgs  = [6.39339e-05,  3.74600e-05,  2.25360e-05,  3.13214e-06,  1.11103e-06,  4.28659e-07,  6.41719e-08,  2.65660e-08,  1.03879e-08,  2.95485e-09,  3.79828e-11,  4.95422e-13,0,0,0,0]

#make dictionaries
waves = dict(zip(bandname, wave))
zerofluxes = dict(zip(bandname,zeroflux))
#zerocgss = dict(zip(bandname,zerocgs))
clight = 2.99e8

#Wavelengths to plot: WISE and Spitzer
slams =np.array([1.235,  1.662, 2.159,  3.6,   4.5,  5.0, 8.0, 24.0, 70.])
wlams =np.array([3.4,4.6,12.0,22.0])

def makemag(fluxes,band):
    return -2.5 * np.log10(fluxes*1e-3 / zerofluxes[band])

def makeergs(fluxes,band):
    return fluxes/waves[band]**2 * 1.e-3 * (clight*100 * 1.e-19)

def makefluxes(mags,band):
    return zerofluxes[band]*10**(-mags/2.5)

def makefluxunc(mags,munc, band):
    return makefluxes(mags,band)*np.log(10)/2.5 * munc

# F(mJy) -> Mag   (Cohen et al. 2003)
#jm = np.where(((cat['J_flux_c'] != -99.) & (cat['J_D_flux_c'] != 0)),makemag(cat['J_flux_c'],'J'),0)

# F(mJy) -> F(erg/cm2/s/mu) 
jff =  np.where(((scat['J_flux_c'] != -99.) & (scat['J_D_flux_c'] != 0)),makeergs(scat['J_flux_c'],'J'),-99.)
hff =  np.where(((scat['H_flux_c'] != -99.) & (scat['H_D_flux_c'] != 0)),makeergs(scat['H_flux_c'],'H'),-99.)
kff =  np.where(((scat['Ks_flux_c'] != -99.) & (scat['Ks_D_flux_c'] != 0)),makeergs(scat['Ks_flux_c'],'K'),-99.)
i1ff =  np.where(((scat['IR1_flux_c'] != -99.) & (scat['IR1_D_flux_c'] != 0)),makeergs(scat['IR1_flux_c'],'I1'),-99.)
i2ff =  np.where(((scat['IR2_flux_c'] != -99.) & (scat['IR2_D_flux_c'] != 0)),makeergs(scat['IR2_flux_c'],'I2'),-99.)
i3ff =  np.where(((scat['IR3_flux_c'] != -99.) & (scat['IR3_D_flux_c'] != 0)),makeergs(scat['IR3_flux_c'],'I3'),-99.)
i4ff =  np.where(((scat['IR4_flux_c'] != -99.) & (scat['IR4_D_flux_c'] != 0)),makeergs(scat['IR4_flux_c'],'I4'),-99.)
m1ff =  np.where(((scat['MP1_flux_c'] != -99.) & (scat['MP1_D_flux_c'] != 0)),makeergs(scat['MP1_flux_c'],'M1'),-99.)
m2ff =  np.where(((scat['MP2_flux_c'] != -99.) & (scat['MP2_D_flux_c'] != 0)),makeergs(scat['MP2_flux_c'],'M2'),-99.)
#Fix M1 and M2 flux for V1121 Oph using Gaia aperture flux
#print 'V1121 24 micron flux', scat['MP1_flux_c'][4]
m1ff[4] = makeergs(4010,'M1')
m2ff[4] = makeergs(2400,'M2')
#WISE
w1ff =  np.where(((wcat['w1mpro'] != -99.) & (wcat['w1sigmpro'] != 0)),makeergs(makefluxes(wcat['w1mpro'],'W1'),'W1'),-99.)
w2ff =  np.where(((wcat['w2mpro'] != -99.) & (wcat['w2sigmpro'] != 0)),makeergs(makefluxes(wcat['w2mpro'],'W2'),'W2'),-99.)
w3ff =  np.where(((wcat['w3mpro'] != -99.) & (wcat['w3sigmpro'] != 0)),makeergs(makefluxes(wcat['w3mpro'],'W3'),'W3'),-99.)
w4ff =  np.where(((wcat['w4mpro'] != -99.) & (wcat['w4sigmpro'] != 0)),makeergs(makefluxes(wcat['w4mpro'],'W4'),'W4'),-99.)
w1f = makefluxes(wcat['w1mpro'],'W1')
w1df = makefluxunc(wcat['w1mpro'],wcat['w1sigmpro'],'W1')
w2f = makefluxes(wcat['w2mpro'],'W2')
w2df = makefluxunc(wcat['w2mpro'],wcat['w2sigmpro'],'W2')
w3f = makefluxes(wcat['w3mpro'],'W3')
w3df = makefluxunc(wcat['w3mpro'],wcat['w3sigmpro'],'W3')
w4f = makefluxes(wcat['w4mpro'],'W4')
w4df = makefluxunc(wcat['w4mpro'],wcat['w4sigmpro'],'W4')

fig1 = plt.figure(1)
fig1.set_size_inches((9,3))
fig1.subplots_adjust(wspace=0.001)
nx=5;ny=1;
fig1.clf()

wstars = [5,7,10,11,12]
sstars = [1,2,3,4,5]
for i in np.arange(5):
    wstar = wstars[i]-1
    sstar = sstars[i]-1
#    print wcat['designation'][wstar], scat['2MASS_name'][sstar]

    wfluxes_star = np.array([w1ff[wstar],w2ff[wstar],w3ff[wstar],w4ff[wstar]])    
    wflam_star=wfluxes_star*wlams
#    print wstar, wfluxes_star, wflam_star

    sfluxes_star=np.array([jff[i],hff[i],kff[i],i1ff[i],i2ff[i],i3ff[i],i4ff[i],m1ff[i],m2ff[i]])
    sflam_star=np.where(sfluxes_star>0,sfluxes_star*slams,0)
    #print sflam_star,slams
    
    plt = fig1.add_subplot(ny,nx,i+1)

    plt.loglog(wlams,wflam_star,\
           color='red',marker='o',linestyle='None',markersize=3)
    plt.loglog(slams,sflam_star,\
           color='white',marker='o',linestyle='None',markersize=3)
#           label='IRAS F16544$-$1604')
#plt.figtext(0.4,0.83,'WISE')
#print plt.gca()
    axes=fig1.gca()
    for label in axes.get_xticklabels() + axes.get_yticklabels():
        label.set_fontsize('x-small')

    if (i < nx):
        plt.set_xlabel(r'$\lambda$ ($\mu$m)',size='x-small')
    if (np.mod(i,nx)==0):
        plt.set_ylabel(r'$\lambda F_\lambda$ (erg/s/cm$^{2}$)',size='x-small')
    else:
        plt.set_yticklabels('',visible='False',size='x-small')
    plt.axis([0.4,150,1e-14,1e-6])
#    print 'i = ', i
    plt.text(0.7,5e-8,'M%d'%(i+1),size='small')
    print "%12s &%15s &\$%4.1f\pm%4.1f\$ &\$%4.1f\pm%4.1f\$ &\$%4.1f\pm%4.1f\$ &\$%4.1f\pm%4.1f\$ &\$%4.1f\pm%4.1f\$ &\$%4.2f\pm%4.2f\$ &\$%4.2f\pm%4.2f\$ &\$%4.2f\pm%4.2f\$ &\$%4.2f\pm%4.2f\$(%s%s) &\$%4.2f\pm%4.2f\$(%s%s)\\\\"%\
    (scat['2MASS_name'][sstar],scat['object_type'][sstar],\
    scat['J_flux_c'][sstar],scat['J_D_flux_c'][sstar],\
    scat['H_flux_c'][sstar],scat['H_D_flux_c'][sstar],\
    scat['Ks_flux_c'][sstar],scat['Ks_D_flux_c'][sstar],\
    w1f[wstar],w1df[wstar],\
    w2f[wstar],w2df[wstar],\
    w3f[wstar],w3df[wstar],\
    w4f[wstar],w4df[wstar],\
    scat['MP1_flux_c'][sstar],scat['MP1_D_flux_c'][sstar],\
    scat['MP2_flux_c'][sstar],scat['MP2_D_flux_c'][sstar],\
    scat['MP2_Q_det_c'][sstar],scat['MP2_imtype'][sstar],\
    scat['MP3_flux_c'][sstar],scat['MP3_D_flux_c'][sstar],\
    scat['MP3_Q_det_c'][sstar],scat['MP3_imtype'][sstar])
fig1.savefig('PS/plotwiseseds.eps', format='eps')

# if (__name__=='__main__'):
#     main()
