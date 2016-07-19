#################################################################################
########### import modules  #############
#Standard modules
import numpy as np
import astropy.io.fits as pyfits
import sys
import os
import string
import mpfit
import matplotlib.pyplot as plt
import copy
from pylab import *
from scipy import stats
import time


#################################################################################
########## Bulk script #############

tempK = '/data/damian/temp_MKR/output/kernel/map/SerpensMain_20141223-auto/SerpensMain_20141223-autotemperature.fits'
tempDB = '/data/damian/temp_MKR_v10/output//map/SerpensMain_20141223_IR2extmask_s2_cal_JypixJH/SerpensMain_20141223_IR2extmask_s2_cal_JypixJH450850temp%5.fits'

label = ['KERNEL','DUAL-BEAM']

#open FITS files - KERNEL
hdulistK = pyfits.open(tempK)
dataK = hdulistK[0].data
#print data
locK = np.where((dataK<70.) & (dataK>5.) & (dataK!=0.0))
#print loc
DataK = dataK[locK]

#open FITS files - DUAL-BEAM
hdulistDB = pyfits.open(tempDB)
dataDB = hdulistDB[0].data
#print data
locDB = np.where((dataDB<70.) & (dataDB>5.) & (dataDB!=0.0))
#print loc
DataDB = dataDB[locDB]

n1, bin_edges1, patches1 = hist(DataK,bins=50,normed=True, histtype='step',label=str(label[0]),color='g')
n2, bin_edges2, patches2 = hist(DataDB,bins=50,normed=True, histtype='step',label=str(label[1]),color='k')

ylabel('Frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([0,70])
ylim([0,0.125])

print len(bin_edges1),np.median(DataK)
print len(bin_edges2),np.median(DataDB)

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('/data/damian/plots/ch5/%s_Tempmethod_histogram.pdf'%(date))

show()
