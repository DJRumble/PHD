#20160516 
#Damian Rumble, UoE

#This is a script for producing comparable temp.histograms



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
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'
cupidr = '/star/bin/cupid'

#################################################################################
########## Bulk script #############

#SEE LATER HALF OR RUN21_FINDBACK.py for original version

CO = 'Aquila_noco_original_Sratio_MSK.fits'#'Aquila_extS2nosm-4am_Sratio_MSK.fits'
noCO = 'Aquila_noco_s2nosm-4am_Sratio_MSK.fits'

#open FITS files
hdulist1 = pyfits.open(CO)
data1 = hdulist1[0].data
#print data
loc = np.where((data1<15.) & (data1>0.) & (data1!=0.0))
CO_data = data1[loc]

#open FITS files
hdulist2 = pyfits.open(noCO)
data2 = hdulist2[0].data
#print data
loc = np.where((data2<15.) & (data2>0.) & (data2!=0.0))
noCO_data = data2[loc]

n1, bin_edges1, patches1 = hist(CO_data,bins=100,normed=False, histtype='step',color='k',cumulative=0)
n2, bin_edges2, patches2 = hist(noCO_data,bins=100,normed=False, histtype='step',color='g',cumulative=0)

ylabel('Normalised frequency')
xlabel('Flux ratio')

xlim([3,14])
ylim([0,1200])

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('%s_COratio_hist.pdf'%(date))

show()

