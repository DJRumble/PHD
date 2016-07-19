#20151103
#Damian Rumble, UoE

#this script produces histograms of clump mass

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

#################################################################################
########## Bulk script ############

#ORIGINAL data
data =  np.loadtxt('data/SMM_master.tab',dtype='string',comments="#")
data16 =  np.loadtxt('data/SMM_15K.tab',dtype='string',comments="#")
#Starless core
data_starless =  np.loadtxt('data/SMM_master_noYSO.tab',dtype='string',comments="#")
#Protostars
data_protostar =  np.loadtxt('data/SMM_master_YSO.tab',dtype='string',comments="#")
data_heat =  np.loadtxt('data/SMM_master_OBheated.tab',dtype='string',comments="#")
#Protostars
data_noheat =  np.loadtxt('data/SMM_master_noOBheated.tab',dtype='string',comments="#")



#mass

#DATA
DATA = map(float, data[:,4])
DATA16 = map(float, data16[:,4])
DATA_starless = map(float, data_starless[:,4])
DATA_proto = map(float, data_protostar[:,4])
DATA_heat = map(float, data_heat[:,4])
DATA_noheat = map(float, data_noheat[:,4])
ERR = map(float, data[:,5])
ERR_starless = map(float, data_starless[:,5])
#ERR_proto = map(float, data_protostar[:5])
#HISTOGRAM
bins = np.logspace(start=-2,stop=2,num=1000,base=10)
#print DATA

fig = plt.figure()


###PLOT###
fig.add_subplot(131)
n, bins, patches = hist(DATA, bins,normed=False,label=('Real temp.'),color='r',histtype='step',cumulative=-1)
n, bins, patches = hist(DATA16, bins,normed=False,label=('Temp. = 15K'),color='b',histtype='step',cumulative=-1)
#Labels
xlabel('Mass (M$_{\odot}$)')
ylabel('Number')

#Plot Salpetre
Y = [((0.005*i)**(-1.35)) for i in bins]
plot(bins,Y,'k--')

yscale('log')
xscale('log')
legend(loc=3)

ylim([1,1000])
xlim([0.01,200])

fig.add_subplot(132)
n, bins, patches = hist(DATA_starless, bins,normed=False,label=('Pre'),color='g',histtype='step',cumulative=-1)
n, bins, patches = hist(DATA_proto, bins,normed=False,label=('Proto'),color='y',histtype='step',cumulative=-1)
#Labels
xlabel('Mass (M$_{\odot}$)')
ylabel('Number')

#Plot Salpetre
Y = [((0.005*i)**(-1.35)) for i in bins]
plot(bins,Y,'k--')

yscale('log')
xscale('log')
legend(loc=3)

ylim([1,500])
xlim([0.01,200])

fig.add_subplot(133)
n, bins, patches = hist(DATA_heat, bins,normed=True,label=('OB heating'),color='magenta',histtype='step',cumulative=-1)
n, bins, patches = hist(DATA_noheat, bins,normed=True,label=('No OB heating'),color='cyan',histtype='step',cumulative=-1)
#HIST errors
#errorbar(bins,n,xerr=0,yerr=0,fmt='k+')

#Labels
xlabel('Mass (M$_{\odot}$)')
ylabel('Normalised frequency')
#label.set_fontsize('x-small')

#Plot Salpetre
Y = [((0.3*i)**(-1.35)) for i in bins]
plot(bins,Y,'k--')

yscale('log')
xscale('log')
legend(loc=3)

#plt.axvline(x=1.865,linestyle=':',color='k')
#plt.axhline(y=0.5,linestyle=':',color='k')

ylim([0.01,1])
xlim([0.01,200])

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_CMF.pdf'%(date))

show()



