#20151105

#A script for producing scatter plots for mass vs radiu

import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt
import time

################################################################################
########## Bulk script ############

#ORIGINAL data
data_YSOHOT =  np.loadtxt('data/SMM_master_OBheated_YSO.tab',dtype='string',comments="#")
#Protostars
data_YSOCOLD =  np.loadtxt('data/SMM_master_real_noOBheated_YSO.tab',dtype='string',comments="#")
data_noYSOHOT =  np.loadtxt('data/SMM_master_OBheated_noYSO.tab',dtype='string',comments="#")
#Protostars
data_noYSOCOLD =  np.loadtxt('data/SMM_master_real_noOBheated_noYSO.tab',dtype='string',comments="#")
#constant temp 15K
data15_noYSOHOT =  np.loadtxt('data/SMM_master_OBheated_noYSO+15Kclumps.tab',dtype='string',comments="#")
data15_YSOHOT =  np.loadtxt('data/SMM_master_OBheated_YSO+15Kclumps.tab',dtype='string',comments="#")

#MASS
mass_ysohot = map(float, data_YSOHOT[:,4])
dmass_ysohot = map(float, data_YSOHOT[:,5])
mass_noysohot = map(float, data_noYSOHOT[:,4])
dmass_noysohot = map(float, data_noYSOHOT[:,5])
mass_ysocold = map(float, data_YSOCOLD[:,4])
dmass_ysocold = map(float, data_YSOCOLD[:,5])
mass_noysocold = map(float, data_noYSOCOLD[:,4])
dmass_noysocold = map(float, data_noYSOCOLD[:,5])
mass15_ysohot = map(float, data15_YSOHOT[:,25])
dmass15_ysohot = map(float, data15_YSOHOT[:,26])
mass15_noysohot = map(float, data15_noYSOHOT[:,25])
dmass15_noysohot = map(float, data15_noYSOHOT[:,26])
#RADIUS
d_ysohot = map(float, data_YSOHOT[:,11])
d_noysohot = map(float, data_noYSOHOT[:,11])
d_ysocold = map(float, data_YSOCOLD[:,11])
d_noysocold = map(float, data_noYSOCOLD[:,11])
d15_ysohot = map(float, data15_YSOHOT[:,11])
d15_noysohot = map(float, data15_noYSOHOT[:,11])

radius_ysohot = [i/2. for i in d_ysohot]
radius_noysohot = [i/2. for i in d_noysohot]
radius_ysocold = [i/2. for i in d_ysocold]
radius_noysocold = [i/2. for i in d_noysocold]
radius15_ysohot = [i/2. for i in d15_ysohot]
radius15_noysohot = [i/2. for i in d15_noysohot]

fig1 = plt.figure()

fig1.add_subplot(121)

ylabel('Clump mass (M$_{\odot}$)')
xlabel('Clump radius (pc)')

errorbar(radius_ysocold,mass_ysocold,yerr=dmass_ysocold,xerr=0,fmt='bo',ecolor='b',label=('No OB heating'))
#errorbar(radius15_ysohot,mass15_ysohot,yerr=dmass15_ysohot,xerr=0,fmt='ko',ecolor='r',label=('15K heating'))
errorbar(radius_ysohot,mass_ysohot,yerr=dmass_ysohot,xerr=0,fmt='ro',ecolor='r',label=('OB heating'))
errorbar(radius_noysocold,mass_noysocold,yerr=dmass_noysocold,xerr=0,fmt='wo',ecolor='b',label=('No OB heating'))
#errorbar(radius15_noysohot,mass15_noysohot,yerr=dmass15_noysohot,xerr=0,fmt='ko',ecolor='k',label=('15K heating'))
errorbar(radius_noysohot,mass_noysohot,yerr=dmass_noysohot,xerr=0,fmt='wo',ecolor='r',label=('OB heating'))


yscale('log')
xscale('log')

legend(loc=4)

X = [0.005,0.25]

#PLOT MODELS
#T1 = 6
T2 = 15
T3 = 22
#MJ1 = [2.71*T1*i for i in X]
MJ2 = [2.71*T2*i for i in X]
MJ3 = [2.71*T3*i for i in X]
#ML = [i**3. for i in radius_clumps]

#plot(X,MJ1,'k-')
#plot(X,MJ2,'k--')
#plot(X,MJ3,'k:')
#plot(radius_clumps,ML,'r-')

axvline(x=0.025,color='k',linestyle='--')

xlim([0.005,0.25])
ylim([0.005,200])

fig1.add_subplot(122)

xlabel('Clump radius (pc)')

#errorbar(radius_noysocold,mass_noysocold,yerr=dmass_noysocold,xerr=0,fmt='wo',ecolor='b',label=('No OB heating'))
errorbar(radius15_noysohot,mass15_noysohot,yerr=dmass15_noysohot,xerr=0,fmt='ko',ecolor='k',label=('15K heating'))
errorbar(radius_noysohot,mass_noysohot,yerr=dmass_noysohot,xerr=0,fmt='wo',ecolor='r',label=('OB heating'))


yscale('log')
xscale('log')

legend(loc=4)

X = [0.005,0.25]

#PLOT MODELS
#T1 = 6
T2 = 15
T3 = 22
#MJ1 = [2.71*T1*i for i in X]
MJ2 = [2.71*T2*i for i in X]
MJ3 = [2.71*T3*i for i in X]
#ML = [i**3. for i in radius_clumps]

#plot(X,MJ1,'k-')
plot(X,MJ2,'k--')
plot(X,MJ3,'k-')
#plot(radius_clumps,ML,'r-')

axvline(x=0.025,color='k',linestyle='--')

xlim([0.005,0.25])
ylim([0.005,200])


#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_envelope.pdf'%(date))


show()
