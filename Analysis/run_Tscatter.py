#20160713

#Damian Rumble, UoE

#This script produces a scatter plot of stability vs T

#I also aim to introduce a box and whisker diagram on to the plot for heated clumps at 15K. 



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

##### DEFFINE c0nstants ##########

au = 149597871000
p = 3.08567758E16
M_x = 1.989E30
h = 6.626068E-34
c = 2.99792458E8
k = 1.3806488E-23
m_h2 = 3.34745E-27 #g
G = 6.67384E-11
pi = 3.14159265359
mu = 2.333333  #ratio of H2 to He (5:1)

A = 5.03517447875E28 #m2
aN_h2 = 1E21 #cm-2
kappa = 0.00012 #m-2/kg

#################################################################################
########## functions ############

def BB(T):
    l = 850E-6
    nu = c / l
    #equation of a blackbody function
    B_a = ((2.0*h*(nu**3.0))/(c**2.0))
    B_b = (1.0/(np.exp((h*nu)/(k*T))-1.0))
    B = B_a*B_b
    return B

def stability(T):
    S = 10E-27 #Jy
    d = 250*p
    M = ((S*(d**2.))/(BB(T)*kappa))/M_x
    MJ = 1.9*(T/10.)*(0.07/0.07)
    return M/MJ


#################################################################################
########## Bulk script ############


#No 15K clumps
datanoOBnoYSO =  np.loadtxt('data/SMM_master_real_noOBheated_noYSO.tab',dtype='string',comments="#") 
datanoOBYSO =  np.loadtxt('data/SMM_master_real_noOBheated_YSO.tab',dtype='string',comments="#")
dataOBnoYSO =  np.loadtxt('data/SMM_master_OBheated_noYSO.tab',dtype='string',comments="#")
dataOBYSO = np.loadtxt('data/SMM_master_OBheated_YSO.tab',dtype='string',comments="#")

dataOBnoYSO15K =  np.loadtxt('data/SMM_master_OBheated_noYSO+15Kclumps.tab',dtype='string',comments="#")
dataOBYSO15K = np.loadtxt('data/SMM_master_OBheated_YSO+15Kclumps.tab',dtype='string',comments="#")

#Temp
T1a = map(float, datanoOBnoYSO[:,6])
dT1a = map(float, datanoOBnoYSO[:,7])
T1b = map(float, datanoOBYSO[:,6])
dT1b = map(float, datanoOBYSO[:,7])
T2a = map(float, dataOBnoYSO[:,6])
dT2a = map(float, dataOBnoYSO[:,7])
T2b = map(float, dataOBYSO[:,6])
dT2b = map(float, dataOBYSO[:,7])

mdT1a = np.mean(dT1a)
mdT1b = np.mean(dT1b)
mdT2a = np.mean(dT2a)
mdT2b = np.mean(dT2b)

mdT1 = (mdT1a + mdT1b)/2.
mdT2 = (mdT2a + mdT2b)/2.

#Stability
J1a = map(float, datanoOBnoYSO[:,17])
dJ1a = map(float, datanoOBnoYSO[:,18])
J1b = map(float, datanoOBYSO[:,17])
dJ1b = map(float, datanoOBYSO[:,18])
J2a = map(float, dataOBnoYSO[:,17])
dJ2a = map(float, dataOBnoYSO[:,18])
J2b = map(float, dataOBYSO[:,17])
dJ2b = map(float, dataOBYSO[:,18])

mdJ1a = np.mean(dJ1a)
mdJ1b = np.mean(dJ1b)
mdJ2a = np.mean(dJ2a)
mdJ2b = np.mean(dJ2b)

mdJ1 = (mdJ1a + mdJ1b)/2.
mdJ2 = (mdJ2a + mdJ2b)/2.

#corrected Stability
J15a = map(float, dataOBnoYSO15K[:,38])
dJ15a = map(float, dataOBnoYSO15K[:,39])
J15b = map(float, dataOBYSO15K[:,38])
dJ15b = map(float, dataOBYSO15K[:,39])

J15 =  J15a + J15b
mJ15 = np.mean(J15)
sJ15 = np.std(J15)
lJ15 = len(J15)

print mJ15, sJ15, lJ15 

norm = np.random.normal(mJ15, sJ15, lJ15)

###PLOT####

fig1 = plt.figure()

### Box plot ###
x = [-15,15]

#bp = plt.boxplot((norm,norm), positions=x, patch_artist=True, showfliers=False)

#for box in bp['boxes']:
##    box.set( color='k',facecolor = 'k', linewidth=2)
#for median in bp['medians']:
#    median.set(color='w', linewidth=2)
#for whisker in bp['whiskers']:
#    whisker.set(color='k', linewidth=2)
#for cap in bp['caps']:
#    cap.set(color='k', linewidth=2)


### Markers ###
errorbar(T1a,J1a,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w',markersize=3)
errorbar(T1b,J1b,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w',markersize=3)
errorbar(T2a,J2a,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w',markersize=3)
errorbar(T2b,J2b,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w',markersize=3)

### Mean standing error bar ###
errorbar(35,1.5,yerr=mdJ1,xerr=mdT1,marker='+',markeredgecolor='r', markerfacecolor='r', ecolor='r', linewidth=2)
errorbar(36,1.6,yerr=mdJ2,xerr=mdT2,marker='+',markeredgecolor='b', markerfacecolor='b', ecolor='b', linewidth=2)


xlim((6,45))
ylim([0,5])

T = [5,45]
T = np.arange(5,45,0.1)

J = [stability(i) for i in T]

plt.plot(T,J,'k', linewidth=2)

#yscale('log')


axhline(y=1,color='k',linestyle='--')
axhline(y=2,color='k',linestyle='-')

#axvline(x=14.8,color='k',linestyle='--')
#axvline(x=17.2,color='k',linestyle=':')
#axvline(x=12.4,color='k',linestyle=':')

xlabel('Temp. (K)')
ylabel('M/M$_{\mathrm{J}}$')

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_temp_scatter_JCMTGBS.pdf'%(date))


show()
