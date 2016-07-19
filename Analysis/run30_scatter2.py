#20160218
#Damian Rumble, UoE

#this script produces scatter plots of clump properties

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

def CD(T):
    S = 10E-27 #Jy
    M = 1*M_x #kg
    d = 250*p
    N = (S*(d**2.))/(M*BB(T)*kappa)
    return N

def stability1(T):
    S = 10E-27 #Jy
    d = 250*p
    M = ((S*(d**2.))/(BB(T)*kappa))/M_x
    MJ = 1.9*(T/10.)*(0.07/0.07)
    return M/MJ

def stability2(N):
    S = 10E-27 #Jy
    d = 250*p
    T = 15.0
    m = 0.5*M_x #kg
    M = (m*N)/M_x
    MJ = 1.9*(T/10.)*(0.07/0.07)
    return M/MJ

    

#################################################################################
########## Bulk script ############


match_prototemp =  np.loadtxt('data/SMM-proto-match.tab',dtype='string',comments="#") 
prototemp =  np.loadtxt('data/GBSprotostar_temperatures_master.tab',dtype='string',comments="#") 


#TEMP
protoT = map(float, prototemp[:,5])
protodT = map(float, prototemp[:,6])
centerT = map(float, prototemp[:,7])
centerdT = map(float, prototemp[:,8])

match_protoT = map(float, match_prototemp[:,5])
match_protodT = map(float, match_prototemp[:,6])
match_clumpT = map(float, match_prototemp[:,17])
match_clumpdT = map(float, match_prototemp[:,18])

#Class
tbol = map(float, prototemp[:,10])
Mtbol = map(float, match_prototemp[:,10])

Class0_pT = []
Class0_cT = []
ClassI_pT = []
ClassI_cT = []

for i in range(len(prototemp)):
    if tbol[i] <= 70.:
        Class0_pT.append(protoT[i])
        Class0_cT.append(centerT[i])
    else:
        ClassI_pT.append(protoT[i])
        ClassI_cT.append(centerT[i])

MClass0_pT = []
MClass0_cT = []
MClassI_pT = []
MClassI_cT = []


for i in range(len(match_prototemp)):
    if Mtbol[i] <= 70.:
        MClass0_pT.append(match_protoT[i])
        MClass0_cT.append(match_clumpT[i])
    else:
        MClassI_pT.append(match_protoT[i])
        MClassI_cT.append(match_clumpT[i])


#### the scatter plots: ####

fig = plt.figure(1, figsize=(16, 8))

### CD vs T ###
fig.add_subplot(122)
errorbar(Class0_pT,Class0_cT,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
errorbar(ClassI_pT,ClassI_cT,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='k', markerfacecolor='None', ecolor='k')

meanx = np.mean(protodT)
meany = np.mean(centerdT)

errorbar(30,10,yerr=meany,xerr=meanx,linestyle='None',marker='+',markeredgecolor='k', markerfacecolor='k', ecolor='k')

xlim([5,35])
ylim([5,35])

xlabel('Envelope Temp. (K)')
ylabel('Centeral Temp. (K)')

T = [5,45]

plt.plot(T,T,'k')

### MJ vs T ###
fig.add_subplot(121)

errorbar(MClass0_pT,MClass0_cT,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
errorbar(MClassI_pT,MClassI_cT,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='k', markerfacecolor='None', ecolor='k')


meanx = np.mean(match_protodT)
meany = np.mean(match_clumpdT)

errorbar(30,10,yerr=meany,xerr=meanx,linestyle='None',marker='+',markeredgecolor='k', markerfacecolor='k', ecolor='k')

xlim([5,35])
ylim([5,35])

xlabel('Envelope Temp. (K)')
ylabel('Clump Temp. (K)')

T = [5,45]

plt.plot(T,T,'k')

show()
