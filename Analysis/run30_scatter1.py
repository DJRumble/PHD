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

def RJ(T):
    l = 850E-6
    nu = c / l
    #equation of a blackbody function
    B = (2*(nu**2.)*k*T)/(c**2.)
    return B

def CD(T):
    S = 10E-27 #Jy
    M = 0.5*M_x #kg
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

#No 15K clumps
datanoOBnoYSO =  np.loadtxt('data/SMM_master_real_noOBheated_noYSO.tab',dtype='string',comments="#") 
datanoOBYSO =  np.loadtxt('data/SMM_master_real_noOBheated_YSO.tab',dtype='string',comments="#")
dataOBnoYSO =  np.loadtxt('data/SMM_master_OBheated_noYSO.tab',dtype='string',comments="#")
dataOBYSO = np.loadtxt('data/SMM_master_OBheated_YSO.tab',dtype='string',comments="#")

dataOBnoYSO15K =  np.loadtxt('data/SMM_master_OBheated_noYSO+15Kclumps.tab',dtype='string',comments="#")
dataOBYSO15K = np.loadtxt('data/SMM_master_OBheated_YSO+15Kclumps.tab',dtype='string',comments="#")



#All clumps
#ALLdatanoOBnoYSO =  np.loadtxt('data/SMM_masterNO-OB-noYSO_full_flt.tab',dtype='float',comments="#") 
#ALLdatanoOBYSO =  np.loadtxt('data/SMM_masterNO-OB-YSO_full_flt.tab',dtype='float',comments="#") 
#ALLdataOBnoYSO =  np.loadtxt('data/SMM_masterOB-noYSO_full_flt.tab',dtype='float',comments="#")
#ALLdataOBYSO =  np.loadtxt('data/SMM_masterOB-YSO_full_flt.tab',dtype='float',comments="#")

#### KEY ####

#1 = noOBnoYSO
#2 = noOBYSO
#3 = OBnoYSO
#4 = OBYSO

#Radius
##R1 = np.array(ALLdatanoOBnoYSO[:,10]).tolist() 
##R2 = np.array(ALLdatanoOBYSO[:,10]).tolist() 
#R3 = np.array(ALLdataOBnoYSO[:,10]).tolist() 
#R4 = np.array(ALLdataOBYSO[:,10]).tolist() 

#MASS
#M1 = map(float, data[:,4])
#dM1 = map(float, data[:,4])
#M2 = map(float, data[:,4])
#dM2 = map(float, data[:,4])
#M3 = np.array(dataOBnoYSO[:,3]).tolist() 
#dM3 = np.array(dataOBnoYSO[:,4]).tolist() 
#M4 = np.array(dataOBYSO[:,3]).tolist() 
#dM4 = np.array(dataOBYSO[:,4]).tolist() 
#allM1 = np.array(ALLdatanoOBnoYSO[:,3]).tolist() 
#alldM1 = np.array(ALLdatanoOBnoYSO[:,4]).tolist() 
#allM2 = np.array(ALLdatanoOBYSO[:,3]).tolist() 
#alldM2 = np.array(ALLdatanoOBYSO[:,4]).tolist() 
#allM3 = np.array(ALLdataOBnoYSO[:,3]).tolist() 
#alldM3 = np.array(ALLdataOBnoYSO[:,4]).tolist() 
#allM4 = np.array(ALLdataOBYSO[:,3]).tolist() 
#alldM4 = np.array(ALLdataOBYSO[:,4]).tolist() 

#TEMP
T1 = map(float, datanoOBnoYSO[:,6])
dT1 = map(float, datanoOBnoYSO[:,7])
T2 = map(float, datanoOBYSO[:,6])
dT2 = map(float, datanoOBYSO[:,7])
T3 = map(float, dataOBnoYSO[:,6])
dT3 = map(float, dataOBnoYSO[:,7])
T4 = map(float, dataOBYSO[:,6])
dT4 = map(float, dataOBYSO[:,7])

#CD
CD1 = map(float, datanoOBnoYSO[:,9])
dCD1 = map(float, datanoOBnoYSO[:,10])
CD2 = map(float, datanoOBYSO[:,9])
dCD2 = map(float, datanoOBYSO[:,10])
CD3 = map(float, dataOBnoYSO[:,9])
dCD3 = map(float, dataOBnoYSO[:,10])
CD4 = map(float, dataOBYSO[:,9])
dCD4 = map(float, dataOBYSO[:,10])
CD5 = map(float, dataOBnoYSO15K[:,30])
dCD5 = map(float, dataOBnoYSO15K[:,31])
CD6 = map(float, dataOBYSO15K[:,30])
dCD6 = map(float, dataOBYSO15K[:,31])

#MJ
J1 = map(float, datanoOBnoYSO[:,17])
dJ1 = map(float, datanoOBnoYSO[:,18])
J2 = map(float, datanoOBYSO[:,17])
dJ2 = map(float, datanoOBYSO[:,18])
J3 = map(float, dataOBnoYSO[:,17])
dJ3 = map(float, dataOBnoYSO[:,18])
J4 = map(float, dataOBYSO[:,17])
dJ4 = map(float, dataOBYSO[:,18])
J5 = map(float, dataOBnoYSO15K[:,38])
dJ5 = map(float, dataOBnoYSO15K[:,39])
J6 = map(float, dataOBYSO15K[:,38])
dJ6 = map(float, dataOBYSO15K[:,39])

#### the scatter plots: ####

fig = plt.figure(1, figsize=(16, 8))
'''
### CD vs T ###
fig.add_subplot(121)
errorbar(T1,CD1,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
errorbar(T2,CD2,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w')
errorbar(T3,CD3,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='None', ecolor='w')
errorbar(T4,CD4,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w')

xlim([5,45])
ylim([5,2500])

xlabel('Temp. (K)')
ylabel('Peak column density (10$^{21}$ H$_{2}$ cm$^{-2}$)')

T = [5,45]
T = np.arange(5,45,0.1)

N = [CD(i) for i in T]

plt.plot(T,N,'magenta')

yscale('log')

axhline(y=20,color='k',linestyle='--')
axhline(y=200,color='k',linestyle='-')

#axvline(x=14.8,color='k',linestyle='--')
#axvline(x=17.2,color='k',linestyle=':')
#axvline(x=12.4,color='k',linestyle=':')

### MJ vs T ###
fig.add_subplot(122)

errorbar(T1,J1,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
errorbar(T2,J2,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w')
errorbar(T3,J3,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='None', ecolor='w')
errorbar(T4,J4,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w')

xlim((6,45))
ylim((0.05,11))

T = [5,45]
T = np.arange(5,45,0.1)

J = [stability1(i) for i in T]

plt.plot(T,J,'magenta')

yscale('log')

axhline(y=1,color='k',linestyle='--')
axhline(y=2,color='k',linestyle='-')

#axvline(x=14.8,color='k',linestyle='--')
#axvline(x=17.2,color='k',linestyle=':')
#axvline(x=12.4,color='k',linestyle=':')

xlabel('Temp. (K)')
ylabel('M/M$_{\mathrm{J}}$')

### CD vs Stability ###
'''
plt.figure(1, figsize=(8, 8))

#errorbar(CD1,J1,xerr=0,yerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
#errorbar(CD2,J2,xerr=0,yerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w')
errorbar(CD3,J3,xerr=0,yerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='None', ecolor='w')
errorbar(CD4,J4,xerr=0,yerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w')
errorbar(CD5,J5,xerr=0,yerr=0,linestyle='None',marker='o',markeredgecolor='k', markerfacecolor='None', ecolor='w')
errorbar(CD6,J6,xerr=0,yerr=0,linestyle='None',marker='o',markeredgecolor='k', markerfacecolor='k', ecolor='w')

#print CD5
#print J5

CD_pre = CD1 + CD3
CD_pro = CD2 + CD4
CD = CD1 + CD3 + CD2 + CD4

bins = np.logspace(start=0.75,stop=3.75,num=10,base=10)

#n, bin_edges, patches = hist(CD_pre,bins,cumulative=0)

#print 'PRE'
#print n
#print bin_edges
#n, bin_edges, patches = hist(CD_pro,bins,cumulative=0)

#print 'PRO'
#print n
#print bin_edges
#n, bin_edges, patches = hist(CD,bins,cumulative=0)

#print 'TOTAL'
#print n
#print bin_edges

mCD1 = np.mean(CD1)
mCD2 = np.mean(CD2)
mCD3 = np.mean(CD3)
mCD4 = np.mean(CD4)
mCD5 = (mCD1+mCD2)/2.
mCD6 = (mCD3+mCD4)/2.
mJ1 = np.mean(J1)
mJ2 = np.mean(J2)
mJ3 = np.mean(J3)
mJ4 = np.mean(J4)
#print mJ1

#errorbar(150,0.05,xerr=mCD5,yerr=0,linestyle='None',marker='+',markeredgecolor='b', markerfacecolor='None', ecolor='b')
#errorbar(700,0.05,xerr=mCD6,yerr=0,linestyle='None',marker='+',markeredgecolor='r', markerfacecolor='None', ecolor='r')


ylim((0.03,13))
xlim([5,3000])

xlabel('Peak column density (10$^{21}$ H$_{2}$ cm$^{-2}$)')
ylabel('M/M$_{\mathrm{J}}$')

#N = [5,3000]
N = np.arange(5,3000,1)

J = [stability2(i) for i in N]

#print J

plt.plot(N,J,'magenta')

yscale('log')
xscale('log')

axhline(y=1,color='k',linestyle='--')
axhline(y=2,color='k',linestyle='-')
axvline(x=20,color='k',linestyle='--')
axvline(x=200,color='k',linestyle='-')

#'''
show()


#Archive

### MJ vs M ###
#errorbar(M1,MJ1,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
#errorbar(M2,MJ2,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w')
#errorbar(M3,MJ3,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='None', ecolor='w')
#errorbar(M4,MJ4,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w')

#yscale('log')

#xlim((-0.5,90))
#ylim((0,11))

#axScatter2.axhline(y=1,color='k',linestyle='--')
#axScatter2.axhline(y=2,color='k',linestyle='-')

#xlabel('Mass (M$_{\odot}$)')
#ylabel('M/M$_{\mathrm{J}}$')

### M vs CD ###
#errorbar(M1,CD1,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='None', ecolor='w')
#errorbar(M2,CD2,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='b', markerfacecolor='b', ecolor='w')
#errorbar(M3,CD3,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='None', ecolor='w')
#errorbar(M4,CD4,yerr=0,xerr=0,linestyle='None',marker='o',markeredgecolor='r', markerfacecolor='r', ecolor='w')

#yscale('log')
#xscale('log')

#xlim((-0.5,90))
#ylim((0,2500))

#xlabel('Mass (M$_{\odot}$)')
#ylabel('Peak column density (10$^{21}$ H$_{2}$ cm$^{-2}$)')

