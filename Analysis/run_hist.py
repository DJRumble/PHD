#20160713
#Damian Rumble, UoE

#This is a script produces two panel histograms for OB temperatures.

#### NOTE #####

# To include Orion A in the top plot go to 'list.txt' and uncomment the Orion A line. Extent the histogram limits (line 65) to 99 (likewise for x-axis on both upper and lower plots). 

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

import slfit

#################################################################################
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'
cupidr = '/star/bin/cupid'

#################################################################################
########## Bulk script #############

#PREP DATA UPPER PANEL#

maps = np.loadtxt("list.txt",dtype='string')

label_U = ['With OB stars','Without OB stars','Orion A']

OB = []
noOB = []
OA = []

for i in maps:
    #split INPUT
    temp,OB_LOGIC = string.split(i,',')
    name = string.split(temp,'/')[5]
    region = string.split(name,'_')[0]
    print region
    
    #convert to FITS
    prefix = string.split(temp,'.sdf')[0]
    fits = prefix+'.fits'
    if (os.path.exists(fits)):
	   os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s noprohis QUIET'%(convdir,temp,fits)
    os.system(cmd)

    #open FITS files
    hdulist = pyfits.open(fits)
    data = hdulist[0].data
    #print data
    loc = np.where((data<50.) & (data>8.) & (data!=0.0))
    #print loc
    Data = data[loc]
    #print Data
    
    if OB_LOGIC == 'TRUE':
        for i in Data:
            #print i
            OB.append(i)
    elif OB_LOGIC == 'FALSE':
        for i in Data:
            #print i
            noOB.append(i)
    elif OB_LOGIC == 'ORION':
        for i in Data:
            #print i
            OA.append(i)


print np.median(OB)
print np.median(noOB)
print np.mean(OB)
print np.mean(noOB)

#PREP DATA LOWER PANEL#

#DATA - sets based on stellar type
dataO =  np.loadtxt('data/SMM_master_real_Otype.tab',dtype='string',comments="#")
dataB1 =  np.loadtxt('data/SMM_master_real_eBtype.tab',dtype='string',comments="#")
dataB2 =  np.loadtxt('data/SMM_master_real_lBtype.tab',dtype='string',comments="#")

TempO = [float(i) for i in dataO[:,6]]
ErrtempO = [float(i) for i in dataO[:,7]]
TempB1 = [float(i) for i in dataB1[:,6]]
ErrtempB1 = [float(i) for i in dataB1[:,7]]
TempB2 = [float(i) for i in dataB2[:,6]]
ErrtempB2 = [float(i) for i in dataB2[:,7]]

DistanceO = [np.log10(float(i)) for i in dataO[:,19]]
DistanceB1 = [np.log10(float(i)) for i in dataB1[:,19]]
DistanceB2 = [np.log10(float(i)) for i in dataB2[:,19]]

distanceO = []
distanceB1 = []
distanceB2 = []
distanceOcut = []
distanceB1cut = []
distanceB2cut = []

tempO = []
errtempO = []
tempB1 = []
errtempB1 = []
tempB2 = []
errtempB2 = []
tempOcut = []
errtempOcut = []
tempB1cut = []
errtempB1cut = []
tempB2cut = []
errtempB2cut = []

#Deffine a distance cut point (equivelent to 3pc)
Cut = 0.48

for i in range(len(DistanceO)):
    D = DistanceO[i]
    if (D>=-1.5) & (D<=1.5):
        distanceO.append(D)
        tempO.append(TempO[i])
        errtempO.append(ErrtempO[i])
    if (D>=-1.5) & (D<=Cut): #trim for better fit
        distanceOcut.append(D)
        tempOcut.append(TempO[i])
        errtempOcut.append(ErrtempO[i])
for i in range(len(DistanceB1)):
    D = DistanceB1[i]
    if (D>=-1.5) & (D<=1.5): #general
        distanceB1.append(D)
        tempB1.append(TempB1[i])
        errtempB1.append(ErrtempB1[i])
    if (D>=-1.5) & (D<=Cut): #trim for better fit
        distanceB1cut.append(D)
        tempB1cut.append(TempB1[i])
        errtempB1cut.append(ErrtempB1[i])
for i in range(len(DistanceB2)):
    D = DistanceB2[i]
    if (D>=-1.5) & (D<=1.5):#general
        distanceB2.append(D)
        tempB2.append(TempB2[i])
        errtempB2.append(ErrtempB2[i])
    if (D>=-1.5) & (D<=Cut): #trim for better fit
        distanceB2cut.append(D)
        tempB2cut.append(TempB2[i])
        errtempB2cut.append(ErrtempB2[i])

label_L = ['O type','Early B type','Late B type']

###PLOT FIGURE####
fig1 = plt.figure()

### UPPER ###
fig1.add_subplot(211)

n1, bin_edges1, patches1 = hist(OB,bins=50,normed=True, histtype='step',label=str(label_U[0]),color='r')
n2, bin_edges2, patches2 = hist(noOB,bins=50,normed=True, histtype='step',label=str(label_U[1]),color='k')
#n3, bin_edges3, patches3 = hist(OA,bins=50,normed=True, histtype='step',label=str(label[2]),color='b')

ylabel('Normalised frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([5,50])
ylim([0,0.1])

### LOWER ###
fig1.add_subplot(212)

nO, bin_edgesO, patchesO = hist(tempO,bins=10,normed=True, histtype='step',label=str(label_L[0]),color='r')
nB1, bin_edgesB1, patchesB1 = hist(tempB1,bins=10,normed=True, histtype='step',label=str(label_L[1]),color='y')
nB2, bin_edgesB2, patchesB2 = hist(tempB2,bins=10,normed=True, histtype='step',label=str(label_L[2]),color='b')

D,P = stats.ks_2samp(nO,nB1)
print 'KS-stats: D = ',D,' P = ',P
D,P = stats.ks_2samp(nB1,nB2)
print 'KS-stats: D = ',D,' P = ',P

ylabel('Normalised frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([5,50])
ylim([0,0.125])

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_temphist_2xHIST_JCMTGBS.pdf'%(date))

show()

