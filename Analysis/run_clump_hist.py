#13/07/2016

#This script produces histrograms of clump populations 

import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt

import slfit

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

label = ['O type','Early B type','Late B type']

fig1 = plt.figure()

nO, bin_edgesO, patchesO = hist(tempO,bins=10,normed=True, histtype='step',label=str(label[0]),color='r')
nB1, bin_edgesB1, patchesB1 = hist(tempB1,bins=10,normed=True, histtype='step',label=str(label[1]),color='y')
nB2, bin_edgesB2, patchesB2 = hist(tempB2,bins=10,normed=True, histtype='step',label=str(label[2]),color='b')

D,P = stats.ks_2samp(nO,nB1)
print 'KS-stats: D = ',D,' P = ',P
D,P = stats.ks_2samp(nB1,nB2)
print 'KS-stats: D = ',D,' P = ',P

ylabel('Normalised frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([5,47])
ylim([0,0.125])

show()
