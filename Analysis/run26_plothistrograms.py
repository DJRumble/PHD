#20150630
#Damian Rumble, UoE

#this script produces histograms of clump properties

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

#################################################################################
########## Bulk script ############

data =  np.loadtxt('SMM/SMM-Aq-ff-W40.tab',dtype='float',comments="#")
data15 =  np.loadtxt('SMM/SMM-Aq-ff-W4015.tab',dtype='float',comments="#")

fig1 = plt.figure()

#mass
fig1.add_subplot(221)
xlabel('Mass (M$_{\odot}$)')
ylabel('Norm. frequency')
#label.set_fontsize('x-small')
DATA = np.array(data[:,6]).tolist()
DATA15 = np.array(data15[:,6]).tolist()
bins = 80
n, bins, patches = hist(DATA, bins, normed=True,label=('mass'),color='r')
xlim([0,65])
n, bins, patches = hist(DATA15, bins, normed=True,label=('mass-15'),color='b')
xlim([0,65])
#CD
fig1.add_subplot(222)
xlabel('Column Density (10$^{25}$ H$_{2}$ cm$^{-2}$)')
ylabel('Norm. frequency')
bins = 80
DATA = np.array(data[:,10]).tolist()
DATA15 = np.array(data15[:,10]).tolist()
n, bins, patches = hist(DATA, bins, normed=True,label=('CD'),color='r')
n, bins, patches = hist(DATA15, bins, normed=True,label=('CD-15'),color='b')
xlim([0,3.3])
#Temp
fig1.add_subplot(223)
xlabel('Temp. (K)')
ylabel('Norm. frequency')
bins = 30
DATA = np.array(data[:,8]).tolist()
DATA15 = np.array(data15[:,8]).tolist()
n, bins, patches = hist(DATA, bins, normed=True,label=('Temp.'),color='r')
n, bins, patches = hist(DATA15, bins, normed=True,label=('Temp.-15'),color='b')
xlim([5,29])
#M/Mj
fig1.add_subplot(224)
xlabel('Jeans stability')
ylabel('Norm. frequency')
bins = 40
DATA = np.array(data[:,16]).tolist()
DATA15 = np.array(data15[:,16]).tolist()
n, bins, patches = hist(DATA, bins, normed=True,label=('Jeans.'),color='r')
n, bins, patches = hist(DATA15, bins, normed=True,label=('Jeans.-15'),color='b')
xlim([0,7.8])
axvline(x=1,c='k',linestyle='dashed')



#legend()
show()
