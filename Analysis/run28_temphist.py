#20151002
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

#################################################################################
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'
cupidr = '/star/bin/cupid'

#################################################################################
########## Bulk script #############

#list.txt = W40
#list2.txt = Serpens Aquila

list = np.loadtxt("list.txt",dtype='string')

print list

label = ['Dust Arc','W40-N','W40-S']
#label = ['W40','SerpS','SerpE','SerpMain','SerpNH3','MWC297','SerpN']
l = len(label)

colours = ['k','r','lime']
j = 0

for i in range(l):
    input = list[i]
    prefix = string.split(input,'.sdf')[0]
    #print 'Converting %s to fits...'%(input)
    fits = prefix+'.fits'
    #print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(fits)):
	   os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s noprohis QUIET'%(convdir,input,fits)
    os.system(cmd)

    hdulist = pyfits.open(fits)
    data = hdulist[0].data

    loc = np.where((data<70.) & (data>0.) & (data!=0.0))
    Data = data[loc]
    bins=50

    n, bin_edges, patches  = hist(Data, bins,normed=True, histtype='step',label=str(label[i]),color=[colours[j]])
    j = j +1

ylabel('Frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([0,75])
ylim([0,0.2])

show()
