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
import time

#################################################################################
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'
cupidr = '/star/bin/cupid'

#################################################################################
########## Bulk script #############

maps = np.loadtxt("list.txt",dtype='string')

label = ['With OB stars','Without OB stars','Orion A']

OB = []
noOB = []
OA = []

for i in maps:
    #split INPUT
    temp,OB_LOGIC = string.split(i,',')
    name = string.split(temp,'/')[1]
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
    loc = np.where((data<70.) & (data>8.) & (data!=0.0))
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

n1, bin_edges1, patches1 = hist(OB,bins=50,normed=True, histtype='step',label=str(label[0]),color='r')
n2, bin_edges2, patches2 = hist(noOB,bins=50,normed=True, histtype='step',label=str(label[1]),color='k')
n3, bin_edges3, patches3 = hist(OA,bins=50,normed=True, histtype='step',label=str(label[2]),color='b')

ylabel('Normalised frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([5,70])
ylim([0,0.1])

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_temphist_JCMTGBS_OB.pdf'%(date))

show()

