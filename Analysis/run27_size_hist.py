#20160218version1/
#Damian Rumble, UoE

#this script produces plots of clump properties

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

#No 15K clumps
data =  np.loadtxt('SMM/SMM_master.tab',dtype='string')

#### KEY ####

#1 = noOBnoYSO
#2 = noOBYSO
#3 = OBnoYSO
#4 = OBYSO


size = [float(i) for i in data[:,11]]

#Radius
R1 = np.array(size).tolist() 
#R2 = np.array(data[:,10]).tolist() 

######## BULK SCRIPT ########

print np.median(size)
print np.median(R1)


#ALL
n1,b1,p1 = hist(R1,(36),histtype='bar',color='k',rwidth=0.8,normed=False)
#n2,b2,p2 = hist(OBT,(10),histtype='step',color='r',normed=True)

xlabel('Flux weighted clump diameter (pc)')
ylabel('Frequency')

#D,P = stats.ks_2samp(n1,n2)
#print 'KS-stats: D = ',D,' P = ',P

#label1 = 'P value = '+str(round(P,2))
#label2 = 'D value = '+str(round(D,2))

#ylim([0,0.18])

#annotate(label1,xy=(32,0.16))
#annotate(label2,xy=(32,0.15))

show()


