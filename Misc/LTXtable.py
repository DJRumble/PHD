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
########## Bulk script ############

data =  np.loadtxt('GBSprotostar_temperatures_master_HEX.tab',dtype='string')

print data[0]

L = len(data[:,0])

SMM = open("JGPC_master_LTX.tab","w")

for i in range(L):
    JGCCID = data[i][0]
    RA = str(data[i][1])
    Dec = str(data[i][2])
    temp_e = str(round(float(data[i][3]),1))+'$\pm$'+str(round(float(data[i][4]),1))
    temp_c = str(round(float(data[i][5]),1))+'$\pm$'+str(round(float(data[i][6]),1))
    a = str(round(float(data[i][7]),1))
    Td = str(round(float(data[i][8]),0))


    SMM.write('%s\t&\t%s\t&\t%s\t&\t%s\t&\t%s\t&\t%s\t&\t%s\taaa\n'%(JGCCID,RA,Dec,temp_e,temp_c,a,Td))
