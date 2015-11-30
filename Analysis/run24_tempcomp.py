#Damian Rumble, UoE
#20150106
#TvT.py

#this script takes two temp maps for the same region made using different methods and plots Ta vs Tb to compare the differences (and errors). 

####################################################
#import maths, ploting and astrophysical packages

import random
import numpy as np
import os
import subprocess
import commands 
import astropy.io.fits as pyfits
import string
from pylab import *
import matplotlib.pyplot as plt
from scipy import stats

########### set CSH Commands directories #############

kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

################################################################################
########### Functions
def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)

    a = np.loadtxt('parameter.txt')
    os.remove('parameter.txt')
    return a

def NDF2FITS(input):
    #convert input maps into fits format
    prefix = string.split(input,'.sdf')[0]
    print 'Converting %s to fits...'%(input)
    fits = prefix+'.fits'
    #print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(fits)):
	   os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s QUIET'%(convdir,input,fits)
    os.system(cmd)

########### Bulk Code #############
#Method A = Dual Beam
#Method B = Kernel

inA = 'ratio/Aquila_IR2_s2_DB_Sratio_mask.sdf' 
inA_fits = 'ratio/Aquila_IR2_s2_DB_Sratio_mask.fits' 
inB = 'ratio/Aquila_IR2_s2_K_Sratio_mask.sdf' 
inB_fits = 'ratio/Aquila_IR2_s2_K_Sratio_mask.fits'
inC = 'ratio/Aquila-4-30am_SratioBLmask.sdf'
inC_fits = 'ratio/Aquila-4-30am_SratioBLmask.fits'

erA = 'ratio/Aquila_IR2_s2_DBSratio_error_mask.sdf'
erA_fits = 'ratio/Aquila_IR2_s2_DBSratio_error_mask.fits'
erB = 'ratio/Aquila_IR2_s2_KSratio_error_mask.sdf'
erB_fits = 'ratio/Aquila_IR2_s2_KSratio_error_mask.fits'
erC = 'ratio/Aquila-4-30amSratioBL_errormask.sdf'
erC_fits = 'ratio/Aquila-4-30amSratioBL_errormask.fits'

'''
inA = 'ratio/mwc297_DBtest_Sratio_SMM1.sdf' 
inA_fits = 'ratio/mwc297_DBtest_Sratio_SMM1.fits'
inB = 'ratio/mwc297_Ktest_Sratio_SMM1.sdf' 
inB_fits = 'ratio/mwc297_Ktest_Sratio_SMM1.fits'

erA = 'ratio/mwc297_DBtestSratio_error_SMM1.sdf'
erA_fits = 'ratio/mwc297_DBtestSratio_error_SMM1.fits'
erB = 'ratio/mwc297_KtestSratio_error_SMM1.sdf'
erB_fits = 'ratio/mwc297_KtestSratio_error_SMM1.fits'

inA = 'ratio/mwc297_DBtest_Sratio.sdf' 
inA_fits = 'ratio/mwc297_DBtest_Sratio.fits'
inB = 'ratio/mwc297_Ktest_Sratio.sdf' 
inB_fits = 'ratio/mwc297_Ktest_Sratio.fits'

erA = 'ratio/mwc297_DBtestSratio_error.sdf'
erA_fits = 'ratio/mwc297_DBtestSratio_error.fits'
erB = 'ratio/mwc297_KtestSratio_error.sdf'
erB_fits = 'ratio/mwc297_KtestSratio_error.fits'

inA = 'temp/mwc297_DBtesttemperature.sdf' 
inA_fits = 'temp/mwc297_DBtesttemperature.fits'
inB = 'temp/mwc297_Ktesttemperature.sdf' 
inB_fits = 'temp/mwc297_Ktesttemperature.fits'

erA = 'temp/mwc297_DBtesttemp_error.sdf'
erA_fits = 'temp/mwc297_DBtesttemp_error.fits'
erB = 'temp/mwc297_Ktesttemp_error.sdf'
erB_fits = 'temp/mwc297_Ktesttemp_error.fits'
'''

if (os.path.exists(inC_fits)):
    os.unlink(inC_fits)
if (os.path.exists(inB_fits)):
    os.unlink(inB_fits)
if (os.path.exists(erC_fits)):
    os.unlink(erC_fits)
if (os.path.exists(erB_fits)):
    os.unlink(erB_fits)

NDF2FITS(inC)
NDF2FITS(inB)
NDF2FITS(erC)
NDF2FITS(erB)

imageA = pyfits.open(str(inC_fits), mode = 'update')
imageB = pyfits.open(str(inB_fits), mode = 'update')
errorA = pyfits.open(str(erC_fits), mode = 'update')
errorB = pyfits.open(str(erB_fits), mode = 'update')

dimA = PARGET(inA,'dims','ndftrace')
dimB = PARGET(inB,'dims','ndftrace')

rowA = int(dimA[1])    #rows of the table from NDFTRACE
columnA = int(dimA[0])    #columns of the table from NDFTRACE
rowB = int(dimB[1])    #rows of the table from NDFTRACE
columnB = int(dimB[0])    #columns of the table from NDFTRACE

if (rowA == rowB) & (columnA == columnB):
    print 'Maps are the same region: GOOD'
else:
    print 'Maps are not the same region: BAD'
    sys.exit()


dataA = imageA[0].data
dataB = imageB[0].data
errordataA = errorA[0].data
errordataB = errorB[0].data

print 'row = ' + str(rowA)
print 'column =' + str(columnA)
print 'imageA =' + str(imageA)
print 'imageB =' + str(imageB)

TA = []
TB = []
EA = []
EB = []

for i in range(0, rowA):
    #print i
    for j in range(0, columnA):
        #print j
        if (dataA[i][j] > 0.) & (dataB[i][j] > 0.):
            TA.append(dataA[i][j])
            TB.append(dataB[i][j])
            EA.append(errordataA[i][j]/2.)
            EB.append(errordataB[i][j]/2.) 

l = len(TA)

for i in range(l):
    #print '=============================='
    #print 'Wavelength start ', data[i][3]
    xi = TA 
    yi = TB
    
    A = array([xi, ones(9)])

    slope, intercept, r_value, p_value, slope_std_err = stats.linregress(xi,yi)

print 'intercept', intercept
print 'slope', slope
print 'standard deviation', slope_std_err
print 'R value',r_value
print 'P value',p_value

X = np.arange(0,20,1)

line = (slope*X)+intercept

plt.scatter(TA,TB,s=1)
plot(X,line,'-')
#plt.errorbar(TA,TB,xerr=EA,yerr=EB,fmt='+')
plt.plot([0,20],[0,20],'r-',lw=2)
#plt.xlabel('Ratio by Dual Beam Method')
#plt.ylabel('Ratio by Kernel Method')
plt.xlabel('Ratio by bilinear method')
plt.ylabel('Ratio by nearest Method')
ylim([3,18])
xlim([3,18])

if (os.path.exists(inA_fits)):
    os.unlink(inA_fits)
if (os.path.exists(inB_fits)):
    os.unlink(inB_fits)
if (os.path.exists(erA_fits)):
    os.unlink(erA_fits)
if (os.path.exists(erB_fits)):
    os.unlink(erB_fits)

plt.show()


