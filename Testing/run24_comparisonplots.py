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
import ndf2fits
from scipy.stats import norm

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

def K(T): #rat

    beta = 1.8
    #print (np.exp((h*c)/(k*850*1E-6*10))-1.) / (np.exp((h*c)/(k*450*1E-6*10))-1.)*(850.**(3.+beta)) / (450**(3.+beta))
    #print (np.exp((h*c)/(k*850*1E-6*20))-1.) / (np.exp((h*c)/(k*450*1E-6*20))-1.)*(850.**(3.+beta)) / (450**(3.+beta))
   

    l2 = 1200
    
    return (np.exp((h*c)/(k*850*1E-6*20))-1.) / (np.exp((h*c)/(k*450*1E-6*20))-1.)*(850.**(3.+beta)) / (450**(3.+beta))


########### Bulk Code #############
#Method A = RUMBLE
#Method B = CHEN

inA_fits = 'PerseusWest_20150317-autotemperatureWCSALN_SMT.fits'
errA_fits = 'PerseusWest_20150317-autotemp_errorWCSALN_SMT.fits'

inB_fits = 'tempMap_Her+850_clean.fits'
errB_fits = 'TempMapErr_Her+850_clean.fits'

Beta_fits = 'betaMap_Her+850_clean.fits'
errBeta_fits = 'betaMapErr_Her+850_clean.fits'

inA = ndf2fits.fits2ndf(inA_fits)
inB =  ndf2fits.fits2ndf(inB_fits)
errA = ndf2fits.fits2ndf(errA_fits)
errB =  ndf2fits.fits2ndf(errB_fits)

k = 0

TA = []
TB = []
EA = []
EB = []
Beta = []
eBeta = []

count_lo = 0
count_hi = 0
count_lower = 0
count_upper = 0
count_lower_sub = 0
count_upper_sub = 0
print 'a'
imageA = pyfits.open(str(inA_fits), mode = 'update') #RUMBLE
imageB = pyfits.open(str(inB_fits), mode = 'update') #CHEN
errorA = pyfits.open(str(errA_fits), mode = 'update') #RUMBLE
errorB = pyfits.open(str(errB_fits), mode = 'update') #CHEN
beta = pyfits.open(str(Beta_fits), mode = 'update') #RUMBLE

#print beta

errbeta = pyfits.open(str(errBeta_fits), mode = 'update') #CHEN
print 'b'
dimA = PARGET(inA,'dims','ndftrace')
dimB = PARGET(inB,'dims','ndftrace')
print 'c'
rowA = int(dimA[1])    #rows of the table from NDFTRACE
columnA = int(dimA[0])    #columns of the table from NDFTRACE
rowB = int(dimB[1])    #rows of the table from NDFTRACE
columnB = int(dimB[0])    #columns of the table from NDFTRACE

if (rowA == rowB) & (columnA == columnB):
    print 'Maps are the same region: GOOD'
    print 'row = ' + str(rowA) +':'+ str(rowB)
    print 'column =' + str(columnA) +':'+ str(columnB)
else:
    print 'Maps are not the same region: BAD'
    print 'row = ' + str(rowA) +':'+ str(rowB)
    print 'column =' + str(columnA) +':'+ str(columnB)
    sys.exit()


dataA = imageA[0].data
dataB = imageB[0].data
dataEA = errorA[0].data
dataEB = errorB[0].data
dataBeta = beta[0].data
dataeBeta = errbeta[0].data

print 'imageA =' + str(imageA)
print 'imageB =' + str(imageB)

for i in range(0, rowA):
#print i
    for j in range(0, columnA):
#print j
        if (dataA[i][j] > 0.) & (dataA[i][j] < 120.5) & (dataB[i][j] > 0.) & (dataB[i][j] < 120.5):
            TA.append(dataA[i][j]) #RUMBLE
            TB.append(dataB[i][j]) #CHEN
            EA.append(dataEA[i][j]/2.)
            EB.append(dataEB[i][j]/2.) 
            Beta.append(dataBeta[i][j]) #CHEN
            eBeta.append(dataeBeta[i][j]) #CHEN

            Tr = dataA[i][j]/dataB[i][j]

            if Tr > 1.:
                count_lo = count_lo + 1
            if Tr <= 1.:
                count_hi = count_hi + 1





l = len(TA)

print 'CHEN > RUMBLE: ',count_lo
print 'RUMBLE < CHEN: ',count_hi
print 'TOTAL: ',count_lo+count_hi

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

fig = plt.figure(1, figsize=(16, 8))

fig.add_subplot(121)

R16err = []
C16err = []
Berr = []
Bsub = []
Tsub = []



#print EA

for i in range(len(EA)):
    if Beta[i] > 2.1:
        count_lower = count_lower + 1
    if Beta[i] < 1.5:
        count_upper = count_upper + 1
    if str(EA[i]) <> 'nan':
        R16err.append(EA[i])
        C16err.append(EB[i])
        Berr.append(eBeta[i])
        Bsub.append(Beta[i])
        Tsub.append(TB[i])

        if Beta[i] < 1.5:
            count_lower_sub = count_lower_sub + 1
        if Beta[i] > 2.1:
            count_upper_sub = count_upper_sub + 1

print 'count upper: ',str(count_upper_sub),' out of ',str(count_upper)
print 'count lower: ',str(count_lower_sub),' out of ',str(count_lower)

meanx = np.mean(R16err)
meany = np.mean(C16err)
meanB = np.mean(Berr)

print meanB

plt.scatter(TA,TB,s=1)
plt.plot([0,100],[0,100],'r-',lw=2)
errorbar(15,40,yerr=meany,xerr=meanx,linestyle='None',marker='+',markeredgecolor='k', markerfacecolor='k', ecolor='k')
plt.xlabel('Temperature by SCUBA-2 flux ratio(K)')
plt.ylabel('Temperature by SED fitting (K)')
ylim([5,50])
xlim([5,50])

fig.add_subplot(122)

#hist,bin_edges = np.histogram(Beta,bins=20, normed=True)
#mu, std = norm.fit(Beta)

#p = norm.pdf(bin_edges, mu, std)

plt.scatter(Beta,TB,s=1,color='r')
plt.scatter(Bsub,Tsub,s=1,color='k')
errorbar(3,40,yerr=meany,xerr=meanB,linestyle='None',marker='+',markeredgecolor='k', markerfacecolor='k', ecolor='k')
plt.xlabel(r'$\beta$ by SED fitting (K)')
plt.ylabel('Temperature by SED fitting (K)')

#plt2 = plt.twinx()
#
#plt2.plot(bin_edges, p, 'b', linewidth=2)

#print p
#n1,b1,p1 = hist(Beta,(15),histtype='step',color='r', normed=True)
axvline(x=1.8,color='k',linewidth=3,linestyle='-')
axvline(x=1.5,color='k',linewidth=1,linestyle='--')
axvline(x=2.1,color='k',linewidth=1,linestyle='--')

ylim([5,50])
xlim([1,3.5])

#,histtype='step'

#plt.scatter(onpointeight,Beta,s=1)
#plt2.set_ylabel('Frequency')





#if (os.path.exists(inA_fits)):
#    os.unlink(inA_fits)
#if (os.path.exists(inB_fits)):
#    os.unlink(inB_fits)
#if (os.path.exists(erA_fits)):
#    os.unlink(erA_fits)
#if (os.path.exists(erB_fits)):
#    os.unlink(erB_fits)

plt.show()


