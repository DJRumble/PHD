#20141203
#Damian Rumble, UoE

#This script is a prototype designed to extract pixels from a temperature map and display them as a histogram. 


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

#################################################################################
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

#################################################################################
########## Bulk script #############

def ndf2fits(file):

#convert input maps into fits format
    prefix = string.split(file,'.sdf')[0]
    print 'Converting %s to fits...'%(file)
    fits = prefix+'.fits'
#print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(fits)):
        os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s QUIET'%(convdir,file,fits)
    os.system(cmd)
    return fits

def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)
    a = np.loadtxt('parameter.txt')
    os.remove('parameter.txt')
    return a


#################################################################################
########## Bulk script #############

in_fileA = str(sys.argv[1])
#in_fileB = str(sys.argv[2])

fits_fileA = ndf2fits(in_fileA)
#fits_fileB = ndf2fits(in_fileB)

numpixA = PARGET(in_fileA,'NUMGOOD','stats')
#numpixB = PARGET(in_fileB,'NUMPIX','stats')

print numpixA
#print numpixB

#Open file, extract data, build histogram, fit model using MPFIT
#open data and variance
hdulistA = pyfits.open(fits_fileA)
#hdulistB = pyfits.open(fits_fileB)
dataA = hdulistA[0].data
#dataB = hdulistB[0].data

locA = np.where((dataA!=0.0)&(dataA>5)&(dataA<999))
#locB = np.where((dataB!=0.0)&(dataB>5)&(dataB<999))

bins=100

histA, bin_edgesA = np.histogram(dataA[locA], bins = bins)
#histB, bin_edgesB = np.histogram(dataB[locB], bins = bins)

#normalised hist
Xa = histA/numpixA
#Xb = histB/numpixB

print "Plotting..."

#fig1 = plt.figure(1)
#fig1.set_size_inches((2,10))

#fig1.subplots_adjust(hspace=0.001)

#nx=1,ny=5;

#fig1.clf()

#j = 0

plt.clf()

plt.scatter(bin_edgesA[:-1], Xa,c=u'r')
#plt.scatter(bin_edgesB[:-1], Xb,c=u'b')

plt.xlabel('Temperature (K)')
plt.ylabel('Normalised Number of pixels')

plt.ylim([0,0.1])
#plt.ylim([0,300])
plt.xlim([5,60])

plt.text(30,0.08,str(in_fileA),size='small')#,relative='True')

plt.show()



if (os.path.exists(fits_fileA)):
    os.unlink(fits_fileA)
#if (os.path.exists(fits_fileB)):
#    os.unlink(fits_fileB)
