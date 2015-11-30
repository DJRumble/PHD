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

def hist(file):

    fits_fileA = ndf2fits(file)

    numpixA = PARGET(file,'NUMPIX','stats')
    print numpixA

#Open file, extract data, build histogram, fit model using MPFIT
#open data and variance
    hdulistA = pyfits.open(fits_fileA)
    dataA = hdulistA[0].data
    locA = np.where((dataA!=0.0)&(dataA>5)&(dataA<999))
    bins=100
    histA, bin_edgesA = np.histogram(dataA[locA], bins = bins)

#normalised hist
    Xa = histA/numpixA

    return bin_edgesA,Xa


#################################################################################
########## Bulk script #############

fileA = 'Aquila_extmask_s2temperature_SerpSouth.sdf'
fileB = 'Aquila_extmask_s2temperature_W40.sdf'
fileC = 'SerpensMain_IR1_Kerneltemperature_NH2.sdf'
fileD = 'SerpensMWC297_IR1_freefree+jettemperature_section.sdf'
fileE = 'SerpensMain_IR1_Kerneltemperature_main.sdf'

files = [fileA,fileB,fileC,fileD,fileE]
colour = ['r','b','g','y','k']

print "Plotting..."

fig1 = plt.figure(1)
fig1.set_size_inches((2,10))

fig1.subplots_adjust(hspace=0.001)

nx=1;ny=5;

fig1.clf()

i = 0
j = 0

for i in range(len(files)):
    plt = fig1.add_subplot(ny,nx,j+1)
    histogram = hist(files[i])
    bin_edges=histogram[0]
    x = histogram[1]
    C = colour

    plt.scatter(bin_edges[:-1],x,c=u'r')

    axes=fig1.gca()
    
    plt.set_xticklabels('',visible='False',size='small')
    plt.set_yticklabels('',visible='True',size='small')
    
    if (j==4):
        plt.set_xlabel(r'Temperature (K)',size='small')
        plt.set_xticklabels('',visible='True',size='small')
        
    plt.axis([5,60,0,1.1E-3])
    
    plt.set_ylabel('Normalised Number of pixels',size='small')

    j = j + 1

fig1.show()
