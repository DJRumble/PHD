#Damian Rumble, UoE
#20141215
#error-testing.py

#this script takes a flat test map of a constant flu and creates a new maps where the pixels are distributed about a mean by a given STDV.

####################################################
#import maths, ploting and astrophysical packages

import random
import numpy as np
import os
import subprocess
import commands 
import astropy.io.fits as pyfits
import string

########### set CSH Commands directories #############

kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

#################################################################################
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
### distribution parameters ### 

#Asign new value to each pixel on the new map - FITS
#convert to FITS

#450
in450 = 'tests/test_80_450.sdf'

mean450 = PARGET(in450,'MEAN','stats')
STDV450 = 0.0165

print 'Map has a mean flux of '+str(mean450)+' and STDV of '+str(STDV450)

NDF2FITS(in450)

fits450 = 'tests/test_80_450.fits'

dim450 = PARGET(in450,'dims','ndftrace')

image450 = pyfits.open(str(fits450), mode = 'update')

row450 = int(dim450[1])    #rows of the table from NDFTRACE
column450 = int(dim450[0])    #rows of the table from NDFTRACE

data450 = image450[0].data

print data450[0][1][0]

print 'row = ' + str(row450)
print 'column =' + str(column450)
print 'image =' + str(image450)

for i in range(0, row450):
    #print i
    for j in range(0, column450):
        #print j
        newval = random.normalvariate(mean450,STDV450)
        #print newval
        #print data450[0][i][j]
        data450[0][i][j] = newval
        #print data450[0][i][j]

image450.flush()
image450.close()

#850

in850 = 'tests/smm11_test_s850.sdf'

mean850 = PARGET(in850,'MEAN','stats')

STDV850 = 0.0022

print 'Map has a mean flux of '+str(mean850)+' and STDV of '+str(STDV850)

NDF2FITS(in850)

fits850 = 'tests/smm11_test_s850.fits'

dim850 = PARGET(in850,'dims','ndftrace')

image850 = pyfits.open(str(fits850), mode = 'update')

row850 = int(dim850[1])    #rows of the table from NDFTRACE
column850 = int(dim850[0])    #rows of the table from NDFTRACE

data850 = image850[0].data

print data850[0][1][0]

print 'row = ' + str(row850)
print 'column =' + str(column850)
print 'image =' + str(image850)

for i in range(0, row850):
    #print i
    for j in range(0, column850):
        #print j
        newval = random.normalvariate(mean850,STDV850)
        #print newval
        #print data850[0][i][j]
        data850[0][i][j] = newval
        #print data850[0][i][j]

image850.flush()
image850.close()
