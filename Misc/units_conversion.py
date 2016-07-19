import os
import numpy as np
import string
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
import mpfit
import copy
import ndf2fits

convdir = '/stardev/bin/convert'
kapdir = '/stardev/bin/kappa'

def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)

    a = np.loadtxt('parameter.txt')
    os.remove('parameter.txt')
    return a

files = np.loadtxt('maps.txt',dtype='string',comments="#")

wave = '850'


for myfile in files:
    if '850' in myfile:
        print "=============="
        print '850 micron map'
        print myfile
        p850 = PARGET(myfile,'fpixscale','ndftrace')
        print "Pixel sizes are "+str(p850[0])
        p850 = p850[0]
        cunits = (p850**2.)/1000
    elif '450' in myfile:
        print "=============="
        print '450 micron map'
        print myfile
        p450 = PARGET(myfile,'fpixscale','ndftrace')
    #p850 = PARGET(input_dir+'/'+file850+'.sdf','fpixscale','ndftrace')
        print "Pixel sizes are "+str(p450[0])
        p450 = p450[0]
        cunits = (p450**2.)/1000

    old_units = '.sdf'
    nu_units = 'Jypix.sdf'
    col_units = 'Jypix_col.sdf'

    prefix = string.split(myfile,old_units)[0]
    print 'Converting %s to %s...'%(old_units,nu_units)

    nu_file = prefix+nu_units
    nu_2dfile = prefix+col_units

    cmd = '%s/cmult in=%s scalar=%f out=%s QUIET'%(kapdir,myfile,cunits,nu_file)
    os.system(cmd)

    cmd = '%s/setunits ndf=%s units=%s'%(kapdir,nu_file,'Jy/pix')
    os.system(cmd)

    cmd = '%s/collapse in=%s out=%s axis=%s'%(kapdir,nu_file,nu_2dfile,3)
    os.system(cmd)
    
    ndf2fits.ndf2fits(nu_2dfile)
