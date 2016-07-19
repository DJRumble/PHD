"""Takes inout function and converts from .sdf to .fits file format"""

import os
import string

#################################################################################
########### set CSH Commands directories #############
convdir = '/star/bin/convert'

####################
#funtction to run ndf2fits
def ndf2fits(input):
    """Takes inout function and converts from .sdf to .fits file format"""
   #convert input maps into fits format
    prefix = string.split(input,'.sdf')[0]
    #print 'Converting %s to fits...'%(input)
    fits = prefix+'.fits'
    #print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(fits)):
	   os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s noprohis QUIET'%(convdir,input,fits)
    os.system(cmd)
    return fits

def fits2ndf(input):
    """Takes inout function and converts from .sdf to .fits file format"""
   #convert input maps into fits format
    prefix = string.split(input,'.fits')[0]
    #print 'Converting %s to sdf...'%(input)
    sdf = prefix+'.sdf'
    #print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(sdf)):
	   os.unlink(sdf)
    cmd = '%s/fits2ndf in=%s out=%s QUIET'%(convdir,input,sdf)
    os.system(cmd)
    return sdf

if __name__ == "__main__":
    ndf2fits(input)
    fits2ndf(input)
