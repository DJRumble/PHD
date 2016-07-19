#################################################################################
########### import modules  #############
#Standard modules
import numpy as np
import astropy.io.fits as pyfits
#import pyfits
import sys
import os
import string
import mpfit
import copy
import re

import ndf2fits

#################################################################################
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

def FWHM(value,LOGIC):
    #LOGIC is 'TRUE' or 'FALSE'.
    #LOGIC == TRUE. Transform FWHM into Sigma
    #LOGIC == FALSE. Transform Sigma into FWHM
    if LOGIC == 'TRUE':
        out = value/(2.*(2.*np.log(2))**0.5)
    elif LOGIC == 'FALSE':
        out = value * (2.*(2.*np.log(2))**0.5)
    else:
        print 'Please enter a LOGIC value as "TRUE" or "False"'
    return out

def FourComponentDualBeam(p450,p850,fwhm450MB,fwhm850MB,fwhm450SB,fwhm850SB,in450,in850,output_dir,file):
#Primary beam convolution using FWHM taken from Demspey2013.
    method = 'beam'
    
    print "primary beam convolution"
    smooth450MB = fwhm850MB/p450
    smooth850MB = fwhm450MB/p850
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in450, smooth450MB, output_dir+'/'+method+'/process/'+file+'/s450convolveMB.sdf')
    os.system(cmd)
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in850, smooth850MB, output_dir+'/'+method+'/process/'+file+'/s850convolveMB.sdf')
    os.system(cmd)

    print "seconary beam convolution"
    smooth450SB = fwhm850SB/p450
    smooth850SB = fwhm450SB/p850
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in450, smooth450SB, output_dir+'/'+method+'/process/'+file+'/s450convolveSB.sdf')
    os.system(cmd)
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in850, smooth850SB, output_dir+'/'+method+'/process/'+file+'/s850convolveSB.sdf')
    os.system(cmd)

    sig450MB = FWHM(fwhm450MB,'TRUE')
    sig850MB = FWHM(fwhm850MB,'TRUE')
    sig450SB = FWHM(fwhm450SB,'TRUE')
    sig850SB = FWHM(fwhm850SB,'TRUE')

#Deffine sigma per pix size for use in normalisation
    SIG_84_MB = sig850MB/p450
    SIG_48_MB = sig450MB/p850
    SIG_84_SB = sig850SB/p450
    SIG_48_SB = sig450SB/p850

 #A (and B) = 1/(2*pi*(a4*(FM4**2)+b4*(FS4**2))) or 1/(2*pi*(a8*(FM8**2)+b8*(FS8**2)))
    A = 1/(2.*np.pi*(alpha450*(SIG_48_MB**2.)+beta450*(SIG_48_SB**2.))) 
    B = 1/(2.*np.pi*(alpha850*(SIG_84_MB**2.)+beta850*(SIG_84_SB**2.))) 

#set numbers for multiplication with convovle maps - in general they are of the form 2pi*(sig/pix)^2*A*alpha
    A_a = 2.*np.pi*alpha450*A*(SIG_48_MB**2.)
    A_b = 2.*np.pi*beta450*A*(SIG_48_SB**2.)
    B_a = 2.*np.pi*alpha850*B*(SIG_84_MB**2.)
    B_b = 2.*np.pi*beta850*B*(SIG_84_SB**2.)

    if (float('%.f'%(A_a+A_b)) == float('%.f'%(B_a+B_b))) == 1.0:
        print 'Normalisation OK'
        print 'A = ',A
        print 'B = ',B
    elif(float('%.f'%(A)) <> float('%.f'%(B))):
        print 'Normalistion FAIL'
        print 'A = ',A
        print 'B = ',B

    print 'printing normalisation constants'
    print '450: ',A_a,'+',A_b,' = ',A_a+A_b
    print '850: ',B_a,'+',B_b,' = ',B_a+B_b

    cmd1 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/'+method+'/process/'+file+'/s450convolveMB.sdf',B_a, output_dir+'/'+method+'/process/'+file+'/s450normMB.sdf')
    cmd2 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/'+method+'/process/'+file+'/s450convolveSB.sdf',B_b, output_dir+'/'+method+'/process/'+file+'/s450normSB.sdf')
    cmd3 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/'+method+'/process/'+file+'/s850convolveMB.sdf',A_a, output_dir+'/'+method+'/process/'+file+'/s850normMB.sdf')
    cmd4 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, output_dir+'/'+method+'/process/'+file+'/s850convolveSB.sdf',A_b, output_dir+'/'+method+'/process/'+file+'/s850normSB.sdf')
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    cmd1 = '%s/add in1=%s in2=%s out=%s'%(kapdir,output_dir+'/'+method+'/process/'+file+'/s450normMB.sdf',output_dir+'/'+method+'/process/'+file+'/s450normSB.sdf', output_dir+'/'+method+'/map/'+file+'/s450convolve.sdf')
    cmd2 = '%s/add in1=%s in2=%s out=%s'%(kapdir,output_dir+'/'+method+'/process/'+file+'/s850normMB.sdf',output_dir+'/'+method+'/process/'+file+'/s850normSB.sdf', output_dir+'/'+method+'/map/'+file+'/s850convolve.sdf')
    os.system(cmd1)
    os.system(cmd2)

s850 = 'temp/temp0.sdf'

#9 sample pixels (from MWC 297)
RaH = [277.06,277.059166667,277.058333333,277.06,277.059166667,277.058333333,277.06,277.059166667,277.058333333]
DecD = [-3.73369444444,-3.73369444444,-3.73369444444,-3.73452777778,-3.73452777778,-3.73452777778,-3.73536111111,-3.73536111111,-3.73536111111]

#Produce grid with 9 pixels on it
for i in range(len(RaH)):
    RaH_i = RaH[i]/15
    DecD_i = DecD[i]
    coords = "%10.7f %s"%(RaH_i, DecD_i)

    cmd = "%s/wcstran ndf=%s posin=\"\'%s\'\" framein=SKY frameout=PIXEL > /dev/null"%(kapdir,s850,coords)
    os.system(cmd)

    posout = (os.popen(kapdir+"/parget posout wcstran").readlines())
    (x,y,junk) = re.split("\s+",posout[0],2)
                
    X = int(round(float(x),1))
    Y = int(round(float(y),1))
            
#Implant point sources into null maps
    section = "\"'%i,%i'\""%(X,Y)

    print section
    inmap = 'temp/temp%s.sdf'%(str(i))
    out = 'temp/temp%s.sdf'%(str(i+1))

    if i == 4:
        n = 1000
    else:
        n = 100

    chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,inmap,out,section,n)
    os.system(chpix1)  



#carry about 2CB convolution
fwhm450MB = 7.9 
fwhm850MB = 13.0 
fwhm450SB = 25.0 
fwhm850SB = 48.0 

alpha450 = 0.94
alpha850 = 0.98
beta450 = 0.06
beta850 = 0.02

output_dir = 'output'
file = '9pix'
method = 'beam'

if not os.path.exists(output_dir+'/'+method):
    os.makedirs(output_dir+'/'+method)
if not os.path.exists(output_dir+'/'+method+'/map'):
    os.makedirs(output_dir+'/'+method+'/map')
if not os.path.exists(output_dir+'/'+method+'/map/'+file):
    os.makedirs(output_dir+'/'+method+'/map/'+file)
if not os.path.exists(output_dir+'/'+method+'/process'):
    os.makedirs(output_dir+'/'+method+'/process')
if not os.path.exists(output_dir+'/'+method+'/process/'+file):
    os.makedirs(output_dir+'/'+method+'/process/'+file)


p850 = PARGET(out,'fpixscale','ndftrace')

FourComponentDualBeam(p850[0],p850[0],fwhm450MB,fwhm850MB,fwhm450SB,fwhm850SB,out,out,output_dir,file)
