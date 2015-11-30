#DJR UoE

#27/03/2015 - run25_W40freefree_sub.py

#This script is designed to subtract freefree point sources from W40 region SCUBA-2 data. Takes an input of positions, with variable alphas, convolves up to the JCMT beam using the convolution Kernel and then substracts from the original SCUBA-2 data. 

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
    return fits

def FourComponentDualBeam(p450,p850,fwhm450MB,fwhm850MB,fwhm450SB,fwhm850SB,in450,in850):
    #THIS FUNCTION TURNS a 450 point into a 450 beam (and vv)

#Primary beam convolution using FWHM taken from Demspey2013.
    print "primary beam convolution"
    smooth450MB = fwhm450MB/p450
    smooth850MB = fwhm850MB/p850

    print smooth450MB

    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in450, smooth450MB,'maps/temp/s450convolveMB.sdf')
    os.system(cmd)
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in850, smooth850MB, 'maps/temp/s850convolveMB.sdf')
    os.system(cmd)

    print "seconary beam convolution"
    smooth450SB = fwhm450SB/p450
    smooth850SB = fwhm850SB/p850
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in450, smooth450SB, 'maps/temp/s450convolveSB.sdf')
    os.system(cmd)
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir, in850, smooth850SB, 'maps/temp/s850convolveSB.sdf')
    os.system(cmd)

    sig450MB = FWHM(fwhm450MB,'TRUE')
    sig850MB = FWHM(fwhm850MB,'TRUE')
    sig450SB = FWHM(fwhm450SB,'TRUE')
    sig850SB = FWHM(fwhm850SB,'TRUE')

#Deffine sigma per pix size for use in normalisation
    SIG_8_MB = sig850MB/p850
    SIG_4_MB = sig450MB/p450
    SIG_8_SB = sig850SB/p850
    SIG_4_SB = sig450SB/p450

 #A (and B) = 1/(2*pi*(a4*(FM4**2)+b4*(FS4**2))) or 1/(2*pi*(a8*(FM8**2)+b8*(FS8**2)))
    A = 1/(2.*np.pi*(alpha450*(SIG_4_MB**2.)+beta450*(SIG_4_SB**2.))) 
    B = 1/(2.*np.pi*(alpha850*(SIG_8_MB**2.)+beta850*(SIG_8_SB**2.))) 

#set numbers for multiplication with convovle maps - in general they are of the form 2pi*(sig/pix)^2*A*alpha
    A_a = 2.*np.pi*alpha450*A*(SIG_4_MB**2.)
    A_b = 2.*np.pi*beta450*A*(SIG_4_SB**2.)
    B_a = 2.*np.pi*alpha850*B*(SIG_8_MB**2.)
    B_b = 2.*np.pi*beta850*B*(SIG_8_SB**2.)

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

    cmd1 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, 'maps/temp/s450convolveMB.sdf',B_a, 'maps/temp/s450normMB.sdf')
    cmd2 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, 'maps/temp/s450convolveSB.sdf',B_b, 'maps/temp/s450normSB.sdf')
    cmd3 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, 'maps/temp/s850convolveMB.sdf',A_a, 'maps/temp/s850normMB.sdf')
    cmd4 = '%s/cmult in=%s scalar=%s out=%s'%(kapdir, 'maps/temp/s850convolveSB.sdf',A_b, 'maps/temp/s850normSB.sdf')
    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    os.system(cmd4)
    cmd1 = '%s/add in1=%s in2=%s out=%s'%(kapdir,'maps/temp/s450normMB.sdf','maps/temp/s450normSB.sdf', 'maps/s450freefree_comb.sdf')
    cmd2 = '%s/add in1=%s in2=%s out=%s'%(kapdir,'maps/temp/s850normMB.sdf','maps/temp/s850normSB.sdf', 'maps/s850freefree_comb.sdf')
    os.system(cmd1)
    os.system(cmd2)


###################################################################
#Start of bulk script
###################################################################
################ 1 -- Create  free-free point sources  #############

#s450 = '/data/damian/maps/serpens/MWC297/IR2/SerpensMWC297_20141219_s450_IR2extmask_s2_cal_mJysqaJH.sdf'
s450 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_extS2_450-4-30am.sdf'
#s850 = '/data/damian/maps/serpens/MWC297/IR2/SerpensMWC297_20141219_s850_IR2extmask_s2_cal_mJysqaJH.sdf'
s850 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_extS2_850-4-30am.sdf'

#make directories
if not os.path.exists('maps/temp'):
    os.makedirs('maps/temp')

#prepare null maps
s450collapse = 'maps/temp/collapse450.sdf'
s850collapse = 'maps/temp/collapse850.sdf'

null450 = "'maps/temp/map450v1.sdf'"
null850 = "'maps/temp/map850v1.sdf'"


thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s450,null450,0,0,0,0)
thresh2 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null850,0,0,0,0)
os.system(thresh1)
os.system(thresh2)

#Calculate sections
pos = np.loadtxt('W40_3_6cm_UCHII.tab') #0,1,3
#pos = np.loadtxt('UCHIIc_FW_W40.tab') #1,2,11

#declare map origin
origx,origy = 277.894708333,-1.90202777778
#calculate displacement of source from the origin
delx,dely = (origx-pos[:,0]),(origy-pos[:,1])
#declare degrees per pixel per wavelength.
pix450 = [5.55E-4,-5.55E-4]
pix850 = [8.33E-4,-8.33E-4]
#list pixel displacements
x450,y450 = (delx/pix450[0]),(dely/pix450[1])
x850,y850 = (delx/pix850[0]),(dely/pix850[1])

for i in range (len(x450)):
    #deffine the source flux at 3.6cm
    s3_6 = pos[i][3]
    alpha = pos[i][7]
    nu_crit = pos[i][8]
    #scale up the flux to SCUBA-2 wavelengths
    s450pix = nu_crit*s3_6*((666.67/8.33)**alpha)
    s850pix = nu_crit*s3_6*((352.9/8.33)**alpha)

    print '3.6cm, 850um, 450um flux = ',s3_6, s850pix, s450pix

    in450 = 'maps/temp/map450v%i.sdf'%(i+1)
    out450 = 'maps/temp/map450v%i.sdf'%(i+2)
    in850 = 'maps/temp/map850v%i.sdf'%(i+1)
    out850 = 'maps/temp/map850v%i.sdf'%(i+2)

    #Implant point sources into null maps
    section = "\"'%i,%i'\""%(int(x450[i]+1),int(y450[i]))
    chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f '%(kapdir,in450,out450,section,s450pix)
    section = "\"'%i,%i'\""%(int(x850[i]+1),int(y850[i]))
    chpix2 = '%s/chpix in=%s out=%s section=%s newval=%f '%(kapdir,in850,out850,section,s850pix)
    os.system(chpix1)
    os.system(chpix2)

#mv and rename point source map.
mv1 = 'mv maps/temp/map450v%s.sdf maps/pointsourcemap450comb.sdf'%(i+2)
mv2 = 'mv maps/temp/map850v%s.sdf maps/pointsourcemap850comb.sdf'%(i+2)
os.system(mv1)
os.system(mv2)

########### Fixed variables
#Set FWHM for SCUBA-2 maps - MB (Main Beam) S (Secondary Beam)
fwhm450MB = 7.9 
fwhm850MB = 13.0 
fwhm450SB = 25.0 
fwhm850SB = 48.0 

#set normalisation constants
alpha450 = 0.94
alpha850 = 0.98
beta450 = 0.06
beta850 = 0.02

in450 = 'maps/pointsourcemap450comb.sdf'
in850 = 'maps/pointsourcemap850comb.sdf'

p450 = PARGET(in450,'fpixscale','ndftrace')
p850 = PARGET(in850,'fpixscale','ndftrace')

print p450[0],p850[0]

#Convolve each point source map with the JCMT beam at that respective wavelength.
FourComponentDualBeam(p450[0],p850[0],fwhm450MB,fwhm850MB,fwhm450SB,fwhm850SB,in450,in850)

sub1='%s/sub in1=%s in2=%s out=%s'%(kapdir,s450,'maps/s450freefree_comb.sdf','maps/Aquila_extS2_450-4-30am_fullfreefree.sdf')
sub2='%s/sub in1=%s in2=%s out=%s'%(kapdir,s850,'maps/s850freefree_comb.sdf','maps/Aquila_extS2_850-4-30am_fullfreefree.sdf')
os.system(sub1)
os.system(sub2)
