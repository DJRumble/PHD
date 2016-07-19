#DJR UoE

#30/04/2015 - run25_W40freefree_sub_large.py

#This script is designed to subtract freefree large structure from W40 region SCUBA-2 data. Takes an input of positions, with variable alphas, convolves up to the JCMT beam using the convolution Kernel and then substracts from the original SCUBA-2 data. 

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
cupdir = '/star/bin/cupid'
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

###################################################################
#Start of bulk script
###################################################################

#input maps
s21 = 'maps/s21/aquila_VLAss_21cm_Jypix.sdf'

#s450 = 'maps/comb/Aquila_extS2_450-4-30am_fullfreefree.sdf'
#s850 = 'maps/comb/Aquila_extS2_850-4-30am_fullfreefree.sdf'
s450 = 'maps/Aquila_extS2nosm_s450-4am-fullfreefree.sdf'
s850 = 'maps/Aquila_noco_extS2nosm_s850-4am-fullfreefree.sdf '

#temp maps
sc450 = 'maps/temp/Aquila_S21at450.sdf'
sc850 = 'maps/temp/Aquila_S21at850.sdf'
a450 = 'maps/temp/align_S21at450.sdf'
a850 = 'maps/temp/align_S21at850.sdf'
aref450 = 'maps/temp/align450.sdf'
aref850 = 'maps/temp/align850.sdf'

#final maps
con450 = 'maps/s21/noCO+FF/Aquila_450-4-30am_ats21specs.sdf'
con850 = 'maps/s21/noCO+FF/Aquila_850-4-30am_ats21specs.sdf'
f450 = 'maps/s21/noCO+FF/aquila_VLAss_21cm-5am_at450.sdf'
f850 = 'maps/s21/noCO+FF/aquila_VLAss_21cm-5am_at850.sdf'
sub450 = 'maps/s21/noCO+FF/Aquila_450freefree21.sdf'
sub850 = 'maps/s21/noCO+FF/Aquila_850freefree21.sdf'
mask450 = 'maps/s21/noCO+FF/Aquila_450freefree21_masked.sdf'
mask850 = 'maps/s21/noCO+FF/Aquila_850freefree21_masked.sdf'

#Scale up the 21cm to SCUBA-2 wavelenghts
scalar450 = (21.E-2/450.E-6)**(-0.1)
scalar850 = (21.E-2/850.E-6)**(-0.1)

cmult1 = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,s21,scalar450,sc450)
cmult2 = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,s21,scalar850,sc850)
os.system(cmult1)
os.system(cmult2)

#Align SCUBA-2 maps with the 21cm maps
align450 = '%s/wcsalign in=%s ref=%s out=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,s450,s21,aref450)
align850 = '%s/wcsalign in=%s ref=%s out=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,s850,s21,aref850)
os.system(align450)
os.system(align850)

#convolve the SCUBA-2 maps to VLA resolution
p21 = PARGET(s21,'fpixscale','ndftrace')
FWHM21 = 45
smooth21 = FWHM21/p21[0]

print smooth21

gsmooth1 = '%s/gausmooth in=%s fwhm=%s out=%s'%(kapdir, aref450, smooth21, con450)
gsmooth2 = '%s/gausmooth in=%s fwhm=%s out=%s'%(kapdir, aref850, smooth21, con850)
os.system(gsmooth1)
os.system(gsmooth2)

#remove structure greater than 5' from the 21cm maps
box = (5*60)/p21[0]

findback1 = '%s/findback in=%s out=%s box=%s sub=TRUE RMS=!'%(cupdir,sc450,f450,int(box))
findback2 = '%s/findback in=%s out=%s box=%s sub=TRUE RMS=!'%(cupdir,sc850,f850,int(box))
os.system(findback1)
os.system(findback2)

#sub freefree contamination
sub1 = '%s/sub in1=%s in2=%s out=%s'%(kapdir,con450,f450,sub450)
sub2 = '%s/sub in1=%s in2=%s out=%s'%(kapdir,con850,f850,sub850)
os.system(sub1)
os.system(sub2)

#masking 
mask = 'maps/s21/ratiomask.sdf'

mask1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,sub450,mask,mask450)
mask2 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,sub850,mask,mask850)
os.system(mask1)
os.system(mask2)
