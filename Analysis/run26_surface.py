#DRJ UoE
#20150521 - run26_surface.py

#this script will take a blank map of a JCMT GBS region and a catalogue of YSOs and produce a surface density plot of some resolution 'r'.

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

###################################################################
#Start of bulk script
###################################################################

#Produce YSO map

#Start with 450 (2'' pixel map)
s450 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_extS2_450-4-30am.sdf'
s850 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_extS2_850-4-30am.sdf'

#if not os.path.exists('/temp'):
#    os.makedirs('/temp')

null = "'temp/YSOv1.sdf'"

#Produce null map
thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null,0,0,0,0)
os.system(thresh1)

A = 1
### r is deffined here.
r = 300. #arcsec
### r is deffined here

#create surface map by convolving with gaussian 'r'
#Load and calculate YSO positions in map
if A == 1:
    pos = np.loadtxt('W40_YSOs_biglist.tab') 
    YSO = 'output/YSOmap850.sdf'
    surface = 'output/Aquila850_YSOsurface_r%s.sdf'%(int(r))
    surfaceunits = 'output/Aquila850_YSOsurface_r%s_units.sdf'%(int(r))
elif A == 2:
    pos = np.loadtxt('W40_YSOs_protostars.tab') 
    YSO = 'output/protostar850map.sdf'
    surface = 'output/Aquila850_protostarsurface_r%s.sdf'%(int(r))
    surfaceunits = 'output/Aquila850_protostarsurface_r%s_units.sdf'%(int(r))
elif A == 3:
    pos = np.loadtxt('W40_YSOs_PMS.tab') 
    YSO = 'output/PMS850map.sdf'
    surface = 'output/Aquila850_PMSsurface_r%s.sdf'%(int(r))
    surfaceunits = 'output/Aquila850_PMSsurface_r%s.sdf'%(int(r))
else:
    sys.exit()

#declare map origin
#origx,origy = 277.894708333,-1.90202777778 #450
origx,origy = 277.895,-1.90147222222 #850

#calculate displacement of source from the origin
delx,dely = (origx-pos[:,0]),(origy-pos[:,1])
#declare degrees per pixel per wavelength.
#pix450 = [5.55E-4,-5.55E-4]
pix850 = [8.33E-4,-8.33E-4]
#list pixel displacements
x,y = (delx/pix850[0]),(dely/pix850[1])

for i in range (len(x)):

    print 'Loop ',i,' : YSO ',int(x[i])+1,':',int(y[i]),' class:',pos[i][2] 
    n = 1

    inyso = 'temp/YSOv%i.sdf'%(i+1)
    outyso = 'temp/YSOv%i.sdf'%(i+2)

    #Implant point sources into null maps
    section = "\"'%i,%i'\""%(int(x[i]+1),int(y[i]))
    chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,inyso,outyso,section,n)
    os.system(chpix1)

#mv and rename point source map.
mv1 = 'mv temp/YSOv%i.sdf %s'%(i+2,YSO)
os.system(mv1)

sig = FWHM(r,'TRUE')
p = PARGET(null,'fpixscale','ndftrace')
smooth = sig/p[0]

cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,YSO,smooth,surface)
os.system(cmd)

#scale up units 
scale450 = 42545
scale850 = 18909

cmult = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,surface,scale850,surfaceunits)
os.system(cmult)

#clean up
#cmd1 = 'rm -r %s/temp/'%(output_dir)
#os.system(cmd1)
