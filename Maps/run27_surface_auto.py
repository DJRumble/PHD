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
        print 'Please enter a LOGIC value as "TRUE" or "FALSE"'
    return out

###################################################################
#Start of bulk script
###################################################################

#Produce YSO map
maps = np.loadtxt('maps_central_pix.txt',dtype='string')

#loop through maps
for i in maps:
    #split INPUT file - #declare map origin
    s850,raO,decO, = string.split(i,',')
    prefix = string.split(s850,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    #deffine region name
    name = string.split(s850,'/')[1]
    region = string.split(name,'_20')[0]

    print region

    null = "'YSO/temp/YSOv1.sdf'"
    nullp = "'YSO/temp/protov1.sdf'"
    nullpms = "'YSO/temp/pmsv1.sdf'"

    #Produce null map
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null,0,0,0,0)
    os.system(thresh1)
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,nullp,0,0,0,0)
    os.system(thresh1)
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,nullpms,0,0,0,0)
    os.system(thresh1)

    ### r is deffined here.
    r = 300. #arcsec
    ### r is deffined here

    #create surface map by convolving with gaussian 'r'
    #Load and calculate YSO positions in map
    pos = np.loadtxt('data/GBS_YSO_master.tab',dtype='string') 

    YSO = 'YSO/points_v2/'+str(region)+'_YSO.sdf'
    YSOsurface = 'YSO/surface_v2/'+str(region)+'_YSOsurface_r%s.sdf'%(int(r))
    YSOsurfaceunits = 'YSO/surface_v2/'+str(region)+'_YSOsurface_r%s_units.sdf'%(int(r))
    proto = 'YSO/points_v2/'+str(region)+'_protostar.sdf'
    protosurface = 'YSO/surface_v2/'+str(region)+'_proto_surface_r%s.sdf'%(int(r))
    protosurfaceunits = 'YSO/surface_v2/'+str(region)+'_proto_surface_r%s_units.sdf'%(int(r))
    PMS = 'YSO/points_v2/'+str(region)+'_PMS.sdf'
    PMSsurface = 'YSO/surface_v2/'+str(region)+'_PMS_surface_r%s.sdf'%(int(r))
    PMSsurfaceunits = 'YSO/surface_v2/'+str(region)+'_PMS_surface_r%s_units.sdf'%(int(r))

    #declare degrees per pixel per wavelength.
    #pix850 = [8.33E-4,-8.33E-4]
    pix850 = [1.02083333333E-3,-8.3333333333E-4]
    n = 1
    k = 0
    l = 0
    m = 0

    #limit catalogue to local region to speed up process. 
    limitra,limitdec = 5,5 #units arcmin

    ralo = float(raO)-float(limitra)
    rahi = float(raO)+float(limitra)
    declo = float(decO)-float(limitdec)
    dechi = float(decO)+float(limitdec)
    print ralo, rahi, declo, dechi

    for j in pos:

        if (float(j[2]) > ralo) & (float(j[2]) < rahi) & (float(j[3]) > declo) & (float(j[3]) < dechi):
            #calculate displacement of source from the origin
            delx,dely = (float(raO)-float(j[2])),(float(decO)-float(j[3]))
            #list pixel displacements
            x,y = (delx/pix850[0]),(dely/pix850[1])
            
            #produce YSO map
            inyso = 'YSO/temp/YSOv%i.sdf'%(l+1)
            outyso = 'YSO/temp/YSOv%i.sdf'%(l+2)
            #Implant point sources into null maps
            section = "\"'%i,%i'\""%(int(x),int(y))
            chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,inyso,outyso,section,n)
            os.system(chpix1)  
            #clean up
            cmd1 = 'rm %s'%(inyso)
            os.system(cmd1)
            l = l+1
            #produce protostar map
            if float(j[4]) > -0.3:
                #print x,y
                inyso = 'YSO/temp/protov%i.sdf'%(k+1)
                outyso = 'YSO/temp/protov%i.sdf'%(k+2)
        #Implant point sources into null maps
                section = "\"'%i,%i'\""%(int(x),int(y))
                chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,inyso,outyso,section,n)
                os.system(chpix1) 
                #clean up
                cmd1 = 'rm %s'%(inyso)
                os.system(cmd1) 
                k = k+1

            #produce PMSstar map
            if float(j[4]) < -0.3:
                #print x,y
                inyso = 'YSO/temp/pmsv%i.sdf'%(m+1)
                outyso = 'YSO/temp/pmsv%i.sdf'%(m+2)
        #Implant point sources into null maps
                section = "\"'%i,%i'\""%(int(x),int(y))
                chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,inyso,outyso,section,n)
                os.system(chpix1) 
                #clean up
                cmd1 = 'rm %s'%(inyso)
                os.system(cmd1) 
                m = m+1

#mv and rename point source map.
    mv1 = 'mv YSO/temp/protov%i.sdf %s'%(k+1,proto)
    os.system(mv1)
#mv and rename point source map.
    mv1 = 'mv YSO/temp/YSOv%i.sdf %s'%(l+1,YSO)
    os.system(mv1)
#mv and rename point source map.
    mv1 = 'mv YSO/temp/pmsv%i.sdf %s'%(m+1,PMS)
    os.system(mv1)

    #print YSO,proto,PMS

    sig = FWHM(r,'TRUE')
    p = PARGET(YSO,'fpixscale','ndftrace')
    smooth = sig/p[0]
    
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,YSO,smooth,YSOsurface)
    os.system(cmd)
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,proto,smooth,protosurface)
    os.system(cmd)
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,PMS,smooth,PMSsurface)
    os.system(cmd)

    #scale up units 
    scale850 = 18909
    cmult = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,YSOsurface,scale850,YSOsurfaceunits)
    os.system(cmult)
    cmult = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,protosurface,scale850,protosurfaceunits)
    os.system(cmult)
    cmult = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,PMSsurface,scale850,PMSsurfaceunits)
    os.system(cmult)

    #convert to fits file
    ndf2fits.ndf2fits(PMSsurfaceunits)
    ndf2fits.ndf2fits(protosurfaceunits)
    ndf2fits.ndf2fits(YSOsurfaceunits)
    #print surfaceunits
    
    #clean up
    cmd1 = 'rm  YSO/temp/*'
    #os.system(cmd1)

    #break
