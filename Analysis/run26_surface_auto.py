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
        print 'Please enter a LOGIC value as "TRUE" or "FALSE"'
    return out

###################################################################
#Start of bulk script
###################################################################

#Produce YSO map
maps = np.loadtxt('maps_surface.txt',dtype='string')

#loop through maps
for i in maps:
    #split INPUT file - #declare map origin
    s850,raO,decO, = string.split(i,',')
    prefix = string.split(s850,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    #deffine region name
    name = string.split(s850,'/')[1]
    region = string.split(name,'_')[0]

    print region

    null = "'YSO/temp/YSOv1.sdf'"
    nullp = "'YSO/temp/protov1.sdf'"

    #Produce null map
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null,0,0,0,0)
    os.system(thresh1)
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,nullp,0,0,0,0)
    os.system(thresh1)

    ### r is deffined here.
    r = 300. #arcsec
    ### r is deffined here

    #create surface map by convolving with gaussian 'r'
    #Load and calculate YSO positions in map
    pos = np.loadtxt('data/c2d+GB_YSOs+DJR.txt',dtype='string') 

    YSO = 'YSO/'+str(region)+'_YSO.sdf'
    surface = 'YSO/'+str(region)+'_YSOsurface_r%s.sdf'%(int(r))
    surfaceunits = 'YSO/'+str(region)+'_YSOsurface_r%s_units.sdf'%(int(r))

    proto = 'YSO/'+str(region)+'_protostar.sdf'

    #declare degrees per pixel per wavelength.
    pix850 = [8.33E-4,-8.33E-4]
    n = 1
    k = 0
    l = 0

    limitra,limitdec = 1,1

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
            if float(j[5]) > -0.3:
            #print x,y,float(j[5])
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

#mv and rename point source map.
    mv1 = 'mv YSO/temp/protov%i.sdf %s'%(k+1,proto)
    os.system(mv1)
#mv and rename point source map.
    mv1 = 'mv YSO/temp/YSOv%i.sdf %s'%(l+1,YSO)
    os.system(mv1)

    print YSO,proto

    sig = FWHM(r,'TRUE')
    p = PARGET(YSO,'fpixscale','ndftrace')
    smooth = sig/p[0]
    
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,YSO,smooth,surface)
    os.system(cmd)

    #scale up units 
    scale850 = 18909
    cmult = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,surface,scale850,surfaceunits)
    os.system(cmult)
    
#clean up
#cmd1 = 'rm -r %s/temp/'%(output_dir)
#os.system(cmd1)

    break
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
        print 'Please enter a LOGIC value as "TRUE" or "FALSE"'
    return out

###################################################################
#Start of bulk script
###################################################################

#Produce YSO map
maps = np.loadtxt('maps_surface.txt',dtype='string')

#loop through maps
for i in maps:
    #split INPUT file - #declare map origin
    s850,raO,decO, = string.split(i,',')
    prefix = string.split(s850,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    #deffine region name
    name = string.split(s850,'/')[1]
    region = string.split(name,'_')[0]

    print region

    null = "'YSO/temp/YSOv1.sdf'"
    nullp = "'YSO/temp/protov1.sdf'"

    #Produce null map
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null,0,0,0,0)
    os.system(thresh1)
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,nullp,0,0,0,0)
    os.system(thresh1)

    ### r is deffined here.
    r = 300. #arcsec
    ### r is deffined here

    #create surface map by convolving with gaussian 'r'
    #Load and calculate YSO positions in map
    pos = np.loadtxt('data/c2d+GB_YSOs+DJR.txt',dtype='string') 

    YSO = 'YSO/'+str(region)+'_YSO.sdf'
    surface = 'YSO/'+str(region)+'_YSOsurface_r%s.sdf'%(int(r))
    surfaceunits = 'YSO/'+str(region)+'_YSOsurface_r%s_units.sdf'%(int(r))

    proto = 'YSO/'+str(region)+'_protostar.sdf'

    #declare degrees per pixel per wavelength.
    pix850 = [8.33E-4,-8.33E-4]
    n = 1
    k = 0
    l = 0

    limitra,limitdec = 1,1

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
            if float(j[5]) > -0.3:
            #print x,y,float(j[5])
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

#mv and rename point source map.
    mv1 = 'mv YSO/temp/protov%i.sdf %s'%(k+1,proto)
    os.system(mv1)
#mv and rename point source map.
    mv1 = 'mv YSO/temp/YSOv%i.sdf %s'%(l+1,YSO)
    os.system(mv1)

    print YSO,proto

    sig = FWHM(r,'TRUE')
    p = PARGET(YSO,'fpixscale','ndftrace')
    smooth = sig/p[0]
    
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,YSO,smooth,surface)
    os.system(cmd)

    #scale up units 
    scale850 = 18909
    cmult = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,surface,scale850,surfaceunits)
    os.system(cmult)
    
#clean up
#cmd1 = 'rm -r %s/temp/'%(output_dir)
#os.system(cmd1)

    break
