
import os
import numpy as np
import string
from astropy import units as u 
from astropy.coordinates import SkyCoord 

kapdir = '/star/bin/kappa'

def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)

    a = np.loadtxt('parameter.txt',dtype='string')
    os.remove('parameter.txt')
    return a

maps = np.loadtxt('maps_surface.txt',dtype='string')

file = open("maps_central_pix.txt","w")
file = open("maps_central_pix.txt","a")

#loop through maps
for i in maps:
    #split INPUT file
    s850 = string.split(i,',')[0]
    #deffine region name
    name = string.split(s850,'/')[1]
    region = string.split(name,'_20')[0]

    out = 'input_s850/center_pix/'+str(region)+'_centerpix.sdf'
    section = "\"'%i,%i'\""%(int(0),int(0))
    newval = 10000


    cmd = '%s/chpix in=%s out=%s section=%s newval=%s'%(kapdir,s850,out,section,newval)
    os.system(cmd)

    print '====='
    print region
    
    pix_wcs =  PARGET(out,'MAXWCS','stats')

    pix_wcs[0] = string.split(pix_wcs[0],',')[0]
    pix_wcs[1] = string.split(pix_wcs[1],',')[0]

    string_wcs = str(pix_wcs[0])+':'+str(pix_wcs[1])
    
    RAh,RAm,RAs,Dech,Decm,Decs = string.split(string_wcs,':')

    string_wcs = RAh+' '+RAm+' '+RAs+' '+Dech+' '+Decm+' '+Decs
    
    print string_wcs
    
    c = SkyCoord('%s'%(string_wcs), unit=(u.hourangle, u.deg))
    print c.ra.degree, c.dec.degree

    file.write(str(s850)+','+str(c.ra.degree)+','+str(c.dec.degree)+'\n')
    
    

    
    
