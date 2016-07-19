## ColumnD_maps.py
#DJR UoE 20150512

#this script takes a map of flux and temperature and produces a map of column density

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
import shutil

#my modules
import mass
import ndf2fits 
import noise

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

def massF(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""
    #Equation of Mass in solar masses, per pixel
    if T == 0:
        M = 0
    else:
        M = 1.55*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/500.0)**2.0)
    return M


#### Bulk Code ####

#Load maps from master file, order: s850,temp,tempr_err,alpha,clumps,clumpcat
maps = np.loadtxt('maps4fluxextractor_master.txt',dtype='string')
OBstars = np.loadtxt('OBstars.txt',dtype='string')

#deffine constants
kappa = 0.012 #cm^2 g^-1

k = 0

#loop through maps
for i in maps:
    #Make directory
    if not os.path.exists('CDmaps/CDpart'):
        os.makedirs('CDmaps/CDpart')

    #split INPUT file
    s850,temp,dtemp,alpha,proto,YSO,clumps,clumpcat = string.split(i,',')
    prefix1 = string.split(s850,'.sdf')[0]
    prefix2 = string.split(temp,'.sdf')[0]
    file = string.split(prefix1,'/')[0]
    #deffine region name
    name = string.split(clumps,'/')[2]
    region = string.split(name,'_2')[0]

    #Convert maps to fits
    ndf2fits.ndf2fits(s850)
    ndf2fits.ndf2fits(temp)
    s850FITS = prefix1+'.fits'
    tempFITS = prefix2+'.fits'
    massFITS = 'CDmaps/mass/'+region+'_mass.fits'
    mass15FITS = 'CDmaps/mass/'+region+'_mass15K.fits'

    #print massFITS
    print temp

    cmd = 'cp %s %s'%(s850FITS,massFITS)
    os.system(cmd)
    cmd = 'cp %s %s'%(s850FITS,mass15FITS)
    os.system(cmd)

    #Open files
    image_s850 = pyfits.open(str(s850FITS)) 
    image_temp = pyfits.open(str(tempFITS))
    image_mass = pyfits.open(str(massFITS), mode = 'update')
    image_mass15 = pyfits.open(str(mass15FITS), mode = 'update')

    print 'flux open'
    s850_data = image_s850[0].data  #flux data in a table
    print 'temp open'
    temp_data = image_temp[0].data  #temp data in a table
    print 'mass open'
    Mass_data = image_mass[0].data  #CD data in a table
    Mass15_data = image_mass15[0].data  #CD data in a table

    #Calculate dim numbers from ndftrace 
    dim = PARGET(s850,'dims','ndftrace')
    
    #find pixel size
    p850 = PARGET(s850,'fpixscale','ndftrace')

    a = p850[0]

    row = int(dim[1])    #rows of the table from NDFTRACE
    column = int(dim[0]) #columns of the table from NDFTRACE  
   
    print 'row = ' + str(row) 
    print 'column = ' + str(column)   

    #DISTANCE
    d = float(OBstars[k][4])
    print 'Distance = ',d,' pc'
    k = k + 1

    #pass through table, completing the Column density caculation
    for i in range(0, row):
        for j in range(0, column):
            T = float(temp_data[i][j])
            S = float(s850_data[i][j])
            Mass_data[i][j]  = massF(S,T,kappa,d)
            Mass15_data[i][j]  = massF(S,20,kappa,d)
        #if Mass[i][j] > 0:
            #print Mass[i][j]
  
    #output new values back to FITS file and save
    image_mass.flush() 
    image_mass.close() 
    image_mass15.flush() 
    image_mass15.close() 

    ndf2fits.fits2ndf(massFITS)
    ndf2fits.fits2ndf(mass15FITS)

    #CALCULATE Column density
    #Convert map of mass (in solar masses) into column density (in H2 cm-2) and Extinction (mag)
    M_x = 1.989E33 #g
    N = 1
    au = 14959787100000 #cm
    m_h = 1.67262178E-24 #g
    mu = 2.8 #333  #ratio of H2 to He (Kauffmann et al. 2008)

    #pixel area
    A = ((((a*d)*au)**2.0)*N) #in cm^2 
    f = M_x/(mu*m_h*A) #column density per cm^2 

    mass = 'CDmaps/mass/'+region+'_mass.sdf'
    CD = 'CDmaps/CDpart/'+region+'_CD.sdf'    
    mass15 = 'CDmaps/mass/'+region+'_mass15K.sdf'
    CD15 = 'CDmaps/CDpart/'+region+'_CD15K.sdf'    

    print 'Making CD maps'
    cmd = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,mass,f,CD)
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,mass15,f,CD15)
    os.system(cmd)

    #### MAKE COMPOSIT MAPS ####

    CDmagic = 'CDmaps/CDpart/'+region+'_CDmagic.sdf'
    CD15magic = 'CDmaps/CDpart/'+region+'_CD15Kmagic.sdf'
    CDmask = 'CDmaps/CDpart/'+region+'_CD_MSK.sdf'
    CD15mask = 'CDmaps/CDpart/'+region+'_CD15K_MSK.sdf' 
    maskbad = 'CDmaps/CDpart/'+region+'_MSKbad.sdf'

    sigma = noise.noise_by_data(s850,'FALSE')

    SIGMA = 3.*sigma

    cmd = '%s/nomagic in=%s out=%s repval=0 QUIET'%(kapdir,CD,CDmagic)
    os.system(cmd)
    cmd = '%s/nomagic in=%s out=%s repval=0 QUIET'%(kapdir,CD15,CD15magic)
    os.system(cmd)

    cmd = '%s/thresh in=%s out=%s thrlo=0 thrhi=0 newlo=0 newhi=1 QUIET'%(kapdir,CDmagic,CDmask)
    os.system(cmd)
    cmd = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=0 newhi=1 QUIET'%(kapdir,s850,CD15mask,SIGMA,SIGMA)
    os.system(cmd)
    cmd = '%s/thresh in=%s out=%s thrlo=%f thrhi=%f newlo=%s newhi=1 QUIET'%(kapdir,s850,maskbad,SIGMA,SIGMA,'bad')
    os.system(cmd)

    mask = 'CDmaps/CDpart/'+region+'_MSK.sdf'

    cmd = '%s/sub in1=%s in2=%s out=%s'%(kapdir,CD15mask,CDmask,mask)
    os.system(cmd)

    CDpart = 'CDmaps/CDpart/'+region+'_CDpart.sdf' 

    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,CD15,mask,CDpart)
    os.system(cmd)

    CDcomp = 'CDmaps/'+region+'_CDcomb.sdf' 
    CDcomp0 = 'CDmaps/CDpart/'+region+'_CDcomb.sdf' 

    cmd = '%s/add in1=%s in2=%s out=%s'%(kapdir,CDmagic,CDpart,CDcomp0)
    os.system(cmd)

    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,CDcomp0,maskbad,CDcomp)
    os.system(cmd)

    print CDcomp
    print '================================================'

    ndf2fits.ndf2fits(CDcomp)

    shutil.rmtree('CDmaps/CDpart/')

    break


