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

#my modules
import mass

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

def ColumnD(S,T,kappa,a,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), Number of pixels in distance (d in pcs)"""
    #constants
    au = 149597871000 #m
    M_x = 1.989E30 #kg
    m_h2 = 3.34745E-24 #g
    mu = 2.333333  #ratio of H2 to He (5:1)
    N = 1

    #each pixel has an area in SI units and each apature contains N pixels
    A = ((((a*d)*au)**2.0)*N)*10000 #in cm
    M = mass.djr(S,T,kappa,d)

    m = M*M_x*1000.0 #A = 1 density in g per pixel

    n = m/(mu*m_h2*A) #column density per cm^2 
    N = m/A #density in g per cm^2
    #print 'mass per pixel in solar masses:', M
    #print 'column density of apature in per cm^2:', n
    #print 'density of apature in g per cm^2:', N
    return N

#### Bulk Code ####

#deffine input maps
s850FITS = 'input/Aquila_extS2_850-4-30am_NM.fits'
tempFITS = 'input/Aquila-4-30amtemperature.fits'
CDFITS = 'input/Aquila-4-30amColumnD.fits'
CDsdf = 'input/Aquila-4-30amColumnD.sdf'
CDatoms = 'input/Aquila-4-30amColumnDatoms.sdf'
Av = 'input/Aquila-4-30am_Av.sdf'

###VLA files - manually enter pixel size of 15
#s850FITS = 'input/s21/Aquila_850freefree21_masked.fits'
#tempFITS = 'input/s21/Aquila_FFtempVLA.fits'
#CDFITS = 'input/s21/Aquila-FF-ColumnD.fits'
#CDsdf = 'input/s21/Aquila-FF-ColumnD.sdf'
#CDatoms = 'input/s21/Aquila-FF-ColumnDatoms.sdf'
#Av = 'input/s21/Aquila-FF-Av.sdf'

#deffine constants
kappa = 0.012 #cm^2 g^-1
d = 500 #pc

#open maps
image_s850 = pyfits.open(str(s850FITS)) #opens FITS file created from NDF
image_temp = pyfits.open(str(tempFITS))
image_CD = pyfits.open(str(CDFITS), mode = 'update')

s850 = image_s850[0].data  #flux data in a table
temp = image_temp[0].data  #temp data in a table
CD = image_CD[0].data  #CD data in a table

#Calculate dim numbers from ndftrace 
dim = PARGET(s850FITS,'dims','ndftrace')
#find pixel size
#p850 = PARGET(s850FITS,'fpixscale','ndftrace')

row = int(dim[1])    #rows of the table from NDFTRACE
column = int(dim[0]) #columns of the table from NDFTRACE  
   
print 'row = ' + str(row) 
print 'column =' + str(column)   
print 'temp image =' + str(image_temp)
print 's850 image =' + str(image_s850)

#pass through table, completing the Column density caculation
for i in range(0, row):
    for j in range(0, column):
        S = float(s850[i][j])
        T = float(temp[i][j])
        #CD[i][j] = S/T
        CD[i][j]  = ColumnD(S,T,kappa,15,d)
        #print CD[i][j]
  
#output new values back to FITS file and save
image_CD.flush() 
image_CD.close() 
cmd = '%s/fits2ndf in=%s out=%s'%(convdir,CDFITS,CDsdf)
os.system(cmd)

