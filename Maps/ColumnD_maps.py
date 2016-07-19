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

def massF(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""
    #Equation of Mass in solar masses, per pixel
    if T == 0:
        M = 0
    else:
        M = 1.55*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/500.0)**2.0)
    return M

#### Bulk Code ####

#deffine input maps
#s850FITS = 'input/S2-FF-CO/Aquila_noco_extS2nosm_s850-4am-fullfreefreeColNM.fits'
#s850 = 'input/S2-FF-CO/Aquila_noco_extS2nosm_s850-4am-fullfreefreeColNM.sdf'
s850FITS = 'input/OphSco_Main_20150318_850_IR2_ext_JypixHKCOL.fits'
s850 = 'input/OphSco_Main_20150318_850_IR2_ext_JypixHKCOL.sdf'

#tempFITS = 'input/S2-FF-CO/Aquila_noco_s2nosm-4am-FFtemperature.fits'
#tempsdf = 'input/S2-FF-CO/Aquila_noco_s2nosm-4am-FFtemperature.sdf'
tempFITS = 'input/OphSco_Main-autotemperature.fits'
tempsdf = 'input/OphSco_Main-autotemperature.sdf'
#tempFITS = 'input/S2-FF-CO/T15.fits'
#tempsdf = 'input/S2-FF-CO/T15.sdf'

massFITS = 'input/OphSco_Main-automass.fits'
masssdf = 'input/OphSco_Main-automass.sdf'
#massFITS = 'input/S2-FF-CO/Aquila-noco4am-massT15.fits'
#masssdf = 'input/S2-FF-CO/Aquila-noco4am-massT15.sdf'

#CDsdf = 'input/S2-FF-CO/Aquila-noco4am-ColumnD.sdf'
#Av = 'input/S2-FF-CO/Aquila-noco4am_Av.sdf'
#CDsdf = 'input/S2-FF-CO/Aquila-noco4am-ColumnDT15.sdf'
#Av = 'input/S2-FF-CO/Aquila-noco4am_AvT15.sdf'
CDsdf = 'input/OphSco_Main-autoCD.sdf'
Av = 'input/OphSco_Main-autoAv.sdf'

###VLA files - manually enter pixel size of 15
#s850FITS = 'input/s21/Aquila_850freefree21_masked.fits'
#tempFITS = 'input/s21/Aquila_FFtempVLA.fits'
#CDFITS = 'input/s21/Aquila-FF-ColumnD.fits'
#CDsdf = 'input/s21/Aquila-FF-ColumnD.sdf'
#CDatoms = 'input/s21/Aquila-FF-ColumnDatoms.sdf'
#Av = 'input/s21/Aquila-FF-Av.sdf'

#deffine constants
kappa = 0.012 #cm^2 g^-1
d = 139 #pc 500

#open maps
image_s850 = pyfits.open(str(s850FITS)) #opens FITS file created from NDF
image_temp = pyfits.open(str(tempFITS))
image_mass = pyfits.open(str(massFITS), mode = 'update')

s850 = image_s850[0].data  #flux data in a table
temp = image_temp[0].data  #temp data in a table
Mass = image_mass[0].data  #CD data in a table

#Calculate dim numbers from ndftrace 
dim = PARGET(tempsdf,'dims','ndftrace')
#find pixel size
p850 = PARGET(tempsdf,'fpixscale','ndftrace')
print str(p850[0])

a = p850[0]

row = int(dim[1])    #rows of the table from NDFTRACE
column = int(dim[0]) #columns of the table from NDFTRACE  
   
print 'row = ' + str(row) 
print 'column =' + str(column)   
print 'temp image =' + str(image_temp)
print 's850 image =' + str(image_s850)


#pass through table, completing the Column density caculation
for i in range(0, row):
    for j in range(0, column):
        T = float(temp[i][j])
        S = float(s850[i][j])
        Mass[i][j]  = massF(S,T,kappa,d)
        #if Mass[i][j] > 0:
            #print Mass[i][j]
  
#output new values back to FITS file and save
image_mass.flush() 
image_mass.close() 

cmd = '%s/fits2ndf in=%s out=%s '%(convdir,massFITS,masssdf)
os.system(cmd)


#Convert map of mass (in solar masses) into column density (in H2 cm-2) and Extinction (mag)
M_x = 1.989E33 #g
N = 1
au = 14959787100000 #cm
m_h = 1.67262178E-24 #g
mu = 2.8 #333  #ratio of H2 to He (Kauffmann et al. 2008)

#pixel area
A = ((((a*d)*au)**2.0)*N) #in cm^2 
f = M_x/(mu*m_h*A) #column density per cm^2 

cmd = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,masssdf,f,CDsdf)
os.system(cmd)
av = 0.9E21
cmd = '%s/cdiv in=%s scalar=%f out=%s'%(kapdir,CDsdf,av,Av)
os.system(cmd)
