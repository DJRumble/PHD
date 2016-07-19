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
import ndf2fits

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

def zoominratio(i,j,min, max, count,ratio,lookup):
    mid = (min+max)/2  #sets middle
    count += 1         #increments flag
    if(count < 14):
        #if value in top half of search, search again with top half of section of array
        if(ratio[i][j] > lookup[mid][1]):
            zoominratio(i,j,int(mid)-1, max, count,ratio,lookup)
        #if in bottom half, call bottom half of array again
        elif(ratio[i][j] < lookup[mid][1]): 
            zoominratio(i,j,min, int(mid)+1, count,ratio,lookup)
    #once zoomed into small enough section of array, find best fit
    else: 
        result = 1                                                             
        difference = 100
        for k in range(min, max):
            #check if closest match or not
            if (np.abs(lookup[k][1] - ratio[i][j]) < difference):
                result = k #set array index if smallest so far
                difference = np.abs(lookup[k][1] - ratio[i][j]) #update smallest difference 
        ratio[i][j] = lookup[result][0] 
    #output correct temperature once best fit found
    return

def temperature(sdf,fits,var,beta,LOGIC):
    #NOTE see earlier version for parts of script that can calculate variance arrays. They ahve been removed here. 

    #Calculate dim numbers from ndftrace 
    dim = PARGET(sdf,'dims','ndftrace')
    #number of dimensions, this item contains to numbers, typically 571 columns 536 row
##### DEFFINE Parrameters ##########
    image = pyfits.open(str(fits), mode = 'update') #opens FITS file created from NDF
    ratio = image[0].data  #Ratio flux data in a table
    if (LOGIC == 'TRUE'):
        variance =  image[1].data  #opens variance array directly from the ndf array
    elif (LOGIC == 'FALSE'):
        imvariance = pyfits.open(str(var), mode = 'update') #opens Variance FITS file 
        variance =  imvariance[0].data #variance data in a table
    row = int(dim[1])    #rows of the table from NDFTRACE
    column = int(dim[0]) #columns of the table from NDFTRACE   
    lookup = np.zeros((99501, 3), float)         
    upperratio = np.zeros((row, column), float)
    lowerratio = np.zeros((row, column), float)
    T = 5.0 #minimum temperature 
##### Main Program ################
    print 'Beta = ' + str(beta)
    print 'row = ' + str(row)
    print 'column =' + str(column)
    print 'image =' + str(image)
#Look up table is created for range of temp. from 5.0K going up in increments of 0.01K.
    for i in range(0, 99501):
        lookup[i][0] = float(T)
        T += 0.01
        lookup[i][1] = (17**(3+beta)) / (9**(3+beta)) * (np.exp(16.93/lookup[i][0])-1) / (np.exp(31.97/lookup[i][0])-1) #Temperature Equation#
#compare table of ratios with real data and use this to build new temp maps - based on varience
    for i in range(0, row):
        for j in range(0, column):
        #if variance flag set, find upper and lower limits based on ratio variance
            if (ratio[i][j] > 0 and LOGIC == "TRUE"):
                upperratio[i][j] = ratio[i][j] + variance[i][j]
                lowerratio[i][j] = ratio[i][j] - variance[i][j]
    for i in range(0, row):
        for j in range(0, column):
            if ratio[i][j] > 0:
                zoominratio(i,j,0, 99500, 0,ratio,lookup) #call functions that find best temperature fit for ratios

#output new values (temperature from ratio) back to FITS file and save
    image.flush() 
    image.close() 

###################################################################
#Start of bulk script
###################################################################

#input maps
s21 = 'maps/s21/aquila_VLAss_21cm_Jypix.sdf'

s450FF = 'maps/Aquila_extS2nosm_s450-4am-fullfreefree.sdf'
s850FF = 'maps/Aquila_noco_extS2nosm_s850-4am-fullfreefree.sdf'
s450 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_extS2nosm_s450-4am.sdf'
s850 = '/data/damian/maps/aquila/serpens_south/IR2/Aquila_noco_extS2nosm_s850-4am.sdf'

#temp maps
a450 = 'maps/temp/align_S21at450.sdf'
a850 = 'maps/temp/align_S21at850.sdf'
aref450 = 'maps/temp/align450.sdf'
aref850 = 'maps/temp/align850.sdf'
aref450FF = 'maps/temp/align450ff.sdf'
aref850FF = 'maps/temp/align850ff.sdf'
f21 = 'maps/temp/aquila_VLAss-5am_21cm.sdf'
sub450 = 'maps/s21/Aquila_450freefree21.sdf'
sub850 = 'maps/s21/Aquila_850freefree21.sdf'

#final maps
###SCUBA-2 @ VLA
#original
con450 = 'maps/s21/noFF/Aquila-450ats21specs.sdf'
con850 = 'maps/s21/noFF/Aquila-850ats21specs.sdf'
con850_fits = 'maps/s21/noFF/Aquila-850ats21specs.fits'
mask450 = 'maps/s21/noFF/Aquila-450ats21_masked.sdf'
mask850 = 'maps/s21/noFF/Aquila-850ats21_masked.sdf'
#smallscale free-free
con450ff = 'maps/s21/FF/Aquila-ff-450resoVLA.sdf'
con850ff = 'maps/s21/FF/Aquila-ff-850resoVLA.sdf'
con850ff_fits = 'maps/s21/FF/Aquila-ff-850resoVLA.fits'
###VLA @ SCUBA-2
sc450 = 'maps/s21/FF/Aquila_S21scaled450.sdf'
sc850 = 'maps/s21/FF/Aquila_S21scaled850.sdf'
###free-free subtracted maps + masks
mask450FF = 'maps/s21/FF/Aquila_450freefree21_masked.sdf'
mask850FF = 'maps/s21/FF/Aquila_850freefree21_masked.sdf'
#ratio 
ratio = 'maps/s21/noFF/Aquila_ratioVLA.sdf'
ratio_fits = 'maps/s21/noFF/Aquila_ratioVLA.fits'
temp_fits = 'maps/s21/noFF/Aquila_tempVLA.fits'
temp = 'maps/s21/noFF/Aquila_tempVLA.sdf'
ratioFF = 'maps/s21/FF/Aquila_FFratioVLA.sdf'
ratioFF_fits = 'maps/s21/FF/Aquila_FFratioVLA.fits'
tempFF_fits = 'maps/s21/FF/Aquila_FFtempVLA.fits'
tempFF = 'maps/s21/FF/Aquila_FFtempVLA.sdf'

########################################
#remove structure greater than 5' from the 21cm maps 
p21 = PARGET(s21,'fpixscale','ndftrace')
box = (5*60)/p21[0]
findback1 = '%s/findback in=%s out=%s box=%s sub=TRUE RMS=!'%(cupdir,s21,f21,int(box))
os.system(findback1)

#Scale up the 21cm to SCUBA-2 wavelenghts
scalar450 = (21.E-2/450.E-6)**(-0.1)
scalar850 = (21.E-2/850.E-6)**(-0.1)

cmult1 = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,f21,scalar450,sc450)
cmult2 = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,f21,scalar850,sc850)
os.system(cmult1)
os.system(cmult2)

#####################################
#Align SCUBA-2 maps with the 21cm maps
#orignal
align450 = '%s/wcsalign in=%s ref=%s out=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,s450,s21,aref450)
align850 = '%s/wcsalign in=%s ref=%s out=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,s850,s21,aref850)
os.system(align450)
os.system(align850)
#small scale free-free
align450 = '%s/wcsalign in=%s ref=%s out=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,s450FF,s21,aref450FF)
align850 = '%s/wcsalign in=%s ref=%s out=%s method=nearest conserve=TRUE accept QUIET'%(kapdir,s850FF,s21,aref850FF)
os.system(align450)
os.system(align850)

#####################################
#convolve the SCUBA-2 maps to VLA resolution

FWHM450to21 = 43.9
FWHM850to21 = 42.6
smooth450to21 = FWHM450to21/p21[0]
smooth850to21 = FWHM850to21/p21[0]

print smooth450to21
print smooth850to21
#original
gsmooth1 = '%s/gausmooth in=%s fwhm=%s out=%s'%(kapdir, aref450, smooth450to21, con450)
gsmooth2 = '%s/gausmooth in=%s fwhm=%s out=%s'%(kapdir, aref850, smooth850to21, con850)
os.system(gsmooth1)
os.system(gsmooth2)
#small scale free-free
gsmooth1 = '%s/gausmooth in=%s fwhm=%s out=%s'%(kapdir, aref450FF, smooth450to21, con450ff)
gsmooth2 = '%s/gausmooth in=%s fwhm=%s out=%s'%(kapdir, aref850FF, smooth850to21, con850ff)
os.system(gsmooth1)
os.system(gsmooth2)

#sub freefree contamination
sub1 = '%s/sub in1=%s in2=%s out=%s'%(kapdir,con450ff,sc450,sub450)
sub2 = '%s/sub in1=%s in2=%s out=%s'%(kapdir,con850ff,sc850,sub850)
os.system(sub1)
os.system(sub2)

#masking 
mask = 'maps/s21/ratiomask.sdf'

#original
mask1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,con450,mask,mask450)
mask2 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,con850,mask,mask850)
os.system(mask1)
os.system(mask2)
#all free-free
mask1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,sub450,mask,mask450FF)
mask2 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,sub850,mask,mask850FF)
os.system(mask1)
os.system(mask2)

############################
#make ratio and temp maps 
#original
div1 = '%s/div in1=%s in2=%s out=%s'%(kapdir,mask450,mask850,ratio)
os.system(div1)
#freefree
div2 = '%s/div in1=%s in2=%s out=%s'%(kapdir,mask450FF,mask850FF,ratioFF)
os.system(div2)

ndf2fits.ndf2fits(ratio)
ndf2fits.ndf2fits(ratioFF)

temperature(ratio,ratio_fits,ratio_fits,1.8,'FALSE')
temperature(ratioFF,ratioFF_fits,ratioFF_fits,1.8,'FALSE')

mv1 = 'mv %s %s'%(ratio_fits,temp_fits)
mv2 = 'mv %s %s'%(ratioFF_fits,tempFF_fits)
os.system(mv1)
os.system(mv2)

ndf2fits.fits2ndf(temp_fits)
ndf2fits.fits2ndf(tempFF_fits)

#compare temp maps
submap = 'maps/s21/submap.sdf'
sub = '%s/sub in1=%s in2=%s out=%s'%(kapdir,tempFF,temp,submap)
os.system(sub)

'''
################################
#### calculate columndensity ###  #see run26
################################

#deffine constants
kappa = 0.012 #cm^2 g^-1
d = 500 #pc

ndf2fits.ndf2fits(con850)
ndf2fits.ndf2fits(con850ff)

CDFITS = 'maps/s21/noFF/CD.fits'
CDFITSFF = 'maps/s21/FF/CDFF.fits'

#open maps
#orignal
image_s850 = pyfits.open(str(con850_fits)) #opens FITS file created from NDF
image_temp = pyfits.open(str(temp_fits))
image_CD = pyfits.open(str(CDFITS), mode = 'update')
##freefree
image_s850FF = pyfits.open(str(con850ff_fits)) #opens FITS file created from NDF
image_tempFF = pyfits.open(str(tempFF_fits))
image_CDFF = pyfits.open(str(CDFITSFF), mode = 'update')

s850 = image_s850[0].data  #flux data in a table
temp = image_temp[0].data  #temp data in a table
CD = image_CD[0].data  #CD data in a table
s850FF = image_s850FF[0].data  #flux data in a table
tempFF = image_tempFF[0].data  #temp data in a table
CDFF = image_CDFF[0].data  #CD data in a table


print s850[1][0]

#Calculate dim numbers from ndftrace 
dim = PARGET(con850_fits,'dims','ndftrace')
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
        Sff = float(s850FF[i][j])
        Tff = float(tempFF[i][j])
        #CD[i][j] = S/T
        CD[i][j]  = ColumnD(S,T,kappa,15,d)
        CDFF[i][j]  = ColumnD(Sff,Tff,kappa,15,d)
        #print CD[i][j]
  
#output new values back to FITS file and save
image_CD.flush() 
image_CD.close() 
image_CDFF.flush() 
image_CDFF.close() 

ndf2fits.fits2ndf(CD)
ndf2fits.fits2ndf(CDFF)

#compare  CD maps
submap = 'maps/s21/submapCD.sdf'
sub = '%s/sub in1=%s in2=%s out=%s'%(kapdir,CDFF,CD,submap)
os.system(sub)
'''
