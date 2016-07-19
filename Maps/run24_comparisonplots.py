#Damian Rumble, UoE
#20150106
#TvT.py

#this script takes two temp maps for the same region made using different methods and plots Ta vs Tb to compare the differences (and errors). 

####################################################
#import maths, ploting and astrophysical packages

import random
import numpy as np
import os
import subprocess
import commands 
import astropy.io.fits as pyfits
import string
from pylab import *
import matplotlib.pyplot as plt
from scipy import stats

########### set CSH Commands directories #############

kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

################################################################################
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

########### Bulk Code #############
#Method A = Dual Beam
#Method B = Kernel

#### 17/06/2016 - This script has been updated to calculate the STD of two alignment methods.


inA = 'ratio/Aquila_IR2_s2_DB_Sratio_mask.sdf' 
inA_fits = 'ratio/Aquila_IR2_s2_DB_Sratio_mask.fits' 
inB = 'ratio/Aquila_IR2_s2_K_Sratio_mask.sdf' 
inB_fits = 'ratio/Aquila_IR2_s2_K_Sratio_mask.fits'
inC = 'ratio/Aquila-4-30am_SratioBLmask.sdf'
inC_fits = 'ratio/Aquila-4-30am_SratioBLmask.fits'

in_sincsinc = 'ratio/Aquila_noco-auto_Sratio_sincsincMSK.fits'
in_nearest = 'ratio/Aquila_noco-auto_Sratio_nearestMSK.fits'
in_bilinear = 'ratio/Aquila_noco-auto_Sratio_bilinearMSK.fits'

main_sb = 'ratio_methods/SerpensMain_20141223-auto_Sratio_MB_MSKTH.fits'
main_db = 'ratio_methods/SerpensMain_20141223-auto_SratioMSKTH.fits'
aquila_sb = 'ratio_methods/Aquila_noco-auto_SratioSB_MSKTH.fits'
aquila_db = 'ratio_methods/Aquila_noco-auto_Sratio_DB_MSKTH.fits'
east_sb = 'ratio_methods/SerpensE_20141219-auto_Sratio_SB_MSKTH.fits'
east_db = 'ratio_methods/SerpensE_20141219-auto_Sratio_DB_MSKTH.fits'

Tmain_sb = 'input_temp/compPLOTdir/SerpensMain_20141223-autotemperature_SB_MSK.fits'
Tmain_db = 'input_temp/compPLOTdir/SerpensMain_20141223-autotemperature_DB_MSK.fits'
Taquila_sb = 'input_temp/compPLOTdir/Aquila_noco-autotemperature_SB_MSK.fits'
Taquila_db = 'input_temp/compPLOTdir/Aquila_noco-autotemperature_DB_MSK.fits'
Teast_sb = 'input_temp/compPLOTdir/SerpensE_20141219-autotemperature_SB_MSK.fits'
Teast_db = 'input_temp/compPLOTdir/SerpensE_20141219-autotemperature_DB_MSK.fits'





#DO NOT INCLUDE NORTH
#north_sb = 'ratio_methods/SerpensN_20141219-auto_SratioSB_MSKTH.fits'
#north_db = 'ratio_methods/SerpensN_20141219-auto_SratioSB_MSKTH.fits'

erA = 'ratio/Aquila_IR2_s2_DBSratio_error_mask.sdf'
erA_fits = 'ratio/Aquila_IR2_s2_DBSratio_error_mask.fits'
erB = 'ratio/Aquila_IR2_s2_KSratio_error_mask.sdf'
erB_fits = 'ratio/Aquila_IR2_s2_KSratio_error_mask.fits'
erC = 'ratio/Aquila-4-30amSratioBL_errormask.sdf'
erC_fits = 'ratio/Aquila-4-30amSratioBL_errormask.fits'

'''
inA = 'ratio/mwc297_DBtest_Sratio_SMM1.sdf' 
inA_fits = 'ratio/mwc297_DBtest_Sratio_SMM1.fits'
inB = 'ratio/mwc297_Ktest_Sratio_SMM1.sdf' 
inB_fits = 'ratio/mwc297_Ktest_Sratio_SMM1.fits'

erA = 'ratio/mwc297_DBtestSratio_error_SMM1.sdf'
erA_fits = 'ratio/mwc297_DBtestSratio_error_SMM1.fits'
erB = 'ratio/mwc297_KtestSratio_error_SMM1.sdf'
erB_fits = 'ratio/mwc297_KtestSratio_error_SMM1.fits'

inA = 'ratio/mwc297_DBtest_Sratio.sdf' 
inA_fits = 'ratio/mwc297_DBtest_Sratio.fits'
inB = 'ratio/mwc297_Ktest_Sratio.sdf' 
inB_fits = 'ratio/mwc297_Ktest_Sratio.fits'

erA = 'ratio/mwc297_DBtestSratio_error.sdf'
erA_fits = 'ratio/mwc297_DBtestSratio_error.fits'
erB = 'ratio/mwc297_KtestSratio_error.sdf'
erB_fits = 'ratio/mwc297_KtestSratio_error.fits'

inA = 'temp/mwc297_DBtesttemperature.sdf' 
inA_fits = 'temp/mwc297_DBtesttemperature.fits'
inB = 'temp/mwc297_Ktesttemperature.sdf' 
inB_fits = 'temp/mwc297_Ktesttemperature.fits'

erA = 'temp/mwc297_DBtesttemp_error.sdf'
erA_fits = 'temp/mwc297_DBtesttemp_error.fits'
erB = 'temp/mwc297_Ktesttemp_error.sdf'
erB_fits = 'temp/mwc297_Ktesttemp_error.fits'
'''

if (os.path.exists(inC_fits)):
    os.unlink(inC_fits)
if (os.path.exists(inB_fits)):
    os.unlink(inB_fits)
if (os.path.exists(erC_fits)):
    os.unlink(erC_fits)
if (os.path.exists(erB_fits)):
    os.unlink(erB_fits)

#NDF2FITS(in_bilinear)
#NDF2FITS(in_nearest)
#NDF2FITS(erC)
#NDF2FITS(erB)

SBlist = [Tmain_sb,Taquila_sb,Teast_sb]
DBlist = [Tmain_db,Taquila_db,Teast_db]

methodlist = [in_sincsinc,in_bilinear]

k = 0

TA = []
TB = []
EA = []
EB = []

count_lo = 0
count_hi = 0


for k in range(len(methodlist)):

    M = methodlist[k]
    M2 = methodlist[k+1]
    print 'MAP:',M
    print 'MAP:',M2
    imageA = pyfits.open(str(M), mode = 'update')
    imageB = pyfits.open(str(M2), mode = 'update')
#errorA = pyfits.open(str(erC_fits), mode = 'update')
#errorB = pyfits.open(str(erB_fits), mode = 'update')

    dimA = PARGET(M,'dims','ndftrace')
    dimB = PARGET(M2,'dims','ndftrace')

    rowA = int(dimA[1])    #rows of the table from NDFTRACE
    columnA = int(dimA[0])    #columns of the table from NDFTRACE
    rowB = int(dimB[1])    #rows of the table from NDFTRACE
    columnB = int(dimB[0])    #columns of the table from NDFTRACE

    if (rowA == rowB) & (columnA == columnB):
        print 'Maps are the same region: GOOD'
        print 'row = ' + str(rowA) +':'+ str(rowB)
        print 'column =' + str(columnA) +':'+ str(columnB)
    else:
        print 'Maps are not the same region: BAD'
        print 'row = ' + str(rowA) +':'+ str(rowB)
        print 'column =' + str(columnA) +':'+ str(columnB)
        sys.exit()

    dataA = imageA[0].data
    dataB = imageB[0].data
#errordataA = errorA[0].data
#errordataB = errorB[0].data

    print 'imageA =' + str(imageA)
    print 'imageB =' + str(imageB)

    for i in range(0, rowA):
    #print i
        for j in range(0, columnA):
        #print j
            if (dataA[i][j] > 0.) & (dataA[i][j] < 120.5) & (dataB[i][j] > 0.) & (dataB[i][j] < 120.5):
                TA.append(dataA[i][j])
                TB.append(dataB[i][j])
            #EA.append(errordataA[i][j]/2.)
            #EB.append(errordataB[i][j]/2.) 

                Tr = dataA[i][j]/dataB[i][j]

                if Tr > 1.:
                    count_lo = count_lo + 1
                if Tr <= 1.:
                    count_hi = count_hi + 1

    break
    print len(TA)

l = len(TA)

ratio = []


for i in range(len(TA)):
    ratio.append(TA[i]/TB[i])

print 'RATIO - STD:',np.std(ratio)


print 'single > dual: ',count_lo
print 'single < dual: ',count_hi
print 'TOTAL: ',count_lo+count_hi

for i in range(l):
    #print '=============================='
    #print 'Wavelength start ', data[i][3]
    xi = TA 
    yi = TB
    
    A = array([xi, ones(9)])

    slope, intercept, r_value, p_value, slope_std_err = stats.linregress(xi,yi)

print 'intercept', intercept
print 'slope', slope
print 'standard deviation', slope_std_err
print 'R value',r_value
print 'P value',p_value

X = np.arange(0,20,1)

line = (slope*X)+intercept

plt.scatter(TA,TB,s=1)
#plot(X,line,'-')
#plt.errorbar(TA,TB,xerr=EA,yerr=EB,fmt='+')
plt.plot([0,100],[0,100],'r-',lw=2)
#plt.xlabel('Ratio by Dual Beam Method')
#plt.ylabel('Ratio by Kernel Method')
plt.xlabel('Temperature by Single-Beam method (K)')
plt.ylabel('Temperature by Dual-Beam method (K)')
#ylim([1.5,12.5])
#xlim([1.5,12.5])
ylim([5,75])
xlim([5,75])

if (os.path.exists(inA_fits)):
    os.unlink(inA_fits)
if (os.path.exists(inB_fits)):
    os.unlink(inB_fits)
#if (os.path.exists(erA_fits)):
#    os.unlink(erA_fits)
#if (os.path.exists(erB_fits)):
#    os.unlink(erB_fits)

#plt.show()


