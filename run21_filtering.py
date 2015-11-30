#20150803
#Damian Rumble, UoE

#This is replacement script for the that which produces ratio histograms for pixels in the W40-N region

#################################################################################
########### import modules  #############
#Standard modules
import numpy as np
import astropy.io.fits as pyfits
import sys
import os
import string
import mpfit
import matplotlib.pyplot as plt
import copy
from pylab import *
from scipy import stats

#################################################################################
########### set CSH Commands directories #############
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'
cupidr = '/star/bin/cupid'

#################################################################################
########## Bulk script #############


box450 = [30,60,90,120,135,150,180]
box850 = [20,40,60,80,90,100,120]
label = ['1am','2am','3am','4am','4-30am','5am','6am']

in450 = 'input/Aquila_20150316_s450_extmask_s2nosm_cal_Jypix_jh.sdf'
in850 = 'input/Aquila_20150731_noco_extS2nosm_s850_IR2_Jypix_DJR.sdf'

l = len(label)

findback = 'False'

if findback == 'TRUE': 
    for i in range(l):
        BOX450 = box450[i]
        BOX850 = box850[i]
        LABEL = label[i]
        print 'running map ',LABEL
        out450 = 'filtered_noco/Aquila_extS2nosm_s450-'+LABEL+'.sdf'
        out850 = 'filtered_noco/Aquila_noco_extS2nosm_s850-'+LABEL+'.sdf'
        cmd450 = '%s/findback in=%s out=%s box=%i sub=True RMS=!'%(cupidr,in450,out450,BOX450)
        cmd850 = '%s/findback in=%s out=%s box=%i sub=True RMS=!'%(cupidr,in850,out850,BOX850)
        os.system(cmd450)
        os.system(cmd850)

label = ['3am','4am','4-30am','5am','6am']

l = len(label)

ratio = 'FALSE'

if ratio == 'TRUE':
    for i in range(l):
        LABEL = label[i]
        print '============='
        print 'RATIO MAP ',LABEL
        print '============='
        cmd = 'python ratiomap450K850_v21.py Aquila_noco-%s Aquila_extS2nosm_s450-%s Aquila_noco_extS2nosm_s850-%s K'%(LABEL,LABEL,LABEL)
        os.system(cmd)

#remove map edges
list = []
thresh  = 'TRUE'

label = ['3am','4am','5am','6am']
l = len(label)

if thresh == 'TRUE':
    for i in range(l):
        LABEL = label[i]
        infile = '/output/kernel/map/Aquila_noco-%s/Aquila_noco-%s_Sratio.sdf'%(LABEL,LABEL)
        outfile = 'filtered_noco/Aquila_noco-%s_SratioE.sdf'%(LABEL)
        cmd = "%s/mult in1=%s in2=Aquila_s850_edgeless.sdf out=%s"%(kapdir,infile,outfile)
        #os.system(cmd)
  
        list.append(outfile)

print list

for i in range(l):
    input = list[i]
    prefix = string.split(input,'.sdf')[0]
    #print 'Converting %s to fits...'%(input)
    fits = prefix+'.fits'
    #print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(fits)):
	   os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s noprohis QUIET'%(convdir,input,fits)
    os.system(cmd)

    hdulist = pyfits.open(fits)
    data = hdulist[0].data

    loc = np.where((data<14.) & (data>0.) & (data!=0.0))
    Data = data[loc]
    bins=75

    n, bin_edges, patches  = hist(Data, bins, histtype='step',label=str(label[i]))

original = 'input/Aquila_noco_s2nosm_SratioMSK.fits'

hdulist = pyfits.open(original)
data = hdulist[0].data

loc = np.where((data<14.) & (data>0.) & (data!=0.0))
Data = data[loc]
bins=75
    
n, bin_edges, patches  = hist(Data, bins, histtype='step',label='none')

ylabel('Frequency')
legend()
xlim([3,14])
X = [9.5,9.5]
Y = [0,700]
plot(X,Y,'--',color='k')

show()

label = ['S2-4am','S2-CO-4am']

files = ['filtered_noco/Aquila_noco-4am_SratioE.fits','output/map/Aquila_extS2-4as/Aquila_extS2-4as_Sratio_edgeless.fits']

l = len(label)
list = []

for i in range(l):
    fits = files[i]
    hdulist = pyfits.open(fits)
    data = hdulist[0].data
    loc = np.where((data<14.) & (data>0.) & (data!=0.0))
    Data = data[loc]
    bins=75
    n, bin_edges, patches  = hist(Data, bins, histtype='step',label=str(label[i]))
    list.append(n)

D,P = stats.ks_2samp(list[0],list[1])
print 'KS-stats: D = ',D,' P = ',P

ylabel('Frequency')
legend()
xlim([3,14])

show()
