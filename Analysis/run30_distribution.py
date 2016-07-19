#20160218version1/
#Damian Rumble, UoE

#this script produces plots of clump properties

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

#################################################################################
########## Bulk script ############


#No 15K clumps
#datanoOBnoYSO =  np.loadtxt('data/SMM_master_real_noOBheated_nYSOcut60_noYSO_H.tab',dtype='string',comments="#") 
#datanoOBYSO =  np.loadtxt('data/SMM_master_real_noOBheated_nYSOcut60_YSO_H.tab',dtype='string',comments="#") 
#dataOBnoYSO =  np.loadtxt('data/Version2/SMM_master-real-OB-noYSO_FLT.tab',dtype='float',comments="#")
#dataOBYSO =  np.loadtxt('data/Version2/SMM_master-real-OB-YSO_FLT.tab',dtype='float',comments="#")

#All clumps
ALLdataOBnoYSO =  np.loadtxt('data/Version2/SMM_masterOB_noYSO_FLT.tab',dtype='float',comments="#")
ALLdataOBYSO =  np.loadtxt('data/Version2/SMM_masterOB_YSO_FLT.tab',dtype='float',comments="#")

noOB = np.loadtxt('data/Version2/SMM_master-real-noOB_FLT.tab',dtype='float',comments="#")
OB = np.loadtxt('data/Version2/SMM_master-real-OB_FLT.tab',dtype='float',comments="#")

#### PROTOSTAR HISTGRAM #####
datanoOBnoYSO =  np.loadtxt('data/SMM_master_real_noOB_noYSO.tab',dtype='string',comments="#") 
datanoOBYSO =  np.loadtxt('data/SMM_master_real_noOB_YSO.tab',dtype='string',comments="#") 
datanoOBnoYSOH =  np.loadtxt('data/SMM_master_real_noOB_noYSO_H.tab',dtype='string',comments="#") 
datanoOBYSOH =  np.loadtxt('data/SMM_master_real_noOB_YSO_H.tab',dtype='string',comments="#") 


'''
#VERSION1
#No 15K clumps
datanoOBnoYSO =  np.loadtxt('data/version1/version1/SMM_masterNO-OB-noYSO_flt.tab',dtype='float',comments="#") 
datanoOBYSO =  np.loadtxt('data/version1/version1/SMM_masterNO-OB-YSO_flt.tab',dtype='float',comments="#") 
dataOBnoYSO =  np.loadtxt('data/version1/version1/SMM_masterOB-noYSO_flt.tab',dtype='float',comments="#")
dataOBYSO =  np.loadtxt('data/version1/version1/SMM_masterOB-YSO_flt.tab',dtype='float',comments="#")

#All clumps
ALLdatanoOBnoYSO =  np.loadtxt('data/version1/version1/SMM_masterNO-OB-noYSO_full_flt.tab',dtype='float',comments="#") 
ALLdatanoOBYSO =  np.loadtxt('data/version1/version1/SMM_masterNO-OB-YSO_full_flt.tab',dtype='float',comments="#") 
ALLdataOBnoYSO =  np.loadtxt('data/version1/version1/SMM_masterOB-noYSO_full_flt.tab',dtype='float',comments="#")
ALLdataOBYSO =  np.loadtxt('data/version1/version1/SMM_masterOB-YSO_full_flt.tab',dtype='float',comments="#")

noOB = np.loadtxt('data/version1/version1/SMM_masterNO-OB-real_flt.tab',dtype='float',comments="#")
OB = np.loadtxt('data/version1/version1/SMM_masterOB-real_flt.tab',dtype='float',comments="#")
'''
#### KEY ####

#1 = noOBnoYSO
#2 = noOBYSO
#3 = OBnoYSO
#4 = OBYSO

#Radius
'''
R1 = np.array(ALLdatanoOBnoYSO[:,10]).tolist() 
R2 = np.array(ALLdatanoOBYSO[:,10]).tolist() 
R3 = np.array(ALLdataOBnoYSO[:,10]).tolist() 
R4 = np.array(ALLdataOBYSO[:,10]).tolist() 


#MASS
M1 = np.array(datanoOBnoYSO[:,3]).tolist() 
dM1 = np.array(datanoOBnoYSO[:,4]).tolist() 
M2 = np.array(datanoOBYSO[:,3]).tolist() 
dM2 = np.array(datanoOBYSO[:,4]).tolist() 
M3 = np.array(dataOBnoYSO[:,3]).tolist() 
dM3 = np.array(dataOBnoYSO[:,4]).tolist() 
M4 = np.array(dataOBYSO[:,3]).tolist() 
dM4 = np.array(dataOBYSO[:,4]).tolist() 

allM1 = np.array(ALLdatanoOBnoYSO[:,3]).tolist() 
alldM1 = np.array(ALLdatanoOBnoYSO[:,4]).tolist() 
allM2 = np.array(ALLdatanoOBYSO[:,3]).tolist() 
alldM2 = np.array(ALLdatanoOBYSO[:,4]).tolist() 
allM3 = np.array(ALLdataOBnoYSO[:,3]).tolist() 
alldM3 = np.array(ALLdataOBnoYSO[:,4]).tolist() 
allM4 = np.array(ALLdataOBYSO[:,3]).tolist() 
alldM4 = np.array(ALLdataOBYSO[:,4]).tolist()
''' 

#TEMP
T1 =  [float(i) for i in datanoOBnoYSO[:,6]]
#dT1 = np.array(datanoOBnoYSO[:,6]).tolist() 
T2 =  [float(i) for i in datanoOBYSO[:,6]]
#dT2 = np.array(datanoOBYSO[:,6]).tolist() 

T1_H =  [float(i) for i in datanoOBnoYSOH[:,6]]
#dT1 = np.array(datanoOBnoYSO[:,6]).tolist() 
T2_H =  [float(i) for i in datanoOBYSOH[:,6]]

#T3 = np.array(dataOBnoYSO[:,5]).tolist() 
#dT3 = np.array(dataOBnoYSO[:,6]).tolist() 
#T4 = np.array(dataOBYSO[:,5]).tolist() 
#dT4 = np.array(dataOBYSO[:,6]).tolist()

#noOBT = np.array(noOB[:,5]).tolist() 
#OBT = np.array(OB[:,5]).tolist() 
'''
allT1 = np.array(ALLdatanoOBnoYSO[:,5]).tolist() 
alldT1 = np.array(ALLdatanoOBnoYSO[:,6]).tolist() 
allT2 = np.array(ALLdatanoOBYSO[:,5]).tolist() 
alldT2 = np.array(ALLdatanoOBYSO[:,6]).tolist() 
allT3 = np.array(ALLdataOBnoYSO[:,5]).tolist() 
alldT3 = np.array(ALLdataOBnoYSO[:,6]).tolist() 
allT4 = np.array(ALLdataOBYSO[:,5]).tolist() 
alldT4 = np.array(ALLdataOBYSO[:,6]).tolist()

#CD
CD1 = np.array(datanoOBnoYSO[:,8]).tolist() 
dCD1 = np.array(datanoOBnoYSO[:,9]).tolist() 
CD2 = np.array(datanoOBYSO[:,8]).tolist() 
dCD2 = np.array(datanoOBYSO[:,9]).tolist() 
CD3 = np.array(dataOBnoYSO[:,8]).tolist() 
dCD3 = np.array(dataOBnoYSO[:,9]).tolist() 
CD4 = np.array(dataOBYSO[:,8]).tolist() 
dCD4 = np.array(dataOBYSO[:,9]).tolist()

#MJ
MJ1 = np.array(datanoOBnoYSO[:,16]).tolist() 
dMJ1 = np.array(datanoOBnoYSO[:,17]).tolist() 
MJ2 = np.array(datanoOBYSO[:,16]).tolist() 
dMJ2 = np.array(datanoOBYSO[:,17]).tolist() 
MJ3 = np.array(dataOBnoYSO[:,16]).tolist() 
dMJ3 = np.array(dataOBnoYSO[:,17]).tolist() 
MJ4 = np.array(dataOBYSO[:,16]).tolist() 
dMJ4 = np.array(dataOBYSO[:,17]).tolist()
'''
######## BULK SCRIPT ########


#print len(allT1),len(allT2),len(T1),len(T2)

#NoOB
#n1,b1,p1 = hist(allT1,(15),histtype='step',color='r', normed=True)
#n2,b2,p2 = hist(allT2,(10),histtype='stepfilled',color='r', normed=True)
n3,b3,p3 = hist(T1_H,(30),histtype='step',color='g', normed=True,alpha=1)
n4,b4,p4 = hist(T2_H,(30),histtype='step',color='y', normed=True,alpha=1)
#hist(product3,(20),histtype='step')
#hist(product4,(20),histtype='step')

#ALL
#n1,b1,p1 = hist(noOBT,(15),histtype='step',color='b',normed=True)
#n2,b2,p2 = hist(OBT,(10),histtype='step',color='r',normed=True)

xlabel('Temp. (K)')
ylabel('Normalsied frequency')

D,P = stats.ks_2samp(n3,n4)
print 'KS-stats: D = ',D,' P = ',P

label1 = 'P value = '+str(round(P,2))
label2 = 'D value = '+str(round(D,2))

#ylim([0,0.15])
xlim([5,25])

#annotate(label1,xy=(20,0.13))
#annotate(label2,xy=(20,0.14))

show()


