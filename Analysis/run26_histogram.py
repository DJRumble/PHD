#20150630
#Damian Rumble, UoE

#this script produces histograms of clump properties

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

data =  np.loadtxt('data/version1/SMM-orion.tab',dtype='float',comments="#")
data15 =  np.loadtxt('data/version1/SMM_realtemp_full.tab',dtype='float',comments="#")
dataYSO =  np.loadtxt('data/version1/SMM_realtemp_YSO_full.tab',dtype='float',comments="#")

#MASS
MASS = np.array(data[:,3]).tolist() 
dMASS = np.array(data[:,4]).tolist() 
MASSYSO = np.array(dataYSO[:,3]).tolist() 
dMASSYSO = np.array(dataYSO[:,4]).tolist() 
MASS15 = np.array(data15[:,3]).tolist()
dMASS15 = np.array(data15[:,4]).tolist()
#TEMP
TEMP = np.array(data[:,5]).tolist()
TEMP15 = np.array(data15[:,5]).tolist()
TEMPYSO = np.array(dataYSO[:,5]).tolist()
dTemp = np.array(data[:,6]).tolist()
dTempYSO = np.array(dataYSO[:,6]).tolist()
dTEMP15 = np.array(data15[:,6]).tolist()
#CD
CD = np.array(data[:,8]).tolist()
CD15 = np.array(data15[:,8]).tolist()
CDYSO = np.array(dataYSO[:,8]).tolist()
dCD = np.array(data[:,9]).tolist()
dCDYSO = np.array(dataYSO[:,9]).tolist()
dCD15 = np.array(data15[:,9]).tolist()
#MJ
MJ = np.array(data[:,16]).tolist()
MJ15 = np.array(data15[:,16]).tolist()
MJYSO = np.array(dataYSO[:,16]).tolist()
dMJ = np.array(data[:,17]).tolist()
dMJYSO = np.array(dataYSO[:,17]).tolist()
dMJ15 = np.array(data15[:,17]).tolist()

nullfmt = NullFormatter()         # no labels

# definitions for the axes
xH1,xH2,xS1,xS2,xS3 = .07,.18,.21,.58,.96
yH1,yH2,yS1,yS2,yS3 = .07,.18,.21,.58,.96
wS,hS,wH,hH = 0.37,0.37,0.10,0.10

rect_scatter1 = [xS1, yS1, wS, hS]
rect_scatter2 = [xS1, yS2, wS, hS]
rect_scatter3 = [xS2, yS1, wS, hS]
rect_scatter4 = [xS2, yS2, wS, hS]
rect_hist1= [xS1, yH1, wS, hH]
rect_hist2= [xS2, yH1, wS, hH]
rect_hist3= [xH1, yS1, wH, hS]
rect_hist4= [xH1, yS2, wH, hS]


# start with a rectangular Figure
plt.figure(1, figsize=(8, 8))

axScatter1 = plt.axes(rect_scatter1)
axScatter2 = plt.axes(rect_scatter2)
axScatter3 = plt.axes(rect_scatter3)
axScatter4 = plt.axes(rect_scatter4)
axHist1 = plt.axes(rect_hist1)
axHist2 = plt.axes(rect_hist2)
axHist3 = plt.axes(rect_hist3)
axHist4 = plt.axes(rect_hist4)

#### the scatter plots: ####

#CD vs temp
#axScatter2.errorbar(TEMP15,CD15,0,0,'wo',ecolor='k')
axScatter2.errorbar(TEMP15,dCD15,0,0,'wo',ecolor='k')
axScatter2.errorbar(TEMPYSO,CDYSO,0,0,'ko')
axScatter2.xaxis.set_major_formatter(nullfmt)
axScatter2.set_xlim((6,38))
axScatter2.set_ylim((0,2500))

#CD vs Mass
axScatter4.errorbar(MASS,CD,0,0,'wo',ecolor='k')
#axScatter4.errorbar(MASS15,CD15,0,0,'wo',ecolor='k')
axScatter4.errorbar(MASSYSO,CDYSO,0,0,'ko')
axScatter4.xaxis.set_major_formatter(nullfmt)
axScatter4.yaxis.set_major_formatter(nullfmt)
axScatter4.set_xlim((-0.5,90))
axScatter4.set_ylim((0,2500))

#Stability vs Temp
axScatter1.errorbar(TEMP,MJ,0,0,'wo',ecolor='k')
#axScatter1.errorbar(TEMP15,MJ15,0,0,'wo',ecolor='k')
axScatter1.errorbar(TEMPYSO,MJYSO,0,0,'ko')
axScatter1.set_xlim((6,38))
axScatter1.set_ylim((0,11))
axScatter1.axhline(y=1,color='k',linestyle='--')

#Stability vs Mass
axScatter3.errorbar(MASS,MJ,0,0,'wo',ecolor='k')
#axScatter3.errorbar(MASS15,MJ15,0,0,'wo',ecolor='k')
axScatter3.errorbar(MASSYSO,MJYSO,0,0,'ko')
axScatter3.yaxis.set_major_formatter(nullfmt)
axScatter3.set_xlim((-0.5,90))
axScatter3.set_ylim((0,11))
axScatter3.axhline(y=1,color='k',linestyle='--')

#### the histograms ####

#MASS
bins = 100
axHist2.yaxis.set_major_formatter(nullfmt)
axHist2.set_xlabel('Mass (M$_{\odot}$)')
axHist2b = axHist2.twinx()
axHist2b.hist(MASS, bins=bins,color='r')
axHist2b.hist(MASS15, bins=bins,color='b')
axHist2b.hist(MASSYSO, bins=bins,color='y')
axHist2b.set_xlim((-0.5,90))
#axHist2b.set_yscale('log')

#TEMP
bins = 150
axHist1.hist(TEMP, bins=bins,color='r')
axHist1.hist(TEMP15, bins=bins,color='b')
axHist1.hist(TEMPYSO, bins=bins,color='y')
axHist1.set_xlim((6,38))
axHist1.set_xlabel('Temp. (K)')
#axHist1.set_yscale('log')

#STABILITY
bins = 43
axHist3.hist(MJ, bins=bins, orientation='horizontal',color='r')
axHist3.hist(MJ15, bins=bins, orientation='horizontal',color='b')
axHist3.hist(MJYSO, bins=bins, orientation='horizontal',color='y')
axHist3.set_ylim((0,11))
axHist3.set_ylabel('M/M$_{\mathrm{J}}$')
#axHist3.set_xscale('log')

#CD
bins = 25
axHist4.xaxis.set_major_formatter(nullfmt)
axHist4.set_ylabel('Peak column density (10$^{21}$ H$_{2}$ cm$^{-2}$)')
axHist4b = axHist4.twiny()
axHist4b.hist(CD, bins=bins, orientation='horizontal',color='r')
axHist4b.hist(CD15, bins=bins, orientation='horizontal',color='b')
axHist4b.hist(CDYSO, bins=bins, orientation='horizontal',color='y')
axHist4b.set_ylim((0,2500))
#axHist4b.set_xscale('log')

show()
