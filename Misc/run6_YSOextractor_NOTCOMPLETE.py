#Damian Rumble, UoE
#05/03/2013
#run6.py 

#This script will *hopefully* be able to extract extranous YSO from collated catlogs and store them in a new file for investigation.

#!!!!! NOT COMPLETED !!!!!!!

#!!!!! ARCHIVED !!!!!!!!


#####################################################
#import maths, ploting and astrophysical packages

import numpy 
import aplpy
import matplotlib.pyplot as plt
import atpy
import c2dtable

####################################################
#Load various masks. 

#serpens_south/aquila_s2_mask_cont.fits
#MWC297/serpensmwc297_s2+h2+ard_mask_cont.fits
serpensmain/serpensmain_s2+h_mask_cont.fits
#serpens_E1/serpensE1_s2+h+ard_mask_cont.fits

####################################################
#load GBS catalogue 

tbl = atpy.Table()
c2dtable.c2d_read('catalog-Serp-YSOc+head.tbl',tbl,3)
ra=tbl.data['RA'] + 200.0
dec=tbl.data['DEC']
a = tbl.data['alpha']. 

####################################################
#Load various data catalogs

1 = numpy.loadtxt('20130211c10_c_0_I_south_deg.txt')
2 = numpy.loadtxt('20130207_G09_C0_alpha_mwc297.txt')
3 = numpy.loadtxt('20130207_G09_CI_alpha_mwc297.txt')
4 = numpy.loadtxt('20130207_G09_CI_mwc297.txt')
5 = numpy.loadtxt('20130206Gutermuth09_yso_I_0_serpens_deg.txt')
6 = numpy.loadtxt('20130207_G09_C0_alpha_serpens.txt')
7 = numpy.loadtxt('20130207_G09_CI_alpha_serpens.txt')
8 = numpy.loadtxt('20130131Evans09_yso_c1_serpens_deg.txt') #note these are columns 1,2
9 = numpy.loadtxt('20130131Evans09_yso_FS_serpens_deg.txt') #note these are columns 1,2
10 = numpy.loadtxt('w07_yso_0_1_wcs_deg_only.txt')
11 = numpy.loadtxt('w07_yso_FS_wcs_deg_only.txt')

####################################################
#Deffine data range - there should be away to read this off the file but for now I will just define the range by eye from GAIA

MWC297 = (ra>227.4)&(ra<276.75)&(dec>-4.07)&(dec<-3.5)



####################################################
#Create array

YSO = []
YSO.append(i[:,0],i[:,1])

####################################################
#iterate through data

