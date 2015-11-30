#Damian Rumble, UoE
#24/10/2013
#run13_jeansmass.py

#this is an advanced code, based on 'mass_calc.py' that converts input 850flux and temp. maps into column density and ulitmately Jeans Mass maps for various clumps, as deffined by ClumpFind FW. 

#####################################################
#import maths, ploting and astrophysical packages

import os
import numpy as np
import aplpy
import matplotlib.pyplot as plt
import atpy
from pylab import *

#My Modules
import BlackBody
import Johnstone
import mass
import flux
import kirk
import rumble_mass
import mass_maps
import jeans

########### set CSH Commands directories #############

kapdir = '/stardev/bin/kappa'
kapdir2 = '/star/bin/kappa'
convdir = '/stardev/bin/convert'

#####################################################
#Constants

au = 149597871000
p = 3.08567758E16
M_x = 1.989E30
h = 6.626068E-34
c = 2.99792458E8
k = 1.3806488E-23
m_h2 = 3.34745E-24 #g
G = 6.67384E-11
pi = 3.14159265359
mu = 2.333333  #ratio of H2 to He (5:1)

#####################################################
#Function to deffine Column Density - Johnstone00 (850 only)

def ColumnD(S,T,kappa,N,d):
    #each pixel has an area in SI units and each apature contains N pixels
    a = 5.99987
    A = ((((a*d)*au)**2.0)*N)*10000 #in cm
    t = 17.0/T
    M = 0.39*S*(np.exp(t)-1.0)*((kappa/0.012)**(-1.0))*((d/250.0)**2.0) #in solarmasses per pixel
    m = M*M_x*1000.00 #A = 1 density in g per pixel

    n = m/(mu*m_h2*A) #column density per cm^2 
    N = m/A #density in g per cm^2
    print 'mass per pixel in solar masses:', M
    print 'column density of apature in per cm^2:', n
    print 'density of apature in g per cm^2:', N
    return


#####################################################
#parameters

kappa = 0.012 #g/cm2
d = 250.0 #parsecs

#####################################################
#bulk code - this code makes maps of jeams mass - current of the whole map.

s850 = 'SerpensMWC297_20130918_232648_s850_DJR_cal_DJR_clumps.sdf'
temp = 'SerpensMWC297_20130916_233448_DJR_cal450850temp1_8snr5%5_DJR_clumps.sdf'
Mtemp = np.loadtxt('SMM_temps.txt', dtype='float')

#s850 = 'mwc297s_test.sdf'
#temp = 'mwc297T_test.sdf'

#s850 = 'johnstone_test.sdf'
#temp = 'johnstoneT_test.sdf'

#####################################################
#read data

signal = s850 #flux per pixel
#T = Mtemp[:,0]
T = temp
t = 40.0
M = np.asarray([0.1,0.5,1.0,2.0,5.0,10.0])
N = Mtemp[:,5]
#print M[0]

print 'making maps with flux map\n'+str(signal)+'\n and temperature map \n'+ str(T)

mass_maps.mass_maps(T,signal,kappa,d)

print "we're done here"

#########################
