#Damian Rumble, UoE
#20140211
#Run17

#A script for producing a large number of temperature maps in a normal distrabution as a test.

####################################################
#import maths, ploting and astrophysical packages

import random
import numpy as np
import os
import subprocess
import commands 

########### set CSH Commands directories #############

kapdir = '/star/bin/kappa'

########### Bulk Code #############

### distribution parameters ### 

mean_450 = 4
mean_850 = 2
sigma = 0.1

irange = 5

### distribution parameters ### 

s450 = 'SerpensMWC297_20130425_s450_IR1_JH.sdf'
s850 = 'SerpensMWC297_20121221_s850_IR1_JH.sdf '

s450_out = '/scratch/damian/run17/input/s450_out.sdf'
s850_out = '/scratch/damian/run17/input/s850_out.sdf'


#temp_code = ''

i = 1

T = []
R = []
A = []
B = []


while i  in range(irange):



    a = random.normalvariate(mean_450,sigma)
    cmd = '%s/thresh in=%s out=%s thrlo=0 thrhi=0 newlo=%d newhi=%d > /dev/null '%(kapdir,s450,s450_out,a,a)
    os.system(cmd)
    print 'flux 450: '+str(a)
    b = random.normalvariate(mean_850,sigma)
    cmd = '%s/thresh in=%s out=%s thrlo=0 thrhi=0 newlo=%d newhi=%d %s '%(kapdir,s850,s850_out,b,b,'> /dev/null')
    os.system(cmd)  
    print 'flux 450: '+str(b)
    csh_command = 'source ./ratiomap450850_test.csh input/ output/ s450_out s850_out 99 mwc297 > /dev/null'
    subprocess.check_call(['/bin/csh', '-c', csh_command])

    ratio = 'output/map/s450_out/s450_out450850.sdf'
    temp = 'output/map/s450_out/s450_out450850temp%99.sdf' 

    cmd = '%s/stats ndf=%s %s'%(kapdir,ratio,'> /dev/null')
    os.system(cmd)
    cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
    status, output = commands.getstatusoutput(cmd)
    ratio_i = output
    print 'flux ratio: '+str(ratio_i)
    R.append(float(ratio_i))
    
    cmd = '%s/stats ndf=%s %s'%(kapdir,temp,'> /dev/null')
    os.system(cmd)
    cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
    status, output = commands.getstatusoutput(cmd)
    temp_i = output
    print 'flux temp: '+str(temp_i)
    T.append(float(temp_i))

    A.append(a)
    B.append(b)

    i = i + 1

    print 'process: '+str(i)+' of '+str(irange)

print 'input 450 fluxes: ',A
print 'input 850 fluxes: ',B
print 'output ratios: ',R
print 'output temperatures: ',T
