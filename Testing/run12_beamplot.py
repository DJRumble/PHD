#DJR UoE
#16/10/2013 - run12.py
#
#This code test the Gausmooth and produces plots of Primary and Secondary beam 

###########  Import modules ######################################

import os
import numpy as np
import string
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
import copy

########### set CSH Commands directories #############

kapdir = '/stardev/bin/kappa'
convdir = '/stardev/bin/convert'

###################################################################
#Start of free parameters
###################################################################



###################################################################
#Start of fixed parameters
##################################################################

#null maps - move
s450 = 'null_s450.sdf'
s850 = 'null_s850.sdf'

#normalisation constants
A_a = 0.6118852617 #450_MB
A_b = 0.3911282854 #450_S
B_a = 0.7825471052 #850_MB
B_b = 0.2177259426 #450_S

#alpha/beta
a4 = 0.94
b4 = 0.06
a8 = 0.98
b8 = 0.02

#fwhm
M4 = 7.9 #fwhm450_MB
M8 = 13.0 #fwhm850_MB
S4 = 25.0 #fwhm450_S
S8 = 48.0 #fwhm850_S

#fwhm conversion
ln = 2.354820045 #2SQRT(2ln2)

#pixel sizes
p450 = 4.0
p850 = 6.0

###################################################################
#Start of bulk script
###################################################################

def exp(M,S,A,B,p,null,v):
    
    out = 'MATHS'+str(v)+'MBS.sdf'
    MB = 2*(((M/p)/ln)**2.0)
    S = 2*(((S/p)/ln)**2.0)
    print MB, S
    exp_MB = "'exp(-(((xa+0.5)**2)+((xb+0.5)**2))/"+str(MB)+")'"
    exp_S = "'exp(-(((xa+0.5)**2)+((xb+0.5)**2))/"+str(S)+")'"
    print exp_MB
    print exp_S

    out_MB = 'temp'+str(v)+'_out_MB.sdf'
    out_S = 'temp'+str(v)+'_out_S.sdf'

    mult_MB = 'temp'+str(v)+'_mult_MB.sdf'
    mult_S = 'temp'+str(v)+'_mult_S.sdf'

    cmd = '%s/maths exp=%s out=%s like=%s'%(kapdir,exp_MB,out_MB,null)
    os.system(cmd)
    cmd = '%s/maths exp=%s out=%s like=%s'%(kapdir,exp_S,out_S,null)
    os.system(cmd)

    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,out_MB,A,mult_MB)
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,out_S,B,mult_S)
    os.system(cmd)

    cmd = '%s/add in1=%s in2=%s out=%s'%(kapdir,mult_MB,mult_S,out)
    os.system(cmd)
    return 


# MBS 450
MBS_450 = exp(M8,S8,a4,b4,p450,s450,450)

# MBS 850
MBS_850 = exp(M4,S4,a8,b8,p850,s850,850)

device = 'xwi'
clear = 'false'



#make LINPLOT 
#MATHS450
incpy = "MATHS450MBS'(-15:15,0:0)'"
outcpy = 'MATHS_lin450MBS.sdf'

cmd = '%s/ndfcopy in=%s out=%s'%(kapdir,incpy,outcpy) 
os.system(cmd)
cmd = '%s/linplot ndf=%s device=%s'%(kapdir,outcpy,device)
os.system(cmd)

#make LINPLOT 
#MATHS850
incpy = "MATHS850MBS'(-10:10,0:0)'"
outcpy = 'MATHS_lin850MBS.sdf'


cmd = '%s/ndfcopy in=%s out=%s'%(kapdir,incpy,outcpy) 
os.system(cmd)
cmd = '%s/linplot ndf=%s device=%s clear=%s'%(kapdir,outcpy,device,clear)
os.system(cmd)


#make lineplot
#Convolve450
outcpy = 'convolve/line_s450convolve.sdf'
cmd = '%s/linplot ndf=%s device=%s clear=%s'%(kapdir,outcpy,device,clear)
os.system(cmd)

#make lineplot
#Convolve850
outcpy = 'convolve/line_s850convolve.sdf'
cmd = '%s/linplot ndf=%s device=%s clear=%s'%(kapdir,outcpy,device,clear)
os.system(cmd)


print 'plots should appear now...'
