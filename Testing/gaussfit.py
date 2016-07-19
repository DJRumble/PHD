#Damian Rumble, UoE
#20141215
#gaussfit.py

#this script takes a map of mean value with statisitical noise anc calculates statistics on it. 

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
import matplotlib.pyplot as plt
import copy

#################################################################################
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

def myfunct(p, fjac=None, x=None, y=None, err=None):
    model = func(p,x)
    # Non-negative status value means MPFIT should
    # continue, negative means stop the calculation.
    status = 0
    return [status, (y-model)/err]


########### Bulk Code #############

print 'Select map to test (.sdf format)'
input = str(sys.argv[1])

mean = PARGET(input,'MEAN','stats')

flux_ulim = mean*1.1
flux_llim = mean*0.9

#convert input maps into fits format
prefix = string.split(input,'.sdf')[0]
print 'Converting %s to fits...'%(input)
fits = prefix+'.fits'
#print 'Converting %s to %s'%(input,fits)
if (os.path.exists(fits)):
    os.unlink(fits)
cmd = '%s/ndf2fits in=%s out=%s QUIET'%(convdir,input,fits)
os.system(cmd)

#Model in the form of a gaussian
func = lambda p,x: p[0]*np.exp(-(x/(np.sqrt(2)*p[1]))**2)

#Set parameters for model fitting limits
p0 = [1000,0.08]
bins = 1000

param_base={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
param_info=[]
for i in range(len(p0)):
    param_info.append(copy.deepcopy(param_base))
for i in range(len(p0)): 
    param_info[i]['value']=p0[i]

param_info[0]['limited'] = [1,1]       # Set the lower bound of parameter 0 to 'True' (i.e. 1)
param_info[0]['limits']  = [0.0,100000.0]   # Set the lower limit of parameter 0 to 0.0

param_info[1]['limited'] = [1,1]       # Set the lower bound of parameter 1 to 'True' (i.e. 1)
param_info[1]['limits']  = [0.0,2.0] 

hdulist = pyfits.open(fits)
data = hdulist[0].data


#calculate histogram
loc = np.where((data!=0.0) & (data<flux_ulim) & (data>flux_llim))
hist, bin_edges = np.histogram((data[loc]), bins = bins)

print hist
print bin_edges

#format data/bins for use in MPFIT
fa = {'x':bin_edges[:-1], 'y':hist, 'err':(np.zeros(np.shape(hist))+1)}

print 'fitting gaussian model to data'
m = mpfit.mpfit(myfunct, p0, parinfo=param_info,functkw=fa,quiet=True)
print 'STDV on the noise = '+str(m.params[1])+' Jy/pixel'

logic = 'TRUE'

if logic == 'TRUE':
    print "Plotting..."
    plt.clf()
    plt.scatter(bin_edges[:-1], hist)
    plt.plot(bin_edges[:-1],func(m.params,bin_edges[:-1]),c='r')
    plt.xlabel('arb. unuts')
    plt.ylabel('Number of pixels')
    plt.show()
else:
    print 'no plot'
