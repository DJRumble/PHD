#!/usr/bin/python
# Carl Salji's script
# JH edited filename conversion to fits
import os
import numpy as np
import string
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
import mpfit
import copy

convdir = '/stardev/bin/convert'
work = '/data/damian/maps/'
lists = '/data/damian/maps/'
#mosdir = '/data/hatchell/SCUBA2/GBS/mosaics/extmask_serpens/'
mosdir = '/data/damian/maps/'

files = np.loadtxt(lists+'list.txt',dtype='string',comments="#")
#['']
print files

plot_stuff = False
flux_lim = [-0.02,0.02]
snr_cut = 2.0
plt.ion()

func = lambda p,x: p[0]*np.exp(-(x/(np.sqrt(2)*p[1]))**2)

def myfunct(p, fjac=None, x=None, y=None, err=None):
	model = func(p,x)
	# Non-negative status value means MPFIT should
	# continue, negative means stop the calculation.
	status = 0
	return [status, (y-model)/err]

p0 = [1000,0.04]

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

for myfile in files:
    print myfile
 
    prefix = string.split(myfile,'.sdf')[0]
    #print 'Converting %s to fits...'%(myfile)
    sdf = mosdir+myfile
    fits = work+prefix+'.fits'
    #print 'Converting %s to %s'%(sdf,fits)
    if (os.path.exists(fits)):
	    os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s'%(convdir,sdf,fits)
    os.system(cmd)
    hdulist = pyfits.open(fits)
    data = hdulist[0].data
    variance = hdulist[1].data
    if len(np.shape(data))>2:
        data = np.reshape(data,[np.shape(data)[1],np.shape(data)[2]])
        variance = np.reshape(variance,[np.shape(variance)[1],np.shape(variance)[2]])
    snr = np.abs(data/np.sqrt(variance))
    mask = np.where(snr<snr_cut,1,np.nan)
    loc = np.where((snr<snr_cut) & (data<flux_lim[1]) & (data>flux_lim[0]) & (data!=0.0))
    bins=10000
    hist, bin_edges = np.histogram((data[loc]), bins = bins)
    fa = {'x':bin_edges[:-1], 'y':hist, 'err':(np.zeros(np.shape(hist))+1)}
    m = mpfit.mpfit(myfunct, p0, parinfo=param_info,functkw=fa,quiet=True)
    print 'Noise estimate for %s %s +/- %s mJy/pixel'%(sdf,m.params[1]*1000.,m.perror[1]*1000.)
    os.unlink(fits)

    if plot_stuff==True:
        print "Plotting..."
	plt.clf()
        plt.figure(1)
        plt.imshow(data)
        plt.colorbar()
        #plt.figure(2)
        #plt.imshow(variance)
        #plt.colorbar()
        #plt.figure(3)
        #plt.imshow(data*mask)
        #plt.colorbar()
        #plt.figure(4)
        plt.scatter(bin_edges[:-1], hist)
        plt.plot(bin_edges[:-1],func(m.params,bin_edges[:-1]),c='r')
        plt.xlabel('Flux Density')
        plt.ylabel('Number of pixels')
	plt.show()




