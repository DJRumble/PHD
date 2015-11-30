"""A module for calculating noise given and input map. This method masks all data above 2sigma and then fits a gaussian to the resulting histogram"""

#################################################################################
########### import modules  #############
#Standard modules
import numpy as np
#import astropy.io.fits as pyfits
import pyfits
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

def noise_by_data(input,logic): #this function calculates noise

    flux_lim = [-0.02,0.02]
    snr_cut = 2.0

    #convert input maps into fits format
    prefix = string.split(input,'.sdf')[0]
    print 'Converting %s to fits...'%(input)
    fits = prefix+'.fits'
    #print 'Converting %s to %s'%(input,fits)
    if (os.path.exists(fits)):
	   os.unlink(fits)
    cmd = '%s/ndf2fits in=%s out=%s QUIET'%(convdir,input,fits)
    os.system(cmd)
    
    #set parameters for gaussian model
    func = lambda p,x: p[0]*np.exp(-(x/(np.sqrt(2)*p[1]))**2)

    def myfunct(p, fjac=None, x=None, y=None, err=None):
        model = func(p,x)
        # Non-negative status value means MPFIT should
        # continue, negative means stop the calculation.
        status = 0
        return [status, (y-model)/err]

#Set parameters for model fitting limits
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

    #Open file, extract data, build histogram, fit model using MPFIT
    #open data and variance
    print 'Producing histogram of data below '+str(snr_cut)+' sigma'
    hdulist = pyfits.open(fits)
    data = hdulist[0].data
    variance = hdulist[1].data
    #reduce dimensions if necessary
    if len(np.shape(data))>2:
        data = np.reshape(data,[np.shape(data)[1],np.shape(data)[2]])
        variance = np.reshape(variance,[np.shape(variance)[1],np.shape(variance)[2]])
    #mask data so that only noise dominated sections remain
    snr = np.abs(data/np.sqrt(variance))
    mask = np.where(snr<snr_cut,1,np.nan)#?????
    #deffine data/bins for use in histogram
    loc = np.where((snr<snr_cut) & (data<flux_lim[1]) & (data>flux_lim[0]) & (data!=0.0))
    bins=10000
    #calculate histogram
    hist, bin_edges = np.histogram((data[loc]), bins = bins)
    #format data/bins for use in MPFIT
    fa = {'x':bin_edges[:-1], 'y':hist, 'err':(np.zeros(np.shape(hist))+1)}
    #Fit model to data
    print 'fitting gaussian model to data'
    m = mpfit.mpfit(myfunct, p0, parinfo=param_info,functkw=fa,quiet=True)
    print 'STDV on the noise = '+str(m.params[1])+' Jy/pixel'
    if logic == 'TRUE':
        print "Plotting..."
        plt.clf()
        plt.scatter(bin_edges[:-1], hist)
        plt.plot(bin_edges[:-1],func(m.params,bin_edges[:-1]),c='r')
        plt.xlabel('Flux Density')
        plt.ylabel('Number of pixels')
        plt.show()
    else:
        print 'no plot'

    return m.params[1] #returns the STDV across all data of a 2sigma or below detection. 

if __name__ == "__main__":
    noise_by_data(input,logic)
