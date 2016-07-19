"""A module for calculating noise given and input map. This method masks all data above 2sigma and then fits a gaussian to the resulting histogram"""

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

def PARGET(file,parameter,kappa):
    #Function takes a file, opens it in either STATS or NDFTRACE and pulls out a specific parameter, returns it and bins the tempoary file.
    cmd = '%s/%s ndf=%s QUIET'%(kapdir,kappa,file)
    os.system(cmd)
    
    cmd = '%s/parget %s %s > parameter.txt'%(kapdir,parameter,kappa)
    os.system(cmd)

    a = np.loadtxt('parameter.txt')
    os.remove('parameter.txt')
    return a

def noise_by_data(input,logic): #this function calculates noise

    flux_lim =  [-0.5,0.5]#[-0.1,0.1]#[-0.02,0.02]
    snr_cut = 2.0

    #convert input maps into fits format
    prefix = string.split(input,'.sdf')[0]
    #print 'Converting %s to fits...'%(input)
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
 
#Set parameters for model fitting limits
    peak = hist.max()
    num = PARGET(input,'NUMPIX','stats')
    hist_limit = num*5E-4   #  /10000 
    
    histcut = []
    bincut = []
    for i in range(len(hist)):
        if hist[i] > hist_limit:
            bincut.append(bin_edges[i])
            histcut.append(hist[i])
    bincutnp = np.array(bincut)
    STDV = bincutnp.std()*0.5

    """
#IMPORTANT: This setup is only 'proven' to work with Aquila and MWC297. Other maps may require altering of the guess_sig values and how they vary between 450um and 850um.
    if wave == '450':
        guess_sig = 0.01
        p0 = [peak,guess_sig]
    elif wave == '850':
        guess_sig = 0.01
        p0 = [peak,guess_sig/3]
   """

    p0 = [peak,STDV]
    #print bincutnp

    #Initial guess of model {peak,mean,STDV,offset from x axis}
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


    #format data/bins for use in MPFIT
    fa = {'x':bin_edges[:-1], 'y':hist, 'err':(np.zeros(np.shape(hist))+1)}
    #Fit model to data
    print 'fitting gaussian model to data'
    m = mpfit.mpfit(myfunct, p0, parinfo=param_info,functkw=fa)
    print 'STDV on the noise = '+str(round(m.params[1],3))+' Jy/pixel'
    if logic == 'TRUE':
        print "Plotting..."
        plt.clf()
        plt.scatter(bin_edges[:-1], hist, color='b')
        plt.scatter(bincut, histcut,edgecolor='r',facecolor='r')
        plt.plot(bin_edges[:-1],func(m.params,bin_edges[:-1]),c='r')
        plt.xlabel('Flux Density (Jy/pix)')
        plt.ylabel('Number of pixels')
        plt.xlim([-0.05,0.05])#plt.xlim([-0.05,0.05])
        plt.ylim([0,p0[0]+1000])
        plt.show()
    else:
        print 'no plot'

    return m.params[1] #returns the STDV across all data of a 2sigma or below detection. 

def noise_by_snr(input,logic): #this function calculates noise for SNR maps (autoFW.py)

    flux_lim = [-10,10]#[-0.1,0.1]#[-0.02,0.02]
    snr_cut = 10.0

    #convert input maps into fits format
    prefix = string.split(input,'.sdf')[0]
    #print 'Converting %s to fits...'%(input)
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
 
#Set parameters for model fitting limits
    peak = hist.max()
    num = PARGET(input,'NUMPIX','stats')
    hist_limit = 5
    print hist_limit

    histcut = []
    bincut = []
    for i in range(len(hist)):
        if hist[i] > hist_limit:
            bincut.append(bin_edges[i])
            histcut.append(hist[i])
    bincutnp = np.array(bincut)
    STDV = bincutnp.std()*0.5

    """
#IMPORTANT: This setup is only 'proven' to work with Aquila and MWC297. Other maps may require altering of the guess_sig values and how they vary between 450um and 850um.
    if wave == '450':
        guess_sig = 0.01
        p0 = [peak,guess_sig]
    elif wave == '850':
        guess_sig = 0.01
        p0 = [peak,guess_sig/3]
   """

    p0 = [peak,STDV]
    #print bincutnp
    print p0

    #Initial guess of model {peak,mean,STDV,offset from x axis}
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


    #format data/bins for use in MPFIT
    fa = {'x':bin_edges[:-1], 'y':hist, 'err':(np.zeros(np.shape(hist))+1)}
    #Fit model to data
    print 'fitting gaussian model to data'
    m = mpfit.mpfit(myfunct, p0, parinfo=param_info,functkw=fa)
    print 'STDV on the noise = '+str(round(m.params[1],3))+' Jy/pixel'
    if logic == 'TRUE':
        print "Plotting..."
        plt.clf()
        plt.scatter(bin_edges[:-1], hist)
        plt.scatter(bincut, histcut,edgecolor='r',facecolor='r')
        plt.plot(bin_edges[:-1],func(m.params,bin_edges[:-1]),c='r')
        plt.xlabel('Flux Density')
        plt.ylabel('Number of pixels')
        #plt.xlim([-0.05,0.05])
        #plt.ylim([0,p0[0]+1000])
        plt.show()
    else:
        print 'no plot'

    return m.params[1] #returns the STDV across all data of a 2sigma or below detection. 


if __name__ == "__main__":
    noise_by_data(input,logic)
    noise_by_snr(input,logic)
