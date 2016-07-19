import astropy.io.fits as pyfits
import numpy as np
import random
import copy
import mpfit
import matplotlib.pyplot as plt
    

func = lambda p,x: p[0]*np.exp(-(x/(np.sqrt(2)*p[1]))**2)

def myfunct(p, fjac=None, x=None, y=None, err=None):
    model = func(p,x)
        # Non-negative status value means MPFIT should
        # continue, negative means stop the calculation.
    status = 0
    return [status, (y-model)/err]

#random dist generator.
data = []
for i in range(10000):
    k = random.normalvariate(1,0.1)
    data.append(k)

p0 = [1,0.1]

#param_base={'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
#param_info=[]
#for i in range(len(p0)):
#    param_info.append(copy.deepcopy(param_base))
#for i in range(len(p0)): 
#    param_info[i]['value']=p0[i]
#print param_info

#param_info[0]['limited'] = [1,1]       # Set the lower bound of parameter 0 to 'True' (i.e. 1)
#param_info[0]['limits']  = [0.0,100000.0]   # Set the lower limit of parameter 0 to 0.0
#param_info[1]['limited'] = [1,1]       # Set the lower bound of parameter 1 to 'True' (i.e. 1)
#param_info[1]['limits']  = [0.0,2.0] 

#loc = np.where((snr<snr_cut) & (data<flux_lim[1]) & (data>flux_lim[0]) & (data!=0.0))

bins=20
y, x = np.histogram(data, bins)

guessp = {1500,1,0.1,0}
fa = {'x':x[:-1], 'y':y, 'err':(np.zeros(np.shape(y))+1)}

m = mpfit.mpfit(myfunct, p0, parinfo=guessp,functkw=fa)
print 'show parameters:',m.params

#print 'STDV on the noise = '+str(m.params[1])+' Jy/pixel'
print "Plotting..."
plt.clf()
plt.scatter(x[:-1], y)
#plt.plot(x[:-1],func(m.params,x[:-1]),c='r')
plt.ylabel('freq')
plt.show()
