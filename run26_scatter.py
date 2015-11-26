#20150721

#A script for producing scatter plots for paper2

import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt

def line(intercept,slope,x):
    X = np.array(x)
    y = (slope*X)+intercept
    return y

################################################################################
########## Bulk script ############
#NOTE this data set should NOT have the null values used in the histogram.
data =  np.loadtxt('SMM/SMM-Aq-ff-W4015.tab',dtype='float',comments="#")
#NOTE this data set should NOT have the null values used in the histogram.

#files
temp = np.array(data[:,8]).tolist()
errtemp = np.array(data[:,9]/2.).tolist()
YSO = np.array(data[:,12]).tolist()
errYSO = np.array(data[:,13]/2.).tolist()
#H70 = np.array(data[:,5]).tolist()
distance = np.array(data[:,18]).tolist()

#log
LOGYSO = np.log10(YSO)
LOGerrYSO = np.log10(errYSO)
#LOGH70 = np.log10(H70)

fig1 = plt.figure()

#Plot 1 distance
fig1.add_subplot(211)
xlabel('Distance (pc)')
ylabel('Temp. (K)')
x,y = distance,temp
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
plot(x,line(intercept,slope,x),'r')
errorbar(x,y,xerr=0,yerr=errtemp,fmt='b+')
ylim([5,28])
fig1.text(4,25,'R = %s'%(str(round(r_value,2))),size='x-small')
fig1.text(4,20,'P = %s'%(str(round(p_value,3))),size='x-small')
print round(r_value,2),round(p_value,3)

#Plot 2 H70
'''
fig1.add_subplot(312)
xlabel('Log 70$\mathrm{\mu}$m Flux density (MJy/Sr)')
ylabel('Temp. (K)')
x,y = LOGH70,temp
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
plot(x,line(intercept,slope,x),'r')
errorbar(LOGH70,temp,xerr=0,yerr=errtemp,fmt='b+')
ylim([5,28])
fig1.text(0.25,25,'R = %s'%(str(round(r_value,2))),size='x-small')
fig1.text(0.25,20,'P = %s'%(str(round(p_value,3))),size='x-small')
print round(r_value,2),round(p_value,3)
'''

#Plot 2 H70
fig1.add_subplot(212)
xlabel('Log YSO density (YSOs per pc$^{2}$)')
ylabel('Temp. (K)')
x,y,errx,erry = LOGYSO,temp,errtemp,LOGerrYSO
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
plot(x,line(intercept,slope,x),'r')
errorbar(LOGYSO,temp,xerr=errx,yerr=erry,fmt='b+')
ylim([5,28])
fig1.text(3,15,'R = %s'%(str(round(r_value,2))),size='x-small')
fig1.text(3,10,'P = %s'%(str(round(p_value,3))),size='x-small')
print round(r_value,2),round(p_value,3)

show()
