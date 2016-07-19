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
data =  np.loadtxt('data/SMM_master_real.tab',dtype='string',comments="#")
#data0 =  np.loadtxt('SMM/S2-FF-CO/SMM-15-0.tab',dtype='float',comments="#") #YSO plot
#dataYSO =  np.loadtxt('SMM/S2-FF-CO/SMM-15-YSO.tab',dtype='float',comments="#")


#NOTE this data set should NOT have the null values used in the histogram.

#files
#p1
temp = np.array(float(data[:,6])).tolist()
errtemp = np.array(float(data[:,7])/1.).tolist()

distance = np.array(data[:,19]).tolist()
distanceY = np.array(dataYSO[:,19]).tolist()

#p2
#temp
#tempY
#errtemp
#errtempY
H70 = np.array(data[:,5]).tolist()
LOGH70 = np.log10(H70)

#p3
#temp
#tempY
#errtemp
#errtempY
S21 = np.array(data[:,6]).tolist()
S21Y = np.array(dataYSO[:,6]).tolist()
LOGS21 = np.log10(S21)

#p4
#temp
#tempY
YSO = np.array(data0[:,13]).tolist()
YSOY = np.array(dataYSO[:,13]).tolist()
temp0 = np.array(data0[:,9]).tolist()
errtemp0 = np.array(data0[:,10]/1.).tolist()
LOGYSO = np.log10(YSO)
LOGYSOY = np.log10(YSOY)


fig1 = plt.figure()

#Plot 1 distance
#fig1.add_subplot(311)
xlabel('Distance (pc)')
ylabel('Temp. (K)')
x,y,xY,yY = distance,temp,distanceY,tempY
errorbar(x,y,xerr=0,yerr=errtemp,fmt='b+')
errorbar(xY,yY,xerr=0,yerr=errtempY,fmt='bo')
ylim([5,40])
xlim([0,5])
plt.axvline(x=1.2,linestyle='dashed',color='red')

#Plot 2 H70
'''
fig1.add_subplot(412)
xlabel('Log 70$\mathrm{\mu}$m Flux density (MJy/Sr)')
ylabel('Temp. (K)')
x,y = LOGH70,temp
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#plot(x,line(intercept,slope,x),'r')
errorbar(LOGH70,temp,xerr=0,yerr=errtemp,fmt='b+')
ylim([5,40])
xlim([0,4])
fig1.text(0.25,25,'R = %s'%(str(round(r_value,2))),size='x-small')
fig1.text(0.25,20,'P = %s'%(str(round(p_value,3))),size='x-small')
print round(r_value,2),round(p_value,6)

#Plot 3 s21

fig1.add_subplot(211)
xlabel('21cm Flux density (Jy/pix)')
ylabel('Temp. (K)')
x,y,xY,yY = S21,temp,S21Y,tempY
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#plot(x,line(intercept,slope,x),'r')
errorbar(S21,temp,xerr=0,yerr=errtemp,fmt='wo',ecolor='k')
errorbar(xY,yY,xerr=0,yerr=errtempY,fmt='bo')
ylim([5,40])
#xlim([-4.2,-0.7])
print round(r_value,2),round(p_value,6)
plt.axvline(x=0.01,linestyle='dashed',color='red')

#Plot 4 YSO
fig1.add_subplot(212)
xlabel('YSO density (YSOs per pc$^{2}$)')
ylabel('Temp. (K)')
x,y,erry,erryY = YSO,temp0,errtemp0,errtempY
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#plot(x,line(intercept,slope,x),'r')
errorbar(x,y,yerr=erry,fmt='wo',ecolor='k')
errorbar(YSOY,tempY,yerr=erryY,fmt='bo')
ylim([5,40])
#xlim([0.4,2.5])
print round(r_value,2),round(p_value,6)
plt.axvline(x=45,linestyle='dashed',color='red')
'''
show()
