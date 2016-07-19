import numpy as np
from numpy import arange,array,ones, loadtxt#,random,linalg
from pylab import plot,show, xlabel, ylabel, errorbar, ylim, xlim
from scipy import stats
import sys

data = loadtxt('GBSprotostar_temperatures_master_FLT.tab')

l = len(data[:,0])

Class0_ET = []
Class0_eET = []
Class0_PT = []
Class0_ePT = []
ClassI_ET = []
ClassI_eET = []
ClassI_PT = []
ClassI_ePT = []

x = []
y = []
ex = []
ey = []

for j in data:
    if j[6] > 0:
        x.append(j[4])
        y.append(j[6])
        ex.append(j[5])
        ey.append(j[7])
        if j[9] < 70:
            Class0_ET.append(j[4])
            Class0_eET.append(j[5])
            Class0_PT.append(j[6])
            Class0_ePT.append(j[7])
        elif j[9] > 70:
            ClassI_ET.append(j[4])
            ClassI_eET.append(j[5])
            ClassI_PT.append(j[6])
            ClassI_ePT.append(j[7])


#for i in range(l):
#    print '=============================='
##    xi =  data[i:,4] 
#    yi = data[i:,6]
    
#    A = array([xi, ones(9)])

Wx = x * np.sqrt(ex)
Wy = y * np.sqrt(ey)

slope, intercept, r_value, p_value, slope_std_err = stats.linregress(x,y)
Wslope, Wintercept, Wr_value, Wp_value, Wslope_std_err = stats.linregress(Wx,Wy)

print 'intercept', intercept
print 'slope', slope
print 'standard deviation', slope_std_err
print 'Wintercept', Wintercept
print 'Wslope', Wslope
print 'Wstandard deviation', Wslope_std_err

line = [slope*i+intercept for i in x]
Wline = [Wslope*i+Wintercept for i in x]

plot(x,line,'k-')
plot(x,Wline,'k--')

X = [0,40]

plot(X,X,'r-',X,X,'+')

errorbar(ClassI_ET,ClassI_PT,xerr=ClassI_eET,yerr=ClassI_ePT,fmt='go')
errorbar(Class0_ET,Class0_PT,xerr=Class0_eET,yerr=Class0_ePT,fmt='bo')

ylim([5,40])
xlim([5,40])

xlabel('Envelope Temp. (K)')
ylabel('Protostar Temp. (K)')
show()

