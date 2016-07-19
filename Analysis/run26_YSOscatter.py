
import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt

import slfit

def line(intercept,slope,x):
    X = np.array(x)
    y = (slope*X)+intercept
    return y

dataO =  np.loadtxt('data/SMM_master_real_noOB.tab',dtype='string',comments="#")

datalo20 =  np.loadtxt('data/SMM_master_real_noOB_lo20.tab',dtype='string',comments="#")

datahi20 =  np.loadtxt('data/SMM_master_real_noOB_hi20.tab',dtype='string',comments="#")

datalo60 =  np.loadtxt('data/SMM_master_real_noOB_lo60.tab',dtype='string',comments="#")

datahi60 =  np.loadtxt('data/SMM_master_real_noOB_hi60.tab',dtype='string',comments="#")

#files
#p1
TempO = [float(i) for i in dataO[:,6]]
ErrtempO = [float(i) for i in dataO[:,7]]
TempB2 = [float(i) for i in datalo20[:,6]]
ErrtempB2 = [float(i) for i in datalo20[:,7]]
TempB4 = [float(i) for i in datahi20[:,6]]
ErrtempB4 = [float(i) for i in datahi20[:,7]]
#TempB2 = [float(i) for i in datalo60[:,6]]
#ErrtempB2 = [float(i) for i in datalo60[:,7]]
#TempB4 = [float(i) for i in datahi60[:,6]]
#ErrtempB4 = [float(i) for i in datahi60[:,7]]

logTempO = [np.log10(float(i)) for i in dataO[:,6]]
propO_a = [(float(i)) for i in dataO[:,7]]
propO_b = [i*np.log(10.) for i in propO_a]
logErrtempO = [a/b for a,b in zip(TempO,propO_b)]

YSOO = [float(i) for i in dataO[:,14]]

#print YSOO

ysoO = []

tempO = []
errtempO = []
tempB1 = []
errtempB1 = []
tempB2 = []
errtempB2 = []

for i in range(len(YSOO)):
    D = YSOO[i]
    if (D>0) & (D<500):
        #print D
        ysoO.append(np.log10(D))
        tempO.append(TempO[i])
        errtempO.append(ErrtempO[i])

#label = ['n$_{\mathrm{YSO}}$ > 20 YSO pc$^{-2}$',r' n$_{\mathrm{YSO}}$ < 20 YSO pc$^{-2}$']
label = [r' n$_{\mathrm{YSO}}$ > 60 YSO pc$^{-2}$',r' n$_{\mathrm{YSO}}$ < 60 YSO pc$^{-2}$']

fig1 = plt.figure()

b = 25.

fig1.add_subplot(121)
nB4, bin_edgesB4, patchesB4 = hist(TempB4,bins=b/2.,normed=True, histtype='step',label=str(label[0]),color='r')
#nB1, bin_edgesB1, patchesB1 = hist(TempB1,bins=b,normed=True, histtype='step',label=str(label[1]),color='g')
nB2, bin_edgesB2, patchesB2 = hist(TempB2,bins=b,normed=True, histtype='step',label=str(label[1]),color='b')

D,P = stats.ks_2samp(nB4,nB2)
print 'KS-stats: D = ',D,' P = ',P

ylabel('Normalised frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([5,30])
#ylim([0,0.17])

fig1.add_subplot(122)
ylabel('Temperature (K)')
xlabel(r'Log YSO density (YSOs per pc$^2$)')
#X = [0,200]
X = [-2,3]
x,y = ysoO,tempO
errorbar(x,y,xerr=0,yerr=errtempO,fmt='ko',ecolor='k')

#slopeO, interceptO, r_value1, p_value1, std_err1 = stats.linregress(x,y)
#plot(X,line(interceptO,slopeO,X),'r')

x,y,ey = np.asarray(ysoO),np.asarray(tempO),np.asarray(errtempO)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
#plot(X,line(c,m,X),'r')


ylim([5,30])
xlim([-3,3])

plt.axhline(y=15,linestyle='dashed',color='k')
plt.axhline(y=13,linestyle=':',color='k')
plt.axhline(y=17,linestyle=':',color='k')
plt.axvline(x=1.3,linestyle='dashed',color='k')
plt.axvline(x=1.8,linestyle='dashed',color='k')
#plt.axvline(x=1.69,linestyle='dashed',color='k')
legend()

show()
