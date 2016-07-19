
import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt

import slfit

def line(intercept,slope,x):
    X = np.array(x)
    y = (slope*X)+intercept
    return y

def dust_temp(a,d,Lx):
    Go = 2.1E4*(Lx/(10**4.))*((0.1/d)**2.0)
    Td = 33.5*((1.E-6/a)**0.2)*((Go/(10**4.))**0.2)
    return Td

def sil_temp(a,d,Lx):
    Go = 2.1E4*(Lx/(10**4.))*((0.1/d)**2.0)
    Td = 50.*((1.E-6/a)**0.06)*((Go/(10**4.))**(1./6.))
    return Td

###################################################
#Plotting

d = np.arange(0.01,100,0.01) #Distance in parsecs



Lx_O5 = 200000
Lx_O9 = 55000
Lx_B0 = 24000
Lx_B1 = 5500
Lx_B3 = 1060
Lx_B4 = 650
Lx_B5 = 380
Lx_B9 = 42

a1 = 0.25E-6
a2 = 0.1E-6
a3 = 0.01E-6


dataO =  np.loadtxt('data/SMM_master_real_Otype.tab',dtype='string',comments="#")
dataB1 =  np.loadtxt('data/SMM_master_real_eBtype.tab',dtype='string',comments="#")
dataB2 =  np.loadtxt('data/SMM_master_real_lBtype.tab',dtype='string',comments="#")

#files
#p1
TempO = [np.log10(float(i)) for i in dataO[:,6]]
propO_a = [(float(i)) for i in dataO[:,7]]
propO_b = [i*np.log(10.) for i in propO_a]
ErrtempO = [a/b for a,b in zip(TempO,propO_b)]

TempB1 = [np.log10(float(i)) for i in dataB1[:,6]]
propB1_a = [(float(i)) for i in dataB1[:,7]]
propB1_b = [i*np.log(10.) for i in propB1_a]
ErrtempB1 = [a/b for a,b in zip(TempB1,propB1_b)]

TempB2 = [np.log10(float(i)) for i in dataB2[:,6]]
propB2_a = [(float(i)) for i in dataB2[:,7]]
propB2_b = [i*np.log(10.) for i in propO_b]
ErrtempB2 = [a/b for a,b in zip(TempB2,propB2_b)]

TempB = TempB1 + TempB2
ErrtempB = ErrtempB1 + ErrtempB2

DistanceO = [np.log10(float(i)) for i in dataO[:,19]]
DistanceB1 = [np.log10(float(i)) for i in dataB1[:,19]]
DistanceB2 = [np.log10(float(i)) for i in dataB2[:,19]]

DistanceB  = DistanceB1 + DistanceB2


distanceO = []
distanceB1 = []
distanceB2 = []
distanceB = []

tempO = []
errtempO = []
tempB1 = []
errtempB1 = []
tempB2 = []
errtempB2 = []
tempB = []
errtempB = []

for i in range(len(DistanceO)):
    D = DistanceO[i]
    if (D>=-1.5) & (D<=1.5):
        distanceO.append(D)
        tempO.append(TempO[i])
        errtempO.append(ErrtempO[i])
for i in range(len(DistanceB1)):
    D = DistanceB1[i]
    if (D>=-1.5) & (D<=1.5):
        distanceB1.append(D)
        tempB1.append(TempB1[i])
        errtempB1.append(ErrtempB1[i])
for i in range(len(DistanceB2)):
    D = DistanceB2[i]
    if (D>=-1.5) & (D<=1.5):
        distanceB2.append(D)
        tempB2.append(TempB2[i])
        errtempB2.append(ErrtempB2[i])
for i in range(len(DistanceB)):
    D = DistanceB[i]
    if (D>=-1.5) & (D<=1.5):
        distanceB.append(D)
        tempB.append(TempB[i])
        errtempB.append(ErrtempB[i])

#p2
#mtemp_yso = np.array(datatemp[:,3]).tolist()
#errmtemp_yso = np.array(datatemp[:,4]).tolist()
#mtemp_noyso = np.array(datatemp[:,1]).tolist()
#errmtemp_noyso = np.array(datatemp[:,2]).tolist()
#distance_temp = np.array(datatemp[:,0]).tolist()

label = ['O type','Early B type','Late B type']

#fig1 = plt.figure()

###########################
#fig1.add_subplot(121)
'''
ylabel('Temperature (K)')
xlabel('Log distance (pc)')
X = [-2,2]

x,y,ey = np.asarray(distanceB),np.asarray(tempB),np.asarray(errtempB)
errorbar(distanceB,tempB,xerr=0,yerr=errtempB,fmt='bo',ecolor='b')
#slopeB2, interceptB2, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'b')

a1 = 0.25E-6
a2 = 0.1E-6
a3 = 0.01E-6

Td = dust_temp(a2,d,Lx_B3)
logd = np.log10(d)
logTd = np.log10(Td)

plt.plot(logd,logTd,label='B3-0.1um',color='b',linestyle='--',linewidth=2)

ylim([0.5,2.0])
xlim([-2,2])
#plt.axvline(x=0.08,linestyle='--',color='k')

plt.axhline(y=1.176,linestyle='dashed',color='k')
legend()
'''
###########################
#fig1.add_subplot(122)
ylabel('Log temperature (K)')
xlabel('Log distance (pc)')
X = [-2,2]

x,y,ey = np.asarray(distanceO),np.asarray(tempO),np.asarray(errtempO)
errorbar(distanceO,tempO,xerr=0,yerr=errtempO,fmt='ro',ecolor='r')
#slopeO, interceptO, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'k')

a1 = 0.25E-6
a2 = 0.1E-6
a3 = 0.5E-6

Td1 = sil_temp(a1,d,Lx_O9)
Td2 = sil_temp(a2,d,Lx_O9)
Td3 = sil_temp(a3,d,Lx_O9)
logd = np.log10(d)
logTd1 = np.log10(Td1)
logTd2 = np.log10(Td2)
logTd3 = np.log10(Td3)

plt.plot(logd,logTd2,label='O9-0.1um',color='red',linestyle='-.',linewidth=2)
plt.plot(logd,logTd1,label='O9-0.25um',color='red',linestyle='--',linewidth=2)
plt.plot(logd,logTd3,label='O9-0.5um',color='red',linestyle=':',linewidth=2)


ylim([0.5,2.0])
xlim([-1.5,1.2])
#plt.axvline(x=0.52,linestyle=':',color='r')
#plt.axvline(x=0.08,linestyle='--',color='k')

plt.axhline(y=1.176,linestyle='dashed',color='k')
legend(loc=3)

show()
