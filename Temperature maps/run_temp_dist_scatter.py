#20160713

#Damian Rumble, UoE,

#This script plots clump temperature as a function of distance to the major OB star in the region (pre deffined in the clump catalogue compilation). Clumps are split into three sets, O type, Early B type and late B type. Weighted lines are fit to each sample, with the sample cut at distances greater than 3pc.

import numpy as np
from pylab import *
from scipy import stats
import matplotlib.pyplot as plt
import time

import slfit

#model fitting parameters (not used)

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

def line(intercept,slope,x):
    X = np.array(x)
    y = (slope*X)+intercept
    return y

def dust_temp(a,d,Lx):
    #Model fitting (not used)
    Go = 2.1E4*(Lx/(10**4.))*((0.1/d)**2.0)
    Td = 33.5*((1.E-6/a)**0.2)*((Go/(10**4.))**0.2)
    return Td

def intercept(m,c,dm,dc):
    y = 15
    i = (y-c)/m
    di = np.sqrt((((-1./m)**2.)*(dc**2.))+((((y-c)/(m**2.))**2.)*(dm**2.)))

    I = 10.**(i)
    dI = np.sqrt(I*2.30258*(di**2.))

    print 'Intercept = '+str(round(I,2))+'pm'+str(round(dI,2))+' pc'

### DATA ###
dataO =  np.loadtxt('data/SMM_master_real_Otype.tab',dtype='string',comments="#")
dataB1 =  np.loadtxt('data/SMM_master_real_eBtype.tab',dtype='string',comments="#")
dataB2 =  np.loadtxt('data/SMM_master_real_lBtype.tab',dtype='string',comments="#")



#files
TempO = [float(i) for i in dataO[:,6]]
ErrtempO = [float(i) for i in dataO[:,7]]
TempB1 = [float(i) for i in dataB1[:,6]]
ErrtempB1 = [float(i) for i in dataB1[:,7]]
TempB2 = [float(i) for i in dataB2[:,6]]
ErrtempB2 = [float(i) for i in dataB2[:,7]]

DistanceO = [np.log10(float(i)) for i in dataO[:,19]]
DistanceB1 = [np.log10(float(i)) for i in dataB1[:,19]]
DistanceB2 = [np.log10(float(i)) for i in dataB2[:,19]]

distanceO = []
distanceB1 = []
distanceB2 = []
distanceOcut = []
distanceB1cut = []
distanceB2cut = []

tempO = []
errtempO = []
tempB1 = []
errtempB1 = []
tempB2 = []
errtempB2 = []
tempOcut = []
errtempOcut = []
tempB1cut = []
errtempB1cut = []
tempB2cut = []
errtempB2cut = []

Cut = 0.48

for i in range(len(DistanceO)):
    D = DistanceO[i]
    if (D>=-1.5) & (D<=1.5):
        distanceO.append(D)
        tempO.append(TempO[i])
        errtempO.append(ErrtempO[i])
    if (D>=-1.5) & (D<=Cut): #trim for better fit
        distanceOcut.append(D)
        tempOcut.append(TempO[i])
        errtempOcut.append(ErrtempO[i])
for i in range(len(DistanceB1)):
    D = DistanceB1[i]
    if (D>=-1.5) & (D<=1.5): #general
        distanceB1.append(D)
        tempB1.append(TempB1[i])
        errtempB1.append(ErrtempB1[i])
    if (D>=-1.5) & (D<=Cut): #trim for better fit
        distanceB1cut.append(D)
        tempB1cut.append(TempB1[i])
        errtempB1cut.append(ErrtempB1[i])
for i in range(len(DistanceB2)):
    D = DistanceB2[i]
    if (D>=-1.5) & (D<=1.5):#general
        distanceB2.append(D)
        tempB2.append(TempB2[i])
        errtempB2.append(ErrtempB2[i])
    if (D>=-1.5) & (D<=Cut): #trim for better fit
        distanceB2cut.append(D)
        tempB2cut.append(TempB2[i])
        errtempB2cut.append(ErrtempB2[i])

fig1 = plt.figure()

ylabel('Temperature (K)')
xlabel('Log distance (pc)')
X = [-2,2]

###O type stars###
x,y,ey = np.asarray(distanceOcut),np.asarray(tempOcut),np.asarray(errtempOcut)
errorbar(distanceO,tempO,xerr=0,yerr=errtempO,fmt='ro',ecolor='r',label=('O type'))

m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'r')

print 'O STARS'
print 'm = '+str(round(m,2))+' pm '+str(round(sigm,2))
print 'c = '+str(round(c,2))+' pm '+str(round(sigc,2))
print 'cov = '+str(round(cov,2))
intercept(m,c,sigm,sigc)

#d = np.arange(0.01,100,0.01) #Distance in parsecs
#Td = dust_temp(0.5E-6,d,Lx_O9)
#plt.plot(np.log10(d),Td,label=r'O9-0.5$\mu m$',color='r',linestyle='-.',linewidth=2)

###Early B type stars###
x,y,ey = np.asarray(distanceB1cut),np.asarray(tempB1cut),np.asarray(errtempB1cut)
errorbar(distanceB1,tempB1,xerr=0,yerr=errtempB1,fmt='yo',ecolor='y',label=('Early B type'))
#slopeB1, interceptB1, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'y')

print 'eB STARS'
print 'm = '+str(round(m,2))+' pm '+str(round(sigm,2))
print 'c = '+str(round(c,2))+' pm '+str(round(sigc,2))
print 'cov = '+str(round(cov,2))
intercept(m,c,sigm,sigc)

#d = np.arange(0.01,100,0.01) #Distance in parsecs
#Td = dust_temp(0.1E-6,d,Lx_B4)
#plt.plot(np.log10(d),Td,label='B4-0.1um',color='y',linestyle='--',linewidth=2)

###Late B type stars###
x,y,ey = np.asarray(distanceB2cut),np.asarray(tempB2cut),np.asarray(errtempB2cut)
errorbar(distanceB2,tempB2,xerr=0,yerr=errtempB2,fmt='bo',ecolor='b',label=('Late B type'))
#slopeB2, interceptB2, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'b')

print 'lB STARS'
print 'm = '+str(round(m,2))+' pm '+str(round(sigm,2))
print 'c = '+str(round(c,2))+' pm '+str(round(sigc,2))
print 'cov = '+str(round(cov,2))
intercept(m,c,sigm,sigc)

#d = np.arange(0.01,100,0.01) #Distance in parsecs
#Td = dust_temp(0.1E-6,d,Lx_B9)
#plt.plot(np.log10(d),Td,label='B9-0.1um',color='b',linestyle='--',linewidth=2)

#Plot format#

ylim([5,50])
xlim([-1.7,1.5])

plt.axvline(x=0.55,linestyle=':',color='r')
plt.axvline(x=-0.54,linestyle=':',color='y')
plt.axvline(x=-0.90,linestyle=':',color='b')
plt.axvline(x=Cut,linestyle='--',color='k')

plt.axhline(y=15,linestyle='dashed',color='k')
legend()

#Save fig specific date and time
date = str(time.strftime("%Y%m%d"))
savefig('plots/%s_temp_dist_JCMTGBS.pdf'%(date))

show()
