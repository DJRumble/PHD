
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

def intercept(m,c,dm,dc):
    y = 15
    i = (y-c)/m
    di = np.sqrt((((-1./m)**2.)*(dc**2.))+((((y-c)/(m**2.))**2.)*(dm**2.)))

    I = 10.**(i)
    dI = np.sqrt(I*2.30258*(di**2.))

    print 'Intercept = '+str(I)+'pm'+str(dI)+' pc'

dataO =  np.loadtxt('data/SMM_master_real_Otype.tab',dtype='string',comments="#")
dataB1 =  np.loadtxt('data/SMM_master_real_eBtype.tab',dtype='string',comments="#")
dataB2 =  np.loadtxt('data/SMM_master_real_lBtype.tab',dtype='string',comments="#")

#files
#p1
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

#p2
#mtemp_yso = np.array(datatemp[:,3]).tolist()
#errmtemp_yso = np.array(datatemp[:,4]).tolist()
#mtemp_noyso = np.array(datatemp[:,1]).tolist()
#errmtemp_noyso = np.array(datatemp[:,2]).tolist()
#distance_temp = np.array(datatemp[:,0]).tolist()

label = ['O type','Early B type','Late B type']

fig1 = plt.figure()

fig1.add_subplot(211)
nO, bin_edgesO, patchesO = hist(tempO,bins=10,normed=True, histtype='step',label=str(label[0]),color='r')
nB1, bin_edgesB1, patchesB1 = hist(tempB1,bins=10,normed=True, histtype='step',label=str(label[1]),color='y')
nB2, bin_edgesB2, patchesB2 = hist(tempB2,bins=10,normed=True, histtype='step',label=str(label[2]),color='b')

D,P = stats.ks_2samp(nO,nB1)
print 'KS-stats: D = ',D,' P = ',P
D,P = stats.ks_2samp(nB1,nB2)
print 'KS-stats: D = ',D,' P = ',P



ylabel('Normalised frequency')
xlabel('Pixel temperature (K)')
legend()
xlim([5,47])
ylim([0,0.125])

fig1.add_subplot(212)
ylabel('Temperature (K)')
xlabel('Log distance (pc)')
X = [-2,2]

x,y,ey = np.asarray(distanceOcut),np.asarray(tempOcut),np.asarray(errtempOcut)
errorbar(distanceO,tempO,xerr=0,yerr=errtempO,fmt='ro',ecolor='r')
#slopeO, interceptO, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'r')

print 'O STARS'
print 'm = '+str(m)+' pm '+str(sigm)
print 'c = '+str(c)+' pm '+str(sigc)
print 'cov = '+str(cov)
intercept(m,c,sigm,sigc)

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

d = np.arange(0.01,100,0.01) #Distance in parsecs
Td = dust_temp(0.5E-6,d,Lx_O9)
plt.plot(np.log10(d),Td,label=r'O9-0.5$\mu m$',color='r',linestyle='-.',linewidth=2)

x,y,ey = np.asarray(distanceB1cut),np.asarray(tempB1cut),np.asarray(errtempB1cut)
errorbar(distanceB1,tempB1,xerr=0,yerr=errtempB1,fmt='yo',ecolor='y')
#slopeB1, interceptB1, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'y')

print 'O STARS'
print 'm = '+str(m)+' pm '+str(sigm)
print 'c = '+str(c)+' pm '+str(sigc)
print 'cov = '+str(cov)
intercept(m,c,sigm,sigc)

d = np.arange(0.01,100,0.01) #Distance in parsecs
Td = dust_temp(0.1E-6,d,Lx_B4)
#plt.plot(np.log10(d),Td,label='B4-0.1um',color='y',linestyle='--',linewidth=2)

x,y,ey = np.asarray(distanceB2cut),np.asarray(tempB2cut),np.asarray(errtempB2cut)
errorbar(distanceB2,tempB2,xerr=0,yerr=errtempB2,fmt='bo',ecolor='b')
#slopeB2, interceptB2, r_value1, p_value1, std_err1 = stats.linregress(x,y)
m,c,sigm,sigc,cov = slfit.slfit(x,y,ey) #Flux Wieghting Module
plot(X,line(c,m,X),'b')

print 'O STARS'
print 'm = '+str(m)+' pm '+str(sigm)
print 'c = '+str(c)+' pm '+str(sigc)
print 'cov = '+str(cov)
intercept(m,c,sigm,sigc)

d = np.arange(0.01,100,0.01) #Distance in parsecs
Td = dust_temp(0.1E-6,d,Lx_B9)
#plt.plot(np.log10(d),Td,label='B9-0.1um',color='b',linestyle='--',linewidth=2)

ylim([5,47])
xlim([-1.7,1.5])
plt.axvline(x=0.55,linestyle=':',color='r')
plt.axvline(x=-0.54,linestyle=':',color='y')
plt.axvline(x=-0.90,linestyle=':',color='b')
plt.axvline(x=Cut,linestyle='--',color='k')

plt.axhline(y=15,linestyle='dashed',color='k')
legend()

show()
