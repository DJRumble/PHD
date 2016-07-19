#A script for plotting greybody model SEDs

from mpl_toolkits.axes_grid.inset_locator import inset_axes
import numpy as np
from pylab import *
import matplotlib.pyplot as plt

##### DEFFINE c0nstants ##########

au = 149597871000
p = 3.08567758E16
M_x = 1.989E30
h = 6.626068E-34
c = 2.99792458E8
k = 1.3806488E-23
m_h2 = 3.34745E-27 #g
G = 6.67384E-11
pi = 3.14159265359
mu = 2.333333  #ratio of H2 to He (5:1)

A = 5.03517447875E28 #m2
aN_h2 = 1E21 #cm-2

##### DEFFINE data - from equation  #########

def GB(l,T,kappa,N):
    #equation of a Greybody function
    d = 250.*p
    M = 2.38888003057e+38
    B_a = ((2.0*h*(c**2.0))/(l**5.0))
    B_b = (1.0/(np.exp((h*c)/(k*T*l))-1.0))
    B = kappa*B_a*B_b*M*N*(1/(d**2.))
    return B

def BB2(l,T,N):
    #equation of a Greybody function
    d = 250.*p
    M = 2.38888003057e+38
    B_a = ((2.0*h*(c**2.0))/(l**5.0))
    B_b = (1.0/(np.exp((h*c)/(k*T*l))-1.0))
    B =B_a*B_b*M*N*(1/(d**2.))
    return B

def BB(l,T):
    #equation of a blackbody function
    B_a = ((2.0*h*(c**2.0))/(l**5.0))
    B_b = (1.0/(np.exp((h*c)/(k*T*l))-1.0))
    B = B_a*B_b
    return B

def alpha(x1,x2,y1,y2):
    top = (np.log10(y2)-np.log10(y1))
    bottom = (np.log10(c/x2))-(np.log10(c/x1))
    a = top/bottom
    return a

def mass(T,l,S,kappa):
    d = 250*p
    return (S*(d**2.))/(BB(l,T)*kappa)

def comb(T,f):
    d = 250*p
    S450o = 6.0 #Jy
    S850o = 1.0 #Jy
    M450 = mass(15,l450,S450o,kappa_4)
    M850 = mass(15,l850,S850o,kappa_8)

    #print M850
    ###15K cloud
    S15_450 = (BB(l450,15)*M450*kappa_4)/(d**2.)
    S15_850 = (BB(l850,15)*M850*kappa_8)/(d**2.)
    ### TK cloud with fraction f
    STf_450 = (BB(l450,T)*M450*f*kappa_4)/(d**2.)
    STf_850 = (BB(l850,T)*M850*f*kappa_8)/(d**2.)
    #Sum
    c450 = S15_450 + STf_450
    c850 = S15_850 + STf_850
    #print round(c450,1),round(c850,1)
    #ratio 
    R = c450/c850
    #fractional change in ratio
    Rf = ((R/6.)-1)
    return R #round(Rf,1)


#range of x
OH5 = np.loadtxt('oh5.tab',dtype='float')

lambda_si = []
kappa = []

#Envelope
S_10K = []
S_15Ka = []
S_15Kb = []
S_20K = []

#Photosphere
S_200K = []
S_300K = []
S_2000K = []

sumA = []
sumB = []
sumC = []



for i in range(len(OH5)):
    lambda_si.append((OH5[i][0])*1E-6)
    kappa.append(OH5[i][1])
    
    S_10K.append(GB(lambda_si[i],15,kappa[i],1))
    S_15Ka.append(GB(lambda_si[i],15,kappa[i],1))
    S_15Kb.append(GB(lambda_si[i],15,kappa[i],0.1))
    S_20K.append(GB(lambda_si[i],20,kappa[i],0.001))

    S_200K.append(BB2(lambda_si[i],150,0.00001))
    S_300K.append(BB2(lambda_si[i],300,0.00001))
    S_2000K.append(BB2(lambda_si[i],800,0.00001))

    sumA.append(S_15Ka[i]+S_200K[i])
    sumB.append(S_15Kb[i]+S_300K[i])
    sumC.append(S_20K[i]+S_2000K[i])

l850 = 850E-6
l450 = 450E-6

#Deffine a 1D Plot using matplotlib
fig = plt.figure(figsize=(22, 6))

fig.add_subplot(141)

plt.loglog( lambda_si,S_10K,'k',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,sum_15_20_N10,'red',linestyle='--',linewidth=1)
#plt.loglog( lambda_si,S_50K_N10,'red',linestyle='-',linewidth=1)

plt.axvline(x=(l450),linestyle='-.',color='k')
plt.axvline(x=(l850),linestyle='-.',color='k')
plt.axvline(x=(2E-6),linestyle=':',color='k')
plt.axvline(x=(24E-6),linestyle=':',color='k')


plt.xlim([1E-6,3E-3])
#plt.ylim([1E-2,5E8])
plt.ylim([8E-2,9E4])

plt.annotate('Starless',color='k',xy=(3E-6,3E4))


#Apply labels
plt.ylabel('Relative flux density (arbitrary units)')
plt.xlabel('Wavelength (m)')

fig.add_subplot(142)

plt.loglog( lambda_si,S_15Ka,'k',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,sumA,'red',linestyle='--',linewidth=1)
#plt.loglog( lambda_si,S_200K,'red',linestyle='-',linewidth=1)

plt.axvline(x=(l450),linestyle='-.',color='k')
plt.axvline(x=(l850),linestyle='-.',color='k')
plt.axvline(x=(2E-6),linestyle=':',color='k')
plt.axvline(x=(24E-6),linestyle=':',color='k')

plt.xlim([1E-6,3E-3])
#plt.ylim([1E-2,5E8])
plt.ylim([8E-2,9E4])

plt.annotate('Class 0',color='k',xy=(3E-6,3E4))

#Apply labels
#plt.ylabel('Relative flux density (arbitrary units)')
plt.xlabel('Wavelength (m)')

fig.add_subplot(143)


plt.loglog( lambda_si,S_15Kb,'k',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,sumB,'red',linestyle='--',linewidth=1)
#plt.loglog( lambda_si,S_300K,'red',linestyle='-',linewidth=1)

plt.axvline(x=(l450),linestyle='-.',color='k')
plt.axvline(x=(l850),linestyle='-.',color='k')
plt.axvline(x=(2E-6),linestyle=':',color='k')
plt.axvline(x=(24E-6),linestyle=':',color='k')

plt.xlim([1E-6,3E-3])
#plt.ylim([1E-2,5E8])
plt.ylim([8E-2,9E4])

plt.annotate('Class I/FS',color='k',xy=(3E-6,3E4))

#Apply labels
#plt.ylabel('Relative flux density (arbitrary units)')
plt.xlabel('Wavelength (m)')

fig.add_subplot(144)

plt.loglog( lambda_si,S_20K,'k',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,sumC,'red',linestyle='--',linewidth=1)
#plt.loglog( lambda_si,S_2000K,'red',linestyle='-',linewidth=1)

plt.axvline(x=(l450),linestyle='-.',color='k')
plt.axvline(x=(l850),linestyle='-.',color='k')
plt.axvline(x=(2E-6),linestyle=':',color='k')
plt.axvline(x=(24E-6),linestyle=':',color='k')

plt.xlim([1E-6,3E-3])
#plt.ylim([1E-2,5E8])
plt.ylim([8E-2,9E4])

plt.annotate('Class II/III',color='k',xy=(3E-6,3E4))

#Apply labels
#plt.ylabel('Relative flux density (arbitrary units)')
plt.xlabel('Wavelength (m)')

plt.show()
