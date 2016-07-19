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
S_15K = []
S_17K = []
S_20K = []
S_35K = []
S_50K = []
S_100K = []

sum_15_50_N3 = []
sum_15_20_N10 = []

S_50K_N1 = []
S_50K_N3 = []
S_50K_N10 = []

for i in range(len(OH5)):
    lambda_si.append((OH5[i][0])*1E-6)
    kappa.append(OH5[i][1])
    
    S_15K.append(GB(lambda_si[i],15,kappa[i],1))
    S_17K.append(GB(lambda_si[i],17,kappa[i],1))
    S_20K.append(GB(lambda_si[i],20,kappa[i],1))
    S_35K.append(GB(lambda_si[i],35,kappa[i],1))
    S_50K.append(GB(lambda_si[i],50,kappa[i],1))
    S_100K.append(GB(lambda_si[i],100,kappa[i],1))

    S_50K_N1.append(GB(lambda_si[i],20,kappa[i],0.01))
    S_50K_N3.append(GB(lambda_si[i],20,kappa[i],0.03))
    S_50K_N10.append(GB(lambda_si[i],20,kappa[i],0.10))

    sum_15_50_N3.append(GB(lambda_si[i],15,kappa[i],1)+GB(lambda_si[i],20,kappa[i],0.03))
    sum_15_20_N10.append(GB(lambda_si[i],15,kappa[i],1)+GB(lambda_si[i],20,kappa[i],0.1))

l850 = 850E-6
l450 = 450E-6
l1200 = 1.2E-3

l250 = 250E-6
l350 = 350E-6
l160 = 160E-6
l70 = 70E-6


#Dictionary
items = {str(lambda_si[i]):kappa[i] for i in range(len(OH5))}

kappa0 = items['0.0005']

kappa_8 = kappa0*((l850/500E-6)**(-1.8))
kappa_8_d = kappa0*((l850/500E-6)**(-1.0))
kappa_4 = kappa0*((l450/500E-6)**(-1.8))
kappa_4_d = kappa0*((l450/500E-6)**(-1.0))


print kappa_4, kappa_8
print kappa_4_d, kappa_8_d
print kappa_4/kappa_4_d, kappa_8/kappa_8_d


L = [l450,l850]

S450 = GB(L[0],15,kappa_4,1)+GB(L[0],20,kappa_4,0.03)
S850 = GB(L[1],15,kappa_8,1)+GB(L[1],20,kappa_8,0.03)

S450o = GB(L[0],15,kappa_4,1)+GB(L[0],50,kappa_4,0.00)
S850o = GB(L[1],15,kappa_8,1)+GB(L[1],50,kappa_8,0.00)

#print S450o,S850o,S450o/S850o

print "RATIOS"
print "0%:",comb(20,0.0)
print "3%:",comb(20,0.03)

#Deffine a 1D Plot using matplotlib
fig = plt.figure()

#plt.loglog( lambda_si,S_15K,'blue',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_20K,'green',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_35K,'yellow',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_50K,'orange',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_100K,'red',linestyle='-',linewidth=1)

plt.loglog( lambda_si,S_15K,'k',linestyle='-',linewidth=1)
plt.loglog( lambda_si,sum_15_20_N10,'red',linestyle='--',linewidth=1)
#plt.loglog( lambda_si,S_17K,'k',linestyle='--',linewidth=1)

#plt.loglog( lambda_si,kappa,'red',linestyle='-',linewidth=1)

#plt.loglog( lambda_si,S_50K_N1,'green',linestyle='-',linewidth=1)
plt.loglog( lambda_si,S_50K_N10,'red',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_50K_N10,'cyan',linestyle='-',linewidth=1)

plt.axvline(x=(l450),linestyle='-.',color='k')
plt.axvline(x=(l850),linestyle='-.',color='k')
plt.axvline(x=(l1200),linestyle=':',color='k')
plt.axvline(x=(l350),linestyle=':',color='k')
plt.axvline(x=(l250),linestyle=':',color='k')
plt.axvline(x=(l160),linestyle=':',color='k')
plt.axvline(x=(l70),linestyle=':',color='k')

plt.xlim([1E-5,3E-3])
#plt.ylim([1E-2,5E8])
plt.ylim([8E-2,9E4])

#Apply labels
plt.ylabel('Relative flux density (arbitrary units)')
#plt.ylabel('Dust opactiy (cm$^{2}$ per g gas+dust)')
#plt.xlabel('Frequency (Hz)')
plt.xlabel('Wavelength (m)')

#plt.loglog( lambda_si,S_50K_N1,'green',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_50K_N3,'red',linestyle='-',linewidth=1)
#plt.loglog( lambda_si,S_50K_N10,'cyan',linestyle='-',linewidth=1)

plt.annotate('15K',color='k',xy=(1.5E-3,1E-1))
#plt.annotate('~17K',color='b',xy=(7E-5,2E3))
#plt.annotate('50K @ 3%',color='r',xy=(1.5E-5,5E4))
plt.annotate('20K @ 10%',color='r',xy=(4.0E-4,1E-1))

#plt.axvline(x=(l450),linestyle='.',color='k')
#plt.axvline(x=(l850),linestyle='.',color='k')
#plt.axvline(x=(l1200),linestyle='--',color='k')
#plt.axvline(x=(l250),linestyle='--',color='k')

plt.show()
