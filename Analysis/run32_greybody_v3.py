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

def temp(T):
    beta = 1.8
    Sr = (850.**(3.+beta)) / (450.**(3.+beta)) * (np.exp((h*c)/(k*850*1E-6*T))-1.) / (np.exp((h*c)/(k*450*1E-6*T))-1.)
    return Sr

def GB(l,T,kappa,N):
    #equation of a Greybody function
    #B_a = ((2.0*h*(c**2.0))/(l**5.0))
    #B_b = (1.0/(np.exp((h*c)/(k*T*l))-1.0))
    #B = kappa*B_a*B_b*m_h2*mu*aN_h2*N

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

def kappa_test(T,l,S,M):
    d = 250*p
    return (S*(d**2.))/(BB(l,T)*M)

def flux(T,l,kappa,M):
    d = 250*p
 
    return (BB(l,T)*M*kappa)/(d**2.)

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
    return Rf #round(Rf,1)


#range of x
OH5 = np.loadtxt('oh5.tab',dtype='float')

lambda_si = []
kappa = []
sum_FLUX = []
L_MAX = []

Mf = np.arange(0,0.5,0.001)


for j in Mf:
    for i in range(len(OH5)):
        lambda_si.append((OH5[i][0])*1E-6)
        kappa.append(OH5[i][1])

        sum_FLUX.append(GB(lambda_si[i],15,kappa[i],1)+GB(lambda_si[i],50,kappa[i],j))
    
    items2 = {str(round(sum_FLUX[k],-7)):lambda_si[k] for k in range(len(OH5))}
    
    #print sum_FLUX
    AMX = max(sum_FLUX)
    #print AMX
    lambda_max = items2[str(round(AMX,-7))]

    L_MAX.append(AMX)
    
    
l850 = 850E-6
l450 = 450E-6
l1200 = 1.2E-3
l250 = 250E-6

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

T = 17.

print 'T = %s: RATIO= '%(T),temp(T)

f = np.arange(0.00,1,0.0001)

SR20 = [comb(20,i) for i in f]
SR25 = [comb(25,i) for i in f]
SR35 = [comb(35,i) for i in f]
SR50 = [comb(50,i) for i in f]
SR100 = [comb(100,i) for i in f]
SR500 = [comb(500,i) for i in f]
SR2000 = [comb(2000,i) for i in f]

#Deffine a 1D Plot using matplotlib
fig = plt.figure()

#Fraction

plt.loglog(f,SR20,'k',linestyle='-',linewidth=1)
#plt.plot(f,SR25,'k',linestyle='-',linewidth=1)
#plt.plot(f,SR35,'k',linestyle='-',linewidth=2)
plt.loglog(f,SR100,'k',linestyle='-',linewidth=2)
plt.loglog(f,SR500,'k',linestyle='-',linewidth=3)
plt.loglog(f,SR2000,'k',linestyle='-',linewidth=4)

plt.annotate('20K',xy=(0.2,0.1))
#plt.annotate('25K',xy=(0.15,0.05))
#plt.annotate('35K',xy=(0.15,0.14))
#plt.annotate('50K',xy=(0.0005,0.1))
plt.annotate('100K',xy=(0.005,0.1))
plt.annotate('500K',xy=(0.0009,0.1))
plt.annotate('2000K',xy=(0.0002,0.1))

plt.annotate('RATIO = 6 (~15K)',xy=(0.00012,0.8))

plt.axvline(x=(0.44),linestyle=':',color='k')
plt.axvline(x=(0.011),linestyle=':',color='k')
plt.axvline(x=(0.0017),linestyle=':',color='k')
plt.axvline(x=(0.00041),linestyle=':',color='k')
plt.axhline(y=(0.083),linestyle=':',color='k')

xlim([0.0001,1])
ylim([0.01,1])

#Apply labels
plt.ylabel(r'Fractional change in 450$\mu m$/850$\mu m$ flux ratio')
plt.xlabel(r"'Hot' to 'cold' cloud mass ratio")

#Peak Lambda
#plt.plot(Mf,L_MAX,'k',linestyle='-',linewidth=1)

#Apply labels
#plt.ylabel('Peak SED wavelength (m)')
#plt.xlabel('Hot to cold cloud mass fraction')

plt.show()
