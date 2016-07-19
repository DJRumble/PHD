#A script for plotting greybody model SEDs

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
kappa = 0.01 #g/cm2
A = 5.03517447875E28 #m2

##### DEFFINE data - from equation  #########

def BB(nu,T):
    #equation of a blackbody function
    B_a = ((2.0*h*(nu**3.0))/(c**2.0))
    B_b = (1.0/(np.exp((h*nu)/(k*T))-1.0))
    B = B_a*B_b
    return B

def alpha(x1,x2,y1,y2):
    top = (np.log10(y2)-np.log10(y1))
    bottom = (np.log10(c/x2))-(np.log10(c/x1))
    a = top/bottom
    return a

#range of x

lambdamin = 1E-6
lambdamax = 5E-2

l850 = 850E-6
l450 = 450E-6
l1200 = 1.2E-3
l250 = 250E-6

T5 = [BB(c/l250,5),BB(c/l450,5),BB(c/l850,5),BB(c/l1200,5)]
T10 = [BB(c/l250,10),BB(c/l450,10),BB(c/l850,10),BB(c/l1200,10)]
T15 = [BB(c/l250,15),BB(c/l450,15),BB(c/l850,15),BB(c/l1200,15)]
T25 = [BB(c/l250,25),BB(c/l450,25),BB(c/l850,25),BB(c/l1200,25)]
T35 = [BB(c/l250,35),BB(c/l450,35),BB(c/l850,35),BB(c/l1200,35)]
T50 = [BB(c/l250,50),BB(c/l450,50),BB(c/l850,50),BB(c/l1200,50)]

#print alpha(l850,l1200,BB(c/l850,1000),BB(c/l1200,1000))

#print alpha(l850,l1200,T5[2],T5[3])
#print alpha(l850,l1200,T10[2],T10[3])
#print alpha(l850,l1200,T15[2],T15[3])
#print alpha(l850,l1200,T25[2],T25[3])
#print alpha(l850,l1200,T35[2],T35[3])
#print alpha(l850,l1200,T50[2],T50[3])


#print T5
#print T10
#print T15
#print T25
#print T35
#print T50


#print (c/lambdamax),(c/lambdamin)

x = np.arange((c/lambdamax),(c/lambdamin),1E10)

#Deffine a 1D Plot using matplotlib
plt.figure()

#plt.loglog(x, BB(x,5),'blue',linestyle='-',linewidth=1)
#plt.loglog(x, BB(x,10),'green',linestyle='-',linewidth=1)
plt.loglog(x, BB(x,15),'blue',linestyle='-',linewidth=1)
#plt.loglog(x, BB(x,25),'yellow',linestyle='-',linewidth=1)
plt.loglog(x, (0.1*BB(x,25)),'cyan',linestyle='-',linewidth=1)
plt.loglog(x, (0.05*BB(x,25)),'red',linestyle='-',linewidth=1)
plt.loglog(x, (0.01*BB(x,25)),'green',linestyle='-',linewidth=1)

print "fraction 'flux'at 450:"
print '10%'
K1a = ((0.1*BB(c/l450,35))+BB(c/l450,15))/BB(c/l450,15)
print K1a
print '5%'
K2a = ((0.05*BB(c/l450,35))+BB(c/l450,15))/BB(c/l450,15)
print K2a
print '1%'
K3a = ((0.01*BB(c/l450,35))+BB(c/l450,15))/BB(c/l450,15)
print K3a
print "TOTAL 'flux'at 850:"
print '10%'
K1b = ((0.1*BB(c/l850,35))+BB(c/l850,15))/BB(c/l850,15)
print K1b
print '5%'
K2b = ((0.05*BB(c/l850,35))+BB(c/l850,15))/BB(c/l850,15)
print K2b
print '1%'
K3b = ((0.01*BB(c/l850,35))+BB(c/l850,15))/BB(c/l850,15)
print K3b

K1 = [((0.1*BB(c/l450,35))+BB(c/l450,15)),((0.1*BB(c/l850,35))+BB(c/l850,15))]
K2 = [((0.07*BB(c/l450,35))+BB(c/l450,15)),((0.07*BB(c/l850,35))+BB(c/l850,15))]
K3 = [((0.02*BB(c/l450,35))+BB(c/l450,15)),((0.02*BB(c/l850,35))+BB(c/l850,15))]
X = [(c/l450),(c/l850)]

print X
print K2

plt.scatter(X,K1,color='cyan',marker='x')
plt.scatter(X,K2,color='red',marker='x')
plt.scatter(X,K3,color='green',marker='x')

#plt.loglog(x, (BB(x,20)),'b',linestyle=':',linewidth=1)
#plt.loglog(x, (BB(x,19)),'b',linestyle=':',linewidth=1)
#plt.loglog(x, (BB(x,18)),'b',linestyle=':',linewidth=1)
plt.loglog(x, (BB(x,17)),'b',linestyle=':',linewidth=1)
#plt.loglog(x, (BB(x,16)),'b',linestyle=':',linewidth=1)



plt.axvline(x=(c/l450),linestyle=':',color='k')
plt.axvline(x=(c/l850),linestyle=':',color='k')
plt.axvline(x=(c/l1200),linestyle='--',color='k')
plt.axvline(x=(c/l250),linestyle='--',color='k')

plt.xlim([2E11,8E12])
plt.ylim([5E-17,3E-15])

#plt.annotate('5K', xy=(3E11,1.5E-17),color='b')
#plt.annotate('10K', xy=(6E11,1E-16),color='g')
plt.annotate('15K', xy=(1E12,4E-16),color='blue')
plt.annotate('10%*50K', xy=(3E12,1.5E-15),color='cyan')
plt.annotate('5%*50K', xy=(3E12,7E-16),color='red')
plt.annotate('1%*50K', xy=(3E12,1.5E-16),color='green')
#plt.annotate('35K', xy=(2.5E12,6E-15),color='r')

plt.gca().invert_xaxis()

#Apply labels
plt.ylabel('Opacity modified Blackbody (W Sr-1 m-2 Hz-1)')
plt.xlabel('Frequency (Hz)')

plt.show()
