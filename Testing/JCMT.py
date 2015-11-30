#Damian Rumble, UoE
#06/01/2014
#equ_plotter.py

#A sand box for plotting results - from EQUATIONS

import numpy as np
#from pylab import *
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

##### DEFFINE labels ##########



##### DEFFINE data - from equation  #########

def BB(T):
    l = 850E-6
    nu = c / l
    #equation of a blackbody function
    B_a = ((2.0*h*(nu**3.0))/(c**2.0))
    B_b = (1.0/(np.exp((h*nu)/(k*T))-1.0))
    B = B_a*B_b
    return B


def f(t):

    return ((((pi**2.0)*k)/(6*G*mu*m_h2))*t*2.53199195516E14)/M_x

def g(s,t):
    
    return 0.23*s*(np.exp(17/t)-1.0)*((kappa/0.02)**(-1.0))*((250.0/250.0)**2.0)

def H(T,M):
    B = BB(T)
    return (B*M*M_x*kappa*10.0/((250.0*p)**2.0))*1E27

def I(T,beta): #temp (given beta


#range of x
#x = np.arange(0,100,1)

#Deffine a 1D Plot using matplotlib
#plt.figure()

#Tile plots here
#plt.subplot(211)


#Plot equations here (from function)
#plt.plot(x, I(x,1.6),'blue',linestyle=':',linewidth=4)
#plt.plot(x, I(x,2.0),'red',linestyle=':',linewidth=4)
#plt.plot(EQU1[3],EQU1[4],'red',linestyle='-')
#plt.plot(EQU2[3],EQU2[4],'blue',linestyle='-')
#plt.plot(EQU1[3],EQU1[5],'red',linestyle='--')
#plt.plot(EQU2[3],EQU2[5],'blue',linestyle='--')
#Annotate plots here
#plt.annotate('Beta = 1.0', xy=(40, 5))

#Apply labels
#plt.ylabel('ratio')
#plt.xlabel('temperature')

    return (17.**(3.+beta)) / (9.**(3.+beta)) * (np.exp(16.93/T)-1.) / (np.exp(31.97/T)-1.)




def J(S,T): #Mass error
    #return fmt: [variable,xaxis,yaxis,save
    dS450 = 0.3
    dS850 = 0.02
    A = (0.39**2.)*(((np.exp(17./T)-1)**2.)*(dS450**2.))
    B = (0.39**2.)*(289.*(S**2.)/(T**4.))*np.exp(34./T)*((0.05*T)**2.)
    
    M = A + B
    dM = np.sqrt(M)
    print M

    #labels
    X1 = 'flux(Jy)'
    Y1 = 'mass err'
    
    #save file
    save = "/images/20140106_masserr1.pdf"

    return dM,X1,Y1,save

def beam1D(wave,x):
    #return fmt: [variable,xaxis,yaxis,x-axis .txt data,y-axis .txt data]
    y = 0

    beam = np.loadtxt('beam.txt')
    xtxt = beam[:,0]

    #print wave

    if wave == '450':
        altpix = 6.
        sigM = 7.9/(2.3548200450309493*altpix)
        sigS = 25.0/(2.3548200450309493*altpix)
        a = 0.94
        b = 0.06
        norm = 1 /(2.*np.pi*((a*(sigM**2.))+(b*(sigS**2.))))
        ytxt = beam[:,1]
        ytxt_old = beam[:,3]
        Io = beam[15][1]
    elif wave == '850':
        altpix = 4.
        sigM = 13.0/(2.3548200450309493*altpix)
        sigS = 48.0/(2.3548200450309493*altpix)
        a=0.98
        b=0.02
        norm = 1 /(2.*np.pi*((a*(sigM**2.))+(b*(sigS**2.))))
        ytxt = beam[:,2]
        ytxt_old = beam[:,4]
        Io = beam[15][2]
    else:
        sigM = 1.
        sigS = 1.
        altpix = 1.
        print 'Not valid WAVE'

    partM = (-((x**2.)+(y**2.))/(2.*((sigM/altpix)**2.)))
    partS = (-((x**2.)+(y**2.))/(2.*((sigS/altpix)**2.)))

    F = norm*((a*np.exp(partM))+(b*np.exp(partS)))

    #labels
    X1 = 'x'
    Y1 = 'intensity'

    #save file
    save = "/images/20140714beam"+str(wave)+".pdf"
    

    return F,X1,Y1,xtxt,ytxt,ytxt_old,Io

###################################################

#range of x
x = np.arange(-5,5,0.1)

#Deffine a 1D Plot using matplotlib
plt.figure()

#Tile plots here
#plt.subplot(211)

#Deffine which equation I am using
EQU1 = beam1D('450',x)
EQU2 = beam1D('850',x)

print 'effective beam:'
print 'EQU:'
print '450 = ',(np.sqrt(1/(2.*np.pi*beam1D('450',0)[0])))*2.3548200450309493*6
print '850 = ',(np.sqrt(1/(2.*np.pi*beam1D('850',0)[0])))*2.3548200450309493*4
print 'SCR:'
print '450 = ',(np.sqrt(1/(2.*np.pi*EQU2[6])))*2.3548200450309493*6
print '850 = ',(np.sqrt(1/(2.*np.pi*EQU1[6])))*2.3548200450309493*4

#Plot equations here (from function)
plt.plot(x*6, EQU1[0],'blue',linestyle=':',linewidth=4)
plt.plot(x*4, EQU2[0],'red',linestyle=':',linewidth=4)
plt.plot(EQU1[3],EQU1[4],'red',linestyle='-')
plt.plot(EQU2[3],EQU2[4],'blue',linestyle='-')
#plt.plot(EQU1[3],EQU1[5],'red',linestyle='--')
#plt.plot(EQU2[3],EQU2[5],'blue',linestyle='--')

#plt.title(T)
plt.grid(True)

# Display plots
plt.show()
#plt.savefig(EQU1[3])
