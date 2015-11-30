#05/01/13
#Damian Rumble, UoE
#This code is used to calcualte the the mass of a gas from its submil flux due to radiative transfer using:

# Mass = Flux x beam correction x d^2 / Planck function x oppacity

import math             #import relevant math packages 


#------------------------------------------------------------------------
#Global Parameters
#------------------------------------------------------------------------

h = 6.626068E-34
c_c = 299792458 
k = 1.3806503E-23
psc = 3.08567758E16
M_s = 1.989E30

#------------------------------------------------------------------------
#Planck Function
#------------------------------------------------------------------------

def planck(h,nu,t):
    B = (2.0*h*(nu**3.0))/((c_c**2.0)*(1.0/(math.exp((h*nu)/(k*t))-1.0)))
    return B

#------------------------------------------------------------------------
#Beam Correction Function
#------------------------------------------------------------------------

def correction(v,omega):
      A = ((v*3600.0)**2.0) / omega
      return A

#------------------------------------------------------------------------
#Parameters
#------------------------------------------------------------------------

l = 850E-6            #m
d = 130.0             #Pcs
omega = 238.23        #arc secs ^2
t = 10.0              #K
kappa = 0.00012       #cm^2/g
v = 0.0016666666666   #CDEL2 from JCMT website degs

#deffine input FLUX in Watts per m^2 per Hz per beam

S = 13E-26 

#------------------------------------------------------------------------
#Basic calculations
#------------------------------------------------------------------------

#find frequency - nu in Hz

nu = c_c / l

#find distance - D in m

D = d * psc;

#Corrected flux

s = S*correction(v,omega)

#Planck function in WHz^-1m^2

b = planck(h,nu,t)

#------------------------------------------------------------------------
#printing
#------------------------------------------------------------------------

print "flux (pre correction) = ",S," W m^-2 Hz^-1 per beam"
print "flux (post correction) = ",s," W m^-2 Hz^-1"
#print "correction = ",correction(v,omega)
print "nu = ",nu
#print "D = ",D
print "Planck function = ",b," W m^2 Hz^-1"
print "Oppacity = ",kappa," M^2/Kg"

#------------------------------------------------------------------------
#Mass calculations
#------------------------------------------------------------------------

m = (s*(D**2.0))/(b*kappa)
M = m/M_s

print "mass in kg = ",m
print "mass in solar masses = ",M
