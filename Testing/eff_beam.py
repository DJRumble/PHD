#this code calculates the effective beam for a 4 component beam

import numpy as np

def FWHM(fwhm):
    return fwhm/(2.*np.sqrt(2.*np.log(2.)))

def quad(n1,n2,s1,s2):
    return (n1*n2*(s1**2)*(s2**2))/((s1**2)+(s2**2))

def quad_part(n1,n2,s1,s2):
    return (n1*n2*(s1**2)*(s2**2))


fwhm450_MB = 7.9
fwhm850_MB = 13.0 
fwhm450_S = 25.
fwhm850_S = 48.

sig4M = FWHM(fwhm450_MB)
sig8M = FWHM(fwhm850_MB)
sig4S = FWHM(fwhm450_S)
sig8S = FWHM(fwhm850_S)

print sig4M, sig4S
print sig8M, sig8S

#print sig4M, sig8M, sig4S, sig8S

a4 = 0.94
b4 = 0.06

a8 = 0.98
b8 = 0.02

#FWHM

#A = 1/(2.*np.pi*((a4*(fwhm450_MB**2.))+(b4*(fwhm450_S**2.))))
#B = 1/(2.*np.pi*((a8*(fwhm850_MB**2.))+(b8*(fwhm850_S**2.))))

#F = 4*(np.pi**2.)*A*B*(quad_part(a4,a8,fwhm450_MB,fwhm850_MB) + quad_part(a8,b4,fwhm450_S,fwhm850_MB) + quad_part(a4,b8,fwhm450_MB,fwhm850_S) + quad_part(b4,b8,fwhm450_S,fwhm850_S))

#I = 2*np.pi*A*B*(quad(a4,a8,fwhm450_MB,fwhm850_MB) + quad(a8,b4,fwhm450_S,fwhm850_MB) + quad(a4,b8,fwhm450_MB,fwhm850_S) + quad(b4,b8,fwhm450_S,fwhm850_S))

#SIGAM

A = 1/(2.*np.pi*(a4*(sig4M**2.)+b4*(sig4S**2.)))
B = 1/(2.*np.pi*(a8*(sig8M**2.)+b8*(sig8S**2.)))

F = 4*(np.pi**2.)*A*B*(quad_part(a4,a8,sig4M,sig8M) + quad_part(a8,b4,sig4S,sig8M) + quad_part(a4,b8,sig4M,sig8S) + quad_part(b4,b8,sig4S,sig8S))

I = 2*np.pi*A*B*(quad(a4,a8,sig4M,sig8M) + quad(a8,b4,sig4S,sig8M) + quad(a4,b8,sig4M,sig8S) + quad(b4,b8,sig4S,sig8S))


print 'Total flux at 450um ', F

print 'Intensity is ',I

print 'beam = ', np.sqrt(1/(2.*np.pi*I))*(2.*np.sqrt(2.*np.log(2.)))
