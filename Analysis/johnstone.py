"""Script for caculating mass as calculated by Johnstone00b"""

import numpy as np

#####################################################
#Function to deffine Mass - Johnstone00 (850 only)
def johnstone(S,T,d,N):
    """Takes in puts of flux(S in Jys), temperature (T in Kelvin), distance (d in pcs) and Number of pixels in area (N)"""
    #constants
    au = 149597871000 #m
    kappa = 0.01 #g/cm2
    #each pixel has an area in SI units and each apature contains N pixels
    a = 5.99987
    A = (((a*d)*au)**2.0)*N
    #Equation of Mass in solar masses, per pixel
    M = 0.19*S*(np.exp(17.0/T)-1.0)*((kappa/0.01)**(-1.0))
    print 'mass of apature in solar masses:', M
    return

if __name__ == "__main__":
    johnstone(S,T,d,N)
