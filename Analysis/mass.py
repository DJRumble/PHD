"""Script for caculating mass as calculated by:
- Johnstone00b
- Kirk06
- Sadavoy10
- Rumble13"""

import numpy as np
import jeans

#####################################################
#Function to deffine Mass - Johnstone00 (850 only)
def johnstone(S,T,kappa,d):
    """Takes in puts of flux(S in Jys), temperature (T in Kelvin),  oppacity (kappa in g cm-2) and distance (d in pcs)"""
    #each pixel has an area in SI units and each apature contains N pixels
    a = 5.99987
    #Equation of Mass in solar masses, per pixel
    M = 0.19*S*(np.exp(17.0/T)-1.0)*((kappa/0.01)**(-1.0))
    print 'mass of apature in solar masses:', M
    return M

#####################################################
#Function to deffine Mass - Kirk06 (850 only)
def kirk(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""

    #Equation of Mass in solar masses, per pixel
    M = 0.23*S*(np.exp(17.0/T)-1.0)*((kappa/0.02)**(-1.0))*((d/250.0)**2.0)
    print 'mass of apature in solar masses:', M
    return M

#####################################################
#Function to deffine Mass - djr13 (850 only)
def djr(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""

    if (T == 0.):
        M = 0
    else:
        #Equation of Mass in solar masses, per pixel
        M = 0.39*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/250.0)**2.0)
    #print 'mass of apature in solar masses:', M
    return M


#####################################################
#Function to deffine Mass - djr13 (850 only)
def djr_err(S,dS,T,dT,N):
    """Produces error in M/Mj = Takes inputs of; Flux (S in Jy), Flux error (S in Jy), Temperature (T in Kelvin) and Temperature error (T in Kelvin)"""
    
    A =(((np.exp(34./T)-(2.*np.exp(17./T))+1.)*(dS**2.)))
    B = (((17.*(S/T**2.))**2.)*(np.exp(34./T)*(dT**2.)))
    C = (0.39**2)*(A+B)
    dM = np.sqrt(C)
    M = djr(S,T,0.012,250)
    Mj = jeans.sadavoy(T,250,N)

    frac_m = dM/M
    frac_mj = dT/T
    print A, B, C
    print 'frac error in M =', frac_m, 'and frac error in Mj =', frac_mj

    
    R = M/Mj

    dR = R*np.sqrt((frac_m**2)+(frac_mj**2))
    print 'R =', R
    print 'error in R =', dR
    print 'Frac error in R =', dR/R

    return dM

#####################################################
#Function to deffine Mass - sadavoy10 (850 only)
def sadavoy(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""

    #Equation of Mass in solar masses, per pixel
    M = 0.074*S*(np.exp(17.0/T)-1.0)*((kappa/0.01)**(-1.0))*((d/100.0)**2.0)
    print 'mass of apature in solar masses:', M
    return M


if __name__ == "__main__":
    johnstone(S,T,d,N)
    djr(S,T,kappa,d)
    kirk(S,T,kappa,d)
    sadavoy(S,T,kappa,d)
