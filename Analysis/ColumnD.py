"""Script for calculating the column density and number density - using rumble_mass"""

import numpy as np
import mass

#####################################################
#Function to deffine Column Density - Johnstone00 (850 only)

def massF(S,T):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""
    #Equation of Mass in solar masses, per pixel
    kappa = 0.012 #cm^2 g^-1
    d = 500 #pc
    if T == 0:
        T = 15
    M = 1.55*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/500.0)**2.0)
    print 'Mass = ',round(M,2),' Mdot'
    return M

def ColumnD(S,T,N):
    ##### THIS METHOD NOT KNOWN TO WORK #####
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), Number of pixels in distance (d in pcs)"""
    #constants
    au = 149597871000 #m
    M_x = 1.989E30 #kg
    m_h2 = 3.34745E-24 #g
    m_h = 1.67262178E-24 #g
    mu = 2.333333  #ratio of H2 to He (5:1)
    av = 1.8E21
    kappa = 0.012 #cm^2 g^-1
    d = 500 #pc

    #each pixel has an area in SI units and each apature contains N pixels
    a = 2.99987
    #a = 5.99987
    A = ((((a*d)*au)**2.0)*N)*10000 #in cm
    #M = mass.djr(S,T,kappa,d)
    M = massF(S,T)
    #print M
    m = M*M_x*1000.0 #A = 1 density in g per pixel

    n = m/(mu*m_h*A) #column density per cm^2 
    N = m/A #density in g per cm^2

    Av = n/av

    print 'mass per pixel in solar masses:', round(M,1)
    print 'column density of apature in per cm^2:', n
    print 'density of apature in g per cm^2:', round(N,4)
    print 'Extinction in mag.:', round(Av,0)
    return

def konyves(S,T):
    ##### THIS METHOD NOT KNOWN TO WORK #####
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin)"""
    #constants
    au = 149597871000 #m
    m_h = 1.67262178E-27 #kg
    mu = 2.8  #ratio of H2 to He (Kauffmann et al. 2008)
    h = 6.62606957E-34 #m2 kg / s
    c = 299792458 #m/s
    Kb = 1.3806488E-23 #m2 kg s-2 K-1
    av = 1.8E21 #H2 cm-2 mag-1    
    pc = 3.08567758E16 #m
    Jy = 1E-26 #Si
    #fixed variables
    nu = c/(850E-6) #m
    kappa = 0.012*((100**2.)/1000) #m^2 kg^-1
    d = 500 #pc
    a = 2.9999997 #arcsec
    N = 1
    #minor terms
    A = ((((a*d)*au)**2.0)*N)
    B = ((2.*h*(nu**3.))/c**2.)*(1/((np.exp((h*nu)/(Kb*T)))-1))
    s = S*Jy 
    #major terms
    N = (s*((d*pc)**2.))/(B*kappa*mu*m_h*A) #H2 m-2
    n = (s*((d*pc)**2.))/(B*kappa*A) #kg m-2

    print 'Column density [H2 cm-2] ',N/(100**2) #H2 cm-2
    print 'Column density [g cm-2] ',n/((100**2)/1000) #g cm-2
    print 'Extinction [mag] ',N/((100**2)*av)
    
if __name__ == "__main__":
    ColumnD(S,T,N)
    konyves(S,T)
    massF(S,T)
