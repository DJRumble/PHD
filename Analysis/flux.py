"""Script for calculating the fluxmass in simplest form"""

import numpy as np
import BlackBody

#####################################################
#Function to deffine Flux 
def flux(T,M,l,kappa,d):
    """Takes inputs of; Temperature (T in Kelvin), Mass (M in solaremasses),  wavelength (l in m), oppacity (kappa in g cm-2), distance (d in pcs)"""
    #constants
    M_x = 1.989E30 #solar mass
    p = 3.08567758E16 #parsec

    # Temp is an array, mass is single float
    B = BlackBody.BB(l,T)
    comp = (B*M_x*kappa*0.1/((d*p)**2.0))*1E26
    s850 = comp*M
    print 'For flux at a single mass', M
    print 'flux at at temp.',T,' equals ',s850
    return

if __name__ == "__main__":
    flux(T,M,l,kappa,d)
