"""Script for calculating the Mass using Kirk's method at 850microns at 250pc at kappa = 0.02"""

import numpy as np

#####################################################
#Function to deffine Mass - Kirk06 (850 only)
def kirk(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""

    #Equation of Mass in solar masses, per pixel
    M = 0.23*S*(np.exp(17.0/T)-1.0)*((kappa/0.02)**(-1.0))*((d/250.0)**2.0)
    print 'mass of apature in solar masses:', M
    return

if __name__ == "__main__":
    kirk(S,T,kappa,d)
