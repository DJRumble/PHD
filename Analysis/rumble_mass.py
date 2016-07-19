"""Script for calculating the Mass using Rumble's method at 850microns at 250pc at kappa = 0.012"""

import numpy as np

#####################################################
#Function to deffine Mass - Kirk06 (850 only)
def djr(S,T,kappa,d):
    """Takes inputs of; Flux (S in Jk), Temperature (T in Kelvin), oppacity (kappa in g cm-2), distance (d in pcs)"""

    #Equation of Mass in solar masses, per pixel
    M = 0.39*S*(np.exp(17.0/T)-1.0)*((kappa/0.012)**(-1.0))*((d/250.0)**2.0)
    print 'mass of apature in solar masses:', M
    return M

if __name__ == "__main__":
    djr(S,T,kappa,d)
