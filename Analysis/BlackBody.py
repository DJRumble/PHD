"""Script for calculating the BlackBody function"""

import numpy as np

#####################################################
#Function to deffine Black Body Constant

def BB(l,T):
    """Takes inputs of; wavelength (l in m) & Temperature (T in Kelvin)"""
    #constants
    h = 6.626068E-34
    c = 2.99792458E8
    k = 1.3806488E-23

    #calculate frequency
    nu = c / l
    #equation of a blackbody function
    B_a = ((2.0*h*(nu**3.0))/(c**2.0))
    B_b = (1.0/(np.exp((h*nu)/(k*T))-1.0))
    B = B_a*B_b
    print 'Blackbody function = ', B
    return B

if __name__ == "__main__":
    BB(l,T)
