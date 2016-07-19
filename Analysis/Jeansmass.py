"""Script for caculating mass as calculated by Jeans Mass"""

#####################################################
#Calculates the jeans mass of a clump. require list of average temps and total pixels

def jeans(T,d,N): 
    """Takes inputs of  temperature (T in Kelvin), distance (d in pcs) and Number of pixels in area (N)"""
    #constants
    au = 149597871000 #m
    M_x = 1.989E30 #kg
    pi = 3.14159265359
    G = 6.67384E-11
    m_h2 = 3.34745E-24 #g
    mu = 2.333333  #ratio of H2 to He (5:1)
    k = 1.3806488E-23
    
    #each pixel has an area in SI units and each apature contains N pixels
    a = 5.99987
    print 'n = ',N
    #total area of the clump.
    A = ((((a*d)*au)**2.0))*N 
    
    r = (A/pi)**(0.5)
    lambdaj = 2.0*r

    print 'jeans lengths are',lambdaj/au,' in au'
    # calculating jeans mass using mean Temp. for the whole clump... assuming spherical sysmetry
    Mj = (((((pi**2.0)*k)/(6.0*G*mu*m_h2))*lambdaj)*1000/M_x)*T
    print 'jeans masses are',Mj,' in solar masses'
    return 

if __name__ == "__main__":
    jeans(T,d,N)
