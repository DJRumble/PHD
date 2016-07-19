"""Script for caculating Jeans mass and length for a clump containing N pixels - assuming spherical symmetry"""

#####################################################
#Calculates the jeans mass of a clump. require list of average temps and total pixels

def jeans(T,d,N): 
    """Takes inputs of  temperature (T in Kelvin), distance (d in pcs) and Number of pixels in area (N)"""
    #constants
    au = 149597871000 #m
    pc = 3.08567758E+16
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
    R = (N**(0.5))*5.06E+11*d #m
    r = (((((a*d)*au)**2.0))*N/pi)**(0.5) #m

    lambdaj = 2.0*r

    print 'jeans lengths are',lambdaj/au,' in au'
    # calculating jeans mass using mean Temp. for the whole clump... assuming spherical sysmetry
    Mj = (((((pi**2.0)*k)/(6.0*G*mu*m_h2))*lambdaj)*1000/M_x)*T
    print 'jeans masses are',Mj,' in solar masses'
    return Mj

def sadavoy(T,d,N):
    """Takes inputs of  temperature (T in Kelvin), distance (d in pcs) and Number of pixels in area (N)"""
    psc = 3.08567758E+16 #m
    M_x = 1.989E30 #kg
    au = 149597871000 #m
    pi = 3.14159265359

    #each pixel has an area in SI units and each apature contains N pixels
    a = 5.99987
    #print 'n = ',N
    #total area of the clump.
    r = (N**(0.5))*5.06E+11*d #m

    lambdaj = 2.0*r

    #Uses SADAVOY's method of calculating jeans mass 
    Mj = 1.9*(T/10)*((r/psc)/0.07)
    #print 'jeans masses are',Mj,' in solar masses'
    return Mj

if __name__ == "__main__":
    jeans(T,d,N)
    sadavoy(T,d,N)
