"""Script for creating maps of mass, Column density & number density from input maps of flux & temperature - V2 is specifically for run14 which requires use in a for loop"""

import numpy as np
import os

########### set CSH Commands directories #############

kapdir = '/stardev/bin/kappa'
kapdir2 = '/star/bin/kappa'
convdir = '/stardev/bin/convert'

#####################################################
#Function should take inputs of a the temperautre map and 850 flux map and return maps of mass and Column densityColumn Density - rumble (850 only)

def mass_map(S,T,kappa,d,i,l,body): 
    """input variables; flux map (S in JY .sdf format), temperature map (T in Kelvin .sdf format), oppacity (kappa in 0.012), distance (d in parsecs), i is an iterable, l is the wavelength & body is clumps or cores used (both the latter are only used the resulting name."""
    #Constants
    au = 149597871000 #m
    M_x = 1.989E30 #kg
    m_h2 = 3.34745E-24 #g
    mu = 2.333333  #ratio of H2 to He (5:1)

    #each pixel has an area in SI units and each apature contains N pixel
    if l == 450:
        a = 1.99987
    elif l == 850:
        a = 2.99987
    else:
        a = 1
    A = ((((a*d)*au)**2.0))*10000 #in cm
    
    # makes a single scalar out of all of the other constants in equations
    masscomp =  0.39*((kappa/0.012)**(-1.0))*((d/250.0)**2.0) # - puts things into grammes

    # various constants in the 'temp' dependant part of the mass equation
    one = 1.0
    sun_g = M_x * 1000
    base = 'Natural'
    
    #MATHS expressions
    exp = "'exp(ia/ib)'"

    #temporary maps
    s17 = '17.sdf'
    massmap_a = 'temp/massmap_a.sdf'
    massmap_b = 'temp/massmap_b.sdf'
    massmap_c = 'temp/massmap_c.sdf'
    mass = body+'/mass/mass_'+str(l) + str(body) + str(i)+'.sdf'

    #Kappa process
    cmd = '%s/maths exp=%s out=%s ia=%s ib=%s'%(kapdir,exp,massmap_a,s17,T)
    #print 'IMPORTANT: right now part a does not work, massmap_a has to be made manually in the terminal\n'
    os.system(cmd)
    cmd = '%s/csub in=%s scalar=%s out=%s'%(kapdir,massmap_a,one,massmap_b) #exp(17/T)-1
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,massmap_b,masscomp,massmap_c)#(exp(17/T)-1)*0.19*((kappa/0.01)(-1.0))
    os.system(cmd)
    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,S,massmap_c,mass)#(exp(17/T)-1)*0.19*((kappa/0.01)**(-1.0))*S
    os.system(cmd)

    print 'Mass calculation complete'
    return mass


def CD_map(S,T,kappa,d,i,l,body):
    """input variables; flux map (S in JY .sdf format), temperature map (T in Kelvin .sdf format), oppacity (kappa in 0.012), distance (d in parsecs), i is an iterable, l is the wavelength & body is clumps or cores used (both the latter are only used the resulting name."""
    #Constants
    au = 149597871000 #m
    M_x = 1.989E30 #kg
    m_h2 = 3.34745E-24 #g
    mu = 2.333333  #ratio of H2 to He (5:1)

    #each pixel has an area in SI units and each apature contains N pixels 
    if l == 450:
        a = 3.99987
    elif l == 850:
        a = 5.99987
    else:
        a = 1

    A = ((((a*d)*au)**2.0))*10000 #in cm


    #import mass 
    mass = mass_map(S,T,kappa,d,i,l,body)


    #Density components - makes a single scalar out of all of the other constants in equations
    ColumnDcomp = 1.0/(mu*A*m_h2)
    Densitycomp = 1.0/A
    sun_g = M_x * 1000

    #temporary files
    massmap_d = 'temp/massmap_d.sdf'

    #Named files
    ColumnD = body+'/mass/SerpensMWC297_20131025_CD'+str(l) + str(body) + str(i)+'.sdf'
    density = body+'/mass/SerpensMWC297_20131025_D'+str(l) + str(body) + str(i)+'.sdf'

    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,mass,sun_g,massmap_d) #puts mass solar masses map into g
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,massmap_d,ColumnDcomp,ColumnD) #mass*(1/A*(100^2))*(1.0/m_h2) 
    os.system(cmd)   
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,massmap_d,Densitycomp,density) #mass*(1/A*(100^2)) 
    os.system(cmd)

    print 'Density calculations complete'
    return 


if __name__ == "__main__":
    mass_map(S,T,kappa,d,i,l,body)
    CD_map(S,T,kappa,d,i,l,body)
