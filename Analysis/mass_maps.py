"""Script for creating maps of mass, Column density & number density from input maps of flux & temperature"""

import numpy as np
import os

########### set CSH Commands directories #############

kapdir = '/stardev/bin/kappa'
kapdir2 = '/star/bin/kappa'
convdir = '/stardev/bin/convert'

#####################################################
#Function should take inputs of a the temperautre map and 850 flux map and return maps of mass and Column densityColumn Density - rumble (850 only)

def mass_map(S,T,d,l): 
    """input variables; flux map (S in JY .sdf format), temperature map (T in Kelvin .sdf format), oppacity (kappa in 0.012), distance (d in parsecs) - 850 wavelengths only, l is string wavelength"""
    #Constants
    au = 149597871000 #m
    pc = 3.08567758E16 #m
    M_x = 1.989E30 #kg
    m_h2 = 3.34745E-24 #g
    mu = 2.333333  #ratio of H2 to He (5:1)

    #each pixel has an area in SI units and each apature contains N pixels
    if l == '450':
        a = 3.99987 #450
        kappa = 0.048
    elif l == '850':
        a = 5.99987 #850
        kappa = 0.012
    else:
        a = 1.

    print a
    A = ((((a*d)*au)**2.0))*10000 #in cm
    V = A*(pc*250*100) #Volume in cm 
    # makes a single scalar out of all of the other constants in equations
    masscomp =  0.39*((kappa/0.012)**(-1.0))*((d/250.0)**2.0) # - puts things into grammes
    ColumnDcomp = 1.0/(mu*V*m_h2)
    Densitycomp = 1.0/A

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
    massmap_d = 'temp/massmap_d.sdf'
    mass = 'SerpensMWC297_s'+l+'_20140115_M.sdf'
    ColumnD = 'SerpensMWC297_s'+l+'_20140115_CD.sdf'
    density = 'SerpensMWC297_s'+l+'_20140115_D.sdf'

    #Kappa process
    cmd = '%s/maths exp=%s out=%s ia=%s ib=%s'%(kapdir,exp,massmap_a,s17,T)
    os.system(cmd)
    cmd = '%s/csub in=%s scalar=%s out=%s'%(kapdir,massmap_a,one,massmap_b) #exp(17/T)-1
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,massmap_b,masscomp,massmap_c)#(exp(17/T)-1)*0.19*((kappa/0.01)(-1.0))
    os.system(cmd)
    cmd = '%s/mult in1=%s in2=%s out=%s'%(kapdir,S,massmap_c,mass)#(exp(17/T)-1)*0.19*((kappa/0.01)**(-1.0))*S
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,mass,sun_g,massmap_d) #puts mass solar masses map into g
    os.system(cmd)
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,massmap_d,ColumnDcomp,ColumnD) #mass*(1/A*(100^2))*(1.0/m_h2) 
    os.system(cmd)   
    cmd = '%s/cmult in=%s scalar=%s out=%s'%(kapdir,massmap_d,Densitycomp,density) #mass*(1/A*(100^2)) 
    os.system(cmd)
    print 'Mass calculation complete'
    return

if __name__ == "__main__":
    mass_map(S,T,d,l)
