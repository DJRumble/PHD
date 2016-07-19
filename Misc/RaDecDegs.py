#Damian Rumble, UoE
#20131212

#This is a module for converting between Ra Dec and degrees. Both ways.

import math
import cmath                                    # these import necessary maths programs to deal with complex numbers

#------------------------------------------------------------
#Class - where appropriate (I or II)
#------------------------------------------------------------


#------------------------------------------------------------
#degree calculations
#------------------------------------------------------------

def deg(A,B,C,D,E,F):
    #Takes hh:mm:ss and hh:mm:ss of Ra and Dec
    a = A*15
    b = B*(1.0/4.0)
    c = C*(1.0/240.0)

    d = D
    e = E*(1.0/60.0)
    f = F*(1.0/3600.0)

    x = a+b+c
    y = d+e+f

    print "\n"
    print "right ascention ", str(A),":",str(B),":",str(C)
    print "Declination ", str(D),":",str(E),":",str(F)

    print "degrees"

    print str(x),  str(y)
    ra = str(x)
    dec =  str(y)
    print ra, dec

    return
    
#------------------------------------------------------------
#time calculations
#------------------------------------------------------------

def timeRa(Ra,hours):
    #Takes Ra  in degrees - returns in hours format
    #----Ra-------#

    A = Ra - 270

    b = A*4

    i, d = divmod(b, 1)

    B = i
    C = d*60

    if hours == 'TRUE':
     #hours format
        Ra_h = 18 + (B/60) + (C/3600)
    else:
        print "\n"
        print "right ascention ", str(18),":",str(B),":",str(round(C,1))
        Ra_h = str(18) +' '+ str(int(B)) +' '+ str(round(C,1))
        
    return Ra_h

def timeDec(Dec,hours):
    #Takes Ra and Dec in degrees - returns in hours format
    #----Dec-------#

    D = Dec + 3

    e = D*60

    i, d = divmod(e, 1)

    E = i
    F = d*60

    if hours == 'TRUE':
    #hours format
        Dec_h = -3 + (E/60) + (F/3600)
    else:
        print "Declination ", str(-3),":",str(E),":",str(round(F,0))
        Dec_h = str(-3) +' '+ str(int(-E)) +' '+ str((int(F)))
    return Dec_h

if __name__ == "__main__":
    deg(A,B,C,D,E,F)
    timeRa(Ra,hours)
    timeDec(Dec,hours)
