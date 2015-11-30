

import math
import cmath                                    # these import necessary maths programs to deal with complex numbers

#------------------------------------------------------------
#Class - where appropriate (I or II)
#------------------------------------------------------------


#------------------------------------------------------------
#degree calculations
#------------------------------------------------------------

def deg(A,B,C,D,E,F):

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

def time(Ra,Dec):

    #----Ra-------#

    A = Ra - 270

    b = A*4

    i, d = divmod(b, 1)

    B = i
    C = d*60

    #----Dec-------#

    D = Dec + 3

    e = D*60

    i, d = divmod(e, 1)

    E = i
    F = d*60

    print "\n"
    print "right ascention ", str(18),":",str(B),":",str(C)
    print "Declination ", str(-3),":",str(E),":",str(F)

    #hours format
    Ra_h = 18 + (B/60) + (C/3600)
    Dec_h = -3 + (E/60) + (F/3600)

    print Ra_h, Dec_h
    return

if __name__ == "__main__":
    deg(A,B,C,D,E,F)
    time(Ra,Dec)
