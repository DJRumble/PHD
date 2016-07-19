""" this script runs the Kolmogorov-Smirnov test for one and two data sets"""

import numpy as np

#####################################################
#TWO

#This function calculates the KS-test for two unique distributions. 

def two(data1,data2):
    """Inputs: Data set 1, data set 2. Returns K-S statistic 'D' and the significance level 'prob' for the null hypothesis (i.e. the two distributions are identical. The smaller the valie of 'prob', the greater the divergence between data sets"""
    
    n1 = len(data1)
    n2 = len(data2)

    i1 = 1
    i2 = 1
    Dnull = 0.0

    #Sort data into ascending order
    Data1 = data1.sort() 
    Data2 = data2.sort()

    while ((i1 <= n1) && (i2 <= n2)):
        if ((Data1[i1]) <= (Data2[i2])):
            fn1 = (i1+1)/n1
        if ((Data2[i2]) <= (Data1[i1])):
            fn2 = (i2+1)/n2
        if ((abs(fn2-fn1))> Dnull):
            Dnull = abs(fn2-fn1)
    
    prob = proKS(np.sqrt(n1*n2/(n1+n2))*Dnull)

    return Dnull, prob

def one(data1,func):
    """Inputs: Data set and a specific function. Returns K-S statistic 'D' and the significance level 'prob' for the null hypothesis (i.e. the two distributions are identical. The smaller the valie of 'prob', the greater the divergence between the data and the model"""

    return

def probKS(alam):
    
    EPS1 = 0.001
    EPS2 = 1.0E-8

    a2 = -2.0*(alam**2.)
    
    termbf = 0.0
    sum = 0.0
    fac = 2.0
    
    j = 1

    for (j <=100):
        term = fac*np.exp(a2*(j**2.))
        sum = sum + term
        if (abs(term) <= EPS1*termbf) || (abs(term) <= EPS2.sum):
            print 'converged'
        fac = -fac
        termbf=abs(term)

    return sum

if __name__ == "__main__":
    two(data1,data2)
