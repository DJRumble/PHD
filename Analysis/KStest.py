""" this script runs the Kolmogorov-Smirnov test for one and two data sets"""

import numpy as np

#####################################################
#TWO

#This function calculates the KS-test for two unique distributions. 

def two(data1,data2):
    """Inputs: Data set 1, data set 2. Returns K-S statistic 'D' and the significance level 'prob' for the null hypothesis (i.e. the two distributions are identical. The smaller the valie of 'prob', the greater the divergence between data sets"""
    
    n1 = len(data1)
    n2 = len(data2)

    i1 = 0
    i2 = 0
    Dnull = 0.0
    fn1 = 0.0
    fn2 = 0.0

    #Sort data into ascending order
    data1.sort() 
    data2.sort()

    print n1,data1
    print n2,data2

    #my version
    D = []

    for i1 in range (n1):
        Dnull = abs(data1[i1]-data2[i2])
        i1 = i1 + 1
        i2 = i2 + 1
        D.append(Dnull)

    print D

    D.sort(reverse=1)
    
    print D

    #numerical recipies version
    #while (i1 <= n1) & (i2 <= n2):
    #    print i1, i2
    #    if ((data1[i1]) <= (data2[i2])):
    #        fn1 = (float(i1))/float(n1)
    #        i1 = i1 + 1
    #        #print 'fn1',fn1
    #    if ((data2[i2]) <= (data1[i1])):
    #        fn2 = (float(i2))/float(n2)
    #        i2 = i2 + 1
     #       #print 'fn2',fn2
    #    if ((abs(fn2-fn1))> Dnull):
    #        Dnull = abs(fn2-fn1)
    #        print Dnull

    D = np.sqrt(n1*n2/(n1+n2))*D[0]

    print D

    prob = prob

KS(D,data1)

    print prob

    return Dnull, prob

def one(data1,func):
    """Inputs: Data set and a specific function. Returns K-S statistic 'D' and the significance level 'prob' for the null hypothesis (i.e. the two distributions are identical. The smaller the valie of 'prob', the greater the divergence between the data and the model"""

    return

def probKS(alam,data1):
    
    N = len(data1)

    i = 0
    Q = 0.0
    
    for i in range(N):
        a2 = -2.0*(alam**2.)*(i**2.0)
        q = ((-1.)**(i-1))*np.exp(a2)
        Q = Q + q

    prob = 1-(2*Q)
    
    return prob
    
    #termbf = 0.0
    #sum = 0.0
    #fac = 2.0
    
    #j = 1

    #for j in range(0,100):
    #    term = fac*np.exp(a2*(j**2.))
    #    sum = sum + term
        #if (abs(term) <= EPS1*termbf) || (abs(term) <= EPS2.sum):
        #fac = -fac
        #termbf=abs(term)
    #return sum

if __name__ == "__main__":
    two(data1,data2)
