#!/usr/bin/python
# slfit.py
# weighted straight line fit adapted from Numerical Recipes by J. Hatchell
# (originally in bayes.py)

import numpy as np

def slfit(x,y,sigma):
    """ Least-chisquared straight line fit from Numerical Recipes Chapter 14
    x, y arrays of data abcissa and value
    sigma = assumed Gaussian uncertainties on y
    Returns:
    m gradient
    c y-intercept
    sigm uncertainty on gradient
    sigc uncertainty on y-intercept
    cov covariance between m and c"""
    #print "x = ", x
    #print "y = ", y
    #print "sigma = ", sigma
    # Calculate weights: fix large where values of sigma are 0
    W = np.where(sigma>0,[1./i for i in sigma] ,100*max(sigma.max(),1e-10))
    S = np.sum(W*W)
    Sx = np.sum(x*W*W)
    Sy = np.sum(y*W*W)
    # Weighted sum NR 14.2.15
    t = W*(x-Sx/S)
    # Weighted sum NR 14.2.16
    Stt = np.sum(t*t)
    #Sy = np.sum(y*W*W)
    #Sxx = np.sum(x*x*W*W)
    #Syy = np.sum(y*y*W*W)
    #Sxy = np.sum(x*y*W*W)
    #print "S*Sxx =", S*Sxx
    #print "Sx*Sx = ", Sx*Sx
    #delta = (S*Sxx-Sx*Sx)
    #c = (Sxx*Sy-Sx*Sxy)/delta
    #m = (S*Sxy-Sx*Sy)/delta
    m = (1/Stt)*np.sum(t*y*W)
    c = (Sy - Sx*m)/S
    sigc = np.sqrt(1/S*(1+Sx*Sx/S/Stt))
    sigm = np.sqrt(1/Stt)

    cov = -Sx/S/Stt

    return(m,c,sigm,sigc,cov)
    
if __name__ == "__main__":
    slfit(x,y,sigma)
