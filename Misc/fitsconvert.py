import ndf2fits
import numpy as np 

data = np.loadtxt('maps.txt',dtype='string')

for i in data:
    ndf2fits.ndf2fits(i)
