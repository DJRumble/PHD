import numpy as np
import os
import string

maps = np.loadtxt('maps4fluxextractor_master.txt',dtype='string')

kapdir = '/star/bin/kappa'


#s850,temp,dtemp,alpha,proto,YSO,clumps,clumpcat

for i in maps:
    #split INPUT file
    files = string.split(i,',')

    print files[0]

    for j in files:
        cmd = '%s/stats ndf=%s QUIET'%(kapdir,j)
        os.system(cmd)

    
