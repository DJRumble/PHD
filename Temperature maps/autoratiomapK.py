#20150924 - autoratiomapK.py
#DJR, UoE
#A wrapper for the automatic production of temp. maps from a list, automaplist.txt.

import numpy as np
import string
import sys
import os
import noise

#Load files- specific order 450, then 850.
files = np.loadtxt('automaplist.txt',dtype='string',comments="#")

#NOI = open("Noise.tab","w")
#NOI = open("Noise.tab","a")
#NOI.write('#Region\ts450\ts850\n')

for i in files:
    input = string.split(i,',')
    A = string.split(input[2],'_')[0]
    B = string.split(input[2],'_')[1]
    outdir = A+'_'+B+'-auto'
    print 'Making temperature map for ',outdir
    
    cmd = 'python ratiomap450K850-AUTO.py %s %s %s %s K'%(input[0],outdir,input[1],input[2])
    os.system(cmd)

    #TEST NOISE
    
    #in450 = '/data/damian/maps/'+str(input[0])+str(input[1])+'.sdf'
    #in850 = '/data/damian/maps/'+str(input[0])+str(input[2])+'.sdf'

    #n450 = 5*noise.noise_by_data(in450,'false')
    #n850 = 5*noise.noise_by_data(in850,'false')

    #print 'TEST NOISE: ',outdir
    #print 'NOISE 450: ',n450, ' NOISE 850: ',n850

    #NOI.write(outdir+'\t'+str(round(n450,4))+'\t'+str(round(n850,4))+'\n')

    break




