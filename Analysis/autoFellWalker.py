#Damian Rumble, UoE
#20151020
#autoFellwalker

#A script for producing FellWalker clump maps from multiple inpouts

#####################################################
#import maths, ploting and astrophysical packages
import os
import numpy as np
import string
import noise

########### set CSH Commands directories #############

cupdir = '/star/bin/cupid'
kapdir = '/star/bin/kappa'

########### Bulk Code #############'

maps = np.loadtxt('mapsFW.txt',dtype='string')

j = 0 
for i in maps:
    prefix = string.split(i,'.sdf')[0]
    #print prefix
    #if j <= 6:
    #    J = 7
    #else:
    #    J = 6

    J = 1
    file = string.split(prefix,'/')[J]
    region = string.split(file,'850')[0]
    print region
    
    #PART 1 #### SNR + FIND BACK 
    print '==============================================='
    print 'Finding SNR + FB in: '+file

    snr = 'FW/temp/'+file+'SNR.sdf'
    cmd = '%s/makesnr in=%s out=%s minvar=0'%(kapdir,i,snr)
    os.system(cmd)

    #box = 120
    #print box
    #snrFB = 'FW/temp/'+file+'snrFB.sdf'
    #cmd = '%s/findback in=%s out=%s box=%s sub=TRUE %s'%(cupdir,snr,snrFB,box,'> /dev/null')
    #os.system(cmd)

    #PART 2 #### FIND clumps
    print '==============================================='
    print 'Finding clumps in: '+snr
 
    sigma = noise.noise_by_snr(snr,'FALSE')

    clumps = 'FW/clumps/'+region+'clumps.sdf'
    cat = 'FW/cat/'+region+'clump_cat.fit'

    config = '^inparams.txt'

    cmd = '%s/findclumps in=%s out=%s outcat=%s method=FellWalker deconv=True config=%s RMS=%s '%(cupdir,snr,clumps,cat,config,sigma)
    os.system(cmd)
    
    print 'CLUMP MAP:'+clumps

    j = j + 1
    #break

tempfolder = 'FW/temp/*'
cmd = 'rm %s'%(tempfolder)
os.system(cmd)
