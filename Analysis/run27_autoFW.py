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


for i in maps:
    map,mask = string.split(i,',')
    prefix = string.split(map,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    region = string.split(file,'850')[0]
    
    '''
    #PART 1 #### FIND BACK 
    box = 80
    FB = 'FW/temp/'+file+'FB.sdf'
    print 'Findback box size ',box
    cmd = '%s/findback in=%s out=%s box=%s sub=TRUE %s'%(cupdir,map,FB,box,'> /dev/null')
    os.system(cmd)
    
    ### mask 850 maps to remove fringes
    MSK = 'FW/temp/'+file+'FB_MSK.sdf'
 
    mult = '%s/mult in1=%s in2=%s out=%s'%(kapdir,FB,mask,MSK)
    os.system(mult)

    print 'MANUAL SECTION: '+str(MSK)+' NOW'
    '''
#Now manually section FB maps, adding 'section' to the end of each file name in 'FW/temp/'

    #PART 2 #### FINDclumps
    print '==============================================='
    print 'Finding clumps in: '+file
    #Base 1sigma level
    sigma = noise.noise_by_data(map,'FALSE')
    print 'sigma noise: '+str(sigma)

    inmap = 'FW/temp/'+file+'FB_MSKsection.sdf'

    clumps = 'FW/'+region+'clumps.sdf'
    cat = 'FW/cat/'+region+'clump_cat.fit'


    config = '^inparams.txt'

    cmd = '%s/findclumps in=%s out=%s outcat=%s method=FellWalker deconv=True config=%s RMS=%s '%(cupdir,inmap,clumps,cat,config,sigma)
    os.system(cmd)
    


