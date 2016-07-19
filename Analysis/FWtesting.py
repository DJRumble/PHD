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

map = '/data/damian/maps/aquila/serpens_south/IR2/W40_noco_extS2nosm_s850-4am-fullfreefree.sdf'
mask = 'SNRmask/Aquila_noco-autoSNRmask.sdf'
#map = '/data/damian/maps/lupus/IR2/Lupus_20150318_850_IR2_ext_HKJypix_col.sdf'
#mask = '/data/damian/temp_MKR_kernel/output/kernel/map/Lupus_20150318-auto/Lupus_20150318-autoSNRmask.sdf'

prefix = string.split(map,'.sdf')[0]
file = string.split(prefix,'/')[7]
region = string.split(file,'_')[0]

config = '^inparams.txt'



#TEST 1 #### FIND BACK 
box = 80    
FB = 'FW/temp/'+file+'FB.sdf'
print '1 FIND BACK '
cmd = '%s/findback in=%s out=%s box=%s sub=TRUE %s'%(cupdir,map,FB,box,'> /dev/null')
#os.system(cmd)
MSK = 'FW/temp/'+file+'FB_MSK.sdf'
mult = '%s/mult in1=%s in2=%s out=%s'%(kapdir,FB,mask,MSK)
os.system(mult)
clumps1 = 'FW/'+region+'clumpsFB.sdf'
cat1 = 'FW/cat/'+region+'clump_catFB.fit'
sigma1 = noise.noise_by_data(FB,'FALSE')
print clumps1
cmd = '%s/findclumps in=%s out=%s outcat=%s method=FellWalker deconv=True config=%s RMS=%s '%(cupdir,MSK,clumps1,cat1,config,sigma1)
os.system(cmd)

#TEST 2 #### FINDBACK + MAKESNR
print '2 FINDBACK + MAKESNR'
FBsnr = 'FW/temp/'+file+'FBSNR.sdf'
cmd = '%s/makesnr in=%s out=%s'%(kapdir,FB,FBsnr)
os.system(cmd)
#snrMSK = 'FW/temp/'+file+'SNR_MSK.sdf'
#mult = '%s/mult in1=%s in2=%s out=%s'%(kapdir,FBsnr,mask,snrMSK)
#os.system(mult)
clumps2 = 'FW/'+region+'clumpsSNR.sdf'
cat2 = 'FW/cat/'+region+'clump_catSNR.fit'
sigma2 = noise.noise_by_snr(FBsnr,'FALSE')
cmd = '%s/findclumps in=%s out=%s outcat=%s method=FellWalker deconv=True config=%s RMS=%s '%(cupdir,FBsnr,clumps2,cat2,config,sigma2)
os.system(cmd)
print clumps2

#TEST 3 #### MAKESNR + FINDBACK
print '3 MAKESNR + FINDBACK'
snr = 'FW/temp/'+file+'SNR.sdf'
cmd = '%s/makesnr in=%s out=%s'%(kapdir,map,snr)
#os.system(cmd)
snrFB = 'FW/temp/'+file+'SNRFB.sdf'
print 'Findback box size ',box
cmd = '%s/findback in=%s out=%s box=%s sub=TRUE %s'%(cupdir,snr,snrFB,box,'> /dev/null')
#os.system(cmd)
#snrFBMSK = 'FW/temp/'+file+'SNRFB_MSK.sdf'
#mult = '%s/mult in1=%s in2=%s out=%s'%(kapdir,snrFB,mask,snrFBMSK)
#os.system(mult)
clumps3 = 'FW/clumps/'+region+'clumpsFBSNR.sdf'
cat3 = 'FW/cat/'+region+'clump_catFBSNR.fit'
sigma3 = noise.noise_by_snr(snrFB,'FALSE')
print clumps3
cmd = '%s/findclumps in=%s out=%s outcat=%s method=FellWalker deconv=True config=%s RMS=%s '%(cupdir,snrFB,clumps3,cat3,config,sigma3)
os.system(cmd)
