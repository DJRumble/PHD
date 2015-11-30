#Damian Rumble, UoE
#20140117
#Run16

#A script for producing FellWalker clump maps

#####################################################
#import maths, ploting and astrophysical packages

import os

########### set CSH Commands directories #############

cupdir = '/star/bin/cupid'
kapdir = '/star/bin/kappa'

########### Bulk Code #############'

#s450 = 'SerpensMWC297_20140325_IR1_s450_freefree.sdf'
s850 = 'SerpensMWC297_20140325_IR1_s850_freefree.sdf'
temp = 'SerpensMWC297_20140325_IR1_freefree450850temp%5.sdf'

area = '(~380,~150)'

section = "'"+str(s850)+str(area)+"'"
print section

inmap = 'SerpensMWC297_IR1_section.sdf'

cmd = '%s/ndfcopy in=%s out=%s'%(kapdir,section,inmap)
os.system(cmd)

box = 20

FB = 'SerpensMWC297_IR1_box'+str(box)+'.sdf'

cmd = '%s/findback in=%s out=%s box=%s sub=TRUE %s'%(cupdir,inmap,FB,box,'> /dev/null')
os.system(cmd)

#Base 1sigma level
sigma = 0.0016814645435

#Other parameters
minpix = 4
maxbad = 100
FWHMbeam = 0
Mindip = 1
Maxjump = 1
minhieght = 3
noise = 3
allowedge = 1

"""
params = open("inparams.txt","w")
params = open("inparams.txt","a")

params.write('FELLWALKER.ALLOWEDGE='+str(allowedge)+'\n')
params.write('FELLWALKER.FWHMBEAM='+str(FWHMbeam)+'\n')
params.write('FELLWALKER.CLEANITER=1\n')
params.write('FELLWALKER.FLATSLOPE=1*RMS\n')
params.write('FELLWALKER.MAXJUMP='+str(Maxjump)+'\n')
params.write('FELLWALKER.MINDIP='+str(Mindip)+'*RMS\n')
params.write('FELLWALKER.MINHEIGHT='+str(minhieght)+'*RMS\n')
params.write('FELLWALKER.MINPIX='+str(minpix)+'\n')
params.write('FELLWALKER.MAXBAD='+str(maxbad)+'\n')
params.write('FELLWALKER.NOISE='+str(noise)+'*RMS\n')
"""

clumps = 'SerpensMWC297_IR1_clumps'+str(noise)+'.sdf'
cat = 'SerpensMWC297_IR1_clumps'+str(noise)+'_params.fit'
config = '^inparams.txt'

cmd = '%s/findclumps in=%s out=%s outcat=%s method=FellWalker config=%s RMS=%s '%(cupdir,FB,clumps,cat,config,sigma)
os.system(cmd)
