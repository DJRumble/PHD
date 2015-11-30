#Damian Rumble, UoE
#20131212

#This code runs Autophotometery for CORES, preselected by FellWalker, in MCW297 region

#####################################################
#import modules

import coordconverssion
import mass_maps2
import RaDecDegs

import os
import commands
import pyfits
import numpy as np


########### set CSH Commands directories #############

photomdir = '/star/bin/photom'
kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

#####################################################
#contants

kappa_8 = 0.012 #opacity
kappa_4 = 0.048 #opacity
d = 250 #distance (pcs)
errd = 50 #distance error
#each pixel has an area in SI units and each apature contains N pixels
a = 5.99987
au = 149597871000 #m
pi = 3.14159265359
psc = 3.08567758E+16 #m

########### Functions #############

def flux_450(mask,s450,c450):
    mul1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,s450,c450)
    os.system(mul1)
    return c450
    
def flux_850(mask,s850,c850):
    mul2 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,s850,c850)
    os.system(mul2)
    return c850

def temp_45(mask,temp,ctemp):
    T = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,temp,ctemp)
    os.system(T)
    return ctemp

def mass_err(S,T,dS): #Mass error
    S = float(S)
    print S
    print dS
    print T
    A = (0.39**2.)*(((np.exp(17./T)-1)**2.)*(dS**2.))
    B = (0.39**2.)*(289.*(S**2.)/(T**4.))*np.exp(34./T)*((0.05*T)**2.)
    
    M = A + B
    dM = np.sqrt(M)

    return dM


########### Bulk Code #############'

#set parameters for coord conversion

s450 = 'SerpensMWC297_20140325_IR1_s450_freefree.sdf'
s850 = 'SerpensMWC297_20140325_IR1_s850_freefree.sdf'
temp = 'SerpensMWC297_20140325_IR1_freefree450850temp%5.sdf'     

params = 'SerpensMWC297_IR1_clumps3_params.fit'

data = pyfits.getdata(params)

#produce the X Y positions in grid coordinates
XYarray_850 = coordconverssion.coord(s850,params)
XYarray_450 = coordconverssion.coord(s450,params)

#write AUTO PHOTOM input paramete file. Format is:
#         1  86.250000 -58.980000     0    0            0           0          OK   5.000000   0.000000   0.000000   regions     circle

inparam4 = open("inparams_450_photom.txt","w")
inparam4 = open("inparams_450_photom.txt","a")

inparam8 = open("inparams_850_photom.txt","w")
inparam8 = open("inparams_850_photom.txt","a")

i = 0

for i in range(len(XYarray_450)):
    #write parameter file
    #450
    inparam4.write('\t'+str(i+1)+'\t')
    inparam4.write(str(XYarray_450[i][0])+'\t')
    inparam4.write(str(XYarray_450[i][1])+'\t')
    inparam4.write('0\t0\t0\t0\tOK\t5.0\t0\t0\tregions\tcircle\n')
    #850
    inparam8.write('\t'+str(i+1)+'\t')
    inparam8.write(str(XYarray_850[i][0])+'\t')
    inparam8.write(str(XYarray_850[i][1])+'\t')
    inparam8.write('0\t0\t0\t0\tOK\t5.0\t0\t0\tregions\tcircle\n')    

    i = i + 1

#### this part doesn't work in the code but does work when run seperately in python shell. As I only need to use it once I'm going to let this slide ###
'''
out='photom_850_cores.txt'

sig_850 = 0.00218242691874
cmd = '%s/autophotom IN=%s INFILE=%s OUTFILE=%s CENTRO=FALSE photon=4 sky=0 skyest=4 skysig=%s usemags=false PADU=1  '%(photomdir, 'SerpensMWC297_20130918_232648_s850_DJR_cal.sdf', 'inparams_850_photom.txt', 'photom_850cores.txt',sig_850)
os.system(cmd)

sig_450 = 0.0168632057269
out='photom_450_cores.txt'
cmd = '%s/autophotom IN=%s INFILE=%s OUTFILE=%s CENTRO=FALSE photon=4 sky=0 skyest=4 skysig=%s usemags=false PADU=1  '%(photomdir, 'SerpensMWC297_20130916_233448_s450_DJR_cal.sdf', 'inparams_450_photom.txt', 'photom_450cores.txt',sig_450)
os.system(cmd)
'''
### end ###

i = 0

core_LTX = open("core_LTX.txt","w")
core_LTX = open("core_LTX.txt","a")

core_LTX.write('index\t&Flux 450\t&\t&Flux 850\t&Mass 850\t&Mean Temp.\t\n')

for i in range(len(XYarray_450)):

    #label Ra and Dec at row
    pX_450 =-1* float(XYarray_450[i][0])
    pY_450 =-1* float(XYarray_450[i][1])
    pX_850 =-1* float(XYarray_850[i][0])
    pY_850 =-1* float(XYarray_850[i][1])

    i = i + 1

    print 'Starting processing of core#'+str(i)
    #create a flux mask for each core postion. This is done by making a simple quadratic plot and then threshing at a given clump size. Clump size is 0.07pcs which varys appature size depending on pixe size

   #450
    if pX_450 > 0:
        pX_450 = '+'+str(pX_450)
    else:
        pX_450 = str(pX_450)

    if pY_450  > 0:
        pY_450 = '+'+str(pY_450)
    else:
        pY_450 = str(pY_450)

    #850
    if pX_850 > 0:
        pX_850 = '+'+str(pX_850)
    else:
        pX_850 = str(pX_850)

    if pY_850  > 0:
        pY_850 = '+'+str(pY_850)
    else:
        pY_850 = str(pY_850)
    
    print pX_450, pY_450, pX_850, pY_850

### FILE SETUP ####

    #set up core files
    c450 = 'cores/s450_core'+str(i)+'.sdf'
    c850 = 'cores/s850_core'+str(i)+'.sdf'
    ctemp = 'cores/temp_core'+str(i)+'.sdf'

#Making mask for all points - fits files
    q450 = 'temp/450/quad'+str(i)+'.sdf'
    q850 = 'temp/850/quad'+str(i)+'.sdf'
    m450 = 'temp/450/mask'+str(i)+'.sdf'
    mq450 = 'temp/450/maskq'+str(i)+'.sdf'
    m850 = 'temp/850/mask'+str(i)+'.sdf'
    masks = 'temp/850/masks/mask'+str(i)+'.fits'
    magic = 'temp/850/magic'+str(i)+'.sdf'

#set up temp file  names
    o450 = 'temp/450/s450_core'+str(i)+'.sdf'
    o850 = 'temp/850/s850_core'+str(i)+'.sdf'
    otemp = 'temp/850/temp_core'+str(i)+'.sdf'
    mtemp = 'temp/meantemp'+str(i)+'.sdf'

### END FILE SETUP ####

    #MATHS quadratic expression, of set to grid coord.
    exp_450 = "'(((xa+0.5)'"+pX_450+"')**2)+(((xb+0.5)'"+pY_450+"')**2)'"
    exp_850 = "'(((xa+0.5)'"+pX_850+"')**2)+(((xb+0.5)'"+pY_850+"')**2)'"

    if  i == 24 or i == 26 or i == 29 :
        print 'Not a real object'
    elif i == 5 or i == 6 or i == 8 or i == 9 or i == 10 or i == 13 or i == 14 or i == 15 or i == 16 or i == 18 or i == 20 or i == 21 or i == 22 or i == 23  or i == 25 or i == 27 or i == 28 :
        print 'Identified as a Clump'
    else:

        cmd = "%s/maths exp=%s out=%s like=%s"%(kapdir,exp_450,q450,s450)
        os.system(cmd)
        cmd = "%s/maths exp=%s out=%s like=%s"%(kapdir,exp_850,q850,s850)
        os.system(cmd)
    
    #THRESH Quadratic maps to creating masks
    #Set aperture levels for cores of 0.05pcs (40'' diameter)
        ap_450 = 26
        ap_850 = 11.2
        
    #make aperture masks and convert to .FITS
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=1 newhi=bad > /dev/null'%(kapdir,q450,mq450,ap_450,ap_450)
        os.system(cmd)
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=1 newhi=bad > /dev/null'%(kapdir,q850,m850,ap_850,ap_850)
        os.system(cmd)

        #Eliminate residual pixels in the mask - see 26/02/2014
        cmd = '%s/substitute in=%s out=%s oldval=26 newval=1 > /dev/null'%(kapdir, mq450, m450)
        os.system(cmd)
    
        cmd = '%s/nomagic in=%s out=%s repval=0  > /dev/null'%(kapdir,m850,magic)
        os.system(cmd)
        cmd = '%s/ndf2fits in=%s out=%s > /dev/null'%(convdir,magic,masks)
        os.system(cmd)

    #Make core maps in Flux 450, 850 & temp. use MULT
    #masking orginal maps, new saved to temp files
        flux_450(m450,s450,c450)
        flux_850(m850,s850,c850)
        temp_45(m850,temp,ctemp)

    #run stats and extract TOTAL value of pixels for 850 and append to files
        cmd = '%s/stats ndf=%s %s'%(kapdir,c450,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c450_i = output
        c450i = float(c450_i)
        print 'full 450 flux = '+str(c450_i)

        cmd = '%s/stats ndf=%s %s'%(kapdir,c850,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c850_i = output
        c850i = float(c850_i)
        print 'full 850 flux = '+str(c850_i)

    #Correction - specific to SMM5/SMM6 due to over lapping aperture. Method: Halve aperture size - scale down full SMM5/SMM6 with respect to other sources. 
        if i == 4:
            print 'Scaling'
            scale_450 = 0.761
            scale_850 = 0.947
            c450i = scale_450*float(c450_i)
            c850i = scale_850*float(c850_i)
            print 'corrected 450 flux = '+str(c450i)
            print 'corrected 850 flux = '+str(c850i)
        elif i == 12:
            print 'Scaling'
            scale_450 = 0.744
            scale_850 = 0.917
            c450i = scale_450*float(c450_i)
            c850i = scale_850*float(c850_i)
            print 'corrected 450 flux = '+str(c450i)
            print 'corrected 850 flux = '+str(c850i)

    #replace 'bad' values in temp map with the MEAN temperature of pixels in the clump
    #Run Stats and extract MEAN value of tempature
        cmd = '%s/stats ndf=%s %s'%(kapdir,ctemp,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'mean','stats') 
        status, output = commands.getstatusoutput(cmd)
        temp_mi = output
        temp_mi = float(temp_mi)

        if temp_mi > 0:
            temp_mi = temp_mi
        else:
            temp_mi = 15.0

        cmd = '%s/nomagic in=%s out=%s repval=%s > /dev/null '%(kapdir,ctemp,mtemp,temp_mi)
        os.system(cmd)

    #Make maps of mass
        m450 = mass_maps2.mass_map(c450,mtemp,kappa_4,d,i,450,'cores')
        m850 = mass_maps2.mass_map(c850,mtemp,kappa_8,d,i,850,'cores')
    
    #Get mass error. - based on mean temp. and final mass 

        dMi_450 = mass_err(c450i,temp_mi,0.3)
        dMi_850 = mass_err(c850i,temp_mi,0.02)

    #Extract core mass as number
        cmd = '%s/stats ndf=%s %s'%(kapdir,m450,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        m450i = output
        cmd = '%s/stats ndf=%s %s'%(kapdir,m850,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        m850i = output
        print 'mass of clump 850 =', m850i

    #Convert Ra/Dec from degrees to Time format.
        Ra = float(data[i-1][1])
        Dec = float(data[i-1][2])

        Ra_T = RaDecDegs.timeRa(Ra,'FALSE')
        Dec_T = RaDecDegs.timeDec(Dec,'FALSE')

    #Write to files #LaTex data file
        core_LTX.write(str(i)+'\t&\t')
        #core_LTX.write(Ra_T+'\t&\t')
        #core_LTX.write(Dec_T+'\t&\t')
        core_LTX.write(str(round(float(c450i),1))+'\t&\t')
        core_LTX.write(str(round(float(c850i),2))+'\t&\t')
        #core_LTX.write(str(round(float(m450i),1))+'('+str(round(float(dMi_450),1))+')\t&\t')
        core_LTX.write(str(round(float(m850i),2))+'('+str(round(float(dMi_850),2))+')\t&\t')
        core_LTX.write(str(round(float(temp_mi),1))+'('+str(round(float(0.05*temp_mi),1))+') \n')
