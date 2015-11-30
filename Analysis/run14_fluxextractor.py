#Damian Rumble, UoE
#22/11/2013
#run14_fluxextractor.py

#This code extracts total flux of a clump, as identified by FellWalker, and the temperature for the same region and stores it in a file for calculation as a mass.

#Flux typically has more pixels than temp. Where there is overlap, bad temp pixels are labelled with the mean temp. of the remaining group. Where no temp. data exists, pixels are labelled with a constant temperature of 15K 

#Jeans Mass is calculated in parrallel 

#All data is written to a file in a suitable format for Latex Tables

#####################################################
#import maths, ploting and astrophysical packages

import numpy as np
import os
import subprocess
import commands
import pyfits

#my modules
import mass_maps2
import jeans
import RaDecDegs
import mass


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


########### set CSH Commands directories #############

kapdir = '/star/bin/kappa'
convdir = '/star/bin/convert'

########### Functions #############

def flux_450(mask,s450,c450):
    mul1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,s450,c450)
    os.system(mul1)
    print 'Flux 450 masked'
    return c450
    
def flux_850(mask,s850,c850):
    mul2 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,s850,c850)
    os.system(mul2)
    print 'Flux 450 masked'
    return c850

def temp_45(mask,temp,out):
    T = '%s/mult in1=%s in2=%s out=%s'%(kapdir,mask,temp,out)
    os.system(T)
    print 'Temp. masked'
    return out

def mass_err(S,T,dS): #Mass error
    S = float(S)
    A = (0.39**2.)*(((np.exp(17./T)-1)**2.)*(dS**2.))
    B = (0.39**2.)*(289.*(S**2.)/(T**4.))*np.exp(34./T)*((0.05*T)**2.)
    
    M = A + B
    dM = np.sqrt(M)

    return dM

########### Bulk Code #############

clumps = 'SerpensMWC297_IR1_clumps3.sdf'

s450 = 'SerpensMWC297_20140325_IR1_s450_freefree.sdf'
s850 = 'SerpensMWC297_20140325_IR1_s850_freefree.sdf'
temp = 'SerpensMWC297_20140325_IR1_freefree450850temp%5.sdf'

#For 450 clumps, the 450 map needs to be regridded onto 850 scale to match the clumps (also 850 scale). Then the total flux devided by 2.25 to regain 450 levels. This requires collapsing first.

d2_450 = 'temp/s450_2d.sdf'
d2_850 = 'temp/s850_2d.sdf'
align = 'temp/SerpensMWC297_20130916_s450_aign_850.sdf'

cmd = '%s/collapse in=%s out=%s axis=%s %s'%(kapdir,s450,d2_450,'3','> /dev/null')
os.system(cmd)
cmd = '%s/collapse in=%s out=%s axis=%s %s'%(kapdir,s850,d2_850,'3','> /dev/null')
os.system(cmd)
cmd = '%s/wcsalign in=%s out=%s ref=%s method=%s conserve=%s %s'%(kapdir,d2_450,align,d2_850,'nearest','TRUE','accept')
os.system(cmd)

params = 'SerpensMWC297_IR1_clumps3_params.fit'

data = pyfits.getdata(params)

flux450 = []
flux850 = []

i = 1

#open files for writing
SMM = open("SMM.txt","w")
SMM = open("SMM.txt","a")

SMM_LTX = open("SMM_LTX.txt","w")
SMM_LTX = open("SMM_LTX.txt","a")

SMM.write('#index\tRa\tDec\tFlux 850\tMass 850\tMean Temp.\tPixel#\teffRad(pc)\tErr_effRad(pc)\tJeans Mass\terr_Jeans Mass\tM/Mj\terrR\n')
SMM_LTX.write('index\t&Ra\t&Dec\t&\tFlux 850\t&\tMass 850\t&\tMean Temp.\t&\tPixel#\t&\teffRad(au)\t&\tJeans Mass\t&\terr_Jeans Mass\tM/Mj\n')

for i in range(1,30):

    print 'Starting processing of clump#'+str(i)

    if i == 24 or i == 26 or i == 27 :
        print 'Not a real clump'
    else:

    #Set up temp. file names
        thresh = 'temp/thresh'+str(i)+'.sdf'
        magic = 'temp/magic'+str(i)+'.sdf'
        mask = 'temp/mask'+str(i)+'.sdf'
        maskbad = 'temp/maskbad'+str(i)+'.sdf'
        mtemp = 'temp/meantemp'+str(i)+'.sdf'
        ctempbad = 'temp/temp_clump_bad'+str(i)+'.sdf'

    #set up clump file names
        c450 = 'clumps/450/SerpensMWC297_s450_clump'+str(i)+'.sdf'
        c850 = 'clumps/850/SerpensMWC297_s850_clump'+str(i)+'.sdf'
        ctemp = 'clumps/temp/SerpensMWC297_temp_clump'+str(i)+'.sdf'
        masktemp = 'clumps/temp/SerpensMWC297_MMtemp_clump'+str(i)+'.sdf'
    
    #process 1 - isolate clump i
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,clumps,thresh,i,i,0,0,'> /dev/null')
        os.system(cmd)

    #process 2 - remove 0value data
        cmd = '%s/setmagic in=%s out=%s repval=%d %s'%(kapdir,thresh,magic,0,'> /dev/null')
        os.system(cmd)

    #recreate mask of clump i
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d %s'%(kapdir,magic,mask,0,0,0,1,'> /dev/null')
        os.system(cmd)

    #recreate mask of clump i
        cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%d %s'%(kapdir,magic,maskbad,0,0,'bad',1,'> /dev/null')
        os.system(cmd)

   #masking orginal maps, new saved to temp files
        flux_450(mask,align,c450)
        flux_850(mask,s850,c850)
        temp_45(maskbad,temp,ctempbad)
        temp_45(mask,temp,ctemp)
        
    #FLUXES
        #Extract total fluxes
        cmd = '%s/stats ndf=%s %s'%(kapdir,c450,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c450i = output 

        cmd = '%s/stats ndf=%s %s'%(kapdir,c850,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        c850i = output
        c850i = float(c850i)

        print '850 flux = ', c850i

    #TEMPERATURE
        #replace 'bad' values in temp map with the MEAN temperature of pixels in the clump
        #Run Stats and extract MEAN value of tempature
        cmd = '%s/stats ndf=%s %s'%(kapdir,ctempbad,'> /dev/null')
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

        temp_45(mask,mtemp,masktemp)

        print 'mean temperature is ', temp_mi, 'K'

    #MASS - 
        #calculate masses
        mass450 = mass_maps2.mass_map(c450,masktemp,kappa_4,d,i,450,'clumps')
        mass850 = mass_maps2.mass_map(c850,masktemp,kappa_8,d,i,850,'clumps')

        #Get mass error. - based on mean temp. and final mass 
        dMi_450 = mass_err(c450i,temp_mi,0.3)
        dMi_850 = mass_err(c850i,temp_mi,0.02)

        #Extract clump mass as number
        cmd = '%s/stats ndf=%s %s'%(kapdir,mass450,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        m450i = output
        cmd = '%s/stats ndf=%s %s'%(kapdir,mass850,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
        status, output = commands.getstatusoutput(cmd)
        m850i = output
  
        print '850 mass by mean = ',mass.djr(float(c850i),float(temp_mi),kappa_4,d)
        print '850 mass by pixel = ',m850i

    #Run Stats and extract NUMBER of pixels
        cmd = '%s/stats ndf=%s %s'%(kapdir,magic,'> /dev/null')
        os.system(cmd)
        cmd = '%s/parget parname=%s applic=%s'%(kapdir,'numgood','stats') 
        status, output = commands.getstatusoutput(cmd)
        Ni = output

    #Calculate effetive radii and append to file
        Ni = float(Ni)

        xxyy = a*d*(au/psc)
        A = Ni*(xxyy**2.0)
        ri = (A**0.5)/pi
        err_ri = ri*(errd/d)

        #Ni = float(Ni)
        #ri = ((Ni**(0.5))*5.06E+11*d)/au #au
        #err_ri = (5.06E+11 * (errd) * (Ni**(0.5)))/au #au 

    #run module for jeans mass
        temp_mi = float(temp_mi)
        Ji = jeans.sadavoy(temp_mi,d,Ni)
        err_Ji = Ji * (0.05) #frac. error in temp. 
        if float(Ji) == 0:
            equ = 0
        else:
            equ = float(m850i) / float(Ji)

    #Error in mass ratio:
        frac = ((dMi_850/float(m850i))**2.)+((err_Ji/float(Ji))**2.)
        err_R = equ * ((frac)**0.5)

    #Convert Ra/Dec from degrees to Time format.
        Ra = float(data[i-1][1])
        Dec = float(data[i-1][2])

        Ra_T = RaDecDegs.timeRa(Ra,'FALSE')
        Dec_T = RaDecDegs.timeDec(Dec,'FALSE')

    #Write to files
    #open files for appending
        SMM = open("SMM.txt","a")
        SMM_LTX = open("SMM_LTX.txt","a")

    #RAW data file
        SMM.write(str(i)+'\t')
        SMM.write(str(round(float(data[i-1][1]),2))+'\t')
        SMM.write(str(round(float(data[i-1][2]),4))+'\t')
        SMM.write(str(c450i)+'\t')
        SMM.write(str(c850i)+'\t')
        SMM.write(str(m450i)+'\t')
        SMM.write(str(m850i)+'\t')
        SMM.write(str(dMi_850)+'\t')
        SMM.write(str(temp_mi)+'\t')
        SMM.write(str(0.05*temp_mi)+'\t')
        SMM.write(str(Ni)+'\t')
        SMM.write(str(ri)+'\t')
        SMM.write(str(err_ri)+'\t')
        SMM.write(str(Ji)+'\t')
        SMM.write(str(err_Ji)+'\t')
        SMM.write(str(equ)+'\t')
        SMM.write(str(err_R)+'\n')

    #LaTex data file
        SMM_LTX.write(str(i)+'\t&\t')
        SMM_LTX.write(Ra_T+'\t&\t')
        SMM_LTX.write(Dec_T+'\t&\t')
        SMM_LTX.write(str(round(float(c450i),1))+'\t&\t')
        SMM_LTX.write(str(round(float(c850i),1))+'\t&\t')
    #SMM_LTX.write(str(round(float(m450i),1))+'\t&\t')
        SMM_LTX.write(str(round(float(m850i),2))+'('+str(round(float(dMi_850),2))+')\t&\t')
        SMM_LTX.write(str(round(float(temp_mi),1))+'('+str(round(float(0.05*temp_mi),1))+')\t&\t')
        SMM_LTX.write(str(round(float(Ni),-0))+'\t&\t')
        #SMM_LTX.write(str(round(float(ri),2))+'\t&\t')
        SMM_LTX.write(str(round(float(Ji),2))+'('+(str(round(float(err_Ji),2))+')\t&\t'))
        SMM_LTX.write(str(round(equ,2))+'('+(str(round(float(err_R),2)))+")\\\n")

    #increment
    i = i + 1
    
print 'complete'
