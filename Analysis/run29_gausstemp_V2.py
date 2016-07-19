#DJR UoE
#30/10/2015 - run29.py
#
#code is designed to calculate the gassian wieghted temp of a core. 

###########  Import modules ######################################

import os
import numpy as np
import string
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
import copy
import commands
from astropy import units as u
from astropy.coordinates import SkyCoord
import re

########### set CSH Commands directories #############

kapdir = '/stardev/bin/kappa'
convdir = '/stardev/bin/convert'

###################################################################
#Start of bulk script
###################################################################

def exp(x,y,FWHM,out,null):
    ln = 2.354820045 #2SQRT(2ln2)
    p = 3
    sig = 2*(((FWHM/p)/ln)**2.0)
    exp = "'exp(-(((xa+0.5"+x+")**2)+((xb+0.5"+y+")**2))/"+str(sig)+")'"

    #norm = 1/(sig*np.sqrt(2*np.pi))
    #print norm
    
    cmd = '%s/maths exp=%s out=%s like=%s'%(kapdir,exp,out,null)
    os.system(cmd)
    #cmd = '%s/cmult in=%s scalar=%f out=%s'%(kapdir,out,norm,Nout)
    #os.system(cmd) 
    return 

def guasmooth(x,y,fwhm):
    pixel = "'Gausstests/pixel.sdf'"

    #Implant point sources into null maps
    section = "\"'%i,%i'\""%(int(x),int(y))
    chpix1 = '%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,null,pixel,section,n)
    os.system(chpix1)  
    sig = fwhm/(2.*(2.*np.log(2))**0.5)
    smooth = sig/3
    cmd = '%s/gausmooth in=%s fwhm=%f out=%s'%(kapdir,pixel,smooth,out)
    os.system(cmd)
    return 

def fmask(maskmap,inmap,outmap):
    #create mask and mask input map
    maskbad = 'MSK.sdf'  
    cmd = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%s newhi=%d %s'%(kapdir,maskmap,maskbad,0,0,'bad',1,'> /dev/null')
    os.system(cmd)
    mul1 = '%s/mult in1=%s in2=%s out=%s'%(kapdir,maskbad,inmap,outmap)
    os.system(mul1)
    #print 'masked'
    return outmap

def Gauss_weight(X,Y,FWHM):
    #Produce gaussian 
    Gmask = 'Gmask.sdf'
    guasmooth(X,Y,FWHM)
    #exp(X,Y,FWHM,out,null)
    #mask gaussian
    fmask(temp,out,Gmask)

    #SUM gauss component
    cmd = '%s/stats ndf=%s %s'%(kapdir,Gmask,'> /dev/null')
    os.system(cmd)
    cmd = '%s/parget parname=%s applic=%s'%(kapdir,'total','stats') 
    status, output = commands.getstatusoutput(cmd)
    gw = float(output)

    #Record gaussian position 
    cmd = '%s/stats ndf=%s %s'%(kapdir,out,'> /dev/null')
    os.system(cmd)
    cmd = '%s/parget parname=%s applic=%s'%(kapdir,'MAXWCS','stats') 
    status, output = commands.getstatusoutput(cmd)
    ysoloc = str(output)
    ysoloc = string.split(ysoloc,',')
    YSOloc = ysoloc[0]+ysoloc[1]
    c = SkyCoord('%s'%(YSOloc), unit=(u.hourangle, u.deg))
    ra = c.ra.degree
    dec = c.dec.degree

    #if gaussian covers a real area, continue 
    if gw > 1E-2:
        #print X,Y
        ###TEMP
        Gtemp = 'Gtemp.sdf'   
        #RUN gaussian wieghting and mean calculation.
        os.system('%s/mult in1=%s in2=%s out=%s'%(kapdir,out,temp,Gtemp))
        os.system('%s/stats ndf=%s %s'%(kapdir,Gtemp,'> /dev/null'))
        status, output = commands.getstatusoutput('%s/parget parname=%s applic=%s'%(kapdir,'total','stats') )
        gwT = float(output)
        ### ERROR
        Gerr = 'Gerr.sdf'   
        #RUN gaussian wieghting and mean calculation.
        os.system('%s/mult in1=%s in2=%s out=%s'%(kapdir,out,temperr,Gerr))
        os.system('%s/stats ndf=%s %s'%(kapdir,Gerr,'> /dev/null'))
        status, output = commands.getstatusoutput('%s/parget parname=%s applic=%s'%(kapdir,'total','stats'))
        gwTe = float(output)
        #print 'gw:',round(gw,1)
        #print 'gwT:',(gwT)
                    
        T = gwT/gw
        dT = gwTe/gw
        #print 'T:',str(round(T,1))+'pm'+str(round(dT,1))+' K'
        #if X == '-5':
                    
        #Find central temperature
        #TEMP
        out2 = 'Gausstests/point.sdf'
        section = "\"'%i,%i'\""%(int(X),int(Y))
        #print chpix1
        os.system('%s/chpix in=%s out=%s section=%s newval=%f'%(kapdir,null,out2,section,1))  
        Ptemp = 'Gausstests/Ptemp.sdf'
        os.system('%s/mult in1=%s in2=%s out=%s'%(kapdir,out2,temp,Ptemp))
        os.system('%s/stats ndf=%s %s'%(kapdir,Ptemp,'> /dev/null'))
        status, output = commands.getstatusoutput('%s/parget parname=%s applic=%s'%(kapdir,'total','stats'))
        ptemp = float(output)
                    
        #ERROR
        Ptemperr = 'Gausstests/Ptemp_err.sdf'
        os.system('%s/mult in1=%s in2=%s out=%s'%(kapdir,out2,temperr,Ptemperr)) 
        os.system('%s/stats ndf=%s %s'%(kapdir,Ptemperr,'> /dev/null'))
        status, output = commands.getstatusoutput('%s/parget parname=%s applic=%s'%(kapdir,'total','stats') )
        ptemperr = float(output)

        #print 'T:',str(round(float(ptemp),1))+'pm'+str(round(float(ptemperr),1))+' K'

    else:
        T, dT, ptemp, ptemperr  = 0,0,0,0
    return T, dT, ptemp, ptemperr, ra, dec

def FWHM(core,d):
    FWHM = (((core*3.08567758E16)/1.49597871E11)/float(d))
    return FWHM

#Produce YSO map
maps = np.loadtxt('maps_gausstemp.txt',dtype='string')

#Open master results file for writing
table = open("GBSprotostar_temperaturesc.tab","w")
table = open("GBSprotostar_temperaturesc.tab","a")
table.write('#region\tra(deg)\tdec(deg)\tGra(deg)\tGdec(deg)\ttemp(K)\ttemperr(K)\tCentral_temp(K)\tCentral_temperr(K)\talphaE\tTbolE\n')

table2 = open("FWHMprotostar_temperaturesc.tab","w")
table2 = open("FWHMprotostar_temperaturesc.tab","a")
table2.write('#0.025pc\t0.050pc\t0.0750\t0.0100\t0.0125\t0.0150\n')

temp_core1 = []
temp_core2 = []
temp_core3 = []
temp_core4 = []
temp_core5 = []
temp_core6 = []

errtemp_core1 = []
errtemp_core2 = []
errtemp_core3 = []
errtemp_core4 = []
errtemp_core5 = []
errtemp_core6 = []

N = 0

out = 'out.sdf'

#loop through maps
for i in maps:
    #LOAD protostars
    pos = np.loadtxt('GBS_YSO_master_protostars.tab',dtype='string')

    #split INPUT file - #declare map origin
    s850,temp,temperr,raO,decO,d = string.split(i,',')
    prefix = string.split(s850,'.sdf')[0]
    file = string.split(prefix,'/')[1]
    #deffine region name
    name = string.split(s850,'/')[1]
    region = string.split(name,'_2')[0]
    print region
    
    MASS = 'CDmaps/mass/'+region+'_mass.sdf'

    #Produce null map
    null = 'null.sdf'
    thresh1 = '%s/thresh in=%s out=%s thrlo=%d thrhi=%d newlo=%d newhi=%d QUIET'%(kapdir,s850,null,0,0,0,0)
    os.system(thresh1)

    #deffine FWHM
    core = [0.025,0.050,0.075,0.100,0.125,0.15] #pc
    
    #print 'FWHM is '+str(FWHM)

    n = 1
    k = 0
    l = 0

    #deffine limits on region search area
    limitra,limitdec = 1,1

    ralo = float(raO)-float(limitra)
    rahi = float(raO)+float(limitra)
    declo = float(decO)-float(limitdec)
    dechi = float(decO)+float(limitdec)
    #print ralo, rahi, declo, dechi

    for j in pos:
        
        #limit catalogue to the specific region
        if (float(j[2]) > ralo) & (float(j[2]) < rahi) & (float(j[3]) > declo) & (float(j[3]) < dechi):
            if float(j[4]) > -0.3:

                #print 'INPUT:'+str(j[2])+','+str(j[3])

                #PIXEL COORD METHOD
                RaH = float(j[2])/15
                DecD = float(j[3])
                coords = "%10.7f %s"%(RaH, DecD)

                cmd = "%s/wcstran ndf=%s posin=\"\'%s\'\" framein=SKY frameout=PIXEL > /dev/null"%(kapdir,s850,coords)
                os.system(cmd)
                posout = (os.popen(kapdir+"/parget posout wcstran").readlines())
                (x,y,junk) = re.split("\s+",posout[0],2)
                
                X = str(int(round(float(x),1)))
                Y = str(int(round(float(y),1)))

                #print X,Y
                i = core[0]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                if TC > 0:
                    temp_core1.append(TC)
                    errtemp_core1.append(eTC)
                i = core[3]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                if TC > 0:
                    temp_core4.append(TC)
                    errtemp_core4.append(eTC)
                i = core[1]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                table.write(str(region)+'\t'+str(j[2])+'\t'+str(j[3])+'\t'+str(RA)+'\t'+str(DEC)+'\t'+str(TC)+'\t'+str(eTC)+'\t'+str(float(pT))+'\t'+str(float(epT))+'\t'+str(j[6])+'\t'+str(j[7])+'\n')
                if TC > 0:
                    temp_core2.append(TC)
                    errtemp_core2.append(eTC)
                    print str(round(TC,1))+'K'
                i = core[2]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                if TC > 0:
                    temp_core3.append(TC)
                    errtemp_core3.append(eTC)
                i = core[3]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                if TC > 0:
                    temp_core4.append(TC)
                    errtemp_core4.append(eTC)
                i = core[4]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                if TC > 0:
                    temp_core5.append(TC)
                    errtemp_core5.append(eTC)
                i = core[5]
                fwhm = FWHM(i,float(d)) #arcsecs
                TC,eTC,pT,epT,RA,DEC = Gauss_weight(X,Y,fwhm)
                if TC > 0:
                    temp_core6.append(TC)
                    errtemp_core6.append(eTC)

                    
    #print temp_core2
    print np.mean(temp_core2)
    table2.write(str(np.mean(temp_core1))+'\t'+str(np.mean(temp_core2))+'\t'+str(np.mean(temp_core3))+'\t'+str(np.mean(temp_core4))+'\t'+str(np.mean(temp_core5))+'\t'+str(np.mean(temp_core6)))

    break


